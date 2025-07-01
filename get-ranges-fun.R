#' Get thermal ranges for species based on occurrence data
#'
#' @param occ_source one of 'api' (to use the OBIS API), 'full' (when using the
#' OBIS full export) or 'speciesgrids' (to use the OBIS speciesgrids product;
#' recommended). Note that 'api' is intended for tests only, and is limited to
#' 10 species.
#' @param occ_path in case of `occ_source` 'full' or 'speciesgrids', the path
#' to the parquet dataset
#' @param species 'all', to run for all species, an AphiaID (numeric), 
#' a csv/txt/tsv file with species list, or a folder containing csv/txt/tsv files. 
#' For the lists, it should contain columns aphiaID and species (not case sensitive)
#' You can also set this to NULL. In that case, you should pass an `aoi` to get
#' all species within the area of interest.
#' @param min_records minimum records available to proceed and calculate the thermal ranges
#' @param min_year minimum year to get the records
#' @param aoi a WKT string or a path to a shapefile containing an area of interest
#' at which to get the species list for calculating the thermal ranges
#' @param output_dir the output directory. If non existent, it will be created.
#' @param output_filename the name of the output file
#' @param output_format the output format. One of 'csv' (default), 'txt', 'tsv' or 'parquet'
#'
#' @return
#' @export
#'
#' @examples
get_thermal_ranges <- function(
    occ_source,
    occ_path = NULL,
    species,
    min_records = 10,
    min_year = 1950,
    aoi = NULL,
    output_dir = "results",
    output_filename = paste0("thermal_ranges_", format(Sys.Date(), "%Y%m%d")),
    output_format = "csv"
) {

    require(dplyr)

    # Initial checking
    if (!occ_source %in% c("api", "speciesgrids", "full")) {
        cli::cli_abort('{.var occ_source} should be one of `api`, `speciesgrids` or`full`')
    } else if (occ_source %in% c("speciesgrids", "full") && is.null(occ_path)) {
        cli::cli_abort("For `speciesgrids` or `full` you need to pass {.var occ_path}")
    }
    if (!output_format %in% c("csv", "parquet", "txt", "tsv")) {
        cli::cli_abort('{.var output_format} should be one of {.arg {c("csv", "parquet", "txt", "tsv")}}')
    }

    # Load species list
    if (is.null(species)) {
        if (is.null(aoi)) {
            cli::cli_abort("If `species` is `NULL` you should pass an `aoi`")
        } else {
            species <- .retrieve_names_from(
                occ_source, occ_path, min_records, aoi
            )
        }
    } else if (is.numeric(species)) {
        species <- data.frame(aphiaID = species, species = worrms::wm_id2name(species), group = 1)
    } else if (length(species) == 1 && species == "all") {
        cli::cli_alert_info("Calculating thermal ranges for all available species")
        species <- .retrieve_names_from(occ_source, occ_path, min_records)
    } else if (tools::file_ext(species) == "") {
        if (!dir.exists(species) || length(list.files(species)) < 1) {
            cli::cli_abort("Directory passed to `species` does not exists or is empty.")
        } else {
            files <- list.files(species, full.names = T) |>
                        lapply(data.table::fread) |>
                        (\(x){do.call("rbind", x)})()
            colnames(files) <- tolower(colnames(files))
            if (any(!c("species", "aphiaid") %in% colnames(files))) {
                cli::cli_abort("Supplied table does not contain a field `species` or `aphiaid`")
            }
            species <- files |> 
                select(aphiaID = aphiaid, species) |>
                distinct(aphiaID, .keep_all = TRUE)
            counts <- .get_counts(species$aphiaID, occ_source, occ_path)
            species <- left_join(species, counts, by = "aphiaID")
            species <- .add_group(species)
        }
    } else if (tools::file_ext(species) %in% c("csv", "tsv", "txt")) {
        species <- data.table::fread(species)
        colnames(files) <- tolower(colnames(files))
        if (any(!c("species", "aphiaid") %in% colnames(files))) {
            cli::cli_abort("Supplied table does not contain a field `species` or `aphiaid`")
        }
        species <- files |> 
            select(aphiaID = aphiaid, species) |>
            distinct(aphiaID, .keep_all = TRUE)
        counts <- .get_counts(species$aphiaID, occ_source, occ_path)
        species <- left_join(species, counts, by = "aphiaID")
        species <- .add_group(species)
    } else {
        cli::cli_abort("Could not load species list")
    }

    # Get temperature data
    sp_temps <- .get_temperature(species, occ_source, occ_path, min_year)

    # Save results
    fs::dir_create(output_dir)
    output_filename <- gsub("\\..*", "", output_filename)
    output_filename <- switch(output_format,
        csv = file.path(output_dir, paste0(output_filename, ".csv")),
        txt = file.path(output_dir, paste0(output_filename, ".txt")),
        tsv = file.path(output_dir, paste0(output_filename, ".tsv")),
        parquet = file.path(output_dir, paste0(output_filename, ".parquet"))
    )
    switch(output_format,
        csv = write.csv(sp_temps, output_filename, row.names = F),
        txt = write.table(sp_temps, output_filename, row.names = F),
        tsv = write.table(sp_temps, output_filename, row.names = F, quote = FALSE, sep = "\t"),
        parquet = arrow::write_parquet(sp_temps, output_filename)
    )

    cli::cli_alert_success(
        "Thermal ranges retrieved and file saved: {.path {output_filename}}"
    )

    return(invisible(NULL))

}


# Extract temperature data
#
# @param species species data frame containing aphiaID, species and group
# @param occ_source where to get occurrence data
# @param occ_path path, if `occ_source` is 'speciesgrids' or 'full'
# @param mult_mode if there are NA values, it will try to get from 8 adjacent cells.
# @param min_year minimum year to get the records
# If mult_mode is distance and more than one adjacent point is valid,
# than it get the closest of the valid points. Otherwise, it samples one point.
#
# @return data.frame
#' @export
.get_temperature <- function(species, occ_source, occ_path, min_year = 1950, mult_mode = "distance") {

    cli::cli_alert_info("Loading environmental layers")
    sst <- .load_temperature()

    if (occ_source == "speciesgrids") {
        if (tools::file_ext(occ_path) != "parquet") occ_path <- file.path(occ_path, "*")
    }

    cli::cli_alert_info("Extracting temperature data")
    total <- length(unique(species$group))
    if (total > 2) {
        pb <- progress::progress_bar$new(total = total)
    }
    results <- vector("list", total)
    for (id in unique(species$group)) {
        if (total > 2) pb$tick()

        group_species <- species$aphiaID[species$group == id]

        if (occ_source == "api") {
            occurrences <- robis::occurrence(
                taxonid = group_species,
                startdate = as.Date(paste0(min_year, "-01-01"))
            )
            occurrences <- occurrences |>
                select(species = scientificName, AphiaID = aphiaID, decimalLongitude, decimalLatitude) |>
                as.data.frame()
        } else if (occ_source == "full") {
            # TODO
        } else if (occ_source == "speciesgrids") {
            con <- DBI::dbConnect(duckdb::duckdb())
            DBI::dbSendQuery(con, "install h3 from community; load h3;")
            occurrences <- DBI::dbGetQuery(con, glue::glue(
                "select species, AphiaID, h3_cell_to_lng(cell) as lng, h3_cell_to_lat(cell) as lat
                from read_parquet('{occ_path}')
                where AphiaID in ({paste(group_species, collapse = ',')}) and min_year >= {min_year}
                "
            ))
            DBI::dbDisconnect(con)
            occurrences <- occurrences |>
                select(species, AphiaID, decimalLongitude = lng, decimalLatitude = lat)
        }

        extracted_sst <- terra::extract(sst, occurrences[,c("decimalLongitude", "decimalLatitude")], 
                                        ID = FALSE, cell = T)

        unique_cells <- extracted_sst |>
            distinct(cell, .keep_all = T)

        # Fill NAssummarise
        is_na <- which(is.na(unique_cells$surface))
        if (length(is_na) > 0) {

            na_done <- data.frame(cell = 0, new_cell = 0)[0,]
            if (exists("pre_proc")) {
                na_done <- pre_proc[pre_proc$cell %in% unique_cells$cell[is_na], ]
                is_na <- is_na[!unique_cells$cell[is_na] %in% pre_proc$cell]
            }

            if (length(is_na) > 0) {
                na_adjs <- terra::adjacent(sst$surface, unique_cells$cell[is_na], directions = "queen")
                na_adjs <- as.vector(t(na_adjs))
                na_adjs_v <- terra::extract(sst$surface, na_adjs)[,1]
                na_adjs_v <- matrix(data = na_adjs_v, ncol = 8, nrow = length(is_na), byrow = T)
                na_adjs <- matrix(data = na_adjs, ncol = 8, nrow = length(is_na), byrow = T)

                unique_new <- lapply(seq_len(nrow(na_adjs_v)), \(nr){
                    adj_v <- na_adjs_v[nr,]
                    if (sum(!is.na(adj_v)) == 1) {
                        new_cell <- na_adjs[nr,][!is.na(adj_v)]
                        cbind(cell = unique_cells$cell[is_na[nr]], new_cell = new_cell)
                    } else if (sum(!is.na(adj_v)) > 1) {
                        new_cell <- na_adjs[nr,][!is.na(adj_v)]
                        if (mult_mode == "distance") {
                            dists <- terra::distance(
                                terra::vect(terra::xyFromCell(sst, unique_cells$cell[is_na[nr]]), crs = "EPSG:4326"),
                                terra::vect(terra::xyFromCell(sst, new_cell), crs = "EPSG:4326")
                            )
                            cbind(
                                cell = unique_cells$cell[is_na[nr]],
                                new_cell = new_cell[which.min(dists)]
                            )
                        } else {
                            new_cell <- sample(new_cell, 1)
                            cbind(
                                cell = unique_cells$cell[is_na[nr]],
                                new_cell = new_cell
                            )
                        }
                    } else {
                        cbind(cell = unique_cells$cell[is_na[nr]], new_cell = unique_cells$cell[is_na[nr]])
                    }
                })
                unique_new <- do.call("rbind", unique_new)
                unique_new <- as.data.frame(unique_new) |>
                    bind_rows(na_done)

                # Input new cells values on old
                unique_new <- unique_new |> 
                    distinct(cell, .keep_all = T)
                if (exists("pre_proc")) {
                    pre_proc <- pre_proc |>
                        bind_rows(unique_new) |>
                        distinct(cell, .keep_all = T)
                } else {
                    pre_proc <- unique_new
                }
            } else {
                unique_new <- na_done
            }

            extracted_sst <- dplyr::left_join(extracted_sst, unique_new, by = "cell")

            extracted_sst$new_cell[is.na(extracted_sst$new_cell)] <- extracted_sst$cell[is.na(extracted_sst$new_cell)]

            extracted_sst <- terra::extract(sst, extracted_sst$new_cell)
        }

        occurrences <- dplyr::bind_cols(occurrences, extracted_sst)

        surface <- occurrences |>
            select(aphiaID = AphiaID, species, temperature = surface) |>
            .prepare_tr_table() |>
            .add_confidence(layer = sst, target = "surface")
        colnames(surface)[2:length(surface)] <- paste0("surface_", colnames(surface)[2:length(surface)])

        bottom <- occurrences |>
            select(aphiaID = AphiaID, species, temperature = bottom) |>
            .prepare_tr_table() |>
            .add_confidence(layer = sst, target = "bottom")
        colnames(bottom)[2:length(bottom)] <- paste0("bottom_", colnames(bottom)[2:length(bottom)])

        unique_aphias <- species |>
            filter(group == id) |>
            select(aphiaID, species)

        results[[id]] <- unique_aphias |>
            left_join(surface, by = "aphiaID") |>
            left_join(bottom, by = "aphiaID")

    }

    return(bind_rows(results))
}


# Calculate thermal range
#
# @param occ_data the data.frame returned by .get_temperature 
#
# @return data.frame with summaries
#' @export
.prepare_tr_table <- function(occ_data) {
    occ_data |>
        group_by(aphiaID) |>
        summarise(
            min = ifelse(all(is.na(temperature)), NA, min(temperature, na.rm = TRUE)),
            q01 = quantile(temperature, 0.01, na.rm = TRUE),
            q05 = quantile(temperature, 0.05, na.rm = TRUE),
            q50 = quantile(temperature, 0.5, na.rm = TRUE),
            q95 = quantile(temperature, 0.95, na.rm = TRUE),
            q99 = quantile(temperature, 0.99, na.rm = TRUE),
            max = ifelse(all(is.na(temperature)), NA, max(temperature, na.rm = TRUE)),
            mean = ifelse(all(is.na(temperature)), NA, mean(temperature, na.rm = TRUE)),
            sd = sd(temperature, na.rm = TRUE),
            na_count = sum(is.na(temperature)),
            valid_count = sum(!is.na(temperature)),
            .groups = "drop"
        )
}


# Load temperature data from Bio-ORACLE
#
# @param layers which layers to load (surface, bottom)
#
# @return
#' @export
.load_temperature <- function(layers = c("surface", "bottom")) {
    if ("surface" %in% layers) {
        surface <- biooracler::download_layers(
            dataset_id = c("thetao_baseline_2000_2019_depthsurf"),
            constraints = list(time = c('2001-01-01T00:00:00Z', '2010-01-01T00:00:00Z')),
            variables = c("thetao_mean")
        )
        surface <- terra::mean(surface)
        names(surface) <- "surface"
    }
    if ("bottom" %in% layers) {
        bottom <- biooracler::download_layers(
            dataset_id = c("thetao_baseline_2000_2019_depthmean"),
            constraints = list(time = c('2001-01-01T00:00:00Z', '2010-01-01T00:00:00Z')),
            variables = c("thetao_mean")
        )
        bottom <- terra::mean(bottom)
        names(bottom) <- "bottom"
    }
    if (length(layers) > 1) {
        return(c(surface, bottom))
    } else if (layers == "surface") {
        return(surface)
    } else {
        return(bottom)
    }
}


# Retrieve species list
#
# @param occ_source the source to retrieve the list
# @param occ_path the path, if occ_source is 'speciesgrids' or 'full'
# @param min_records minimum number of records
# @param aoi an optional area of interest
#
# @return species list
#' @export
.retrieve_names_from <- function(occ_source, occ_path, min_records = 10, aoi = NULL) {

    if (occ_source == "api" & is.null(aoi)) {
        cli::cli_abort("When `occ_source = 'api'` you should supply an `aoi` to get a species list.")
    }

    if (occ_source == "api") {
        species_list <- robis::checklist(
            geometry = aoi
        )
        species_list <- species_list |>
            filter(records >= min_records) |>
            select(aphiaID = taxonID, species = scientificName, records) |>
            distinct(aphiaID, .keep_all = T)
    } else if (occ_source == "full") {
        # TODO
    } else if (occ_source == "speciesgrids") {
        if (tools::file_ext(occ_path) != "parquet") occ_path <- file.path(occ_path, "*")
        con <- DBI::dbConnect(duckdb::duckdb())
        DBI::dbSendQuery(con, "install spatial; load spatial;")
        occ_data <- DBI::dbGetQuery(con, ifelse(
            is.null(aoi),
            glue::glue(
                    "select species, AphiaID, count(*) as total
                    from read_parquet('{occ_path}')
                    group by species, AphiaID
                    having total >= {min_records}
                    "
                ),
            glue::glue(
                    "select distinct species, AphiaID, count(*) as total
                    from read_parquet('{occ_path}')
                    where ST_Intersects (geometry, '{aoi}')
                    group by species, AphiaID
                    having total >= {min_records}
                    "
                )
        ))
        DBI::dbDisconnect(con)

        species_list <- occ_data |>
            select(aphiaID = AphiaID, species, records = total) |>
            distinct(aphiaID, .keep_all = T)
    } else {
        cli::cli_abort("`occ_source` invalid")
    }

    species_list <- species_list[order(-species_list$records), ]

    cum_sum <- cumsum(species_list$records)

    batch_id <- cumsum(c(1, diff(cum_sum %/% 100000)) > 0)

    species_list$group <- batch_id

    return(species_list)

}


# Get number of records
#
# @param aphiaid a vector of AphiaIDs
# @param occ_source the source
# @param occ_path the path, if occ_source is 'speciesgrids' or 'full'
#
# @return the list with records
#' @export
.get_counts <- function(aphiaid, occ_source, occ_path) {
    if (occ_source == "api") {
        robis::checklist(taxonid = aphiaid) |>
            select(aphiaID = taxonID, species = scientificName, records)
    } else if (occ_source == "full") {
        # TODO
    } else if (occ_source == "speciesgrids") {
        if (tools::file_ext(occ_path) != "parquet") occ_path <- file.path(occ_path, "*")
        con <- DBI::dbConnect(duckdb::duckdb())
        occurrences <- DBI::dbGetQuery(con, glue::glue(
            "select AphiaID as aphiaID, species, count(*) as records
            from read_parquet('{occ_path}')
            where AphiaID in ({paste(aphiaid, collapse = ',')})
            group by AphiaID, species
            "
        ))
        DBI::dbDisconnect(con)
        return(occurrences)
    } else {
        cli::cli_abort("Invalid `occ_source`")
    }
}


# Divide list in subgroups (batches) with approximately 100,000 records
#
# @param species_list the species list
#
# @return
#' @export
.add_group <- function(species_list) {
    species_list <- species_list[order(-species_list$records), ]
    cum_sum <- cumsum(species_list$records)
    batch_id <- cumsum(c(1, diff(cum_sum %/% 100000)) > 0)
    species_list$group <- batch_id
    return(species_list)
}


# Add a confidence note
#
# @param species_list the species list
# @param layer SST SpatRaster
# @param target to which layer to get confidence
#
# @return
#' @export
.add_confidence <- function(species_list, layer, target) {
    layer <- terra::minmax(layer[[target]])
    layer <- layer["max",]

    species_list <- species_list |>
        mutate(
            upper_lim = layer - max
        ) |>
        mutate(
            confidence = ifelse(upper_lim == 0, "close-to-limit-0",
                ifelse(upper_lim <= 0.5, "close-to-limit-0.5",
                    ifelse(upper_lim <= 1, "close-to-limit-1", "ok")
                )
            )
        ) |>
        select(-upper_lim)

    return(species_list)
}
