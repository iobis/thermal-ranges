get_thermal_ranges <- function(
    occ_source,
    occ_path = NULL,
    species,
    aoi,
    output_dir,
    output_filename,
    output_format,
    future_temperature,
    future_sites,
    future_mode,
    future_period,
    future_scenarios
) {

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
    if (is.numeric(species)) {
        species <- data.frame(aphiaID = species, species = worrms::wm_id2name(species))
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
            species <- species |> select(aphiaID = aphiaid, species)
        }
    } else if (tools::file_ext(species) %in% c("csv", "tsv", "txt")) {
        species <- data.table::fread(species)
        colnames(files) <- tolower(colnames(files))
        if (any(!c("species", "aphiaid") %in% colnames(files))) {
            cli::cli_abort("Supplied table does not contain a field `species` or `aphiaid`")
        }
        species <- species |> select(aphiaID = aphiaid, species)
    } else {
        cli::cli_abort("Could not load species list")
    }

    # Check if AOI is not null
    if (!is.null(aoi)) {
        if (tools::file_ext(aoi) != "") {
            aoi <- aoi |>
                sf::read_sf() |>
                sf::st_as_sfc() |>
                sf::st_as_text()
        }
        if (nchar(aoi) > 1500) {
            cli::cli_alert_warning("`aoi` is very complex, it may fail!")
        }
    }

    # Get species data
    if (occ_source == "api") {
        occ_data <- robis::occurrence(taxonid = species$taxonID)

        if (fast_mode) {
            colnames(occ_data)[colnames(occ_data) == "sst"] <- "temperature"
            return(.prepare_tr_table(occ_data))
        }
    } else {
        con <- DBI::dbConnect(duckdb::duckdb())
        DBI::dbSendQuery(con, "install h3 from community; load h3;")
        DBI::dbSendQuery(con, "install spatial; load spatial;")

        if (occ_source == "speciesgrids") {
            if (tools::file_ext(occ_path) != "parquet") occ_path <- file.path(occ_path, "*")
            if (!is.null(aoi)) {
                occ_data <- DBI::dbGetQuery(con, glue::glue(
                    "select species, AphiaID, records, cell
                    from read_parquet('{occ_path}')
                    where AphiaID in ({paste(species, collapse = ',')}) and ST_Intersects (geometry, '{aoi}')
                    "
                ))
            } else {
                occ_data <- DBI::dbGetQuery(con, glue::glue(
                    "select species, AphiaID, records, cell
                    from read_parquet('{occ_path}')
                    where AphiaID in ({paste(species, collapse = ',')})
                    "
                ))
            }
            occ_data_coords <- h3jsr::cell_to_point(occ_data$cell) |>
                sf::st_coordinates() |>
                as.data.frame() |>
                rename(longitude = X, latitude = Y)
            occ_data <- bind_cols(occ_data, occ_data_coords)
        } else if (occ_source == "full") {
            # TODO
        }

        DBI::dbDisconnect(con)

        temp <- .load_temperature()

        occ_data$temperature <- terra::extract(temp, occ_data[,c("longitude", "latitude")], ID = F)[,1]

        occ_data <- occ_data |>
            select(aphiaID = AphiaID, species, temperature)

        return(.prepare_tr_table(occ_data))

    }

}

r <- data.table::fread("https://github.com/iobis/whomp-species-lists/raw/refs/heads/main/output/biosphere_marine.tsv")
head(r)


.prepare_tr_table <- function(occ_data) {
    occ_data |>
        group_by(aphiaID, species) |>
        summarise(
            q01 = quantile(temperature, 0.01, na.rm = TRUE),
            q05 = quantile(temperature, 0.05, na.rm = TRUE),
            q50 = quantile(temperature, 0.5, na.rm = TRUE),
            q95 = quantile(temperature, 0.95, na.rm = TRUE),
            q99 = quantile(temperature, 0.99, na.rm = TRUE),
            mean = mean(temperature, na.rm = TRUE),
            sd = sd(temperature, na.rm = TRUE),
            na_count = sum(is.na(temperature))
        )
}


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


.get_nearby <- function() {

}