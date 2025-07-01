get_temperature_sites <- function(shapefiles,
                                  output_folder,
                                  periods = NULL,
                                  scenarios = NULL) {

    if (length(shapefiles) == 1) {
        if (tools::file_ext(shapefiles) == "" & dir.exists(shapefiles)) {
            shapefiles <- list.files(shapefiles, full.names = T)
            shapefiles <- shapefiles[grepl("shp$|gpkg$|json$", shapefiles)]
        }
    }

    cli::cli_alert_info("Processing {length(shapefiles)} shapefile{?s}")

    sst <- .load_temperature_f(periods = periods, scenarios = scenarios)

    for (sh in seq_along(shapefiles)) {
        sel_sh <- shapefiles[sh]
        cli::cli_progress_step("Processing {basename(sel_sh)}")

        geoms <- terra::vect(sel_sh)

        geoms <- terra::project(geoms, "EPSG:4326")

        sst_sites <- terra::extract(sst, geoms, fun = mean, na.rm = TRUE)

        geoms <- cbind(geoms, sst_sites[,-1])

        geoms <- terra::as.data.frame(geoms)

        write.csv(geoms, 
                  file.path(output_folder, 
                            paste0(gsub("\\.shp|\\.gpkg|\\.json", "", basename(sel_sh)), "_exttemp.csv")), row.names = F)

    }

    cli::cli_progress_done()
    cli::cli_alert_success("All files saved at folder {.path {output_folder}}")
}



# Load temperature data from Bio-ORACLE including future
#
# @param layers which layers to load (surface, bottom)
# @param periods decade or decades for which to get future conditions. It always use a period of 20 years.
# Thus, if you supply 2050, it will actually get the period of 2030-2050.
# @param scenarios which SSP scenarios to use. You should write in the format "sspXXX". Available: "ssp119", "ssp126",
# "ssp245", "ssp370", "ssp460", and "ssp585"
#
# @return
#' @export
.load_temperature_f <- function(layers = c("surface", "bottom"), periods = NULL, scenarios = NULL) {

    if (!is.null(periods)) {
        if (is.null(scenarios)) {
            cli::cli_abort("When `periods` is supplied, you should also supply `scenarios`")
        } else {
            do_future <- TRUE
            future_matrix <- expand.grid(periods = periods, scenarios = scenarios)
        }
    } else {
        do_future <- FALSE
    }

    if ("surface" %in% layers) {
        surface <- biooracler::download_layers(
            dataset_id = c("thetao_baseline_2000_2019_depthsurf"),
            constraints = list(time = c('2001-01-01T00:00:00Z', '2010-01-01T00:00:00Z')),
            variables = c("thetao_mean")
        )
        surface <- terra::mean(surface)
        names(surface) <- "surface"

        if (do_future) {
            for (sc in seq_len(nrow(future_matrix))) {
                sel_p <- future_matrix$periods[sc]
                sel_s <- future_matrix$scenarios[sc]

                if (sel_p < 2030) {
                    cli::cli_abort("Minimum period date is 2030")
                } else if (sel_p == 2030) {
                    time_const <- c('2020-01-01T00:00:00Z', '2020-01-01T00:00:00Z')
                } else {
                    time_const <- c(paste0(sel_p - 20, '-01-01T00:00:00Z'), paste0(sel_p - 10, '-01-01T00:00:00Z'))
                }

                temp <- biooracler::download_layers(
                    dataset_id = paste0("thetao_", sel_s, "_2020_2100_depthsurf"),
                    constraints = list(time = time_const),
                    variables = c("thetao_mean")
                )
                temp <- terra::mean(temp)
                names(temp) <- paste0("surface_", sel_s, "_", sel_p)

                surface <- c(surface, temp)
            }
        }
    }
    if ("bottom" %in% layers) {
        bottom <- biooracler::download_layers(
            dataset_id = c("thetao_baseline_2000_2019_depthmean"),
            constraints = list(time = c('2001-01-01T00:00:00Z', '2010-01-01T00:00:00Z')),
            variables = c("thetao_mean")
        )
        bottom <- terra::mean(bottom)
        names(bottom) <- "bottom"

        if (do_future) {
            for (sc in seq_len(nrow(future_matrix))) {
                sel_p <- future_matrix$periods[sc]
                sel_s <- future_matrix$scenarios[sc]

                if (sel_p < 2030) {
                    cli::cli_abort("Minimum period date is 2030")
                } else if (sel_p == 2030) {
                    time_const <- c('2020-01-01T00:00:00Z', '2020-01-01T00:00:00Z')
                } else {
                    time_const <- c(paste0(sel_p - 20, '-01-01T00:00:00Z'), paste0(sel_p - 10, '-01-01T00:00:00Z'))
                }

                temp <- biooracler::download_layers(
                    dataset_id = paste0("thetao_", sel_s, "_2020_2100_depthmean"),
                    constraints = list(time = time_const),
                    variables = c("thetao_mean")
                )
                temp <- terra::mean(temp)
                names(temp) <- paste0("bottom_", sel_s, "_", sel_p)

                bottom <- c(bottom, temp)
            }
        }
    }
    if (length(layers) > 1) {
        return(c(surface, bottom))
    } else if (layers == "surface") {
        return(surface)
    } else {
        return(bottom)
    }
}
