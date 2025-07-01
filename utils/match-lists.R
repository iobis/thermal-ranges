# Match the results of a species thermal ranges file with lists of species
match_lists <- function(thermal_ranges, lists, output_folder) {

    fs::dir_create(output_folder)
    tr <- arrow::open_dataset(thermal_ranges)

    if (length(lists) == 1) {
        if (tools::file_ext(lists) == "" & dir.exists(lists)) {
            lists <- list.files(lists, full.names = T)
            lists <- lists[grepl("csv$|tsv$|txt$|parquet$", lists)]
        }
    }

    cli::cli_alert_info("Processing {length(lists)} list{?s}")

    for (l in seq_along(lists)) {
        sel_f <- lists[l]
        cli::cli_progress_step("Processing {basename(sel_f)}")
        list_content <- switch(
            tools::file_ext(sel_f),
            csv = data.table::fread(sel_f),
            tsv = data.table::fread(sel_f),
            txt = data.table::fread(sel_f),
            parquet = arrow::read_parquet(sel_f)
        )
        colnames(list_content)[grepl("aphiaid", tolower(colnames(list_content)))] <- "aphiaID"

        tr_l <- tr |>
            dplyr::filter(aphiaID %in% list_content$aphiaID) |>
            dplyr::collect()

        list_content <- dplyr::left_join(
            list_content, tr_l |> dplyr::select(-species),
            by = "aphiaID"
        )

        data.table::fwrite(
            list_content, file.path(output_folder, basename(sel_f))
        )
    }

    cli::cli_progress_done()
    cli::cli_alert_success("All files saved at folder {.path {output_folder}}")

    return(invisible(NULL))
}
