# Get thermal ranges for all species within the speciesgrids product
library(dplyr)
source("get-ranges-fun.R")
source("utils/match-lists.R")
source("utils/get-temperature.R")

# Download the species grids if not existent
grids_path <- "/Volumes/OBIS2/speciesgrids/h3_7"

get_thermal_ranges(
    occ_source = "speciesgrids",
    occ_path = grids_path,
    species = "all",
    min_records = 10,
    output_format = "parquet"
)

# Match lists
base_url <- "https://github.com/iobis/whomp-species-lists/raw/refs/heads/main/output"
whs_lists <- c(
    "biosphere_marine.tsv",
    "geoparks_marine.tsv",
    "wh_cult_marine.tsv",
    "wh_nat_marine.tsv"
)
whs_lists <- file.path(base_url, whs_lists)

match_lists("results/thermal_ranges_20250701.parquet", whs_lists, "whs_sites")

# Extract temperature at each site and add to file
# Rename original files to match lists
sites <- list.files("whs_shapes")
sites

sites_new <- sites |>
    (\(x) {gsub("WH_Cultural_v2_marine", "wh_cult_marine", x)})() |>
    (\(x) {gsub("WH_Natural_v2_marine", "wh_nat_marine", x)})() |>
    (\(x) {gsub("BRmay25_v2_marine", "biosphere_marine", x)})() |>
    (\(x) {gsub("UGGP_v2_marine", "geoparks_marine", x)})()

file.rename(
    file.path("whs_shapes", sites),
    file.path("whs_shapes", sites_new)
)

get_temperature_sites(
    "whs_shapes", "whs_sites",
    periods = c(2050, 2100),
    scenarios = c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585")
)

# Match with lists to create sensitivity score
lists <- c(
    "biosphere_marine",
    "geoparks_marine",
    "wh_cult_marine",
    "wh_nat_marine"
)

for (i in seq_along(lists)) {
    sel_site <- lists[i]

    species <- data.table::fread(file.path("whs_sites", paste0(sel_site, ".tsv")))

    species <- species |>
        select(1:10, s_tmax = surface_q95,
               b_tmax = bottom_q95,
               s_n = surface_valid_count,
               b_n = bottom_valid_count,
               s_conf = surface_confidence,
               b_conf = bottom_confidence) |>
        mutate(s_tmax = round(s_tmax, 2), b_tmax = round(b_tmax, 2))

    site_temp <- read.csv(file.path("whs_sites", paste0(sel_site, "_exttemp.csv")))
    site_temp <- site_temp |> rename(area_id = ID) |>
        select(c("area_id", starts_with("surface"), starts_with("bottom")))

    species_s <- left_join(species, site_temp, by = join_by(area_id))

    species_s <- species_s |>
        mutate(across(starts_with("surface"), ~ round(s_tmax - .x, 2))) |>
        mutate(across(starts_with("bottom"), ~ round(b_tmax - .x, 2))) |>
        rename(surface_current = surface, bottom_current = bottom) |>
        rename(surface_tmax = s_tmax, bottom_tmax = b_tmax, surface_n = s_n, bottom_n = b_n,
               surface_confidence = s_conf, bottom_confidence = b_conf)
    
    data.table::fwrite(species_s, file.path("whs_sites", paste0(sel_site, "_sensitivity.tsv")))
}

# On final files, a negative value means that the species is out of the thermal limit
# by that amount. Positive values means that the species is below the thermal limit
# by that amount.