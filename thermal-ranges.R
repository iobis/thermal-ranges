# Get thermal ranges for all species within the speciesgrids product
source("get-ranges-fun.R")

# Download the species grids if not existent
grids_path <- "data/speciesgrids/h3_7"
if (dir.exists(grids_path)) {
    # TODO: add download code
}

get_thermal_ranges(
    occ_source = "speciesgrids",
    occ_path = grids_path,
    species = "all",
    min_records = 10
)