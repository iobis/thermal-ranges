

# Load configuration file
config <- yaml::read_yaml("config.yaml")

get_thermal_ranges(
    occ_source = config$occ_source,
    species = config$species,
    aoi = config$aoi,
    output_dir = config$output_dir,
    output_filename = config$output_filename,
    output_format = config$output_format,
    future_temperature = config$future_temperature,
    future_sites = config$future_sites,
    future_mode = config$future_mode,
    future_period = config$future_period,
    future_scenarios = config$future_scenarios
)