# Create a function to check if is installed
is_package_installed <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
}

# Check which ones are not installed and install if needed:
req_packages <- read.table("requirements.txt")[, 1]
for (i in seq_along(req_packages)) {
    if (!is_package_installed(req_packages[i])) {
        install.packages(req_packages[i])
    }
}
