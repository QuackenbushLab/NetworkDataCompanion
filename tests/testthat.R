library(testthat)

# Tests are files in the testthat directory
# They will be run alphabetically

test_out = test_check(package = "NetworkDataCompanion",
           reporter = "summary")
