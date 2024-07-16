context("[NetZooDataCompanion] Testing mapUUIDtoTCGA function ... ")

test_that("Testing mapUUIDtoTCGA",{

  my_friend = NetZooDataCompanion::CreateNetZooDataCompanionObject()
  # Example from
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAutils/inst/doc/TCGAutils.html#uuid-history-lookup
  # Visited 08/24/23
  expect_equal(
    as.character(my_friend$mapUUIDtoTCGA("b4bce3ff-7fdc-4849-880b-56f2b348ceac")['submitter_id']),
    "TCGA-B0-5094-11A"
    )
})
