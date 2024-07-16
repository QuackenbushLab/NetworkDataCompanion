context("[NetZooDataCompanion] Testing extract functions ... ")

test_that("extract functions correctly extract identifiers from barcode",{

  my_friend = NetZooDataCompanion::CreateNetZooDataCompanionObject()
  barcode1 = "TCGA-A1-1234-01A"
  barcode2 = "TCGA-A1-1234-11A"
  barcode3 = "TCGA-A1-1234-21A"
  barcode4 = "TCGA-A1-1234-01A"
  barcodes = c(barcode1, barcode2, barcode3, barcode4)

  samples <- my_friend$extractSampleOnly(barcodes)
  expect_equal(as.character(samples),rep("TCGA-A1-1234", 4))

  samples_and_type <- my_friend$extractSampleAndType(barcodes)
  expect_equal(as.character(samples_and_type), gsub('.{1}$', '', barcodes))

  vial <- my_friend$extractSampleAndTypeAndVial(barcodes)
  expect_equal(as.character(vial), barcodes)

  vial <- my_friend$extractVialOnly(barcodes)
  expect_equal(as.character(vial), rep("A", 4))

  sample_type <- my_friend$extractSampleType(barcodes)
  expect_equal(as.character(sample_type), c("01", "11", "21", "01"))
})
