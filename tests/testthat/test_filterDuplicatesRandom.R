context("[NetworkDataCompanion] Testing filterDuplicatesRandom function ... ")

test_that("Testing filterDuplicatesRandom",{
  
  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()
  barcode1 = "TCGA-A1-1234-01A"
  barcode2 = "TCGA-A1-1234-11A"
  barcode3 = "TCGA-A1-1234-21A"
  barcode4 = "TCGA-A1-1234-01B"
  barcode5 = "TCGA-A1-1234-21B"
  barcode6 = "TCGA-A1-1234-21C"
  barcodes = c(barcode1, barcode2, barcode3, barcode4, barcode5, barcode6)
  set.seed(1989)
  permuted_barcodes = sample(my_friend$extractSampleAndType(barcodes))
  keep_barcodes_manual = names(permuted_barcodes[!duplicated(permuted_barcodes)])
  keep_barcodes_function = my_friend$filterDuplicatesRandom(barcodes,seed=1989)
  expect_equal(keep_barcodes_manual, keep_barcodes_function)

})
