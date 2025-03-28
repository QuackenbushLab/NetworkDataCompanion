context("[NetworkDataCompanion] Testing filterDuplicatesSeqDepthOther function ... ")

test_that("Testing filterDuplicatesSeqDepthOther",{

  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()
  barcode1 = "TCGA-A1-1234-01A"
  barcode2 = "TCGA-A1-1234-11A"
  barcode3 = "TCGA-A1-1234-21A"
  barcode4 = "TCGA-A1-1234-01B"
  barcode5 = "TCGA-A1-1234-21B"
  barcode6 = "TCGA-A1-1234-21C"
  barcodes = c(barcode1, barcode2, barcode3, barcode4, barcode5, barcode6)

  expr <- data.frame(rbind(c(1, 2, 3, 4, 5, 6), c(rep(1, 6))))
  colnames(expr) <- barcodes
  expect_equal(my_friend$filterDuplicatesSeqDepthOther(expr,barcodes), c(6))
})
