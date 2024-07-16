context("[NetZooDataCompanion] Testing filterSampleType function ... ")

test_that("filterSampleType correctly selects samples matching input types",{

  my_friend = NetZooDataCompanion::CreateNetZooDataCompanionObject()
  barcode1 = "TCGA-A1-1234-01A"
  barcode2 = "TCGA-A1-1234-11A"
  barcode3 = "TCGA-A1-1234-21A"
  barcode4 = "TCGA-A1-1234-01A"
  barcodes = c(barcode1, barcode2, barcode3, barcode4)

  idx_01 = my_friend$filterSampleType(barcodes,types_of_samples = c("01"))
  expect_equal(idx_01$index,c(1,4))

  idx_11 = my_friend$filterSampleType(barcodes,types_of_samples = c("11"))
  expect_equal(idx_11$index,c(2))

  idx_01_11 = my_friend$filterSampleType(barcodes,types_of_samples = c("01","11"))
  expect_equal(idx_01_11$sample_type,c("01","11","01"))
  expect_equal(idx_01_11$index,c(1,2,4))

  barcode0 = "TCGA-A1-1234-99A" # push something with an invalid sample type on top of vector
  barcodes = c(barcode0,barcodes)
  idx_01_with_99 = my_friend$filterSampleType(barcodes,types_of_samples = c("01"))
  expect_equal(idx_01_with_99$index,c(2,5))

  })

test_that("filterNormalSamples, filterTumorSamples, filterControlSamples correctly select samples",
{
  my_friend = NetZooDataCompanion::CreateNetZooDataCompanionObject()
  barcode1 = "TCGA-A1-1234-11A"
  barcode2 = "TCGA-A1-1234-21A"
  barcode3 = "TCGA-A1-1234-11B"
  barcode4 = "TCGA-A1-1234-01A"
  barcode5 = "TCGA-A1-1234-06B"
  barcode6 = "TCGA-A1-1234-**B" # something that's not a valid code, make sure it is not spuriously picked up
  barcodes = c(barcode1, barcode2, barcode3,barcode4,barcode5,barcode6)

  idx_normal = my_friend$filterNormalSamples(barcodes)
  expect_equal(idx_normal, c(1,3))
  idx_tumor = my_friend$filterTumorSamples(barcodes)
  expect_equal(idx_tumor, c(4,5))
  idx_control = my_friend$filterControlSamples(barcodes)
  expect_equal(idx_control,c(2))

  })
