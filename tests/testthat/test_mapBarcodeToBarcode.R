context("[NetSciDataCompanion] Testing mapBarcodeToBarcode function ... ")

test_that("mapBarcodeToBarcode function correctly extract information for a pair of barcodes",{
  
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  
  # Basic test
  bc1 <- c("a", "b", "c")
  bc2 <- c("b", "e")
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter1,c(FALSE, TRUE, FALSE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs1,c(NA, 1, NA))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter2,c(TRUE, FALSE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs2,c(2, NA))
  
  # Test empty list
  bc1 <- c("a")
  bc2 <- c("")
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter1,c(FALSE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs1,c(NA))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter2,c(FALSE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs2,c(NA))
  
  # Test subset
  bc1 <- c("a", "b")
  bc2 <- c("aa", "a", "a1", "b")
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter1,c(TRUE, TRUE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs1,c(2, 4))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter2,c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs2,c(NA, 1, NA, 2))
  
  # Test superset
  bc1 <- c("aa", "a", "a1", "b")
  bc2 <- c("a", "b")
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter1,c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs1,c(NA, 1, NA, 2))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter2,c(TRUE, TRUE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs2,c(2, 4))
  
  # Test equal lists
  bc1 <- c("a", "b")
  bc2 <- c("a", "b")
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter1,c(TRUE, TRUE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs1,c(1, 2))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter2,c(TRUE, TRUE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs2,c(1, 2))
  
  # Test permuted
  bc1 <- c("a", "b")
  bc2 <- c("b", "a")
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter1,c(TRUE, TRUE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs1,c(2, 1))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$is_inter2,c(TRUE, TRUE))
  expect_equal(my_friend$mapBarcodeToBarcode(bc1, bc2)$idcs2,c(2, 1))
})