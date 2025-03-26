context("[NetworkDataCompanion] Testing getTissueType function ... ")

numcodes <- c(
  "01",
  "02",
  "03",
  "04",
  "05",
  "06",
  "07",
  "08",
  "09",
  "10",
  "11",
  "12",
  "13",
  "14",
  "15",
  "16",
  "20",
  "40",
  "50",
  "60",
  "61",
  "99"
)

descriptions <- c(
  "Primary Solid Tumor",
  "Recurrent Solid Tumor",
  "Primary Blood Derived Cancer - Peripheral Blood",
  "Recurrent Blood Derived Cancer - Bone Marrow",
  "Additional - New Primary",
  "Metastatic",
  "Additional Metastatic",
  "Human Tumor Original Cells",
  "Primary Blood Derived Cancer - Bone Marrow",
  "Blood Derived Normal",
  "Solid Tissue Normal",
  "Buccal Cell Normal",
  "EBV Immortalized Normal",
  "Bone Marrow Normal",
  "sample type 15",
  "sample type 16",
  "Control Analyte",
  "Recurrent Blood Derived Cancer - Peripheral Blood",
  "Cell Lines",
  "Primary Xenograft Tissue",
  "Cell Line Derived Xenograft Tissue",
  "sample type 99"
)
my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()

test_that("Testing getTissueType for a single TCGA barcode",{

  barcode = "TCGA-A1-1234-"
  for(numcode in numcodes){
   expect_equal(as.character(my_friend$getTissueType(paste0(barcode, numcode))['description']), descriptions[numcodes == numcode])
  }
  expect_true(is.null(my_friend$getTissueType(paste0(barcode, "53"))))
})

test_that("Testing getTissueType vectorizes properly",{
  
  barcode1 = "TCGA-A1-1234-01A"
  barcode2 = "TCGA-A1-1234-10A"
  barcode3 = "TCGA-A1-1234-10B"
  barcode4 = "TCGA-A1-1234-02A"
  barcodes = c(barcode1,
               barcode2,
               barcode3,
               barcode4)
  
  expect_equal(my_friend$getTissueType(c(barcode1,
                                         barcode2,
                                         barcode3,
                                         barcode4))$description,
                                       c(descriptions[1],
                                         descriptions[10],
                                         descriptions[10],
                                         descriptions[2]))
})

test_that("Testing getTissueType is not bamboozled by extra suffix items",{
  
  barcode = "TCGA-A1-1234-01A-1989"
  expect_equal(my_friend$getTissueType(barcode)$description,descriptions[1])
  
})
