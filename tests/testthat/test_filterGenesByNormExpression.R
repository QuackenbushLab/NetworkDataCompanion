context("[NetworkDataCompanion] Testing filterGenesByNormExpression function ... ")

test_that("Testing filterGenesByNormExpression",{

  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()
  expr <- rbind(c(1, 2, 3, 4), c(0, 0, 1, 1), c(0, 1, 2, 1), c(0, 0, 10, 10))

  # Test illegal inputs
  expect_error(my_friend$filterGenesByNormExpression(expr, 1, -1))
  expect_error(my_friend$filterGenesByNormExpression(expr, 1, 1.1))
  expect_no_error(my_friend$filterGenesByNormExpression(expr, 1, 0.5))
  expect_error(my_friend$filterGenesByNormExpression(expr, -1, 0.5))
  expect_error(my_friend$filterGenesByNormExpression(expr, "error", 0.5))

  # Test functionality
  expect_equal(my_friend$filterGenesByNormExpression(expr, 2, 0.5), c(1, 4))
  expect_equal(my_friend$filterGenesByNormExpression(expr, 2, 0.75), c(1))
  expect_equal(my_friend$filterGenesByNormExpression(expr, 10, 0.4), c(4))
})
