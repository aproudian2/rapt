library(rapt)

#### readResult ####
test_that("TAPSim results read correctly", {
  res.3225 <- readResult('results_3225.00000001')
  expect_named(res.3225, c("index","id","number","voltage",
                           "startX","startY","startZ",
                           "stopX","stopY","stopZ","tof","probability",
                           "potentialBefore",
                           "fieldBeforeX","fieldBeforeY","fieldBeforeZ",
                           "potentialAfter",
                           "fieldAfterX","fieldAfterY","fieldAfterZ",
                           "normalX","normalY","normalZ",
                           "apexX","apexY","apexZ"))
  expect_type(res.3225, "list")
  expect_s3_class(res.3225, "data.frame")

  res.557 <- readResult('results_557.00000001')
  expect_named(res.557, c("index","id","number","voltage",
                          "startX","startY","startZ",
                          "stopX","stopY","stopZ","tof","probability",
                          "normalX","normalY","normalZ",
                          "apexX","apexY","apexZ"))
  expect_type(res.557, "list")
  expect_s3_class(res.557, "data.frame")
})
