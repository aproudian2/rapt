library(rapt)

#### readPOS ####
test_that("POS reads correctly", {
  pos <- readPOS("example-pos.pos")
  expect_named(pos, c("x", "y", "z", "mass"))
  expect_type(pos, "list")
  expect_s3_class(pos, "data.frame")
  pos.meta <- attr(pos, "metaData")
  expect_named(pos.meta, "name")
  expect_identical(pos.meta$name, "example.pos")
})

#### readATO ####
test_that("ATO reads correctly", {
  ato <- readATO("example-ato.ato")
  expect_named(ato, c("x", "y", "z", "mass", "clusID", "pIndex", "Vdc",
                      "TOF", "dx", "dy", "Vp", "shank", "FouR", "FouI"))
  expect_type(ato, "list")
  expect_s3_class(ato, "data.frame")
  ato.meta <- attr(ato, "metaData")
  expect_named(ato.meta, "name")
  expect_identical(ato.meta$name, "example.ato")
})

#### readRRNG ####
test_that("RRNG reads correctly", {
  rng <- readRRNG("example-rrng.rrng")
  expect_named(rng, c("start", "end", "volume", "name", "formula", "color"))
  expect_type(rng, "list")
  expect_s3_class(rng, "data.frame")
  expect_type(rng$start, "double")
  expect_type(rng$end, "double")
  expect_type(rng$volume, "double")
  expect_type(rng$name, "character")
  expect_type(rng$formula, "character")
  expect_match(rng$color, "^#[[:xdigit:]]{6}$")
})
