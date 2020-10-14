library(rapt)

#### readPOS ####
test_that("POS reads correctly", {
  pos <- readPOS("example-pos.pos")

  expect_named(pos, c("x", "y", "z", "mass"))
  expect_type(pos, "list")
  expect_s3_class(pos, "data.frame")
  pos.meta <- attr(pos, "metaData")
  expect_named(pos.meta, "name")
  expect_identical(pos.meta$name, "example-pos")
})
pos <- readPOS("example-pos.pos")

#### readATO ####
test_that("ATO reads correctly", {
  ato <- readATO("example-ato.ato")

  expect_named(ato, c("x", "y", "z", "mass", "clusID", "pIndex", "Vdc",
                      "TOF", "dx", "dy", "Vp", "shank", "FouR", "FouI"))
  expect_type(ato, "list")
  expect_s3_class(ato, "data.frame")
  ato.meta <- attr(ato, "metaData")
  expect_named(ato.meta, "name")
  expect_identical(ato.meta$name, "example-ato")
})
ato <- readATO("example-ato.ato")

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
  rng.meta <- attr(rng, "metaData")
  expect_named(rng.meta, "name")
  expect_identical(rng.meta$name, "example-rrng")
})

#### read.rcp ####
# scaleUp = TRUE should be tested with scaleRCP function
test_that("read.rcp reads correctly", {
  rcp <- read.rcp("FinalConfig1", "system1")

  expect_type(rcp, "list")
  expect_s3_class(rcp, "ppx")
  expect_s3_class(rcp, "pp3")
  expect_equal(signif(min(nndist(rcp)), digits = 2), 0.05)
  expect_equal(round(rcp$domain$xrange), c(0,1))
  expect_equal(round(rcp$domain$yrange), c(0,1))
  expect_equal(round(rcp$domain$zrange), c(0,1))
})

#### createSpat ####
test_that("createSpat converts POS and ATO data.frames to pp3", {
  pp.pos <- createSpat(pos)

  expect_type(pp.pos, "list")
  expect_s3_class(pp.pos, "ppx")
  expect_s3_class(pp.pos, "pp3")
  pos.meta <- attr(pp.pos, "metaData")
  expect_named(pos.meta, "name")
  expect_identical(pos.meta$name, "example-pos")
  expect_equal(round(pp.pos$domain$xrange), c(-20,20))
  expect_equal(round(pp.pos$domain$yrange), c(-20,20))
  expect_equal(round(pp.pos$domain$zrange), c(0,51))

  pp.clip <- createSpat(pos, win = box3())
  expect_equal(domain(pp.clip),box3())
  expect_equal(pp.clip$data, pp.pos[inside.boxx(pp.pos, w = box3())]$data)

  pp.ato <- createSpat(ato)

  expect_type(pp.ato, "list")
  expect_s3_class(pp.ato, "ppx")
  expect_s3_class(pp.ato, "pp3")
  ato.meta <- attr(pp.ato, "metaData")
  expect_named(ato.meta, "name")
  expect_identical(ato.meta$name, "example-ato")
  expect_equal(round(pp.ato$domain$xrange), c(-202,204))
  expect_equal(round(pp.ato$domain$yrange), c(-204,205))
  expect_equal(round(pp.ato$domain$zrange), c(0,507))

})

#### createSpec ####
test_that("createSpec converts POS and ATO data.frames to MassSpectrum", {
  pos.ms <- createSpec(pos)
  expect_type(pos.ms, "S4")
  expect_s4_class(pos.ms, "MassSpectrum")
  pos.meta <- attr(pos.ms, "metaData")
  expect_named(pos.meta, "name")
  expect_identical(pos.meta$name, "example-pos")

  pos.ms <- createSpec(pos, clip = c(10,50))
  expect_type(pos.ms, "S4")
  expect_s4_class(pos.ms, "MassSpectrum")

  ato.ms <- createSpec(ato)
  expect_type(ato.ms, "S4")
  expect_s4_class(ato.ms, "MassSpectrum")
  ato.meta <- attr(ato.ms, "metaData")
  expect_named(ato.meta, "name")
  expect_identical(ato.meta$name, "example-ato")

  ato.ms <- createSpec(ato, clip = c(10,50))
  expect_type(ato.ms, "S4")
  expect_s4_class(ato.ms, "MassSpectrum")
  ato.meta <- attr(ato.ms, "metaData")
  expect_named(ato.meta, "name")
  expect_identical(ato.meta$name, "example-ato")
})
