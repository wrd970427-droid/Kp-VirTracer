test_that("stage-1 config example exists", {
  cfg <- system.file("extdata", "config.yaml.example", package = "KpVirTracer")
  expect_true(file.exists(cfg))
})
