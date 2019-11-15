context("Gene id parsing")
library(exprsex)

test_that("parse id from a multi-col works", {

  expect_equal(sum(.parse_multi_col(toy_gpl@table, 8, "Hs\\.[0-9]+")$gene_col!=""), 40)
  expect_equal(nrow(.parse_multi_col(toy_gpl@table, 12, "[0-9]+")), 78)
})

test_that("find loc of id within column works", {
  expect_equal(nrow(.find_col_loc(toy_gpl@table, 13, "N[R|M][_][\\d]+[_.-]*[\\w\\d]*")), 254)
})
