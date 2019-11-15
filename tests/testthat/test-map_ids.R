context("ID mapping")
library(exprsex)

test_that("refseq mapping works", {
  my.tmp <- tempdir()
  save(ref_dat_sm,
       file=sprintf("%s/%s_%s.RData", my.tmp, "human", "gene_map"))
  mapped <- .map_from_refseq(toy_gpl@table, 13, "human", ref_dir = my.tmp)
  expect_equal(nrow(mapped),62)
})


test_that("hgnc mapping works", {
  my.tmp <- tempdir()
  save(ref_dat_sm,
       file=sprintf("%s/%s_%s.RData", my.tmp, "human", "gene_map"))
  alt_toy <- toy_gpl@table[,1:11]
  alt_toy[5,11] <- "PAX8 /// GAPDH"
  mapped <- .map_from_hgnc(alt_toy[1:50,], 11, "human", ref_dir=my.tmp)
  expect_equal(nrow(mapped),50)
})


test_that("ensembl mapping works", {
  my.tmp <- tempdir()
  save(ref_dat_sm,
       file=sprintf("%s/%s_%s.RData", my.tmp, "human", "gene_map"))
  alt_toy <- toy_gpl@table[,1:8]
  set.seed(1115)
  alt_toy[,6] <- ref_dat_sm$ensembl_gene_id[sample(1:nrow(ref_dat_sm), nrow(alt_toy))]
  mapped <- .map_from_ensembl(alt_toy, 6, "human", ref_dir=my.tmp)
  expect_equal(nrow(mapped),61)
})


test_that("unigene mapping works", {
  my.tmp <- tempdir()
  save(unigene_sm,
       file=sprintf("%s/%s_%s.RData", my.tmp, "human", "unigene"))
  mapped <- .map_from_unigene(toy_gpl@table, 8, "human", ref_dir = my.tmp)
  expect_equal(nrow(mapped), 5)
})


test_that("genbank mapping works", {
  my.tmp <- tempdir()
  save(genbank_sm,
       file=sprintf("%s/%s_%s.RData", my.tmp, "human", "genbank"))
  alt_toy <- toy_gpl@table[,1:8]
  alt_toy$GB_ACC[5] <- "AAA52518"
  mapped <- .map_from_genbank(alt_toy, 2, "human", ref_dir = my.tmp)
  expect_equal(nrow(mapped), 52)
})
