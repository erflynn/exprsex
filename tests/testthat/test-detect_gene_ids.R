context("Gene ids found in GPL")
library(exprsex)


test_that("find entrez ids works", {

  # check find entrez name (all caps)
  expect_equal(.detect_entrez_ids_name(toy_gpl@columns), 12)
  toy_gpl_cols2 <- toy_gpl@columns

  # find entrez name (lower/mixed case) and multiple
  toy_gpl_cols2$Column <- tolower(toy_gpl_cols2$Column)
  toy_gpl_cols2$Column[7] <- "My Entrez "
  expect_equal(.detect_entrez_ids_name(toy_gpl_cols2), c(7,12))

  # we can detect integer only columns
  expect_equal(.detect_int_cols(toy_gpl@table, MIN.OVERLAP=30), 12)
})


test_that("entrez overlap works", {
    # create a toy list of entrez ids
    toy.entrez.col <- toy_gpl@table[,12]
    gene.list <- unique(sapply(unlist(sapply(toy.entrez.col,function(x)
      strsplit(x, "///"))), stringr::str_trim))
    names(gene.list) <- NULL
    my.tmp <- tempdir()
    save(gene.list, file=sprintf("%s/%s_%s.RData", my.tmp, "human", "entrezids"))

    # check that this works for an empty column and a column where all
    # are present but sometimes multiple per line
    overlap.l <- .check_entrez_overlap(toy_gpl@table, c(7,12), "human", my.tmp)
    expect_equal( names(overlap.l), c("7", "12"));
    expect_equal(overlap.l[[1]], 0);
    expect_equal(overlap.l[[2]], 77);

 })

 test_that("find genbank ids works", {

   # by column name
   expect_equal(.detect_genbank_cols_name(toy_gpl@columns), 2)

   # by column contents
   alt_toy <- toy_gpl@table[,1:8]
   alt_toy$GB_ACC[5] <- "AAA52518"
   expect_equal(.detect_genbank_cols(alt_toy, "human"), 2)

   # returns nothing if not present
   expect_equal(length(.detect_genbank_cols(toy_gpl@table[,1:7],
                                           "human")), 0)

 })



test_that("find refseq ids works", {
  # problem - may be not refseq! could be genbank -- but still works!
  expect_equal(.detect_refseq_cols(toy_gpl@table,
                                   MIN.OVERLAP = 30), 13)

  # returns nothing if not present
  expect_equal(length(.detect_refseq_cols(toy_gpl@table[,1:7],
                                          "human")), 0)

})

test_that("find HGNC ids works", {

  alt_toy <- toy_gpl@table[,1:11]
  alt_toy[5,11] <- "PAX8 /// GAPDH"
  expect_equal(.detect_hgnc_cols(alt_toy), 11)

  # return nothing if not present
  expect_equal(length(.detect_hgnc_cols(alt_toy[,1:10])), 0)
})


test_that("find ensembl ids works", {

  # finds the column with the most ensembl ids
  alt_toy <- toy_gpl@table[,1:11]
  alt_toy[5,11] <- "ENSG1000"
  alt_toy[7,11] <- "ENSG1002"
  alt_toy[4,6] <- "ENSG1003"
  expect_equal(
    .detect_ensembl_cols(alt_toy, "human", MIN.OVERLAP=1), 11)

  # return nothing if not present
  expect_equal(
   length(.detect_ensembl_cols(toy_gpl@table, "human", MIN.OVERLAP=1)), 0)
})

test_that("find unigene ids works", {

  # finds the column with the most ensembl ids
  expect_equal(
    .detect_unigene_cols(toy_gpl@table, MIN.OVERLAP=1), 8)

  # return nothing if not present
  expect_equal(
    length(.detect_unigene_cols(toy_gpl@table[,1:5], MIN.OVERLAP=1)), 0)
})


