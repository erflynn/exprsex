# Code for examining an example dataset


### with MetaIntegrator
require('MetaIntegrator')
gse <- getGEOData("GSE55668")
exp_rank <- expDataToRanks(gse$originalData$GSE55668)
sex_lab <- gse$originalData$GSE55668$pheno$`Sex:ch1`
map_list <- list("Female"=0, "Male"=1)
key_mapping <- gse$originalData$GSE55668$keys
sex_lab_reform <- unlist(sapply(sex_lab, function(x) map_list[[x]]))
names(sex_lab_reform) <- colnames(exp_rank)

fit_ex <- trainSexLab(exp_rank, sex_lab_reform)
sex_lab_f <- predSexLab(fit_ex, exp_rank)

print(sex_lab)


### with GEOquery
require('GEOquery')
gse2 <- getGEO("GSE55668")
expData <- exprs(gse2$GSE55668_series_matrix.txt.gz)
fData <- fData(gse2$GSE55668_series_matrix.txt.gz)
mapToGenes <- fData[,c("ID", "GENE_SYMBOL")]
gene_symbol_list <- mapToGenes$GENE_SYMBOL
# TODO - does this handle multiple symbols

names(gene_symbol_list) <- mapToGenes$ID
sex_lab <- pData(gse2$GSE55668_series_matrix.txt.gz)$`Sex:ch1`
map_list <- list("Female"=0, "Male"=1)

sex_lab_reform <- unlist(sapply(sex_lab, function(x) map_list[[x]]))
names(sex_lab_reform) <- colnames(expData)
exp_rank2 <- expDataToRanks(expData, gene_symbol_list)
fit_ex2 <- trainSexLab(exp_rank2, sex_lab_reform)
sex_lab_f2 <- predSexLab(fit_ex2, exp_rank2)

