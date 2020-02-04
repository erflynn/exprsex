#' Default genes used for sex labeling.
#'
#' A list containing two attributes: "f" and "m" which are the female and male-specific genes respectively.
#' These genes were identified using meta-analysis, and are used to sex label data unless other genes
#' are provided.
"sex_lab_genes"

#' Default fit for sex labeling.
#'
#' A list containing the "f" and "m" sex_lab genes and the default threshold.
#' This was trained with using mixed sex data from a wide range of platforms.
"default_fit"
