#' Helper function to generate a predicted probability df for
#' a given data frame and fit that contains the scores and predictions
#'
#' @param fit a sex labeling fit object
#' @param df the expression matrix, rows are genes and samples are columns
.predDf <- function(fit, df){
  require(dplyr)
  pred_df <- predSexLab(fit, df, scores=TRUE)
  pred_df <- pred_df %>%
    dplyr::rename(pred_sex=sex) %>%
    dplyr::mutate(diff_score=(score_m-score_f))
}

#' Fit a logistic curve to the sex labeling fit so that we can get
#' accurate proba scores for the output.
#' Note that it is possible (and encouraged!) to use a representative
#' but different dataset to fit the curve as used to generate the SL object.
#'
#' @param fit a sex labeling fit object
#' @param df the expression matrix, rows are genes and samples are columns
#' @param labels sex labels in numeric (0/1) form for the expression mat
trainProba <- function(fit=NULL, df, labels){
  # // TODO input checks
  pred_df <- .predDf(fit, df)
  pred_df$true_sex <- labels
  lm.fit <- glm(true_sex ~ diff_score, data=pred_df, family="binomial")
  return(list("lm.fit"=lm.fit, "fit"="fit"))
}

#' Return the predicted probabilities and sex labels for an expression matrix.
#'
#' @param fit a sex labeling fit object with a logisitic fit
#' @param df the expression matrix, rows are genes and samples are columns
#' @return a data frame with a column for samples, sex labels, and proba scores
predProba <- function(fit=NULL, df){
  # // TODO input checks
  pred_df <- .predDf(fit$fit, df)
  preds <- predict(fit$lm.fit, pred_df, type="response")
  pred_class <- ifelse(preds >= 0.5, 1, 0)
  return(data.frame(cbind("sample"=colnames(df), "sex"=pred_class, "prob"=preds)))
}
