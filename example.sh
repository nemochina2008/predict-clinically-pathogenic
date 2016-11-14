## Train a model based on the labels in data-2016-11-09/train.txt and
## save it to data-2016-11-09/train.txt.RData
Rscript train_model.R data-2016-11-09/train.txt

## Read the model saved in data-2016-11-09/train.txt.RData and use it
## to compute predicted scores for data-2016-11-09/test1.txt -- the
## scores and predictions with un-modified input features are in
## data-2016-11-09/test1.txt_predictions.txt, and the modified input
## features with univariate model predictions for each data point (no
## missing values) are in data-2016-11-09/test1.txt_predictions.RData
Rscript predict_scores.R data-2016-11-09/train.txt.RData data-2016-11-09/test1.txt

## Read the predicted scores in
## data-2016-11-09/test1.txt_predictions.RData and plot ROC curves.
Rscript plot_ROC.R data-2016-11-09/test1.txt_predictions.RData
