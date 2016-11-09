## Train a model based on the labels in data-2016-11-09/train.txt and
## save it to data-2016-11-09/train.txt.RData
Rscript train_model.R data-2016-11-09/train.txt

## Read the model saved in data-2016-11-09/train.txt.RData and use it
## to compute predicted scores for data-2016-11-09/test1.txt -- the
## predicted scores are in data-2016-11-09/test1.txt_predictions.txt 
Rscript predict_scores.R data-2016-11-09/train.txt.RData data-2016-11-09/test1.txt

## Read the predicted scores in data-2016-11-09/test1.txt_predictions.txt and plot ROC curves to 
Rscript plot_ROC.R data-2016-11-09/test1.txt_predictions.txt data-2016-11-09/test1.txt_predictions.txt.png
