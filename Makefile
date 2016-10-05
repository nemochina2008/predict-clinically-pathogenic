figure-BP-other-models/index.html: figure-BP-other-models.R BP.other.models.RData trainData.RData
	R --no-save < $<
BP.other.models.RData: BP.other.models.R BP.other.RData
	R --no-save < $<
BP.other.RData: BP.other.R
	R --no-save < $<
testData.RData: testData.R 
	R --no-save < $<
testData.folds.RData: testData.folds.R testData.RData
	R --no-save < $<
figure-some-test-error.png: figure-some-test-error.R some.models.RData
	R --no-save < $<
figure-two-features.png: figure-two-features.R
	R --no-save < $<
figure-test-error.png: figure-test-error.R models.RData
	R --no-save < $<
models.RData: models.R folds.RData trainData.RData  
	R --no-save < $<
some.models.RData: some.models.R folds.RData trainData.RData  
	R --no-save < $<
folds.RData: folds.R trainData.RData
	R --no-save < $<
trainData.RData: trainData.R
	R --no-save < $<
