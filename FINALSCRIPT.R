setwd("/Users/kevine/Google Drive/36-402/FINAL")
setwd("C:\\Users\\Kevin\\Google Drive\\36-402\\FINAL")

library(glmnet)
load("cancerxxx.Rdata")

#---------------------------------------------------------
#Variable screening

#Find the number of patients in each class. Store count in nk.
nk = numeric(14)
for(i in 1:14){
	nk[i] = length(which(ytrain == i))
}


#Find mean of patients in class k along gene. Store in mxkj.
#The row of mxkj is the class. The column is the gene.
mxkj = matrix(0,14,16063)
for(k in 1: 14){
	for(j in 1: 16063){
		sum = 0
		for(i in which(ytrain == k)){
			sum = sum + xtrain[i,j]
		}	
		sum = sum/nk[k]
		mxkj[k,j] = sum
	}
}

#Find the overall mean of patients along gene j. Store in mxj.
mxj = numeric(16063)
for(j in 1: 16063){
	mxj[j] = mean(xtrain[,j])
}

#Find between class variation. Store in v.
v = numeric(16063)
for(j in 1: 16063){
	sum = 0
	for(k in 1: 14){
		sum = sum + nk[k] * ( mxkj[k,j] - mxj[j])^2
	}
	v[j] = sum
}

#Screen the set of variables.
findMaxV = v
findMaxIndex = numeric(50)
for(i in 1: 50){
	findMaxIndex[i] = which.max(findMaxV)
	findMaxV[which.max(findMaxV)] = -1
}
#Plot screened genes along with variance scores.
plot(seq(1:50),v[findMaxIndex], ylab = "Variance Score", cex.lab = 1.3)

scrndIndex = findMaxIndex
scrndValues = v[scrndIndex]

scrndValues[which(diff(scrndValues) == 0)] = -1
filter = which(scrndValues > -1)

scrndIndexF = scrndIndex[filter]
scrndValuesF = scrndValues[filter]

#4 pefectly correlated variables so 50 - 4 = 46 genes to consider.

#Find largest correlation
cm = cor(xtrain[,scrndIndexF])
cm = cm - diag(46)
cm = abs(cm)
which(cm == max(cm), arr.ind = TRUE)
#largest element in position (20,14) of matrix cm


#---------------------------------------------------------
#Linear Regression of Indicators

#Function returns coefficient matrix of linear regression model.
#Two arguments: response of training set and measurement of training set
regIndicatorCoeff = function(response, predictor){
	#response set is given in the columns of yClass
	yClass = matrix(0,length(response),14)
	for(i in 1: 14){
	index = which(response == i)
	yClass[index,i] = 1
	}

	#the kth column of coeffMat contains beta^k 
	#46+1 = 46 predictors + 1 intercept
	coeffMat = matrix(0,46+1,14)

	for(k in 1: 14){
	reg.mod = lm(yClass[,k]~predictor)
	coeffMat[,k] = reg.mod$coefficients
	}
	return(coeffMat)
}

#Function gives response of indicator model
#Two arguments: Coefficient matrix of model and values to evaluate model at
predictIndicator = function(coeffMat, predictors){
	yResponse = numeric(length(predictors[,1]))
	for(i in 1: length(predictors[,1])){
		yOut = c(1,predictors[i,])%*%coeffMat
		yResponse[i] = which.max(yOut)
	}
	return(yResponse)
}



#Model fit on the entire training set
regIndicatorCoeff(ytrain, xtrain[,scrndIndexF])
p = predictIndicator(regIndicatorCoeff(ytrain, xtrain[,scrndIndexF]),xtrain[,scrndIndexF])

#Number of misclassified responses
length(which((ytrain - p) != 0)) 


#create cross validation folds
case.folds = matrix(0,4,36)
f1 = numeric(0)
f2 = numeric(0)
f3 = numeric(0)
f4 = numeric(0)
for(i in 1: 14){
	folds = rep(1:4, length.out = length(which(ytrain == i)))
	folds = sample(folds)
	classIndex = which(ytrain == i)

	f1 = c(f1,classIndex[which(folds == 1)])
	f2 = c(f2,classIndex[which(folds == 2)])
	f3 = c(f3,classIndex[which(folds == 3)])
	f4 = c(f4,classIndex[which(folds == 4)])
}

case.folds[1,] = f1
case.folds[2,] = f2
case.folds[3,] = f3
case.folds[4,] = f4


#check proportions. No response means it worked
for(i in 1:4){
	for(j in 1:14){
		if(length(which(ytrain[case.folds[i,]] == j))/length(case.folds[i,]) != nk[j]/sum(nk))
			print("We have a problem")
	}
}


scrndPredictors = xtrain[,scrndIndexF]
#The response vector for the ith fold is given in the ith column of yResp
yResp = matrix(0,length(case.folds[1,]), 4)
for(i in 1: 4){
	train.setIndex = numeric(0)
	for(j in 1: 4){
		if(j != i){
			train.setIndex = c(train.setIndex,case.folds[j,])
		}
	}
	m = regIndicatorCoeff(ytrain[train.setIndex],scrndPredictors[train.setIndex,])
	yResp[,i] = predictIndicator(m, scrndPredictors[case.folds[i,],])	
}

#calculate average misclassification rate over the four folds
misClassTot = numeric(4)
for(i in 1: 4){
	misClassTot[i] = length(which((ytrain[case.folds[i,]] - yResp[,i]) != 0))
}
mean(misClassTot)
mean(misClassTot/36)


#---------------------------------------------------------
#Lasso Method


bal.folds = numeric(144)
bal.folds[case.folds[1,]] = 1
bal.folds[case.folds[2,]] = 2
bal.folds[case.folds[3,]] = 3
bal.folds[case.folds[4,]] = 4

las.cv = cv.glmnet(xtrain,ytrain,alpha=1,family="multinomial",foldid= bal.folds)
las.mod = las.cv$glmnet.fit
lasClass.cv = cv.glmnet(xtrain,ytrain,lambda = las.cv$lambda, alpha=1,family="multinomial",type.measure="class",foldid= bal.folds)

pdf(file = "4a.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
plot(las.cv,cex.lab = 1.3)
plot(lasClass.cv,cex.lab = 1.3)
dev.off()

indexMin = which(las.cv$lambda == las.cv$lambda.min)
index1Se = which(las.cv$lambda == las.cv$lambda.1se)

lasClass.cv$cvm[indexMin]
lasClass.cv$cvm[index1Se]

m1 = coef(las.mod,s = las.cv$lambda.min)
m2 = coef(las.mod,s = las.cv$lambda.1se)


m1List = list()
m2List = list()
m1List = m1
m2List = m2
totalm1 = numeric(0)
totalm2 = numeric(0)
for(i in 1: 14){
	totalm1 = c(totalm1, as.matrix(summary(m1List[[i]]))[,1])
	totalm2 = c(totalm2, as.matrix(summary(m2List[[i]]))[,1])
}
#total number of different coefficients. Subtract one for the intercept
#to get total number of gene variables used.
length(unique(totalm1))-1
length(unique(totalm2))-1

#find shared genes between two models.
#Note unique(totalm2)-1 is used to correct for gene index. They are off by one because of intercept.
intersect(unique(totalm2)-1,scrndIndexF)

#---------------------------------------------------------
#Test set

#Model fit on the entire training set
p2 = predictIndicator(regIndicatorCoeff(ytrain, xtrain[,scrndIndexF]),xtest[,scrndIndexF])

#Number of misclassified responses for indicator model
length(which((ytest - p2) != 0)) 

yLogClass = numeric(length(xtest[,1]))
temp = predict(las.mod, xtest, s=las.cv$lambda.1se, "response")
for(i in 1:length(xtest[,1])){
	p3 = as.matrix(temp[i,,])	
	p3 = t(p3)
	yLogClass[i] = which.max(p3)
}
#Number of misclassified responses for lasso model
length(which((ytest - yLogClass) != 0)) 

#---------------------------------------------------------
#Conclusion

#get proportions on types of cancer in full data set
nk2 = numeric(14)
for(i in 1:14){
	nk2[i] = length(which(ytest == i))
}
nktot = nk + nk2
nktot/sum(nktot)

#Get total number of misclassifications for lasso method on training set
yLogClass = numeric(length(xtest[,1]))
temp = predict(las.mod, xtrain, s=las.cv$lambda.1se, "response")
for(i in 1:length(xtest[,1])){
	p3 = as.matrix(temp[i,,])	
	p3 = t(p3)
	yLogClass[i] = which.max(p3)
}
#Number of misclassified responses for lasso model
length(which((ytest - yLogClass) != 0)) 








