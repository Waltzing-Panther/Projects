#Initialize
rm(list=ls(all=TRUE)) 
setwd("C:/Users/Kevin/Google Drive/36-401/TEST REPORT 2")
mydata <- read.csv(file="forestfires.csv",head=TRUE,sep=",")
attach(mydata)
library(MASS)


newdata = mydata
for(i in  0:516){
	if(mydata$area[517 - i] >= 300){
		newdata = newdata[-(517 - i),]
	}
	
}
attach(newdata)


newdata = mydata
for(i in  0:516){
	if(mydata$area[517 - i] >= 300){
		newdata = newdata[-(517 - i),]
	}
	
}
attach(newdata)



#UDA
Month = NULL
Month[length(month)] = NA
Month[month == "jan"] = 1
Month[month == "feb"] = 2
Month[month == "mar"] = 3
Month[month == "apr"] = 4
Month[month == "may"] = 5
Month[month == "jun"] = 6
Month[month == "jul"] = 7
Month[month == "aug"] = 8
Month[month == "sep"] = 9
Month[month == "oct"] = 10
Month[month == "nov"] = 11
Month[month == "dec"] = 12
month = Month

Day = NULL
Day[length(day)] = NA
Day[day == "mon"] = 1
Day[day == "tue"] = 2
Day[day == "wed"] = 3
Day[day == "thu"] = 4
Day[day == "fri"] = 5
Day[day == "sat"] = 6
Day[day == "sun"] = 7
day = Day

par(mfrow=c(4,4))
vars = cbind(names(mydata))

for(i in 1:8){
	hist(get(vars[i]),breaks = 10,xlab = vars[i],main = paste("Histogram of", vars[i], sep = " "))
}
i =9
hist(get(vars[i]),breaks = 10,xlab = "temp (C)",main = paste("Histogram of", vars[i], sep = " "))
i = 10
hist(get(vars[i]),breaks = 10,xlab = "RH (%)",main = paste("Histogram of", vars[i], sep = " "))
i = 11
hist(get(vars[i]),breaks = 10,xlab = "wind (km/h)",main = paste("Histogram of", vars[i], sep = " "))
i = 12
hist(get(vars[i]),breaks = 10,xlab = "rain (mm/m^2)",main = paste("Histogram of", vars[i], sep = " "))
i = 13
hist(get(vars[i]),breaks = 10,xlab = "area (ha)",main = paste("Histogram of", vars[i], sep = " "))


par(mfrow=c(2,2))
boxplot(area~X,xlab = "X", ylab = "Area", ylim = range(0:40), main = "Area vs. X")
boxplot(area~Y,xlab = "Y", ylab = "Area",ylim = range(0:40), main = "Area vs. Y")
boxplot(area~month,xlab = "Month", ylab = "Area",ylim = range(0:40), main = "Area vs. Month")
boxplot(area~day,xlab = "Day", ylab = "Area",ylim = range(0:40), main = "Area vs. Day")




#Table Generator
M = matrix(nrow = length(names(mydata)),ncol = 5)
for(i in 1:length(names(mydata))){
	M[i,1] = mean(get(vars[i]))
	M[i,2] = median(get(vars[i]))
	M[i,3] = max(get(vars[i])) - min(get(vars[i]))
	M[i,4] = IQR(get(vars[i]))
	M[i,5] = var(get(vars[i]))
}


reg.line = lm(area~ X + Y + summerFallMonth + day + FFMC + DMC + DC + ISI + wind + Rain + RH + temp + FFMC*DMC + FFMC*DC + FFMC*ISI + DMC*DC+DMC*ISI +DC*ISI)

reg.line = lm(area~  summerFallMonth  + FFMC + DMC + DC + ISI + wind + rain + RH)

summerMonth = ifelse(Month >= 6 & Month <= 9,1,0)
summerFallMonth = ifelse(Month >= 6 & Month <= 12,1,0)
Rain = ifelse(rain > 0,1,0)
Ynew = ifelse(Y > 6,1,0)

par(mfrow=c(4,2))
reg.line = lm(1/sqrt(area+1)~ summerFallMonth + FFMC + DMC + DC + ISI + wind + RH)
F = log(1/FFMC^(2))
reg.line = lm(area~ summerFallMonth   + F + DMC + DC + ISI + wind + rain)
reg.line = lm(area~ DMC  + RH)

reg.line = lm(area~ F)
reg.line = lm(1/sqrt(area+1)~ FFMC)

plot(F,reg.line$res)
plot(FFMC,reg.line$res)
plot(DMC,reg.line$res)
plot(DC,reg.line$res)
plot(ISI,reg.line$res)
plot(wind,reg.line$res)
plot(rain,reg.line$res)
plot(RH,reg.line$res)


qqnorm(reg.line$res,pch=16,main="Linear Model Residuals: QQ Plot")
qqline(reg.line$res,lwd = 2,col = 2)


reg.line = lm(asin(sqrt(area))~ISI)
plot(reg.line$fit, reg.line$res , xlab = "Fitted Value" , ylab = "Residual", pch = 16)

reg.line = lm(1/sqrt(area+1)~ISI)
plot(reg.line$fit, reg.line$res , xlab = "Fitted Value" , ylab = "Residual", pch = 16)


reg.line = lm(area~ X + Y + month + day + FFMC + DMC + DC + ISI + wind + rain)
plot(reg.line$fit, reg.line$res , xlab = "Fitted Value" , ylab = "Residual", pch = 16)
qqnorm(reg.line$res,pch=16,main="Linear Model Residuals: QQ Plot")
qqline(reg.line$res,lwd = 2,col = 2)








par(mfrow=c(4,2))

reg.line = lm(1/sqrt(area+1)~ FFMC  + DC + ISI + wind + summerFallMonth + Rain)
reduced = lm(1/sqrt(area+1)~ wind + summerFallMonth + temp)
plot(FFMC,reg.line$res)
plot(DMC,reg.line$res)
plot(DC,reg.line$res)
plot(ISI,reg.line$res)
plot(wind,reg.line$res)
plot(Rain,reg.line$res)


Temp = 1/temp^2
summerMonth2 = 1/summerMonth^2

reg.line = lm(1/sqrt(area+1)~ wind + Rain + summerMonth*temp)
reg.line = lm(1/sqrt(area+1)~ wind)

reg.line = lm(area~ exp(ISI))
Area = 1/sqrt(area+1)
reg.line = lm(Area ~ ISI)
reg.line = lm(area~ X + Y + FFMC + DMC + DC + ISI + wind + Rain + summerMonth*temp)
reg.line = lm(area~ I(1/(X^2)) + I(1/(Y^2)) + I(1/(FFMC^2)) + I(1/(DMC^2)) + I(1/(DC^2))  + I(1/(wind^2)) + Rain + summerMonth*Temp)
reg.line = lm(Area~Wind)
reg.line = lm(area~I(1/(wind^2)))
reg.line = lm(1/sqrt(Area)~Wind)
reg.line = lm(log(area)~ FFMC + DMC + DC + ISI + wind + Rain + summerMonth*temp)
reg.line = lm(log(area)~ FFMC + DMC + DC + ISI + RH + summerMonth)
reduced = lm(log(area)~ FFMC + DMC + ISI + wind)
reduced = lm(log(area)~ DC)
reg.line = lm(log(area)~ X + Y + FFMC + DMC + DC + ISI + wind + Rain + summerMonth*temp + FFMC*DMC + FFMC*DC + FFMC*ISI + DMC*DC+DMC*ISI +DC*ISI)
reg.line = lm(log(area)~ Ynew + X + FFMC + DMC + DC + ISI + summerMonth + X*Ynew)

plot(reg.line$fit, reg.line$res , xlab = "Fitted Value" , ylab = "Residual", pch = 16)
qqnorm(reg.line$res,pch=16,main="Linear Model Residuals: QQ Plot")
qqline(reg.line$res,lwd = 2,col = 2)
summary(reg.line)
plot(ISI,area)
plot(FFMC,area)





reg.line = lm(area~ summerMonth*temp)
qqnorm(reg.line$res,pch=16,main="Linear Model Residuals: QQ Plot")
qqline(reg.line$res,lwd = 2,col = 2)


table(FFMC)
barplot(table(FFMC))
table(DMC)
barplot(table(DMC))
round(table(FFMC)/length(FFMC),3)


source("panelfxns.R")
rm(vars)
#vars = cbind(X,Y,month,day,temp,RH,wind,rain)
vars = cbind(area,FFMC, DMC, DC ,ISI, temp,RH,wind,rain)
#vars = cbind(area, X,Y,month,day,FFMC, DMC, DC ,ISI, temp,RH,wind,rain)

pairs(vars);cov(vars)
pairs(vars,lower.panel=panel.cor)
