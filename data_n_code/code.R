#### A.3.1. Gathering Data ####
# load the dataset
df = read.csv('Data.csv', header = FALSE)

# based on the "Data Set Description" from the Website http://archive.ics.uci.edu/ml/datasets/Ozone+Level+Detection
# the last column, including 0 and 1, indicates normal day and ozone day respectively
# therefore, drop the last column
df <- df[,-c(74)]

# loading header
header = read.delim('Data_Details.txt', skip = 1, header = FALSE, sep = ":")

# keep the column name only
header <- header[, c(1)]

# naming the columns in the dataset
names(df) <- header
str(df)
View(df)





#### A.3.2 Data Preprocessing ####

# A.3.2.1
# replacing all "?" by missing values NA
df[df=="?"] <- NA
# converting all columns into numeric columnss, except from Date
for(i in 2:ncol(df)){
  df[[i]] <- as.numeric(df[[i]])
}

# A.3.2.2
# examine the correlations
library(corrplot)
cor.matrix <- cor(df[,-c(1)], use = 'pairwise.complete.obs')
corrplot(cor.matrix, tl.pos = 'n')

# A.3.2.3. reporting missing values
library(VIM)
library(mice)
# visualising missing values
missing_report <- aggr(df[,-c(1)], prop = F)
missing_report
length(missing_report$missings$Variable)

# A.3.2.4 imputing missing values
mice.df <- mice(df[,-c(1)], m=5, meth="pmm")
summary(mice.df)
complete_df <- complete(mice.df, 1)
complete_df$Date <- df$Date
# reorder the columns
df <- complete_df[,c(73, 1:72)]
str(df)

# A.3.2.5 turn df into time series dataframe
library(xts)
df <- xts(df[,-1], order.by = as.Date(df[,1], "%m/%d/%Y"))
str(df)
View(df)

# A.3.2.6 creating time series
var1 <- ts(df$WSR_PK, start = c(1998,1,1), frequency = 365.25)
var2 <- ts(df$T_PK, start = c(1998,1,1), frequency = 365.25)
var3 <- ts(df$T_AV, start = c(1998,1,1), frequency = 365.25)
var4 <- ts(df$T85, start = c(1998,1,1), frequency = 365.25)
var5 <- ts(df$RH85, start = c(1998,1,1), frequency = 365.25)
var6 <- ts(df$HT85, start = c(1998,1,1), frequency = 365.25)
var7 <- ts(df$T70, start = c(1998,1,1), frequency = 365.25)
var8 <- ts(df$KI, start = c(1998,1,1), frequency = 365.25)
var9 <- ts(df$TT, start = c(1998,1,1), frequency = 365.25)
var10 <- ts(df$SLP, start = c(1998,1,1), frequency = 365.25)
var11 <- ts(df$SLP_, start = c(1998,1,1), frequency = 365.25)



#### A.4. Results ####
library(astsa)
library(forecast)
library(ggplot2)
library(skimr)
library(summarytools)
farima <- function(x, h) {
  forecast(auto.arima(x), h = h)
}
#### A.4.1. WSR_PK (var1) ####
# A.4.1.1. Descriptive Statistics #
# summary
descr(var1)
# histogram
hist(var1,
     main = 'Histogram of WSR_PK', 
     xlab = "variable's unit",
     probability = T,
     xlim = c(min(var1) - sd(var1), max(var1) + sd(var1)),
     col = 'grey93')
lines(density(var1), col='red')
# time series plot
plot(var1, main = 'Time Series WSR_PK', ylab = "variable's unit")
abline(tslm(var1~time(var1)), col='red')
legend('topleft', 
       legend=c("WSR_PK", "Fitted Regression Line"),
       lty=1,
       col = c('black', 'red'),
       cex=0.8)

# A.4.1.2. Checking add or multi / Decomposing time series 
# Checking additive or Multiplicative
sum_squared_corr_add <- sum((acf(decompose(var1, type = 'additive')$random, na.action = na.omit)$acf)^2)
sum_squared_corr_add
sum_squared_corr_multi <- sum((acf(decompose(var1, type = 'multiplicative')$random, na.action = na.omit)$acf)^2)
sum_squared_corr_multi
if(sum_squared_corr_add < sum_squared_corr_multi){
  type <- 'additive'
}else{
  type <- 'multiplicative'
}
type
# Decomposing time series
plot(decompose(var1, type))

random <- decompose(var1, type = 'multiplicative')$random
Box.test(random,type = "Ljung-Box")


# A.4.1.3. Autocorrelation
ggAcf(var1, lag = length(var1)) + theme_classic() + ggtitle('WSR_PK Autocorrelation')


# A.4.1.4. Simple Exponential Smoothing
var1_ses <- ses(var1, h=365)
var1_ses$model
autoplot(var1_ses) + 
  autolayer(fitted(var1_ses)) + 
  theme_classic()
checkresiduals(var1_ses)
var1_ses_error <- tsCV(df$WSR_PK, ses, h =1) # df$WSR_PK is faster than var1 but the same result              
mean(var1_ses_error^2, na.rm = TRUE)

# A.4.1.5.  holt's linear method
var1_holt <- holt(var1, h=365)
var1_holt$model
autoplot(var1_holt) + autolayer(fitted(var1_holt)) + theme_classic()
checkresiduals(var1_holt)
var1_holt_error <- tsCV(df$WSR_PK, holt, h =1)             # df$WSR_PK is faster than var1 but the same result
mean(var1_holt_error^2, na.rm = TRUE)

# A.4.1.6.  damped trend
var1_dam <- holt(var1, h=365, damped = T)
var1_dam$model
autoplot(var1_dam) + autolayer(fitted(var1_dam)) + theme_classic()
checkresiduals(var1_dam)
var1_dam_error <- tsCV(df$WSR_PK, holt, damped=T, h =1)    # df$WSR_PK is faster than var1 but the same result
mean(var1_dam_error^2, na.rm = TRUE)

# A.4.1.7. arima
var1_arima <- auto.arima(var1)
var1_arima
autoplot(forecast(var1_arima, h=365)) + autolayer(fitted(var1_arima)) + theme_classic()
checkresiduals(forecast(var1_arima, h=365))
var1_arima_error <- tsCV(df$WSR_PK, farima, h = 1)        # df$WSR_PK is faster than var1 but the same result
mean(var1_arima_error^2, na.rm = TRUE)




#### A.4.2. GATHERING RESULTS ####

#### IMPORTANT: TO CHECK THE RESULTS OF ANY VARIABLE FROM var1, var2,... var11
#### replace the variable of interest below 
#### var <- ... 
var <- var11
var_name <- dimnames(var)[[2]]
var_name


# A.4.2.1. Descriptive Statistics #
# summary
descr(var)
# histogram
hist(var,
     probability = T,
     main = paste('Histogram of ', var_name),
     xlim = c(min(var) - sd(var), max(var) + sd(var)),
     col = 'grey93')
lines(density(var), col='red')
# time series plot
plot(var, main = paste('Time Series ', var_name), ylab = "variable's unit")
abline(tslm(var~time(var)), col='red')
legend('topleft', 
       legend=c(var_name, "Fitted Regression Line"),
       lty=1,
       col = c('black', 'red'),
       cex=0.8)

# A.4.1.2. Checking add or multi / Decomposing time series 
# Checking additive or Multiplicative
sum_squared_corr_add <- sum((acf(decompose(var, type = 'additive')$random, na.action = na.omit)$acf)^2)
sum_squared_corr_add
sum_squared_corr_multi <- sum((acf(decompose(var, type = 'multiplicative')$random, na.action = na.omit)$acf)^2)
sum_squared_corr_multi
if(sum_squared_corr_add < sum_squared_corr_multi){
  type <- 'additive'
}else{
  type <- 'multiplicative'
}
type
# Decomposing time series
plot(decompose(var, type))


# A.4.1.3. Autocorrelation
# Plotting ACF
ggAcf(var, lag = length(var)) + theme_classic() + ggtitle(paste(var_name,'Autocorrelation'))


# A.4.1.4. Simple Exponential Smoothing
# fitting the simple exponential smoothing model with 1-year prediction 
var_ses <- ses(var, h=365)
# getting the model's parameters
var_ses$model
# plotting the model's prediction and its fitted values 
autoplot(var_ses) + 
  autolayer(fitted(var_ses)) + 
  theme_classic()
# checking model's residuals - Ljung-Box test
checkresiduals(var_ses)
# 1-step cross validation 
var_ses_error <- tsCV(df[,colnames(df) == var_name], ses, h =1) 
# cross-validated mean squared error
mean(var_ses_error^2, na.rm = TRUE)


# A.4.1.5.  holt's linear method
# fitting the holt's linear method with 1-year prediction 
var_holt <- holt(var, h=365)
# getting the model's parameters
var_holt$model
# plotting the model's prediction and its fitted values 
autoplot(var_holt) + autolayer(fitted(var_holt)) + theme_classic()
# checking model's residuals - Ljung-Box test
checkresiduals(var_holt)
# 1-step cross validation
var_holt_error <- tsCV(df[,colnames(df) == var_name], holt, h =1) 
# cross-validated mean squared error
mean(var_holt_error^2, na.rm = TRUE)


# A.4.1.6.  damped trend method
# fitting damped trend method with 1-year prediction 
var_dam <- holt(var, h=365, damped = T)
# getting the model's parameters
var_dam$model
# plotting the model's prediction and its fitted values 
autoplot(var_dam) + autolayer(fitted(var_dam)) + theme_classic()
# checking model's residuals - Ljung-Box test
checkresiduals(var_dam)
# 1-step cross validation
var_dam_error <- tsCV(df[,colnames(df) == var_name], holt, damped=T, h =1)
# cross-validated mean squared error
mean(var_dam_error^2, na.rm = TRUE)

# A.4.1.7. arima
# fitting arima model
var_arima <- auto.arima(var)
# getting the model's parameters
var_arima
# plotting the model's 1-year prediction and its fitted values
autoplot(forecast(var_arima, h=365)) + autolayer(fitted(var_arima)) + theme_classic()
# checking model's residuals - Ljung-Box test
checkresiduals(forecast(var_arima, h=365))
# 1-step cross validation
var_arima_error <- tsCV(df[,colnames(df) == var_name], farima, h = 1) 
# cross-validated mean squared error
mean(var_arima_error^2, na.rm = TRUE)




