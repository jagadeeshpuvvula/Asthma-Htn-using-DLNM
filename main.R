library(dlnm)
library(splines)
library(zoo)
library(dplyr)
library(mgcv)
library(lubridate)
library(dplyr)
library(naniar)
library(mice)
library(pcaMethods)
dat <- read.csv ("C:/Users/jagad/Desktop/douglas_cnty_asth_htn_ER/criteria_poll/analytic_dataset.csv", 
                 header = T,fileEncoding="UTF-8-BOM")

#temporal variables
dat$date <- as.Date(dat$date, format = "%m/%d/%Y")
dat$month<- as.factor(month(dat$date))
dat$year<- as.factor(format(dat$date, '%Y'))
dat$dow <- wday(as.Date(dat$date, format = "%m/%d/%Y"))
weekdays1 <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday')
dat$wDay <- factor((weekdays(dat$date) %in% weekdays1), 
                   levels=c(FALSE, TRUE), labels=c('weekend', 'weekday'))
dat$week_num<- week(ymd(dat$date))
#year to quarter conversion use zoo library for as.nearmon function
yq <- as.yearqtr(as.yearmon(dat$date, "%Y %m/%d/") + 1/12)
dat$season <- factor(format(yq, "%q"), levels = 1:4,
                     labels = c("winter", "spring", "summer", "fall"))

#convert all numeric variables to numeric
con.names = dat %>% select_if(is.numeric) %>% colnames()
dat[,con.names] = data.frame(apply(dat[con.names], 2, as.numeric))

dat$dow<- as.factor(dat$dow)
dat$dow<- as.factor(dat$dow)

str(dat)
vis_miss(dat_max)

dat_max<- subset(dat, select = c(asthma_ped, asthma_all, htn_all,
                                      Co.1Hr.Max, No.Max,Noy.Max, Noyno.Max, 
                                      Ozone.1Hr.Max, Pm1025Lc.Max,
                                      pm25_Max, so2_1hr_max,so2_5min_max, month,
                                      year, dow, wDay, week_num, season, pres_burn, 
                                      tmax, date))

#imputing missing values using MCAR (missing completely at random) assumption
#pmm : predictive mean matching
#m=5 creates five imputations
tempData <- mice(dat_max,m=5,maxit=50,meth='pmm',seed=500)
#view imputations
tempData$imp$Ozone.1Hr.Max

completedData <- complete(tempData,3)
vis_miss(completedData)

write.csv(completedData, "C:/Users/jagad/Desktop/douglas_cnty_asth_htn_ER/criteria_poll/mice_imputation.csv", 
          row.names = F)


#PPCA vs MICE imputation: https://www.sciencedirect.com/science/article/pii/S2352914819302783
pc <- pca(dat_max, nPcs=3, method="ppca")
pca_imputed <- completeObs(pc)

write.csv(pca_imputed, "C:/Users/jagad/Desktop/douglas_cnty_asth_htn_ER/criteria_poll/pca_imputed.csv", 
          row.names = F)

write.csv(dat_max, "C:/Users/jagad/Desktop/douglas_cnty_asth_htn_ER/criteria_poll/max_pollu.csv", 
          row.names = F)


##### Working with imputed data
dat <- read.csv ("C:/Users/jagad/Desktop/douglas_cnty_asth_htn_ER/criteria_poll/pca_imputed.csv", 
                 header = T,fileEncoding="UTF-8-BOM")


dat$date <- as.Date(dat$date, format = "%m/%d/%Y")
dat$month<- as.factor(month(dat$date))
dat$year<- as.factor(format(dat$date, '%Y'))
dat$dow <- wday(as.Date(dat$date, format = "%m/%d/%Y"))
weekdays1 <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday')
dat$wDay <- factor((weekdays(dat$date) %in% weekdays1), 
                   levels=c(FALSE, TRUE), labels=c('weekend', 'weekday'))
dat$week_num<- week(ymd(dat$date))
#year to quarter conversion use zoo library for as.nearmon function
yq <- as.yearqtr(as.yearmon(dat$date, "%Y %m/%d/") + 1/12)
dat$season <- factor(format(yq, "%q"), levels = 1:4,
                     labels = c("winter", "spring", "summer", "fall"))

dat$season<- as.factor (dat$season)
dat$wDay<- as.factor (dat$wDay)
dat$dow<- as.numeric(dat$dow)

str(dat)

#convert all numeric variables to numeric
con.names = dat %>% select_if(is.numeric) %>% colnames()
dat[,con.names] = data.frame(apply(dat[con.names], 2, as.numeric))

vis_miss(dat)
summary(dat)

#DLNM MODEL for time searies data

#continuous TS
cb1.pm <- crossbasis(dat$pm25_Max, lag=7, argvar=list(fun="lin"),
                     arglag=list(fun="poly",degree=4))
cb1.temp <- crossbasis(dat$tmax, lag=3, argvar=list(df=5),
                       rglag=list(fun="strata",breaks=1))

model1 <- glm(asthma_ped ~ cb1.pm + cb1.temp + ns(date, 7*4) + dow,
              family=quasipoisson(), dat)

pred1.pm <- crosspred(cb1.pm, model1, at=0:20, bylag=0.2, cumul=TRUE)

plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
     main="Association with a 10-unit increase in PM10")
plot(pred1.pm, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR",
     main="Cumulative association with a 10-unit increase in PM10")

#stratified analysis by burn season

dat_pres <- subset(dat, pres_burn %in% 1)

cb2.o3 <- crossbasis(dat_pres$pm25_Max, lag=5,
                     argvar=list(fun="thr",thr=15), arglag=list(fun="integer"),
                     group=dat_pres$year)
cb2.temp <- crossbasis(dat_pres$tmax, lag=10,argvar=list(fun="thr",thr=c(15,25)), 
                       arglag=list(fun="strata",breaks=c(2,6)),
                         group=dat_pres$year)

#ns(doy, 4) + removed (creat a variables for day of year)
model2 <- glm(asthma_ped ~ cb2.o3 + cb2.temp + ns(date,3) + dow + season + year,
              family=quasipoisson(), dat_pres)

pred2.o3 <- crosspred(cb2.o3, model2, at=c(0:25,15,17))

plot(pred2.o3, "slices", var=17, ci="bars", type="p", col=2, pch=19,
       ci.level=0.80, main="Lag-response a 10-unit increase above threshold (80CI)")

plot(pred2.o3,"overall",xlab="Ozone", ci="l", col=3, ylim=c(0.9,1.3), lwd=2,
       ci.arg=list(col=1,lty=3), main="Overall cumulative association for 5 lags")

#REF 6-A figure
#http://www.ag-myresearch.com/uploads/1/3/8/6/13864925/2019_vicedo-cabrera_epidem.pdf
#https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata/blob/master/07DemChangesAdapt.R