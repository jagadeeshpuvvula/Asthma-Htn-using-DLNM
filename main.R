library(dlnm)
library(splines)
library(zoo)
library(dplyr)
library(mgcv)
library(lubridate)
library(dplyr)

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