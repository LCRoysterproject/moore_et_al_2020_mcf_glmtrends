###########################################
## Oyster Transects - Data Management    ##
## Jennifer Moore, Bill Pine, Mel Moreno ##
## Updated: July 2019####                ##
###########################################

##a really nice feature is that this is the file for use in reducing the length of the transects
##that have the -999 for specific segments (2.5 m intervals) that result from the transect flooding
##before counts can be completed.

##note also this file averages the counts from the two observers

##at the bottom is a small "other" summary that matches the white board summaries of 
#transects conducted


#load packages
library("lubridate") #format dates

#read in data
#transect production file
#could change this to read directly from excel workbook
tr <- read.csv("data/transect_data.csv", header = T)

#format date column as a date object
tr$date<-mdy(tr$date)

#create season column
tr$season<-ifelse(tr$month ==1 | tr$month ==2 |tr$month ==3 | tr$month == 10 | tr$month == 11 | tr$month == 12, "winter", "summer")

#remove rows with NA for count_live
tr <- tr[!is.na(tr$count_live),]

#remove rows with NA for tran_length
tr <-tr[!is.na(tr$tran_length),]

###need to create periods of time when the sampling can be collapsed
###as an example "winter 2018/2019" would be period X
#I've examined all sampling and here are the periods

#here is the pattern for the periods.  This pattern would extend through 2019
#period 1 = year =2010 and months 5-9
#period 2 = year = 2010 months 10-12 and year = 2011 months 1-3
#period 3 = year = 2011 months 4-9
#period 4 = year =2011 months = 10-12 and year = 2012 months 1-3 
#currently 17 periods total

#create a period column based on the information above
tr$period <- NA
years <- sort(rep(unique(tr$year), times = 1, each = 2))
for(i in unique(years)){
  print(i)
  y <- i #year
  p <- which(years == i) #period number - 2010 = 1 and 2, 2011 = 3 and 4, and so forth.
  for(j in 1:nrow(tr)){
    if(tr$year[j] == y & tr$month[j] > 3 & tr$month[j] < 10) tr$period[j] = p[1] #year i months 4-9
    if(tr$year[j] == y & tr$month[j] > 9) tr$period[j] = p[2] #year i months 10-12
    if(tr$year[j] == y+1 & tr$month[j] < 4) tr$period[j] = p[2] #year i+1 months 1-3
  }
}




#only need to use locality == HB, CR, LC, CK for epoch 1-3 analyses

#only locality = LC, HB, CR, CK
tr <- tr[tr$locality != "BT",]
tr <- tr[tr$locality != "LT",]
tr <- tr[tr$locality != "NN",]
tr$locality = droplevels(tr$locality) #this removes those levels from the dataset


table(tr$station,tr$period)

#use BT, LC, LT, NN for epoch 3 analyses

# #only locality = BT, LC, LT, NN
# tr <- tr[tr$locality != "CK",]
# tr <- tr[tr$locality != "CR",]
# tr <- tr[tr$locality != "HB",]
# tr$locality = droplevels(tr$locality) #this removes those levels from the dataset



#only years after 2017
#tr <- tr[tr$year > 2017,]

#this removes the odd data from Jan 2018 with the grub box
tr<- tr[tr$date != "2018-01-30",]

#for 2018/2019 there are 2 passes (pass 1, pass 2)
#average values together for this - can change this later if we calculate detection probabilities

#keep day month year here because you want it to be averaged
#on the same day the two counts were done
dta0=aggregate(count_live~date+day+month+year+season+period+treatment+
                 locality+site+bar+station+transect+tran_length,data = tr, FUN = "mean")

#aggregate live count data for each transect
#first remove all rows with -999 then sum live count for each transect segment
dta0.2 <- dta0[dta0$count_live > -1,]

#oyster live counts by transect


#make sure year is here, maybe month too, but some counts go between months, so drop month for now
live=aggregate(count_live~year+season+period+treatment+locality+site+bar+station,data=dta0.2,sum)


#aggregate transect length for each transect
#for each row with -999 reduce transect length by 2.5

#max length for each transect
dta2=aggregate(tran_length~year+season+period+treatment+locality+site+bar+station+transect,data=dta0,max)

#find rows with -999
miss <- which(dta0$count_live < -1)
for(i in 1:length(miss)){
  ind <- miss[i]
  station <- dta0$station[ind]
  transect <- dta0$transect[ind]
  #subtract 2.5 from the transect length for this particular station/transect that has a missing value
  dta2$tran_length[dta2$station == station & dta2$transect == transect] <- dta2$tran_length[dta2$station == station & dta2$transect == transect] - 2.5
}

#sum over all transects
tranlength=aggregate(tran_length~year+season+period+treatment+locality+site+bar+station,data=dta2,sum)

#merge live count total data frame with the tran_length total data frame
#dta3=merge(live,tranlength,by=c("day","month","year","season","treatment","locality","site","bar","station"))
dta3=merge(live,tranlength,by=c("year","season","period","treatment","locality","site","bar","station"))


#calculate density
dta3$area = dta3$tran_length*.1524
dta3$density = dta3$count_live/dta3$area

dta3$area <- round(dta3$area,digits=2)
dta3$density <- round(dta3$density,digits=2)

#read in strata file
st <- read.csv("data/archive/strata.csv", header = T)
#attach strata data to final file
for(i in 1:nrow(dta3)){
  station <- as.character(dta3$station[i])
  strata <- as.character(st$strata[st$station == station])
  ifelse(length(strata) == 0, dta3$strata[i] <- NA, dta3$strata[i] <- strata)
}

#sort data by period
sort.dta3<-dta3[order(dta3$period),]



#write cleaned production file to .csv
write.csv(sort.dta3,file="data/all_epoch_transect_data.csv")

###
##here is a summary table to match the white board totals by strata
dtax=aggregate(count_live~day+month+year+season+treatment+locality+site+bar+station+transect,data=dta0,sum)
st <- read.csv("data/archive/strata.csv", header = T)
for(i in 1:nrow(dtax)){
  station <- as.character(dtax$station[i])
  print(station)
  strata <- as.character(st$strata[st$station == station])
  print(strata)
  ifelse(length(strata) == 0, dtax$strata[i] <- NA, dtax$strata[i] <- strata)
}

whiteboard<-table(dtax$strata)
write.table(whiteboard,file="data_output/whiteboard.txt", sep = ",", quote = FALSE, row.names = F)



#look at each strata individually
dtax_reduced<- dtax[dtax$strata == "N_SM",]

