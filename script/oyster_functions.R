#############################
## Oyster Data Management  ##
## Output Summary Tables   ##
## Plots                   ##
## Analyses                ##
## FUNCTIONS               ##
##                         ##
## Jennifer Moore          ##
## October 29, 2019        ##
#############################

#organizeData()
#Function to organize raw data
#formats dates
#changes missing data to NA
#add columns for season, period, strata, rocksBinary, rocks, harvest
organizeData <- function(data){
  
  #load packages
  library("lubridate") #format dates
  
  #format date column as a date object
  data$date <- mdy(data$date)
  
  #create season column and make it a factor
  data$season <- ifelse(data$month ==1 | data$month ==2 |data$month ==3 | data$month == 10 | data$month == 11 | data$month == 12, "winter", "summer")
  data$season <- as.factor(data$season)
  
  #remove rows with -999 for count_live
  data <- subset(data, data$count_live > -1)
  
  #remove rows with -999 for tran_length
  data <- subset(data, data$tran_length > -1)
  
  #remove time columns because it is missing alot
  data <- data[,-c(6,7)]
  
  #add period column
  #here is the pattern for the periods.  This pattern would extend through the end of sampling
  #period 1 = year = 2010 and months 4-9
  #period 2 = year = 2010 months 10-12 and year = 2011 months 1-3
  #period 3 = year = 2011 months 4-9
  #period 4 = year = 2011 months = 10-12 and year = 2012 months 1-3 
  #period 5 = year = 2012 months 4-9
  #period 6 = year = 2012 months = 10-12 and year = 2013 months 1-3
  
  data$period <- NA
  firstyear <- min(data$year)
  endyear <- max(data$year)
  years <- sort(rep(firstyear:endyear, times = 1, each = 2))
  for(i in unique(years)){
    y <- i #year
    p <- which(years == i) #period number - 2010 = 1 and 2, 2011 = 3 and 4, and so forth.
    for(j in 1:nrow(data)){
      if(data$year[j] == y & data$month[j] > 3 & data$month[j] < 10) data$period[j] = p[1] #year i months 4-9
      if(data$year[j] == y & data$month[j] > 9) data$period[j] = p[2] #year i months 10-12
      if(data$year[j] == y+1 & data$month[j] < 4) data$period[j] = p[2] #year i+1 months 1-3
    }
  }
  
  #this removes the odd data from Jan 2018 with the grub box
  data <- data[data$date != "2018-01-30",]
  
  #read in strata file
  st <- read.csv("data/strata_period.csv", header = T)
  
  #attach strata as a new column
  data$strata <- NA
  for(i in 1:nrow(data)){
    station <- as.character(data$station[i])
    period <- data$period[i]
    strata <- as.character(st$strata[st$station == station & st$period == period])
    ifelse(length(strata) == 0, data$strata[i] <- NA, data$strata[i] <- strata)
  }
  
  #add rocks column (second half of strata) - either LG, SM, or NA
  data$rocks <- substr(data$strata,3,4)
  
  #add harvest column (first half of strata) - either Y or N
  data$harvest <- substr(data$strata,1,1)
  
  return(data)
}

#calculateCountsDensity()
#Function to calculate live counts and density by transect
#aggregates over each transect
#averages 2 passes (if there are two)
#add columns for live count, transect length, area, and density
calculateCountsDensity <- function(data, data2) {
  
  #average the two passes 
  dta1=aggregate(count_live~date+day+month+year+season+period+treatment+
                   locality+site+bar+station+transect+tran_length+strata+rocks+harvest,data = data2, FUN = "mean")
  
  #first remove all rows with -999 
  dta2 <- dta1[dta1$count_live > -1,]
  
  #sum live counts for each transect
  live=aggregate(count_live~season+period+treatment+locality+site+bar+station+strata+rocks+harvest,data=dta2,sum)
  
  #aggregate transect length for each transect
  #max length for each transect
  length=aggregate(tran_length~season+period+treatment+locality+site+bar+station+transect+strata+rocks+harvest,data=dta2,max)
  
  #add period in the data file so that we can compare with the new file
  data$period <- NA
  firstyear <- min(data$year)
  endyear <- max(data$year)
  years <- sort(rep(firstyear:endyear, times = 1, each = 2))
  for(i in unique(years)){
    y <- i #year
    p <- which(years == i) #period number - 2010 = 1 and 2, 2011 = 3 and 4, and so forth.
    for(j in 1:nrow(data)){
      if(data$year[j] == y & data$month[j] > 3 & data$month[j] < 10) data$period[j] = p[1] #year i months 4-9
      if(data$year[j] == y & data$month[j] > 9) data$period[j] = p[2] #year i months 10-12
      if(data$year[j] == y+1 & data$month[j] < 4) data$period[j] = p[2] #year i+1 months 1-3
    }
  }
  #remove rows with tran_length = -999
  data <- subset(data, data$tran_length > -1)
  
  #find rows with -999
  miss <- which(data$count_live < -1)
  for(i in 1:length(miss)){
    ind <- miss[i]
    period <- data$period[ind]
    station <- data$station[ind]
    transect <- data$transect[ind]
    #subtract 2.5 from the transect length for this particular station/transect that has a missing value
    length$tran_length[length$station == station & length$transect == transect & length$period == period] <- length$tran_length[length$station == station & length$transect == transect & length$period == period] - 2.5
  }
  
  #sum over all transects
  tranlength=aggregate(tran_length~season+period+treatment+locality+site+bar+station+strata+rocks+harvest,data=length,sum)
  
  #merge live count total data frame with the tran_length total data frame
  dta3=merge(live,tranlength,by=c("season","period","treatment","locality","site","bar","station","strata","rocks","harvest"))
  
  #calculate density
  dta3$area = dta3$tran_length*.1524
  dta3$density = dta3$count_live/dta3$area
  
  dta3$area <- round(dta3$area,digits=2)
  dta3$density <- round(dta3$density,digits=2)
  
  return(dta3)
}

#calculates counts and densities separately for each pass
doublePasses <- function(data, data2){
  
  #first remove all rows with -999 
  dta1 <- data2[data2$count_live > -1,]
  
  #sum live counts for each transect
  live=aggregate(count_live~season+period+treatment+locality+site+bar+station+strata+rocks+harvest+pass,data=dta1,sum)
  
  #aggregate transect length for each transect
  #max length for each transect
  length=aggregate(tran_length~season+period+treatment+locality+site+bar+station+transect+strata+rocks+harvest+pass,data=dta1,max)
  
  #add period in the data file so that we can compare with the new file
  data$period <- NA
  firstyear <- min(data$year)
  endyear <- max(data$year)
  years <- sort(rep(firstyear:endyear, times = 1, each = 2))
  for(i in unique(years)){
    y <- i #year
    p <- which(years == i) #period number - 2010 = 1 and 2, 2011 = 3 and 4, and so forth.
    for(j in 1:nrow(data)){
      if(data$year[j] == y & data$month[j] > 3 & data$month[j] < 10) data$period[j] = p[1] #year i months 4-9
      if(data$year[j] == y & data$month[j] > 9) data$period[j] = p[2] #year i months 10-12
      if(data$year[j] == y+1 & data$month[j] < 4) data$period[j] = p[2] #year i+1 months 1-3
    }
  }
  #remove rows with tran_length = -999
  data <- subset(data, data$tran_length > -1)
  
  #find rows with -999
  miss <- which(data$count_live < -1)
  for(i in 1:length(miss)){
    ind <- miss[i]
    period <- data$period[ind]
    station <- data$station[ind]
    transect <- data$transect[ind]
    pass <- data$pass[ind]
    #subtract 2.5 from the transect length for this particular station/transect that has a missing value
    length$tran_length[length$station == station & length$transect == transect & length$period == period & length$pass == pass] <- length$tran_length[length$station == station & length$transect == transect & length$period == period] - 2.5
  }
  
  #sum over all transects
  tranlength=aggregate(tran_length~season+period+treatment+locality+site+bar+station+strata+rocks+harvest+pass,data=length,sum)
  
  #merge live count total data frame with the tran_length total data frame
  dta3=merge(live,tranlength,by=c("season","period","treatment","locality","site","bar","station","strata","rocks","harvest","pass"))
  
  #calculate density
  dta3$area = dta3$tran_length*.1524
  dta3$density = dta3$count_live/dta3$area
  
  dta3$area <- round(dta3$area,digits=2)
  dta3$density <- round(dta3$density,digits=2)
  
  return(dta3)
}

#display summary tables
#total number of transects walked per station, per strata, and per period
summaryEffort <- function(data){
  cat("Effort by Locality")
  line <- readline()
  a= (aggregate(count_live ~ locality, data = data, FUN = 'length'))
  a2<- aggregate(tran_length ~ locality, data= data, FUN = 'sum')
  a3<- merge (a, a2, by= "locality")
  colnames(a3) = c("Locality", "Number of Transects", "Total Length (m)")
  print(a3, row.names = FALSE)
  line <- readline()
  
  cat("Effort by Strata")
  line <- readline()
  b= (aggregate(count_live ~ strata, data = data, FUN = 'length'))
  b2<- aggregate(tran_length ~ strata, data= data, FUN = 'sum')
  b3<- merge (b, b2, by= "strata")
  colnames(b3) = c("Strata", "Number of Transects", "Total Length (m)")
  print(b3, row.names = FALSE)
  line <- readline()
  
  cat("Effort by Period")
  line <- readline()
  c=(aggregate(count_live ~ period, data = data, FUN = 'length'))
  c2<- aggregate(tran_length ~ period, data= data, FUN = 'sum')
  c3<- merge (c, c2, by= "period")
  colnames(c3) = c("Period",  "Number of Transects", "Total Length (m)" )
  print(c3, row.names = FALSE)
  line <- readline()
  
  cat("Effort by Locality and Period")
  line <- readline()
  d=(aggregate(count_live ~ period+locality, data = data, FUN = 'length'))
  d2<- aggregate(tran_length ~ period + locality, data= data, FUN = 'sum')
  d3<- merge (d, d2)
  colnames(d3) = c("Period", "Locality", "Number of Transects", "Total Length (m)" )
  print(d3, row.names = FALSE)
  line <- readline()
  
  cat("Effort by Strata and Period")
  line <- readline()
  e=(aggregate(count_live ~ period+strata, data = data, FUN = 'length'))
  e2<- aggregate(tran_length ~ period+strata, data= data, FUN = 'sum')
  e3<- merge (e, e2)
  colnames(e3) = c("Period", "Strata", "Number of Transects", "Total Length (m)" )
  print(e3, row.names = FALSE)
  line <- readline()
}

#display summary tables 
#total counts per station, per strata, and per period
summaryCounts <- function(data){
  cat("Total Counts by Locality")
  line <- readline()
  a <- aggregate(count_live ~ locality, data = data, FUN = function (x) sumstats(x))
  a2 <- as.data.frame(a$count_live)
  a2$Locality <- a$locality
  a2 <- a2[,c(12, 1:11)]
  colnames(a2) <- c("Locality", "Mean", "Median", "SD", "Var","CV","SE", "L95", "U95", "Bstrap_Mean", "L95_Bstrap", "U95_Bstrap")
  print(a2, row.names = FALSE)
  line <- readline()
  
  cat("Total Counts by Strata")
  line <- readline()
  b <- aggregate(count_live ~ strata, data = data, FUN = function (x) sumstats(x))
  b2 <- as.data.frame(b$count_live)
  b2$Strata <- b$strata
  b2 <- b2[,c(12, 1:11)]
  colnames(b2) <- c("Strata", "Mean", "Median", "SD", "Var","CV","SE", "L95", "U95", "Bstrap_Mean", "L95_Bstrap", "U95_Bstrap")
  print(b2, row.names = FALSE)
  line <- readline()
  
  cat("Total Counts by Period")
  line <- readline()
  c <- (aggregate(count_live ~ period, data = data, FUN = function (x) sumstats(x)))
  c2 <- as.data.frame(c$count_live)
  c2$Period <- c$period
  c2 <- c2[,c(12, 1:11)]
  colnames(c2) <- c("Period", "Mean", "Median", "SD", "Var","CV","SE", "L95", "U95", "Bstrap_Mean", "L95_Bstrap", "U95_Bstrap")
  print(c2, row.names = FALSE)
  line <- readline()
}

#display summary tables 
#average density per station, per strata, and per period
summaryDensity <- function(data){
  cat("Density by Locality")
  line <- readline()
  a <- aggregate(density ~ locality, data = data, FUN = function (x) sumstats(x))
  a2 <- as.data.frame(a$density)
  a2$Locality <- a$locality
  a2 <- a2[,c(12, 1:11)]
  colnames(a2) <- c("Locality", "Mean", "Median", "SD", "Var","CV","SE", "L95", "U95", "Bstrap_Mean", "L95_Bstrap", "U95_Bstrap")
  print(a2, row.names = FALSE)
  line <- readline()
  
  cat("Density by Strata")
  line <- readline()
  b <- aggregate(density ~ strata, data = data, FUN = function (x) sumstats(x))
  b2 <- as.data.frame(b$density)
  b2$Strata <- b$strata
  b2 <- b2[,c(12, 1:11)]
  colnames(b2) <- c("Strata", "Mean", "Median", "SD", "Var","CV","SE", "L95", "U95", "Bstrap_Mean", "L95_Bstrap", "U95_Bstrap")
  print(b2, row.names = FALSE)
  line <- readline()
  line <- readline()
  
  cat("Density by Period")
  line <- readline()
  c <- (aggregate(density ~ period, data = data, FUN = function (x) sumstats(x)))
  c2 <- as.data.frame(c$density)
  c2$Period <- c$period
  c2 <- c2[,c(12, 1:11)]
  colnames(c2) <- c("Period", "Mean", "Median", "SD", "Var","CV","SE", "L95", "U95", "Bstrap_Mean", "L95_Bstrap", "U95_Bstrap")
  print(c2, row.names = FALSE)
  line <- readline()
}

#display plots
#frequency plot for density by locality
#boxplots for average density per strata and per period
plotsDensity <- function(data){
  print(ggplot(cal_den, aes(x=density, fill = locality)) +
          geom_density(size = 2)+
          xlab("Density")+
          ylab("Probability Density Function")+
          labs(title = "Oyster Density by Locality", fill = "Locality",
               caption= "Figure- Calculated oyster density by locality for all periods including period 20 (current period).")+
          facet_wrap(~locality, scales = 'free', ncol = 2)+
          scale_x_continuous(limits=c(0,1500)) +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  cat('\r\n\r\n')
  
  print(ggplot(cal_den, aes(x=density, fill = strata)) +
          geom_density(size = 2)+
          xlab("Density")+
          ylab("Probability Density Function")+
          labs(title = "Oyster Density by Strata", fill = "Strata",
               caption = "Figure- Calculated oyster density by strata for all periods including period 20 (current period).")+
          facet_wrap(~strata, scales = 'free', ncol = 2)+
          scale_x_continuous(limits=c(0,1500)) +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  cat('\r\n\r\n')
  
  print(ggplot(cal_den, aes(x=density, fill = as.factor(period))) +
          geom_density(size = 2)+
          xlab("Density")+
          ylab("Probability Density Function")+
          labs(title = "Oyster Density by Period", fill = "Period",
               caption = "Figure- Calculated oyster density for all periods including period 20 (current period) using a probability density function.")+
          facet_wrap(~period, scales = 'free', ncol = 2)+
          scale_x_continuous(limits=c(0,1500)) +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  cat('\r\n\r\n')
  
  print(ggplot(data, aes(locality, density))+
          geom_boxplot()+ coord_flip()+
          labs(title = "Oyster Density by Locality", 
               caption = "Figure- Box plot depicting density by locality for all periods including period 20 (current period).")+
          xlab("Locality") +
          ylab ("Oyster density per m^2") +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  cat('\r\n\r\n')
  
  print(ggplot(data, aes(strata, density))+
          geom_boxplot()+ coord_flip()+
          labs(title = "Oyster Density by Strata", 
               caption = "Figure- Box plot depicting density by strata for all periods including period 20 (current period).")+
          xlab("Strata") +
          ylab ("Oyster density per m^2") +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  cat('\r\n\r\n')  
  
  print(ggplot(data, aes(as.factor(period), density))+
          geom_boxplot()+ coord_flip()+
          labs(title = "Oyster Density by Period", 
               caption = "Figure- Box plot depicting density by period for all periods including period 20 (current period).")+
          xlab("Period") +
          ylab ("Oyster density per m^2") +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  cat('\r\n\r\n')  
  
  print(ggplot(data, aes(density, period, shape=locality, colour=locality))+
          geom_point(size=5, alpha=0.5)+
          scale_shape_manual(values = c(15,16,17,18,19,3,8))+
          scale_color_manual(values = c("#E69F00", "#000000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
          labs(title = "Oyster Density by Locality and Period", shape = "Locality", colour="Locality", caption = "Figure - Oyster density by locality and period for all periods including period 20 (current period). ")+
          scale_y_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")) +
          ylab("Period") +
          xlab ("Oyster density per m^2") +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  cat('\r\n\r\n')  
  
  print(ggplot(data, aes(density, period, shape=strata, colour=strata))+
          geom_point(size=5, alpha=0.5) +
          scale_shape_manual(values = c(15,16,17,18,3,8))+
          scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"))+
          labs(title = "Oyster Density by Strata and Period", shape = "Strata", colour="Strata", caption = "Figure - Oyster density by strata and period for all periods including period 20 (current period). ") +
          scale_y_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")) +
          ylab("Period") +
          xlab ("Oyster density per m^2") +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
}

#pilot sites
#LCO10B, LCO11A, LCO8B, LCO9A
#plots for oyster counts and density before and after
pilotSites <- function(data){
  data2 <- subset(data, data$station == "LCO10B" | data$station == "LCO11A" |
                    data$station == "LCO8B" | data$station == "LCO9A")
  
  print(ggplot(data2, aes(density, period, shape=station, colour=station))+
          geom_point(size=5, alpha=0.5)+
          scale_shape_manual(values = c(15,16,17,18,19,3,8))+
          scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
          labs(title = "Average Density by Station and Period", shape = "Station", colour="Station", caption = " Figure - Average density comparison by period for all stations that were sampled during the pilot study.")+
          ylab("Period") +
          xlab ("Oyster density per m^2") +
          scale_y_continuous(breaks = seq(7, 16, by = 1)) +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0))
  )
  line <- readline()
  
  
  print(ggplot(data2, aes(density, strata, shape=station, colour=station))+
          geom_point(size=5, alpha=0.5)+
          scale_shape_manual(values = c(15,16,17,18,19,3,8))+
          scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
          labs(title = "Average Density by Station and Strata", shape = "Station", colour = "Station", caption = "Figure - Average density comparison by strata and period for all stations that were sampled during the pilot study.")+
          ylab("Strata") +
          xlab("Oyster density per m^2") +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  line <- readline()
  
}


effortPlot<- function(data) {
  
  a= (aggregate(count_live ~ locality, data = data, FUN = 'length'))
  a2<- aggregate(tran_length ~ locality, data= data, FUN = 'sum')
  a3<- merge (a, a2, by= "locality")
  colnames(a3) = c("locality", "num_transects", "total_length")
  
  print(ggplot(a3) +
          geom_bar(aes(y= total_length, x= locality),stat = "identity") +
          theme_bw() +
          ylab("Total Length of Transects (m)") +
          xlab("Locality") +
          labs(title = "Total Transect Length Sampled by Locality", caption = "Figure - Bar plot of total transect length in meters sampled by locality for all periods.") +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  line <- readline()
  
  
  b= (aggregate(count_live ~ strata, data = data, FUN = 'length'))
  b2<- aggregate(tran_length ~ strata, data= data, FUN = 'sum')
  b3<- merge (b, b2, by= "strata")
  colnames(b3) = c("strata", "num_transects", "total_length")
  line <- readline()
  
  
  print(ggplot(b3) +
          geom_bar(aes(y= total_length, x= strata),stat = "identity") +
          theme_bw() +
          ylab("Total Length of Transects (m)") +
          xlab ("Strata") +
          labs(title = "Total Transect Length Sampled by Strata", 
               caption= "Figure - Bar plot of total transect length in meters sampled by strata for all periods.")+
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
  line <- readline()
  
  
  c=(aggregate(count_live ~ period, data = data, FUN = 'length'))
  c2<- aggregate(tran_length ~ period, data= data, FUN = 'sum')
  c3<- merge (c, c2, by= "period")
  colnames(c3) = c("period",  "num_transects", "total_length" )
  
  print(ggplot(c3) +
          geom_bar(aes(y= total_length, x= period), stat = "identity") +
          theme_bw() +
          ylab("Total Length of Transects (m)") +
          xlab ("Period")+
          labs(title = "Total Transect Length Sampled by Period", 
               caption= "Figure- Bar plot of total transect length in meters sampled by period for all periods.") +
          scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")) +
          theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA,size=1, linetype="solid"), plot.caption = element_text(hjust = 0)))
}


#progress bar
#number of proposed sites
#number of completed sites as of last date in data file
#input organized data file - not aggregated
progress <- function(data){
  s <- subset(data, data$period == 20)
  f <- function(x){length(unique(x))}
  s2 <- aggregate(transect ~ station+treatment, data = s, FUN = f)
  s3 <- aggregate(transect ~ treatment, data = s2, FUN = 'sum')
  
  totalRocks <- 66
  totalControl <- 40
  
  rockBar <- (s3$transect[s3$treatment == "rocks"] / totalRocks) * 100
  controlBar <- (s3$transect[s3$treatment == "control"] / totalControl) * 100
  
  plot(c(0,100), c(0,2), type = 'n', xlab = 'Percentage Complete', ylab = '', yaxt = 'n', mar=c(3,3,3,3))
  #plot rock sites progress
  rect(0, 0.1+1-1, 100, 0.9+1-1)
  rect(0, 0.1+1-1, rockBar, 0.9+1-1, col = 'blue')
  text(40, 0.5+1-1, paste('Rock Sites: ', round(rockBar,2), '%', sep=''), adj = 0, col = 'black')
  #plot controls ites progress
  rect(0, 0.1+2-1, 100, 0.9+2-1)
  rect(0, 0.1+2-1, controlBar, 0.9+2-1, col = 'orange')
  text(40, 0.5+2-1, paste('Control Sites: ', round(controlBar,2), '%', sep=''), adj = 0, col = 'black')
  title('Field Work Progress')
}

#compare most recent years (this year to last year)
#right now this is period 20 to period 19 (winter 2019-2020 to summer 2019)
summaryRecent <- function(data){
  data2 <- subset(data, data$period > 18)
  print("Total Counts per Period")
  print(aggregate(count_live ~ period, FUN = "sum", data = data2))
  print("Total Counts per Period per Strata")
  print(aggregate(count_live ~ period + strata, FUN = 'sum', data = data2))
  print("Average Density per Period")
  print(aggregate(density ~ period, FUN = 'mean', data = data2))
  print("Average Density per Period per Strata")
  print(aggregate(density ~ period + strata, FUN = 'mean', data = data2))
}

#subset NFWF
#LC, NN, BT - any year
subsetNFWF <- function(data){
  data2 <- subset(data, data$locality == "LC" | data$locality == "NN" | data$locality == "BT")
  data2$locality <- droplevels(data2$locality)
  return(data2)
}

#subset repeated measures sites  
#LC09C, LCO10A, LCO11B, LCO12 (control)
#LCO10B, LCO11A, LCO8B, LCO9A (pilot rocks)
subsetRepeatedMeasures <- function(data){
  data2 <- subset(data, data$station == "LCO9C" | 
                    data$station == "LCO10A" | data$station == "LCO11B" | 
                    data$station == "LCO12" | data$station == "LCO10B" |
                    data$station == "LCO11A" | data$station == "LCO8B" |
                    data$station == "LCO9A")
  data3 <- subset(data2, data2$tran_length > -1)
  data4 <- subset(data3, data3$count_live > -1)
  return(data2)
}


## A function to find the start date of dates
start_date = function(year) {
  # Finds start state for this calendar year
  jan1 = as.Date(paste(year, '-01-01', sep=''))
  wday = as.numeric(MMWRweekday(jan1))
  jan1 - (wday-1) + 7*(wday>4)
}


### Summary Stats

sumstats = function(x){ 
  y=na.omit(x[x>0])
  bstrap <- c()
  for (i in 1:1000){
    bstrap <- c(bstrap, mean(sample(y,(length(y)),replace=T), na.rm = T))}
  c(
    Mean=mean(y), 
    Median=median(y),
    SD=sd(y), 
    Var=var(y),
    CV=sd(y)/mean(y),
    SE=sd(y)/sqrt(length(y)),
    L95SE=mean(y)-1.96*(sd(y)/sqrt(length(y))),
    U95SE=mean(y)+1.96*(sd(y)/sqrt(length(y))),
    BSMEAN = mean(bstrap),
    L95BS = quantile(bstrap,.025),
    U95BS= quantile(bstrap,.975))
}
