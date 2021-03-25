############
## GLM Models using glmmTMB package - nbinom2
## Jennifer Moore
## Big Bend Trends Paper
## August 2019
############

#import libraries
library(glmmTMB)
library(sjPlot)
library(waterData)
library(reshape)
library(ggplot2)
library(lubridate)
library(AICcmodavg)
library(DMwR)
library(ggeffects)
library(cowplot)

#import data
d <- read.csv("data/all_epoch_transect_data.csv", header = T)

#turn oyster counts into integers
#b/c we averaged 2 observers not all values are integers
d$count_live<-round(d$count_live, digits=0)

#only keep stations that have been sampled more than once
t <- data.frame(count = rowSums(table(d$station, d$period)))
a <- which(t$count == 1)
b <- row.names(t)[a]
d2 <- subset(d, !(d$station %in% b))
d2$station <- droplevels(d2$station)




################
#basic plots

#hist (distribution) of count data overall
hist(d2$count_live, xlab ="Oyster Counts", main = NULL)
#boxplots per locality
boxplot(count_live ~ locality, data = d2, ylab = "Oyster Counts")
#boxplot per period
boxplot(count_live ~ period, data = d2, ylab = "Oyster Counts", xlab = "Period")
#mean oyster count per period
agg <- aggregate(count_live ~ period, data = d2, FUN = 'mean')
plot(agg$period, agg$count_live, type = 'l', xlab = "Period", ylab = "Mean Oyster Count")

#histogram of oyster counts by density
h <- ggplot(d2, aes(count_live))+
  geom_histogram(aes(x = count_live, y = ..density..),bins=55, fill = 'grey', col = 'black')+
  theme_classic()+
  ylim (0,0.0025) +
  ylab("Density") +
  xlab ("Oyster Counts")

#plot of live counts by period
ggplot(d2, aes(period, count_live, colour=site))+
  geom_point(size=3)+
  labs(title = "Lone Cabbage oyster count by site and period")+
  ylab("Live oyster count") +
  xlab ("Period")+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))

# Count Histogram, for assessing GLM family fit to count data - fit neg bin model
theta=c(mu=1,k=0.2)
nb_LL=function(theta)
{
  -sum(dnbinom(d$count_live,mu=theta[1],size=theta[2],log=TRUE))
}
fit_nb=optim(theta,nb_LL)
out = dnbinom(seq(0,5350,25),mu=fit_nb$par[1],size=fit_nb$par[2])
h + geom_line(aes(x = seq(0,5350,25), y = out), size = 1.5, col = 'red')



#######################
#basic GLMM w/o covariates

glmm1 <- glmmTMB(count_live ~ period + offset(log(tran_length)), data = d2, family = 'nbinom2') 
glmm2 <- glmmTMB(count_live ~ period + locality + offset(log(tran_length)), data = d2, family = 'nbinom2') 
glmm3 <- glmmTMB(count_live ~ period * locality + offset(log(tran_length)), data = d2, family = 'nbinom2')
glmm4 <- glmmTMB(count_live ~ period + site + offset(log(tran_length)), data = d2, family = 'nbinom2') 
glmm5 <- glmmTMB(count_live ~ period * site + offset(log(tran_length)), data = d2, family = 'nbinom2') 
glmm6 <- glmmTMB(count_live ~ period + locality + site + offset(log(tran_length)), data = d2, family = 'nbinom2')
glmm7 <- glmmTMB(count_live ~ period + locality * site + offset(log(tran_length)), data = d2, family = 'nbinom2') 
glmm8 <- glmmTMB(count_live ~ period * locality + site + offset(log(tran_length)), data = d2, family = 'nbinom2') 
glmm9 <- glmmTMB(count_live ~ period * locality * site + offset(log(tran_length)), data = d2, family = 'nbinom2') 
glmm10 <- glmmTMB(count_live ~ period * site + locality + offset(log(tran_length)), data = d2, family = 'nbinom2') 


cand.set = list(glmm1, glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, glmm8, glmm9, glmm10)
modnames = c("period", "period + locality", "period * locality", "period + site", "period * site", "period + locality + site", "period + locality * site", "period * locality + site", "period * locality * site", "period * site + locality")
aictab(cand.set, modnames, second.ord = FALSE) #model selection table with AIC

#top  model is glmm8

#now create a 4 panel plot for all the localities
#first 4 panel plot just shows raw data
p1 = ggplot(d2[d2$locality == "CK",], aes(period, count_live, colour=site))+
  geom_point(size=3)+
  labs(title = "CK")+
  ylab("Live oyster count") +
  xlab ("Period")+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))+
  ylim(0,5000)
p2 = ggplot(d2[d2$locality == "CR",], aes(period, count_live, colour=site))+
  geom_point(size=3)+
  labs(title = "CR")+
  ylab("Live oyster count") +
  xlab ("Period")+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))+
  ylim(0,5000)
p3 = ggplot(d2[d2$locality == "HB",], aes(period, count_live, colour=site))+
  geom_point(size=3)+
  labs(title = "HB")+
  ylab("Live oyster count") +
  xlab ("Period")+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))+
  ylim(0,5000)
p4 = ggplot(d2[d2$locality == "LC",], aes(period, count_live, colour=site))+
  geom_point(size=3)+
  labs(title = "LC")+
  ylab("Live oyster count") +
  xlab ("Period")+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))+
  ylim(0,5000)
#arrange all 4 plots on one page w/o legend
p <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2
)
#add legend back in
# extract the legend from one of the plots
legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(p, legend, rel_widths = c(2.5, .4))


##################
#second 4 panel plot shows raw data
#and predicted model line
testCK = ggpredict(glmm8, terms = c("period", "site", "locality[CK]", "tran_length[33.23]"), type = c("fe"), data = d2)
newCK = data.frame(x = testCK$x, predicted = testCK$predicted, site = rep(c("I", "N", "O"), each = 1, times = 10), conf.low = testCK$conf.low, conf.high = testCK$conf.high)
testCR = ggpredict(glmm8, terms = c("period", "site", "locality[CR]", "tran_length[33.23]"), type = c("fe"), data = d2)
newCR = data.frame(x = testCR$x, predicted = testCR$predicted, site = rep(c("I", "N", "O"), each = 1, times = 10), conf.low = testCR$conf.low, conf.high = testCR$conf.high)
testHB = ggpredict(glmm8, terms = c("period", "site", "locality[HB]", "tran_length[33.23]"), type = c("fe"), data = d2)
newHB = data.frame(x = testHB$x, predicted = testHB$predicted, site = rep(c("I", "N", "O"), each = 1, times = 10), conf.low = testHB$conf.low, conf.high = testHB$conf.high)
testLC = ggpredict(glmm8, terms = c("period", "site", "locality[LC]", "tran_length[33.23]"), type = c("fe"), data = d2)
newLC = data.frame(x = testLC$x, predicted = testLC$predicted, site = rep(c("I", "N", "O"), each = 1, times = 10), conf.low = testLC$conf.low, conf.high = testLC$conf.high)

pm1 = ggplot(newCK, aes(x, predicted, group=site, col=site))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = site), alpha = .5) +
  geom_point(data = d2[d2$locality == "CK",], mapping = aes(period, count_live), size = 2) +
  labs(title = "CK")+
  ylab("Live oyster count") +
  xlab ("Period")+
  theme(panel.grid = element_blank(), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face='bold'),
        plot.title=element_text(size=16,face='bold'))+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))+
  ylim(0,12000)
pm2 = ggplot(newCR, aes(x, predicted, group=site, col=site))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill=site), alpha = .5) +
  geom_point(data = d2[d2$locality == "CR",], mapping = aes(period, count_live),size=2) +
  labs(title = "CR")+
  ylab("Live oyster count") +
  xlab ("Period")+
  theme(panel.grid = element_blank(), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face='bold'),
        plot.title=element_text(size=16,face='bold'))+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))+
  ylim(0,12000)
pm3 = ggplot(newHB, aes(x, predicted, group = site, col=site))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = site), alpha = .5) +
  geom_point(data = d2[d2$locality == "HB",], mapping = aes(period, count_live),size=2) +
  labs(title = "HB")+
  ylab("Live oyster count") +
  xlab ("Period")+
  theme(panel.grid = element_blank(), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face='bold'),
        plot.title=element_text(size=16,face='bold'))+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))+
  ylim(0,12000)
pm4 = ggplot(newLC, aes(x, predicted, group=site, col=site))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill=site), alpha = .5) +
  geom_point(data = d2[d2$locality == "LC",], mapping = aes(period, count_live),size=2) +
  labs(title = "LC")+
  ylab("Live oyster count") +
  xlab ("Period")+
  theme(panel.grid = element_blank(), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face='bold'),
        plot.title=element_text(size=16,face='bold'))+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks=seq(1,17,1))+
  ylim(0,12000)
#arrange all 4 plots on one page w/o legend
pm <- plot_grid(
  pm1 + theme(legend.position="none"),
  pm2 + theme(legend.position="none"),
  pm3 + theme(legend.position="none"),
  pm4 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2
)
#add legend back in
# extract the legend from one of the plots
legend <- get_legend(
  pm1 + theme(legend.box.margin = margin(0, 0, 0, 12),
              legend.title = element_text(size=14),
              legend.text = element_text(size=12))
)
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(pm, legend, rel_widths = c(2.5,.4))








#check dispersion
E4 <- resid(glmm8, type = 'pearson')
N <- nrow(d2)
p <- length(fixef(glmm10))+1
dispersion <- sum(E4^2) / (N-p)
dispersion # 1.25



##########################################
#GLMM w/ covariates
#look at trend in covariate and then run model

#pull in covariates
#river discharge data - suwannee
#station to analyze
station = '02323500'   
#get site name to use in plot titles and such
stinfo  = siteInfo(station)
#read entire time series
dis = importDVs(staid=station,code='00060',stat='00003', sdate= "1941-01-01") 
dis$year = as.numeric(strftime(dis$dates,format="%Y"))
#Naming columns, using the Diver sensors, collects date, pressure, temp, conductivity
colnames(dis) <- c("StaID", "Discharge", "Date", "QualCode", "Year")

#pull in the wacasassa data
disW <- read.csv("data/WaccasassaRiverDischarge.csv", header = T)
disW$Date <- as.Date(disW$Date, "%m/%d/%Y")
disW$year    = as.numeric(strftime(disW$Date,format="%Y"))
disW <- subset(disW, select = c("Station.ID", "Value", "Date", "X.1", "year"))
colnames(disW) <- c("StaID", "Discharge", "Date", "QualCode", "Year")


dis_mean_year <- aggregate(Discharge ~ Year, dis, mean)
dis_mean_year$YearLag1 <- dis_mean_year$Year + 1
dis_mean_year$YearLag2 <- dis_mean_year$Year + 2
dis_mean_yearW <- aggregate(Discharge ~ Year, disW, mean)
dis_mean_yearW$YearLag1 <- dis_mean_yearW$Year + 1
dis_mean_yearW$YearLag2 <- dis_mean_yearW$Year + 2
dis_total_year <- aggregate(Discharge ~ Year, dis, sum)
dis_total_year$YearLag1 <- dis_total_year$Year + 1
dis_total_year$YearLag2 <- dis_total_year$Year + 2
dis_total_yearW <- aggregate(Discharge ~ Year, disW, sum)
dis_total_yearW$YearLag1 <- dis_total_yearW$Year + 1
dis_total_yearW$YearLag2 <- dis_total_yearW$Year + 2
colnames(dis_mean_yearW) <- c("YearW", "DischargeW", "YearLag1W", "YearLag2W")
colnames(dis_total_yearW) <- c("YearW", "DischargeW", "YearLag1W", "YearLag2W")

#add covariate data to original data set - first line adds discharge data for suwannee
#second line adds discharge data for wacasassa
d3 <- merge(d2, dis_mean_year[,c("Year", "Discharge")], by.x = "year", by.y = "Year", all.x = TRUE)
names(d3)[length(names(d3))] <- "AD"
d3 <- merge(d3, dis_mean_yearW[,c("YearW", "DischargeW")], by.x = "year", by.y = "YearW", all.x = TRUE)
names(d3)[length(names(d3))] <- "AD_W"

d4 <- merge(d3, dis_mean_year[,c("YearLag1", "Discharge")], by.x = "year", by.y = "YearLag1", all.x = TRUE)
names(d4)[length(names(d4))] <- "ADLag1"
d4 <- merge(d4, dis_mean_yearW[,c("YearLag1W", "DischargeW")], by.x = "year", by.y = "YearLag1W", all.x = TRUE)
names(d4)[length(names(d4))] <- "ADLag1_W"

d5 <- merge(d4, dis_mean_year[,c("YearLag2", "Discharge")], by.x = "year", by.y = "YearLag2", all.x = TRUE)
names(d5)[length(names(d5))] <- "ADLag2"
d5 <- merge(d5, dis_mean_yearW[,c("YearLag2W", "DischargeW")], by.x = "year", by.y = "YearLag2W", all.x = TRUE)
names(d5)[length(names(d5))] <- "ADLag2_W"

d6 <- merge(d5, dis_total_year[,c("Year", "Discharge")], by.x = "year", by.y = "Year", all.x = TRUE)
names(d6)[length(names(d6))] <- "TD"
d6 <- merge(d6, dis_total_yearW[,c("YearW", "DischargeW")], by.x = "year", by.y = "YearW", all.x = TRUE)
names(d6)[length(names(d6))] <- "TD_W"

d7 <- merge(d6, dis_total_year[,c("YearLag1", "Discharge")], by.x = "year", by.y = "YearLag1", all.x = TRUE)
names(d7)[length(names(d7))] <- "TDLag1"
d7 <- merge(d7, dis_total_yearW[,c("YearLag1W", "DischargeW")], by.x = "year", by.y = "YearLag1W", all.x = TRUE)
names(d7)[length(names(d7))] <- "TDLag1_W"

d8 <- merge(d7, dis_total_year[,c("YearLag2", "Discharge")], by.x = "year", by.y = "YearLag2", all.x = TRUE)
names(d8)[length(names(d8))] <- "TDLag2"
d8 <- merge(d8, dis_total_yearW[,c("YearLag2W", "DischargeW")], by.x = "year", by.y = "YearLag2W", all.x = TRUE)
names(d8)[length(names(d8))] <- "TDLag2_W"

d8$ADsc <- scale(d8$AD)
d8$ADLag1sc <- scale(d8$ADLag1)
d8$ADLag2sc <- scale(d8$ADLag2)
d8$TDsc <- scale(d8$TD)
d8$TDLag1sc <- scale(d8$TDLag1)
d8$TDLag2sc <- scale(d8$TDLag2)
d8$ADsc_W <- scale(d8$AD_W)
d8$ADLag1sc_W <- scale(d8$ADLag1_W)
d8$ADLag2sc_W <- scale(d8$ADLag2_W)
d8$TDsc_W <- scale(d8$TD_W)
d8$TDLag1sc_W <- scale(d8$TDLag1_W)
d8$TDLag2sc_W <- scale(d8$TDLag2_W)

#oyster landing data
ol <- read.csv("data/suwannee_data.csv", header = T)
landing <- subset(ol, ol$value == "LANDINGS (KG)")
landing$YearLag1 <- landing$Year + 1
landing$YearLag2 <- landing$Year + 2
trips <- subset(ol, ol$value == "TRIPS")
trips$YearLag1 <- trips$Year + 1
trips$YearLag2 <- trips$Year + 2
d9 <- merge(d8, landing[,c("Year", "measurement")], by.x = "year", by.y = "Year", all.x = TRUE)
names(d9)[length(names(d9))] <- "Landings"
d10 <- merge(d9, landing[,c("YearLag1", "measurement")], by.x = "year", by.y = "YearLag1", all.x = TRUE)
names(d10)[length(names(d10))] <- "LandingsLag1"
d11 <- merge(d10, landing[,c("YearLag2", "measurement")], by.x = "year", by.y = "YearLag2", all.x = TRUE)
names(d11)[length(names(d11))] <- "LandingsLag2"
d12 <- merge(d11, trips[,c("Year", "measurement")], by.x = "year", by.y = "Year", all.x = TRUE)
names(d12)[length(names(d12))] <- "Trips"
d13 <- merge(d12, trips[,c("YearLag1", "measurement")], by.x = "year", by.y = "YearLag1", all.x = TRUE)
names(d13)[length(names(d13))] <- "TripsLag1"
d14 <- merge(d13, trips[,c("YearLag2", "measurement")], by.x = "year", by.y = "YearLag2", all.x = TRUE)
names(d14)[length(names(d14))] <- "TripsLag2"

d14$Landingssc <- scale(d14$Landings)
d14$LandingsLag1sc <- scale(d14$LandingsLag1)
d14$LandingsLag2sc <- scale(d14$LandingsLag2)
d14$Tripssc <- scale(d14$Trips)
d14$TripsLag1sc <- scale(d14$TripsLag1)
d14$TripsLag2sc <- scale(d14$TripsLag2)

#fishing (yes/no) based on strata column
d14$fish <- substr(as.character(d14$strata), 1, 1)



#############
## Covariate models
############
#start with top model glmm8
#glmm8: count_live ~ period * locality + site + offset(log(tran_length)) 


#river discharge year of sampling (yearly mean)
glmm11 <- glmmTMB(count_live ~ period * locality + site + ADsc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#river discharge 1 year lag (yearly mean)
glmm12 <- glmmTMB(count_live ~ period * locality + site + ADLag1sc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#river discharge 2 year lag (yearly mean)
glmm13 <- glmmTMB(count_live ~ period * locality + site + ADLag2sc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#river discharge year of sampling (yearly total)
glmm14 <- glmmTMB(count_live ~ period * locality + site + TDsc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#river discharge 1 year lag (yearly total)
glmm15 <- glmmTMB(count_live ~ period * locality + site + TDLag1sc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#river discharge 2 year lag (yearly total)
glmm16 <- glmmTMB(count_live ~ period * locality + site + TDLag2sc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#oyster meat weight (landings) year of sampling
glmm17 <- glmmTMB(count_live ~ period * locality + site + Landingssc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#oyster meat weight 1 year lag
glmm18 <- glmmTMB(count_live ~ period * locality + site + LandingsLag1sc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#oyster meat weight 2 year lag
glmm19 <- glmmTMB(count_live ~ period * locality + site + LandingsLag2sc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#oyster fishing trips year of sampling
glmm20 <- glmmTMB(count_live ~ period * locality + site + Tripssc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#oyster fishing trips 1 year lag
glmm21 <- glmmTMB(count_live ~ period * locality + site + TripsLag1sc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#oyster fishing trips 2 year lag
glmm22 <- glmmTMB(count_live ~ period * locality + site + TripsLag2sc + offset(log(tran_length)), data = d14, family = 'nbinom2') 
#oyster fishing (yes/no) 
glmm23 <- glmmTMB(count_live ~ period * locality + site + fish + offset(log(tran_length)), data = d14, family = 'nbinom2')

cand.set = list(glmm11, glmm12, glmm13, glmm14, glmm15, glmm16, glmm17, glmm18, glmm19, glmm20, glmm21, glmm22, glmm23)
modnames = c("AD", "ADLag1", "ADLag2", "TD", "TDLag1", "TDLag2", "Landings", "LandingsLag1", "LandingsLag2", "Trips", "TripsLag1", "TripsLag2", "Harvest")
aictab(cand.set, modnames, second.ord = FALSE) #model selection table with AIC

#top model glmm12

s <- function(x) {(x*sd(d14$ADLag1))+mean(d14$ADLag1)}

plot_model(glmm12, type = 'pred', terms = c('ADLag1sc[s]'))
plot_model(glmm12, type = 'pred', terms = c('period', 'site', 'ADLag1sc'))
plot_model(glmm15, type = 'pred', terms = 'TDLag1sc')

#plot model 12 for ADLag1 unscaled 
library(effects)
library(ggeffects)
library(ggthemes)

# I redid this cause I was having an issue with ggpredict not finding the name
# "ADLag1sc" (I think this is b/c its a funky format in d14, see str(d14))
new.dat = data.frame(count_live = d14$count_live,
                     period = d14$period,
                     site = d14$site,
                     locality = d14$locality,
                     ADLag1sc = d14$ADLag1sc,
                     tran_length = log(d14$tran_length))

glmm12_new <- glmmTMB(count_live ~ period * locality + site + ADLag1sc + offset(tran_length), data = new.dat, family = 'nbinom1') 

# see https://urldefense.proofpoint.com/v2/url?u=https-3A__strengejacke.github.io_ggeffects_reference_ggpredict.html&d=DwIGAg&c=sJ6xIWYx-zLMB3EPkvcnVg&r=0N2A3Co00ZnRVDd_w9mp4Q&m=vg5HvjIe8pL7XypTtg9b5_8orUTjE3VJgISgC4C4z5U&s=VSs2A0265UFXr_d08pSBSlbkRF1aSaz1XoRh84RKhBE&e= 
# see the text talking about how these are marginal effects (i.e., at the mean
# value of everything else in the model, or a reference level, when you only
# specify a subset of terms, I think thats what is showing up in the output of
# 'test', as 'Adjusted for:')
test = ggpredict(glmm12_new, terms = c("ADLag1sc", "tran_length[exp=1]"), type = c("fe"), data = new.dat)
test2 = ggpredict(glmm12_new, terms = c("ADLag1sc", "tran_length[exp=0]"), type = c("fe"), data = new.dat)

# predictions for the scaled x values
p = ggplot(test, aes(x = x, y = predicted)) +
  geom_line()

p

# see link: https://urldefense.proofpoint.com/v2/url?u=https-3A__stats.stackexchange.com_questions_209784_rescale-2Dpredictions-2Dof-2Dregression-2Dmodel-2Dfitted-2Don-2Dscaled-2Dpredictors&d=DwIGAg&c=sJ6xIWYx-zLMB3EPkvcnVg&r=0N2A3Co00ZnRVDd_w9mp4Q&m=vg5HvjIe8pL7XypTtg9b5_8orUTjE3VJgISgC4C4z5U&s=c8W0XZEIH7NyBS6Ekhkv1RRIX1Zdtm6BBTvDW5Ye79s&e= 
scaled_lag <- scale(d8$ADLag1)
scale_list = list(scale = attr(scaled_lag, "scaled:scale"),
                  center = attr(scaled_lag, "scaled:center"))

test$x_unscaled = test$x * scale_list$scale + scale_list$center
test2$x_unscaled = test$x * scale_list$scale + scale_list$center

# predictions back on the natural scale of x 
p2 = ggplot(test, aes(x = x_unscaled, y = predicted)) +
  geom_line()

p2



# i think this means, more water, more oysters?!
p3 = ggplot(test, aes(x = x_unscaled, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "salmon", alpha = .5) +
  xlim(3500,12500) +
  labs(y = 'Oyster Counts (per 1 m)', x = 'Mean Annual Daily Discharge (cf/s)') +
  theme_classic()


p3


##########################
## Simulations
##########################

#simulate data (nsim = number of times) from the live count data (sample command)
#data simulated for each period, site, locality

period <- NA
pval <- NA
new <- data.frame()
sim <- list()
nsim = 500

#loop over each unique period, site and locality
#pull out the data from the full data set for each combination of period,site,locality
#resample from these data with replacement

#add the simulated counts back to the original dataset
#sim[[1]] -> simulation 1 counts
#sim[[2]] -> simulation 2 counts etc.
for(a in 1:nsim){
  for(i in unique(d14$period)){
    for(j in unique(d14$site)){
      for(k in unique(d14$locality)){
        temp <- subset(d14, d14$period == i & d14$site == j & d14$locality == k)
        temp$count2 <- sample(temp$count_live, nrow(temp), replace = T)
        new <- rbind(new, temp)
      }
    }
  }
  endrow <- a*215
  startrow <- endrow-214
  sim[[a]] <- new[startrow:endrow,]
}

saveRDS(mod,"simdata.rds") #save the output of the simulation list as rds
simdata<-readRDS("simdata.rds") #this is how your get it back


#use the simulated counts to rerun the model
#pull out the beta coefficient for period and the associated p value
#create a plot of the trajectories for counts over periods

mod <- list()
#test <- ggpredict(glmm10, terms = c('period'), type = c('fe'), data = d14)
plot(test$x, test$predicted, type = 'l', xlab = "Period", ylab = "Live Oyster Counts", col = 'blue', ylim = c(0,500))

for(i in 1:nsim){
  m <- glmmTMB(count2 ~ period * site + locality + offset(log(tran_length)), family = 'nbinom1', data = sim[[i]])
  period[i] <- as.vector(getME(m, 'beta'))[2]
  pval[i] <- summary(m)[[6]]$cond[2,4]
  mod[[i]] <- ggpredict(m, terms = c('period'), type = c('fe'), data = sim[[i]])
  lines(mod[[i]]$x, mod[[i]]$predicted)
  print(i)
}

lines(test$x, test$predicted, type = 'l',lwd=13, col = 'white')


#create ggplot version of the graph - draw a line for each simulation
p = ggplot(test, aes(x = x, y = predicted)) +
  geom_line() +  
  labs(y = 'Oyster Counts', x = 'Period') +
  theme_classic()
i = 1
while(i <= nsim){
  df <- mod[[i]]
  p <- p + geom_line(data = df, aes(x, predicted))
  i <- i+1
}
p <- p + geom_line(data = test, aes(x,predicted), size = 1.2, col = "blue") #add observed data line on top in blue
p #display plot


#saveRDS(mod,"simoutput.rds") #save the output of the simulation list as rds
#simout<-readRDS("simoutput.rds") #this is how your get it back
save.image("simdata.RData") #save output as rData file
load("simdata.RData") #reload data

#ok need to know how many beta are negative

sum(period < 0)
#all 1000

#how many pvalue < 0.05
sum(pval < 0.05)
#952

hist(pval)

mean(pval)

zz<-density(pval)
plot(zz,type = 'l', xlab = "p-value", ylab = "Density", col = 'blue', ylim = c(0,150),main="")

#convert density plot for pvalues into ggplot
den <- data.frame(zz$x, zz$y)#input to ggplot must be a dataframe
ggplot(den, aes(x = zz.x, y = zz.y)) +
  geom_line() +  
  labs(y = 'Density', x = 'p-value') +
  theme_classic()