
library(multimode)
library(tidyverse)
library(plotrix)
library(readxl)
library(ggplot2)


ARD_V6_1_ <- read_excel("ARD_V6 (1).xls")
data <- ARD_V6_1_

#Whole Sample
  #DIP-Test whole sample
dip.test(data$ifhpol)

  #Silverman Tests Whole Sample
modetest(data$ifhpol, mod0 = 1, method = "SI")
modetest(data$ifhpol, mod0 = 2, method = "SI")
modetest(data$ifhpol, mod0 = 3, method = "SI")
modetest(data$ifhpol, mod0 = 4, method = "SI")

  #Mode location whole sample
locmodes(data$ifhpol, mod0 = 1)
locmodes(data$ifhpol, mod0 = 2)
locmodes(data$ifhpol, mod0 = 3)
    #Mod0 = 3 is just here for posterity

#Start and end of period
dm1972 <- filter(data, year==1972)
dm2014 <- filter(data, year==2014)

#Dip Test for Start and End of Period
  #1972
dip.test(dm1972$ifhpol)

  #2014
dip.test(dm2014$ifhpol)

#Silverman Tests for Start and End of period
  #1972
modetest(dm1972$ifhpol, mod0 = 1, method = "SI")
modetest(dm1972$ifhpol, mod0 = 2, method = "SI")
modetest(dm1972$ifhpol, mod0 = 3, method = "SI")
modetest(dm1972$ifhpol, mod0 = 4, method = "SI")

  #2014
modetest(dm2014$ifhpol, mod0 = 1, method = "SI")
modetest(dm2014$ifhpol, mod0 = 2, method = "SI")
modetest(dm2014$ifhpol, mod0 = 3, method = "SI")
modetest(dm2014$ifhpol, mod0 = 4, method = "SI")

  #Density functions for start and end
kd72u <- density(dm1972$ifhpol, na.rm = TRUE)
kd14u <- density(dm2014$ifhpol, na.rm = TRUE)

#Kernel Density Plots for Start and end of Period
plot(kd72u, main = "Kernel Density Plot 1972")
plot(kd14u, main = "Kernel Density Plot 2014")

#Combined Plot
startend <- subset(data, year == 1972 | year == 2014)


# create a kernel density plot with separate curves for each year
ggplot(startend, aes(x = ifhpol, fill = factor(year), group = factor(year))) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("red", "blue"), name = "Year") +
  labs(x = "ifhpol", y = "Density") +
  theme_bw() +
  theme(legend.position = "top",
        legend.box.spacing = unit(0.2, "cm"),
        legend.margin = margin(l = 0, r = 0, t = 0, b = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(1.2, "lines"))

#Graphs on a decade-by-decade basis (2012 omitted due to (temporal) proximity to 2014)
  #1972
plot(kd72u, , main = "Density in 1972")

  #1982
dm1982 <- filter(data, year==1982)
kd82u <- density(dm1982$ifhpol, na.rm = TRUE)
plot(kd82u, main = "Density in 1982")

  #1992
dm1992 <- filter(data, year==1992)
kd92u <- density(dm1992$ifhpol, na.rm = TRUE)
plot(kd92u, main = "Density in 1992")

  #2002
dm2002 <- filter(data, year==2002)
kd02u <- density(dm2002$ifhpol, na.rm = TRUE)
plot(kd02u, main = "Density in 2002")

  #2014
plot(kd14u, main = "Density in 2014")

#Potential years of special interest (large movements)
  #1988
dm1988 <- filter(data, year == 1988)
kd88u <- density(dm1988$ifhpol, na.rm = TRUE)
plot(kd88u)

  #1989
dm1989 <- filter(data, year == 1989)
kd89u <- density(dm1989$ifhpol, na.rm = TRUE)
plot(kd89u)

  #1990
dm1990 <- filter(data, year == 1990)
kd90u <- density(dm1990$ifhpol, na.rm = TRUE)
plot(kd90u)

  #1991
dm1991 <- filter(data, year == 1991)
kd91u <- density(dm1991$ifhpol, na.rm = TRUE)
plot(kd91u)

#Note to self: Inconsistent y-limits on plots, will specify if/when necessary
##############################################################################
#Regime Survival
  #Filtering for the different regime types
    #Should find a more efficient way to do this
ARD_V6 <- data
dem <- dplyr::filter(ARD_V6, ARD_V6$regimeny == 'democracy')
dem <- dem[!duplicated(dem[,c('totdurny2','country')]),]
summary(dem$totdurny2)

limp <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='limited multiparty')
limp <- limp[!duplicated(limp[,c('totdurny2','country')]),]
summary(limp$totdurny2)

nop <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='no-party')
nop <- nop[!duplicated(nop[,c('totdurny2','country')]),]
summary(nop$totdurny2)

mil <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='military trad')
mil <- mil[!duplicated(mil[,c('totdurny2','country')]),]
summary(mil$totdurny2)

milnp <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='military no-party')
milnp <- milnp[!duplicated(milnp[,c('totdurny2','country')]),]
summary(milnp$totdurny2)

milmp <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='military multiparty')
milmp <- milmp[!duplicated(milmp[,c('totdurny2','country')]),]
summary(milmp$totdurny2)

milop <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='military one-party')
milop <- milop[!duplicated(milop[,c('totdurny2','country')]),]
summary(milop$totdurny2)

op <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='one-party')
op <- op[!duplicated(op[,c('totdurny2','country')]),]
summary(op$totdurny2)

oth <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='other')
oth <- oth[!duplicated(oth[,c('totdurny2','country')]),]
summary(oth$totdurny2)

opmon <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='one-party monarchy')
opmon <- opmon[!duplicated(opmon[,c('totdurny2','country')]),]
summary(opmon$totdurny2)

mon <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='monarchy')
mon <- mon[!duplicated(mon[,c('totdurny2','country')]),]
summary(mon$totdurny2)

theo <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='theocracy')
theo <- theo[!duplicated(theo[,c('totdurny2','country')]),]
summary(theo$totdurny2)

nopmon <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='no-party monarchy')
nopmon <- nopmon[!duplicated(nopmon[,c('totdurny2','country')]),]
summary(nopmon$totdurny2)

mpmon <- dplyr:: filter(ARD_V6, ARD_V6$regimeny=='multiparty monarchy')
mpmon <- mpmon[!duplicated(mpmon[,c('totdurny2','country')]),]
summary(mpmon$totdurny2)



  #Average duration for different Regime types
ggplot(NULL, aes(x=regimeny)) +
  geom_bar(data=limp, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=oth, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=mil, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=milmp, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=milnp, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=milop, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=mon, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=nopmon, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=mpmon, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=op, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  geom_bar(data=dem, aes(y=totdurny2), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  labs(title='Average regime duration', y='years', x = 'regime') +
  coord_flip()

#Average Regime duration by score
df6 <- dplyr::filter(ARD_V6, ARD_V6$totdurny2 != 'NA',
                     ARD_V6$ifhpol != 'NA')
df6$cat <- plyr::round_any(df6$ifhpol, 0.5)

ggplot(NULL, aes(y=totdurny2 , x=cat), group=year) +
  geom_bar(df6,mapping = aes(y=totdurny2 , x=cat), position = "dodge", stat = "summary", fun.y = "mean", fill="#E69F00") +
  labs(title='Average regime duration by freedom score', y='years', x = 'score') + 
  scale_x_continuous(breaks = seq(0,10, by = 1))

  #Setup for regression, removing non-relevant regime types
df7 <- dplyr::filter(ARD_V6, ARD_V6$regimeny != 'civil war')
df7 <- dplyr::filter(df7, df7$regimeny != 'rebel regime')
df7 <- dplyr::filter(df7, df7$regimeny != 'occupation')
df7 <- dplyr::filter(df7, df7$regimeny != 'transitional')

ggplot(df7, aes(x=ifhpol, y=totdurny2))+
  geom_point(data=df7, aes(x=round(ifhpol,digits=1), y=totdurny2, colour= regimeny)) +
  scale_color_manual(values=c('deepskyblue3', 'blueviolet','green','darkseagreen4','chartreuse','darkolivegreen','darkorange','darkgoldenrod1','brown1','brown4', 'azure4', 'darksalmon', "black", "darkblue", "purple")) +
  labs(title='Distribution', y='years', x = 'score') +
  geom_smooth(
    method = "lm",
    formula = y~I(x^3)+I(x^2)+x,
    se = TRUE,
    fullrange = TRUE) 

#Linear "vs" non-linear
linreg <- lm(totdurny2 ~ ifhpol, data = data)
summary(linreg)
AIC(linreg)

nonreg <- lm(totdurny2 ~ I(ifhpol^3)+I(ifhpol^2)+ifhpol, data = data)
summary(nonreg)
AIC(nonreg)

#Higher r2 for non-linear, implies better fit. Lower AIC-"score", implies that the 
  #higher r2 is not (purely) a result of additional dependent variables

#Reducing scope of dataset for Sirnes
data2 <- data %>% select(country, year, ifhpol, regimeny, totdurny2)
saveRDS(data2, file = "data.Rda")

write.csv(data2, "C:\\Users\\CEPTER\\Clean Data\\datasimplified.csv")
write.csv(data, "C:\\Users\\CEPTER\\Clean Data\\datafull.csv")
