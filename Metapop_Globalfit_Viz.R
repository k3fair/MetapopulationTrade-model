#     Metapop_Globalfit_Viz.R visualizes the output from Metapop_Globalfit.m
#     Copyright (C) 2021 Kathyrn R Fair
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.


library(ggplot2) #for plotting
library(ggpubr) #for creating a multipanel figure
library(RColorBrewer) #for ColorBrewer palettes

### Select yield scenario (low or high)
yield<-"low"

#Read in files with timseries data (generated using Metapop_Globalfit.m)
df.low <-read.csv(sprintf("GlobalTrajectoryMODEL_%syield.csv", yield), header=FALSE, stringsAsFactors = FALSE)
df.low$tag<-"Model"

df.fao<-read.csv("GlobalTrajectoryFAO.csv", header=FALSE, stringsAsFactors = FALSE)
df.fao$tag<-"FAO"

df.tot<-rbind(df.fao, df.low)
colnames(df.tot)<-c("Year", "Population", "Ag", "Food", "Yield", "Source")

#Create custom palette based on the RdYlBu scheme
rdylbu.pal<-palette(brewer.pal(5,"RdYlBu"))
my.pal<-c(rdylbu.pal[1], rdylbu.pal[5])

#Generate plots
p.pop<-ggplot(df.tot, aes(x=Year, y=Population, group=Source, colour=Source)) + 
  geom_line() +
  scale_colour_manual(values=my.pal) +
  ylab(expression("Population (x "*10^{9}*")")) +
  theme_bw() +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

p.ag<-ggplot(df.tot, aes_string(x="Year", y="Ag", group="Source", colour="Source")) + 
  geom_line() + 
  scale_colour_manual(values=my.pal) +
  ylab(expression("Agricultural land (x"*10^{9}*" hectares)")) +
  theme_bw() +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

p.food<-ggplot(df.tot, aes_string(x="Year", y="Food", group="Source", colour="Source")) + 
  geom_line() + 
  scale_colour_manual(values=my.pal) +
  ylab(expression("Food supply (x"*10^{9}*" tonnes)")) +
  theme_bw() +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

p.yield<-ggplot(df.tot, aes_string(x="Year", y="Yield", group="Source", colour="Source")) + 
  geom_line() + 
  scale_colour_manual(values=my.pal) +
  theme_bw() +
  ylab("Yield (tonnes per hectare)") +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

p.total<-ggarrange(p.pop + theme(axis.title.x = element_blank()), p.ag + theme(axis.title.x = element_blank()), p.food, p.yield,
                     labels = c("a", "b", "c", "d"), ncol=2, nrow=2, common.legend =TRUE, legend="right")

png(sprintf("GlobalTrajectoryMODELplusFAO_%syield.png", yield), width=6.5, height=5.75, units="in", res=500)
print(p.total)
dev.off()

