#     Metapop_HeterogeneityExperiment_betaprop_Viz.R visualizes the output from Metapop_HeterogeneityImportThreshold.m
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

rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggpubr) # to make figures with subplots
library(scales) #to rescale colours
library(reldist) #for gini coeff
library(data.table) #for dealing with large datasets
library(dplyr) #for IDing runs
library(moments) # for skewness and kurtosis
library(ggforce) #to draw hulls, etc
library(ggpmisc) #to add equations for lm fit
library(RColorBrewer) #to make custom palette
library(patchwork)  #for multipanel figures

#pick parameters to compare (leave as-is)
param.list<-c("betaA", "betaI", "gammaA", "gammaI")
param1<-param.list[1]
param2<-param.list[2]
outparam1<-param.list[3]
outparam2<-param.list[4]

path = getwd() #Set to wherever you have the files from NetworkStatistics_forHeterogeneityExperiment folder stored

netsize<-c(100)

for (h in 1:length(netsize))
{
  file.names <- dir(path, pattern =sprintf("WSnetstats_1_%iN_26_0.5_r*",netsize[h]))
  
  for (i in 1:length(file.names))
  {
    stats<-fread(file.names[i],header=TRUE,stringsAsFactors = FALSE);
    splits<-strsplit(file.names[i], "_")[[1]]
    stats.annotated=cbind(length(stats$node), as.numeric(splits[[4]]), as.numeric(splits[[5]]), as.numeric(splits[[7]]),  stats)
    
    if (h==1 & i==1)
    {
      stats.comb.0=stats.annotated
    } else {
      stats.comb.0=rbind(stats.annotated, stats.comb.0)
    }
  }
}
colnames(stats.comb.0)[1:4]<-c("N", "neigh", "rewire", "realization")
stats.comb.tidy <- subset(stats.comb.0, select = -c(nodes))

#clean up from stats
rm("stats.annotated", "stats.comb.0", "stats", "splits", "file.names")

#read in ts data
experiment.tag<-"beta_netexpHETERO"
scenario.tag<-"lowyield"
eperc.tag<-16

samples<-seq(0,1,0.05)

for (i in 1:length(samples)) {
  df.temp0<-fread(sprintf("PPlane_NOurban_%s_%s_WSnet_N100_%s%.2f_%s_export%i_new.csv", param1, param2, experiment.tag, samples[i],scenario.tag,eperc.tag),
                  header=FALSE,stringsAsFactors = FALSE);
  colnames(df.temp0)<-c("N", "neigh", "rewire", "gammaI", "betaI", "gammaA", "betaA", "node", "ts", "yield", "pop", "ag", "food", "realization", "highprop", "import" )
  
  df.temp<-df.temp0[df.temp0$ts==100,]
  
  if (i==1)
  {
    df.comb<-df.temp
  } else {
    df.comb<-rbind(df.comb,df.temp)
  }
}

#clean up from reading in data
rm("df.temp0", "df.temp")

DT.0<-merge(df.comb, stats.comb.tidy, by=c("node","N", "neigh", "rewire", "realization"))

#clean up from combining data
rm("df.comb", "stats.comb.tidy")

#ID rows by the simulation they correspond to
DT.tidy <- data.table(DT.0, key=c("N", "neigh", "rewire", "realization", "highprop"))
DT.tidy[, runID:=.GRP, by=key(DT.tidy)]

#clean up after tidying
rm("DT.0") 

#add in column for nonag land
DT.tidy$nonag<-(9.625323203/DT.tidy$N)-DT.tidy$ag
#add in column for nonag land excluding urban area # new added may 5 2021 #
DT.tidy$natland<-(9.625323203/DT.tidy$N) - DT.tidy$ag - (DT.tidy$pop*0.06)

#calculate gini index
for (i in 1:length(unique(DT.tidy$runID)))
{
  for (j in 1:length(unique(DT.tidy$ts)))
  {
    sub<-DT.tidy[(DT.tidy$runID==unique(DT.tidy$runID)[i] & DT.tidy$ts==unique(DT.tidy$ts)[j]),]
    sub$gini.food<-gini(sub$food/sub$pop, sub$pop)
    sub$gini.land<-gini(sub$nonag/sub$pop, sub$pop)
    sub$gini.import<-gini(sub$food/sub$import, sub$pop)
    
    if (i==1 && j==1)
    {
      DT.new<-sub
    } else
    {
      DT.new<-rbind(DT.new, sub)
    }
  }
}

rm("sub")

### Visualize global behaviour
df.sum<-DT.new %>%
  group_by(realization,ts,highprop,sw.measure.TJHBL,density,diameter,av.pl,av.deg,edges,avg.local.cc,avg.local.transitivity,n.components,runID,gini.food,gini.land,gini.import) %>%
  summarize(import.dep=weighted.mean(import/food, pop), pop = sum(pop), ag=sum(ag), natland=sum(natland), food=sum(food), fpc=weighted.mean(food/pop, pop), nonag=sum(nonag), lpc=weighted.mean(nonag/pop, pop))

### Plot global outcomes as a function of highprop

#global network comparison (with summary stats)

var.sample<-c("pop", "ag", "food", "import.dep", "fpc", "lpc", "gini.food", "gini.land", "natland")
varts.labels<-c(expression("Population (x "*10^{9}*")"), expression("Agricultural land (x"*10^{9}*" hectares)"), expression("Food supply (x"*10^{9}*" tonnes)"), 
                paste('Mean import dependency', sep=''), paste("Mean food per capita (tonnes)", sep=''), paste("Mean non-ag. land per capita (hectares)", sep=''), 
                paste('Gini index - food per capita', sep=''), paste('Gini index - non-ag. land per capita', sep=''),expression("Natural land-states (x"*10^{9}*" hectares)"))

p.globalnetcompare.list<-lapply(1:length(var.sample),function(i)
{
  
  ggplot(data=df.sum, aes_string(x="1-highprop", y=var.sample[i])) +
    stat_summary(geom="errorbar", fun.min=min, fun.max=max,  width=0.025)+
    stat_summary(geom="point", fun=mean)+
    xlab("Proportion of patches with high import demand threshold") +
    ylab(varts.labels[i]) +
    scale_x_continuous(expand=c(0,0)) +
    theme_bw()
})

p.netcompare<-wrap_plots(p.globalnetcompare.list[[1]] + theme(axis.title.x = element_blank()) + guides(colour=FALSE),
                         p.globalnetcompare.list[[2]] + theme(axis.title.x = element_blank())  + ylim(c(5.5,6)) + guides(colour=FALSE),
                         p.globalnetcompare.list[[3]]  + guides(colour=FALSE),
                         p.globalnetcompare.list[[4]],
                         nrow = 2) +
  plot_annotation(tag_levels = 'a')

png(sprintf("HIGHPROPnetcompare_NOurban_%s_%s_WSnet_%s_%s_export%i_new2.png", param1, param2, experiment.tag, scenario.tag, eperc.tag), width=17.6, height=17.6, units="cm", res=500)
print(p.netcompare)
dev.off()

#create custom palette based on the RdBu scheme
rdylbu.pal<-palette(brewer.pal(5,"RdYlBu"))
my.pal<-c(rdylbu.pal[1], rdylbu.pal[5])

###Generate spectrum of node-level results
varnode.sample<-c("pop", "ag", "food", "fpc", "lpc", "import.dep", "natland") 
varnodets.labels<-c(expression("Population (x "*10^{9}*")"), expression("Agricultural land (x"*10^{9}*" hectares)"), expression("Food supply (x"*10^{9}*" tonnes)"), 
                    paste("Food per capita (tonnes)", sep=''), paste("Non-agricultural land per capita (hectares)", sep=''), paste('Import dependency', sep=''),
                    expression("Natural land-states (x"*10^{9}*" hectares)"))

DT.new$fpc<-DT.new$food/DT.new$pop
DT.new$lpc<-DT.new$nonag/DT.new$pop
DT.new$import.dep<-DT.new$import/DT.new$food

df.beta.sum<-DT.new %>%
  group_by(highprop, betaI) %>%
  summarise_at(c("pop", "ag", "food", "fpc", "lpc", "import.dep", "natland"), list(~min(.), ~max(.), ~mean(.)))

p.globalnodecompare.list<-lapply(1:length(varnode.sample),function(i)
{
  
  ggplot(data=df.beta.sum, aes_string(x="1-highprop", y=paste0(varnode.sample[i], "_mean"), colour="as.factor(betaI)", group="1-highprop")) +
    geom_errorbar(aes_string(ymin=paste0(varnode.sample[i], "_min"), ymax=paste0(varnode.sample[i], "_max")), width=.025) +
    geom_point() +
    xlab("Proportion of patches with high import demand threshold") +
    scale_colour_manual("Patch import demand threshold", values=my.pal, labels=c(expression("high ("*beta^{I}*"=0.5)"),expression("low ("*beta^{I}*"=1)"))) +
    scale_x_continuous(expand=c(0,0)) +
    ylab(varnodets.labels[i]) +
    theme_bw()  +
    guides(colour = guide_legend(override.aes = list(alpha=1)))
  
})

p.nodecompare<-wrap_plots(p.globalnodecompare.list[[1]] + theme(axis.title.x = element_blank()) + guides(colour=FALSE),
                          p.globalnodecompare.list[[3]] + theme(axis.title.x = element_blank()) + guides(colour=FALSE),
                          p.globalnodecompare.list[[6]], nrow = 3) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')  & theme(legend.position = 'bottom', plot.margin = margin(2, 4, 2, 0, "pt"))

png(sprintf("HIGHPROPnodecompare_NOurban_%s_%s_WSnet_%s_%s_export%i_new2.png", param1, param2, experiment.tag, scenario.tag, eperc.tag), width=17.6, height=23.4, units="cm", res=500)
print(p.nodecompare)
dev.off()

#################

p.fin1<-wrap_plots(p.globalnetcompare.list[[7]] + theme(axis.title.x = element_blank()) + guides(colour=FALSE) + ggtitle("Global outcomes") ,
                  p.globalnetcompare.list[[8]] + theme(axis.title.x = element_blank()),
                  p.globalnetcompare.list[[2]] + theme(axis.title.x = element_blank()),
                  p.globalnetcompare.list[[9]] + theme(axis.title.x = element_blank()),
                  nrow = 2) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')  & theme(legend.position = 'bottom', plot.margin = margin(2, 4, 2, 1, "pt"))

p.fin2<-wrap_plots(p.globalnodecompare.list[[4]] + theme(axis.title.x = element_blank()) + guides(colour=FALSE)  + ggtitle("Patch-level outcomes") , 
                   p.globalnodecompare.list[[5]] + theme(axis.title.x = element_blank()), 
                   p.globalnodecompare.list[[2]], 
                   p.globalnodecompare.list[[7]],
                   nrow = 2) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')  & theme(legend.position = 'bottom', plot.margin = margin(2, 4, 2, 1, "pt"))

p.fin<-wrap_plots(p.fin1,p.fin2,nrow = 2) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect', heights = c(1, 1))  & theme(legend.position = 'bottom', plot.margin = margin(0, 6, 0, 6, "pt"))


ggsave(sprintf("MainBetaVarExp_NOurban_%s_%s_WSnet_%s_%s_export%i_new2.eps", param1, param2, experiment.tag, scenario.tag, eperc.tag), p.fin, width=24, height=32, units="cm", dpi=800)
