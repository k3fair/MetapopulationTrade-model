#     Metapop_HeterogeneityExperiment_betaprop_Viz.R visualizes the output from Metapop_HeterogeneityImportThreshold.m or Metapop_HeterogeneityResponsiveness.m 
#     Copyright (C) 2021 Kathyrn R Fair
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

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
library(patchwork) #for multipanel figures

#pick parameters to compare (leave as-is)
param.list<-c("betaA", "betaI", "gammaA", "gammaI")
param1<-param.list[1]
param2<-param.list[2]
outparam1<-param.list[3]
outparam2<-param.list[4]

#Select experiment; hetergeneous responsiveness ("gamma_netexpHETERO") or heterogeneous import demand threshold ("beta_netexpHETERO")
experiment.tag<-"beta_netexpHETERO"

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
scenario.tag<-"lowyield"
eperc.tag<-16

if (experiment.tag=="gamma_netexpHETERO")
{
  node.tag<-"as.factor(gammaA)"
  node.labels<-c(expression("low ("*gamma*"=2.5)"),expression("high ("*gamma*"=7.5)"))
  node.title<-"Patch demand responsiveness"
  net.title<-"Patch demand \nresponsiveness distribution"
  group.tag="interaction(highprop,gammaA)"
  net.labels<- c(expression("100% low ("*gamma*"=2.5)"),expression("50% low ("*gamma*"=2.5), 50% high ("*gamma*"=7.5)"),expression("100% high ("*gamma*"=7.5)"))
  combo.title<-paste("Patch demand responsiveness \n(responsiveness distribution)", sep='')
  combo.labels<-c(paste("Low \n(100% low)", sep=''), paste("Low \n(50% low, 50% high)", sep=''), paste("High \n(50% low, 50% high)", sep=''), paste("High \n(100% high)", sep=''))
} else {
  node.tag<-"as.factor(betaI)"
  node.labels<-c(expression("High ("*beta^{I}*"=0.5)"),expression("Low ("*beta^{I}*"=1)"))
  node.title<-"Patch import demand threshold"
  net.title<-paste("Patch import demand \nthreshold distribution", sep='')
  group.tag="interaction(highprop,betaI)"
  net.labels<-c(expression("100% high ("*beta^{I}*"=0.5)"),expression("50% high ("*beta^{I}*"=0.5), 50% low ("*beta^{I}*"=1)"), expression("100% low ("*beta^{I}*"=1)"))
  combo.title<-paste("Patch import demand threshold \n(threshold distribution)", sep='')
  combo.labels<-c(paste("High \n(100% High)", sep=''), paste("High \n(50% low, 50% high)", sep=''), 
                  paste("Low \n(50% low, 50% high)", sep=''), paste("Low \n(100% low)", sep=''))
}

df.h0<-fread(sprintf("PPlane_NOurban_%s_%s_WSnet_N100_%s0.00_%s_export%i_new.csv", param1, param2, experiment.tag, scenario.tag,eperc.tag),
                 header=FALSE,stringsAsFactors = FALSE);
colnames(df.h0)<-c("N", "neigh", "rewire", "gammaI", "betaI", "gammaA", "betaA", "node", "ts", "yield", "pop", "ag", "food", "realization", "highprop", "import" )

df.h50<-fread(sprintf("PPlane_NOurban_%s_%s_WSnet_N100_%s0.50_%s_export%i_new.csv", param1, param2, experiment.tag, scenario.tag,eperc.tag),
                       header=FALSE,stringsAsFactors = FALSE);
colnames(df.h50)<-c("N", "neigh", "rewire", "gammaI", "betaI", "gammaA", "betaA", "node", "ts", "yield", "pop", "ag", "food", "realization", "highprop", "import" )

df.h100<-fread(sprintf("PPlane_NOurban_%s_%s_WSnet_N100_%s1.00_%s_export%i_new.csv", param1, param2, experiment.tag, scenario.tag,eperc.tag),
              header=FALSE,stringsAsFactors = FALSE);
colnames(df.h100)<-c("N", "neigh", "rewire", "gammaI", "betaI", "gammaA", "betaA", "node", "ts", "yield", "pop", "ag", "food", "realization", "highprop", "import" )

df.comb<-rbind(df.h0, df.h50, df.h100)

#clean up from reading in data
rm("df.h0", "df.h50", "df.h100")

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
  group_by_at(c(5,10,15,22:31,34:36)) %>%
  summarize(import.dep=weighted.mean(import/food, pop), pop = sum(pop), ag=sum(ag), food=sum(food), fpc=weighted.mean(food/pop, pop), nonag=sum(nonag), lpc=weighted.mean(nonag/pop, pop), natland=sum(natland))


#create custom palette based on the RdBu scheme
rdylbu.pal<-palette(brewer.pal(5,"RdYlBu"))
my.pal<-c(rdylbu.pal[1], "darkorange", rdylbu.pal[5])

#global network comparison (with summary stats)
var.sample<-c("pop", "ag", "food", "import.dep", "fpc", "lpc", "gini.food", "gini.land", "natland")
varts.labels<-c(expression("Population (x "*10^{9}*")"), expression("Agricultural land (x"*10^{9}*" hectares)"), expression("Food supply (x"*10^{9}*" tonnes)"), 
                paste('Mean import dependency', sep=''), paste("Mean food per capita (tonnes)", sep=''), paste("Mean non-ag. land per capita (hectares)", sep=''), 
                paste('Gini index - food per capita', sep=''), paste('Gini index - non-ag. land per capita', sep=''), expression("Natural land-states (x"*10^{9}*" hectares)"))

p.globalnetcompare.list<-lapply(1:length(var.sample),function(i)
{
  
    ggplot(data=df.sum, aes_string(x="ts", y=var.sample[i], colour="as.factor(highprop)", fill="as.factor(highprop)", group="highprop")) +
    geom_line(aes(group=runID), alpha=0.1) +
    stat_summary(geom="line", fun=mean, linetype="solid")+
    scale_colour_manual(net.title, values=my.pal, labels=net.labels) +
    scale_fill_manual(net.title, values=my.pal, labels=net.labels) +
    scale_x_continuous(expand=c(0,0)) +
    xlab("Time-step (t)") +
    ylab(varts.labels[i]) +
    theme_bw()
})


p.ts.netcompare<-wrap_plots(p.globalnetcompare.list[[1]] + theme(axis.title.x = element_blank()) + guides(colour = "none"),
                            p.globalnetcompare.list[[3]] + theme(axis.title.x = element_blank()) + guides(colour = "none"),
                            p.globalnetcompare.list[[2]] + theme(axis.title.x = element_blank()) + guides(colour = "none"),
                            p.globalnetcompare.list[[9]] + theme(axis.title.x = element_blank()) + guides(colour = "none"),
                            p.globalnetcompare.list[[5]] + theme(axis.title.x = element_blank()) + guides(colour = "none"),
                            p.globalnetcompare.list[[6]] + theme(axis.title.x = element_blank()) + guides(colour = "none"),
                            p.globalnetcompare.list[[7]] + guides(colour = "none"),
                            p.globalnetcompare.list[[8]],nrow = 4) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')  & theme(legend.box = "vertical", legend.position = 'bottom', plot.margin = margin(2, 4, 0, 0, "pt"))

ggsave(sprintf("TSnetcompare_NOurban_%s_%s_WSnet_%s_%s_export%i_new.eps", param1, param2, experiment.tag, scenario.tag, eperc.tag), p.ts.netcompare, width=20, height=30, units="cm", dpi=800, device=cairo_ps)

#patch-level comparison (with summary stats)
varnode.sample<-c("pop", "ag", "food", "food/pop", "nonag/pop", "import/food", "nonag", "natland")
varnodets.labels<-c(expression("Population (x "*10^{9}*")"), expression("Agricultural land (x"*10^{9}*" hectares)"), expression("Food supply (x"*10^{9}*" tonnes)"), 
                    paste("Food per capita (tonnes)", sep=''), paste("Non-agricultural land per capita (hectares)", sep=''), paste('Import dependency', sep=''),
                    expression("Non-agricultural land (x"*10^{9}*" hectares)"), expression("Natural land-states (x"*10^{9}*" hectares)"))


p.globalnodecompare.list<-lapply(1:length(varnode.sample),function(i)
{
  ggplot(data=DT.tidy, aes_string(x="ts", y=varnode.sample[i], colour=group.tag, fill=group.tag, group="interaction(runID,node)")) +
    stat_summary(data=DT.tidy, geom="ribbon", fun.max=max, fun.min=min, linetype="blank", aes_string(group=group.tag), alpha=0.25)+
    stat_summary(data=DT.tidy, geom="line", fun=mean, size=0.75, aes_string(group=group.tag))+
    scale_colour_brewer(combo.title, palette="RdYlBu", direction=1, labels=combo.labels) +
    scale_fill_brewer(combo.title, palette="RdYlBu", direction=1, labels=combo.labels) +
    scale_x_continuous(expand=c(0,0)) +
    xlab("Time-step (t)") +
    ylab(varnodets.labels[i]) +
    theme_bw()
})

if (experiment.tag=="gamma_netexpHETERO")
{
  ylims1<-c(0.05,0.06)
  ylims2<-c(0.06,0.14)
} else {
  ylims1<-c(0.05,0.06)
  ylims2<-c(0.06,0.16)
}


p.ts.nodecompare<-wrap_plots(p.globalnodecompare.list[[1]] + theme(axis.title.x = element_blank()) + guides(colour = "none", fill= "none", shape = "none"), 
                             p.globalnodecompare.list[[3]]  + ylim(ylims2) + theme(axis.title.x = element_blank()) + guides(colour = "none", fill= "none", shape = "none"), 
                             p.globalnodecompare.list[[2]] + ylim(ylims1) + theme(axis.title.x = element_blank()) + guides(colour = "none", fill= "none", shape = "none"), 
                             p.globalnodecompare.list[[8]] + theme(axis.title.x = element_blank()) + guides(colour = "none", fill= "none", shape = "none"), 
                             p.globalnodecompare.list[[4]] + guides(colour = "none", fill= "none", shape = "none"), 
                             p.globalnodecompare.list[[5]] + guides(
                               color = guide_legend(order = 0),
                               fill = guide_legend(order = 0),
                               shape = guide_legend(order = 1)
                             ), nrow = 3) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')  & theme(legend.box = "vertical", legend.position = 'bottom', plot.margin = margin(2, 4, 0, 0, "pt"))

ggsave(sprintf("TSnodecompare_NOurban_%s_%s_WSnet_%s_%s_export%i_new.eps", param1, param2, experiment.tag, scenario.tag,eperc.tag), p.ts.nodecompare, width=20, height=25, units="cm", dpi=800 , device=cairo_ps)


### node-level network attribute comparison

#functions to normalize node-level metrics, as we don't really care about the actual value just how node compares to others in it's network
normalize <- function(x){
  return(if(min(x) != max(x)) {(x-min(x)) / (max(x)-min(x))} else {1})
}

DT.scaled <- 
  DT.new[DT.new$ts==100,] %>%
  group_by(runID) %>%
  mutate(nor.degree = normalize(degree))


if (experiment.tag=="gamma_netexpHETERO")
{
  w<-7.25
  h<-7.25
  xcoords<-c(0.99, 0.99, 0.01, 0.01)
  ycoords<- c(0.06, 0.01, 0.95, 0.9)
  ylims1<-c(0.0525,0.06)
  ylims2<-c(0.1,0.15)
  ylims3<-c(0.1,0.2)
} else {
  w<-7.25
  h<-7.25
  xcoords<-c(0.99, 0.01, 0.99, 0.01)
  ycoords<- c(0.075, 0.975, 0.025, 0.925)
  ylims1<-c(0.055,0.06)
  ylims2<-c(0.1,0.15)
  ylims3<-c(0.05,0.3)
}


nodevar.sample<-c("pop", "ag", "food", "food/pop", "nonag/pop", "import/food", "nonag", "natland")
nodevar.labels<-c(expression("Population (x "*10^{9}*")"), expression("Agricultural land (x"*10^{9}*" hectares)"), expression("Food supply (x"*10^{9}*" tonnes)"), paste("Food per capita (tonnes)", sep=''), paste("Non-agricultural land per capita (hectares)", sep=''), paste('Import dependency', sep=''),expression("Non-agricultural land (x"*10^{9}*" hectares)"), expression("Natural land-states (x"*10^{9}*" hectares)"))


p.nodevar.list<-lapply(1:length(nodevar.sample),function(i)
{
  ggplot(data=DT.scaled, aes_string(y=nodevar.sample[i], x="nor.degree", colour=group.tag)) +
    geom_point(size=0.25)+
    scale_colour_brewer(combo.title, palette="RdYlBu", direction=1, labels=combo.labels) +
    theme_bw() +
    xlab("Normalized degree centrality") +
    ylab(nodevar.labels[i])  +
    guides(shape=FALSE)
  
})

p.nodevar<-wrap_plots(p.nodevar.list[[1]] + theme(axis.title.x = element_blank()) + guides(colour = "none", fill= "none", shape = "none"), 
                      p.nodevar.list[[3]] + theme(axis.title.x = element_blank()) + ylim(ylims2) + guides(colour = "none", fill= "none", shape = "none"), 
                      nrow = 2) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')

p.nodevar2<-p.nodevar.list[[5]] +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')


p.combo<-wrap_plots(p.nodevar + guides(colour = "none", fill= "none", shape = "none"), p.nodevar2 + guides(color = guide_legend(order = 0, override.aes = list(size = 5))), nrow=2) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect',heights = c(2, 1)) & theme(legend.box = "vertical", legend.position = 'bottom', plot.margin = margin(2, 6, 2, 0, "pt"))

ggsave(sprintf("DegreeCentrality_var_comparer_NOurban_%s_%s_WSnet_%s_%s_export%i_new.eps", param1, param2, experiment.tag, scenario.tag, eperc.tag), p.combo, width=20, height=23.4, units="cm", dpi=800)

DT.scaled.fin<-DT.scaled[DT.scaled$ts==max(DT.scaled$ts),]

bwidth<-0.05 
 
  p.dist.pop<-ggplot(data=DT.scaled.fin, aes_string(x="pop", group=group.tag, colour=group.tag)) +
    scale_colour_brewer(combo.title, palette="RdYlBu", direction=1, labels=combo.labels) +
    geom_freqpoly(alpha=1, binwidth=bwidth) +
    scale_x_continuous(expand=c(0,0)) +
    theme_bw()
  
  p.dist.fpc<-ggplot(data=DT.scaled.fin, aes_string(x="food/pop", group=group.tag, colour=group.tag)) +
    scale_colour_brewer(combo.title, palette="RdYlBu", direction=1, labels=combo.labels) +
    geom_freqpoly(alpha=1, binwidth=bwidth) +
    scale_x_continuous(expand=c(0,0)) +
    xlab("Food per capita (tonnes)") +
    ylab("Number of patches") +
    theme_bw()
  
  p.dist.imp<-ggplot(data=DT.scaled.fin, aes_string(x="import/food", group=group.tag, colour=group.tag)) +
    scale_colour_brewer(combo.title, palette="RdYlBu", direction=1, labels=combo.labels) +
    geom_freqpoly(alpha=1, binwidth=bwidth) +
    scale_x_continuous(expand=c(0,0)) +
    theme_bw()
 
  p1<- wrap_plots(p.globalnetcompare.list[[1]]+ ggtitle("Global trajectories") ,
                 p.globalnetcompare.list[[5]] + ylim(c(0.4,1.2) ) + guides(colour = "none"),
                 p.globalnetcompare.list[[7]] + guides(colour = "none"), nrow=1) +
    plot_layout(guides = 'collect',heights = c(1, 1, 1)) & theme(legend.box = "vertical", legend.position = 'bottom', plot.margin = margin(0, 7, 0, 2, "pt"))

  
  p2<- wrap_plots(p.nodevar.list[[4]] + guides(color = guide_legend(order = 0, override.aes = list(size = 5))) + ggtitle("Patch-level outcomes at t=100"),
                  p.nodevar.list[[2]] + guides(color = guide_legend(order = 0, override.aes = list(size = 5))), 
                  p.nodevar.list[[8]] + guides(color = guide_legend(order = 0, override.aes = list(size = 5))), nrow=1)  +
    plot_layout(guides = 'collect',heights = c(1, 1, 1)) & theme(legend.box = "vertical", legend.position = 'bottom', plot.margin = margin(0, 7, 0, 2, "pt"))
  
  p<-wrap_plots(p1, p2, nrow=2) + plot_annotation(tag_levels = 'a') & theme(plot.title = element_text(size = 9, face = "bold"))

 png(sprintf("combodemo_%s_new.png", experiment.tag), width=20, height=18, units="cm", res=500)
 print(p)
 dev.off()