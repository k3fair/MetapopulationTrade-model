#     Metapop_ParameterPlanes_Viz.R visualizes the output from Metapop_NetworkExperiment.m
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
library(reldist) #for gini coeff
library(data.table) #for dealing with large datasets
library(dplyr) #for IDing runs
library(ggforce) #to draw hulls, etc
library(patchwork) #for multipanel figures
library(RColorBrewer) # for custom palette
library(scales) #For log10 scale

#Sets parameters for filename to be read in (leave as-is) 
param.list<-c("betaA", "betaI", "gammaA", "gammaI")
param1<-param.list[1]
param2<-param.list[3]

path = getwd() #Set to wherever you have the files from NetworkStatistics_forNetworkExperiment folder stored

netsize<-c(100)

for (h in 1:length(netsize))
{
file.names <- dir(path, pattern =sprintf("WSnetstats_1_%iN_*",netsize[h]))

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
experiment.tag<-"WSnetexp"
scenario.tag<-"lowyield"
eperc.tag<-16

df.N100.neigh26<-fread(sprintf("PPlane_NOurban_%s_%s_sWnet_N100_neigh26_%s_%s_export%i_new.csv", param1, param2, experiment.tag, scenario.tag, eperc.tag),
                 header=FALSE,stringsAsFactors = FALSE);
colnames(df.N100.neigh26)<-c("N", "neigh", "rewire", "gammaI", "betaI", "gammaA", "betaA", "node", "ts", "yield", "pop", "ag", "food", "realization", "import"  )

df.N100.neigh17<-fread(sprintf("PPlane_NOurban_%s_%s_sWnet_N100_neigh17_%s_%s_export%i_new.csv", param1, param2, experiment.tag, scenario.tag, eperc.tag),
                       header=FALSE,stringsAsFactors = FALSE);
colnames(df.N100.neigh17)<-c("N", "neigh", "rewire", "gammaI", "betaI", "gammaA", "betaA", "node", "ts", "yield", "pop", "ag", "food", "realization", "import"  )

df.comb<-rbind(df.N100.neigh26, df.N100.neigh17)

#clean up from reading in data
rm("df.N100.neigh17", "df.N100.neigh26")

DT.0<-merge(df.comb, stats.comb.tidy, by=c("node","N", "neigh", "rewire", "realization"))

#clean up from combining data
rm("df.comb", "stats.comb.tidy")

#ID rows by the simulation they correspond to
DT.tidy <- data.table(DT.0, key=c("N", "neigh", "rewire", "realization", param1,param2))
DT.tidy[, runID:=.GRP, by=key(DT.tidy)]

#clean up after tidying
rm("DT.0") 

#add in column for nonag land
DT.tidy$nonag<-(9.625323203/DT.tidy$N)-DT.tidy$ag
#add in column for nonag land excluding urban area # new added may 5 2021 #
DT.tidy$natland<-(9.625323203/DT.tidy$N) - DT.tidy$ag - (DT.tidy$pop*0.06)

#Calculate Gini index values (weighting by population size)
DT.fin.0 <- DT.tidy[(DT.tidy$ts==100),]
system.time(for (i in 1:length(unique(DT.fin.0$runID)))
{
    sub<-DT.fin.0[DT.fin.0$runID==unique(DT.fin.0$runID)[i],]
    sub$gini.food.final<-gini(sub$food/sub$pop, sub$pop)
    sub$gini.land.final<-gini(sub$nonag/sub$pop, sub$pop)
    sub$gini.import.final<-gini(sub$food/sub$import, sub$pop)
    
    if (i==1)
    {
      DT.fin<-sub
    } else
    {
      DT.fin<-rbind(DT.fin, sub)
    }
  #print(i)
})

rm("DT.fin.0", "sub")

### Visualize global behaviour
df.sum<-DT.tidy %>%
  group_by(N,neigh,rewire,realization,gammaI,betaI,gammaA,betaA,ts,sw.measure.TJHBL,density,diameter,av.pl,av.deg,edges,avg.local.cc,avg.local.transitivity,n.components,runID) %>%
  summarize(import.dep=weighted.mean(import/food, pop), pop = sum(pop), ag=sum(ag), food=sum(food), fpc=weighted.mean(food/pop, pop), nonag=sum(nonag), lpc=weighted.mean(nonag/pop, pop),natland=sum(natland))

#create final timestep data frame
df.sum.fin<-merge(df.sum[df.sum$ts==100,], unique(DT.fin[,c("runID", "gini.food.final", "gini.land.final","gini.import.final")]), by="runID")

#create custom palette based on the RdYlBu scheme
rdylbu.pal<-palette(brewer.pal(5,"RdYlBu"))
my.pal<-c(rdylbu.pal[1], rdylbu.pal[5])

### global network rewire comparison
var.sample<-c("pop", "ag", "food", "import.dep",  "fpc", "lpc", "gini.food.final", "gini.land.final", "natland")
labels.var<-c(expression("Population (x "*10^{9}*")"), expression("Agricultural land (x"*10^{9}*" hectares)"), expression("Food supply (x"*10^{9}*" tonnes)"), paste('Mean import \ndependency', sep=''), 
  paste("Mean food per capita (tonnes)", sep=''), paste("Mean non-ag. land per capita (hectares)", sep=''), paste('Gini index - food per capita', sep=''),
  paste('Gini index - non-ag. land per capita', sep=''), expression("Natural land-states (x"*10^{9}*" hectares)"))

df.rewire.sum<-df.sum.fin[,c(4,12,20:29)] %>%
  group_by(rewire,density) %>%
  summarise_at(c("pop", "ag", "food", "fpc", "lpc", "gini.food.final", "gini.land.final", "import.dep", "natland"), list(~min(.), ~max(.), ~mean(.)))

p.rewire.list<-lapply(1:length(var.sample),function(i)
{
  ggplot(data=df.rewire.sum, aes_string(y=paste0(var.sample[i], "_mean"), x="rewire", group="density", colour="as.factor(round(density,3))")) +
    geom_errorbar(aes_string(ymin=paste0(var.sample[i], "_min"), ymax=paste0(var.sample[i], "_max")), width=.15) +
    geom_point() +
    scale_colour_manual("Network density (d)", values=my.pal) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("Rewiring probability (p)") +
    ylab(labels.var[i]) +
    theme_bw() +
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
})

p.rewire<-wrap_plots(p.rewire.list[[1]] + theme(axis.title.x = element_blank()),
                     p.rewire.list[[3]] + theme(axis.title.x = element_blank())+ ylim(c(12.75,13.25)), 
                     p.rewire.list[[2]] + theme(axis.title.x = element_blank()) + ylim(c(5.5,5.75)),
                     p.rewire.list[[9]] + theme(axis.title.x = element_blank()) + ylim(c(3,3.25)),
                     p.rewire.list[[5]] + theme(axis.title.x = element_blank()) + ylim(c(0.75,1)), 
                     p.rewire.list[[6]] + ylim(c(0.2,0.4)), 
                     p.rewire.list[[7]] + ylim(c(0,0.35)), p.rewire.list[[8]] + ylim(c(0,0.35)), guide_area(), nrow = 3) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')  & theme(text = element_text(size=9), legend.background = element_blank(),
                                           plot.margin = margin(2, 2, 2, 0, "pt"), legend.box.margin=margin(-8,2,0,2))

png("Fig3_Rewire_comparer_new.png", width=20, height=20, units="cm", res=300)
print(p.rewire)
dev.off()

#####

p.rewire.list.pres<-lapply(1:length(var.sample),function(i)
{
  ggplot(data=df.rewire.sum, aes_string(y=paste0(var.sample[i], "_mean"), x="rewire", group="density", colour="as.factor(round(density,3))")) +
    geom_errorbar(aes_string(ymin=paste0(var.sample[i], "_min"), ymax=paste0(var.sample[i], "_max")), width=.15) +
    geom_point() +
    scale_colour_manual("Network density", values=my.pal) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab("Rewiring probability") +
    ylab(labels.var[i]) +
    theme_bw() +
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
})

netmetric.sample<-c("sw.measure.TJHBL", "av.pl", "avg.local.cc")
netmetric.labels<-c("Small-world measure", "Average path length", "Average clustering coefficient")

p.netchar.list<-lapply(1:length(netmetric.sample),function(i)
{
  
  ggplot(data=df.sum.fin, aes_string(y=netmetric.sample[i], x="rewire", group="density", colour="as.factor(round(density,3))",fill="as.factor(round(density,3))")) +
    stat_summary(geom="ribbon", fun.min="min", fun.max="max", alpha=0.1) +
    stat_summary(geom="line", fun=mean, linetype="dashed")+
    scale_colour_manual("Network density (d)", values=my.pal) +
    scale_fill_manual("Network density (d)", values=my.pal) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), expand=c(0,0)) +
    scale_shape_manual("Network density", values=c(1,4)) +
    xlab("Rewiring probability (p)") +
    ylab(netmetric.labels[i])+
    theme_bw()
  
})

p.netchar<-wrap_plots(p.netchar.list,nrow = 1) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.background = element_blank(), plot.margin = margin(2, 7, 0, 0, "pt"))

png(sprintf("NetworkMetric_rewire_comparer_NOurban_%s_%s_WSnet_%s_%s_new.png", param1, param2, experiment.tag, scenario.tag), width=20, height=8.5, units="cm", res=300)
print(p.netchar)
dev.off()

df.sum<-df.sum.fin %>% group_by(density,rewire) %>% summarize(mean.sw=mean(sw.measure.TJHBL))

#functions to normalize node-level metrics so we can make relative comparisons
normalize <- function(x){
  return(if(min(x) != max(x)) {(x-min(x)) / (max(x)-min(x))} else {1})
}

DT.scaled <- 
  DT.fin[DT.fin$rewire %in% c(0.0001, 0.01, 1),] %>%
  group_by(runID) %>%
  mutate(degree = normalize(degree), betweenness=normalize(betweenness), closeness=normalize(closeness), eigen=normalize(eigen), local.cc=normalize(local.cc))

nodemetric.sample<-c("degree", "betweenness", "closeness", "eigen", "local.cc")
nodemetric.labels<-c("Degree centrality", "Betweenness centrality", "Closeness centrality", "Eigenvector centrality", "Clustering coefficient")
p.nodechar.list<-lapply(1:length(nodemetric.sample),function(i)
{
  ggplot(data=DT.scaled, aes_string(x=nodemetric.sample[i], y="food/pop")) +
    geom_point(size=0.5) +
    geom_mark_hull(aes(colour = as.factor(round(density,3)), fill = as.factor(round(density,3)),  linetype=as.factor(rewire)), expand = unit(2.5, "mm"), concavity=2.5, alpha=0.1, size=0.75) +
    scale_colour_manual("Network density (d)", values=my.pal) +
    scale_fill_manual("Network density (d)", values=my.pal) +
    scale_linetype_manual("Rewiring probability (p)", values=c("solid","dotted", "blank")) +
    theme_bw() +
    ylab("Food per capita (tonnes)") +
    xlab(nodemetric.labels[i])+
    theme(legend.background = element_rect(fill="gray95", colour="black",
                                            linetype="solid"))

})


p.nodechar<-wrap_plots(p.nodechar.list[[1]] + guides(colour=FALSE, fill=FALSE, linetype=FALSE),
                       p.nodechar.list[[2]] + theme(axis.title.y = element_blank()) + guides(colour=FALSE, fill=FALSE, linetype=FALSE),
                       p.nodechar.list[[3]] + guides(colour=FALSE, fill=FALSE, linetype=FALSE),
                       p.nodechar.list[[4]] + theme(axis.title.y = element_blank()) + guides(colour=FALSE, fill=FALSE, linetype=FALSE),
                       p.nodechar.list[[5]], nrow = 3) +
  guide_area() + 
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect') & theme(legend.position = 'right', plot.margin = margin(2, 4, 0, 0, "pt"))

png(sprintf("NodeMetric_fpc_comparer_NOurban_%s_%s_WSnet_%s_%s_export%i_new.png", param1, param2, experiment.tag, scenario.tag,eperc.tag), width=17.6, height=23.4, units="cm", res=500)
print(p.nodechar)
dev.off()
