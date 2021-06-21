#     Metapop_ParameterPlanes_Viz.R visualizes the output from Metapop_ParameterPlanes.m
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

library(ggplot2) #For plotting
library(reldist) #for gini coeff
library(data.table) #for dealing with large datasets
library(dplyr) #for IDing runs
library(patchwork) #For creating multipanel figures

#### Pick plot to generate: "betagamma" i.e. (beta, gamma) planes or "betabeta" i.e. (beta_A, beta_I) planes
select.plot<-"betagamma"
#### Pick yield scenario to assume (lowyield or highyield)
scenario.tag<-"highyield"

#picks parameters to compare, based on your choice of select.plot
param.list<-c("betaA", "betaI", "gammaA", "gammaI")

if(select.plot=="betabeta"){
param1<-param.list[1];
param2<-param.list[2]; } else {
  param1<-param.list[1];
  param2<-param.list[3];
}
outparam1<-param.list[3]
outparam2<-param.list[4]

if (param1==param.list[1] & param2==param.list[3])
{
  xlabel<- expression("Demand responsiveness ("*gamma*")")
  ylabel<-expression("Demand threshold ("*beta*")")
}

if (param1==param.list[1] & param2==param.list[2])
{
  xlabel<-expression("Import demand threshold ("*beta^{I}*")")
  ylabel<-expression("Agriculture demand threshold ("*beta^{A}*")")
}

#read in ts data
if(select.plot=="betabeta"){
experiment.tag<-"gamma7_5"; } else {
  experiment.tag<-"maxbeta2_5"; 
}

eperc.tag<-16

df.N10<-fread(sprintf("PPlane_NOurban_%s_%s_sWnet_N10_%s_%s_export%i.csv", param1, param2, experiment.tag, scenario.tag, eperc.tag),
                 header=FALSE,stringsAsFactors = FALSE);
colnames(df.N10)<-c("N", "neigh", "rewire", "gammaI", "betaI", "gammaA", "betaA", "node", "ts", "yield", "pop", "ag", "food", "import" )

df.N100<-fread(sprintf("PPlane_NOurban_%s_%s_sWnet_N100_%s_%s_export%i.csv", param1, param2, experiment.tag, scenario.tag, eperc.tag),
                 header=FALSE,stringsAsFactors = FALSE);
colnames(df.N100)<-c("N", "neigh", "rewire", "gammaI", "betaI", "gammaA", "betaA", "node", "ts", "yield", "pop", "ag", "food", "import" )

df.comb<-rbind(df.N10, df.N100)

##clean up from reading in data
rm("df.N10", "df.N100")

#ID rows by the simulation they correspond to
DT.tidy <- data.table(df.comb, key=c("N", "neigh", "rewire",param1,param2))
DT.tidy[, runID:=.GRP, by=key(DT.tidy)]

#clean up after tidying
rm("df.comb")

#add in column for nonag land
DT.tidy$nonag<-(9.625323203/DT.tidy$N)-DT.tidy$ag
#add in column for nonag land excluding urban area # new added may 5 2021 #
DT.tidy$natland<-(9.625323203/DT.tidy$N) - DT.tidy$ag - (DT.tidy$pop*0.06)

#Calculate Gini index values (weighting by population size)
DT.fin.0 <- DT.tidy[(DT.tidy$ts==100),]
for (i in 1:length(unique(DT.fin.0$runID)))
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
}

rm("DT.fin.0", "sub")

### Summarize global behaviour
df.sum<-DT.tidy %>%
  group_by(N,neigh,rewire,gammaI,betaI,gammaA,betaA,ts,runID) %>%
  summarize(import.dep=weighted.mean(import/food, pop), pop = sum(pop), ag=sum(ag), food=sum(food), fpc=weighted.mean(food/pop, pop), nonag=sum(nonag), lpc=weighted.mean(nonag/pop, pop), natland=sum(natland))

#create final timestep
df.sum.fin<-merge(df.sum[df.sum$ts==100,], unique(DT.fin[,c("runID","gini.food.final","gini.land.final","gini.import.final")]), by="runID")
df.sum.fin$importrel<-df.sum.fin$import/df.sum.fin$food

###Generate plots

var.sample.final<-c("pop", "ag", "food", "fpc", "lpc", "gini.food.final", "gini.land.final", "import.dep", "natland")
var.sample.labels<-c(expression("Population (x "*10^{9}*")"), expression("Agricultural land (x"*10^{9}*" hectares)"), expression("Food supply (x"*10^{9}*" tonnes)"),
                     paste("Mean food, \nper capita (tonnes)", sep=''), paste("Mean non-ag. \nland, per capita (hectares)", sep=''), paste('Gini index - food \nper capita', sep=''),
                     paste('Gini index - non-ag. \nland per capita', sep=''), paste('Mean import \ndependency', sep=''), expression("Natural land-states (x"*10^{9}*" hectares)"))

net.sample<-c(10,100)

p.final.list<-lapply(1:length(var.sample.final),function(i)
{
  lapply(1:length(net.sample),function(j) {
      
      ggplot(data=df.sum.fin[(df.sum.fin$N==net.sample[j]),]) +
        geom_point(aes_string(x=param2, y=param1, colour=var.sample.final[i])) +
      scale_colour_distiller(var.sample.labels[i], palette="RdYlBu", direction=1,
                             limits=c(ifelse((var.sample.final[i] %in% c("gini.food.final", "gini.land.final")),min(df.sum.fin[c("gini.food.final", "gini.land.final")], na.rm=TRUE), min(df.sum.fin[[var.sample.final[[i]]]], na.rm=TRUE)),
                                      ifelse((var.sample.final[i] %in% c("gini.food.final", "gini.land.final")),max(df.sum.fin[c("gini.food.final", "gini.land.final")], na.rm=TRUE), max(df.sum.fin[[var.sample.final[[i]]]], na.rm=TRUE))),
                             na.value = "black") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      xlab(xlabel) +
      ylab(ylabel) +
      theme(panel.background = element_rect(fill = "black",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "black"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "black")) +
      guides(colour= guide_colourbar(frame.colour="black", ticks.colour = "black"))
      
  })
  
})


p.pop<-wrap_plots(p.final.list[[1]][[1]] + ggtitle("N=10") + theme(axis.title.x = element_blank()), p.final.list[[1]][[2]] + ggtitle("N=100") + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')
p.ag<-wrap_plots(p.final.list[[2]][[1]] + theme(axis.title.x = element_blank()) , p.final.list[[2]][[2]]  + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')
p.food<-wrap_plots(p.final.list[[3]][[1]]+ ggtitle("N=10") + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), p.final.list[[3]][[2]]  + ggtitle("N=100") + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')
p.fpc<-wrap_plots(p.final.list[[4]][[1]] + theme(axis.title.x = element_blank()), p.final.list[[4]][[2]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')
p.lpc<-wrap_plots(p.final.list[[5]][[1]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), p.final.list[[5]][[2]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')
p.ginifood<-wrap_plots(p.final.list[[6]][[1]], p.final.list[[6]][[2]]+ theme(axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')
p.giniland<-wrap_plots(p.final.list[[7]][[1]]+ theme(axis.title.y = element_blank()), p.final.list[[7]][[2]]+ theme(axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')
p.import<-wrap_plots(p.final.list[[8]][[1]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), p.final.list[[8]][[2]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')
p.natland<-wrap_plots(p.final.list[[9]][[1]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), p.final.list[[9]][[2]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),nrow = 1) +
  plot_layout(guides = 'collect')

p.combined.test<-wrap_plots(p.pop,  p.food, p.ag, p.natland, p.fpc, p.lpc, p.ginifood, p.giniland,  nrow=4) +
  plot_annotation(tag_levels = 'a') +
   plot_layout(ncol=2) & theme(legend.justification = "center", aspect.ratio=1)  & theme(text = element_text(size=8), 
                                                                                        legend.position = 'bottom', plot.margin = margin(0, 1, 0, 1, "pt"),
                                                                                  legend.box.margin=margin(-12,0,-6,0))

if (select.plot=="betabeta" && scenario.tag=="lowyield") {
png("Fig2_PPplaneVARS_betaA_betaI_lowyield.png", width=20, height=24, units="cm", res=500)
print(p.combined.test)
dev.off() } 

if (select.plot=="betagamma" && scenario.tag=="lowyield") {
  png("Fig1_PPplaneVARS_betaA_gammaA_lowyield.png", width=20, height=24, units="cm", res=500)
  print(p.combined.test)
  dev.off()
}

if (select.plot=="betagamma" && scenario.tag=="highyield") {
  png("SIFig_PPplaneVARS_betaA_gammaA_highyield.png", width=20, height=24, units="cm", res=500)
  print(p.combined.test)
  dev.off()
}