rm(list=ls())
source("./Legacy_functions.R")

# ------------------------------------- Main analyses and figures ----
## >> Global map----

d_irregular=read.table("./Data/data_sites_CLIM.csv",sep=";")
d_regular=read.table("./Data/data_sites_Regular_Buxton.csv",sep=";")  

past_biome=read_sf("../Data/Paleo_clim/Desert_past/world_cut.shp")

p=ggplot(past_biome%>%Change_biome_name_LGM(.)%>%
           mutate(., Desert=sapply(1:nrow(.),function(x){
             if (any(grep("desert",.$Save_name[x]))){
               return(.$Save_name[x])
             }else{
               return("Non-desert biomes")
             }
           }))) + 
  geom_sf(aes(fill=as.factor(Desert))) + 
  geom_point(data=tibble(Long=c(d_regular$Longitude,d_irregular$Longitude),
                         Lat=c(d_regular$Latitude,d_irregular$Lat),
                         Type=c(rep("Dataset2: Regular patterns in the Sahel",each=nrow(d_regular)),
                                rep("Dataset1: Global survey",each=nrow(d_irregular)))),
             aes(x=Long,y=Lat,shape=Type)) + 
  scale_fill_manual(values=c("#E8DCBB","#F5D992","#DCB95D","#B7912E","#7D6015","grey"))+
  scale_shape_manual(values=c(17,15))+
  coord_sf()+
  the_theme2+
  guides(fill=guide_legend(nrow=3,byrow=TRUE))+
  ggtitle("LGM biomes")+
  theme( legend.box="vertical")+
  labs(fill="",shape="",x="Longitude",y="Latitude")

ggsave("./Figures/Distribution_sites_regular_irregular.pdf",p,width = 6,height = 7)


## >> BIOCOM partition variable and RF ----

id_metric_kept=c(1,2,3,8,12,14,18)
for (ID_PC in 1:2){
  
  d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
    Closer_to_normality(.)%>%
    Perform_PCA_spatial_struc(.,F)
  
  d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
  d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
  
  d=d%>%
    dplyr::mutate(.,
                  Desert_Current=as.numeric(Desert_Current=="Desert"),
                  Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
    dplyr::select(.,-Biome_LGM,-Biome_Current)
  
  d$Struct1=d[,paste0("Struct",ID_PC)]
  
  #LGM clim
  d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))[c(id_metric_kept,20)]])]
  #Holocene clim
  d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
  
  #Current clim
  d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept,20)]])]
  
  assign(paste0("d_",ID_PC),
         Plot_variance_partition(varpart(d[,c("Struct1")],d_LGM,d_Holo,d_Current),
                                 bar = T,four_groups = F))
  
  combined_d=cbind(cbind(d_LGM),
                   d_Holo,
                   cbind(d_Current),
                   d$Struct1)
  colnames(combined_d)[ncol(combined_d)]="Struct1"
  train_id=sample(1:nrow(combined_d),round(.7*nrow(combined_d)))
  
  #Changing to factors
  if (any(colnames(combined_d)=="Desert_Current")){
    combined_d[,"Desert_Current"]=as.factor(combined_d[,"Desert_Current"])
  }
  
  if (any(colnames(combined_d)=="Desert_LGM")){
    combined_d[,"Desert_LGM"]=as.factor(combined_d[,"Desert_LGM"])
  }
  
  set.seed(123)
  train_id=sample(1:nrow(d),round(.7*nrow(d)))
  
  BRT_model = gbm(Struct1 ~ .,
                      distribution = "gaussian",
                      data = combined_d[train_id,],
                      interaction.depth = 7,
                      n.trees = 500,
                      cv.folds = 10,
                      n.minobsinnode = 4, 
                      shrinkage = .1,
                      train.fraction = .7,
                      n.cores = 1,
                      bag.fraction = .8)
  
  assign(paste0("pRF_",ID_PC),
         ggplot(as.data.frame(summary(BRT_model,plot=F))%>%
                  add_column(., mean_imp=.$rel.inf)%>%
                  Change_name_RF(.,"var")%>%
                  dplyr::arrange(., mean_imp)%>%
                  add_column(., Order_stat=1:nrow(.))%>%
                  dplyr::filter(., Order_stat>(max(Order_stat)-12))%>%
                  mutate(.,Name_stat = fct_reorder(Name_stat, Order_stat)))+
           geom_bar(aes(y=mean_imp,group=interaction(Name_stat,Group_variable),x=Name_stat,
                        fill=Group_variable),stat="identity")+
           
           coord_flip()+
           the_theme2+
           scale_fill_manual(values=c("LGM legacy"="#FFF9BF", "Holocene legacy"="#FFD2A0","Current climate"="pink",
                                      "Shared past clim."="#DAB2BA","Other shared variance"="#E8E8E8"))+
           labs(y="Mean Importance",x="",fill=""))
}  

pvariance=ggarrange(d_1$Figure+ggtitle("PC1: mean patch-size, cover")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "none"),
                    d_2$Figure+ggtitle("PC2: spatial aggregation")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "none"),
                    ncol=2,common.legend = T,legend="none")

PCA_stat=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.,return_PCA = T)

p0=plotPCA(PCA_stat,1,2,F,xLim = c(-3,3),yLim = c(-3,3),colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")

p_RF=ggarrange(
  pRF_1+ggtitle("PC1: mean patch-size, cover")+theme(legend.position = "none"),
  pRF_2+ggtitle("PC2: spatial aggregation")+theme(legend.position = "none"),
  ncol=2,common.legend = T,legend = "none")


p_legend=get_legend(ggplot(tibble(ID=1:5,Name=c("LGM legacy","Holocene legacy",
                                                "Current climate","Shared past clim.",
                                                "Other shared variance")))+
                      geom_bar(aes(x=1,y=1,fill=as.factor(Name)),stat="identity")+
                      scale_fill_manual(values=c("LGM legacy"="#FFF9BF",
                                                 "Holocene legacy"="#FFD2A0",
                                                 "Current climate"="pink",
                                                 "Shared past clim."="#DAB2BA",
                                                 "Other shared variance"="#E8E8E8"))+
                      labs(fill="")+theme(legend.position = "bottom")
)

p_tot=ggarrange(
  ggarrange(p0+the_theme2+theme(legend.position = "none"),
            pvariance,ncol=2,widths = c(.7,1),labels = letters[1:2],align = "hv"),p_legend,nrow=2,heights = c(1,.15),align = "hv")

ggsave("./Figures/Figure_irregular.pdf",p_tot,width = 9,height = 4.5)

p_tot_RF=ggarrange(p_RF,p_legend,nrow=2,heights = c(1,.15),align = "hv")

ggsave("./Figures/SI/Figure_RF_irregular.pdf",p_tot_RF,width = 8,height = 4.5)



save_names_LGM=colnames(d_LGM)[-length(colnames(d_LGM))]
save_names_Holo=colnames(d_Holo)
save_names_Current=colnames(d_Current)[-length(colnames(d_Current))]

## >> Regular partition variable and RF ----

id_metric_kept=c(1,2,3,8,12,14,18)
for (ID_PC in 1:2){
  
  d=read.table("./Data/data_sites_Regular_Buxton.csv",sep=";")%>%
    Closer_to_normality(.)%>%
    Perform_PCA_spatial_struc(.,F)%>%
    dplyr::select_if(function(x) length(unique(x)) > 1 & !any(is.na(x)))%>%
    mutate(across(colnames(dplyr::select_if(., is.numeric)), ~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
  
  d$Struct1=d[,paste0("Struct",ID_PC)]
  
  #LGM clim
  d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))[c(id_metric_kept)]])]
  
  #Holocene clim
  d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
  
  #Current clim
  d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept)]])]
  
  assign(paste0("d_",ID_PC),
         Plot_variance_partition(varpart(d[,c("Struct1")],d_LGM,d_Holo,d_Current),
                                 bar = T,four_groups = F))
  combined_d=cbind(cbind(d_LGM),
                   d_Holo,
                   cbind(d_Current),
                   d$Struct1)
  colnames(combined_d)[ncol(combined_d)]="Struct1"
  train_id=sample(1:nrow(combined_d),round(.7*nrow(combined_d)))
  
  #Changing to factors
  if (any(colnames(combined_d)=="Desert_Current")){
    combined_d[,"Desert_Current"]=as.factor(combined_d[,"Desert_Current"])
  }
  
  if (any(colnames(combined_d)=="Desert_LGM")){
    combined_d[,"Desert_LGM"]=as.factor(combined_d[,"Desert_LGM"])
  }
  
  set.seed(123)
  BRT_model = gbm(Struct1 ~ .,
                  distribution = "gaussian",
                  data = combined_d[train_id,],
                  interaction.depth = 7,
                  n.trees = 500,
                  cv.folds = 10,
                  n.minobsinnode = 4, 
                  shrinkage = .1,
                  train.fraction = .7,
                  n.cores = 1,
                  bag.fraction = .8)
  
  assign(paste0("pRF_",ID_PC),
         ggplot(as.data.frame(summary(BRT_model,plot=F))%>%
                  add_column(., mean_imp=.$rel.inf)%>%
                  Change_name_RF(.,"var")%>%
                  dplyr::arrange(., mean_imp)%>%
                  add_column(., Order_stat=1:nrow(.))%>%
                  dplyr::filter(., Order_stat>(max(Order_stat)-12))%>%
                  mutate(.,Name_stat = fct_reorder(Name_stat, Order_stat)))+
           geom_bar(aes(y=mean_imp,group=interaction(Name_stat,Group_variable),x=Name_stat,
                        fill=Group_variable),stat="identity")+
           
           coord_flip()+
           the_theme2+
           scale_fill_manual(values=c("LGM legacy"="#FFF9BF", "Holocene legacy"="#FFD2A0","Current climate"="pink",
                                      "Shared past clim."="#DAB2BA","Other shared variance"="#E8E8E8"))+
           labs(y="Mean Importance",x="",fill=""))
}  
  
pvariance=ggarrange(d_1$Figure+ggtitle("PC1, spatial structure")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_2$Figure+ggtitle("PC2, spatial structure")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    ncol=2,common.legend = T,legend="none")

PCA_stat=read.table("./Data/data_sites_Regular_Buxton.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.,return_PCA = T)
  
p0=plotPCA(PCA_stat,1,2,F,xLim = c(-3,3),yLim = c(-3,3),
           colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")

p_RF=ggarrange(
  pRF_1+ggtitle("PC1, spatial structure")+theme(legend.position = "none"),
  pRF_2+ggtitle("PC2, spatial structure")+theme(legend.position = "none"),
  ncol=2,common.legend = T,legend = "none")


p_legend=get_legend(ggplot(tibble(ID=1:5,Name=c("LGM legacy","Holocene legacy",
                                                "Current climate","Shared past clim.",
                                                "Other shared variance")))+
                      geom_bar(aes(x=1,y=1,fill=as.factor(Name)),stat="identity")+
                      scale_fill_manual(values=c("LGM legacy"="#FFF9BF",
                                                 "Holocene legacy"="#FFD2A0",
                                                 "Current climate"="pink",
                                                 "Shared past clim."="#DAB2BA",
                                                 "Other shared variance"="#E8E8E8"))+
                      labs(fill="")+theme(legend.position = "bottom")
)

p_tot=ggarrange(
  ggarrange(p0+the_theme2+theme(legend.position = "none"),
            pvariance,ncol=2,widths = c(.7,1),labels = letters[1:2],align = "hv"),p_legend,nrow=2,heights = c(1,.15),align = "hv")

ggsave("./Figures/Figure_regular.pdf",p_tot,width = 9,height = 4.5)

p_tot_RF=ggarrange(p_RF,p_legend,nrow=2,heights = c(1,.15),align = "hv")

ggsave("./Figures/SI/Figure_RF_regular.pdf",p_tot_RF,width = 8,height = 4.5)


## >> SEM indirect effects without desert factors ----

d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.)%>%
  mutate(., Sand=bcPower(.$Sand,6))%>%
  mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
d=d%>%
  dplyr::mutate(.,
                Desert_Current=as.factor(as.numeric(Desert_Current=="Desert")),
                Desert_LGM=as.factor(as.numeric(Desert_LGM=="Desert")))%>%
  dplyr::select(.,-Biome_LGM,-Biome_Current)

library(piecewiseSEM)
d_site=dplyr::distinct(d,Site_ID,.keep_all = T)

mod_MF=lm(MF~Isothermality_Holocene+
              MAP_Holocene+MAP_Current+
              MAP_LGM+MAT_LGM+Isothermality_Current+
              MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
              Isothermality_LGM,
            data = d_site,na.action = na.fail)

mod_sand=lm(Sand~Isothermality_Holocene+
            MAP_Holocene+MAP_Current+
            MAP_LGM+MAT_LGM+Isothermality_Current+
            MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
            Isothermality_LGM,
          data = d_site,na.action = na.fail)

mod_Struct1=lmer(Struct1~Isothermality_Holocene+
                   MAP_Holocene+MAP_Current+
                   MAP_LGM+MAT_LGM+Isothermality_Current+
                   MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
                   Isothermality_LGM+Sand+MF+
                   (1|Site_ID),data = d,na.action = na.fail,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_Struct2=lmer(Struct2~Isothermality_Holocene+
                   MAP_Holocene+MAP_Current+
                   MAP_LGM+MAT_LGM+Isothermality_Current+
                   MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
                   Isothermality_LGM+Sand+MF+
                   (1|Site_ID),data = d,na.action = na.fail,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

d_SEM=psem(
  mod_sand,
  mod_MF,
  mod_Struct1,
  mod_Struct2,
  MF%~~%Sand
)

LMERConvenienceFunctions::mcp.fnc(mod_Struct1) #-> ok
LMERConvenienceFunctions::mcp.fnc(mod_Struct2) #-> ok
DHARMa::simulateResiduals(mod_sand,plot=T) #-> ok
DHARMa::simulateResiduals(mod_MF,plot=T) #-> small quantile deviation but thats ok

SEM_boot = bootEff(d_SEM, R = 1000, seed = 13, parallel = "snow",ran.eff = "Site_ID",ncpus = 40)
saveRDS(SEM_boot,"./Data/SEM_boot_full.rds")

SEM_boot=readRDS("./Data/SEM_boot_full.rds")
Effects = semEff(SEM_boot)

#Getting information about indirect effects through direct and total effects
#Aggregating in a tibble
get_bootstrapped_pval=function(x){
  return(ifelse(length(which(x>0))/length(x)>.5,length(which(x<0))/length(x),length(which(x>0))/length(x)))
}

Direct_effects=tibble()
for (predictors in c("Struct1","Struct2","MF","Sand")){
  effect_pred=Effects$Effects$Bootstrapped[[predictors]]
  Direct_effects=rbind(Direct_effects,
                       tibble(q1=apply(effect_pred$Direct[,-1],2,quantile,.025),
                              q3=apply(effect_pred$Direct[,-1],2,quantile,.975),
                              q2=apply(effect_pred$Direct[,-1],2,median),
                              pval=apply(effect_pred$Direct[,-1],2,get_bootstrapped_pval),
                              Term=names(apply(effect_pred$Direct[,-1],2,quantile,.025)),
                              Response=predictors))
}

effect_on_struct1=Effects$Effects$Bootstrapped$Struct1
effect_on_struct2=Effects$Effects$Bootstrapped$Struct2
effect_on_sand=Effects$Effects$Bootstrapped$Sand
effect_on_MF=Effects$Effects$Bootstrapped$MF

Total_effects=rbind(
  tibble(q1=apply(effect_on_struct1$Total[,-1],2,quantile,.025),
         q3=apply(effect_on_struct1$Total[,-1],2,quantile,.975),
         q12=apply(effect_on_struct1$Total[,-1],2,quantile,.05),
         q32=apply(effect_on_struct1$Total[,-1],2,quantile,.95),
         q2=apply(effect_on_struct1$Total[,-1],2,median),
         pval=apply(effect_on_struct1$Total[,-1],2,get_bootstrapped_pval),
         Term=names(apply(effect_on_struct1$Total[,-1],2,quantile,.025)),
         Response="Struct1"),
  tibble(q1=apply(effect_on_struct2$Total[,-1],2,quantile,.025),
         q3=apply(effect_on_struct2$Total[,-1],2,quantile,.975),
         q12=apply(effect_on_struct2$Total[,-1],2,quantile,.05),
         q32=apply(effect_on_struct2$Total[,-1],2,quantile,.95),
         q2=apply(effect_on_struct2$Total[,-1],2,median),
         pval=apply(effect_on_struct2$Total[,-1],2,get_bootstrapped_pval),
         Term=names(apply(effect_on_struct2$Total[,-1],2,quantile,.025)),
         Response="Struct2"),
  tibble(q1=apply(effect_on_sand$Total[,-1],2,quantile,.025),
         q3=apply(effect_on_sand$Total[,-1],2,quantile,.975),
         q12=apply(effect_on_sand$Total[,-1],2,quantile,.05),
         q32=apply(effect_on_sand$Total[,-1],2,quantile,.95),
         q2=apply(effect_on_sand$Total[,-1],2,median),
         pval=apply(effect_on_sand$Total[,-1],2,get_bootstrapped_pval),
         Term=names(apply(effect_on_sand$Total[,-1],2,quantile,.025)),
         Response="Sand"),
  tibble(q1=apply(effect_on_MF$Total[,-1],2,quantile,.025),
         q3=apply(effect_on_MF$Total[,-1],2,quantile,.975),
         q12=apply(effect_on_MF$Total[,-1],2,quantile,.05),
         q32=apply(effect_on_MF$Total[,-1],2,quantile,.95),
         q2=apply(effect_on_MF$Total[,-1],2,median),
         pval=apply(effect_on_MF$Total[,-1],2,get_bootstrapped_pval),
         Term=names(apply(effect_on_MF$Total[,-1],2,quantile,.025)),
         Response="MF"))%>%
  mutate(., Term=sapply(1:nrow(.),function(x){
    if (any(grep("LGM",.$Term[x]))){
      return(gsub(".L","_L",.$Term[x]))
    }else if (any(grep("Current",.$Term[x]))){
      return(gsub(".C","_C",.$Term[x]))
    }else if (any(grep("Holo",.$Term[x]))){
      return(gsub(".H","_H",.$Term[x]))
    }else{
      return(.$Term[x])
    }
  }))%>%
  mutate(., Variable=sapply(1:nrow(.),function(x){
    if (length(strsplit(.$Term[x],"_")[[1]])==2){
      return(strsplit(.$Term[x],"_")[[1]][1])
    }else {
      return(.$Term[x])
    }
  }))%>%
  mutate(., Class=sapply(1:nrow(.),function(x){
    if (length(strsplit(.$Term[x],"_")[[1]])==2){
      return(strsplit(.$Term[x],"_")[[1]][2])
    }else {
      return("Soil")
    }
  }))%>%
  mutate(., Variable=recode_factor(Variable,
                                   "MAT"="Mean temperature",
                                   "MAP"="Mean precipitation",
                                   "MF"="Multifunctionality",
                                   "Sand"="Sand soil content",
                                   "Desert_Current1"="Desert, Current",
                                   "Desert_LGM"="Desert, LGM",
  ))%>%
  mutate(., Class=recode_factor(Class,
                                "Current1"="Current",
                                "LGM1"="LGM"
  ))%>%
  mutate(., Class=recode_factor(Class,
                                "Holocene"="Holocene legacy",
                                "LGM"="LGM legacy",
                                "Current"= "Current climate"
  ))

Total_effects$Signif2=sapply(1:nrow(Total_effects),function(x){return(is_signif(Total_effects$pval[x]))})

for (k in 1:2){
  assign(paste0("p",k),ggplot(Total_effects%>%dplyr::filter(., Response %in% paste0("Struct",k))%>%
                                dplyr::arrange(.,Response,q2)%>%
                                add_column(.,Order=(1:(nrow(.))))%>%
                                dplyr::group_by(., Response)%>%
                                dplyr::mutate(.,Term=fct_reorder(Term,Order)))+
           geom_linerange(aes(x=Variable,ymin=q1,y=q2,ymax=q3,group=interaction(Variable,Class),color=Class,alpha=Signif2),
                          position = position_dodge2(width = .8, preserve = "single"))+
           geom_linerange(aes(x=Variable,ymin=q12,y=q2,ymax=q32,group=interaction(Variable,Class),color=Class,alpha=Signif2),
                          position = position_dodge2(width = .8, preserve = "single"),lwd=2)+
           geom_point(aes(x=Variable,y=q2,group=interaction(Variable,Class),alpha=Signif2),color="black",
                      position = position_dodge2(width = .8, preserve = "single"))+
           
           #facet_wrap(.~Response,scales="free")+
           the_theme2+
           guides(alpha="none")+
           facet_wrap(.~Response,scales="free")+
           # scale_fill_manual(values=c("Soil"="#A2D086","LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           scale_color_manual(values=c("Soil"="#A2D086","LGM legacy"="#FFE699","Holocene legacy"="#F3A875","Current climate"="#A97858"))+
           scale_alpha_manual(values=c(.2,1,1))+
           # scale_color_manual(values=c("grey","black"))+
           geom_hline(yintercept = 0)+
           labs(y=c("Total effect on PC1: \n mean patch-size and cover","Total effect on PC2: \n spatial aggregation")[k],fill="",color="",x="Predictor")+
           theme(strip.text.x.top = element_blank(),strip.text.y.top = element_blank())+
           coord_flip())
}

for (k in 1:2){
  assign(paste0("p",k+2),ggplot(Total_effects%>%dplyr::filter(., Response %in% c("Sand","MF")[k])%>%
                                  dplyr::arrange(.,Response,q2)%>%
                                  add_column(.,Order=(1:(nrow(.))))%>%
                                  dplyr::group_by(., Response)%>%
                                  dplyr::mutate(.,Term=fct_reorder(Term,Order)))+
           geom_linerange(aes(x=Variable,ymin=q1,y=q2,ymax=q3,group=interaction(Variable,Class),color=Class,alpha=Signif2),
                          position = position_dodge2(width = .8, preserve = "single"))+
           geom_linerange(aes(x=Variable,ymin=q12,y=q2,ymax=q32,group=interaction(Variable,Class),color=Class,alpha=Signif2),
                          position = position_dodge2(width = .8, preserve = "single"),lwd=2)+
           geom_point(aes(x=Variable,y=q2,group=interaction(Variable,Class),alpha=Signif2),color="black",
                      position = position_dodge2(width = .8, preserve = "single"))+
           
           #facet_wrap(.~Response,scales="free")+
           the_theme2+
           guides(alpha="none")+
           facet_wrap(.~Response,scales="free")+
           # scale_fill_manual(values=c("Soil"="#A2D086","LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           scale_color_manual(values=c("Soil"="#A2D086","LGM legacy"="#FFE699","Holocene legacy"="#F3A875","Current climate"="#A97858"))+
           scale_alpha_manual(values=c(.2,1,1))+
           # scale_color_manual(values=c("grey","black"))+
           geom_hline(yintercept = 0)+
           labs(y=c("Effect on sand soil content","Effect on multifunctionality")[k],fill="",color="",x="Predictor")+
           theme(strip.text.x.top = element_blank(),strip.text.y.top = element_blank())+
           coord_flip())
}



p_tot=ggarrange(p3+theme(axis.title.x = element_text(color="#4F6F3B")),
                p4+theme(axis.text.y = element_blank(),
                         axis.title.x = element_text(color="#4F6F3B"),
                         axis.title.y = element_blank(),
                         axis.ticks.y = element_blank()),
                widths = c(1,.75),legend = "none",common.legend = T,labels = letters[1:2],font.label = list(size=15),hjust = c(-7,-2))

p_tot2=ggarrange(p1+theme(axis.title.x = element_text(color="#9446CE")),
                 p2+theme(axis.text.y = element_blank(),
                          axis.title.x = element_text(color="#9446CE"),
                          axis.title.y = element_blank(),
                          axis.ticks.y = element_blank()),
                 widths = c(1,.75),legend = "bottom",common.legend = T,labels = letters[2:3],font.label = list(size=15),hjust = c(-7,-2))

ggsave("./Figures/SEM.pdf",
       ggarrange(ggplot()+theme_void(),#p_tot,
                 p_tot2,nrow=2,labels = c("a",""),heights = c(1,1.2),hjust = c(-7.5),
                 font.label = list(size=15)),width = 7,height = 7)
ggsave("./Figures/SI/SEM_sand_MF.pdf",p_tot,width = 7,height = 4)

## >> Some bivariate plots ----

d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.)%>%
  mutate(., Sand=bcPower(.$Sand,6))%>%
  mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
d=d%>%
  dplyr::mutate(.,
                Desert_Current=as.factor(as.numeric(Desert_Current=="Desert")),
                Desert_LGM=as.factor(as.numeric(Desert_LGM=="Desert")))%>%
  dplyr::select(.,-Biome_LGM,-Biome_Current)

library(piecewiseSEM)

d_site=dplyr::distinct(d,Site_ID,.keep_all = T)

mod_Struct1=lmer(Struct1~Isothermality_Holocene+
                   MAP_Holocene+MAP_Current+
                   MAP_LGM+MAT_LGM+Isothermality_Current+
                   MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
                   Isothermality_LGM+Sand+MF+
                   (1|Site_ID),data = d,na.action = na.fail,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_Struct2=lmer(Struct2~Isothermality_Holocene+
                   MAP_Holocene+MAP_Current+
                   MAP_LGM+MAT_LGM+Isothermality_Current+
                   MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
                   Isothermality_LGM+Sand+MF+
                   (1|Site_ID),data = d,na.action = na.fail,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

d_res=d_points=tibble()
for (k in c("Isothermality_Holocene",
            "Desert_LGM",
            "MAP_Holocene","MAP_Current",
            "MAP_LGM","MAT_LGM","Isothermality_Current",
            "MAT_Current","MAT_Holocene",
            "Isothermality_LGM")){
  
  
  res_reg1=visreg::visreg(mod_Struct1,k,plot = F)
  res_reg1$fit=res_reg1$fit%>%melt(., measure.vars=k)
  d_res=rbind(d_res,tibble(Upper=res_reg1$fit$visregUpr,Lower=res_reg1$fit$visregLwr,
                           Median=res_reg1$fit$visregFit,Driver=res_reg1$fit$value,
                           Response="Struct1",Fullname=k,
                           Climate_var=strsplit(k,"_")[[1]][1],
                           Period=strsplit(k,"_")[[1]][2])) 
  
  res_reg1$res=res_reg1$res%>%melt(., measure.vars=k)
  d_points=rbind(d_points,tibble(Points=res_reg1$res$visregRes,Driver=res_reg1$res$value,
                                 Response="Struct1",Fullname=k,
                                 Climate_var=strsplit(k,"_")[[1]][1],
                                 Period=strsplit(k,"_")[[1]][2])) 
  
  res_reg2=visreg::visreg(mod_Struct2,k,plot = F)
  res_reg2$fit=res_reg2$fit%>%melt(., measure.vars=k)
  d_res=rbind(d_res,tibble(Upper=res_reg2$fit$visregUpr,Lower=res_reg2$fit$visregLwr,
                           Median=res_reg2$fit$visregFit,Driver=res_reg2$fit$value,
                           Response="Struct2",Fullname=k,
                           Climate_var=strsplit(k,"_")[[1]][1],
                           Period=strsplit(k,"_")[[1]][2])) 
  
  res_reg2$res=res_reg2$res%>%melt(., measure.vars=k)
  d_points=rbind(d_points,tibble(Points=res_reg2$res$visregRes,Driver=res_reg2$res$value,
                                 Response="Struct2",Fullname=k,
                                 Climate_var=strsplit(k,"_")[[1]][1],
                                 Period=strsplit(k,"_")[[1]][2])) 
  
}
d_points$Driver=as.numeric(d_points$Driver)
d_res$Driver=as.numeric(d_res$Driver)

for (id in 1:4){
  assign(paste0("p",id),ggplot(d_res%>%
                                 dplyr::filter(., Response=="Struct1", Fullname %in% c("MAP_LGM",
                                                                                       "Isothermality_LGM",
                                                                                       "Isothermality_Holocene",
                                                                                       "MAT_LGM")[id]))+
           geom_point(data=d_points%>%dplyr::filter(., Response=="Struct1", Fullname %in% c("MAP_LGM",
                                                                                            "Isothermality_LGM",
                                                                                            "Isothermality_Holocene",
                                                                                            "MAT_LGM")[id]),
                      aes(x=Driver,y=Points),alpha=.4,color="grey")+
           geom_line(aes(x=Driver,y=Median,color=Period),lwd=1.5)+
           geom_ribbon(aes(x=Driver,y=Median,ymin=Lower,ymax=Upper,fill=Period),alpha=.4)+
           the_theme2+
           labs(x="",y="PC1: mean patch-size\n and cover",fill="")+
           scale_fill_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           scale_color_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           theme(legend.text = element_text(size=12),strip.placement = "outside")+
           guides(color="none",fill="none"))
  
}

p51=Increase_size_axes(ggplot(d_points%>%
            dplyr::filter(., Response=="Struct1", Fullname =="Desert_LGM"))+
  geom_boxplot(aes(x=Driver,y=Points,group=Driver),alpha=.4,color="#000",width = .2,outlier.shape = NA)+
  geom_jitter(aes(x=Driver,y=Points),alpha=.4,color="grey",width = .2)+
  the_theme2+
  labs(x="",y="PC1: mean patch-size\n and cover",fill="")+
  theme(legend.text = element_text(size=12),strip.placement = "outside")+
    scale_x_continuous(breaks = c(0,1),labels = c("No desert \n during LGM","Desert \n during LGM"))+
    guides(color="none",fill="none"))


p_tot1=ggarrange(Increase_size_axes(p1)+labs(x="Precipitation legacy LGM"),
          Increase_size_axes(p2)+theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.title.y = element_blank())+labs(x="Isothermality legacy LGM"),
          Increase_size_axes(p4)+theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.title.y = element_blank())+labs(x="Temperature legacy LGM"),
          Increase_size_axes(p3)+theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.title.y = element_blank())+labs(x="Isothermality legacy Holocene"),
          ncol=4,widths = c(1.25,1,1,1),hjust = -1)

for (id in 1:4){
  assign(paste0("p",id),ggplot(d_res%>%
                                 dplyr::filter(., Response=="Struct2", Fullname %in% c("Isothermality_Holocene",
                                                                                       "MAT_Holocene",
                                                                                       "Isothermality_LGM",
                                                                                       "MAT_LGM")[id]))+
           geom_point(data=d_points%>%dplyr::filter(., Response=="Struct2", Fullname %in% c("Isothermality_Holocene",
                                                                                            "MAT_Holocene",
                                                                                            "Isothermality_LGM",
                                                                                            "MAT_LGM")[id]),
                      aes(x=Driver,y=Points),alpha=.4,color="grey")+
           geom_line(aes(x=Driver,y=Median,color=Period),lwd=1.5)+
           geom_ribbon(aes(x=Driver,y=Median,ymin=Lower,ymax=Upper,fill=Period),alpha=.4)+
           the_theme2+
           labs(x="",y="PC2: spatial\n aggregation",fill="")+
           scale_fill_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           scale_color_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           theme(legend.text = element_text(size=12),strip.placement = "outside")+
           guides(color="none",fill="none"))
  
}

p52=Increase_size_axes(ggplot(d_points%>%
            dplyr::filter(., Response=="Struct2", Fullname =="Desert_LGM"))+
  geom_boxplot(aes(x=Driver,y=Points,group=Driver),alpha=.4,color="#000",width = .2,outlier.shape = NA)+
  geom_jitter(aes(x=Driver,y=Points),alpha=.4,color="grey",width = .2)+
  the_theme2+
  labs(x="",y="PC2: spatial\n aggregation",fill="")+
  theme(legend.text = element_text(size=12),strip.placement = "outside")+
  scale_x_continuous(breaks = c(0,1),labels = c("No desert \n during LGM","Desert \n during LGM"))+
  guides(color="none",fill="none"))

p_tot2=ggarrange(Increase_size_axes(p3)+labs(x="Isothermality legacy LGM"),
                 Increase_size_axes(p4)+theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.title.y = element_blank())+labs(x="Temperature legacy LGM"),
                 Increase_size_axes(p1)+theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.title.y = element_blank())+labs(x="Isothermality legacy Holocene"),
                 Increase_size_axes(p2)+theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.title.y = element_blank())+labs(x="Temperature legacy Holocene"),
          ncol=4,widths = c(1.25,1,1,1))

p_tot=ggarrange(p_tot1,p_tot2,nrow=2,labels = letters[1:2],hjust = c(-2,-1))

ggsave("./Figures/Partial_res.pdf",ggarrange(p_tot,
                                             ggarrange(ggplot()+theme_void(),p51+theme(axis.text.x = element_text(size=10)),p52+theme(axis.text.x = element_text(size=10)),
                                                       ggplot()+theme_void(),ncol=4,
                                                       widths = c(.6,1,1,.6),labels = c("c","","",""),hjust = -11,
                                                       font.label = list(size=15)),nrow=2,heights = c(1,.5)),width = 12,height = 8)



# ------------------------------------- SI analyses and figures ----

## >> Variance partitioning with env, biocom ----

id_metric_kept=c(1,2,3,8,12,14,18)

for (ID_PC in 1:2){
  
  d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
    Closer_to_normality(.)%>%
    Perform_PCA_spatial_struc(.)
  
  d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
  d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
  
  d=d%>%
    dplyr::mutate(.,
                  Desert_Current=as.numeric(Desert_Current=="Desert"),
                  Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
    dplyr::select(.,-Biome_LGM,-Biome_Current)
  
  d$Struct1=d[,paste0("Struct",ID_PC)]
  
  #LGM clim
  d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))])]
  #Holocene clim
  d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))])]
  
  #Current clim
  d_Current=d[,c(colnames(d)[grep("Current",colnames(d))])]
  
  #Environmental var
  d_Env=d[,c("Sand","MF","Slope","Elevation")]
  
  
  assign(paste0("d_",ID_PC),
         Plot_variance_partition(varpart(d[,c("Struct1")],d_LGM,d_Holo,d_Current,d_Env),bar = T))
  
  
  combined_d=cbind(cbind(d_LGM),
                   d_Holo,
                   cbind(d_Current),
                   d_Env,
                   d$Struct1)
  colnames(combined_d)[ncol(combined_d)]="Struct1"
  train_id=sample(1:nrow(combined_d),round(.7*nrow(combined_d)))
  
  # RDA_mod_RF= rda(combined_d$Struct1 ~ ., data = combined_d)
  # Select_var_RDA_Env = Select_variable_VarPart(RDA_all = RDA_mod_RF,
  #                                              data = combined_d,
  #                                              response_var = d$Struct1,
  #                                              forwardsel = F)
  
  #Changing to factors
  if (any(colnames(combined_d)=="Desert_Current")){
    combined_d[,"Desert_Current"]=as.factor(combined_d[,"Desert_Current"])
  }
  
  if (any(colnames(combined_d)=="Desert_LGM")){
    combined_d[,"Desert_LGM"]=as.factor(combined_d[,"Desert_LGM"])
  }
  
  train_id=sample(1:nrow(d),round(.7*nrow(d)))
  
  BRT_model = gbm(Struct1 ~ .,
                  distribution = "gaussian",
                  data = combined_d[train_id,],
                  interaction.depth = 7,
                  n.trees = 500,
                  cv.folds = 10,
                  n.minobsinnode = 4, 
                  shrinkage = .1,
                  train.fraction = .7,
                  n.cores = 1,
                  bag.fraction = .8)
  
  
  assign(paste0("pRF_",ID_PC),
         ggplot(as.data.frame(summary(BRT_model,plot=F))%>%
                  add_column(., mean_imp=.$rel.inf)%>%
                  Change_name_RF(.,"var")%>%
                  dplyr::arrange(., mean_imp)%>%
                  add_column(., Order_stat=1:nrow(.))%>%
                  dplyr::filter(., Order_stat>(max(Order_stat)-12))%>%
                  mutate(.,Name_stat = fct_reorder(Name_stat, Order_stat)))+
           geom_bar(aes(y=mean_imp,group=interaction(Name_stat,Group_variable),x=Name_stat,
                        fill=Group_variable),stat="identity")+
           geom_linerange(aes(x=Name_stat,ymin=ifelse(mean_imp-sd_imp<0,0,mean_imp-sd_imp),
                              ymax=mean_imp+sd_imp,group=interaction(Name_stat,Group_variable),
                              fill=Group_variable),shape=21,size=1)+
           
           coord_flip()+
           the_theme2+
           scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Current"="pink","Env"="#C1E0B8",
                                      "Shared past clim."="#DAB2BA","Other shared variance"="#E8E8E8"))+
           labs(y="Mean Importance",x="",fill=""))
  
}  

pvariance=ggarrange(d_1$Figure+ggtitle("PC1: mean patch-size, cover")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "none"),
                    d_2$Figure+ggtitle("PC2: spatial aggregation")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "none"),
                    ncol=2,common.legend = T,legend="none")

PCA_stat=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.,return_PCA = T)

p0=plotPCA(PCA_stat,1,2,F,xLim = c(-3,3),yLim = c(-3,3),
           colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")


p_legend=get_legend(ggplot(tibble(ID=1:6,Name=c("LGM legacy","Holocene legacy","Env",
                                                "Current climate","Shared past clim.",
                                                "Other shared variance")))+
                      geom_bar(aes(x=1,y=1,fill=as.factor(Name)),stat="identity")+
                      scale_fill_manual(values=c("LGM legacy"="#FFF9BF",
                                                 "Holocene legacy"="#FFD2A0",
                                                 "Env"="#C1E0B8",
                                                 "Current climate"="pink",
                                                 "Shared past clim."="#DAB2BA",
                                                 "Other shared variance"="#E8E8E8"))+
                      labs(fill="")+theme(legend.position = "bottom"))

p_tot=ggarrange(
  ggarrange(p0+the_theme2+theme(legend.position = "none"),
            pvariance,ncol=2,widths = c(.7,1),labels = letters[1:2],align = "hv"),p_legend,nrow=2,heights = c(1,.15),align = "hv")

ggsave("./Figures/SI/Figure_irregular_with_env.pdf",p_tot,width = 9,height = 4.5)

p_tot_RF=ggarrange(p_RF,p_legend,nrow=2,heights = c(1,.15),align = "hv")

ggsave("./Figures/SI/Figure_RF_irregular_with_env.pdf",p_tot_RF,width = 8,height = 4.5)

## >> Metric by metric analysis ----

list_spatial_metric=c("Small_patches","moran_I","mean_psd","fmax_psd","flow_length","core_area","cv_psd")
Name_spatial=c("# small patch","Moran I","Mean PSD","Largest patch","Bare connectivity","Shape patch","CV PSD")

index=1
id_metric_kept=c(1,2,3,8,12,14,18)

for (metric_id in list_spatial_metric){
  
  d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
    Closer_to_normality(.)%>%
    Perform_PCA_spatial_struc(.)
  
  d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
  d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
  
  d=d%>%
    dplyr::mutate(.,
                  Desert_Current=as.numeric(Desert_Current=="Desert"),
                  Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
    dplyr::select(.,-Biome_LGM,-Biome_Current)
  
  d$Struct1=d[,metric_id]
  d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))[c(id_metric_kept,20)]])]
  #Holocene clim
  d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
  
  #Current clim
  d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept,20)]])]
  
  assign(paste0("d_",index),
         Plot_variance_partition(varpart(d[,c("Struct1")],d_LGM,d_Holo,d_Current),
                                 bar = T,four_groups = F))
  
  combined_d=cbind(cbind(d_LGM),
                   d_Holo,
                   cbind(d_Current),
                   d$Struct1)
  colnames(combined_d)[ncol(combined_d)]="Struct1"
  train_id=sample(1:nrow(combined_d),round(.7*nrow(combined_d)))
  
  #Changing to factors
  if (any(colnames(combined_d)=="Desert_Current")){
    combined_d[,"Desert_Current"]=as.factor(combined_d[,"Desert_Current"])
  }
  
  if (any(colnames(combined_d)=="Desert_LGM")){
    combined_d[,"Desert_LGM"]=as.factor(combined_d[,"Desert_LGM"])
  }
  
  
  combined_d=cbind(cbind(Var_LGM),
                   Var_Holo,
                   cbind(Var_Current),
                   # Var_Env,
                   d$Struct1)
  colnames(combined_d)[ncol(combined_d)]="Struct1"
  train_id=sample(1:nrow(combined_d),round(.7*nrow(combined_d)))
  
  if (any(colnames(combined_d)=="Desert_Current")){
    combined_d[,"Desert_Current"]=as.factor(combined_d[,"Desert_Current"])
  }
  
  if (any(colnames(combined_d)=="Desert_LGM")){
    combined_d[,"Desert_LGM"]=as.factor(combined_d[,"Desert_LGM"])
  }
  
  set.seed(123)
  train_id=sample(1:nrow(d),round(.7*nrow(d)))
  
  BRT_model = gbm(Struct1 ~ .,
                  distribution = "gaussian",
                  data = combined_d[train_id,],
                  interaction.depth = 7,
                  n.trees = 500,
                  cv.folds = 10,
                  n.minobsinnode = 4, 
                  shrinkage = .1,
                  train.fraction = .7,
                  n.cores = 1,
                  bag.fraction = .8)
  
  assign(paste0("pRF_",index),
         ggplot(as.data.frame(summary(BRT_model,plot=F))%>%
                  add_column(., mean_imp=.$rel.inf)%>%
                  Change_name_RF(.,"var")%>%
                  dplyr::arrange(., mean_imp)%>%
                  add_column(., Order_stat=1:nrow(.))%>%
                  dplyr::filter(., Order_stat>(max(Order_stat)-12))%>%
                  mutate(.,Name_stat = fct_reorder(Name_stat, Order_stat)))+
           geom_bar(aes(y=mean_imp,group=interaction(Name_stat,Group_variable),x=Name_stat,
                        fill=Group_variable),stat="identity")+
           
           coord_flip()+
           the_theme2+
           scale_fill_manual(values=c("LGM legacy"="#FFF9BF", "Holocene legacy"="#FFD2A0","Current climate"="pink",
                                      "Shared past clim."="#DAB2BA","Other shared variance"="#E8E8E8"))+
           labs(y="Mean Importance",x="",fill=""))
  index=index+1
}  

PCA_stat=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.,F,T)

p0=plotPCA(PCA_stat,1,2,F,xLim = c(-3,3),yLim = c(-3,3),colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")

pvariance=ggarrange(d_1$Figure+ggtitle("# small patch")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_2$Figure+ggtitle("Moran I")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_3$Figure+ggtitle("Mean PSD")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_4$Figure+ggtitle("Largest patch")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_5$Figure+ggtitle("Bare connectivity")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_6$Figure+ggtitle("Shape patch")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_7$Figure+ggtitle("CV PSD")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    p0,
                    ncol=4,nrow=2,common.legend = T,legend="bottom")

ggsave("./Figures/SI/All_variances_individualstats.pdf",pvariance,width = 10,height = 7)

p_RF=ggarrange(
  pRF_1+ggtitle("# small patch")+theme(legend.position = "none"),
  pRF_2+ggtitle("Moran I")+theme(legend.position = "none"),
  pRF_3+ggtitle("Mean PSD")+theme(legend.position = "none"),
  pRF_4+ggtitle("Largest patch")+theme(legend.position = "none"),
  pRF_5+ggtitle("Bare connectivity")+theme(legend.position = "none"),
  pRF_6+ggtitle("Shape patch")+theme(legend.position = "none"),
  pRF_7+ggtitle("CV PSD")+theme(legend.position = "none"),
  ncol=4,nrow=2,common.legend = T,legend = "bottom")

ggsave("./Figures/SI/All_RF_individualstats.pdf",p_RF,width = 14,height = 7)

## >> Histogram image color ----

pdf("./Figures/SI/Distribution_color_image.pdf",width = 15,height = 6)
par(mfrow=c(1,3))
k="00-sentinel2-11-58n-27-94e-sudan-2020-09-28-19-33-14_PROCESSED_RGB"
list_images_years=list.files(paste0("./Data/Regular_images_Buxton/Unzipped/",k))
img=imager::load.image(paste0("./Data/Regular_images_Buxton/Unzipped/",k,"/",list_images_years[length(list_images_years)])) #last image
plot(img)
kmean_img=k_means_RGB(img,2)
cats=get_cut_grayscale_values(2)[[1]]
mat=kmean_img %>% binarize(cats[[1]], cats[[2]])
image(mat,col=c("white","black"))
hist(img,main="Distribution of pixel image color",xlab="Image color RGB",ylab="Number of pixels")
dev.off()

## >> Some images ----


pdf("./Figures/SI/Some_regular_images.pdf",width = 15,height = 6)
par(mfrow=c(1,3))
k=c(2,31,32)
list_folder=list.files(paste0("./Data/Regular_images_Buxton/Unzipped/"))

list_image=list.files(paste0("./Data/Regular_images_Buxton/Unzipped/",list_folder[k[1]]))
img=imager::load.image(paste0("./Data/Regular_images_Buxton/Unzipped/",list_folder[k[1]],
                              "/",list_image[length(list_image)])) #last image
plot(img)

list_image=list.files(paste0("./Data/Regular_images_Buxton/Unzipped/",list_folder[k[2]]))
img=imager::load.image(paste0("./Data/Regular_images_Buxton/Unzipped/",list_folder[k[2]],
                              "/",list_image[length(list_image)])) #last image
plot(img)

list_image=list.files(paste0("./Data/Regular_images_Buxton/Unzipped/",list_folder[k[3]]))
img=imager::load.image(paste0("./Data/Regular_images_Buxton/Unzipped/",list_folder[k[3]],
                              "/",list_image[length(list_image)])) #last image
plot(img)
dev.off()


## >> Correlation climatic variables ----

id_metric_kept=c(1,2,3,8,12,14,18)

d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.,F)

d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")

d=d%>%
  dplyr::mutate(.,
                Desert_Current=as.numeric(Desert_Current=="Desert"),
                Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
  dplyr::select(.,-Biome_LGM,-Biome_Current)

d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))[c(id_metric_kept)]])]
d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept)]])]

p1=Plot_correlation_variables(d_LGM)
p2=Plot_correlation_variables(d_Holo)
p3=Plot_correlation_variables(d_Current)

p_tot=ggarrange(p1+ggtitle("LGM climatic conditions"),
          p2+ggtitle("Mid-Holocene climatic conditions"),
          p3+ggtitle("Current climatic conditions"),nrow = 3,labels = letters[1:3])

ggsave("./Figures/SI/Correlation_climate.pdf",p_tot,width = 5,height = 15)

## >> SEM indirect effects with desert factors----

d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.)%>%
  mutate(., Sand=bcPower(.$Sand,6))%>%
  mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
d=d%>%
  dplyr::mutate(.,
                Desert_Current=as.factor(as.numeric(Desert_Current=="Desert")),
                Desert_LGM=as.factor(as.numeric(Desert_LGM=="Desert")))%>%
  dplyr::select(.,-Biome_LGM,-Biome_Current)

library(piecewiseSEM)
Plot_correlation_variables(d[,c("Isothermality_Holocene","MAP_Holocene","MAP_Current","MAP_LGM","MAT_LGM",#MAT_Holocene+MAP_Seasonality_LGM+MAP_Seasonality_Holocene+
                                "Isothermality_Current","MAT_Current","MAT_Holocene","Isothermality_LGM", 
                                "Sand","MF")])

d_site=dplyr::distinct(d,Site_ID,.keep_all = T)

mod_sand=lm(MF~Isothermality_Holocene+
              MAP_Holocene+MAP_Current+
              MAP_LGM+MAT_LGM+Isothermality_Current+
              MAT_Current+MAT_Holocene+
              Desert_LGM+Desert_Current+
              Isothermality_LGM,
            data = d_site,na.action = na.fail)

mod_MF=lm(Sand~Isothermality_Holocene+
            MAP_Holocene+MAP_Current+
            MAP_LGM+MAT_LGM+Isothermality_Current+
            MAT_Current+MAT_Holocene+
            Desert_LGM+Desert_Current+
            Isothermality_LGM,
          data = d_site,na.action = na.fail)

mod_Struct1=lmer(Struct1~Isothermality_Holocene+
                   MAP_Holocene+MAP_Current+
                   MAP_LGM+MAT_LGM+Isothermality_Current+
                   MAT_Current+MAT_Holocene+
                   Desert_LGM+Desert_Current+
                   Isothermality_LGM+Sand+MF+
                   (1|Site_ID),data = d,na.action = na.fail,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_Struct2=lmer(Struct2~Isothermality_Holocene+
                   MAP_Holocene+MAP_Current+
                   MAP_LGM+MAT_LGM+Isothermality_Current+
                   MAT_Current+MAT_Holocene+
                   Desert_LGM+Desert_Current+
                   Isothermality_LGM+Sand+MF+
                   (1|Site_ID),data = d,na.action = na.fail,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


vif(mod_sand)
vif(mod_MF)
vif(mod_Struct1)
vif(mod_Struct2)

# confint(mod_sand)
# confint(mod_MF)
# confint(mod_Struct1)
# confint(mod_Struct2)
# confint(mod_sand2)
# confint(mod_MF2)
# confint(mod_Struct12)
# confint(mod_Struct22)




d_SEM=psem(
  mod_sand,
  mod_MF,
  mod_Struct1,
  mod_Struct2,
  MF%~~%Sand
)

piecewiseSEM::fisherC(d_SEM)
piecewiseSEM::dSep(d_SEM)


SEM_boot = bootEff(d_SEM, R = 1000, seed = 13, parallel = "snow",ran.eff = "Site_ID",ncpus = 40)
saveRDS(SEM_boot,"./Data/SEM_boot.rds")


SEM_boot=readRDS("./Data/SEM_boot.rds")
Effects = semEff(SEM_boot)

#Getting information about indirect effects through direct and total effects
#Aggregating in a tibble
get_bootstrapped_pval=function(x){
  return(ifelse(length(which(x>0))/length(x)>.5,length(which(x<0))/length(x),length(which(x>0))/length(x)))
}

Direct_effects=tibble()
for (predictors in c("Struct1","Struct2","MF","Sand")){
  effect_pred=Effects$Effects$Bootstrapped[[predictors]]
  Direct_effects=rbind(Direct_effects,
                       tibble(q1=apply(effect_pred$Direct[,-1],2,quantile,.025),
                              q3=apply(effect_pred$Direct[,-1],2,quantile,.975),
                              q2=apply(effect_pred$Direct[,-1],2,median),
                              pval=apply(effect_pred$Direct[,-1],2,get_bootstrapped_pval),
                              Term=names(apply(effect_pred$Direct[,-1],2,quantile,.025)),
                              Response=predictors))
}

effect_on_struct1=Effects$Effects$Bootstrapped$Struct1
effect_on_struct2=Effects$Effects$Bootstrapped$Struct2
effect_on_sand=Effects$Effects$Bootstrapped$Sand
effect_on_MF=Effects$Effects$Bootstrapped$MF

Total_effects=rbind(
  tibble(q1=apply(effect_on_struct1$Total[,-1],2,quantile,.025),
         q3=apply(effect_on_struct1$Total[,-1],2,quantile,.975),
         q2=apply(effect_on_struct1$Total[,-1],2,median),
         pval=apply(effect_on_struct1$Total[,-1],2,get_bootstrapped_pval),
         Term=names(apply(effect_on_struct1$Total[,-1],2,quantile,.025)),
         Response="Struct1"),
  tibble(q1=apply(effect_on_struct2$Total[,-1],2,quantile,.025),
         q3=apply(effect_on_struct2$Total[,-1],2,quantile,.975),
         q2=apply(effect_on_struct2$Total[,-1],2,median),
         pval=apply(effect_on_struct2$Total[,-1],2,get_bootstrapped_pval),
         Term=names(apply(effect_on_struct2$Total[,-1],2,quantile,.025)),
         Response="Struct2"),
  tibble(q1=apply(effect_on_sand$Total[,-1],2,quantile,.025),
         q3=apply(effect_on_sand$Total[,-1],2,quantile,.975),
         q2=apply(effect_on_sand$Total[,-1],2,median),
         pval=apply(effect_on_sand$Total[,-1],2,get_bootstrapped_pval),
         Term=names(apply(effect_on_sand$Total[,-1],2,quantile,.025)),
         Response="Sand"),
  tibble(q1=apply(effect_on_MF$Total[,-1],2,quantile,.025),
         q3=apply(effect_on_MF$Total[,-1],2,quantile,.975),
         q2=apply(effect_on_MF$Total[,-1],2,median),
         pval=apply(effect_on_MF$Total[,-1],2,get_bootstrapped_pval),
         Term=names(apply(effect_on_MF$Total[,-1],2,quantile,.025)),
         Response="MF"))


ggplot(Total_effects%>%dplyr::filter(., Response %in% paste0("Struct",1:2)))+
  geom_pointrange(aes(x=Term,ymin=q1,y=q2,ymax=q3,group=Response,color=Response),
                  position=position_jitterdodge(seed=123))+
  #facet_wrap(.~Response,scales="free")+
  the_theme2+
  geom_hline(yintercept = 0)+
  theme(axis.text.x = element_text(angle = 60,hjust=1))+
  coord_flip()

ggplot(Total_effects%>%dplyr::filter(., Response %!in% paste0("Struct",1:2)))+
  geom_pointrange(aes(x=Term,ymin=q1,y=q2,ymax=q3,group=Response,color=Response),
                  position=position_jitterdodge(seed=123))+
  #facet_wrap(.~Response,scales="free")+
  the_theme2+
  geom_hline(yintercept = 0)+
  theme(axis.text.x = element_text(angle = 60,hjust=1))+
  coord_flip()



## >> PCA clim for SEM ----

id_metric_kept=c(1,2,3,8,12,14,18)
d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.,F)

d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")

d=d%>%
  dplyr::mutate(.,
                Desert_Current=as.numeric(Desert_Current=="Desert"),
                Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
  dplyr::select(.,-Biome_LGM,-Biome_Current)

#LGM clim
d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))[c(id_metric_kept)]])]
#Holocene clim
d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
#Current clim
d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept)]])]

p0=plotPCA(principal(d_LGM,nfactors = 2),1,2,F,
           xLim = c(-3,3),yLim = c(-3,3),
           colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")

p1=plotPCA(principal(d_Holo,nfactors = 2),1,2,F,
           xLim = c(-3,3),yLim = c(-3,3),
           colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")

p2=plotPCA(principal(d_Current,nfactors = 2),1,2,F,
           xLim = c(-3,3),yLim = c(-3,3),
           colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")


p_tot=ggarrange(p0+ggtitle("LGM period"),
          p1+ggtitle("Holocene period"),
          p2+ggtitle("Present-day"),ncol=3)

ggsave("./Figures/SI/PCA_CLIM_periods.pdf",p_tot,width = 12,height = 4)


## >> RF all variables ----

list_spatial_metric=c("Small_patches","moran_I","mean_psd","fmax_psd","flow_length","core_area","cv_psd")
Name_spatial=c("# small patch","Moran I","Mean PSD","Largest patch","Bare connectivity","Shape patch","CV PSD")
id_metric_kept=c(1,2,3,8,12,14,18)


R2=tibble()
index=1
for (metric_id in list_spatial_metric){
  
  d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
    Closer_to_normality(.)%>%
    Perform_PCA_spatial_struc(.,F)
  
  d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
  d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
  
  d=d%>%
    dplyr::mutate(.,
                  Desert_Current=as.numeric(Desert_Current=="Desert"),
                  Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
    dplyr::select(.,-Biome_LGM,-Biome_Current)
  
  d$Struct1=d[,metric_id]
  
  #LGM clim
  d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))[c(id_metric_kept,20)]])]
  #Holocene clim
  d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
  
  #Current clim
  d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept,20)]])]
  
  assign(paste0("d_",index),
         Plot_variance_partition(varpart(d[,c("Struct1")],d_LGM,d_Holo,d_Current),bar = T,four_groups = F))
  
  
  combined_d=cbind(cbind(d_LGM),
                   d_Holo,
                   cbind(d_Current),
                   d$Struct1)
  colnames(combined_d)[ncol(combined_d)]="Struct1"
  train_id=sample(1:nrow(combined_d),round(.7*nrow(combined_d)))
  
  if (any(colnames(combined_d)=="Desert_Current")){
    combined_d[,"Desert_Current"]=as.factor(combined_d[,"Desert_Current"])
  }
  
  if (any(colnames(combined_d)=="Desert_LGM")){
    combined_d[,"Desert_LGM"]=as.factor(combined_d[,"Desert_LGM"])
  }
  
  set.seed(123)
  train_id=sample(1:nrow(d),round(.7*nrow(d)))
  
  BRT_model = gbm(Struct1 ~ .,
                  distribution = "gaussian",
                  data = combined_d[train_id,],
                  interaction.depth = 7,
                  n.trees = 500,
                  cv.folds = 10,
                  n.minobsinnode = 4, 
                  shrinkage = .1,
                  train.fraction = .7,
                  n.cores = 1,
                  bag.fraction = .8)
  
  R2=rbind(R2,tibble(Stat=metric_id,R2=cor(predict(BRT_model,combined_d[-train_id,]),combined_d$Struct1[-train_id])))
  
  assign(paste0("pRF_",index),
         ggplot(as.data.frame(summary(BRT_model,plot=F))%>%
                  add_column(., mean_imp=.$rel.inf)%>%
                  Change_name_RF(.,"var")%>%
                  dplyr::arrange(., mean_imp)%>%
                  add_column(., Order_stat=1:nrow(.))%>%
                  dplyr::filter(., Order_stat>(max(Order_stat)-12))%>%
                  mutate(.,Name_stat = fct_reorder(Name_stat, Order_stat)))+
           geom_bar(aes(y=mean_imp,group=interaction(Name_stat,Group_variable),x=Name_stat,
                        fill=Group_variable),stat="identity")+
           
           coord_flip()+
           the_theme2+
           scale_fill_manual(values=c("LGM legacy"="#FFF9BF", "Holocene legacy"="#FFD2A0","Current climate"="pink",
                                      "Shared past clim."="#DAB2BA","Other shared variance"="#E8E8E8"))+
           labs(y="Mean Importance",x="",fill=""))
  index=index+1
}  

PCA_stat=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.,F,T)

p0=plotPCA(PCA_stat,1,2,F,xLim = c(-3,3),yLim = c(-3,3),
           colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")

pvariance=ggarrange(d_1$Figure+ggtitle("# small patch")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_2$Figure+ggtitle("Moran I")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_3$Figure+ggtitle("Mean PSD")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_4$Figure+ggtitle("Largest patch")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_5$Figure+ggtitle("Bare connectivity")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_6$Figure+ggtitle("Shape patch")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    d_7$Figure+ggtitle("CV PSD")+
                      labs(y="% of shared variance")+
                      theme(legend.position = "none"),
                    p0,
                    ncol=4,nrow=2,common.legend = T,legend="bottom")

ggsave("./Figures/SI/All_variances_individualstats.pdf",pvariance,width = 10,height = 7)

p_RF=ggarrange(
  pRF_1+ggtitle("# small patch, R2 = 0.69")+theme(legend.position = "none"),
  pRF_2+ggtitle("Moran I, R2 = 0.84")+theme(legend.position = "none"),
  pRF_3+ggtitle("Mean PSD, R2 = 0.76")+theme(legend.position = "none"),
  pRF_4+ggtitle("Largest patch, R2 = 0.72")+theme(legend.position = "none"),
  pRF_5+ggtitle("Bare connectivity, R2 = 0.77")+theme(legend.position = "none"),
  pRF_6+ggtitle("Shape patch, R2 = 0.69")+theme(legend.position = "none"),
  pRF_7+ggtitle("CV PSD, R2 = 0.76")+theme(legend.position = "none"),
  ncol=2,nrow=4,common.legend = T,legend = "bottom")

ggsave("./Figures/SI/All_RF_individualstats.pdf",p_RF,width = 7,height = 14)


## >> Some bivariate plots soil ----

d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.)%>%
  mutate(., Sand=bcPower(.$Sand,6))%>%
  mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
d=d%>%
  dplyr::mutate(.,
                Desert_Current=as.factor(as.numeric(Desert_Current=="Desert")),
                Desert_LGM=as.factor(as.numeric(Desert_LGM=="Desert")))%>%
  dplyr::select(.,-Biome_LGM,-Biome_Current)

library(piecewiseSEM)

d_site=dplyr::distinct(d,Site_ID,.keep_all = T)

mod_sand=lm(MF~Isothermality_Holocene+
              MAP_Holocene+MAP_Current+
              MAP_LGM+MAT_LGM+Isothermality_Current+
              MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
              Isothermality_LGM,
            data = d_site,na.action = na.fail)

mod_MF=lm(Sand~Isothermality_Holocene+
            MAP_Holocene+MAP_Current+
            MAP_LGM+MAT_LGM+Isothermality_Current+
            MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
            Isothermality_LGM,
          data = d_site,na.action = na.fail)

mod_Struct1=lmer(Struct1~Isothermality_Holocene+
                   MAP_Holocene+MAP_Current+
                   MAP_LGM+MAT_LGM+Isothermality_Current+
                   MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
                   Isothermality_LGM+Sand+MF+
                   (1|Site_ID),data = d,na.action = na.fail,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

mod_Struct2=lmer(Struct2~Isothermality_Holocene+
                   MAP_Holocene+MAP_Current+
                   MAP_LGM+MAT_LGM+Isothermality_Current+
                   MAT_Current+MAT_Holocene+Desert_Current+Desert_LGM+
                   Isothermality_LGM+Sand+MF+
                   (1|Site_ID),data = d,na.action = na.fail,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

d_res=d_points=tibble()
for (k in c("Isothermality_Holocene",
            "Desert_LGM",
            "MAP_Holocene","MAP_Current",
            "MAP_LGM","MAT_LGM","Isothermality_Current",
            "MAT_Current","MAT_Holocene",
            "Isothermality_LGM")){
  
  
  res_reg1=visreg::visreg(mod_sand,k,plot = F)
  res_reg1$fit=res_reg1$fit%>%melt(., measure.vars=k)
  d_res=rbind(d_res,tibble(Upper=res_reg1$fit$visregUpr,Lower=res_reg1$fit$visregLwr,
                           Median=res_reg1$fit$visregFit,Driver=res_reg1$fit$value,
                           Response="Sand",Fullname=k,
                           Climate_var=strsplit(k,"_")[[1]][1],
                           Period=strsplit(k,"_")[[1]][2])) 
  
  res_reg1$res=res_reg1$res%>%melt(., measure.vars=k)
  d_points=rbind(d_points,tibble(Points=res_reg1$res$visregRes,Driver=res_reg1$res$value,
                                 Response="Sand",Fullname=k,
                                 Climate_var=strsplit(k,"_")[[1]][1],
                                 Period=strsplit(k,"_")[[1]][2])) 
  
  res_reg2=visreg::visreg(mod_MF,k,plot = F)
  res_reg2$fit=res_reg2$fit%>%melt(., measure.vars=k)
  d_res=rbind(d_res,tibble(Upper=res_reg2$fit$visregUpr,Lower=res_reg2$fit$visregLwr,
                           Median=res_reg2$fit$visregFit,Driver=res_reg2$fit$value,
                           Response="MF",Fullname=k,
                           Climate_var=strsplit(k,"_")[[1]][1],
                           Period=strsplit(k,"_")[[1]][2])) 
  
  res_reg2$res=res_reg2$res%>%melt(., measure.vars=k)
  d_points=rbind(d_points,tibble(Points=res_reg2$res$visregRes,Driver=res_reg2$res$value,
                                 Response="MF",Fullname=k,
                                 Climate_var=strsplit(k,"_")[[1]][1],
                                 Period=strsplit(k,"_")[[1]][2])) 
  
}
d_points$Driver=as.numeric(d_points$Driver)
d_res$Driver=as.numeric(d_res$Driver)

for (id in 1:4){
  assign(paste0("p",id),ggplot(d_res%>%
                                 dplyr::filter(., Response=="Struct1", Fullname %in% c("MAP_LGM",
                                                                                       "Isothermality_LGM",
                                                                                       "Isothermality_Holocene",
                                                                                       "MAT_LGM")[id]))+
           geom_point(data=d_points%>%dplyr::filter(., Response=="Struct1", Fullname %in% c("MAP_LGM",
                                                                                            "Isothermality_LGM",
                                                                                            "Isothermality_Holocene",
                                                                                            "MAT_LGM")[id]),
                      aes(x=Driver,y=Points),alpha=.4,color="grey")+
           geom_line(aes(x=Driver,y=Median,color=Period),lwd=1.5)+
           geom_ribbon(aes(x=Driver,y=Median,ymin=Lower,ymax=Upper,fill=Period),alpha=.4)+
           the_theme2+
           labs(x="",y="PC1: mean patch-size and cover",fill="")+
           scale_fill_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           scale_color_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           theme(legend.text = element_text(size=12),strip.placement = "outside")+
           guides(color="none",fill="none"))
  
}

p51=ggplot(d_points%>%
             dplyr::filter(., Response=="Struct1", Fullname =="Desert_LGM"))+
  geom_boxplot(aes(x=Driver,y=Points,group=Driver),alpha=.4,color="#000",width = .2,outlier.shape = NA)+
  geom_jitter(aes(x=Driver,y=Points),alpha=.4,color="grey",width = .2)+
  the_theme2+
  labs(x="",y="PC1: mean patch-size and cover",fill="")+
  theme(legend.text = element_text(size=12),strip.placement = "outside")+
  scale_x_continuous(breaks = c(0,1),labels = c("No desert during LGM","Desert during LGM"))+
  guides(color="none",fill="none")


p_tot1=ggarrange(p1+labs(x="Mean precipitation, LGM"),
                 p2+theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+labs(x="Isothermality, LGM"),
                 p4+theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+labs(x="Mean Temperature, LGM"),
                 p3+theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+labs(x="Isothermality, Holocene"),
                 ncol=4,widths = c(1.25,1,1,1),hjust = -1)

for (id in 1:4){
  assign(paste0("p",id),ggplot(d_res%>%
                                 dplyr::filter(., Response=="Struct2", Fullname %in% c("Isothermality_Holocene",
                                                                                       "MAT_Holocene",
                                                                                       "Isothermality_LGM",
                                                                                       "MAT_LGM")[id]))+
           geom_point(data=d_points%>%dplyr::filter(., Response=="Struct2", Fullname %in% c("Isothermality_Holocene",
                                                                                            "MAT_Holocene",
                                                                                            "Isothermality_LGM",
                                                                                            "MAT_LGM")[id]),
                      aes(x=Driver,y=Points),alpha=.4,color="grey")+
           geom_line(aes(x=Driver,y=Median,color=Period),lwd=1.5)+
           geom_ribbon(aes(x=Driver,y=Median,ymin=Lower,ymax=Upper,fill=Period),alpha=.4)+
           the_theme2+
           labs(x="",y="PC2: spatial aggregation",fill="")+
           scale_fill_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           scale_color_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
           theme(legend.text = element_text(size=12),strip.placement = "outside")+
           guides(color="none",fill="none"))
  
}

p52=ggplot(d_points%>%
             dplyr::filter(., Response=="Struct2", Fullname =="Desert_LGM"))+
  geom_boxplot(aes(x=Driver,y=Points,group=Driver),alpha=.4,color="#000",width = .2,outlier.shape = NA)+
  geom_jitter(aes(x=Driver,y=Points),alpha=.4,color="grey",width = .2)+
  the_theme2+
  labs(x="",y="PC2: spatial aggregation",fill="")+
  theme(legend.text = element_text(size=12),strip.placement = "outside")+
  scale_x_continuous(breaks = c(0,1),labels = c("No desert during LGM","Desert during LGM"))+
  guides(color="none",fill="none")

p_tot2=ggarrange(p3+labs(x="Isothermality, LGM"),
                 p4+theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+labs(x="Mean Temperature, LGM"),
                 p1+theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+labs(x="Isothermality, Holocene"),
                 p2+theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+labs(x="Mean Temperature, Holocene"),
                 ncol=4,widths = c(1.25,1,1,1))

p_tot=ggarrange(p_tot1,p_tot2,nrow=2,labels = letters[1:2],hjust = c(-2,-1))

ggsave("./Figures/Partial_res.pdf",ggarrange(p_tot,
                                             ggarrange(ggplot()+theme_void(),p51+theme(axis.text.x = element_text(size=10)),p52+theme(axis.text.x = element_text(size=10)),
                                                       ggplot()+theme_void(),ncol=4,
                                                       widths = c(.6,1,1,.6),labels = c("c","","",""),hjust = -11,
                                                       font.label = list(size=15)),nrow=2,heights = c(1,.5)),width = 12,height = 8)




## >> Distribution legacy variables -----

d=read.table("./Data/data_sites_CLIM.csv",sep=";")

p=ggplot(d%>%melt(., measure.vars=c("MAT_LGM","MAT_Current","MAT_Holocene",
                                    "MAP_Holocene","MAP_LGM","MAP_Current",
                                    "Isothermality_Current","Isothermality_Holocene",
                                    "Isothermality_LGM"))%>%
           mutate(., Period=sapply(1:nrow(.),function(x){return(str_split(.$variable[x],"_")[[1]][2])})))+
  geom_histogram(aes(x=value,fill=Period))+
  the_theme2+
  facet_wrap(.~variable,scales="free")+
  scale_fill_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
  labs(x="Value",y="Number of sites",fill="")

ggsave("./Figures/SI/Distribution_climate_metrics.pdf",p,width = 7,height = 6)


d=read.table("./Data/data_sites_CLIM_nolegacy.csv",sep=";")


p=ggplot(d%>%melt(., measure.vars=c("MAT_LGM","MAT_Current","MAT_Holocene",
                                    "MAP_Holocene","MAP_LGM","MAP_Current",
                                    "Isothermality_Current","Isothermality_Holocene",
                                    "Isothermality_LGM"))%>%
           mutate(., Period=sapply(1:nrow(.),function(x){return(str_split(.$variable[x],"_")[[1]][2])})))+
  geom_histogram(aes(x=value,fill=Period))+
  the_theme2+
  facet_wrap(.~variable,scales="free")+
  scale_fill_manual(values=c("LGM"="#FFE699","Holocene"="#F3A875","Current"="#A97858"))+
  labs(x="Value",y="Number of sites",fill="")

ggsave("./Figures/SI/Distribution_climate_metrics_nolegacy.pdf",p,width = 7,height = 6)


