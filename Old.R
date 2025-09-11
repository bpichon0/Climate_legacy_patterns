## >> Mixed effect models ----

d=read.table("../Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality_spatialstructure(.)%>%
  Perform_PCA_spatial_struc(.)%>%
  mutate(., Isothermality_Current=bcPower(Isothermality_Current,-2),
         MAT_LGM=bcPower(MAT_LGM+20,-2))%>%
  mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))


d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")

d=d%>%
  dplyr::mutate(.,
                Desert_Current=as.numeric(Desert_Current=="Desert"),
                Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
  dplyr::select(.,-Biome_LGM,-Biome_Current)


Plot_correlation_variables(d[,c("MAT_Current","MAT_range_diurn_Current","MAP_Current",
                                "MAPdry_Current","Isothermality_Holocene",
                                "Desert_Current","Desert_LGM")])

d=d%>%mutate(., Desert_Current=as.factor(Desert_Current),Desert_LGM=as.factor(Desert_LGM))
mod_struct1=lmer(Struct1~Desert_LGM + Desert_Current +
                   MAT_Current+MAT_range_diurn_Current +MAP_Current+
                   MAPdry_Current+Isothermality_Holocene+(1|Site_ID),
                 data = d,na.action = na.fail)
DHARMa::simulateResiduals(mod_struct1,plot = T) #OK

vif(mod_struct1) #<=2 -> OK

model_dredge1 = MuMIn::dredge(global.model = mod_struct1, 
                              options(na.action = "na.fail"),trace = 2 )
MuMIn::sw(MuMIn::get.models(model_dredge1,subset = delta<4))

d_Struct1=Perform_model_averaging(model_dredge1,4)


Plot_correlation_variables(d[,c("Isothermality_Holocene","Isothermality_Current","MAT_Current",
                                "MAT_LGM","MAT_Holocene")])

mod_struct2=lmer(Struct2~Isothermality_Holocene + MAT_LGM +MAT_Current+
                   Isothermality_Current+
                   Desert_LGM + Desert_Current+
                   (1|Site_ID),
                 data = d,na.action = na.fail)
DHARMa::simulateResiduals(mod_struct2,plot = T) #not OK

vif(mod_struct2) #we remove MaxMATcurrent because of VIF ISSUES

model_dredge2 = MuMIn::dredge(global.model = mod_struct2, 
                              options(na.action = "na.fail"),trace = 2 )
MuMIn::sw(MuMIn::get.models(model_dredge2,subset = delta<4))

d_Struct2=Perform_model_averaging(model_dredge2,4)

#Extracting effects
par(mfrow=c(3,3))
partial_res=visreg::visreg(MuMIn::get.models(model_dredge1,subset = delta < 4)[[1]],plot=T)

assembling_point=lapply(partial_res,function(x){
  partial_res_x=x$res%>%
    dplyr::select_if(., function(x) length(unique(x))>1)
  return(tibble(Value=as.numeric(partial_res_x[,1]),Response=partial_res_x[,2],Name=colnames(partial_res_x)[1]))
})%>%bind_rows(.)

assembling_lines=lapply(partial_res,function(x){
  partial_res_x=x$fit%>%
    dplyr::select_if(., function(x) length(unique(x))>1)
  return(tibble(Value=as.numeric(partial_res_x[,1]),Q2=partial_res_x[,2],
                Q1=partial_res_x[,3],Q3=partial_res_x[,4],Name=colnames(partial_res_x)[1]))
})%>%bind_rows(.)

p1_1=ggplot(assembling_point%>%
              dplyr::filter(.,Name!='Desert_LGM')%>%
              Change_name_RF(., name_predictor = "Name"))+
  geom_point(aes(x=Value,y=Response),alpha=.5,color="grey70")+
  geom_line(data=assembling_lines%>%
              dplyr::filter(.,Name!='Desert_LGM')%>%
              Change_name_RF(., name_predictor = "Name"),
            aes(x=Value,y=Q2,group=Name_stat),color="pink")+
  geom_ribbon(data=assembling_lines%>%
                dplyr::filter(.,Name!='Desert_LGM')%>%
                Change_name_RF(., name_predictor = "Name"),
              aes(x=Value,y=Q2,ymin=Q1,ymax=Q3,group=Name_stat),fill="pink",alpha=.5)+
  facet_wrap(.~Name_stat,scales="free",switch = "x",nrow = 2)+
  labs(x="",y="PC1, spatial structure")+
  the_theme2+  theme(strip.placement = "outside")

p1_2=ggplot(assembling_point%>%
              dplyr::filter(.,Name=='Desert_LGM')%>%
              mutate(., Value=as.character(Value))%>%
              mutate(., Value=recode_factor(Value,
                                            "No desert during LGM"=unique(Value)[1],
                                            "Desert during LGM"=unique(.$Value)[2]))%>%
              Change_name_RF(., name_predictor = "Name"))+
  geom_boxplot(aes(x=(Value),y=Response),alpha=.5,color="pink")+
  facet_wrap(.~Name_stat,scales="free",switch = "x",nrow = 2)+
  labs(x="",y="PC1, spatial structure")+
  scale_x_discrete(labels=c("No desert during LGM","Desert during LGM"))+
  the_theme2+  theme(strip.placement = "outside")+
  theme(axis.text.x =  element_text(angle=60,hjust=1))



#Extracting effects
partial_res=visreg::visreg(MuMIn::get.models(model_dredge2,subset = delta < 4)[[1]],plot=T)

assembling_point=lapply(partial_res,function(x){
  partial_res_x=x$res%>%
    dplyr::select_if(., function(x) length(unique(x))>1)
  return(tibble(Value=as.numeric(partial_res_x[,1]),Response=partial_res_x[,2],Name=colnames(partial_res_x)[1]))
})%>%bind_rows(.)

assembling_lines=lapply(partial_res,function(x){
  partial_res_x=x$fit%>%
    dplyr::select_if(., function(x) length(unique(x))>1)
  return(tibble(Value=as.numeric(partial_res_x[,1]),Q2=partial_res_x[,2],
                Q1=partial_res_x[,3],Q3=partial_res_x[,4],Name=colnames(partial_res_x)[1]))
})%>%bind_rows(.)

p2=ggplot(assembling_point%>%
            Change_name_RF(., name_predictor = "Name"))+
  geom_point(aes(x=Value,y=Response),alpha=.5,color="grey70")+
  geom_line(data=assembling_lines%>%Change_name_RF(., name_predictor = "Name"),
            aes(x=Value,y=Q2,group=Name_stat),color="pink")+
  geom_ribbon(data=assembling_lines%>%Change_name_RF(., name_predictor = "Name"),
              aes(x=Value,y=Q2,ymin=Q1,ymax=Q3,group=Name_stat),fill="pink",alpha=.5)+
  facet_wrap(.~Name_stat,scales="free",switch = "x",nrow = 1)+
  labs(x="",y="PC2, spatial structure")+
  the_theme2+  theme(strip.placement = "outside")


ggsave("../Figures/Final_figs/Mixed_models.pdf",
       #ggarrange(
       # ggarrange(d_Struct1$Fig+ggtitle("PC1, spatial structure"),
       #           d_Struct2$Fig+ggtitle("PC2, spatial structure"),
       #           nrow=2,labels = c("a","c"),align="hv"),
       ggarrange(ggarrange(p1_1,ggarrange(ggplot()+theme_void(),
                                          p1_2,ggplot()+theme_void(),
                                          nrow=3,heights = c(.2,1,.1)),ncol=2,widths = c(1,.8)),
                 ggarrange(ggplot()+theme_void(),
                           p2,ggplot()+theme_void(),
                           nrow=3,heights = c(.2,1,.1)),
                 nrow=2,labels = c("b","d"))#,ncol=2,align="hv"
       #)
       ,width = 8,height = 8)



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
  RDA_mod_LGM= rda(d$Struct1 ~ .,d_LGM)
  Select_var_RDA_LGM = Select_variable_VarPart(RDA_all = RDA_mod_LGM,
                                               data = d_LGM,
                                               response_var = d$Struct1,
                                               forwardsel = F)
  #Holocene clim
  d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
  RDA_mod_Holo= rda(d$Struct1 ~ ., data = d_Holo)
  Select_var_RDA_Holo = Select_variable_VarPart(RDA_all = RDA_mod_Holo,
                                                data = d_Holo,
                                                response_var = d$Struct1,
                                                forwardsel = F)
  
  #Current clim
  d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept,20)]])]
  RDA_mod_Current= rda(d$Struct1 ~ ., data = d_Current)
  Select_var_RDA_Current = Select_variable_VarPart(RDA_all = RDA_mod_Current,
                                                   data = d_Current,
                                                   response_var = d$Struct1,
                                                   forwardsel = F)
  
  #Getting subset from RDA
  Var_LGM = subset(d_LGM, select = Select_var_RDA_LGM)
  Var_Holo = subset(d_Holo, select = Select_var_RDA_Holo)
  Var_Current = subset(d_Current, select = Select_var_RDA_Current)
  
  # assign(paste0("d_",ID_PC),
  #        Plot_variance_partition(varpart(d[,c("Struct1")],Var_LGM,Var_Holo,Var_Current),
  #                                bar = T,four_groups = F))
  assign(paste0("d_",ID_PC),
         Plot_variance_partition(varpart(d[,c("Struct1")],d_LGM,d_Holo,d_Current),
                                 bar = T,four_groups = F))
  
  # combined_d=cbind(cbind(Var_LGM),
  #                  Var_Holo,
  #                  cbind(Var_Current),
  #                  d$Struct1)
  combined_d=cbind(cbind(d_LGM),
                   d_Holo,
                   cbind(d_Current),
                   d$Struct1)
  colnames(combined_d)[ncol(combined_d)]="Struct1"
  train_id=sample(1:nrow(combined_d),round(.7*nrow(combined_d)))
  
  
  RDA_mod_RF= rda(combined_d$Struct1 ~ ., data = combined_d)
  Select_var_RDA_Env = Select_variable_VarPart(RDA_all = RDA_mod_RF,
                                               data = combined_d,
                                               response_var = d$Struct1,
                                               forwardsel = F)
  
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
                  data = combined_d[train_id,],#c(Select_var_RDA_Env,"Struct1")],
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
                  # dplyr::group_by(Name_stat)%>%
                  # dplyr::summarise(., mean_imp=mean(Struct1),
                  #                  sd_imp=sd(Struct1))%>%
                  # dplyr::ungroup(.)%>%
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
           scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Current"="pink",
                                      "Shared past clim."="#DAB2BA","Other shared variance"="#E8E8E8"))+
           labs(y="Mean Importance",x="",fill=""))
}  

pvariance=ggarrange(d_1$Figure+ggtitle("PC1, spatial structure")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "none"),
                    d_2$Figure+ggtitle("PC2, spatial structure")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "none"),
                    ncol=2,common.legend = T,legend="none")

PCA_stat=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.,return_PCA = T)

p0=plotPCA(PCA_stat,1,2,F,xLim = c(-3,3),yLim = c(-3,3),colorLow = "lightgrey",colorMid = "grey40",colorHigh = "grey10")

p_RF=ggarrange(
  pRF_1+ggtitle("PC1, spatial structure")+theme(legend.position = "none"),
  pRF_2+ggtitle("PC2, spatial structure")+theme(legend.position = "none"),
  ncol=2,common.legend = T,legend = "none")


p_legend=get_legend(ggplot(tibble(ID=1:5,Name=c("LGM","Holocene",
                                                "Current","Shared past clim.",
                                                "Other shared variance")))+
                      geom_bar(aes(x=1,y=1,fill=as.factor(Name)),stat="identity")+
                      scale_fill_manual(values=c("LGM"="#FFF9BF",
                                                 "Holocene"="#FFD2A0",
                                                 "Current"="pink",
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

for (ID_PC in 1:2){
  
  d=read.table("./Data/data_sites_Regular_Buxton.csv",sep=";")%>%
    Closer_to_normality(.)%>%
    Perform_PCA_spatial_struc(.,F)%>%
    dplyr::select_if(function(x) length(unique(x)) > 1 & !any(is.na(x)))%>%
    # mutate(., Cos_long=cos(Longitude),Sin_long=sin(Longitude))%>%
    mutate(across(colnames(dplyr::select_if(., is.numeric)), ~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
  
  d$Struct1=d[,paste0("Struct",ID_PC)]
  
  #LGM clim
  d_LGM=d[,save_names_LGM[-6]]
  RDA_mod_LGM= rda(d$Struct1 ~ .,d_LGM)
  Select_var_RDA_LGM = Select_variable_VarPart(RDA_all = RDA_mod_LGM,
                                               data = d_LGM,
                                               response_var = d$Struct1,
                                               forwardsel = T)
  #Holocene clim
  d_Holo=d[,save_names_Holo[-6]]
  RDA_mod_Holo= rda(d$Struct1 ~ ., data = d_Holo)
  Select_var_RDA_Holo = Select_variable_VarPart(RDA_all = RDA_mod_Holo,
                                                data = d_Holo,
                                                response_var = d$Struct1,
                                                forwardsel = T)
  
  #Current clim
  d_Current=d[,save_names_Current[-6]]
  RDA_mod_Current= rda(d$Struct1 ~ ., data = d_Current)
  Select_var_RDA_Current = Select_variable_VarPart(RDA_all = RDA_mod_Current,
                                                   data = d_Current,
                                                   response_var = d$Struct1,
                                                   forwardsel = T)
  
  #Getting subset from RDA
  Var_LGM = subset(d_LGM, select = Select_var_RDA_LGM)
  Var_Holo = subset(d_Holo, select = Select_var_RDA_Holo)
  Var_Current = subset(d_Current, select = Select_var_RDA_Current)
  
  
  assign(paste0("d_",ID_PC),
         Plot_variance_partition(varpart(d[,c("Struct1")],d_LGM,d_Holo,d_Current),
                                 bar = T,four_groups = F))
  combined_d=cbind(cbind(d_LGM),
                   d_Holo,
                   cbind(d_Current),
                   d$Struct1)
  colnames(combined_d)[ncol(combined_d)]="Struct1"
  train_id=sample(1:nrow(combined_d),round(.7*nrow(combined_d)))
  
  
  RDA_mod_RF= rda(combined_d$Struct1 ~ ., data = combined_d)
  Select_var_RDA_Env = Select_variable_VarPart(RDA_all = RDA_mod_RF,
                                               data = combined_d,
                                               response_var = d$Struct1,
                                               forwardsel = F)
  
  #Changing to factors
  if (any(colnames(combined_d)=="Desert_Current")){
    combined_d[,"Desert_Current"]=as.factor(combined_d[,"Desert_Current"])
  }
  
  if (any(colnames(combined_d)=="Desert_LGM")){
    combined_d[,"Desert_LGM"]=as.factor(combined_d[,"Desert_LGM"])
  }
  
  BRT_model = gbm(Struct1 ~ .,
                  distribution = "gaussian",
                  data = combined_d[train_id,],#c(Select_var_RDA_Env,"Struct1")],
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
           scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Current"="pink",
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


p_legend=get_legend(ggplot(tibble(ID=1:5,Name=c("LGM","Holocene",
                                                "Current","Shared past clim.",
                                                "Other shared variance")))+
                      geom_bar(aes(x=1,y=1,fill=as.factor(Name)),stat="identity")+
                      scale_fill_manual(values=c("LGM"="#FFF9BF",
                                                 "Holocene"="#FFD2A0",
                                                 "Current"="pink",
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
