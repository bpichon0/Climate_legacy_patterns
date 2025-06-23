rm(list=ls())
source("./Legacy_functions.R")

# --------------------- Step 1: Irregular biocom patterns ----
## >> 1) Recompute spatial statistics ----

library(parallel)
list_landscape=list.files("../Data/Images_irregular/landscapes/")
runing_spatialstats=mclapply(X = 1:length(list_landscape),
                             FUN=function(x){
                               info_biocom=read.table("../Data/data_sites.csv",sep=";")
                               landscape=as.matrix(read.table(paste0("../Data/Images_irregular/landscapes/",list_landscape[x]),sep="\t"))
                               name_x=gsub(".txt","",list_landscape[x])
                               d_sumstat=Get_sumstat(landscape,slope = info_biocom$Slope[which(info_biocom$File_ID==name_x)],compute_KS = F)%>%
                                 add_column(., Name_plot=name_x)
                               return(d_sumstat)
                             },mc.cores = 35)%>%
  bind_rows(.)

write.table(runing_spatialstats,"../Data/Spatial_stats.csv",sep=";")

## >> 2) Get past & current climatic data ----

database=read.table("../Data/data_sites.csv",sep=";")
dryland_sf = st_as_sf(database[,c("Longitude","Lattitude")], coords = c("Longitude", "Lattitude"), crs = 4326)

name_bioclim=paste0(c("MAT","MAT_range_diurn","Isothermality","MAT_Seasonality","MaxMAT","MinMAT",
                      "MAT_range","MATwet","MATdry","MATwarm","MATcol",
                      "MAP","MAPwet","MAPdry","MAP_Seasonality","MAPwetQ","MAPdryQ","MAPwarmQ","MAPcoldQ"),"_")
list_bioclim=c(1:19)

data_clim=data.frame(File_ID=database$File_ID)
for (bioclim_id in list_bioclim){ #ID of the bioclimatic climate variable
  
  for (period_acronyme in c("mid","lgm")){ #Two climatic periods
    
    data_clim_id=tibble()
    
    for (period in list.files("../Data/Paleo_clim/",period_acronyme)){ # for each climatic model of the period
      name_period=ifelse(any(grep("mid",period)),"Holocene","LGM")
      nc_data_aridity=rast(paste0("../Data/Paleo_clim/",period,"/",gsub("_2-5m","",period),bioclim_id,".tif"))
      
      climate_values = raster::extract(nc_data_aridity, dryland_sf,method="bilinear")[,-1]
      data_clim_id=rbind(data_clim_id,climate_values)
    }
    
    data_clim_id=as.data.frame(t(data_clim_id))
    data_clim_id$File_ID=database$File_ID
    data_clim_id$ID=database$ID
    data_clim_id$Site_ID=database$Site_ID

    data_clim_id=data_clim_id%>%
      melt(.,id.vars = c("File_ID","ID","Site_ID"))%>%
      dplyr::group_by(., File_ID,ID,Site_ID)%>%
      dplyr::summarise(.,.groups = "keep",Mean_value=mean(value))%>%
      ungroup(.)
    
    colnames(data_clim_id)[4]=paste0(name_bioclim[bioclim_id],name_period)
    data_clim=cbind(data_clim,data_clim_id[,4])
  }
}

#Current data
data_clim_current=tibble()
for (bioclim_id in list_bioclim){
  name_period="Current"
  nc_data_aridity=rast(paste0("../Data/Paleo_clim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_",bioclim_id,".tif"))
  climate_values = raster::extract(nc_data_aridity, dryland_sf,method="bilinear")[,-1]
  data_clim_current=rbind(data_clim_current,climate_values)
}

data_clim_current=as.data.frame(t(data_clim_current))
colnames(data_clim_current)=paste0(c("MAT","MAT_range_diurn","Isothermality","MAT_Seasonality","MaxMAT","MinMAT",
                                     "MAT_range","MATwet","MATdry","MATwarm","MATcol",
                                     "MAP","MAPwet","MAPdry","MAP_Seasonality","MAPwetQ","MAPdryQ","MAPwarmQ","MAPcoldQ"),"_Current")

data_clim$Site_ID=data_clim_id$Site_ID
data_clim=cbind(data_clim,data_clim_current)

database=cbind(database,data_clim[,-1])

# data_LGM=readxl::read_xls("../Data/Data_past_biome_biocom.xls")%>%
#   dplyr::rename(., Site_ID=`Site ID`)

biomes_irregular = read.table("../Data/Biomes_irregular.csv",sep=";",header = T)
database$Biome_LGM=biomes_irregular$Past_biome
database$Biome_Current=biomes_irregular$Current_biome

# database$pH=sapply(1:nrow(database),function(x){return(data_LGM$pH[which(data_LGM$Site_ID==database$Site_ID[x])])})


d_stats=read.table("../Data/Spatial_stats.csv",sep=";")
d_stats$Site_ID=sapply(1:nrow(d_stats),function(x){return(database$Site_ID[which(database$File_ID==d_stats$Name_plot[x])])})
database=database[,-c(14:24)]
save_database=database

#Computing legacies
for (period in c("LGM","Holocene")){
  for (metric_name in colnames(database)[grep(period,colnames(database))[1:19]]){
    if (any(grep("MAT",metric_name))){
      database[,metric_name]=-(database[,metric_name]/10-database[,gsub(paste0("_",period),"_Current",metric_name)])
    }else{
      database[,metric_name]=-(database[,metric_name]-database[,gsub(paste0("_",period),"_Current",metric_name)])
    }
  }
}

#Taken from https://figshare.com/articles/dataset/Dataset_and_R_code_from_the_article_The_interplay_between_facilitation_and_habitat_type_drive_spatial_vegetation_patterns_in_global_drylands_/5640193?file=11863505
d_berdugo=read.table("../Data/data_berdugo_2019.csv",sep=",",header = T)

database$File_ID=as.numeric(gsub("-b","",gsub("-a","",gsub("-c","",database$File_ID))))

database$CWM_Height=sapply(1:nrow(database),function(x){return(d_berdugo$HEIGHT[which(d_berdugo$plotn==database$File_ID[x])])})
database$Woody=sapply(1:nrow(database),function(x){return(d_berdugo$PercRWC[which(d_berdugo$plotn==database$File_ID[x])])})
database$Facilitation=sapply(1:nrow(database),function(x){return(d_berdugo$Facil[which(d_berdugo$plotn==database$File_ID[x])])})
database$Type_veg=sapply(1:nrow(database),function(x){return(d_berdugo$VEG[which(d_berdugo$plotn==database$File_ID[x])])})

##Adding the stats
write.table(cbind(database,d_stats[,-c(ncol(d_stats)-1,ncol(d_stats))]),"./Data/data_sites_CLIM.csv",sep=";")


#Taken from https://figshare.com/articles/dataset/Dataset_and_R_code_from_the_article_The_interplay_between_facilitation_and_habitat_type_drive_spatial_vegetation_patterns_in_global_drylands_/5640193?file=11863505
d_berdugo=read.table("../Data/data_berdugo_2019.csv",sep=",",header = T)
database=save_database

database$File_ID=as.numeric(gsub("-b","",gsub("-a","",gsub("-c","",database$File_ID))))

database$CWM_Height=sapply(1:nrow(database),function(x){return(d_berdugo$HEIGHT[which(d_berdugo$plotn==database$File_ID[x])])})
database$Woody=sapply(1:nrow(database),function(x){return(d_berdugo$PercRWC[which(d_berdugo$plotn==database$File_ID[x])])})
database$Facilitation=sapply(1:nrow(database),function(x){return(d_berdugo$Facil[which(d_berdugo$plotn==database$File_ID[x])])})
database$Type_veg=sapply(1:nrow(database),function(x){return(d_berdugo$VEG[which(d_berdugo$plotn==database$File_ID[x])])})

##Adding the stats
write.table(cbind(database,d_stats[,-c(ncol(d_stats)-1,ncol(d_stats))]),"./Data/data_sites_CLIM_nolegacy.csv",sep=";")




# --------------------- Step 2: Irregular biodesert patterns ----
## >> 1) Get past & current climatic data ----

database=read.table("../Data/data_biodesert.csv",sep=";")
dryland_sf = st_as_sf(database[,c("Longitude","Lattitude")], coords = c("Longitude", "Lattitude"), crs = 4326)

name_bioclim=paste0(c("MAT","MAT_range_diurn","Isothermality","MAT_Seasonality","MaxMAT","MinMAT",
                      "MAT_range","MATwet","MATdry","MATwarm","MATcol",
                      "MAP","MAPwet","MAPdry","MAP_Seasonality","MAPwetQ","MAPdryQ","MAPwarmQ","MAPcoldQ"),"_")
list_bioclim=c(1:19)

data_clim=data.frame(File_ID=database$ID)
for (bioclim_id in list_bioclim){ #ID of the bioclimatic climate variable
  
  for (period_acronyme in c("mid","lgm")){ #Two climatic periods
    
    data_clim_id=tibble()
    
    for (period in list.files("../Data/Paleo_clim/",period_acronyme)){ # for each climatic model of the period
      name_period=ifelse(any(grep("mid",period)),"Holocene","LGM")
      nc_data_aridity=rast(paste0("../Data/Paleo_clim/",period,"/",gsub("_2-5m","",period),bioclim_id,".tif"))
      
      climate_values = raster::extract(nc_data_aridity, dryland_sf,method="bilinear")[,-1]
      data_clim_id=rbind(data_clim_id,climate_values)
    }
    
    data_clim_id=as.data.frame(t(data_clim_id))
    data_clim_id$File_ID=1:504
    data_clim_id$ID=database$ID

    
    data_clim_id=data_clim_id%>%
      melt(.,id.vars = c("File_ID","ID"))%>%
      dplyr::group_by(., File_ID,ID)%>%
      dplyr::summarise(.,.groups = "keep",Mean_value=mean(value))%>%
      ungroup(.)
    
    colnames(data_clim_id)[3]=paste0(name_bioclim[bioclim_id],name_period)
    data_clim=cbind(data_clim,data_clim_id[,3])
  }
}

#Current data
data_clim_current=tibble()
for (bioclim_id in list_bioclim){
  name_period="Current"
  nc_data_aridity=rast(paste0("../Data/Paleo_clim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_",bioclim_id,".tif"))
  climate_values = raster::extract(nc_data_aridity, dryland_sf,method="bilinear")[,-1]
  data_clim_current=rbind(data_clim_current,climate_values)
}

data_clim_current=as.data.frame(t(data_clim_current))
colnames(data_clim_current)=paste0(c("MAT","MAT_range_diurn","Isothermality","MAT_Seasonality","MaxMAT","MinMAT",
                                     "MAT_range","MATwet","MATdry","MATwarm","MATcol",
                                     "MAP","MAPwet","MAPdry","MAP_Seasonality","MAPwetQ","MAPdryQ","MAPwarmQ","MAPcoldQ"),"_Current")
data_clim=cbind(data_clim,data_clim_current)
database=cbind(database,data_clim[,-1])

#Computing legacies
for (period in c("LGM","Holocene")){
  for (metric_name in colnames(database)[grep(period,colnames(database))[1:19]]){
    if (any(grep("MAT",metric_name))){
      database[,metric_name]=database[,metric_name]/10-database[,gsub(paste0("_",period),"_Current",metric_name)]
    }else{
      database[,metric_name]=database[,metric_name]-database[,gsub(paste0("_",period),"_Current",metric_name)]
    }
  }
}


##Adding the stats
write.table(database,"../Data/data_biodesert_CLIM.csv",sep=";")


# --------------------- Step 3: Regular patterns ----
## >> 1) Aggregating climate and spatial statistics ----

d=read.table("../Data/Images_regular/Regular_patterns.csv",sep=";",header = T)
coords = d[,c("POINT_X","POINT_Y")]
coords=as.data.frame(apply(coords,2,function(x){
  return(as.numeric(gsub(",",".",x)))
}))
colnames(coords)=c("POINT_X","POINT_Y")
coords_sf = st_as_sf(coords, coords = c("POINT_X","POINT_Y"), crs = 32633)
dryland_sf = st_transform(coords_sf, crs = 4326)

name_bioclim=paste0(c("MAT","MAT_range_diurn","Isothermality","MAT_Seasonality","MaxMAT","MinMAT",
                      "MAT_range","MATwet","MATdry","MATwarm","MATcol",
                      "MAP","MAPwet","MAPdry","MAP_Seasonality","MAPwetQ","MAPdryQ","MAPwarmQ","MAPcoldQ"),"_")
list_bioclim=c(1:19)

data_clim=data.frame(File_ID=d$OBJECTID)
for (bioclim_id in list_bioclim){ #ID of the bioclimatic climate variable
  
  for (period_acronyme in c("mid","lgm")){ #Two climatic periods
    
    data_clim_id=tibble()
    
    for (period in list.files("../Data/Paleo_clim/",period_acronyme)){ # for each climatic model of the period
      name_period=ifelse(any(grep("mid",period)),"Holocene","LGM")
      nc_data_aridity=rast(paste0("../Data/Paleo_clim/",period,"/",gsub("_2-5m","",period),bioclim_id,".tif"))
      
      climate_values = raster::extract(nc_data_aridity, dryland_sf,method="bilinear")[,-1]
      data_clim_id=rbind(data_clim_id,climate_values)
    }
    
    data_clim_id=as.data.frame(t(data_clim_id))
    data_clim_id$File_ID=d$FID
    
    
    data_clim_id=data_clim_id%>%
      melt(.,id.vars = c("File_ID"))%>%
      dplyr::group_by(., File_ID)%>%
      dplyr::summarise(.,.groups = "keep",Mean_value=mean(value))%>%
      ungroup(.)
    
    colnames(data_clim_id)[2]=paste0(name_bioclim[bioclim_id],name_period)
    data_clim=cbind(data_clim,data_clim_id[,2])
  }
}

#Current data
data_clim_current=tibble()
for (bioclim_id in list_bioclim){
  name_period="Current"
  nc_data_aridity=rast(paste0("../Data/Paleo_clim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_",bioclim_id,".tif"))
  climate_values = raster::extract(nc_data_aridity, dryland_sf,method="bilinear")[,-1]
  data_clim_current=rbind(data_clim_current,climate_values)
}

data_clim_current=as.data.frame(t(data_clim_current))
colnames(data_clim_current)=paste0(c("MAT","MAT_range_diurn","Isothermality","MAT_Seasonality","MaxMAT","MinMAT",
                                     "MAT_range","MATwet","MATdry","MATwarm","MATcol",
                                     "MAP","MAPwet","MAPdry","MAP_Seasonality","MAPwetQ","MAPdryQ","MAPwarmQ","MAPcoldQ"),"_Current")

data_clim$Site_ID=data_clim_id$File_ID
data_clim=cbind(data_clim,data_clim_current)

database=cbind(data_clim,data.frame(Skewness=as.numeric(gsub(",",".",d$Skewness))))

database=database%>%
  add_column(., Type_pattern=sapply(1:nrow(database),function(x){
    if (.$Skewness[x]< -.5){
      return("Spots")
    }else if (.$Skewness[x] < 0){
      return("Labyrinths")
    }else {
      return("Gaps")
    }
  }))

#Computing legacies
for (period in c("LGM","Holocene")){
  for (metric_name in colnames(database)[grep(period,colnames(database))[1:19]]){
    if (any(grep("MAT",metric_name))){
      database[,metric_name]=database[,metric_name]/10-database[,gsub(paste0("_",period),"_Current",metric_name)]
    }else{
      database[,metric_name]=database[,metric_name]-database[,gsub(paste0("_",period),"_Current",metric_name)]
    }
  }
}

##Adding the stats
write.table(database,"../Data/data_regular.csv",sep=";")



# --------------------- Step 4: Regular patterns Buxton ----

## >> 1) Extract images ----


list_buxton=list.files("../Data/Regular_images_Buxton/Raw/","RGB")
dir.create("../Data/Regular_images_Buxton/Unzipped",showWarnings = F)
for (k in list_buxton){ #extracting .gz
  untar(paste0("../Data/Regular_images_Buxton/Raw/",k),exdir = paste0("../Data/Regular_images_Buxton/Unzipped/",gsub(".tar.gz","",k)))
}


list_sites_regular=list.files("../Data/Regular_images_Buxton/Unzipped/")

pdf("../Figures/All_images_Buxton.pdf",width = 6,height = 6)
for (k in list_sites_regular){
  
  list_images_years=list.files(paste0("../Data/Regular_images_Buxton/Unzipped/",k))
  img=imager::load.image(paste0("../Data/Regular_images_Buxton/Unzipped/",k,"/",list_images_years[length(list_images_years)])) #last image
  plot(img)
  
  #and binarize the kmean output 
  kmean_img=k_means_RGB(img,2)
  cats=get_cut_grayscale_values(2)[[1]]
  mat=kmean_img %>% binarize(cats[[1]], cats[[2]])
  
  write.table(mat,paste0("../Data/Regular_images_Buxton/Binary/Binary_",k,".csv"),sep=";",row.names = F,col.names = F)
  
}

dev.off()


#Binarize images
for (k in list_sites_regular){
  list_images_years=list.files(paste0("../Data/Regular_images_Buxton/Unzipped/",k))
  img=imager::load.image(paste0("../Data/Regular_images_Buxton/Unzipped/",k,"/",list_images_years[length(list_images_years)])) #last image
  # and binarize the kmean output
  kmean_img=k_means_RGB(img,2)
  cats=get_cut_grayscale_values(2)[[1]]
  mat=kmean_img %>% binarize(cats[[1]], cats[[2]])
}

d_all=tibble()
for (k in list.files("../Data/Regular_images_Buxton/Binary/")){
  list_images_years=list.files(paste0("../Data/Regular_images_Buxton/Binary/",k))
  d_sumstat=Get_sumstat(read.table(paste0("../Data/Regular_images_Buxton/Binary/",k),sep=";")%>%as.matrix(.),slope = 0,compute_KS = F)%>%
    add_column(., Name_plot=strsplit(k,"-")[[1]][1])
  d_all=rbind(d_sumstat,d_all)
  print(strsplit(k,"-")[[1]][1])
}
write.table(d_all,'./Spatial_stats_Buxton.csv',sep=";")


## >> 2) Get past & current climatic data ----

#From Buxton et al., GCB 2021
data = data.frame(
  ID = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 18, 20, 21, 23, 25, 26, 27, 28, 29, 30, 31, 48, 49,
         50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66),
  Latitude = c(11.58, 11.12, 10.96, 13.12, 13.17, 15.2, 15.09, 15.8, 15.11, 15.03, 15.34, 14.85,
               14.97, 15.02, 16.19, 16.17, 16.48, 15.95, 15.86, 14.8, 14.94, 15.48, 15.57, 15.58,
               12.58, 12.7, 12.54, 13.12, 11.07, 11.28, 11.27, 11.47, 11.51, 11.22, 11.62, 11.32,
               11.37, 11.6, 11.46, 11.71),
  Longitude = c(27.94, 28.37, 28.2, 2.59, 1.58, -15.2, -15.04, -14.36, -14.53, -0.87, -1.15, -1.43,
                -1.12, -1.35, -1.83, -2.03, -1.87, -1.52, -2.05, -3.38, -3.56, -5.83, -5.92, -13, 
                3.75, 2.63, 2.26, 2.17, 27.93, 27.96, 27.55, 27.97, 27.87, 27.73, 27.86, 27.88, 
                27.68, 27.73, 27.68, 27.91),
  Type = c("Spots", "Labyrinths", "Gaps", "Labyrinths", "Labyrinths", "Labyrinths", "Labyrinths",
           "Gaps", "Gaps", "Spot–labyrinths", "Spot–labyrinths", "Spot–labyrinths", "Spot–labyrinths",
           "Spot–labyrinths", "Spot–labyrinths", "Spot–labyrinths", "Spot–labyrinths", "Spot–labyrinths",
           "Spot–labyrinths", "Labyrinths", "Labyrinths", "Labyrinths", "Labyrinths","Gaps", "Labyrinths",
           "Labyrinths", "Gaps","Labyrinths","Gaps",
           "Gaps", "Spots", "Spots", "Spots", "Spots", "Spots", "Spots", "Spots", "Spots", "Spots","Spots")
)

database=read.table("../Data/Spatial_stats_Buxton.csv",sep=";")%>%
  add_column(., ID=gsub("Binary_","",.$Name_plot))%>%
  dplyr::filter(., as.numeric(ID) %in% data$ID)%>%
  dplyr::arrange(.,ID)%>%
  add_column(., 
             Latitude=data$Latitude,
             Longitude=data$Longitude,
             Type=data$Type)


dryland_sf = st_as_sf(database[,c("Longitude","Latitude")], coords = c("Longitude", "Latitude"), crs = 4326)

name_bioclim=paste0(c("MAT","MAT_range_diurn","Isothermality","MAT_Seasonality","MaxMAT","MinMAT",
                      "MAT_range","MATwet","MATdry","MATwarm","MATcol",
                      "MAP","MAPwet","MAPdry","MAP_Seasonality","MAPwetQ","MAPdryQ","MAPwarmQ","MAPcoldQ"),"_")
list_bioclim=c(1:19)

data_clim=data.frame(ID=database$ID)
for (bioclim_id in list_bioclim){ #ID of the bioclimatic climate variable
  
  for (period_acronyme in c("mid","lgm")){ #Two climatic periods
    
    data_clim_id=tibble()
    
    for (period in list.files("../Data/Paleo_clim/",period_acronyme)){ # for each climatic model of the period
      name_period=ifelse(any(grep("mid",period)),"Holocene","LGM")
      nc_data_aridity=rast(paste0("../Data/Paleo_clim/",period,"/",gsub("_2-5m","",period),bioclim_id,".tif"))
      
      climate_values = raster::extract(nc_data_aridity, dryland_sf,method="bilinear")[,-1]
      data_clim_id=rbind(data_clim_id,climate_values)
    }
    
    data_clim_id=as.data.frame(t(data_clim_id))
    data_clim_id$ID=database$ID
    
    
    data_clim_id=data_clim_id%>%
      melt(.,id.vars = c("ID"))%>%
      dplyr::group_by(., ID)%>%
      dplyr::summarise(.,.groups = "keep",Mean_value=mean(value))%>%
      ungroup(.)
    
    colnames(data_clim_id)[2]=paste0(name_bioclim[bioclim_id],name_period)
    data_clim=cbind(data_clim,data_clim_id[,2])
  }
}

#Current data
data_clim_current=tibble()
for (bioclim_id in list_bioclim){
  name_period="Current"
  nc_data_aridity=rast(paste0("../Data/Paleo_clim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_",bioclim_id,".tif"))
  climate_values = raster::extract(nc_data_aridity, dryland_sf,method="bilinear")[,-1]
  data_clim_current=rbind(data_clim_current,climate_values)
}

data_clim_current=as.data.frame(t(data_clim_current))
colnames(data_clim_current)=paste0(c("MAT","MAT_range_diurn","Isothermality","MAT_Seasonality","MaxMAT","MinMAT",
                                     "MAT_range","MATwet","MATdry","MATwarm","MATcol",
                                     "MAP","MAPwet","MAPdry","MAP_Seasonality","MAPwetQ","MAPdryQ","MAPwarmQ","MAPcoldQ"),"_Current")

data_clim$ID=data_clim_id$ID
data_clim=cbind(data_clim,data_clim_current)
database=cbind(database,data_clim[,-1])

#Computing legacies
for (period in c("LGM","Holocene")){
  for (metric_name in colnames(database)[grep(period,colnames(database))[1:19]]){
    if (any(grep("MAT",metric_name))){
      database[,metric_name]=database[,metric_name]/10-database[,gsub(paste0("_",period),"_Current",metric_name)]
    }else{
      database[,metric_name]=database[,metric_name]-database[,gsub(paste0("_",period),"_Current",metric_name)]
    }
  }
}

biomes_regular = read.table("../Data/Biomes_regular.csv",sep=";",header = T)
database$Biome_LGM=biomes_regular$Past_biome
database$Biome_Current=biomes_regular$Current_biome

##Adding the stats
write.table(database,"../Data/data_sites_Regular_Buxton.csv",sep=";")

