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
write.table(cbind(database,d_stats[,-c(ncol(d_stats)-1,ncol(d_stats))]),"./Data/data_sites_CLIM.csv",sep=";")




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


# --------------------- Step 5: Chelsea database ----
## >> 1) Data collection ----

library(purrr)
library(httr) 
library(parallel)
library(terra)


download_chelsa_file=function(time_step) {
  time_str=sprintf("%04d", time_step)  # e.g., "-0210"
  filename=paste0("CHELSA_TraCE21k_bio01_", time_str, "_V.1.0.tif")
  
  base_url="https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/bioclim/bio01/"
  url=paste0(base_url, filename)
  
  # Local file path
  local_path=file.path("./all_data", filename)
  
  # Skip if already downloaded
  
  # Try to download
  tryCatch({
    message("Downloading: ", filename)
    
    # Use httr for more robust downloading
    response=GET(url, write_disk(local_path, overwrite = TRUE), 
                    timeout(300), progress())
    
    if(response$status_code == 200) {
      message("Success: ", filename)
      return(local_path)
    } else {
      message("Failed (HTTP ", response$status_code, "): ", filename)
      file.remove(local_path)  # Remove partial file
      return(NULL)
    }
  }, error = function(e) {
    message("Error downloading ", filename, ": ", e$message)
    if(file.exists(local_path)) file.remove(local_path)
    return(NULL)
  })
}
mclapply(seq(-200, 20, by = 1),download_chelsa_file,mc.cores = 40)


download_chelsa_file=function(time_step) {
  time_str=sprintf("%04d", time_step)  # e.g., "-0210"
  filename=paste0("CHELSA_TraCE21k_bio03_", time_str, "_V.1.0.tif")
  
  base_url="https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/bioclim/bio03/"
  url=paste0(base_url, filename)
  
  # Local file path
  local_path=file.path("./all_data", filename)
  
  # Skip if already downloaded
  
  # Try to download
  tryCatch({
    message("Downloading: ", filename)
    
    # Use httr for more robust downloading
    response=GET(url, write_disk(local_path, overwrite = TRUE), 
                    timeout(300), progress())
    
    if(response$status_code == 200) {
      message("Success: ", filename)
      return(local_path)
    } else {
      message("Failed (HTTP ", response$status_code, "): ", filename)
      file.remove(local_path)  # Remove partial file
      return(NULL)
    }
  }, error = function(e) {
    message("Error downloading ", filename, ": ", e$message)
    if(file.exists(local_path)) file.remove(local_path)
    return(NULL)
  })
}
mclapply(seq(-200, 20, by = 1),download_chelsa_file,mc.cores = 40)

download_chelsa_file=function(time_step) {
  time_str=sprintf("%04d", time_step)  # e.g., "-0210"
  filename=paste0("CHELSA_TraCE21k_bio12_", time_str, "_V.1.0.tif")
  
  base_url="https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/bioclim/bio12/"
  url=paste0(base_url, filename)
  
  # Local file path
  local_path=file.path("./all_data", filename)
  
  # Skip if already downloaded
  
  # Try to download
  tryCatch({
    message("Downloading: ", filename)
    
    # Use httr for more robust downloading
    response=GET(url, write_disk(local_path, overwrite = TRUE), 
                    timeout(300), progress())
    
    if(response$status_code == 200) {
      message("Success: ", filename)
      return(local_path)
    } else {
      message("Failed (HTTP ", response$status_code, "): ", filename)
      file.remove(local_path)  # Remove partial file
      return(NULL)
    }
  }, error = function(e) {
    message("Error downloading ", filename, ": ", e$message)
    if(file.exists(local_path)) file.remove(local_path)
    return(NULL)
  })
}
mclapply(seq(-200, 20, by = 1),download_chelsa_file,mc.cores = 40)


database=read.table("../Data/data_sites.csv",sep=";")
sites=data.frame(
  site_id = database$File_ID,
  lon = database$Longitude,    
  lat = database$Lat
)

extract_raster_values=function(raster_path, sites_df) {
  tryCatch({
    rast_obj=rast(raster_path)
    crs(rast_obj) = "EPSG:4326" 
    sites_vect=vect(sites_df, geom = c("lon", "lat"), crs = "EPSG:4326")
    values=terra::extract(rast_obj, sites_vect)
    
    result=data.frame(
      site_id = sites_df$site_id,
      value = values[, 2]
    )
    
    return(result)
  }, error = function(e) {
    message("Error processing ", basename(raster_path), ": ", e$message)
    return(data.frame(
      site_id = sites_df$site_id,
      value = NA
    ))
  })
}

parse_chelsa_metadata=function(filename) {
  basename=basename(filename)
  
  # Parse variable
  var_match=regmatches(basename, regexpr("bio[0-9]+", basename))
  variable=ifelse(length(var_match) > 0, var_match[1], NA)
  
  parts=unlist(strsplit(basename, "_"))
  time_number=NA
  
  if (length(parts) >= 4) {
    time_str=parts[4]
    time_number=as.integer(time_str)
  }
  
  if (is.na(time_number)) {
    time_match=regmatches(basename, regexpr("_-?[0-9]+_", basename))
    if (length(time_match) > 0) {
      time_number=as.integer(gsub("_", "", time_match[1]))
    }
  }
  
  if (!is.na(time_number)) {
    calendar_year=time_number * 100
    years_bp=1950 - calendar_year
    
    return(list(
      variable = variable,
      file_number = time_number,
      calendar_year = calendar_year, 
      years_bp = years_bp            
    ))
  }
  
  return(list(
    variable = variable,
    file_number = NA,
    calendar_year = NA,
    years_bp = NA
  ))
}


message("Extracting values for ", nrow(sites), " sites...")

extraction_results=list()
file_list=list.files("./all_data/", full.names = TRUE)

for(file_path in file_list) {
  filename=basename(file_path)
  file_meta=parse_chelsa_metadata(file_path)
  year_bp=file_meta$years_bp
  variable=file_meta$variable
  calendar_year_=file_meta$calendar_year
  
  if(is.na(year_bp) || is.na(variable)) {
    message("Could not parse metadata from: ", filename)
    next
  }
  
  extracted=extract_raster_values(file_path, sites)
  extracted$year_bp=year_bp
  extracted$variable=variable
  extracted$calendar_year=calendar_year_
  
  extraction_results[[file_path]]=extracted
}

all_data=bind_rows(extraction_results)

final_list=list()

for(var in c("bio01", "bio12", "bio03")) {
  var_data=all_data %>%
    filter(variable == var) %>%
    select(site_id, year_bp, calendar_year, value) %>%
    arrange(year_bp)
  
  final_list[[var]]=as.data.frame(var_data)
  write.table(var_data, paste0("chelsa_", var, "_timeseries.csv"), sep=";", row.names=FALSE)
}


# aggregating

d_clim=rbind(read.table("./Data/chelsa_bio01_timeseries.csv",sep=";",header = T)%>%mutate(., Clim_ID="Temperature",value=value-273.15),
             read.table("./Data/chelsa_bio03_timeseries.csv",sep=";",header = T)%>%mutate(., Clim_ID="Isothermality"),
             read.table("./Data/chelsa_bio12_timeseries.csv",sep=";",header = T)%>%mutate(., Clim_ID="Precipitation"))

little_ice_age=c(1300, 1850)    
medieval_warm_period=c(950, 1250)  
LGM=c(-24000, -18000)
Holo=c(-6500, -5500) 

p=ggplot(d_clim%>%dplyr::filter(.,site_id %in% unique(.$site_id)[c(1,75,115)]))+
  geom_line(aes(x=calendar_year,y=value))+
  the_theme2+labs(y="")+
  facet_wrap(Clim_ID~site_id,scales="free",strip.position="left")+
  annotate("rect", 
           xmin = little_ice_age[1], xmax = little_ice_age[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "blue", color = NA) +
  annotate("rect", 
           xmin = medieval_warm_period[1], xmax = medieval_warm_period[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "red", color = NA) +
  annotate("rect", 
           xmin = LGM[1], xmax = LGM[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "brown", color = NA) +
  annotate("rect", 
           xmin = Holo[1], xmax = Holo[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "lightgreen", color = NA)+
  theme(strip.placement = "outside")


ggsave("./Figures/SI/Sites_Chelsa.pdf",p,width = 11,height = 9)

p=ggplot(d_clim%>%
         dplyr::group_by(., calendar_year,Clim_ID)%>%
         dplyr::summarise(., value_mean=mean(value),.groups = "keep"))+
  geom_line(aes(x=calendar_year,y=value_mean))+
  facet_wrap(.~Clim_ID,scales="free",strip.position="left")+
  the_theme2+labs(y="Mean annual temperature")+
  annotate("rect", 
           xmin = little_ice_age[1], xmax = little_ice_age[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "blue", color = NA) +
  annotate("rect", 
           xmin = medieval_warm_period[1], xmax = medieval_warm_period[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "red", color = NA) +
  annotate("rect", 
           xmin = LGM[1], xmax = LGM[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "brown", color = NA) +
  annotate("rect", 
           xmin = Holo[1], xmax = Holo[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "lightgreen", color = NA)+
  theme(strip.placement = "outside")

ggsave("./Figures/SI/Aggregated_sites_Chelsa.pdf",p,width = 12,height = 4)



periods=list(
  "LGM" = c(-21000, -19000),          
  "mid_Holocene" = c(-7000, -5000),   
  "medieval_warm" = c(950, 1350),     
  "little_ice_age" = c(1350, 1850),   
  "current" = c(1950, 2000)           
)

all_data_with_periods=all_data %>%
  mutate(period = case_when(
    calendar_year >= periods$LGM[1] & calendar_year <= periods$LGM[2] ~ "LGM",
    calendar_year >= periods$mid_Holocene[1] & calendar_year <= periods$mid_Holocene[2] ~ "mid_Holocene",
    calendar_year >= periods$medieval_warm[1] & calendar_year <= periods$medieval_warm[2] ~ "medieval_warm",
    calendar_year >= periods$little_ice_age[1] & calendar_year <= periods$little_ice_age[2] ~ "little_ice_age",
    calendar_year >= periods$current[1] & calendar_year <= periods$current[2] ~ "current",
    TRUE ~ "other"
  )
  )


period_means=all_data_with_periods %>%
  filter(!is.na(period)) %>%
  group_by(site_id, variable, period) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

result_wide=period_means %>%
  pivot_wider(
    id_cols = site_id,
    names_from = c(variable, period),
    values_from = mean_value,
    names_sep = "_"
  )
write.table(result_wide,"./Data/All_periods_Chelsea.csv",sep=";")

# adding last interglacial period

database=read.table("../Data/data_sites.csv",sep=";")
sites=data.frame(
  site_id = database$File_ID,
  lon = database$Longitude,    
  lat = database$Lat
)

rast_T=rast("./Data/Interglacial_Temp.tiff")
rast_I=rast("./Data/Interglacial_Iso.tiff")
rast_P=rast("./Data/Interglacial_Precip.tiff")
sites_vect=vect(sites, geom = c("lon", "lat"), crs = "EPSG:4326")
d_interglacial=data.frame(Site_ID=database$File_ID)%>%
  mutate(., 
         MAP_Interglacial=terra::extract(rast_P, sites_vect)[,2],
         MAT_Interglacial=terra::extract(rast_T, sites_vect)[,2],
         Isothermality_Interglacial=terra::extract(rast_I, sites_vect)[,2])
write.table(d_interglacial,"./Data/Values_interglacial.csv",sep=";")


#Correlation


d=read.table("./Data/data_sites_CLIM.csv",sep=";")
d1=read.table("./Data/Values_interglacial.csv",sep=";")
colnames(d1)[1]="File_ID"

d2=read.table("./Data/All_periods_Chelsea.csv",sep=";")%>%
  dplyr::select(., -bio01_other,-bio03_other,-bio12_other)%>%
  dplyr::mutate(.,
                bio01_little_ice_age=bio01_current-bio01_little_ice_age,
                bio01_medieval_warm=bio01_current-bio01_medieval_warm,
                #
                bio03_little_ice_age=bio03_current-bio03_little_ice_age,
                bio03_medieval_warm=bio03_current-bio03_medieval_warm,
                #
                bio12_little_ice_age=bio12_current-bio12_little_ice_age,
                bio12_medieval_warm=bio12_current-bio12_medieval_warm
  )

colnames(d2)[1]="File_ID"

d=merge(d,d2,by="File_ID")%>%
  merge(.,d1,by="File_ID")%>%
  dplyr::mutate(., 
                MAP_Interglacial=MAP_Current-MAP_Interglacial,
                MAT_Interglacial=MAT_Current-MAT_Interglacial/10,
                Isothermality_Interglacial=Isothermality_Current-Isothermality_Interglacial)

colnames(d)[grep("little",colnames(d))]=c("MAT_Little_Ice_age","Isothermality_Little_Ice_age","MAP_Little_Ice_age")

p1=Plot_correlation_variables(d[,c("MAT_Current","MAP_Current","Isothermality_Current","MAT_Holocene","MAP_Holocene","Isothermality_Holocene",
                                "MAT_LGM","MAP_LGM","Isothermality_LGM",
                                "MAT_Little_Ice_age","Isothermality_Little_Ice_age","MAP_Little_Ice_age","MAP_Interglacial","MAT_Interglacial",
                                "Isothermality_Interglacial")])
p2=ggplot(d, aes(x =MAT_Little_Ice_age, y = MAP_Little_Ice_age)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightpink") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  labs(
    x = "Change temperature since \n little age (current - past)",
    y = "Change precipitation since \n little age (current - past)"
  ) +
  the_theme2

p3=ggplot(d, aes(x = MAP_Little_Ice_age, y = MAP_Current)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightpink") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  labs(
    x = "Change precipitation since \n little age (current - past)",
    y = "Current precipitation" 
  ) +
  the_theme2


p_tot=ggarrange(p1+theme(legend.position = "bottom"),ggarrange(p2,p3,nrow=2),ncol=2,widths = c(1,.6),labels = letters[1:2])
ggsave("./Figures/SI/Climate_medieval.pdf",p_tot,width = 16,height = 11)



d=read.table("./Data/data_sites_CLIM.csv",sep=";")

d2=read.table("./Data/All_periods_Chelsea.csv",sep=";")
colnames(d2)[1]="File_ID"

d=merge(d,d2,by="File_ID")
colnames(d)[grep("little",colnames(d))]=c("MAT_Little_Ice_age","Isothermality_Little_Ice_age","MAP_Little_Ice_age")

p1=Plot_correlation_variables(d[,c("MAT_Current","MAP_Current","Isothermality_Current",
                                   "MAT_Little_Ice_age","Isothermality_Little_Ice_age","MAP_Little_Ice_age")])

ggsave("./Figures/SI/Climate_medieval_vs_current.pdf",p1,width = 5,height = 5)



## >> 2) Analysis ----

id_metric_kept=c(1,2,3,8,12,14,18)
for (ID_PC in 1:2){
  
  d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
    Closer_to_normality(.)%>%
    Perform_PCA_spatial_struc(.,F)
  
  d2=read.table("./Data/All_periods_Chelsea.csv",sep=";")%>%
    dplyr::select(., -bio01_other,-bio03_other,-bio12_other,-bio01_LGM,-bio03_LGM,-bio12_LGM,
                  -bio01_mid_Holocene,-bio03_mid_Holocene,-bio12_mid_Holocene)%>%
    dplyr::mutate(.,
                  bio01_little_ice_age=bio01_current-bio01_little_ice_age,
                  bio01_medieval_warm=bio01_current-bio01_medieval_warm,
                  #
                  bio03_little_ice_age=bio03_current-bio03_little_ice_age,
                  bio03_medieval_warm=bio03_current-bio03_medieval_warm,
                  #
                  bio12_little_ice_age=bio12_current-bio12_little_ice_age,
                  bio12_medieval_warm=bio12_current-bio12_medieval_warm
    )
  
  colnames(d2)[1]="File_ID"
  
  d=merge(d,d2,by="File_ID")
  
  d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
  d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
  
  d=d%>%
    dplyr::mutate(.,
                  Desert_Current=as.numeric(Desert_Current=="Desert"),
                  Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
    dplyr::select(.,-Biome_LGM,-Biome_Current)%>%
    dplyr::mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
  
  d$Struct1=d[,paste0("Struct",ID_PC)]
  
  #LGM clim
  d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))[c(id_metric_kept,20)]])]
  #Holocene clim
  d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
  d_medieval=d[,c("bio01_little_ice_age","bio03_little_ice_age","bio12_little_ice_age")]
  
  #Current clim
  d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept,20)]])]
  
  Var_part=varpart(d[,c("Struct1")],d_LGM,d_Holo,d_medieval,d_Current)
  
  d_partition_var=tibble(
    R2=c(Var_part$part$indfract$Adj.R.square[c(1:4)],
         1-Var_part$part$indfract$Adj.R.square[length(Var_part$part$indfract$Adj.R.square)]-sum(Var_part$part$indfract$Adj.R.square[1:4])),
    Name=c("Holocene","LGM","Little Ice Age","Current","Other shared variance"),
  )
  d_partition_var$With_R2=paste0(d_partition_var$Name," = ",round(d_partition_var$R2,2))
  d_partition_var$ID=1
  if (any(round(d_partition_var$R2,2)==0 | d_partition_var$R2<0)){
    d_partition_var=d_partition_var[-which(round(d_partition_var$R2,2)==0 | d_partition_var$R2<0),]
  }
  
  d_partition_var$Order=order(d_partition_var$R2)
  d_partition_var=mutate(d_partition_var,Name=fct_reorder(Name,Order,.desc = T))
  d_partition_var=arrange(d_partition_var,Order)
  d_partition_var$Cumulated_R2=sapply(1:nrow(d_partition_var),
                                      function(x){
                                        if (x==1){
                                          return(d_partition_var$R2[x]/2)
                                        }else if (x==nrow(d_partition_var)){
                                          return(sum(d_partition_var$R2)-d_partition_var$R2[x]/2)
                                        }else{
                                          return(sum(d_partition_var$R2[1:(x-1)])+d_partition_var$R2[x]/2)
                                        }
                                      })
  
  
  p_tot=ggplot(d_partition_var, aes(x = ID, fill = Name,label=round(R2,2),y=R2)) +
    geom_bar(stat="identity")+
    geom_label(aes(x=1,y=Cumulated_R2,label=round(R2,2)))+
    scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Shared past clim."="#DAB2BA",
                               "Current"="pink","Little Ice Age"="#C1E0B8","Other shared variance"="#E8E8E8"))+
    the_theme2+labs(fill="")+
    theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())+
    ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))+
    guides(fill = guide_legend(nrow = 2)) 
  
  
  assign(paste0("p",ID_PC),p_tot)
}  

pvariance1=ggarrange(p1+ggtitle("PC1: mean patch-size, cover")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "none"),
                    p2+ggtitle("PC2: spatial aggregation")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "bottom"),align = "hv")



id_metric_kept=c(1,2,3,8,12,14,18)
for (ID_PC in 1:2){
  
  d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
    Closer_to_normality(.)%>%
    Perform_PCA_spatial_struc(.,F)
  
  d2=read.table("./Data/All_periods_Chelsea.csv",sep=";")%>%
    dplyr::select(., -bio01_other,-bio03_other,-bio12_other,-bio01_LGM,-bio03_LGM,-bio12_LGM,
                  -bio01_mid_Holocene,-bio03_mid_Holocene,-bio12_mid_Holocene)%>%
    dplyr::mutate(.,
                  bio01_little_ice_age=bio01_current-bio01_little_ice_age,
                  bio01_medieval_warm=bio01_current-bio01_medieval_warm,
                  #
                  bio03_little_ice_age=bio03_current-bio03_little_ice_age,
                  bio03_medieval_warm=bio03_current-bio03_medieval_warm,
                  #
                  bio12_little_ice_age=bio12_current-bio12_little_ice_age,
                  bio12_medieval_warm=bio12_current-bio12_medieval_warm
    )
  
  colnames(d2)[1]="File_ID"
  
  d=merge(d,d2,by="File_ID")
  
  d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
  d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
  
  d=d%>%
    dplyr::mutate(.,
                  Desert_Current=as.numeric(Desert_Current=="Desert"),
                  Desert_LGM=as.numeric(Desert_LGM=="Desert"))%>%
    dplyr::select(.,-Biome_LGM,-Biome_Current)%>%
    dplyr::mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
  
  d$Struct1=d[,paste0("Struct",ID_PC)]
  
  #LGM clim
  d_LGM=d[,c(colnames(d)[grep("LGM",colnames(d))[c(id_metric_kept,20)]])]
  #Holocene clim
  d_Holo=d[,c(colnames(d)[grep("Holocene",colnames(d))[c(id_metric_kept)]])]
  d_medieval=d[,c("bio01_medieval_warm","bio03_medieval_warm","bio12_medieval_warm")]
  
  #Current clim
  d_Current=d[,c(colnames(d)[grep("Current",colnames(d))[c(id_metric_kept,20)]])]
  
  Var_part=varpart(d[,c("Struct1")],d_LGM,d_Holo,d_medieval,d_Current)
  
  d_partition_var=tibble(
    R2=c(Var_part$part$indfract$Adj.R.square[c(1:4)],
         1-Var_part$part$indfract$Adj.R.square[length(Var_part$part$indfract$Adj.R.square)]-sum(Var_part$part$indfract$Adj.R.square[1:4])),
    Name=c("Holocene","LGM","Medieval Warm Period","Current","Shared var."),
  )
  d_partition_var$With_R2=paste0(d_partition_var$Name," = ",round(d_partition_var$R2,2))
  d_partition_var$ID=1
  if (any(round(d_partition_var$R2,2)==0 | d_partition_var$R2<0)){
    d_partition_var=d_partition_var[-which(round(d_partition_var$R2,2)==0 | d_partition_var$R2<0),]
  }
  
  d_partition_var$Order=order(d_partition_var$R2)
  d_partition_var=mutate(d_partition_var,Name=fct_reorder(Name,Order,.desc = T))
  d_partition_var=arrange(d_partition_var,Order)
  d_partition_var$Cumulated_R2=sapply(1:nrow(d_partition_var),
                                      function(x){
                                        if (x==1){
                                          return(d_partition_var$R2[x]/2)
                                        }else if (x==nrow(d_partition_var)){
                                          return(sum(d_partition_var$R2)-d_partition_var$R2[x]/2)
                                        }else{
                                          return(sum(d_partition_var$R2[1:(x-1)])+d_partition_var$R2[x]/2)
                                        }
                                      })
  
  
  p_tot=ggplot(d_partition_var, aes(x = ID, fill = Name,label=round(R2,2),y=R2)) +
    geom_bar(stat="identity")+
    geom_label(aes(x=1,y=Cumulated_R2,label=round(R2,2)))+
    scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Shared past clim."="#DAB2BA",
                               "Current"="pink","Medieval Warm Period"="#C1E0B8","Shared var."="#E8E8E8"))+
    the_theme2+labs(fill="")+
    theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())+
    ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))+
    guides(fill = guide_legend(nrow = 2)) 
  
  
  assign(paste0("p",ID_PC),p_tot)
}  

pvariance2=ggarrange(p1+ggtitle("PC1: mean patch-size, cover")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "none"),
                    p2+ggtitle("PC2: spatial aggregation")+
                      labs(y="% of explained variance")+
                      theme(legend.position = "bottom"),align = "hv")
ggsave("./Figures/SI/Medieval_period.pdf",ggarrange(pvariance1+ggtitle("Adding Little Ice Age period"),
                                                pvariance2+ggtitle("Adding Medieval Warm Period"),nrow=2,labels = letters[1:2]),
       width = 8,height = 10)








d=read.table("./Data/data_sites_CLIM.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  Perform_PCA_spatial_struc(.)%>%
  dplyr::mutate(., Sand=bcPower(.$Sand,6))%>%
  dplyr::mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))
d=Is_biome_desert(d,name_variable = "Biome_Current","Desert_Current")
d=Is_biome_desert(d,name_variable = "Biome_LGM","Desert_LGM")
d=d%>%
  dplyr::mutate(.,
                Desert_Current=as.factor(as.numeric(Desert_Current=="Desert")),
                Desert_LGM=as.factor(as.numeric(Desert_LGM=="Desert")))%>%
  dplyr::select(.,-Biome_LGM,-Biome_Current)

d2=read.table("./Data/All_periods_Chelsea.csv",sep=";")%>%
  dplyr::select(., -bio01_other,-bio03_other,-bio12_other)%>%
  dplyr::mutate(.,
                bio01_little_ice_age=bio01_current-bio01_little_ice_age,
                bio01_medieval_warm=bio01_current-bio01_medieval_warm,
                bio01_mid_Holocene=bio01_current-bio01_mid_Holocene,
                bio01_LGM=bio01_current-bio01_LGM,
                #
                bio03_little_ice_age=bio03_current-bio03_little_ice_age,
                bio03_medieval_warm=bio03_current-bio03_medieval_warm,
                bio03_mid_Holocene=bio03_current-bio03_mid_Holocene,
                bio03_LGM=bio03_current-bio03_LGM,
                #
                bio12_little_ice_age=bio12_current-bio12_little_ice_age,
                bio12_medieval_warm=bio12_current-bio12_medieval_warm,
                bio12_mid_Holocene=bio12_current-bio12_mid_Holocene,
                bio12_LGM=bio12_current-bio12_LGM
  )%>%
  dplyr::mutate(across(colnames(dplyr::select_if(., is.numeric)),~ (. - mean(.,na.rm=T)) / sd(.,na.rm = T)))

colnames(d2)[1]="File_ID"

d=merge(d,d2,by="File_ID")

library(piecewiseSEM)


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



