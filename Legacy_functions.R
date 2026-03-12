x=c("tidyverse","reshape2","car","lme4","DHARMa","jtools","raster","terra",
    "sf","missMDA","FactoMineR","factoextra","vegan","randomForestSRC",
    "R.utils","spatialwarnings","ggpubr","treemapify","gbm","psych","semEff")
lapply(x, require, character.only = TRUE)
#lapply(x, install.packages, character.only = TRUE)


`%!in%` = Negate(`%in%`)


the_theme2 = theme_classic() + theme(
  legend.position = "bottom",
  panel.border = element_rect(colour = "black", fill=NA),
  strip.background = element_rect(fill = "transparent",color="transparent"),
  strip.text.y = element_text(size = 10, angle = -90),
  strip.text.x = element_text(size = 10),title = element_text(size=8),
  axis.title.y=element_text(size = 10),
  axis.title.x=element_text(size = 10),
  #legend.box="vertical",
  legend.text = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)
Increase_size_axes=function(p,size_title=12,size_text=12){
  return(p+theme(
    strip.text.y = element_text(size = size_title, angle = -90),
    strip.text.x = element_text(size = size_title),title = element_text(size=8),
    axis.title.y=element_text(size = size_title),
    axis.title.x=element_text(size = size_title),
    axis.text.y=element_text(size = size_text),
    axis.text.x=element_text(size = size_text),
    #legend.box="vertical",
    legend.text = element_text(size = size_title), text = element_text(family = "NewCenturySchoolbook")
  ))
}

is_signif=function(x){
  if (x<=.001){
    return("*")
  }else if (x>.001 & x<=.01){
    return("*")
  }else if (x>.01 & x<.05){
    return("*")
  }else if (x>.05 & x<.1){
    return("°")
  }else {
    return("")
  }
}

Plot_correlation_variables=function(d){
  
  library(psych)
  corr_pred=corr.test(d,use = "pairwise.complete.obs",adjust = "none")
  
  corr_pred$r=round(corr_pred$r,2)
  corr_pred$r[lower.tri(corr_pred$r)]=NA
  diag(corr_pred$r)=NA
  
  
  p=ggplot(corr_pred$r%>%
             melt(.)%>%
             add_column(., pval=sapply(1:nrow(.), function(x){
               if (is.na(.$value[x])){
                 return(NA)
               }else{return(melt(corr_pred$p)$value[x])}
             })))+
    geom_tile(aes(x=Var1,Var2,fill=ifelse(pval<.05,value,0)))+
    theme_classic()+
    geom_text(aes(x=Var1,Var2,label=ifelse(pval<.05,round(value,2),"X")),size=3)+
    scale_fill_gradient2(low="red",mid="white",high = "blue",midpoint = 0,na.value = "white")+
    theme(axis.text.x = element_text(angle=60,hjust=1))+
    labs(x="",y="",fill="")
  return(p)
}

Plot_distributions_variables=function(d){
  
  p=ggplot(d%>%
             melt(.))+
    geom_histogram(aes(x=value),fill="grey")+
    the_theme2+
    facet_wrap(.~variable,scales = "free")
    
  return(p)
}


## Spatial statistics ----

Get_sumstat=function(landscape,log_=T,slope=0,compute_KS=T){
  
  
  if (any(landscape==1 |landscape==T)){
    
    cover = sum(landscape) / (dim(landscape)[1]**2)
    
    # number of neighbors
    #vegetation clustering
    neighbors_mat = simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    mean_nb_neigh = mean(neighbors_mat[which(landscape == 1)]) #mean number of plant neighbors
    mean_clustering = mean_nb_neigh / cover
    spatial_ews = generic_sews(landscape>0,4,moranI_coarse_grain = T)$value
    
    spectral_ratio = as.data.frame(spectral_sews(landscape>0,quiet=T))$value
    
    psd=spatialwarnings::patchdistr_sews(landscape>0)
    max_patchsize=max(psd$psd_obs)
    cv_patch=sd(psd$psd_obs)/mean(psd$psd_obs)
    PLR=spatialwarnings::raw_plrange(landscape>0)
    
    Small_patches=table(patchsizes(landscape>0))[1]
    Small_patches_fraction=table(patchsizes(landscape>0))[1]/sum(table(patchsizes(landscape>0)))
    
    fit_psd=safe_psd_types(psd$psd_obs)
    
    if (log_){
      mean_clustering=log(mean_clustering)
      spectral_ratio=log(spectral_ratio)
      max_patchsize=log(max_patchsize/length(landscape))
    }
    
    if (compute_KS){
      if(length(psd$psd_obs>0)){
        ks_dist  = Get_KS_distance(landscape>0,n_shuffle = 199)
      }else{
        ks_dist  = NA
      }
    }else{
      ks_dist=NA
    }
    
    #flow length
    flow_length=flowlength_sews(landscape>0,slope = slope,
                                cell_size = 50/sqrt(dim(landscape)[1]*dim(landscape)[2]))
    
    #power relationship between area and perimeter
    beta=lsm_c_pafrac(raster(landscape), directions = 8, verbose = TRUE)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #mean perimeter-area ratio 
    mean_perim_area=lsm_c_para_mn(raster(landscape), directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #median, mean, sd of patch size (not in pixel unit but in m? to account for differences in resolutions)
    psd_scaled=lsm_p_area(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    mean_psd=mean(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    median_psd=median(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    sd_psd=sd(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    
    #core area metric
    core_area=lsm_c_cai_mn(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    core_area_land=lsm_c_cpland(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #division of patches
    division=lsm_c_division(raster(landscape), directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #fractal dimension 
    fractal_dim=lsm_c_frac_mn(raster(landscape), directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #contig
    contig=lsm_c_contig_mn(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #shape
    Shape_metric=lsm_c_shape_mn(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #All complexity measures at the landscape scale
    
    complex_land=calculate_lsm(raster(landscape), 
                               what = c("lsm_l_ent", "lsm_l_condent", "lsm_l_joinent","lsm_l_mutinf","lsm_l_relmutinf"),
                               full_name = TRUE,directions = 8,neighbourhood = 4)%>%
      dplyr::select(., value,name)
    
    Hx=lsm_l_ent(raster(landscape),neighbourhood = 4)%>%pull(., value)
    
    
    d=tibble(rho_p=cover,
             nb_neigh=mean_nb_neigh,clustering=mean_clustering,
             skewness=spatial_ews[2],variance=spatial_ews[1],moran_I=spatial_ews[3],
             Spectral_ratio=spectral_ratio,PLR=PLR,PL_expo=fit_psd["slope_best"],cv_psd=cv_patch,
             fmax_psd=max_patchsize,cutoff=fit_psd["cutoff_tpl"],slope_tpl=fit_psd["slope_tpl"],
             flow_length=flow_length$value, #flow length
             perim_area_scaling=beta, #scaling power relationship between area and perimeter 
             mean_perim_area=mean_perim_area, #mean perim/area 
             mean_psd=mean_psd, #mean patch size
             median_psd=median_psd, #median patch size
             Small_patches=Small_patches, #number of smaller patches
             Small_patches_fraction=Small_patches_fraction, #number of smaller patches
             sd_psd=sd_psd, #sd patch size
             KS_dist=ks_dist, #Kolmogorov distance with null expectation
             core_area=core_area, #mean % of core area
             core_area_land=core_area_land,  #same but at the landscape scale
             division=division, #how much patches are divided or constitute big patches
             fractal_dim=fractal_dim, #fractal dimension
             contig=contig, #mean connectedness of cells in patches
             Shape_metric=Shape_metric,#shape of vegetation patches
             Cond_H=complex_land$value[1], #conditional entropy
             Shannon_H=complex_land$value[2], #shannon entropy
             Joint_H=complex_land$value[3], #joint entropy
             mutual_inf=complex_land$value[4], #mutual information
             relat_mutual_inf=complex_land$value[5] #relative mutual information
    )
  }else{
    d=tibble(rho_p=0,
             nb_neigh=0,
             clustering=0,
             skewness=0,variance=0,moran_I=0,
             Spectral_ratio=0,PLR=0,PL_expo=0,cv_psd=0,
             fmax_psd=0,cutoff=0,slope_tpl=0,
             flow_length=0,
             perim_area_scaling=0,
             mean_perim_area=0,
             mean_psd=0,
             median_psd=0,
             Small_patches=0,
             Small_patches_fraction=0, 
             sd_psd=0,
             KS_dist=0,
             core_area=0,
             core_area_land=0,
             division=0,
             fractal_dim=0,
             contig=0,
             Shape_metric=0,
             Cond_H=0,
             Shannon_H=0,
             Joint_H=0,
             mutual_inf=0,
             relat_mutual_inf=0
    )
    
    
  }
  
  return(d)
}

Correct_spatialindic_cover=function(d){
  
  list_spatial_indic=c("moran_I","PL_expo","Small_patches",
                       "cv_psd","fmax_psd","mean_psd",
                       "flow_length","variance")
  
  for (stat_id in list_spatial_indic){
    
    d=d%>%melt(.,measure.vars = stat_id)%>%dplyr::select(., -variable)
    
    residual_value=rep(NA,345)
    residual_value[which(!is.na(d$value))]=residuals(lm(d$value~d$Cover,na.action = na.omit))
    d$value = residual_value
    colnames(d)[which(colnames(d)=="value")]=stat_id
  }
  return(d)
}

Get_KS_distance = function(mat,n_shuffle) {
  obs_psd = spatialwarnings::patchsizes(mat)
  # Compute KS distance with null N_SHUFFLE times
  all_ks = unlist(lapply(seq.int(n_shuffle), function(i) {
    null_mat = matrix(sample(mat), nrow = nrow(mat), ncol = ncol(mat))
    null_psd = spatialwarnings::patchsizes(null_mat)
    ks.test(obs_psd, null_psd)[["statistic"]]
  }))
  
  return( mean(all_ks) )
}

safe_psd_types = function(psd) {
  
  if (length(unique(psd)) > 2) {
    
    # Trying fits for xmin = 1
    pl = safe_fit(pl_fit, psd, 1, F, c("plexpo"))
    tpl = safe_fit(tpl_fit, psd, 1, F, c("plexpo", "cutoff"))
    exp = safe_fit(exp_fit, psd, 1, F, c("cutoff"))
    lnorm = safe_fit(lnorm_fit, psd, 1, F, c("meanlog", "sdlog"))
    
    # detect warnings
    warnings = c(
      pl$warn,
      tpl$warn,
      exp$warn,
      lnorm$warn
    )
    
    bics = c(
      pl$bic,
      tpl$bic,
      exp$bic,
      lnorm$bic
    )
    
    comparison = data.frame(
      cbind(
        warn = warnings,
        bic = bics,
        type = 1:4
      )
    )
    # we want the type of the line without warn and with a normal bic
    best = comparison %>%
      dplyr::filter(!warn & bic != -Inf) %>%
      dplyr::filter(bic == min(bic)) %>%
      dplyr::select(type) %>%
      pull
    
    best_2 = comparison %>%
      slice_head(n = 3) %>% # exclude lnorm from comparison
      dplyr::filter(!warn & bic != -Inf) %>%
      dplyr::filter(bic == min(bic)) %>%
      dplyr::select(type) %>%
      pull
    
    # PSD shape parameters
    # we ensure that they are set to NA if they have warnings:
    # we need that every call returns vectors of the same lengths
    
    if (!tpl$warn) {
      slope_tpl  = tpl$plexpo
      cutoff_tpl  = tpl$cutoff
    } else {
      slope_tpl  = NA
      cutoff_tpl  = NA
    }
    
    if (!pl$warn) {
      slope_pl  = pl$plexpo
    } else {
      slope_pl  = NA
    }
    
    if (!is.na(slope_pl) & !is.na(slope_tpl)) {
      if (pl$bic < tpl$bic) {
        slope_best  = slope_pl
      } else {
        slope_best  = slope_tpl
      }
    } else {
      slope_best  = NA
    }
  } else {
    # if psd is too short, fitting anything does not make any sense:
    # return a NA vector
    best = NA
    best_2 = NA
    slope_pl  = NA
    slope_tpl  = NA
    cutoff_tpl = NA
    slope_best = NA
  }
  
  return(c(
    best = best,
    best_2 = best_2,
    slope_pl = slope_pl,
    slope_tpl = slope_tpl,
    cutoff_tpl = cutoff_tpl,
    slope_best = slope_best
  ))
}

safe_fit = function(fit_fun,psd,xmin = 1,bic_only = FALSE,force_output = NULL) {
  
  # define a compute bic function in internal scope
  # (needed for avoiding namespace issues when this is called by any function 
  # wrapped by spw::compute_indicator())
  
  compute_bic_local = function(fit, psd, xmin) {
    
    psd = psd[psd >= xmin] # fitting was done on x >= xmin
    
    return(fit$npars * log(length(psd)) - 2 * fit$ll)
  }
  
  try = myTryCatch(expr = {
    fit_fun(psd, xmin)
  })
  
  
  if (!length(try$error)) {
    fit = try$value
    warn = ifelse(length(try$warning) > 0, 1, 0) # did optim work well ?
    estimates = fit[names(fit) %in% c("plexpo", "cutoff", "meanlog", "sdlog")]
    
    if (any(sapply(estimates, is.nan))) { 
      # in some cases, some NaNs can be produced without warnings, catch them
      # estimates is a list and we want to keep it like that, thus the sapply()
      
      warn = 1
      estimates[sapply(estimates, is.nan)] = NA
    }
    
    bic = ifelse(
      warn,
      -Inf,
      compute_bic_local(fit, psd, xmin)
    )
  } else {
    warn = 1
    bic = -Inf
    estimates = NULL
  }
  
  out = list(
    bic = bic,
    warn = warn
  )
  
  if (!bic_only) {
    
    if (!length(estimates) & length(force_output)) {
      
      estimates = vector(mode = "list", length = length(force_output))
      names(estimates) = force_output
      
    }
    
    out = c(out, estimates)
  }
  
  return(out)
}

myTryCatch=function(expr) {
  #' Catches errors and warnings in evaluation of expr
  #' 
  #' Particularly useful to detect optim problems in PL fit
  #' 
  #' @export
  
  warn=err=NULL
  
  value=withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <= e
      NULL
    }), warning = function(w) {
      warn <= w
      invokeRestart("muffleWarning")
    })
  
  list(value = value, warning = warn, error = err)
}


## For statistical analyses ----

Perform_PCA_spatial_struc=function(df,plot=F,return_PCA=F){
  save=df
  
  #Performing PCA on the spatial structure metrics 
  struct_variables=c("Small_patches","moran_I","mean_psd","fmax_psd","flow_length","core_area","cv_psd")
  #struct_variables=c("Small_patches","moran_I","PL_expo","mean_psd","fmax_psd","flow_length","core_area")
  
  res.comp=imputePCA(df[,which(colnames(df) %in% struct_variables)],ncp=4,scale = T) 
  
  #colnames(res.comp$completeObs)=c("Moran I","PSD slope","Largest \n patch","Bare connec.","Mean \n PSD","# small patch","Shape patch")
  
  if ("completeObs" %in% names(res.comp)){
    colnames(res.comp$completeObs)=c("Moran I","CV PSD","Largest \n patch","Bare connec.","Mean \n PSD","# small patch","Shape patch")
    
    res.pca=PCA(res.comp$completeObs, ncp = 4,  graph=F)
  }else {
    colnames(res.comp)=c("Moran I","CV PSD","Largest \n patch","Bare connec.","Mean \n PSD","# small patch","Shape patch")
    
    res.pca=PCA(res.comp, ncp = 4,  graph=F)
  }
  fit = principal(res.comp, nfactors=2, rotate="varimax")
  
  
  #ploting the PCA
  if (plot){
    print(factoextra::fviz_pca_var(res.pca,col.var="#2746B1")+
            the_theme2+
            ggtitle("")+
            labs(x=paste0("PC 1 (",round(res.pca$eig[1,2],1),")"),
                 y=paste0("PC 2 (",round(res.pca$eig[2,2],1),")")))
  }
  
  if (return_PCA){
    #We extract the first 2 ones (2/3 of the variance) and add them to the data-frame
    save=save%>%
      add_column(.,Struct1=fit$scores[,1],Struct2=fit$scores[,2])
    return(fit)
    
  }else{
    save=save%>%
      add_column(.,Struct1=fit$scores[,1],Struct2=fit$scores[,2])
    return(save)
    
  }
  
}

Get_PCA_first_axes=function(d,names_variables,names_PC="LGM"){
  
  res.pca=PCA(d[,names_variables], ncp = 4,  graph=F)
  
  d$C1=res.pca$ind$coord[,1]
  d$C2=res.pca$ind$coord[,2]
  d$C3=res.pca$ind$coord[,3]
  d$C4=res.pca$ind$coord[,4]
  
  colnames(d)[(ncol(d)-3):ncol(d)]=c(paste0("PC1_",names_PC),paste0("PC2_",names_PC),
                                     paste0("PC3_",names_PC),paste0("PC4_",names_PC))
  return(d)
}

plotPCA= function(d_PCA, 
                  xIndex, yIndex, classical_PCA=F,
                  xLim=c(-2.5, 2.5), yLim=c(-2.5, 2.7), annotateFactor = 2.6,
                  colorLow = "#FBF2F8", colorHigh = "#7D0D5E", colorMid = "#FFB1E9", midPoint = 0.09,
                  pointSize = 1, pointAlpha = .2, pointColor = "grey85", xTimeSegment = 2.5, 
                  yTimeSegment = 2.5, sizeAnnotation = 3.5) {
  
  if (classical_PCA){
    fitScores = d_PCA$ind$coord
    fitLoadings = d_PCA$var$cor
    fitVaccounted = d_PCA$eig
    Var_ax1 = d_PCA$eig[xIndex,2]
    Var_ax2 = d_PCA$eig[yIndex,2]
    
    
  }else{
    fitScores = d_PCA$scores
    fitLoadings = d_PCA$loadings
    fitVaccounted = d_PCA$Vaccounted
    Var_ax1 = 100*d_PCA$Vaccounted[2,xIndex]
    Var_ax2 = 100*d_PCA$Vaccounted[2,yIndex]
    
  }
  
  
  return(ggplot2::ggplot(data = NULL, ggplot2::aes(x = fitScores[, xIndex], y = fitScores[, yIndex])) +
           ggplot2::theme_classic() +
           ggplot2::stat_density_2d(ggplot2::aes(fill = ..level..), geom = "polygon") +
           ggplot2::geom_point(size = pointSize, alpha = pointAlpha, color = "grey") +
           ggplot2::geom_segment(data = NULL, ggplot2::aes(x = 0, y = 0, xend = (fitLoadings[, xIndex] * xTimeSegment),
                                                           yend = (fitLoadings[, yIndex] * yTimeSegment)),
                                 arrow = ggplot2::arrow(length = ggplot2::unit(1 / 2, units = "picas")),
                                 color = "black") +
           ggplot2::annotate("text", label = rownames(fitLoadings), size = sizeAnnotation,
                             x = (fitLoadings[, xIndex] * annotateFactor), y = (fitLoadings[, yIndex] * annotateFactor)) +
           ggplot2::scale_fill_gradient2(low = colorLow, high = colorHigh, mid = colorMid,  midpoint = midPoint) +
           ggplot2::xlim(xLim) +
           ggplot2::ylim(yLim) +
           ggplot2::xlab(paste0("PC ", xIndex, " (",
                                round(x = Var_ax1, digits = 1), "%)")) +
           ggplot2::ylab(paste0("PC ", yIndex, " (",
                                round(x = Var_ax2, digits = 1), "%)")) +
           ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank()) +
           ggplot2::theme(legend.position = "none") +
           ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 14),
                          axis.title.y = ggplot2::element_text(size = 14)))
}

  
Change_name_RF=function(d,name_predictor="predictor"){
  
  if ("Name_stat" %!in% colnames(d)){
    d=d%>%melt(., measure.vars=name_predictor, value.name = "Name_stat")
  }
  
  d=d%>%
    add_column(., Group_variable=sapply(1:nrow(.),function(x){
      
      if (any(grep("Current",.$Name_stat[x]))){
        return("Current climate")
      }else if (any(grep("LGM",.$Name_stat[x]))){
        return("LGM legacy")
      }else if (any(grep("Holocene",.$Name_stat[x]))){
        return("Holocene legacy")
      }else{
        return("Env")
      }
    }))%>%
    mutate(., Name_stat=recode_factor(Name_stat,
                                      "MAT_range_diurn_LGM"="Diurn range, LGM",
                                      "MAPdryQ_LGM"="MAP dryiest Q",
                                      "MAT_LGM"="MAT, LGM",              
                                      "Desert_LGM"="Desert",
                                      "MAP_Holocene"="MAP Holo",
                                      "MAP_LGM"="MAP LGM",
                                      "MAT_Current"="MAT Current",
                                      "MAP_Current"="MAP Current",
                                      "MAP_Seasonality_LGM"="MAP seasonality LGM",              
                                      "MAP_Seasonality_Current"="MAP seasonality Current",              
                                      "MAP_Seasonality_Holocene"="MAP seasonality Holo",              
                                      "MAT_Seasonality_LGM"="MAT seasonality LGM",              
                                      "MAT_Seasonality_Current"="MAT seasonality Current",              
                                      "MAT_Seasonality_Holocene"="MAT seasonality Holo",              
                                      "MATdry_LGM"="MAT dry, LGM",              
                                      "MAPcoldQ_LGM"="MAP coldest Q, LGM",      
                                      "MinMAT_LGM"="Min MAT, LGM",       
                                      "MAT_range_diurn_Holocene"="Diurn range, Holo", 
                                      "MAPdryQ_Holocene"="MAP dryiest Q, Holo",       
                                      "MinMAT_Holocene"="Min MAT, Holo",         
                                      "MATdry_Current"="MAT dry, Current",    
                                      "MAPdryQ_Current"="MAP dryiest Q, Current",        
                                      "MinMAT_Current"="Min MAT, Current",          
                                      "MAPcoldQ_Current"="MAP coldest Q, Current",        
                                      "MAPdry_Current"="MAP dryiest, Current",           
                                      "Isothermality_Current"="Isotherm., Current",  
                                      "MAPwarmQ_Current"="MAP warmest Q, Current",       
                                      "MAT_range_diurn_Current"="Diurn range, Current", 
                                      "Desert_Current"="Desert, Current",      
                                      "Soil_A"="Soil amelior."
    ))
  return(d)
}

Select_variable_VarPart=function(RDA_all=NULL,data,response_var,forwardsel=T,double=F){
  
  if (double){
    if (forwardsel){
      select_variables=adespatial::forward.sel(response_var$Struct1+response_var$Struct2,data)$variables
      return(select_variables)
    }else{
      Select_var_RDA_LGM = ordiR2step(rda(response_var$Struct1+response_var$Struct2 ~ 1, data = data),
                                      scope = formula(RDA_all), R2scope = RsquareAdj(RDA_all)$adj.r.squared,
                                      direction = "both", trace = FALSE)
      return(names(Select_var_RDA_LGM$terminfo$ordered))
    }
  }else{
    if (forwardsel){
      select_variables=adespatial::forward.sel(response_var,data)$variables
      return(select_variables)
    }else{
      Select_var_RDA_LGM = ordiR2step(rda(response_var ~ 1, data = data),
                                      scope = formula(RDA_all), R2scope = RsquareAdj(RDA_all)$adj.r.squared,
                                      direction = "both", trace = FALSE)
      return(names(Select_var_RDA_LGM$terminfo$ordered))
    }
  }
}  

Get_desert_data=function(d_long_lat){

  #First past biomes
  past_biome=read_sf("../Data/Paleo_clim/Desert_past/world_cut.shp")
  matching_id=read.table("../Data/Paleo_clim/Desert_past/lgm_veg.txt",sep = "\t")
  
  all_biomes=tibble()
  for (k in 1:nrow(d_long_lat)){
    target_point = st_sfc(st_point(c(d_long_lat$Longitude[k], d_long_lat$Lat[k])), crs = st_crs(past_biome))
    #Get local info
    nearest_feature_index = st_nearest_feature(target_point, past_biome)
    nearest_feature = past_biome[nearest_feature_index, ]
    name_biome=matching_id$V2[which(matching_id$V1==nearest_feature$VEG_ID)]
    all_biomes=rbind(all_biomes,tibble(ID=k,Past_biome=name_biome))
  }
  
  current_biome=read_sf("../Data/Paleo_clim/Desert_current/Ecoregions2017_corrected.shp")
  #test=st_make_valid(current_biome)
  #st_write(test, "../Data/Paleo_clim/Desert_current/Ecoregions2017_corrected.shp")
  
  all_biomes$Current_biome=NA
  
  
  library(parallel)
  currentbiome=mclapply(1:nrow(d_long_lat),function(k){
    target_point = st_sfc(st_point(c(d_long_lat$Longitude[k], d_long_lat$Lat[k])), crs = st_crs(current_biome))
    #Get local info
    nearest_feature_index = st_nearest_feature(target_point, current_biome)
    nearest_feature = current_biome[nearest_feature_index, ]
    return(tibble(Current_biome=nearest_feature$BIOME_NAME))
  },mc.cores = 50)%>%bind_rows(.)
  
  all_biomes$Current_biome=currentbiome$Current_biome
  return(all_biomes)
}

Change_biome_name_LGM=function(d){
  
  if (all(is.numeric(d$VEG_ID))){
    
    change_names=read.table("../Data/Paleo_clim/Desert_past/lgm_veg.txt",sep = "\t")
    
    d$VEG_ID=unlist(sapply(1:nrow(d),function(x){
      return(change_names$V2[which(change_names$V1==d$VEG_ID[x])])
    }))
    
    d$Save_name=d$VEG_ID
    d$VEG_ID=gsub(" ","",d$VEG_ID)
    d$VEG_ID=gsub("-","",d$VEG_ID)
    
  }
  
  d=d%>%
    mutate(.,VEG_ID=recode_factor(VEG_ID,
                                      "Foreststeppe"="Temperate Grasslands, Savannas & Shrublands",
                                      "Drysteppe"="Temperate Grasslands, Savannas & Shrublands",
                                      "Temperatesteppegrassland"="Temperate Grasslands, Savannas & Shrublands",
                                      "MontaneMosaic"="Montane Grasslands & Shrublands",
                                      "Semi‐arid temperate woodland or scrub"=" Mediterranean Forests, Woodlands & Scrub",
                                      "Subalpineparkland"="Montane Grasslands & Shrublands",
                                      "Temperatedesert"="Deserts & Xeric Shrublands",
                                      "Tropicalextremedesert"="Deserts & Xeric Shrublands",
                                      "Tropicalgrassland"="Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                      "Tropicalsemidesert"="Deserts & Xeric Shrublands",
                                      "Tropicalthornscrubandscrubwoodland"="Tropical & Subtropical Dry Broadleaf Forests",
                                      "Tundra"="Tundra"))
  
  return(d)
}

Is_biome_desert=function(d,name_variable="Biome_Current",new_name="Desert_Current",numeric=T){
  d$Desert=sapply(1:nrow(d),function(x){
    return(ifelse(any(grep("desert",d[x,name_variable])),"Desert",ifelse(any(grep("Desert",d[x,name_variable])),"Desert","Other biome")))
    # return(ifelse(any(grep("desert",d[x,name_variable])),d$Save_name[x],ifelse(any(grep("Desert",d[x,name_variable])),d$Save_name[x],"Other biome")))
  })
  colnames(d)[ncol(d)]=new_name
  
  #d[,ncol(d)]=as.numeric()
  
  return(d)
}

Plot_variance_partition=function(Var_part,merge_past=F,four_groups=T,bar=F){
  
  if (bar){ #CUmulated bar
    
    if (four_groups){
      if (merge_past==F){
        d_partition_var=tibble(
          R2=c(Var_part$part$indfract$Adj.R.square[c(1:5)],
               1-Var_part$part$indfract$Adj.R.square[16]-sum(Var_part$part$indfract$Adj.R.square[1:5])),
          Name=c("Holocene","LGM","Current","Env","Shared past clim.","Other shared variance"),
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
                                     "Current"="pink","Env"="#C1E0B8","Other shared variance"="#E8E8E8"))+
          the_theme2+labs(fill="")+
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank())+
          ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))
        
        
        
      }else{
        d_partition_var=tibble(
          R2=c(sum(Var_part$part$indfract$Adj.R.square[c(1,2,5)]),
               Var_part$part$indfract$Adj.R.square[c(3,4)],
               1-Var_part$part$indfract$Adj.R.square[16]-sum(Var_part$part$indfract$Adj.R.square[1:5])),
          Name=c("Shared past clim.","Current","Env","Other shared variance"),
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
          scale_fill_manual(values=c("Shared past clim."="#DAB2BA",
                                     "Current"="pink","Env"="#C1E0B8","Other shared variance"="#E8E8E8"))+
          the_theme2+labs(fill="")+
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank())+
          ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))
        
        
      }
      
    }else{ # 3groups
      
      if (merge_past==F){
        d_partition_var=tibble(
          R2=c(Var_part$part$indfract$Adj.R.square[c(1:4)],
               1-Var_part$part$indfract$Adj.R.square[8]-sum(Var_part$part$indfract$Adj.R.square[1:4])),
          Name=c("Holocene","LGM","Current","Shared past clim.","Other shared variance"),
        )
        d_partition_var$With_R2=paste0(d_partition_var$Name," = ",round(d_partition_var$R2,2))
        d_partition_var$ID=1
        if (any(round(d_partition_var$R2,2)==0 | d_partition_var$R2<0)){
          d_partition_var=d_partition_var[-which(round(d_partition_var$R2,2)==0 | d_partition_var$R2<0),]
        }
        
        d_partition_var=arrange(d_partition_var,(R2))
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
        d_partition_var$Order=1:nrow(d_partition_var)
        d_partition_var=mutate(d_partition_var,Name=fct_reorder(Name,Order,.desc = T))
        
        
        p_tot=ggplot(d_partition_var, aes(x = ID, fill = Name,label=round(R2,2),y=R2)) +
          geom_bar(stat="identity")+
          geom_label(aes(x=1,y=Cumulated_R2,label=round(R2,2)))+
          scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Shared past clim."="#DAB2BA",
                                     "Current"="pink","Other shared variance"="#E8E8E8"))+
          the_theme2+labs(fill="")+
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank())+
          ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))
        
        
      }else{
        d_partition_var=tibble(
          R2=c(sum(Var_part$part$indfract$Adj.R.square[c(1,2,4)]),
               Var_part$part$indfract$Adj.R.square[c(3)],
               1-Var_part$part$indfract$Adj.R.square[8]-sum(Var_part$part$indfract$Adj.R.square[1:4])),
          Name=c("Shared past clim.","Current","Other shared variance"),
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
          scale_fill_manual(values=c("Shared past clim."="#DAB2BA",
                                     "Current"="pink","Env"="#C1E0B8","Other shared variance"="#E8E8E8"))+
          the_theme2+labs(fill="")+
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank())+
          ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))
      }
    }
  }else{
    if (four_groups){
      if (merge_past==F){
        d_partition_var=tibble(
          R2=c(Var_part$part$indfract$Adj.R.square[c(1:5)],
               1-Var_part$part$indfract$Adj.R.square[16]-sum(Var_part$part$indfract$Adj.R.square[1:5])),
          Name=c("Holocene","LGM","Current","Env","Shared past clim.","Other shared variance"),
        )
        d_partition_var$With_R2=paste0(d_partition_var$Name," = ",round(d_partition_var$R2,2))
        
        p_tot=ggplot(d_partition_var, aes(area = R2, fill = Name,label=round(R2,2))) +
          geom_treemap()+
          geom_treemap_text(colour = "black", place = "centre",size=12)+
          scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Shared past clim."="#DAB2BA",
                                     "Current"="pink","Env"="#C1E0B8","Other shared variance"="#E8E8E8"))+
          the_theme2+labs(fill="")+
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
          ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))
        
        
        
      }else{
        d_partition_var=tibble(
          R2=c(sum(Var_part$part$indfract$Adj.R.square[c(1,2,5)]),
               Var_part$part$indfract$Adj.R.square[c(3,4)],
               1-Var_part$part$indfract$Adj.R.square[16]-sum(Var_part$part$indfract$Adj.R.square[1:5])),
          Name=c("Shared past clim.","Current","Env","Other shared variance"),
        )
        d_partition_var$With_R2=paste0(d_partition_var$Name," = ",round(d_partition_var$R2,2))
        
        p_tot=ggplot(d_partition_var, aes(area = R2, fill = Name,label=round(R2,2))) +
          geom_treemap()+
          geom_treemap_text(colour = "black", place = "centre",size=12)+
          scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Shared past clim."="#DAB2BA",
                                     "Current"="pink","Env"="#C1E0B8","Other shared variance"="#E8E8E8"))+
          the_theme2+labs(fill="")+
          ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))
        
        
      }
      
    }else{ # 3groups
      
      if (merge_past==F){
        d_partition_var=tibble(
          R2=c(Var_part$part$indfract$Adj.R.square[c(1:4)],
               1-Var_part$part$indfract$Adj.R.square[8]-sum(Var_part$part$indfract$Adj.R.square[1:4])),
          Name=c("Holocene","LGM","Current","Shared past clim.","Other shared variance"),
        )
        d_partition_var$With_R2=paste0(d_partition_var$Name," = ",round(d_partition_var$R2,2))
        
        p_tot=ggplot(d_partition_var, aes(area = R2, fill = Name,label=round(R2,2))) +
          geom_treemap()+
          geom_treemap_text(colour = "black", place = "centre",size=12)+
          scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Shared past clim."="#DAB2BA",
                                     "Current"="pink","Other shared variance"="#E8E8E8"))+
          the_theme2+labs(fill="")+
          theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
          ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))
        
        
        
      }else{
        d_partition_var=tibble(
          R2=c(sum(Var_part$part$indfract$Adj.R.square[c(1,2,4)]),
               Var_part$part$indfract$Adj.R.square[c(3)],
               1-Var_part$part$indfract$Adj.R.square[8]-sum(Var_part$part$indfract$Adj.R.square[1:4])),
          Name=c("Shared past clim.","Current","Other shared variance"),
        )
        d_partition_var$With_R2=paste0(d_partition_var$Name," = ",round(d_partition_var$R2,2))
        
        p_tot=ggplot(d_partition_var, aes(area = R2, fill = Name,label=round(R2,2))) +
          geom_treemap()+
          geom_treemap_text(colour = "black", place = "centre",size=12)+
          scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Shared past clim."="#DAB2BA",
                                     "Current"="pink","Other shared variance"="#E8E8E8"))+
          the_theme2+labs(fill="")+
          ggtitle(paste0("Unexplained variance = ",1-round(max(d_partition_var$Cumulated_R2),2)))
      }
    }
    
  }
  return(list(Figure=p_tot,data=d_partition_var))
}

Perform_model_averaging=function(model_avg,delta_AIC=4,plot=T){
  
  get_models_avg=MuMIn::model.avg(object =model_avg, 
                                  subset = delta < delta_AIC)
  model_averaging=summary(object = get_models_avg)$coefmat.subset[c(-1), c(1:2,5)]
  
  model_averaging_CI=confint(get_models_avg)[-1,]
  colnames(model_averaging_CI)=c("q1_95","q3_95")
  model_averaging_CI=model_averaging_CI[order(rownames(model_averaging_CI)),]

  model_averaging_CI_90=confint(get_models_avg,level = .9)[-1,]
  colnames(model_averaging_CI_90)=c("q1_90","q3_90")
  model_averaging_CI_90=model_averaging_CI_90[order(rownames(model_averaging_CI_90)),]
  
  model_averaging=as.data.frame(model_averaging[order(rownames(model_averaging)), ])
  model_averaging=cbind(model_averaging,model_averaging_CI,model_averaging_CI_90)
  model_averaging=model_averaging[order(model_averaging$Estimate), ]
  model_averaging$label=rownames(model_averaging)
  colnames(model_averaging)[1:3]=c("coeff","stdError","Pval")
  
  q2=confint(get_models_avg,level = 1e-9)[-1,1]
  q2=q2[rownames(model_averaging)]
  
  model_averaging$coeff=q2
  
  model_averaging=arrange(model_averaging,q2)
  model_averaging$Order=1:nrow(model_averaging)
  model_averaging=mutate(model_averaging,label=fct_reorder(label,Order))
  
  model_averaging$Past_clim=NA
  if (any(grep("LGM",model_averaging$label))){
    model_averaging$Past_clim[grep("LGM",model_averaging$label)]="LGM"
  }
  if (any(grep("Holocene",model_averaging$label))){
    model_averaging$Past_clim[grep("Holocene",model_averaging$label)]="Holocene"
  }
  if (any(grep("Current",model_averaging$label))){
    model_averaging$Past_clim[grep("Current",model_averaging$label)]="Current"
  }
  
  fig_mod=ggplot2::ggplot(model_averaging, ggplot2::aes(y=label, x=coeff)) + 
    ggplot2::theme_classic()+
    ggplot2::geom_vline(xintercept = 0,linetype="dotted",size = 1)+
    ggplot2::geom_linerange(aes(xmin=q1_90,xmax=q3_90,fill=Past_clim),lwd=1.5)+
    ggplot2::geom_pointrange(ggplot2::aes(xmin=q1_95, xmax=q3_95,fill=Past_clim),size=1,shape=21,color="black",lwd=1)+
    ggplot2::xlab("Model Coefficient")+
    ggplot2::theme(legend.position="none")+ 
    ggplot2::theme(axis.title.y=ggplot2::element_blank())+
    scale_fill_manual(values=c("LGM"="#FFF9BF", "Holocene"="#FFD2A0","Current"="pink"))
  if (plot){
    print(fig_mod)
  }
  return(list(model_averaging=model_averaging,Fig=fig_mod))
}




## Images ----

Extract_binary_matrix=function(row_id,data){  
  img=readJPEG(data$Full_name_image[row_id]) 
  
  #and binarize the kmean output 
  kmean_img=k_means_RGB(img,Nclust)
  
  cats=get_cut_grayscale_values(Nclust)[[Cut]]
  mat=kmean_img %>% binarize(cats[[1]], cats[[2]])
  return(mat)
}

k_means_RGB = function(img, k) {
  
  imgDm = dim(img)
  
  # separate R,G,B
  
  imgRGB = data.frame(
    x = rep(1:imgDm[2], each = imgDm[1]),
    y = rep(imgDm[1]:1, imgDm[2]),
    R = as.vector(img[, , 1]),
    G = as.vector(img[, , 2]),
    B = as.vector(img[, , 3])
  )
  
  newBW = k_means(imgRGB, k, imgDm)
  
  return(newBW)
}

k_means = function(imgBands, k, dm) {
  
  k = as.integer(k)
  imgDm = dm
  kMeans = kmeans(imgBands[, c("R", "G", "B")], centers = k, nstart = 1)
  
  colors2 = c(1, 0)
  colors3 = c(1, 0.5, 0)
  colors4 = c(1, 0.7, 0.3, 0)
  
  colors = list(colors2, colors3, colors4)
  
  kMeansGray = kMeans$centers %>%
    as.data.frame() %>%
    mutate(num = 1:k, gray = 0.2126 * R + 0.7152 * G + 0.0722 * B) %>% # HDTV gray
    arrange(desc(gray))
  
  imgBands.update = cbind(imgBands, clust = kMeans$cluster) # %>% # add cluster info to img
  # ungroup
  
  imgBW = matrix(imgBands.update$clust,
                 nrow = imgDm[1],
                 ncol = imgDm[2]
  )
  
  # convert clusters numbers to gray scale
  
  newBW = imgBW # separate matrix
  for (l in 1:k) { # for each cluster category
    # replace clusters by ascending color (in gray : 1 -> 0)
    newBW[imgBW == kMeansGray$num[l]] = colors[[k - 1]][l]
  }
  
  return(newBW)
}

binarize = function(mat, cat0, cat1) {
  
  mat2 = mat
  mat2[mat %in% cat1] = 1 # full : black so gray is 0
  mat2[mat %in% cat0] = 0 # empty : white so gray is 1
  
  return(mat2)
}

get_cut_grayscale_values = function(nclust) {
  
  if (nclust == 2) {
    scale= c(1, 0)
  } else if (nclust == 3) {
    scale= c(1, 0.5, 0)
  } else if (nclust == 4) {
    scale=c(1, 0.7, 0.3, 0)
  }
  
  
  if (nclust == 2) {
    cut = list(scale)
  } else if (nclust == 3) {
    cut = list(
      list(scale[1], scale[2:3]),
      list(scale[1:2], scale[3])
    )
  } else if (nclust == 4) {
    cut = list(
      list(scale[1], scale[2:4]),
      list(scale[1:2], scale[3:4]),
      list(scale[1:3], scale[4])
    )
  }
  
  return(cut)
}

Closer_to_normality=function(d){
  
  return(
    
    d%>%mutate(., 
               flow_length=bcPower(flow_length,.2),
               mean_psd=bcPower(mean_psd,-.2),
               PL_expo=bcPower(PL_expo+1,5),
               core_area=bcPower(core_area,-.2),
               cv_psd=log(cv_psd),
               Small_patches=bcPower(Small_patches,-.1),
               MAT_LGM=bcPower(MAT_LGM+15,-3.5),
               Isothermality_Current=bcPower(Isothermality_Current,-1)
               )
  )  
  
}
