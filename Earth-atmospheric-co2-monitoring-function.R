###################################################################################
## constant values
###################################################################################

# aggregate variable 
# grid resolution, 1 degree;
# aggregate duration, 1 stands for 1 month;
# moving window parameter, squre window, [lon=5, lat=4] 
aggregate_var<-data.frame( aggregate_res = 1,
                           dur_factor = 1,
                           window_width_lon = 5/2,
                           window_width_lat = 4/2)

# range for Au map
au_lon <- c(105,160)
au_lon_ticks <- seq(110,160,by=10)
au_lat <- c(-50,-0.5)
au_lat_ticks <- c(seq(-50,0,by=10))
au_lon_rg <- c(au_lon[1]-aggregate_var$aggregate_res/2, au_lon[2]+aggregate_var$aggregate_res/2)
au_lat_rg <- c(au_lat[1]-aggregate_var$aggregate_res/2, au_lat[2]+aggregate_var$aggregate_res/2)

lag_num <- 12  # moving year windown

# #generate the coordinates of the centers of aggregation block
# centers pts for northern extratropical region
centrepoints_raster <- raster(xmn=-180 - (aggregate_var$aggregate_res/2),
                              xmx=180  + (aggregate_var$aggregate_res/2),
                              ymn=23   - (aggregate_var$aggregate_res/2),
                              ymx=87   + (aggregate_var$aggregate_res/2), res=aggregate_var$aggregate_res)
centrepoints_raster_area <- raster::area(centrepoints_raster)
centrepoints_raster_area_df <- as.data.frame(centrepoints_raster_area,xy=TRUE) %>% mutate(rel_layer=rescale(layer, from=c(0, max(layer))))
names(centrepoints_raster_area_df) <- c("x","y","area","rel_area")
centrepoints <- as.data.frame(coordinates(centrepoints_raster))
names(centrepoints) <- c("lon","lat")

# centers pts for au region
au_interested_rg <- c(au_lon_rg,au_lat_rg)
centrepoints_au <- as.data.frame(coordinates(raster(xmn=au_interested_rg[1],
                                                    xmx=au_interested_rg[2],
                                                    ymn=au_interested_rg[3],
                                                    ymx=au_interested_rg[4], res=aggregate_var$aggregate_res)))

names(centrepoints_au) <- c("lon","lat")

# define the sub region for NT, ST, NX, SX, global and au,
# use lat 0.5, i.e., center of the 1x1 grid,
# as a margin to avoid double counting problem in practices,
# when processing the boundaries between sub regions
sub_region<-data.frame(northern_tropics=c(-180,180,0,23.5),
                       southern_tropics=c(-180,180,-23.5,0),
                       northern_extratropics=c(-180,180,23.5,90),
                       southern_extratropics=c(-180,180,-90,-23.5),
                       globalearth=c(-180,180,-90,90),
                       tropics=c(-180,180,-23.5,23.5),
                       australia=au_interested_rg)

###################################################################################
## functions
## NOTE: the ave of the interested variables is used in the paper,
##       but, the median is also calculated in the below functions
###################################################################################
##---------------------------------------------------------------------------------
# functions for time series plots
##---------------------------------------------------------------------------------
# Get the aggregation xco2 in 1x1, perform localloess smooth, 
cell_xco2_impute_summary<-function(center,data,lat_width=2,lon_width=2.5, loess_period = 6, loess_degree = 1, startendcheck = 0, thres = 3, fillmode = 0){
  
  # load zoo for interpolation and extrapolation
  require(zoo)
  require(timetk)
  require(xts)
  
  print(paste0("lon = ",toString(center[1]),", lat = ",toString(center[2])))
  
  
  blk_lon_left <- center[1]-lon_width
  blk_lon_right <- center[1]+lon_width
  
  blk_lat_bottom <- center[2]-lat_width
  blk_lat_up <- center[2]+lat_width
  
  if(blk_lon_left < (-180)){
    blk_lon_warping <- c((blk_lon_left+360),180)
    blk_lon_left <- (-180)
  }else if(blk_lon_right > 180){
    blk_lon_warping <- c(-180, (blk_lon_right-360))
    blk_lon_right <- 180
  }else{blk_lon_warping <- NULL}
  
  block_subset <- data %>% dplyr::filter(lon>=blk_lon_left & lon<=blk_lon_right & lat>= blk_lat_bottom & lat<=blk_lat_up)
  
  if (!is.null(blk_lon_warping)){
    block_subset <- bind_rows(block_subset,
                              data %>% dplyr::filter(lon>=blk_lon_warping[1] & lon<=blk_lon_warping[2] & lat>= blk_lat_bottom & lat<=blk_lat_up))
  }
  # process the data to get the x(s,t)
  if (nrow(block_subset) != 0){
      
    block_subset__mean <- block_subset %>% 
      dplyr::group_by(year,month) %>% 
      mutate(cnt=n()) %>% 
      dplyr::group_by(year,month,cnt) %>%
      summarise_at("xco2",mean,na.rm=TRUE) %>% ungroup()
    span_time <- span_time_df()
    roll_factor = 12
    
    # moving local spline interpolation
    cell_monthly_res <-smooth_xco2_localloess_cell(raw_df = block_subset__mean, time_df = span_time,
                                                       roll_factor = roll_factor, loess_period = loess_period, loess_degree = loess_degree,
                                                       startendcheck = startendcheck, thres = thres, fillmode = fillmode )
      
    cell_monthly_res$lon <- center[1]
    cell_monthly_res$lat <- center[2]
    cell_monthly_res <- cell_monthly_res %>% as.data.frame() %>% dplyr::select(lon,lat,datetime,xco2,xco2_spline,bt,mt,anomaly,cnt,sd)
    
  }else{
    cell_monthly_res <- data.frame(lon=center[1],
                                       lat=center[2],
                                       datetime=NA,
                                       xco2=NA,
                                       xco2_spline=NA,
                                       bt=NA,
                                       mt=NA,
                                       anomaly=NA,
                                       cnt=NA,
                                       sd=NA)
    
  }
  return(cell_monthly_res)
}

span_time_df<-function(start=as.Date("2014-09-05", tz="UTC"),end=as.Date("2022-02-28",tz="UTC")){
  # This function generates the full time steps (monthly) for a given period
  span_date <- data.frame(datetime=seq(start,end, by = "day")) %>% 
    mutate(year=year(datetime),month=month(datetime),halfmonth=LETTERS[-9][2*(month(datetime)-1) + 1 + (mday(datetime) > 15)]) 
  span_date <- span_date %>% dplyr::select(year, month) %>% unique() %>% arrange(year, month)
  
  return(span_date)
}

smooth_xco2_localloess_cell<-function(raw_df, time_df, roll_factor = 12, offset = 0, loess_period = 6, loess_degree = 1, startendcheck = 0, thres = 3, fillmode = 0){
  
  # This function uses localloess to perform the interpolation on the time series{x(t),x(t-1),...,x(t-23)}
  # the x(t) or x(t-23) might be NAs, do the extropolation 
  
  time_width<-roll_factor + offset # as t-0,t-1,...,t-roll_factor-1, use offset to adjust the window size
  loess_span<-loess_period/roll_factor
  
  impute_df <- full_join(raw_df,time_df,by=names(time_df)) 

  impute_df <- impute_df %>% 
    dplyr::arrange(year, month) %>%
    mutate(day=1, datetime=make_datetime(year,month,day))
  
  #embed function put the beginning at the last column
  raw_xco2<-embed(rev(impute_df$xco2),time_width)
  raw_xco2_datetime<-rev(impute_df$datetime)[1:dim(raw_xco2)[1]]
  
  impute_xco2<-apply(raw_xco2,1,loess_predict,loess_span=loess_span,loess_degree = 1,startendcheck = 0, thres = thres, fillmode = fillmode) %>% t()
  
  
  bt<-apply(impute_xco2[,(1+offset):time_width],1,mean,na.rm=TRUE)
  
  # xco2_spline, this is the results from the localloess
  xco2_spline<-unlist(lapply(split(impute_xco2,row(impute_xco2)-col(impute_xco2)),mean,na.rm=TRUE)) %>% as.vector() %>% rev()
  bt_df=data.frame(datetime=raw_xco2_datetime,bt=bt)
  
  # join the results
  impute_df<-full_join(impute_df,bt_df,by="datetime") %>%
    mutate(xco2_spline=xco2_spline,bt=roll_meanr(xco2_spline,n = roll_factor,na.rm = FALSE)) %>% 
    mutate(anomaly=xco2_spline-bt,mt=diff.xts(bt)) %>%
    mutate(sd=1/sqrt(cnt)) %>%
    dplyr::select(datetime,xco2,xco2_spline,anomaly,bt,mt,sd,cnt)
  
  # find the index xco2=NA but xco2_spline has value, interpolate the missing value
  sd_NA_ind <- which(is.na(impute_df$xco2) & (!is.na(impute_df$xco2_spline)))
  xco2_not_NA_ind <- which(!is.na(impute_df$xco2))
  for (get_sd_i in sd_NA_ind){
    left_bound_ind <- max(xco2_not_NA_ind[xco2_not_NA_ind < get_sd_i])
    right_bound_ind <- min(xco2_not_NA_ind[xco2_not_NA_ind > get_sd_i])
    impute_df$sd[get_sd_i] <- sqrt((1/(impute_df$cnt[left_bound_ind]) + 1/(impute_df$cnt[right_bound_ind]))/2)
  }
  
  return(impute_df)
  
}

loess_predict<-function(rawdata, loess_span = 0.25, loess_degree = 1, startendcheck = 0,  thres = 3, fillmode = 0){
  # perform loess smooth
  if (fillmode == 0){
    fill_val = rep(NA, length(rawdata))
  }else{ fill_val = rawdata}
  
  NA_num <- cnt_consecutive_NA(rawdata, startendcheck = startendcheck)
  
  if (NA_num >= thres){
    smoothdata<-fill_val
  }else{
    # do the loess
    
    rawdata_df<-data.frame(ind=seq_along(rawdata),val=rawdata)
    
    rawdata_df_na.trim <- rawdata_df %>% zoo::na.trim()
    
    smoothdata<-tryCatch(predict(loess(val~ind,rawdata_df_na.trim,span=loess_span,degree=loess_degree,control = loess.control(surface = "interpolate")),rawdata_df_na.trim$ind,se=FALSE),
                         error=function(e) rawdata_df_na.trim$val)
    smoothdata<-dplyr::na_if(smoothdata, 0)  # the values of xco2 can't be 0, loess returns all 0 if fails
    rawdata_df_na.trim$smoothdata <- smoothdata
    smoothdata_df<-full_join(rawdata_df,rawdata_df_na.trim %>% dplyr::select(ind,smoothdata), by="ind")
    
    # if the start/end has NA values, but not all NAs
    if((sum(is.na(smoothdata_df$smoothdata)) > 0) && (sum(is.na(smoothdata_df$smoothdata)) != dim(smoothdata_df)[1])){
      if(sum(is.na(smoothdata_df$smoothdata))/dim(smoothdata_df)[1] <= 0.3){ # extropolation, otherwise, return all NA
        NA_ind <- is.na(smoothdata_df$smoothdata)
        smoothdata2<-tryCatch(predict(loess(smoothdata~ind,smoothdata_df,span=loess_span,degree=loess_degree,control = loess.control(surface = "direct")),smoothdata_df$ind,se=FALSE),
                              error=function(e) smoothdata_df$smoothdata)
        smoothdata_df$smoothdata[NA_ind] <- smoothdata2[NA_ind]
      }else{
        smoothdata_df$smoothdata<-fill_val
      }
    }
    smoothdata<-smoothdata_df$smoothdata
  }
  return(smoothdata)
}

cnt_consecutive_NA <- function(rawdata, startendcheck = 0){
  # This function counts consecutive NA
  if ((is.na(rawdata[1]) || is.na(rawdata[length(rawdata)])) && startendcheck == 1){
    consecutive_NA_num <- length(rawdata)
  } else{
    # use run length encoding to find the consecutive NA length
    rawdata.isna <- is.na(rawdata)
    consecutiveNA_length <- rle(rawdata.isna)$lengths[rle(rawdata.isna)$values]
    if(length(consecutiveNA_length) == 0){
      # no NA value, return 0
      consecutive_NA_num <- 0
    } else {
      consecutive_NA_num <- max(consecutiveNA_length)
    }
  }
  return(consecutive_NA_num)
}

# draw spatial map
spatial_map<-function(draw_df,
                      plot_title,
                      fill_lim,
                      legend_breaks,
                      legendtitle=TRUE){
  suppressMessages({library(extrafont)})
  extrafont::loadfonts()
  draw_df <- draw_df[complete.cases(draw_df),]

  draw_plot <-(ggplot(draw_df) +
                 theme_Publication() +
                 geom_rect(data=grayrect, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90",color="grey90") +
                 geom_raster(aes(x=x,y=y,fill=value),interpolate = TRUE) +
                 scale_x_continuous(limits = c(-181,181), breaks=seq(-180,180,by=30),expand = c(0,0)) +
                 scale_y_continuous(limits = c(0,90), breaks=seq(0,90,by=10)) +
                 scale_fill_gradientn(colors=c("#0000FF","#ADD8E6","#008000","#FFFF00","#FFA500","#FF0000","#730005"),
                                      oob=squish,na.value="white",limits=fill_lim,breaks = legend_breaks) +
                 labs(x="", y="",fill="ppm",title=plot_title) +
                 theme(legend.position = "right",legend.margin=margin(0, 0, 0, 0), legend.key.size = unit(0.4, "cm"),
                       legend.direction = "vertical",plot.margin=unit(c(0,1,0,1),"mm"),
                       plot.title = element_markdown(size = 20),
                       legend.key.height=unit(2, 'cm')) +
                 coord_quickmap()) %>% draw_world(inc_border=TRUE)

  
  if (legendtitle == FALSE){
    draw_plot <- draw_plot + theme(legend.title=element_blank()) 
  }
  
  return(draw_plot)
}

# calculate the xco2 temporal summary in the provided region 
# To call this function, mclapply could be used to speed up
# if there are multiple regions, for example:
# sub_region_monthly_summary <- mclapply(sub_region, 
#                                 xco2_temporal_summary,
#                                 xco3_df,
#                                 mc.cores = length(sub_region))
# 
xco2_temporal_summary<-function(region,xco2_df){
  get_subset <- xco2_extract(xco2_df, 
                             lon_left = region[1], 
                             lon_right = region[2],
                             lat_bottom = region[3],
                             lat_up = region[4])
  
  # in case there is data missing, generate the full time stamps
  span<-as.period(as.Date(min(get_subset$datetime)) %--% as.Date(max(get_subset$datetime)))
  totalmonths<-abs(span@year)*12 + abs(span@month) + abs((span@day)/31)
  print(paste0("  processing time span: ",span, " ,total months = ",toString(totalmonths)))
  span_date <- data.frame(datetime=seq(as.Date(min(get_subset$datetime)),as.Date(max(get_subset$datetime)), by = "day")) %>% 
    mutate(year=year(datetime),month=month(datetime)) 
  span_date <- span_date %>% dplyr::select(year, month) %>% unique() %>% arrange(year, month)
  
  get_subset_summary <- get_subset %>% 
    group_by_at(vars(names(span_date))) %>% 
    summarise_at(vars(xco2), list(~ mean(.,na.rm=TRUE),~ median(.,na.rm=TRUE))) 
  
  get_subset_summary <- full_join(get_subset_summary,span_date, by=names(span_date))
  get_subset_summary <- get_subset_summary %>% dplyr::arrange(year, month)

  return(as.data.frame(get_subset_summary))
}

# extract a subset by specific duration and region
xco2_extract<-function(xco2_df, startdate, enddate, lon_left, lon_right, lat_bottom, lat_up){
  # requirement of the input: 
  #   1. xco2_df is a dataframe with variables "datetime" "lon" "lat" "xco2"  
  #   2. startdate and enddate specify the temporal dimension,
  #        if given, the data between [startdate, enddate) will be extracted
  #        NOTE: the startdate and enddate are the Date class, time zone is set to UTC
  #   3. lon_left, lon_right, lat_bottom, lat_up specify the spatial dimension, 
  #        if given, the data between lon <- [lon_left, lon_right],
  #        lat <- [lat_bottom,lat_up] will be extracted
  #      NOTE: this function doesn't deal with the warping situation, e.g. lon <- [-170, -180] ∪ [0, 5]
  #      thus, for the warping case, we can call this function twice and then do the rbind.
  
  data_subset <- xco2_df
  
  if (!missing(startdate)){
    if (is.Date(startdate)){
      data_subset <- subset(data_subset,datetime >= startdate)
    }else{
      print("  startdate is not a valid Date class object, can't constrain the start date.")
    }
  }
  
  
  if (!missing(enddate)){
    if (is.Date(enddate)){
      data_subset <- subset(data_subset,datetime < enddate)
    }else{
      print("  enddate is not a valid Date class object, can't constrain the end date.")
    }
  }
  
  if (!missing(lon_left)){
    if (is.numeric(lon_left)){
      data_subset <- subset(data_subset,lon >= lon_left)
    }else{
      print("  lon_left is not a valid numeric value, can't constrain the minimum longitude.")
    }
  }
  
  if (!missing(lon_right)){
    if (is.numeric(lon_right)){
      data_subset <- subset(data_subset,lon <= lon_right)
    }else{
      print("  lon_right is not a valid numeric value, can't constrain the maximum longitude.")
    }
  }
  
  if (!missing(lat_bottom)){
    if (is.numeric(lat_bottom)){
      data_subset <- subset(data_subset,lat >= lat_bottom)
    }else{
      print("  lat_bottom is not a valid numeric value, can't constrain the minimum latitude")
    }
  }
  
  if (!missing(lat_up)){
    if (is.numeric(lat_up)){
      data_subset <- subset(data_subset,lat <= lat_up)
    }else{
      print("  lat_up is not a valid numeric value, can't constrain the maximum latitude")
    }
  }
  
  if(nrow(data_subset) == 0){ print("No dataset available for specified conditions!!")}
  else{
    print(paste0("  Extracted Dataset temporal dimension: ",min(data_subset$datetime)," to ",max(data_subset$datetime)))
    print(paste0("  Extracted Dataset spatial dimension: lon [",toString(min(data_subset$lon)),",",toString(max(data_subset$lon)),
                 "], lat [",toString(min(data_subset$lat)),",",toString(max(data_subset$lat)),"]"))
  }
  return(data_subset)
}

# linear interpolation for missing value
# because there is a missing value in the dataset, the interpolation is performed
# for example, interpolation the results from “sub_region_monthly_summary” calculated from above example
# sub_region_monthly_summary_impute <- lapply(sub_region_monthly_summary, xco2_missingdata_impute)
xco2_missingdata_impute<-function(xco2_df){
  # D(t) = X(t) - x(t-12)
  # 2X(t) <- (D(t-1) + D(t+1))/2 - (D(t+11) + D(t+13))/2 + x(t -12) + x(t+12)
  #          {missingdata_inter1}   {missingdata_inter2}   {missingdata_inter3}
  
  missingdata_ind <-  which(is.na(xco2_df$mean))

  missingdata_inter3 <- xco2_df[c(missingdata_ind -12 , missingdata_ind + 12),] %>% 
    group_by(month) %>% 
    summarise(across(everything(), sum)) %>% 
    mutate(year=year/2,month=month/2) %>% ungroup() %>% as.data.frame()
  xco2_df_diff <- as.data.frame(diff(as.matrix(xco2_df%>% subset(select=-c(year,month))),lag = 12))
  xco2_df_diff <- cbind(xco2_df[as.numeric(rownames(xco2_df_diff)),c("year","month")],xco2_df_diff)

  
  missingdiff_ind <- which(is.na(xco2_df_diff$mean)) %>% matrix(nrow=2,byrow = TRUE)
  library(zoo)
  missingdata_inter1 <- xco2_df_diff[c(missingdiff_ind[1,] - 1, missingdiff_ind[1,] ,missingdiff_ind[1,] + 1),] %>% 
    unique() 
  missingdata_inter1_NA_ind <- which(is.na(missingdata_inter1$mean))
  missingdata_inter1$mean <- zoo(missingdata_inter1$mean,as.numeric(rownames(missingdata_inter1))) %>% na.approx()
  missingdata_inter1$median <- zoo(missingdata_inter1$median,as.numeric(rownames(missingdata_inter1))) %>% na.approx()
  missingdata_inter1 <- missingdata_inter1[missingdata_inter1_NA_ind,]
  
  missingdata_inter2 <- xco2_df_diff[c(missingdiff_ind[2,] - 1, missingdiff_ind[2,], missingdiff_ind[2,] + 1),] %>% 
    unique()
  missingdata_inter2_NA_ind <- which(is.na(missingdata_inter2$mean))
  missingdata_inter2$mean <- zoo(missingdata_inter2$mean,as.numeric(rownames(missingdata_inter2))) %>% na.approx()
  missingdata_inter2$median <- zoo(missingdata_inter2$median,as.numeric(rownames(missingdata_inter2))) %>% na.approx()
  missingdata_inter2 <- missingdata_inter2[missingdata_inter2_NA_ind,]
  
  joinval <- missingdata_inter3 %>% subset(select=-c(mean,median))
  missingdata_inter1_mat <- missingdata_inter1 %>% subset(select=c(mean,median)) %>% as.matrix()
  missingdata_inter2_mat <- missingdata_inter2 %>% subset(select=c(mean,median)) %>% as.matrix()
  missingdata_inter3_mat <- missingdata_inter3 %>% subset(select=c(mean,median)) %>% as.matrix()
  impute_data <- (missingdata_inter1_mat - missingdata_inter2_mat + missingdata_inter3_mat)/2
  xco2_df[missingdata_ind,c("mean","median")] <- impute_data
  return(xco2_df)
}

# calculate the moving window trend bt
xco2_moving_time_window_trend<-function(xco2_df, n){
  xco2_df2 <- xco2_df[order(nrow(xco2_df):1),]
  xco2_rollsum <- xco2_df2 %>% as.matrix() %>% RcppRoll::roll_mean(n,fill= NA, align = "left") %>% as.data.frame() %>% dplyr::select(mean,median)
  
  xco2_rollsum <-xco2_rollsum[complete.cases(xco2_rollsum),]
  get_date <- xco2_df2[1:nrow(xco2_rollsum),c("year","month")]
  
  
  xco2_rollsum <- cbind(get_date, xco2_rollsum)
  xco2_rollsum <- xco2_rollsum[order(nrow(xco2_rollsum):1),]
  row.names(xco2_rollsum) <- NULL
  return(xco2_rollsum)
}

# diff function to caclulate mt
Calc_diff <- function(xco2_df, n, f = 1){
  if ("halfmonth" %in% names(xco2_df)){xco2_df <- xco2_df %>% subset(select=-c(halfmonth))}
  xco2_diff <- diff(as.matrix(xco2_df),lag = n) %>% as.data.frame() %>% dplyr::select(mean,median)
  xco2_diff <- xco2_diff[complete.cases(xco2_diff),] 
  xco2_diff <- xco2_diff/f
  xco2_diff_date <- tail(xco2_df,dim(xco2_diff)[1]) %>% dplyr::select(year,month)
  xco2_diff <- cbind(xco2_diff_date,xco2_diff)
  row.names(xco2_diff) <- NULL
  return(xco2_diff)
}

# calculate at
xco2_calc_anomaly<-function(rg_ind, xco2_df_list, trend_list, rg_namelist){
  rg_name <- rg_namelist[rg_ind]
  xco2 <- xco2_df_list[[rg_name]]
  trend <- trend_list[[rg_name]]
  anomaly <- full_join(xco2,trend,by=c("year","month")) %>%
    mutate(mean.diff = mean.x - mean.y, median.diff = median.x - median.y)
  anomaly <- anomaly[complete.cases(anomaly),] %>% dplyr::select(year, month, mean.diff, median.diff)
  names(anomaly) <- c("year", "month", "mean", "median")
  return(anomaly)
}

# draw contrast filter
draw_plot_month_compare <- function(draw_df,
                                    title_str = "Anomaly time series {*a*<sub>t</sub>} - {*a*<sub>2018</sub>}",
                                    pt_shape=1,
                                    pt_size=3,
                                    spanyear=c(2015,2022,1),
                                    yscale=c(-0.4,1.2,0.4),
                                    yhline=FASLE,
                                    yhline_type="dotted",
                                    yhline_intercept=0.0,
                                    yhline_size=1.0,
                                    yhline_color="grey70",
                                    tick_margin=c(0,0,0,0)){
  plotlist_at_comparison<-list()
  
  for (mon_i in c(12,1:11)){
    compare_df <- draw_df %>% dplyr::filter(month == mon_i) %>% dplyr::select(year,mean)
    compare_median_val <- median(compare_df$mean,na.rm = TRUE)
    compare_df$mean <- compare_df$mean - compare_median_val
    names(compare_df) <- c("year","diff")
    print(summary(compare_df$diff))
    comparison_plot <- ggplot(compare_df,aes(x=year,y=diff)) +
      theme_timeserials() +
      geom_line() +
      geom_point(shape=pt_shape,size=pt_size) +
      scale_x_continuous(breaks=seq(spanyear[1],spanyear[2],by=spanyear[3]),limits=c(spanyear[1],spanyear[2]),expand = c(0.01,0.01),labels = seq(spanyear[1],spanyear[2],by=spanyear[3])%%100) +
      scale_y_continuous(breaks=seq(yscale[1],yscale[2],yscale[3]),limits=c(yscale[1]-tick_margin[3], yscale[2]+tick_margin[4])) +
      labs(x="",y="",title=month.name[mon_i]) + 
      theme(plot.title = element_markdown(size = 16, lineheight = 3.0, margin=margin(0,0,0,0),vjust = -10),
            plot.margin = unit(c(1,8,1,6), "mm"),
            axis.title.y = element_markdown(size = 16, lineheight = 2.0))
    if(yhline) {
      comparison_plot <- comparison_plot + geom_hline(yintercept=yhline_intercept, linetype=yhline_type,size=yhline_size, color=yhline_color)
    }
    plotlist_at_comparison[[length(plotlist_at_comparison) + 1]] <- comparison_plot
  }
  
  tmpplot<-plot_grid(plotlist_at_comparison[[1]],plotlist_at_comparison[[2]],plotlist_at_comparison[[3]],
                     plotlist_at_comparison[[4]],plotlist_at_comparison[[5]],plotlist_at_comparison[[6]],
                     plotlist_at_comparison[[7]],plotlist_at_comparison[[8]],plotlist_at_comparison[[9]],
                     plotlist_at_comparison[[10]],plotlist_at_comparison[[11]],plotlist_at_comparison[[12]],
                     ncol = 3)
  tmpplot_title <- ggplot() + 
    labs(title = title_str) +
    theme_minimal() +
    theme(plot.title = element_markdown(size = 20,hjust = 0.5),
          plot.margin=unit(c(2,0,0,1),"mm"))
  tmpplot <- plot_grid(tmpplot_title, tmpplot, ncol = 1, rel_heights = c(0.05, 1))
  
  return(tmpplot)
}

