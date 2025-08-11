suppressMessages({
  library(sp)
  library(spacetime)
  library(ggplot2)
  library(dplyr)
  library(DBI)
  library(RSQLite)
  library(lubridate)
  library(mapproj)
  library(FRK)
  library(raster)
  library(rasterVis)
  library(cowplot)
  library(scales)
  library(tidyr)
  #library(conflicted)
  library(stringr)
  library(parallel)
  library(gridExtra)
  library(tidyverse)
  library(ozmaps)
  library(sf)
  library(patchwork)
  library(RcppRoll)
  library(ggtext)
  library(fields) 
})

# load the functions and constant variables
source('Earth-atmospheric-co2-monitoring-function.R')
load("test.RData")
# set the variables for plots
# time duration of the dataset, monthly label
x_label_all <- c("Feb\n2022","Jan\n2022",
                 "Dec\n2021","Nov\n2021","Oct\n2021","Sep\n2021","Aug\n2021","Jul\n2021","Jun\n2021","May\n2021","Apr\n2021","Mar\n2021","Feb\n2021","Jan\n2021",
                 "Dec\n2020","Nov\n2020","Oct\n2020","Sep\n2020","Aug\n2020","Jul\n2020","Jun\n2020","May\n2020","Apr\n2020","Mar\n2020","Feb\n2020","Jan\n2020",
                 "Dec\n2019","Nov\n2019","Oct\n2019","Sep\n2019","Aug\n2019","Jul\n2019","Jun\n2019","May\n2019","Apr\n2019","Mar\n2019","Feb\n2019","Jan\n2019",
                 "Dec\n2018","Nov\n2018","Oct\n2018","Sep\n2018","Aug\n2018","Jul\n2018","Jun\n2018","May\n2018","Apr\n2018","Mar\n2018","Feb\n2018","Jan\n2018",
                 "Dec\n2017","Nov\n2017","Oct\n2017","Sep\n2017","Aug\n2017","Jul\n2017","Jun\n2017","May\n2017","Apr\n2017","Mar\n2017","Feb\n2017","Jan\n2017",
                 "Dec\n2016","Nov\n2016","Oct\n2016","Sep\n2016","Aug\n2016","Jul\n2016","Jun\n2016","May\n2016","Apr\n2016","Mar\n2016","Feb\n2016","Jan\n2016",
                 "Dec\n2015","Nov\n2015","Oct\n2015","Sep\n2015","Aug\n2015","Jul\n2015","Jun\n2015","May\n2015","Apr\n2015","Mar\n2015","Feb\n2015","Jan\n2015",
                 "Dec\n2014","Nov\n2014","Oct\n2014","Sep\n2014")

# grey rect on the NX for drawing plots
grayrect<-data.frame(x1=c(-180,-180), x2=c(180,180), y1=c(0,87+1/2), y2=c(23-1/2,90))

# define the plot theme
# time seriel plots
theme_timeserials <- function(base_size=16, base_family="Helvetica") {
  suppressMessages({
    library(grid)
    library(ggthemes)
  })
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.5), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.35)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = rel(1.2)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_blank(),
            legend.text = element_text(size = rel(1.2)),
            legend.margin=margin(-15, 0, 0, 0),
            plot.margin=unit(c(5,15,5,15),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

# theme for publication
theme_Publication <- function(base_size=14, base_family="Helvetica") {
  suppressMessages({
    library(grid)
    library(ggthemes)
  })
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2, size = rel(1.05)),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = rel(1.0)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(),
            legend.text = element_text(size = rel(1.0)),
            legend.margin=margin(-15, 0, 0, 0),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

# range of fill in the spatial maps
fill_lim<-c(-6,6)
fill_lim_breaks<-seq(-6,6,3)
mt_fill_lim<-c(-0.6,0.6)
mt_fill_lim_breaks<-seq(-6,6,3)/10
lat_upperlim <- 90
#---------------------------------------------------------------------------

# Due to the size of the whole dataset, calculated variables are provided to assist
# the reproduction of the results.
# To run the code snippets of functions related to spatial results. The first step is to download the 
# TCCON-calibrate version 10r (v10r) OCO-2 retrievals and extract the xco2 value with xco2_quality_flag = 0
# For time series results, a variable "sub_region_monthly_summary_impute" of monthly summary (i.e. average)
# is provided, the other variables related to time series results are calculate from this variable. 
# For spatial results, the code snippets are explained. Tu run the codes, the dataset from OCO-2 butof generating the results.  
# Functions used to calculate results are available at "Earth-atmospheric-co2-monitoring-function.R" 
#--------------------------------------------------------------------------------------

##########
## spatial results
##########
## The whole dataset of the TCCON-calibrate version 10r (v10r) OCO-2 retrievals is hard to share on Github.
## Please download the dataset from the NASA Goddard Earth Sciences (GES) Data and Information Services Center (DISC).
## There are many methods to load the dataset,  choose one method that you are familiar with,
## make sure the format is the same as shown in below:
#   datetime      lon      lat     xco2   year month
# 1 2014-09-06 160.5118 54.46574 394.4846 2014     9
# 2 2014-09-06 160.5019 54.48532 392.9321 2014     9
# 3 2014-09-06 160.5380 54.48214 394.0697 2014     9
# 4 2014-09-06 160.5100 54.50333 394.2112 2014     9
# The following code snippets assume the dataset is loaded to variable "xco2_dataset" and show the processing in the NX
if (exists("xco2_dataset")){
  get_xco2_subset<-xco2_dataset %>% dplyr::filter(lat>=21)
  num_cores<-detectCores(logical = FALSE)/4
  # Get the interpolated xco2 using localloess spline, anomaly (i.e.,at), bt and mt (using the interpolated xco2)
  xco2_proc_anomaly_smooth<-data.frame()
  
  for (blk_lat_i in unique(centrepoints$lat)){
    blk_centers<-centrepoints %>% dplyr::filter(lat == blk_lat_i)  
    print(paste0("  processing lat band at ",toString(blk_lat_i)," , cell number = ",dim(blk_centers)[1]))
    blk_centers_list<-unname(as.list(as.data.frame(t(blk_centers))))
    xco2_tmp_proc_lat_band <- get_xco2_subset %>%
      dplyr::filter(lat >= (blk_lat_i - (aggregate_var$window_width_lat + 1)) & lat <= (blk_lat_i + (aggregate_var$window_width_lat + 1)))
    tmp_res_list<-mclapply(blk_centers_list,
                           cell_xco2_impute_summary,
                           data=xco2_tmp_proc_lat_band,
                           lat_width=aggregate_var$window_width_lat,
                           lon_width=aggregate_var$window_width_lon,
                           loess_period = 6, 
                           loess_degree = 1,
                           thres = 3,
                           mc.cores = num_cores)
    for (list_i in 1:length(tmp_res_list)){
      xco2_proc_anomaly_smooth<-rbind(xco2_proc_anomaly_smooth,tmp_res_list[[list_i]])
    }
  }
  
  # convert the datetime from unix timestamp to Date format, check on datetime, replace NaN to NA
  xco2_proc_anomaly_smooth <- xco2_proc_anomaly_smooth %>% 
    dplyr::filter(!is.na(datetime)) %>% 
    mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
    mutate(datetime=as_datetime(datetime)) 
  
  # Get the median map of at and mt
  at_median <- xco2_proc_anomaly_smooth %>% 
    dplyr::select(lon,lat,datetime,anomaly) %>% 
    mutate(month=month(datetime)) %>%
    group_by(lon,lat,month) %>% summarise_at("anomaly",median,na.rm=TRUE) %>%
    ungroup()
  names(at_median) <- c("lon","lat","month","anomalymedian")
  
  mt_median <- xco2_proc_anomaly_smooth %>% 
    dplyr::select(lon,lat,datetime,mt) %>% 
    mutate(month=month(datetime)) %>%
    group_by(lon,lat,month) %>% summarise_at("mt",median,na.rm=TRUE) %>%
    ungroup()
  names(mt_median) <- c("lon","lat","month","mtmedian")
  
  # calculate the anomaly contrast and mt contrast
  at_mt_median <- full_join(at_median,mt_median,by=c("lon","lat","month"))
  xco2_proc_anomaly_smooth <- xco2_proc_anomaly_smooth %>% mutate(year=year(datetime),month=month(datetime),day=day(datetime))
  xco2_proc_anomaly_smooth <- full_join(xco2_proc_anomaly_smooth,at_mt_median,by=c("lon","lat","month"))
  
  # calculate contrast
  xco2_proc_anomaly_smooth <- xco2_proc_anomaly_smooth %>% mutate(anomalycontrast=anomaly-anomalymedian,mtcontrast=mt-mtmedian)
  
  # get at median map for August
  print(paste0("Get the median map for ",month.name[8]))
  at_median_monthly <- at_median %>% dplyr::filter(month==8) %>% dplyr::select(lon,lat,month,anomalymedian)
  names(at_median_monthly) <- c("lon","lat","month","value") 
  
  # get mt contrast for 2021 July
  draw_data_mtcontrast <- xco2_proc_anomaly_smooth %>% 
    dplyr::filter((year == 2021) & (month == 7) & (lat <= lat_upperlim)) %>%
    dplyr::select(lon,lat,mtcontrast)
  names(draw_data_mtcontrast) <- c("lon","lat","value")
  
  # NX spatial CO2 monthly trend average contrast 
  mtcontrast_ave <- xco2_proc_anomaly_smooth %>% 
    dplyr::select(lon,lat,year,month,mtcontrast) %>% 
    group_by(year,month) %>% 
    summarise(avg=mean(mtcontrast,na.rm=TRUE),.groups = "drop") %>%
    mutate(time_label = paste0(month.abb[month],"\n",as.character(year)), ind = row_number()) %>%
    mutate(yearly_avg=roll_sumr(avg,n=12)) %>% 
    dplyr::filter(ind >= 62) # start from Oct.1 2019
  
}
# an example to draw at median map for August
at_median_august <- spatial_map(draw_df=at_median_monthly %>% dplyr::rename(x=lon,y=lat),
            plot_title=paste0("N Extra-tropics ",
                              "median spatial CO<sub>2</sub> anomaly",
                              " for ",
                              month.name[8]),
            fill_lim = fill_lim,
            legend_breaks = fill_lim_breaks)

# an example to draw mt contrast for 2021 July
mtcontrast_2021_July <- spatial_map(draw_df=draw_data_mtcontrast %>% dplyr::rename(x=lon,y=lat),
            plot_title=paste0("N Extra-tropics ",
                              "spatial CO<sub>2</sub> ",
                              "monthly",
                              " trend contrast",
                              ": ",
                              2021,
                              " ",
                              month.abb[7]),
            fill_lim = mt_fill_lim,
            legend_breaks = mt_fill_lim_breaks)

# NX spatial CO2 monthly trend average contrast map
mtcontrast_time_series <- (mtcontrast_ave %>% dplyr::filter(year > 2019) %>%
  ggplot(aes(x=ind,y=avg)) + 
  theme_timeserials() + 
  geom_line() +
  geom_point(shape=1,size=3) +
  scale_x_continuous(breaks=mtcontrast_ave$ind[seq(4,length(mtcontrast_ave$ind),2)],labels=mtcontrast_ave$time_label[seq(4,length(mtcontrast_ave$ind),2)],expand = c(0.01,0.01)) +
  scale_y_continuous(breaks=c(seq(-0.03,0.03,0.01)),limits = c(-0.032,0.032)) +
  labs(x="Month",
       y=paste0("CO<sub>2</sub> ","monthly"," trend contrast (ppm/mo)"),
       title=paste("N Extra-tropics spatial CO<sub>2</sub> ","monthly"," trend average contrast map")) +
  theme(legend.position = "none",
        plot.title = element_markdown(size = 20),
        axis.title.x = element_markdown(size = 18),
        axis.title.y = element_markdown(size = 18)))
##########
## time series results
##########
if (exists("xco2_dataset")){
  # if dataset is available, compute the variables for time series results
  sub_region_monthly_summary<-mclapply(sub_region, xco2_temporal_summary, xco2_dataset, mc.cores = length(sub_region))
  sub_region_monthly_summary_impute <- lapply(sub_region_monthly_summary, xco2_missingdata_impute)
}

## An example to plot the global CO2 using variable sub_region_monthly_summary_impute
x_tick_max<-length(x_label_all)
# as the input label is in reverse order, reverse it
x_label <- rev(x_label_all)
region_notation = "globalearth"
method_str = "mean"
title_str = "Global CO<sub>2</sub> from OCO-2 satellite"
Ymin = 395
Ymax = 416.3
Yinterval = 5
label_interval = 4
ylabel="CO<sub>2</sub> (ppm)"

drawdata <- sub_region_monthly_summary_impute[[region_notation]] 
getnames <- names(drawdata)
drawdata <- drawdata[,method_str]

want_label<-seq(1,length(x_label),label_interval)
x_label[-want_label]=""
x_label_ind <- seq(1,length(x_label),1)
breaks_ind <- 1:x_tick_max
breaks_label <- rep("",x_tick_max)
breaks_label[x_label_ind] <- x_label

if(length(drawdata) < x_tick_max){
  breaks_ind <- (1 + x_tick_max - length(drawdata)):x_tick_max
  breaks_label <- breaks_label[breaks_ind]
  ave_res<-data.frame(ind = (1 + x_tick_max - length(drawdata)):x_tick_max,
                      data = drawdata %>% tail(x_tick_max))
  
}else{
  ave_res<-data.frame(ind = 1:x_tick_max,
                      data = drawdata %>% tail(x_tick_max))
}

Xmax <- x_tick_max
suppressMessages({  library(extrafont)})
extrafont::loadfonts()

# get the color by scales::show_col(scales::hue_pal()(2))
Global_CO2 <- ggplot(ave_res,aes(x=ind,)) +
  theme_timeserials() +
  geom_line(aes(y = data)) +
  geom_point(aes(y = data), shape=1,size=3) +
  scale_x_continuous(breaks=breaks_ind,labels=breaks_label,limits=c(1,length(x_label)),expand = c(0.01,0.01)) +
  scale_y_continuous(breaks=seq(Ymin,Ymax,Yinterval),
                     limits=c(min(c(Ymin,ave_res$south,ave_res$north)),
                              max(c(Ymax,ave_res$south,ave_res$north)))) +
  labs(x="Month",y="XCO2 (PPM)",title=title_str) + theme(legend.position = "none",
                                                         plot.title = element_markdown(size = 20),
                                                         axis.title.x = element_markdown(size = 18),
                                                         axis.title.y = element_markdown(size = 18))
##

# Calculate bt
sub_region_bt <- lapply(sub_region_monthly_summary_impute, xco2_moving_time_window_trend,lag_num)
# Calculate mt 
sub_region_mt <- lapply(sub_region_monthly_summary_impute, Calc_diff, 12, 12)
# Calculate yt for NX, SX, Global and tropics
sub_region_yt <- sub_region_mt[c("globalearth","northern_extratropics","southern_extratropics","tropics")] 
sub_region_yt$globalearth <- sub_region_mt$globalearth %>% mutate(movingmean=roll_sumr(mean,n = 12))
sub_region_yt$northern_extratropics <- sub_region_mt$northern_extratropics %>% mutate(movingmean=roll_sumr(mean,n = 12))
sub_region_yt$southern_extratropics <- sub_region_mt$southern_extratropics %>% mutate(movingmean=roll_sumr(mean,n = 12))
sub_region_yt$tropics <- sub_region_mt$tropics %>% mutate(movingmean=roll_sumr(mean,n = 12))
# Calculate at
sub_region_ind <- as.data.frame(t(seq(1:dim(sub_region)[2])))
names(sub_region_ind) <- names(sub_region)
rg_name <- names(sub_region)
sub_region_at <- lapply(sub_region_ind, xco2_calc_anomaly, sub_region_monthly_summary_impute, sub_region_bt, rg_name)

# To draw the contrast filter in the Figture 8 of paper, please call draw_plot_month_compare function
# Here is an example to draw monthly trend contrast filter of NX
mt_NX <- sub_region_mt$northern_extratropics %>% dplyr::select(year,month,mean)
mt_tilde <- draw_plot_month_compare(mt_NX,
                                    yscale=c(-0.1,0.1,0.05),
                                    spanyear=c(2015,2022,1),
                                    title_str = "N Extra-tropics: Change in CO<sub>2</sub> monthly trend (by month, in ppm)",
                                    yhline=TRUE,
                                    yhline_type="dotted",
                                    yhline_intercept=0.0,
                                    yhline_size=1,
                                    yhline_color="grey70")

# permutation monte carlo is used to infer the extreme 
# The result is rank, p_value and color_tag, see paper for color definition
# Then, you can overlay color on the time series plot.
# This function can also be used on spatial map to get the permuation monte carlo result
#   by feeding spatial values to the function.
set.seed(123)
df_permute_res <- calculate_permutation_ranks(mt_NX$mean)
