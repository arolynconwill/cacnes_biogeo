#####
#####
# Generate data frames and save plots
#####
#####

##
#Load packages
require('magrittr')
require(reshape2)
require(stringr)
require(plyr);  require(dbplyr); require(dplyr)
require(tidyr)
require(scales)
require(gridExtra)
require(splines)
require(data.table)
require(readxl) 
require(ggplot2)  
require(zoo)   
require(ggpubr)
require(RColorBrewer)
require(Hmisc)

## Load color palette
coul <- brewer.pal(12, "Set3")
coul <- colorRampPalette(coul)(20)

#Colorblind Friendly
colorBlind9  <- c("#D55E00", "#56B4E9", "#F0E442", "#000000", "#E69F00",  "#009E73",
                  "#0072B2",  "#CC79A7", "#999999")

##
#Source functions
source("Functions.R")

##
#Import the growth curve data
# Required data format: A data frame with the first column as the time and the remaining columns are for each sample's time course.
#The rows are the measurement at each time point.
growth_time_course <- read.csv(file.choose()) # "~/growth_curve_time_course.csv"
View(growth_time_course)

#Import metadata
metadata <- read.csv(file.choose()) # "~/plate_set_up.csv"
View(metadata)

##
#Plot the growth curve time series by the strain and save plots
plot_odbystrain(growth_time_course, metadata)

##
#Plot the growth curve time series by the lineage and save plots
plot_odbylineage(growth_time_course, metadata)

##
#Plot the growth curve time series by the slst and save plots
plot_odbyslst(growth_time_course, metadata)

##
#Calculate growth rates
#Provide the name of the data frame. 
#Enter the range to consider the max (here from the 6th to the 100th time measurement is used). Note: you will need to exam the growth curves and decided on the range. 
#Provide the spline size - the window over which to fit a linear regression (here the size is 15 time measurements).
growthrates_dataframe = calc_growthrate(growth_time_course, metadata, 6:100, 15) 

#Add metadata
growthrates_dataframe = left_join(growthrates_dataframe, metadata)

##
#Compare the slst types pairwise for the growth rate
compare_means(temp_slope ~ SLST, 
              growthrates_dataframe %>% filter(temp_slope != "NA"), 
              ref.group = ".all.", method = "wilcox.test")

##
#Compare doubling times and plot all isolates from each subject together
ggplot(growthrates_dataframe %>% mutate(Subject = ifelse(Subject == "Subject A", "Subject 1", "Subject 2")), 
       aes(x = factor(Name, level = as.character(unique(growthrates_dataframe[with(growthrates_dataframe, order(SLST, temp_slope)),]$Name))), y = doubling, color=SLST)) +
  ylim(5,18) +
  ggtitle("Growth Rate Comparison by Isolate") +
  xlab("Isolate") + 
  ylab("Doubling time (hrs)") +
  theme_classic() +
  theme(legend.position = "bottom", legend.text = element_text(size=14), legend.title = element_text(size=14)) + 
  theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14, angle = 90),
        axis.title.y  = element_text(size=16),
        axis.title.x  = element_text(size=16),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA)) + 
  stat_summary(fun=mean, geom="point", size=3.75, color="black")  +
  stat_summary(fun=mean, geom="point", size=3, aes(fill="legend"))  +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", color="black", alpha=0.5) + labs(col="Lineage") +
  scale_color_manual(values = c("red2", "darkorange1", "gold1", "lightblue1", "steelblue1")) +
  facet_wrap(~Subject, ncol =1, scales = "free") +
  stat_compare_means(method = "anova", label.y = 16.5, label.x = 1) 

ggsave("GrowthRateANOVAIsolate.pdf", width = 25, height = 30, units = "cm")

