##
## Functions to calculate growth rates, plot the log(OD) vs. time, and compare doubling times
##

#####
#####
#Define function to calculate growth rate and max yield 
#####
#####
#This function will calculate the maximum growth rate and plot the fits
#Required data format: A data frame with the first column as the time and the remaining columns are for each sample's time course.
#The rows are the measurement at each time point.

calc_growthrate = function(growthcurve_dataframe, PlateSetUp, range_to_consider_max, spline_size) { #Manually inspect growth curves and choose the "range_to_consider_max" as the period of exponential growth. Choose a spline_size - the size of the window over which a linear regression is fitted to calculate the growth rate.
  
  #Create directory
  dir.create("growth_rate_fit_plots/")
  
  #Define variables
  time = growthcurve_dataframe[,1] #time column
  `%notin%` <- Negate(`%in%`)
  well_index = names(growthcurve_dataframe)[colnames(growthcurve_dataframe) %notin% metadata[metadata$Name == "Blank","Well"]]
  growthrates_dataframe = data.frame(matrix(ncol=4,nrow=length(well_index) - 1)) #create data frame with 6 columns and rows for each sample
  colnames(growthrates_dataframe) = c('temp_rsquared', 'temp_slope', 'started', 'ended') # We'll store the R-squared, the growth rate, and the start and end of the spline/interval
  rownames(growthrates_dataframe) = well_index[-1] #each sample will be a row
  len_time = length(time) - spline_size

  
  #Plot fits and save files
  for ( i in 2:length(well_index) ) {
    tryCatch({
      #Create growthcurve_dataframe where fits will be stored for each spline
      all_fits_per_spline = data.frame(matrix(ncol=4, nrow= len_time)) #creates temporary dataframe that gets overwritten each loop. This has the 4 columns (as named below) and each row is for the different splines for which the regression is calculated over the time interval
      colnames(all_fits_per_spline) = c('temp_rsquared', 'temp_slope', 'started', 'ended')
      
      #Load data and take the rolling mean every seven measurements
      log_data = log(matrix(as.numeric(unlist(growthcurve_dataframe[well_index][i])),nrow=length(time))) #log transform
      mean_logdata = zoo::rollmean(log_data, k=7, fill = NA, na.rm = TRUE) # rolling mean over a window of 7
      
      #Populate all_fits_per_spline with the fits for the rolling splines
      for (started in 1:length(time)) {
        if (started <= length(time) - spline_size) {
          ended = started + spline_size
          spline_window = started:ended
          fit = lm(mean_logdata[spline_window] ~ as.matrix(growthcurve_dataframe[colnames(growthcurve_dataframe) %notin% metadata[metadata$Name == "Blank","Well"]][spline_window,1])) # linear regression then store all_fits_per_spline into separate columns
          temp_rsquared = summary(fit)$r.squared
          temp_slope = summary(fit)$coefficients[2,1]
          temp_intercept = summary(fit)$coefficients[1,1]
        } 
        all_fits_per_spline[started,] = data.frame(temp_rsquared, temp_slope, started, ended, stringsAsFactors=FALSE)
      }
      
      #Find the max growth rate that meets certain criteria
      #Save the time and window at which the max growth rate occurs for plotting
      max_time = which(all_fits_per_spline[,2] == all_fits_per_spline %>% 
                          filter(temp_rsquared >= 0.99) %>% # Consider only fits that had an R2 > 0.99
                          filter(ended <= max(range_to_consider_max)) %>% #Report the max growth rate over the given interval - "range_to_consider_max"
                          filter(started >= min(range_to_consider_max)) %>%
                          select(temp_slope) %>%
                          max())# record the time for plotting
      
      max_spline_window = seq(max_time, max_time + spline_size, 1) # record the window for plotting

      # Add the max growth rate to a single data frame for all samples. 
      growthrates_dataframe[i-1,] = all_fits_per_spline[max_time, ] 

      #Plot results
      pdf(file = paste("growth_rate_fit_plots/",
                       PlateSetUp[PlateSetUp$Well == row.names(growthrates_dataframe)[i-1],7],
                       "_",
                       PlateSetUp[PlateSetUp$Well == row.names(growthrates_dataframe)[i-1],4],
                       ".pdf",
                       sep=""), width=5, height=15)
      par(mfrow=c(4,1))
      plot(mean_logdata, ylab="ln(OD)", type="p",
           main = paste("Strain:", PlateSetUp[PlateSetUp$Well == row.names(growthrates_dataframe)[i-1],7],
                              "Lineage:", PlateSetUp[PlateSetUp$Well == row.names(growthrates_dataframe)[i-1],5],
                              "SLST:", PlateSetUp[PlateSetUp$Well == row.names(growthrates_dataframe)[i-1],6], sep="  "),
           xlim=c(0,150))
      points(all_fits_per_spline[max_time + spline_size/2,3], mean_logdata[max_time + spline_size/2], col="red", pch=19)
      abline(lm(mean_logdata[max_spline_window] ~ max_spline_window), col= "red")
      plot(all_fits_per_spline[,2], ylab="growth rate", col=ifelse(all_fits_per_spline[,2] == all_fits_per_spline[max_time,2] , "red", "black"), xlim=c(0,150))
      plot(all_fits_per_spline[,1], ylab="r-squared", xlim=c(0,150), ylim=c(0.99,1))
      points(all_fits_per_spline[max_time,3], all_fits_per_spline[max_time,1], col="red")
      plot(diff(zoo::rollmean(all_fits_per_spline[,2], k=3, fill = NA, na.rm = TRUE)), ylab="second derivative", xlim=c(0,150), ylim=c(-0.01,0.01))
      abline(h=0, col="blue")
      points(all_fits_per_spline[max_time,3], diff(zoo::rollmean(all_fits_per_spline[,2], k=3, fill = NA, na.rm = TRUE))[max_time-1], col="red")
      dev.off()
      #The fit for each well is saved as a plot.
      #Inspect the plots and adjust the range_to_consider_max and spline_size as needed.

    }, error=function(e){cat("ERROR :",conditionMessage(e), well_index, "\n")}) #Loop keeps running and prints any errors that occur with each iteration rather than breaking
  }

  while (!is.null(dev.list()))  dev.off() # Prevent files from being overwritten
  
  #Convert to data table
  growthrates_dataframe = setDT(growthrates_dataframe, keep.rownames = TRUE)[]
  names(growthrates_dataframe)[1] = "Well" # Rename the first column 
  
  #Add doubling time ln(2)/r
  growthrates_dataframe$doubling = log(2)/growthrates_dataframe$temp_slope
  
  #Save output as .csv file
  write.csv(growthrates_dataframe, "growthrates_dataframe.csv", row.names = FALSE)
  
  #Store create data frame in global environment
  return(growthrates_dataframe)
  
}

#####
#####

#####
#####
#Plot log10(OD) vs. time for each sample and its replicates
#####
#####
#This function requires the growthcurve_dataframe and a metadata file.
#The metadata file needs to have the format as the "~/plate_set_up.csv" file. 


##Individual OD plots by strain

plot_odbystrain <- function(gcdf, PlateSetUp) { #two inputs required
 
  #Create directory
  dir.create("plot_odbystrain/")
  
  #Create list of strains
  list_strains <- as.character(unique(unlist(PlateSetUp$Name))) #groups wells by the "Names column and plots together based on this
  
  for (i in 1:length(list_strains)) {
    gcdf_i <- gcdf[,names(gcdf) %in% c("time", as.vector(unlist(PlateSetUp[PlateSetUp$Name == list_strains[i], "Well"])))] %>% 
      mutate(time = factor(time, levels = unique(time))) %>%
      gather(sample, value, -time)
    
    ggplot(NULL, aes(time,value)) + 
      geom_line(data = gcdf_i, aes(color = sample, group = sample), alpha=0.75, size=1) +
      geom_smooth(data = gcdf_i, span=0.25, se=FALSE, aes(color = sample, group = sample), alpha=0.75, size=1) +
      theme_minimal() + xlab("Time (hrs)") + ylab("log (OD)") +
      scale_y_continuous(trans='log10', limits=c(0.05, 1.5)) +
      scale_color_manual(values = colorBlind9) +
      scale_x_discrete(breaks = list(gcdf[which.min(abs(gcdf[,1]-0)),1], gcdf[which.min(abs(gcdf[,1]-24)),1],
                                     gcdf[which.min(abs(gcdf[,1]-48)),1], gcdf[which.min(abs(gcdf[,1]-72)),1],
                                     gcdf[which.min(abs(gcdf[,1]-96)),1]),
                       labels = c("0", "24", "48", "72", "96"))  + #manually sets the time on the x-axis with custom label
      theme(axis.text.x = element_text(), axis.text.y = element_text(), axis.ticks = element_line())+
      ggtitle(list_strains[i])+ 
      theme(legend.position = "bottom", legend.title = element_blank()) 
    ggsave(paste("plot_odbystrain/", "plotbystrain_", list_strains[i], ".pdf", sep=""), width = 20, height = 20, units = "cm") # saves plots
  }
}


##Individual OD plots by lineage

plot_odbylineage <- function(gcdf, PlateSetUp) {
  
  #Create directory
  dir.create("plot_odbylineage/")
  
  #Create list of lineages
  list_lineage <- as.character(unique(unlist(PlateSetUp$Lineage)))
  
  #Plot replicates as separate but still plot by lineage
  for (i in 1:length(list_lineage)) {
    gcdf_i <- gcdf[,names(gcdf) %in% c("time", as.vector(unlist(PlateSetUp[PlateSetUp$Lineage == list_lineage[i], "Well"])))] %>% 
      mutate(time = factor(time, levels = unique(time))) %>%
      gather(sample, value, -time)
    
    ggplot(NULL, aes(time,value)) + 
      geom_line(data = gcdf_i, aes(color = sample, group = sample), alpha=0.75, size=1) +
      geom_smooth(data = gcdf_i, span=0.25, se=FALSE, aes(color = sample, group = sample), alpha=0.75, size=1) +
      theme_minimal() + xlab("Time (hrs)") + ylab("log (OD)") +
      scale_color_manual(values = coul) +
      scale_y_continuous(trans='log10', limits=c(0.05, 1.5)) +
      scale_x_discrete(breaks = list(gcdf[which.min(abs(gcdf[,1]-0)),1], gcdf[which.min(abs(gcdf[,1]-24)),1],
                                     gcdf[which.min(abs(gcdf[,1]-48)),1], gcdf[which.min(abs(gcdf[,1]-72)),1],
                                     gcdf[which.min(abs(gcdf[,1]-96)),1]),
                       labels = c("0", "24", "48", "72", "96"))  +
      theme(axis.text.x = element_text(), axis.text.y = element_text(), axis.ticks = element_line())+
      ggtitle(list_lineage[i])+ 
      theme(legend.position = "bottom", legend.title = element_blank())
    ggsave(paste("plot_odbylineage/", "plotbylineage_", list_lineage[i], ".pdf", sep=""), width = 20, height = 20, units = "cm")
  }
}


##Individual OD plots by SLST

plot_odbyslst <- function(gcdf, PlateSetUp) {
  
  #Create directory
  dir.create("plot_odbyslst/")
  
  #Create list of SLSTs
  list_slst <- as.character(unique(unlist(PlateSetUp$SLST)))
  
  #Plot replicates as separate but still plot by slst
  for (i in 1:length(list_slst)) {
    df1 <- gcdf[,names(gcdf) %in% c("time", as.vector(unlist(PlateSetUp[PlateSetUp$SLST == list_slst[i], "Well"])))] %>% 
      mutate(time = factor(time, levels = unique(time))) %>%
      gather(sample, value, -time)   
 
    ggplot(NULL, aes(time,value)) + 
      geom_line(data = df1, aes(time, value, color = sample, group = sample), alpha=0.75, size=1) +
      geom_smooth(data = df1, span=0.25, se=FALSE,  aes(time, value, color = sample, group = sample), alpha=0.75, size=1) +
      theme_minimal() + xlab("Time (hrs)") + ylab("log (OD)") +
      # scale_color_manual(values = coul) +
      scale_y_continuous(trans='log10', limits=c(0.05, 1.5)) +
      scale_x_discrete(breaks = list(gcdf[which.min(abs(gcdf[,1]-0)),1], gcdf[which.min(abs(gcdf[,1]-24)),1],
                                     gcdf[which.min(abs(gcdf[,1]-48)),1], gcdf[which.min(abs(gcdf[,1]-72)),1],
                                     gcdf[which.min(abs(gcdf[,1]-96)),1]),
                       labels = c("0", "24", "48", "72", "96"))  +
      theme(axis.text.x = element_text(), axis.text.y = element_text(), axis.ticks = element_line())+
      ggtitle(list_slst[i])+ 
      theme(legend.position = "bottom", legend.title = element_blank())
    ggsave(paste("plot_odbyslst/", "plotbySLST_", list_slst[i], ".pdf", sep=""), width = 20, height = 20, units = "cm")
  }
}


#####
#####


