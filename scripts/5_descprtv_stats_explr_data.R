################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############    5. Descriptive statistics and data exploration    #############
#############                                                      #############
#############                  By: Zach Laubach                    #############
#############                created: 21 Aug 2024                  #############
#############             last updated: 25 Feb 2024                #############
################################################################################


  ### PURPOSE: Calculate descriptive statistics and explore data 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Descriptive stats 

  


###############################################################################
##############             1.  Configure work space              ##############
###############################################################################

  ### 1.1 Global options
    ## a) clear global environment
      rm(list = ls())
    
    ## b) prevent R from automatically reading character strings as factors
      options(stringsAsFactors = FALSE)
    
  
  ### 1.2 Install and load CRAN packages   
    ## a) Data Manipulation and Descriptive Stats Packages
    # load tidyverse packages
      library ('tidyverse')
      
    # load igraph package
      library ('igraph')  
    
    # load lubridate package
      library ('lubridate')
      
    # load gridExtra packages
      library ('gridExtra')  
      
    # load here package
      library ('here')
    
    
    
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 4.4.2 (2024-10-31)
    # Platform: x86_64-apple-darwin20
    # Running under: macOS Sequoia 15.1.
    
  
  ### 1.4 Set working directory 
    setwd(here())
  
  
  ### 1.5 Source functions
    ## a) Source scripts path
      source_path <- paste("~/WD/Git/source_code/")
    
    ## b) all_char_to_lower function
      source(file = paste0(source_path, "all_char_to_lower.R"))
    
    ## c) format_var_names function
      source(file = paste0(source_path, "format_var_names.R"))
      
      

###############################################################################
##############                  2. Import data                   ##############
###############################################################################    

  ### 2.1 Import sample data files
    ## a) Load mediation data for CHR 2015
      load(here('data/5_6_chr15_mediation_data.RData'))
      
   

###############################################################################
##############               3. Descriptive stats                ##############
###############################################################################

  ### 3.1 Univariate descriptive stats for attribute data
    ## a) Total sample size
        n <- sum(!is.na(chr15_attrib_pre_df$Tag))
        
    ## b) Descriptive stats for sex  
      chr15_attrib_pre_df %>%
        group_by(Sex) %>%
        summarise (n.Sex = sum(!is.na(Sex)))%>%
        mutate(freq = n.Sex / sum(n.Sex)) %>%
        ungroup() 
        
    ## c) Descriptive stats for Age.category  
      chr15_attrib_pre_df %>%
        group_by(Sex, Age.category) %>%
        summarise(n.Age.category = sum(!is.na(Age.category)))%>%
        mutate(freq = n.Age.category / n) %>%
        ungroup()
      
    ## d) Descriptive stats morphology variables
      univar_morph <- chr15_attrib_pre_df %>%
        group_by(Sex) %>%
              # throat average brightness
        summarise (n.T_avg.bright = sum(!is.na(T_avg.bright)),
                   avg.T_avg.bright = round (mean(T_avg.bright, 
                                                  na.rm = T),2),
                   stdev.T_avg.bright = round (sd(T_avg.bright, 
                                                  na.rm = T), 2),
                   med.T_avg.bright = round(median(T_avg.bright,
                                                   na.rm = T), 2),
                   min.T_avg.bright = round(min(T_avg.bright,
                                                na.rm = T), 2),
                   max.T_avg.bright = round(max(T_avg.bright,
                                                na.rm = T), 2),
                   # breast average brightness pre-manip
                   n.R_avg.bright = sum(!is.na(R_avg.bright)),
                   avg.R_avg.bright = round (mean(R_avg.bright, 
                                                  na.rm = T),2),
                   stdev.R_avg.bright = round (sd(R_avg.bright, 
                                                  na.rm = T), 2),
                   med.R_avg.bright = round(median(R_avg.bright,
                                                   na.rm = T), 2),
                   min.R_avg.bright = round(min(R_avg.bright,
                                                na.rm = T), 2),
                   max.R_avg.bright = round(max(R_avg.bright,
                                                na.rm = T), 2),
                   # belly average brightness pre-manip
                   n.B_avg.bright = sum(!is.na(B_avg.bright)),
                   avg.B_avg.bright = round (mean(B_avg.bright, 
                                                  na.rm = T),2),
                   stdev.B_avg.bright = round (sd(B_avg.bright, 
                                                  na.rm = T), 2),
                   med.B_avg.bright = round(median(B_avg.bright,
                                                   na.rm = T), 2),
                   min.B_avg.bright = round(min(B_avg.bright,
                                                na.rm = T), 2),
                   max.B_avg.bright = round(max(B_avg.bright,
                                                na.rm = T), 2)
                   ) %>%
        ungroup()
      
    ## e) transpose the data frame for easier viewing
      univar_morph <- as.data.frame(t(univar_morph))
      univar_morph <- univar_morph %>%
        rename(c(female = V1,
                 male = V2))
      
    ## f) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_morph.pdf'), height = 13, width = 4)
      grid.table(univar_morph)
      dev.off()

    ## g) Descriptive stats morphology variables
      univar_R_brigh_post <- chr15_attrib_post_df %>%
        ungroup() %>%
        filter(Color.manipulation. == 'y') %>%  
      # breast average brightness post-manip
      summarise(n.R.bright.treat.and.orig = sum(!is.na
                                      (R.bright.treat.and.orig)),
      avg.R.bright.treat.and.orig = round(mean
                                           (R.bright.treat.and.orig, 
                                             na.rm = T),2),
      stdev.R.bright.treat.and.orig = round(sd
                                             (R.bright.treat.and.orig, 
                                               na.rm = T), 2),
      med.R.bright.treat.and.orig = round(median
                                          (R.bright.treat.and.orig,
                                            na.rm = T), 2),
      min.R.bright.treat.and.orig = round(min
                                          (R.bright.treat.and.orig,
                                            na.rm = T), 2),
      max.R.bright.treat.and.orig = round(max
                                          (R.bright.treat.and.orig,
                                            na.rm = T), 2))
      
  ### 3.2  Univariate descriptive stats for social intx. data
      ## a) Descriptive stats social intx. variables pre-manip
      univar_soc_intx_pre <- chr15_attrib_pre_df %>%
        group_by(Sex) %>%
        # number of unique social partners
        summarise (n.pre = sum(!is.na(degree.pre)),
                   med.degree.pre = round(median(degree.pre,
                                                   na.rm = T), 2),
                   qt25.degree.pre = round(quantile(degree.pre, 
                                                          0.25, na.rm = T), 2),
                   qt75.degree.pre = round(quantile(degree.pre, 
                                                          0.75, na.rm = T), 2),
                   avg.degree.pre = round(mean(degree.pre, 
                                             na.rm = T),2),
                   stdev.degree.pre = round(sd(degree.pre, 
                                             na.rm = T), 2),
                   min.degree.pre = round(min(degree.pre,
                                           na.rm = T), 2),
                   max.degree.pre = round(max(degree.pre,
                                           na.rm = T), 2), 
         # max number of social interactions with any one soc. partner
                  med.max.cnt.soc.intx.pre = round(median
                                              (max.cnt.soc.intx.pre,
                                                          na.rm = T), 2),
                  qt25.max.cnt.soc.intx.pre = round(quantile
                                              (max.cnt.soc.intx.pre, 
                                                          0.25, na.rm = T), 2),
                  qt75.max.cnt.soc.intx.pre = round(quantile
                                              (max.cnt.soc.intx.pre, 
                                                          0.75, na.rm = T), 2),
                  avg.max.cnt.soc.intx.pre = round (mean
                                              (max.cnt.soc.intx.pre, 
                                                         na.rm = T),2),
                  stdev.max.cnt.soc.intx.pre = round (sd
                                              (max.cnt.soc.intx.pre, 
                                                         na.rm = T), 2),
                  min.max.cnt.soc.intx.pre = round(min
                                              (max.cnt.soc.intx.pre,
                                                       na.rm = T), 2),
                  max.max.cnt.soc.intx.pre = round(max
                                              (max.cnt.soc.intx.pre,
                                                       na.rm = T), 2),
         # total number of social interactions with all soc. partners   
                 med.strength.pre = round(median(strength.pre,
                                                     na.rm = T), 2),
                 qt25.strength.pre = round(quantile(strength.pre, 
                                                        0.25, na.rm = T), 2),
                 qt75.strength.pre = round(quantile(strength.pre, 
                                                        0.75, na.rm = T), 2),
                 avg.strength.pre = round(mean(strength.pre, 
                                                    na.rm = T),2),
                 stdev.strength.pre = round(sd(strength.pre, 
                                                    na.rm = T), 2),
                 min.strength.pre = round(min(strength.pre,
                                                  na.rm = T), 2),
                 max.strength.pre = round(max(strength.pre,
                                                  na.rm = T), 2),
         # total duration of social intx with all social partners
               med.tot.dur.pre = round(median(tot.dur,  na.rm = T), 2),
               qt25.tot.dur.pre = round(quantile(tot.dur, 0.25, na.rm = T), 2),
               qt75.tot.dur.pre = round(quantile(tot.dur, 0.75, na.rm = T), 2),
               avg.tot.dur.pre = round (mean(tot.dur, na.rm = T),2),
               stdev.tot.dur.pre = round (sd(tot.dur, na.rm = T), 2),
               min.tot.dur.pre = round(min (tot.dur, na.rm = T), 2),
               max.tot.dur.pre = round(max(tot.dur, na.rm = T), 2))
      
      
    ## b) transpose the data frame for easier viewing
      univar_soc_intx_pre <- as.data.frame(t(univar_soc_intx_pre))
      univar_soc_intx_pre <- univar_soc_intx_pre %>%
        rename(c(female = V1,
                 male = V2))
      
    ## c) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_soc_intx_pre.pdf'), height = 11, width = 4)
      grid.table(univar_soc_intx_pre)
      dev.off()
      
    ## d) Descriptive stats social intx. pre-manip variables by age class
      univar_soc_intx_age_pre <- chr15_attrib_pre_df %>%
        group_by(Sex, Age.category) %>%
        # number of unique social partners
        summarise (n.pre = sum(!is.na(degree.pre)),
                   med.degree.pre = round(median(degree.pre,
                                                 na.rm = T), 2),
                   qt25.degree.pre = round(quantile(degree.pre, 0.25, 
                                                    na.rm = T), 2),
                   qt75.degree.pre = round(quantile(degree.pre, 0.75, 
                                                    na.rm = T), 2),
                   avg.degree.pre = round(mean(degree.pre,
                                               na.rm = T),2),
                   stdev.degree.pre = round(sd(degree.pre,
                                               na.rm = T), 2),
                   min.degree.pre = round(min(degree.pre,
                                              na.rm = T), 2),
                   max.degree.pre = round(max(degree.pre,
                                              na.rm = T), 2), 
                # max number of social interactions with any one soc. partner
                   med.max.cnt.soc.intx.pre = round(median
                                                (max.cnt.soc.intx.pre,
                                                      na.rm = T), 2),
                   qt25.max.cnt.soc.intx.pre = round(quantile
                                                (max.cnt.soc.intx.pre, 
                                                       0.25, na.rm = T), 2),
                   qt75.max.cnt.soc.intx.pre = round(quantile
                                                (max.cnt.soc.intx.pre, 
                                                       0.75, na.rm = T), 2),
                   avg.max.cnt.soc.intx.pre = round(mean
                                                (max.cnt.soc.intx.pre, 
                                                       na.rm = T),2),
                   stdev.max.cnt.soc.intx.pre = round(sd
                                                (max.cnt.soc.intx.pre, 
                                                         na.rm = T), 2),
                   min.max.cnt.soc.intx.pre = round(min
                                                (max.cnt.soc.intx.pre,
                                                      na.rm = T), 2),
                   max.max.cnt.soc.intx.pre = round(max
                                                (max.cnt.soc.intx.pre,
                                                      na.rm = T), 2),
                # total number of social interactions with all soc. partners   
                   med.strength.pre = round(median(strength.pre,
                                                      na.rm = T), 2),
                   qt25.strength.pre = round(quantile(strength.pre, 
                                                       0.25, na.rm = T), 2),
                   qt75.strength.pre = round(quantile(strength.pre, 
                                                       0.75, na.rm = T), 2),
                   avg.strength.pre = round(mean(strength.pre, 
                                                       na.rm = T),2),
                   stdev.strength.pre = round(sd(strength.pre, 
                                                         na.rm = T), 2),
                   min.strength.pre = round(min(strength.pre,
                                                      na.rm = T), 2),
                   max.strength.pre = round(max(strength.pre,
                                                      na.rm = T), 2),
              # total duration of social intx with all social partners
                   med.tot.dur.pre = round(median(tot.dur,  na.rm = T), 2),
                   qt25.tot.dur.pre = round(quantile(tot.dur, 0.25, na.rm = T), 2),
                   qt75.tot.dur.pre = round(quantile(tot.dur, 0.75, na.rm = T), 2),
                   avg.tot.dur.pre = round (mean(tot.dur, na.rm = T),2),
                   stdev.tot.dur.pre = round (sd(tot.dur, na.rm = T), 2),
                   min.tot.dur.pre = round(min (tot.dur, na.rm = T), 2),
                   max.tot.dur.pre = round(max(tot.dur, na.rm = T), 2))
      
    ## e) transpose the data frame for easier viewing
      univar_soc_intx_age_pre <- as.data.frame(t(univar_soc_intx_age_pre))
      univar_soc_intx_age_pre <- univar_soc_intx_age_pre %>%
        rename(c(female.sy = V1,
                 female.asy = V2,
                 male.sy = V3,
                 male.asy = V4))    
      
    ## f) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_soc_intx_age_pre.pdf'), height = 11, width = 6)
      grid.table(univar_soc_intx_age_pre)
      dev.off()
      
    ## g) Descriptive stats social intx. variables post manip
      univar_soc_intx_post <- chr15_attrib_post_df %>%
        group_by(Sex) %>%
        # number of unique social partners
        summarise (n.post = sum(!is.na(degree.post)),
                   med.degree.post = round(median(degree.post,
                                                          na.rm = T), 2),
                   qt25.degree.post = round(quantile(degree.post, 
                                                           0.25, na.rm = T), 2),
                   qt75.degree.post = round(quantile(degree.post, 
                                                           0.75, na.rm = T), 2),
                   avg.degree.post = round(mean(degree.post, 
                                                           na.rm = T),2),
                   stdev.degree.post = round(sd(degree.post, 
                                                             na.rm = T), 2),
                   min.degree.post = round(min(degree.post,
                                                          na.rm = T), 2),
                   max.degree.post = round(max(degree.post,
                                                          na.rm = T), 2), 
              # max number of social interactions with any one soc. partner
                   med.max.cnt.soc.intx.post = round(median
                                                (max.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                   qt25.max.cnt.soc.intx.post = round(quantile
                                                  (max.cnt.soc.intx.post, 
                                                       0.25, na.rm = T), 2),
                   qt75.max.cnt.soc.intx.post = round(quantile
                                                  (max.cnt.soc.intx.post, 
                                                       0.75, na.rm = T), 2),
                   avg.max.cnt.soc.intx.post = round(mean
                                                  (max.cnt.soc.intx.post, 
                                                       na.rm = T),2),
                   stdev.max.cnt.soc.intx.post = round(sd
                                                  (max.cnt.soc.intx.post, 
                                                         na.rm = T), 2),
                   min.max.cnt.soc.intx.post = round(min
                                                  (max.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                   max.max.cnt.soc.intx.post = round(max
                                                  (max.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                # total number of social interactions with all soc. partners 
                   med.strength.post = round(median(strength.post,
                                                      na.rm = T), 2),
                   qt25.strength.post = round(quantile(strength.post, 
                                                       0.25, na.rm = T), 2),
                   qt75.strength.post = round(quantile
                                                (strength.post, 
                                                       0.75, na.rm = T), 2),
                   avg.strength.post = round(mean(strength.post, 
                                                       na.rm = T),2),
                   stdev.strength.post = round(sd(strength.post, 
                                                         na.rm = T), 2),
                   min.strength.post = round(min(strength.post,
                                                      na.rm = T), 2),
                   max.strength.post = round(max(strength.post,
                                                      na.rm = T), 2),
              # total duration of social intx with all social partners
                   med.tot.dur.post = round(median(tot.dur,  na.rm = T), 2),
                   qt25.tot.dur.post = round(quantile(tot.dur, 0.25, na.rm = T), 2),
                   qt75.tot.dur.post = round(quantile(tot.dur, 0.75, na.rm = T), 2),
                   avg.tot.dur.post = round (mean(tot.dur, na.rm = T),2),
                   stdev.tot.dur.post = round (sd(tot.dur, na.rm = T), 2),
                   min.tot.dur.post = round(min (tot.dur, na.rm = T), 2),
                   max.tot.dur.post = round(max(tot.dur, na.rm = T), 2))
      
    ## h) transpose the data frame for easier viewing
      univar_soc_intx_post <- as.data.frame(t(univar_soc_intx_post))
      univar_soc_intx_post <- univar_soc_intx_post %>%
        rename(c(female = V1,
                 male = V2))
      
    ## i) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_soc_intx_post.pdf'), height = 11, width = 4)
      grid.table(univar_soc_intx_post)
      dev.off()
      
    ## j) Descriptive stats social intx. post manip variables by age class
      univar_soc_intx_age_post <- chr15_attrib_post_df %>%
        group_by(Sex, Age.category) %>%
        # number of unique social partners
        summarise (n.post = sum(!is.na(degree.post)),
                   med.degree.post = round(median(degree.post,
                                                          na.rm = T), 2),
                   qt25.degree.post = round(quantile(degree.post, 
                                                           0.25, na.rm = T), 2),
                   qt75.degree.post = round(quantile(degree.post, 
                                                           0.75, na.rm = T), 2),
                   avg.degree.post = round(mean(degree.post, 
                                                           na.rm = T),2),
                   stdev.degree.post = round(sd(degree.post, 
                                                             na.rm = T), 2),
                   min.degree.post = round(min(degree.post,
                                                          na.rm = T), 2),
                   max.degree.post = round(max(degree.post,
                                                          na.rm = T), 2), 
              # max number of social interactions with any one soc. partner
                   med.max.cnt.soc.intx.post = round(median
                                                (max.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                   qt25.max.cnt.soc.intx.post = round(quantile
                                                (max.cnt.soc.intx.post, 
                                                       0.25, na.rm = T), 2),
                   qt75.max.cnt.soc.intx.post = round(quantile
                                                (max.cnt.soc.intx.post, 
                                                       0.75, na.rm = T), 2),
                   avg.max.cnt.soc.intx.post = round (mean
                                                (max.cnt.soc.intx.post, 
                                                       na.rm = T),2),
                   stdev.max.cnt.soc.intx.post = round (sd
                                                (max.cnt.soc.intx.post, 
                                                         na.rm = T), 2),
                   min.max.cnt.soc.intx.post = round(min
                                                (max.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                   max.max.cnt.soc.intx.post = round(max
                                                (max.cnt.soc.intx.post,
                                                      na.rm = T), 2),
              # total number of social interactions with all soc. partners  
                   med.strength.post = round(median(strength.post,
                                                      na.rm = T), 2),
                   qt25.strength.post = round(quantile(strength.post, 
                                                       0.25, na.rm = T), 2),
                   qt75.strength.post = round(quantile(strength.post, 
                                                       0.75, na.rm = T), 2),
                   avg.strength.post = round(mean(strength.post, 
                                                       na.rm = T),2),
                   stdev.strength.post = round(sd(strength.post, 
                                                         na.rm = T), 2),
                   min.strength.post = round(min(strength.post,
                                                      na.rm = T), 2),
                   max.strength.post = round(max(strength.post,
                                                      na.rm = T), 2),
                # total duration of social intx with all social partners
                   med.tot.dur.post = round(median(tot.dur,  na.rm = T), 2),
                   qt25.tot.dur.post = round(quantile(tot.dur, 0.25, na.rm = T), 2),
                   qt75.tot.dur.post = round(quantile(tot.dur, 0.75, na.rm = T), 2),
                   avg.tot.dur.post = round (mean(tot.dur, na.rm = T),2),
                   stdev.tot.dur.post = round (sd(tot.dur, na.rm = T), 2),
                   min.tot.dur.post = round(min (tot.dur, na.rm = T), 2),
                   max.tot.dur.post = round(max(tot.dur, na.rm = T), 2))
      
    ## k) transpose the data frame for easier viewing
      univar_soc_intx_age_post <- as.data.frame(t
                                        (univar_soc_intx_age_post))
      univar_soc_intx_age_post <- univar_soc_intx_age_post %>%
        rename(c(female.sy = V1,
                 female.asy = V2,
                 male.sy = V3,
                 male.asy = V4))    
      
    ## l) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_soc_intx_age_post.pdf'), height = 11, width = 6)
      grid.table(univar_soc_intx_age_post)
      dev.off()  
      
      
  ### 3.3 Univariate descriptive stats for reproduction data 
    ## a) Frequency distribution of female total fecundity
     f_total_repro_freq_dist_plot <- 
        # raw data plot
        ggplot(data = subset(chr15_attrib_pre_df
                             , Sex == 'f'), aes(x = total.fecundity)) +
        geom_histogram(binwidth = 1, color = 'black', fill = 'steelblue4') +
                # scale the density plot line
        geom_density(aes(y=after_stat(density)/0.025), color = 'red', 
                     linewidth = 3) +
       # Titles, axes, and legends
       labs(title = 'Frequency distribution of female total fecundity') +
       theme(plot.title = element_text(hjust = 0.5, size = 14)) + # center title
       # bold and size title and axes labels
       theme(text = element_text(size=18, face = 'bold')) +
       #theme(legend.position = 'none') +
       theme(axis.ticks = element_blank()) + # remove axis ticks
       # remove background color
       theme(panel.background = element_rect(fill = 'white')) +
       # add major axes
       theme(axis.line = element_line(colour = 'black',
                                      size = 0.5, linetype = 'solid')) +
       
       # change axes font style, color, size, angle, margin, and legend
       theme(axis.text.x = element_text(face='bold', color='black', 
                                        size=18, angle=0,
                                        margin = margin(t = 10, r = 0, 
                                                        b = 10, l = 0)),
             axis.text.y = element_text(face='bold', color='black', 
                                        size=18, angle=0, 
                                        margin = margin(t = 0, r = 10, 
                                                        b = 0, l = 0)),
             legend.title = element_text(size = 16),
             legend.text = element_text(size=14),
             legend.position = 'none',
             legend.key = element_blank()) +
       xlab(expression(bold('female total fecundity (from 3 clutches)'))) +
       ylab(expression(bold('count of individuals'))) 
     
    ## b) View plot
     print(f_total_repro_freq_dist_plot)
     
    ## c) Save Plot
     # use ggsave to save the plot
     ggsave('f_total_repro_freq_dist_plot.pdf', 
            plot = f_total_repro_freq_dist_plot, 
            device = NULL,
            path = paste0(here(),'/output'), 
            scale = 1, width = 7,
            height = 6,
            units = c('in'), dpi = 300, limitsize = TRUE) 
     
    ## d) Frequency distribution of male total paternity
     m_total_repro_freq_dist_plot <- 
       # raw data plot
       ggplot(data = subset(chr15_attrib_pre_df
                            , Sex == 'm'), aes(x = total.paternity)) +
       geom_histogram(binwidth = 1, color = 'black', fill = 'wheat3') +
       # scale the density plot line
       geom_density(aes(y=after_stat(density)/0.025), color = 'red', 
                    linewidth = 3) +
       # Titles, axes, and legends
       labs(title = 'Frequency distribution of male total paternity') +
       theme(plot.title = element_text(hjust = 0.5, size = 14)) + # center title
       # bold and size title and axes labels
       theme(text = element_text(size=18, face = 'bold')) +
       #theme(legend.position = 'none') +
       theme(axis.ticks = element_blank()) + # remove axis ticks
       # remove background color
       theme(panel.background = element_rect(fill = 'white')) +
       # add major axes
       theme(axis.line = element_line(colour = 'black',
                                      size = 0.5, linetype = 'solid')) +
       
       # change axes font style, color, size, angle, margin, and legend
       theme(axis.text.x = element_text(face='bold', color='black', 
                                        size=18, angle=0,
                                        margin = margin(t = 10, r = 0, 
                                                        b = 10, l = 0)),
             axis.text.y = element_text(face='bold', color='black', 
                                        size=18, angle=0, 
                                        margin = margin(t = 0, r = 10, 
                                                        b = 0, l = 0)),
             legend.title = element_text(size = 16),
             legend.text = element_text(size=14),
             legend.position = 'none',
             legend.key = element_blank()) +
       xlab(expression(bold('male total paternity (from 3 clutches)'))) +
       ylab(expression(bold('count of individuals'))) 
     
    ## e) View plot
     print(m_total_repro_freq_dist_plot)
     
    ## f) Save Plot
     # use ggsave to save the plot
     ggsave('m_total_repro_freq_dist_plot.pdf', 
            plot = m_total_repro_freq_dist_plot, 
            device = NULL,
            path = paste0(here(),'/output'), 
            scale = 1, width = 7,
            height = 6,
            units = c('in'), dpi = 300, limitsize = TRUE) 
     
     
                   
  ### 3.4 Descriptive stats for reproductive success measures all ages
    ## a) Summarize female total fecundity from all clutches 
      univar_f_total_offspring <- chr15_attrib_pre_df %>%
        ungroup() %>%
        filter(Sex == 'f') %>%
        # total fecundity
        summarise (n.total.fecundity = sum(!is.na(total.fecundity)),
                   med.total.fecundity = round(median(total.fecundity,
                                                       na.rm = T), 2),
                   qt25.total.fecundity = round(quantile(total.fecundity, 
                                                          0.25, na.rm = T), 2),
                   qt75.total.fecundity = round(quantile(total.fecundity, 
                                                          0.75, na.rm = T), 2),
                   avg.total.fecundity = round (mean(total.fecundity, 
                                                      na.rm = T),2),
                   stdev.total.fecundity = round (sd(total.fecundity, 
                                                      na.rm = T), 2),
                   min.total.fecundity = round(min(total.fecundity,
                                                    na.rm = T), 2),
                   max.total.fecundity = round(max(total.fecundity,
                                                    na.rm = T), 2))
      
    ## b) Summarize nest attempt 1 male total paternity
      univar_m_clutch_1_pat <- chr15_attrib_pre_df %>%
        ungroup() %>%
        filter(Sex == 'm') %>%
        # nest attempt 1 total paternity
        summarise(n.attmpt.1.tot.pat = sum(!is.na(attmpt.1.tot.pat)),
                  med.attmpt.1.tot.pat = round(median(attmpt.1.tot.pat,
                                                   na.rm = T), 2),
                  qt25.attmpt.1.tot.pat = round(quantile(attmpt.1.tot.pat, 
                                                      0.25, na.rm = T), 2),
                  qt75.attmpt.1.tot.pat = round(quantile(attmpt.1.tot.pat, 
                                                      0.75, na.rm = T), 2),
                  avg.attmpt.1.tot.pat = round (mean(attmpt.1.tot.pat, 
                                                  na.rm = T),2),
                  stdev.attmpt.1.tot.pat = round (sd(attmpt.1.tot.pat, 
                                                  na.rm = T), 2),
                  min.attmpt.1.tot.pat = round(min(attmpt.1.tot.pat,
                                                na.rm = T), 2),
                  max.attmpt.1.tot.pat = round(max(attmpt.1.tot.pat,
                                                na.rm = T), 2))
      
    ## c) Summarize nest attempt 2 male total paternity 
      univar_m_clutch_2_pat <- chr15_attrib_post_df %>%
        ungroup() %>%
        filter(Sex == 'm') %>%
        # nest attempt 1 total paternity
        summarise(n.attmpt.2.tot.pat = sum(!is.na(attmpt.2.tot.pat)),
                  med.attmpt.2.tot.pat = round(median(attmpt.2.tot.pat,
                                                   na.rm = T), 2),
                  qt25.attmpt.2.tot.pat = round(quantile(attmpt.2.tot.pat, 
                                                      0.25, na.rm = T), 2),
                  qt75.attmpt.2.tot.pat = round(quantile(attmpt.2.tot.pat, 
                                                      0.75, na.rm = T), 2),
                  avg.attmpt.2.tot.pat = round (mean(attmpt.2.tot.pat, 
                                                  na.rm = T),2),
                  stdev.attmpt.2.tot.pat = round (sd(attmpt.2.tot.pat, 
                                                  na.rm = T), 2),
                  min.attmpt.2.tot.pat = round(min(attmpt.2.tot.pat,
                                                na.rm = T), 2),
                  max.attmpt.2.tot.pat = round(max(attmpt.2.tot.pat,
                                                na.rm = T), 2))
      
    ## d) Summarize total number of male paternity in second clutch
      univar_m_clutch_all_tot_pat <- chr15_attrib_pre_df %>%
        ungroup() %>%
        filter(Sex == 'm') %>%
        # first clutch paternity
        summarise(n.total.paternity = sum(!is.na(total.paternity)),
                  med.total.paternity = round(median(total.paternity,
                                              na.rm = T), 2),
                  qt25.total.paternity = round(quantile(total.paternity, 
                                                 0.25, na.rm = T), 2),
                  qt75.total.paternity = round(quantile(total.paternity, 
                                                 0.75, na.rm = T), 2),
                  avg.total.paternity = round (mean(total.paternity, 
                                             na.rm = T),2),
                  stdev.total.paternity = round (sd(total.paternity, 
                                             na.rm = T), 2),
                  min.total.paternity = round(min(total.paternity,
                                           na.rm = T), 2),
                  max.total.paternity = round(max(total.paternity,
                                           na.rm = T), 2))
      
  ### 3.5 Descriptive stats for reproductive success measures by age category
    ## a) Summarize female total fecundity from all clutches by age
      univar_f_total_offspring_by_age <- chr15_attrib_pre_df %>%
        group_by(Age.category) %>%
        filter(Sex == 'f') %>%
        # total fecundity
        summarise (n.total.fecundity = sum(!is.na(total.fecundity)),
                   med.total.fecundity = round(median(total.fecundity,
                                                      na.rm = T), 2),
                   qt25.total.fecundity = round(quantile(total.fecundity, 
                                                         0.25, na.rm = T), 2),
                   qt75.total.fecundity = round(quantile(total.fecundity, 
                                                         0.75, na.rm = T), 2),
                   avg.total.fecundity = round (mean(total.fecundity, 
                                                     na.rm = T),2),
                   stdev.total.fecundity = round (sd(total.fecundity, 
                                                     na.rm = T), 2),
                   min.total.fecundity = round(min(total.fecundity,
                                                   na.rm = T), 2),
                   max.total.fecundity = round(max(total.fecundity,
                                                   na.rm = T), 2))
      
    ## b) Summarize nest attempt 1 male total paternity by age
      univar_m_clutch_1_pat_by_age <- chr15_attrib_pre_df %>%
        group_by(Age.category) %>%
        filter(Sex == 'm') %>%
        # nest attempt 1 total paternity
        summarise(n.attmpt.1.tot.pat = sum(!is.na(attmpt.1.tot.pat)),
                  med.attmpt.1.tot.pat = round(median(attmpt.1.tot.pat,
                                                      na.rm = T), 2),
                  qt25.attmpt.1.tot.pat = round(quantile(attmpt.1.tot.pat, 
                                                         0.25, na.rm = T), 2),
                  qt75.attmpt.1.tot.pat = round(quantile(attmpt.1.tot.pat, 
                                                         0.75, na.rm = T), 2),
                  avg.attmpt.1.tot.pat = round (mean(attmpt.1.tot.pat, 
                                                     na.rm = T),2),
                  stdev.attmpt.1.tot.pat = round (sd(attmpt.1.tot.pat, 
                                                     na.rm = T), 2),
                  min.attmpt.1.tot.pat = round(min(attmpt.1.tot.pat,
                                                   na.rm = T), 2),
                  max.attmpt.1.tot.pat = round(max(attmpt.1.tot.pat,
                                                   na.rm = T), 2))
      
    ## c) Summarize nest attempt 2 male total paternity by age
      univar_m_clutch_2_pat_by_age <- chr15_attrib_post_df %>%
        group_by(Age.category) %>%
        filter(Sex == 'm') %>%
        # nest attempt 1 total paternity
        summarise(n.attmpt.2.tot.pat = sum(!is.na(attmpt.2.tot.pat)),
                  med.attmpt.2.tot.pat = round(median(attmpt.2.tot.pat,
                                                      na.rm = T), 2),
                  qt25.attmpt.2.tot.pat = round(quantile(attmpt.2.tot.pat, 
                                                         0.25, na.rm = T), 2),
                  qt75.attmpt.2.tot.pat = round(quantile(attmpt.2.tot.pat, 
                                                         0.75, na.rm = T), 2),
                  avg.attmpt.2.tot.pat = round (mean(attmpt.2.tot.pat, 
                                                     na.rm = T),2),
                  stdev.attmpt.2.tot.pat = round (sd(attmpt.2.tot.pat, 
                                                     na.rm = T), 2),
                  min.attmpt.2.tot.pat = round(min(attmpt.2.tot.pat,
                                                   na.rm = T), 2),
                  max.attmpt.2.tot.pat = round(max(attmpt.2.tot.pat,
                                                   na.rm = T), 2))
      
    ## d) Summarize total number of male paternity in second clutch by age
      univar_m_clutch_all_tot_pat_by_age <- chr15_attrib_pre_df %>%
        group_by(Age.category) %>%
        filter(Sex == 'm') %>%
        # first clutch paternity
        summarise(n.total.paternity = sum(!is.na(total.paternity)),
                  med.total.paternity = round(median(total.paternity,
                                                     na.rm = T), 2),
                  qt25.total.paternity = round(quantile(total.paternity, 
                                                        0.25, na.rm = T), 2),
                  qt75.total.paternity = round(quantile(total.paternity, 
                                                        0.75, na.rm = T), 2),
                  avg.total.paternity = round (mean(total.paternity, 
                                                    na.rm = T),2),
                  stdev.total.paternity = round (sd(total.paternity, 
                                                    na.rm = T), 2),
                  min.total.paternity = round(min(total.paternity,
                                                  na.rm = T), 2),
                  max.total.paternity = round(max(total.paternity,
                                                  na.rm = T), 2))
      

    