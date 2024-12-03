################################################################################
#############        The role of age and plumage traits as         #############  
#############       determinants of reproductive success and       #############
#############       the mediating role of social interactions      #############
#############                                                      #############
#############    5. Descriptive statistics and data exploration    #############
#############                                                      #############
#############                  By: Zach Laubach                    #############
#############                created: 21 Aug 2024                  #############
#############              last updated: 3 Dec 2024                #############
################################################################################


  ### PURPOSE: Calculate descriptive statistics and explore data 
  
  # Code Blocks
    # 1: Configure work space
    # 2: Import data 
    # 3: Node level descriptive stats  

  


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
      
    # load assortnet package
      library ('assortnet')  
    
    # # load corrr package
    #   library ('corrr')
    
    # load lubridate package
      library ('lubridate')
    
    # load gridExtra packages
      library ('gridExtra')  
    
    # load smplot2 packages   
      library('smplot2')
      
    # load chisq.posthoc.test package  
      library('chisq.posthoc.test')
      
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
    ## a) Load node and dyad level data for CHR 2015
      load(here('data/5_6_chr15_mediation_data.RData'))
      
   

###############################################################################
##############          3. Node level descriptive stats          ##############
###############################################################################

  ### 3.1 Univariate descriptive stats for node-level attribute data
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
      univar_node_morph <- chr15_attrib_pre_df %>%
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
                   ) %>%
        ungroup()
      
    ## e) transpose the data frame for easier viewing
      univar_node_morph <- as.data.frame(t(univar_node_morph))
      univar_node_morph <- univar_node_morph %>%
        rename(c(female = V1,
                 male = V2))
      
    ## f) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_node_morph.pdf'), height = 13, width = 4)
      grid.table(univar_node_morph)
      dev.off()

    ## g) Descriptive stats morphology variables
      univar_node_R_brigh_post <- chr15_attrib_post_df %>%
        ungroup() %>%
        filter(Color.manipulation. == 'y') %>%  
      # breast average brightness post-manip
      summarise(n.R.bright.treat.and.orig = sum(!is.na
                                      (R.bright.treat.and.orig)),
      avg.R.bright.treat.and.orig = round (mean
                                           (R.bright.treat.and.orig, 
                                             na.rm = T),2),
      stdev.R.bright.treat.and.orig = round (sd
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
      
  ### 3.2  Univariate descriptive stats for node-level social intx. data
      ## a) Descriptive stats social intx. variables pre-manip
      univar_node_soc_intx_pre <- chr15_attrib_pre_df %>%
        group_by(Sex) %>%
        # number of unique social partners
        summarise (n.n.unique.soc.partner.pre = sum(!is.na
                                                    (n.unique.soc.partner.pre)),
                   med.n.unique.soc.partner.pre = round(median
                                                    (n.unique.soc.partner.pre,
                                                   na.rm = T), 2),
                   qt25.n.unique.soc.partner.pre = round(quantile
                                                    (n.unique.soc.partner.pre, 
                                                          0.25, na.rm = T), 2),
                   qt75.n.unique.soc.partner.pre = round(quantile
                                                    (n.unique.soc.partner.pre, 
                                                          0.75, na.rm = T), 2),
                   avg.n.unique.soc.partner.pre = round (mean
                                                    (n.unique.soc.partner.pre, 
                                             na.rm = T),2),
                   stdev.n.unique.soc.partner.pre = round (sd
                                                    (n.unique.soc.partner.pre, 
                                             na.rm = T), 2),
                   min.n.unique.soc.partner.pre = round(min
                                                    (n.unique.soc.partner.pre,
                                           na.rm = T), 2),
                   max.n.unique.soc.partner.pre = round(max
                                                    (n.unique.soc.partner.pre,
                                           na.rm = T), 2), 
         # max number of social interactions with any one soc. partner
                  n.max.cnt.soc.intx.pre = sum(!is.na
                                              (max.cnt.soc.intx.pre)),
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
                 n.tot.cnt.soc.intx.pre = sum(!is.na
                                          (tot.cnt.soc.intx.pre)),
                 med.tot.cnt.soc.intx.pre = round(median
                                          (tot.cnt.soc.intx.pre,
                                                     na.rm = T), 2),
                 qt25.tot.cnt.soc.intx.pre = round(quantile
                                          (tot.cnt.soc.intx.pre, 
                                                        0.25, na.rm = T), 2),
                 qt75.tot.cnt.soc.intx.pre = round(quantile
                                          (tot.cnt.soc.intx.pre, 
                                                        0.75, na.rm = T), 2),
                 avg.tot.cnt.soc.intx.pre = round (mean
                                          (tot.cnt.soc.intx.pre, 
                                                    na.rm = T),2),
                 stdev.tot.cnt.soc.intx.pre = round (sd
                                          (tot.cnt.soc.intx.pre, 
                                                    na.rm = T), 2),
                 min.tot.cnt.soc.intx.pre = round(min
                                          (tot.cnt.soc.intx.pre,
                                                  na.rm = T), 2),
                 max.tot.cnt.soc.intx.pre = round(max
                                          (tot.cnt.soc.intx.pre,
                                                  na.rm = T), 2))
      
      
    ## b) transpose the data frame for easier viewing
      univar_node_soc_intx_pre <- as.data.frame(t(univar_node_soc_intx_pre))
      univar_node_soc_intx_pre <- univar_node_soc_intx_pre %>%
        rename(c(female = V1,
                 male = V2))
      
    ## c) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_node_soc_intx_pre.pdf'), height = 11, width = 4)
      grid.table(univar_node_soc_intx_pre)
      dev.off()
      
    ## d) Descriptive stats social intx. pre-manip variables by age class
      univar_node_soc_intx_age_pre <- chr15_attrib_pre_df %>%
        group_by(Sex, Age.category) %>%
        # number of unique social partners
        summarise (n.n.unique.soc.partner.pre = sum(!is.na
                                                (n.unique.soc.partner.pre)),
                   med.n.unique.soc.partner.pre = round(median
                                                (n.unique.soc.partner.pre,
                                                          na.rm = T), 2),
                   qt25.n.unique.soc.partner.pre = round(quantile
                                                (n.unique.soc.partner.pre, 
                                                           0.25, na.rm = T), 2),
                   qt75.n.unique.soc.partner.pre = round(quantile
                                                (n.unique.soc.partner.pre, 
                                                           0.75, na.rm = T), 2),
                   avg.n.unique.soc.partner.pre = round (mean
                                                (n.unique.soc.partner.pre, 
                                                           na.rm = T),2),
                   stdev.n.unique.soc.partner.pre = round (sd
                                                (n.unique.soc.partner.pre, 
                                                             na.rm = T), 2),
                   min.n.unique.soc.partner.pre = round(min
                                                (n.unique.soc.partner.pre,
                                                          na.rm = T), 2),
                   max.n.unique.soc.partner.pre = round(max
                                                (n.unique.soc.partner.pre,
                                                          na.rm = T), 2), 
                   # max number of social interactions with any one soc. partner
                   n.max.cnt.soc.intx.pre = sum(!is.na
                                                (max.cnt.soc.intx.pre)),
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
                   n.tot.cnt.soc.intx.pre = sum(!is.na
                                                (tot.cnt.soc.intx.pre)),
                   med.tot.cnt.soc.intx.pre = round(median
                                                (tot.cnt.soc.intx.pre,
                                                      na.rm = T), 2),
                   qt25.tot.cnt.soc.intx.pre = round(quantile
                                                (tot.cnt.soc.intx.pre, 
                                                       0.25, na.rm = T), 2),
                   qt75.tot.cnt.soc.intx.pre = round(quantile
                                                (tot.cnt.soc.intx.pre, 
                                                       0.75, na.rm = T), 2),
                   avg.tot.cnt.soc.intx.pre = round (mean
                                                (tot.cnt.soc.intx.pre, 
                                                       na.rm = T),2),
                   stdev.tot.cnt.soc.intx.pre = round (sd
                                                (tot.cnt.soc.intx.pre, 
                                                         na.rm = T), 2),
                   min.tot.cnt.soc.intx.pre = round(min
                                                (tot.cnt.soc.intx.pre,
                                                      na.rm = T), 2),
                   max.tot.cnt.soc.intx.pre = round(max
                                                (tot.cnt.soc.intx.pre,
                                                      na.rm = T), 2))
      
    ## e) transpose the data frame for easier viewing
      univar_node_soc_intx_age_pre <- as.data.frame(t(univar_node_soc_intx_age_pre))
      univar_node_soc_intx_age_pre <- univar_node_soc_intx_age_pre %>%
        rename(c(female.sy = V1,
                 female.asy = V2,
                 male.sy = V3,
                 male.asy = V4))    
      
    ## f) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_node_soc_intx_age_pre.pdf'), height = 11, width = 6)
      grid.table(univar_node_soc_intx_age_pre)
      dev.off()
      
    ## g) Descriptive stats social intx. variables post manip
      univar_node_soc_intx_post <- chr15_attrib_post_df %>%
        group_by(Sex) %>%
        # number of unique social partners
        summarise (n.n.unique.soc.partner.post = sum(!is.na
                                                    (n.unique.soc.partner.post)),
                   med.n.unique.soc.partner.post = round(median
                                                    (n.unique.soc.partner.post,
                                                          na.rm = T), 2),
                   qt25.n.unique.soc.partner.post = round(quantile
                                                    (n.unique.soc.partner.post, 
                                                           0.25, na.rm = T), 2),
                   qt75.n.unique.soc.partner.post = round(quantile
                                                    (n.unique.soc.partner.post, 
                                                           0.75, na.rm = T), 2),
                   avg.n.unique.soc.partner.post = round (mean
                                                    (n.unique.soc.partner.post, 
                                                           na.rm = T),2),
                   stdev.n.unique.soc.partner.post = round (sd
                                                    (n.unique.soc.partner.post, 
                                                             na.rm = T), 2),
                   min.n.unique.soc.partner.post = round(min
                                                    (n.unique.soc.partner.post,
                                                          na.rm = T), 2),
                   max.n.unique.soc.partner.post = round(max
                                                    (n.unique.soc.partner.post,
                                                          na.rm = T), 2), 
                   # max number of social interactions with any one soc. partner
                   n.max.cnt.soc.intx.post = sum(!is.na
                                                (max.cnt.soc.intx.post)),
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
                   n.tot.cnt.soc.intx.post = sum(!is.na
                                                (tot.cnt.soc.intx.post)),
                   med.tot.cnt.soc.intx.post = round(median
                                                (tot.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                   qt25.tot.cnt.soc.intx.post = round(quantile
                                                (tot.cnt.soc.intx.post, 
                                                       0.25, na.rm = T), 2),
                   qt75.tot.cnt.soc.intx.post = round(quantile
                                                (tot.cnt.soc.intx.post, 
                                                       0.75, na.rm = T), 2),
                   avg.tot.cnt.soc.intx.post = round (mean
                                                (tot.cnt.soc.intx.post, 
                                                       na.rm = T),2),
                   stdev.tot.cnt.soc.intx.post = round (sd
                                                (tot.cnt.soc.intx.post, 
                                                         na.rm = T), 2),
                   min.tot.cnt.soc.intx.post = round(min
                                                (tot.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                   max.tot.cnt.soc.intx.post = round(max
                                                (tot.cnt.soc.intx.post,
                                                      na.rm = T), 2))
      
    ## h) transpose the data frame for easier viewing
      univar_node_soc_intx_post <- as.data.frame(t(univar_node_soc_intx_post))
      univar_node_soc_intx_post <- univar_node_soc_intx_post %>%
        rename(c(female = V1,
                 male = V2))
      
    ## i) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_node_soc_intx_post.pdf'), height = 11, width = 4)
      grid.table(univar_node_soc_intx_post)
      dev.off()
      
    ## j) Descriptive stats social intx. post manip variables by age class
      univar_node_soc_intx_age_post <- chr15_attrib_post_df %>%
        group_by(Sex, Age.category) %>%
        # number of unique social partners
        summarise (n.n.unique.soc.partner.post = sum(!is.na
                                                  (n.unique.soc.partner.post)),
                   med.n.unique.soc.partner.post = round(median
                                                  (n.unique.soc.partner.post,
                                                          na.rm = T), 2),
                   qt25.n.unique.soc.partner.post = round(quantile
                                                  (n.unique.soc.partner.post, 
                                                           0.25, na.rm = T), 2),
                   qt75.n.unique.soc.partner.post = round(quantile
                                                  (n.unique.soc.partner.post, 
                                                           0.75, na.rm = T), 2),
                   avg.n.unique.soc.partner.post = round (mean
                                                  (n.unique.soc.partner.post, 
                                                           na.rm = T),2),
                   stdev.n.unique.soc.partner.post = round (sd
                                                  (n.unique.soc.partner.post, 
                                                             na.rm = T), 2),
                   min.n.unique.soc.partner.post = round(min
                                                  (n.unique.soc.partner.post,
                                                          na.rm = T), 2),
                   max.n.unique.soc.partner.post = round(max
                                                  (n.unique.soc.partner.post,
                                                          na.rm = T), 2), 
                   # max number of social interactions with any one soc. partner
                   n.max.cnt.soc.intx.post = sum(!is.na
                                                (max.cnt.soc.intx.post)),
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
                   n.tot.cnt.soc.intx.post = sum(!is.na
                                                (tot.cnt.soc.intx.post)),
                   med.tot.cnt.soc.intx.post = round(median
                                                (tot.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                   qt25.tot.cnt.soc.intx.post = round(quantile
                                                (tot.cnt.soc.intx.post, 
                                                       0.25, na.rm = T), 2),
                   qt75.tot.cnt.soc.intx.post = round(quantile
                                                (tot.cnt.soc.intx.post, 
                                                       0.75, na.rm = T), 2),
                   avg.tot.cnt.soc.intx.post = round (mean
                                                (tot.cnt.soc.intx.post, 
                                                       na.rm = T),2),
                   stdev.tot.cnt.soc.intx.post = round (sd
                                                (tot.cnt.soc.intx.post, 
                                                         na.rm = T), 2),
                   min.tot.cnt.soc.intx.post = round(min
                                                (tot.cnt.soc.intx.post,
                                                      na.rm = T), 2),
                   max.tot.cnt.soc.intx.post = round(max
                                                (tot.cnt.soc.intx.post,
                                                      na.rm = T), 2))
      
    ## k) transpose the data frame for easier viewing
      univar_node_soc_intx_age_post <- as.data.frame(t
                                        (univar_node_soc_intx_age_post))
      univar_node_soc_intx_age_post <- univar_node_soc_intx_age_post %>%
        rename(c(female.sy = V1,
                 female.asy = V2,
                 male.sy = V3,
                 male.asy = V4))    
      
    ## l) save the data frame of summary stats as a pdf into output file
      pdf(here('output/univar_node_soc_intx_age_post.pdf'), height = 11, width = 6)
      grid.table(univar_node_soc_intx_age_post)
      dev.off()  
      
      
  ### 3.3 Univariate descriptive stats for node-level reproduction data 
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
      

    