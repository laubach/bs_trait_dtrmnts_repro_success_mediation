################################################################################
#############       The role of social and phenotypic dyadic       #############  
#############           relationships as determinants of           #############
#############                 reproductive success                 #############
#############                                                      #############
############# 5. Node descriptive statistics and data exploration  #############
#############                                                      #############
#############                  By: Zach Laubach                    #############
#############               created: 21 March 2025                 #############
#############              last updated: 14 Oct 2025               #############
################################################################################


  ### PURPOSE: Calculate descriptive statistics and explore data 
  
  # Code Blocks
    # 1.  Configure work space
    # 2.  Import data 
    # 3.  Descrpt. stats demog., plum., morph. & phys.
    # 4.  Descriptive stats social intx. 
    # 5.  Descriptive stats reproductive success
    # 6.  Female-Male social interaction plots
    # 7.  Female-Female social interaction plots
    # 8.  Male-Male social interaction plots
    # 9.  Female fecundity plots 
    # 10. Male paternity plots

  


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
      
      # load lubridate package
        library ('lubridate')
 
    ## b) Graphing/visualization packages        
      # load geomtextpath package
        library('geomtextpath')
      
      # load gridExtra package 
        library ('gridExtra')  
      
      # load rempsyc to make publication ready tables
      # see https://cran.r-project.org/web/packages/rempsyc/vignettes/table.html
        library('rempsyc')
      
        pkgs <- c('flextable', 'broom', 'report', 'effectsize')
        install_if_not_installed(pkgs)
        
    # ## c) Network packages   
    #   # load igraph package
    #     library ('igraph')  
    #     
    #   # load intergraph package: used in quadratic assignment procedure (QAP)
    #     library ('intergraph')
    #     
    #   # load statnet package: used in QAP
    #     library ('statnet')
      
    ## d) Other packages      
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
      source_path <- paste('~/WD/Git/source_code/')
    
    ## b) all_char_to_lower function
      source(file = paste0(source_path, 'all_char_to_lower.R'))
    
    ## c) format_var_names function
      source(file = paste0(source_path, 'format_var_names.R'))
      
      

###############################################################################
##############                  2. Import data                   ##############
###############################################################################    

  ### 2.1 Import sample data files
      ## a) Load node data for CHR 2015
      load(here('data/5_6_chr_node_df.RData'))
      
   

###############################################################################
##############  3. Descriptive stats for demography and plumage  ##############
###############################################################################

  ### 3.1 Descriptive stats for age and sex
    ## a) Total sample size
        n <- sum(!is.na(chr15_attrib_df$Band.ID))
        
    ## b) Descriptive stats for sex  
      chr15_attrib_df %>%
        group_by(Sex) %>%
        summarise (n.Sex = sum(!is.na(Sex)))%>%
        mutate(freq = n.Sex / sum(n.Sex)) %>%
        ungroup() 
        
    ## c) Descriptive stats for Age.category  
      sex_age_smmry <- chr15_attrib_df %>%
        group_by(Sex, Age.category) %>%
        summarise(n = sum(!is.na(Age.category))) %>%
        ungroup() %>%
        mutate(freq = round(n / sum(n), 2)) 
      
    ## d) transpose the data frame for easier viewing
      sex_age_smmry <- as.data.frame(t(sex_age_smmry))
      
    ## e) tidy sex_age_smmry
      sex_age_smmry <- sex_age_smmry %>%
        # mutate(Sex = 'f') %>%
        # relocate(Sex, .before = Age.category) %>%
        mutate(type = 'demography') %>%
        relocate(type, .before = V1)
      

  ### 3.2 Descriptive stats for plumage traits
    ## a) Descriptive stats plumage variables
      plumage_smmry <- chr15_attrib_df %>%
        group_by(Sex, Age.category) %>%
        # throat average brightness
        summarise (#n.T.avg.bright = sum(!is.na(T_avg.bright)),
                   avg.T.avg.bright = round (mean(T_avg.bright, 
                                                  na.rm = T), 2),
                   stdev.T.avg.bright = round (sd(T_avg.bright, 
                                                  na.rm = T), 2),
                   # med.T.avg.bright = round(median(T_avg.bright,
                   #                                 na.rm = T), 2),
                   # min.T.avg.bright = round(min(T_avg.bright,
                   #                              na.rm = T), 2),
                   # max.T.avg.bright = round(max(T_avg.bright,
                   #                              na.rm = T), 2),
        # breast average brightness pre-manip
                   #n.R.avg.bright = sum(!is.na(R_avg.bright)),
                   avg.R.avg.bright = round (mean(R_avg.bright, 
                                                  na.rm = T), 2),
                   stdev.R.avg.bright = round (sd(R_avg.bright, 
                                                  na.rm = T), 2),
                   # med.R.avg.bright = round(median(R_avg.bright,
                   #                                 na.rm = T), 2),
                   # min.R.avg.bright = round(min(R_avg.bright,
                   #                              na.rm = T), 2),
                   # max.R.avg.bright = round(max(R_avg.bright,
                   #                              na.rm = T), 2),
        # belly average brightness pre-manip
                   #n.B.avg.bright = sum(!is.na(B_avg.bright)),
                   avg.B.avg.bright = round (mean(B_avg.bright, 
                                                  na.rm = T), 2),
                   stdev.B.avg.bright = round (sd(B_avg.bright, 
                                                  na.rm = T), 2),
                   # med.B.avg.bright = round(median(B_avg.bright,
                   #                                 na.rm = T), 2),
                   # min.B.avg.bright = round(min(B_avg.bright,
                   #                              na.rm = T), 2),
                   # max.B.avg.bright = round(max(B_avg.bright,
                   #                              na.rm = T), 2),
        # experimental throat average brightness
                   n.Post.R.avg.bright = sum(!is.na(Post.R_avg.bright)),
                   avg.Post.R.avg.bright = round(mean(Post.R_avg.bright, 
                                                   na.rm = T), 2),
                   stdev.Post.R.avg.bright = round(sd(Post.R_avg.bright, 
                                                     na.rm = T), 2),
                   # med.Post.R.avg.bright = round(median(Post.R_avg.bright,
                   #                                 na.rm = T), 2)
                   # min.Post.R.avg.bright = round(min(Post.R_avg.bright,
                   #                                 na.rm = T), 2),
                   # max.Post.R.avg.bright = round(max(Post.R_avg.bright,
                   #                                 na.rm = T), 2)
        # tail streamer length
                  n.TS = sum(!is.na(Mean.TS)),
                  avg.TS = round (mean(Mean.TS, 
                                       na.rm = T), 2),
                  stdev.TS = round (sd(Mean.TS, 
                                       na.rm = T), 2)
                  # med.TS = round(median(Mean.TS,
                  #                      na.rm = T), 2),
                  # min.TS = round(min(Mean.TS,
                  #                              na.rm = T), 2),
                  # max.TS = round(max(Mean.TS,
                  #                              na.rm = T), 2),
        ) %>%
        ungroup()
      
    ## b) transpose the data frame for easier viewing
      plumage_smmry <- as.data.frame(t(plumage_smmry))
      
    ## c) tidy plumage_smmry
      plumage_smmry <- plumage_smmry %>%
        # mutate(Sex = 'f') %>%
        # relocate(Sex, .before = Age.category) %>%
        mutate(type = 'plumage') %>%
        relocate(type, .before = V1) %>%
        slice(c(-1, -2))
      

  ### 3.3 Format demography and plumage summary table for publication  
    ## a) row bind sex_age_smmry and plumage_smmry
      dem_plum_smmry <- rbind(sex_age_smmry, plumage_smmry)

    ## c) rename variables   
      dem_plum_smmry <- dem_plum_smmry %>%
        rename(c(`female.sy` = V1,
                 `female.asy` = V2,
                 `male.sy` = V3,
                 `male.asy` = V4)) %>%
        slice(c(-1, -2)) %>%
        rownames_to_column(' ')
      
    ## d) tidy dem_plum_smmry
      # rename blank column
      colnames(dem_plum_smmry)[1] <- 'variable'
      
      dem_plum_smmry <- dem_plum_smmry %>%
        # mutate(Sex = 'f') %>%
        # relocate(Sex, .before = Age.category) %>%
       # colnames(dem_plum_morph_phys_smmry)[1] <- 'variable' %>%
        relocate(type, .before = variable)
      
    ## e) format table creating flextable object 
      dem_plum_smmry_frmt <- 
        nice_table(dem_plum_smmry, separate.header = T)
      
    ## f) preview the table as a word doc    
      #print(dem_plum_morph_phys_smmry_frmt, preview = 'docx')
      
    ## f) save dem_plum_morph_phys_smmry_frmt as a word doc
      flextable::save_as_docx(dem_plum_smmry_frmt, 
              path = here('output/dem_plum_smmry_frmt.docx'))      
      

###############################################################################
##############         4. Descriptive stats social intx.         ##############
###############################################################################          
      
           
  ### 4.1 Descriptive stats for female-male node level social interactions   
    ## a) Descriptive stats female-male social interactions
      soc_intx_fxm_smmry <- chr15_attrib_df %>%
        group_by(Sex, Age.category) %>%
                # strength.fxm
        summarise(n = sum(!is.na(strength.fxm)),
                   avg.strength = format(round(mean(strength.fxm, 
                                                na.rm = T), 2), nsmall = 2),
                   stdev.strength = format(round(sd(strength.fxm, 
                                                na.rm = T), 2), nsmall = 2),
                   med.strength = format(round(median(strength.fxm,
                                                na.rm = T), 2), nsmall = 2),
                   # min.strength = format(round(min(strength.fxm,
                   #                            na.rm = T), 2), nsmall = 2),
                   # max.strength = format(round(max(strength.fxm,
                   #                            na.rm = T), 2), nsmall = 2),
                # degree.fxm
                   avg.degree = format(round(mean(degree.fxm, 
                                                na.rm = T), 2), nsmall = 2),
                   stdev.degree = format(round(sd(degree.fxm, 
                                                na.rm = T), 2), nsmall = 2),
                   med.degree = format(round(median(degree.fxm,
                                                na.rm = T), 2), nsmall = 2),
                   # min.degree= format(round(min(degree.fxm,
                   #                            na.rm = T), 2), nsmall = 2),
                   # max.degree = format(round(max(degree.fxm,
                   #                            na.rm = T), 2), nsmall = 2)
               ) %>%
        ungroup()
      
    ## j) tidy soc_intx_fxm_smmry
      soc_intx_fxm_smmry <- soc_intx_fxm_smmry %>%
        mutate(Intx = 'fxm') %>%
        relocate(Intx, .before = Sex)
  

  ### 4.2 Descriptive stats for female-female node level social interactions   
    ## a) Descriptive stats female-female social interactions
      soc_intx_fxf_smmry <- chr15_attrib_df %>%
        group_by(Age.category) %>%
                # strength.fxf
        summarise(n = sum(!is.na(strength.fxf)),
                  avg.strength = format(round(mean(strength.fxf, 
                                              na.rm = T), 2), nsmall = 2),
                  stdev.strength = format(round(sd(strength.fxf, 
                                              na.rm = T), 2), nsmall = 2),
                  med.strength = format(round(median(strength.fxf,
                                              na.rm = T), 2), nsmall = 2),
                  # min.strength = format(round(min(strength.fxf,
                  #                           na.rm = T), 2), nsmall = 2),
                  # max.strength= format(round(max(strength.fxf,
                  #                           na.rm = T), 2), nsmall = 2),
                # degree.fxf
                  avg.degree = format(round(mean(degree.fxf, 
                                              na.rm = T),2), nsmall = 2),
                  stdev.degree = format(round(sd(degree.fxf, 
                                              na.rm = T), 2), nsmall = 2),
                  med.degree = format(round(median(degree.fxf,
                                              na.rm = T), 2), nsmall = 2),
                  # min.degree = format(round(min(degree.fxf,
                  #                           na.rm = T), 2), nsmall = 2),
                  # max.degree. = format(round(max(degree.fxf,
                  #                           na.rm = T), 2), nsmall = 2)
        ) %>%
        ungroup()
      
    ## b) tidy soc_intx_fxf_smmry
      soc_intx_fxf_smmry <- soc_intx_fxf_smmry %>%
        mutate(Sex = 'f') %>%
        relocate(Sex, .before = Age.category) %>%
        mutate(Intx = 'fxf') %>%
        relocate(Intx, .before = Sex)
      
  
  ### 4.3 Descriptive stats for male-male node level social interactions                
    ## a) Descriptive stats male-male social interactions            
      soc_intx_mxm_smmry <- chr15_attrib_df %>%
                  group_by(Age.category) %>%
                # strength.mxm
                  summarise(n = sum(!is.na(strength.mxm)),
                  avg.strength = format(round(mean(strength.mxm, 
                                              na.rm = T), 2), nsmall = 2),
                  stdev.strength = format(round(sd(strength.mxm, 
                                              na.rm = T), 2), nsmall = 2),
                  med.strength = format(round(median(strength.mxm,
                                              na.rm = T), 2), nsmall = 2),
                  # min.strength = format(round(min(strength.mxm,
                  #                           na.rm = T), 2), nsmall = 2),
                  # max.strength = format(round(max(strength.mxm,
                  #                           na.rm = T), 2), nsmall = 2),
                # degree.mxm
                  avg.degree = format(round(mean(degree.mxm, 
                                              na.rm = T),2), nsmall = 2),
                  stdev.degree = format(round(sd(degree.mxm, 
                                              na.rm = T), 2), nsmall = 2),
                  med.degree = format(round(median(degree.mxm,
                                              na.rm = T), 2), nsmall = 2),
                  # min.degree = format(round(min(degree.mxm,
                  #                           na.rm = T), 2), nsmall = 2),
                  # max.degree = format(round(max(degree.mxm,
                  #                           na.rm = T), 2), nsmall = 2)
        ) %>%
        ungroup()

      
    ## b) tidy soc_intx_mxm_smmry
      ## h) tidy soc_intx_fxf_smmry
      soc_intx_mxm_smmry <- soc_intx_mxm_smmry %>%
        mutate(Sex = 'm') %>%
        relocate(Sex, .before = Age.category) %>%
        mutate(Intx = 'mxm') %>%
        relocate(Intx, .before = Sex)
      
  ### 3.4 Format soc intx tables for publication 
    ## a) row bind soc_intx_fxm_smmry, soc_intx_fxf_smmry and soc_intx_mxm_smmry
      soc_intx_smmry <- rbind(soc_intx_fxm_smmry, soc_intx_fxf_smmry, 
                                      soc_intx_mxm_smmry)
            
    ## b) transpose the data frame for easier viewing
      soc_intx_smmry <- as.data.frame(t(soc_intx_smmry))
      
    ## m) rename variables   
      soc_intx_smmry <- soc_intx_smmry %>%
        rename(c(`female-male.female.sy` = V1,
                 `female-male.female.asy` = V2,
                 `female-male.male.sy` = V3,
                 `female-male.male.asy` = V4,
                 `female-female.female.sy` = V5,
                 `female-female.female.asy` = V6,
                 `male-male.male.sy` = V7,
                 `male-male.male.asy` = V8)) %>%
        slice(c(-1, -2, -3)) %>%
        rownames_to_column(' ')
      
    ## n) format table creating flextable object 
      soc_intx_smmry_frmt <- nice_table(soc_intx_smmry, 
                                        separate.header = T)
      
    ## o) preview the table as a word doc    
      #print(soc_intx_smmry_frmt, preview = 'docx')
      
    ## p) save soc_intx_fxf_mxm_smmry_frmt as a word doc
      flextable::save_as_docx(soc_intx_smmry_frmt, 
                        path = here('output/soc_intx_smmry_frmt.docx'))  
      

      
###############################################################################
##############     5. Descriptive stats reproductive success     ##############
###############################################################################   
  
  ### 5.1 Descriptive stats female fecundity 
    ## a) Summarize female total fecundity and from clutches 2 and 3
     fecund_smmry <- chr15_attrib_df %>%
        ungroup() %>%
        group_by(Age.category) %>%
        filter(Sex == 'f') %>%
        summarise(
          # total fecundity
                  n.repro = round(sum(!is.na(total.fecundity)), 0),
                  avg.tot.fecund = format(round(mean(total.fecundity, 
                                           na.rm = T),2), nsmall = 2),
                  stdev.tot.fecund = format(round(sd(total.fecundity, 
                                           na.rm = T), 2), nsmall = 2),
                  # med.tot.fecund = format(round(median(total.fecundity,
                  #                            na.rm = T), 2), nsmall = 2),
                  # qt25.tot.fecund = format(round(quantile(total.fecundity, 
                  #                          0.25, na.rm = T), 2), nsmall = 2),
                  # qt75.tot.fecund = format(round(quantile(total.fecundity, 
                  #                          0.75, na.rm = T), 2), nsmall = 2),
                  # min.tot.fecund = format(round(min(total.fecundity,
                  #                          na.rm = T), 2), nsmall = 2),
                  # max.tot.fecund = format(round(max(total.fecundity,
                  #                          na.rm = T), 2), nsmall = 2)
          # nest 2 and 3 fecundity
                  n.nest.2.3.fecund = round(sum(!is.na(nest.2.3.fecund)), 0),
                  avg.nest.2.3.fecund = format(round(mean(nest.2.3.fecund, 
                                           na.rm = T),2), nsmall = 2),
                  stdev.nest.2.3.fecund = format(round(sd(nest.2.3.fecund, 
                                           na.rm = T), 2), nsmall = 2),
                  # med.nest.2.3.fecund = format(round(median(nest.2.3.fecund,
                  #                            na.rm = T), 2), nsmall = 2),
                  # qt25.nest.2.3.fecund = format(round(quantile(nest.2.3.fecund, 
                  #                          0.25, na.rm = T), 2), nsmall = 2),
                  # qt75.nest.2.3.fecund = format(round(quantile(nest.2.3.fecund, 
                  #                          0.75, na.rm = T), 2), nsmall = 2),
                  # min.nest.2.3.fecund = format(round(min(nest.2.3.fecund,
                  #                          na.rm = T), 2), nsmall = 2),
                  # max.nest.2.3.fecund = format(round(max(nest.2.3.fecund,
                  #                          na.rm = T), 2), nsmall = 2)
        ) %>%
        ungroup()
      
    # ## b) tidy fecund_smmry
    #   fecund_smmry <- fecund_smmry %>%
    #     mutate(Sex = 'f') %>%
    #     relocate(Sex, .before = Age.category) 
     
    # ## c) transpose the data frame for easier viewing
    #   fecund_smmry <- as.data.frame(t(fecund_smmry))
    #   
    # ## d) Reformat values to characters and fename variable
    #   fecund_smmry <- fecund_smmry %>%
    #     mutate(across(everything(), as.character)) %>%
    #     rename(c(`female.sy` = `V1`,
    #              `female.asy` = `V2`))
    #   
      
  ### 5.2 Descriptive stats male paternity from clutches 2 and 3
    ## a) Summarize male total paternity, epp, and spp from nests 2 and 3 
      patern_smmry <- chr15_attrib_df %>%
        ungroup() %>%
        group_by(Age.category) %>%
        filter(Sex == 'm') %>%
          # total paternity nests 2 and 3
        summarise(
                  #n.nest.2.3.tot.pat = round(sum(!is.na(nest.2.3.tot.pat)), 0),
                  avg.nest.2.3.tot.pat = format(round(mean(nest.2.3.tot.pat, 
                                              na.rm = T),2), nsmall = 2),
                  stdev.nest.2.3.tot.pat = format(round(sd(nest.2.3.tot.pat, 
                                              na.rm = T), 2), nsmall = 2),
                  # med.nest.2.3.tot.pat = format(round(median(nest.2.3.tot.pat,
                  #                               na.rm = T), 2), nsmall = 2),
                  # qt25.nest.2.3.tot.pat = format(round(quantile(nest.2.3.tot.pat, 
                  #                           0.25, na.rm = T), 2), nsmall = 2),
                  # qt75.nest.2.3.tot.pat = format(round(quantile(nest.2.3.tot.pat, 
                  #                           0.75, na.rm = T), 2), nsmall = 2),
                  # min.nest.2.3.tot.pat = format(round(min(nest.2.3.tot.pat,
                  #                           na.rm = T), 2), nsmall = 2),
                  # max.nest.2.3.tot.patt = format(round(max(nest.2.3.tot.pat,
                  #                           na.rm = T), 2), nsmall = 2),
          # epp nests 2 and 3
                  #n.nest.2.3.epp = round(sum(!is.na(nest.2.3.epp)), 0),
                  avg.nest.2.3.epp = format(round(mean(nest.2.3.epp, 
                                          na.rm = T), 2), nsmall = 2),
                  stdev.nest.2.3.epp = format(round(sd(nest.2.3.epp, 
                                          na.rm = T), 2), nsmall = 2),
                  # med.nest.2.3.epp = format(round(median(nest.2.3.epp,
                  #                           na.rm = T), 2), nsmall = 2),
                  # qt25.nest.2.3.epp = format(round(quantile(nest.2.3.epp, 
                  #                         0.25, na.rm = T), 2), nsmall = 2),
                  # qt75.nest.2.3.epp = format(round(quantile(nest.2.3.epp, 
                  #                         0.75, na.rm = T), 2), nsmall = 2),
                  # min.nest.2.3.epp = format(round(min(nest.2.3.epp,
                  #                         na.rm = T), 2), nsmall = 2),
                  # max.nest.2.3.epp = format(round(max(nest.2.3.epp,
                  #                         na.rm = T), 2), nsmall = 2),
          # spp nests 2 and 3
                  #n.nest.2.3.spp = round(sum(!is.na(nest.2.3.spp)), 0),
                  avg.nest.2.3.spp = format(round(mean(nest.2.3.spp, 
                                              na.rm = T), 2), nsmall = 2),
                  stdev.nest.2.3.spp = format(round(sd(nest.2.3.spp, 
                                              na.rm = T), 2), nsmall = 2),
                  # med.nest.2.3.spp = format(round(median(nest.2.3.spp,
                  #                               na.rm = T), 2), nsmall = 2)
                  # qt25.nest.2.3.spp = format(round(quantile(nest.2.3.spp, 
                  #                         0.25, na.rm = T), 2), nsmall = 2),
                  # qt75.nest.2.3.spp = format(round(quantile(nest.2.3.spp, 
                  #                         0.75, na.rm = T), 2), nsmall = 2),
                  # min.nest.2.3.spp = format(round(min(nest.2.3.spp,
                  #                         na.rm = T), 2), nsmall = 2),
                  # max.nest.2.3.spp = format(round(max(nest.2.3.spp,
                  #                         na.rm = T), 2), nsmall = 2),
        ) %>%
        ungroup()
    
    ## b) tidy patern_smmry
      # patern_smmry <- patern_smmry %>%
      #   mutate(Sex = 'm') %>%
      #   relocate(Sex, .before = Age.category) 

    # ## c) transpose the data frame for easier viewing
    #   patern_smmry <- as.data.frame(t(patern_smmry))
    # 
    # ## d) Reformat values to characters and rename variable
    #   patern_smmry <- patern_smmry %>%
    #     mutate(across(everything(), as.character)) %>%
    #     rename(c(`male.sy` = `V1`,
    #              `male.asy` = `V2`))
  
          
  ### 5.3 Format reproduction tables for publication 
    ## a) tidy fecund_smmry and patern_smmry
      fecund_smmry <- fecund_smmry %>%
        select(-c(Age.category))
      
      patern_smmry <- patern_smmry %>%
       select(-c(Age.category))
      
    ## b) column bind fecund_smmry and patern_smmry
     repro_smmry <- cbind(fecund_smmry, patern_smmry)
      
    ## c) transpose the data frame for easier viewing
     repro_smmry <- as.data.frame(t(repro_smmry))
     
    ## d) tidy repro_smmry
     repro_smmry <- repro_smmry %>%
       rownames_to_column('variable') %>%
       slice(c(-4)) %>%
       mutate(sex = ifelse(grepl('fecund', variable, fixed = T), 'f', 'm')) %>%
       relocate(sex, .before = V1) 
      
    ## e) rename variables   
     repro_smmry <- repro_smmry %>%
        rename(c(`reproduction.sy` = V1,
                 `reproduction.asy` = V2))
        
      
    ## f) format table creating flextable object 
     repro_smmry_frmt <- nice_table(repro_smmry, 
                                        separate.header = T)
      
    ## g) preview the table as a word doc    
      #print(repro_smmry_frmt, preview = 'docx')
      
    ## h) save repro_smmry_frmt as a word doc
      flextable::save_as_docx(repro_smmry_frmt, 
                              path = here('output/repro_smmry_frmt.docx'))  
      

      
###############################################################################
##############      6. Female-Male social interaction plots      ##############
###############################################################################      
    
  ### 6.1 Plots of DEGREE (female-male) for FEMALE nodes
    ## a) Density plot of CHR 2015 FEMALE degree
      fxm_f_deg_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'f') %>%
        #ggplot(aes(x = degree.fxm, fill = Year)) +
        ggplot(aes(x = degree.fxm, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.3 ) +
        xlim(0, 22) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('palegreen4', 'steelblue4'), 
                          name = 'Age.category',
                          labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. degree
        geom_vline(xintercept = as.numeric(soc_intx_smmry[5, 2]),
                   color = 'palegreen4', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[5, 2]) -0.5), 
                 y=0.02, label='sy avg. degree', angle=90) +
        # add vertical line for mean after second year fem. degree
        geom_vline(xintercept = as.numeric(soc_intx_smmry[5, 3]),
                   color = 'steelblue4', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[5, 3]) +0.4), 
                 y=0.02, label='asy avg. degree', angle=90) +
        labs(title = 'Density plot of female-male node level degree') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('female degree'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(fxm_f_deg_density_plot)
      
    ## c) Save plot
      ggsave('fxm_f_deg_density_plot.pdf', plot = fxm_f_deg_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)   
      
      
  ### 6.2 Plots of DEGREE (female-male) for MALE nodes
    ## a) Density plot of CHR 2015 MALE degree
      fxm_m_deg_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'm') %>%
        #ggplot(aes(x = degree.fxm, fill = Year)) +
        ggplot(aes(x = degree.fxm, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.25 ) +
        xlim(0, 22) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('palegreen3', 'steelblue3'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. degree
        geom_vline(xintercept = as.numeric(soc_intx_smmry[5, 4]),
                   color = 'palegreen3', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[5, 4]) +0.2), 
                 y=0.025, label='sy avg. degree', angle=90) +
        # add vertical line for mean after second year fem. degree
        geom_vline(xintercept = as.numeric(soc_intx_smmry[5, 5]),
                   color = 'steelblue3', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[5, 5]) -0.3), 
                 y=0.025, label='asy avg. degree', angle=90) +
        labs(title = 'Density plot of female-male node level degree') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('male degree'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(fxm_m_deg_density_plot)
      
    ## c) Save plot
      ggsave('fxm_m_deg_density_plott.pdf', plot = fxm_m_deg_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)   
      
     
  ### 6.3 Plots of STRENGTH (female-male) for FEMALE nodes
    ## a) Density plot of CHR 2015 FEMALE strength
      fxm_f_strength_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'f') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = strength.fxm, color = Age.category, 
                   label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.2 ) +
        xlim(0, 101) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('brown4', 'purple4'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. strength
        geom_vline(xintercept = as.numeric(soc_intx_smmry[2, 2]),
                   color = 'brown4', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[2, 2]) -1.9), 
                 y=0.005, label='sy avg. strength', angle=90) +
        # add vertical line for mean after second year fem. strength
        geom_vline(xintercept = as.numeric(soc_intx_smmry[2, 3]),
                   color = 'purple4', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[2, 3]) +1.5), 
                 y=0.005, label='asy avg. strength', angle=90) +
        labs(title = 'Density plot of female-male node level strength') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('female strength'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(fxm_f_strength_density_plot)
      
    ## c) Save plot
      ggsave('fxm_f_strength_density_plot.pdf', 
             plot = fxm_f_strength_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
  
  ### 6.4 Plots of STRENGTH (female-male) for MALE nodes
    ## a) Density plot of CHR 2015 MALE strength
      fxm_m_strength_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'm') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = strength.fxm, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.2 ) +
        xlim(0, 101) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('brown3', 'purple3'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. strength
        geom_vline(xintercept = as.numeric(soc_intx_smmry[2, 4]),
                   color = 'brown3', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[2, 4]) -1.9), 
                 y=0.005, label='sy avg. strength', angle=90) +
        # add vertical line for mean after second year fem. strength
        geom_vline(xintercept = as.numeric(soc_intx_fxm_smmry[2, 5]),
                   color = 'purple3', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_fxm_smmry[2, 5]) +1.6), 
                 y=0.005, label='asy avg. strength', angle=90) +
        labs(title = 'Density plot of female-male node level strength') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('male strength'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(fxm_m_strength_density_plot)
      
    ## c) Save plot
      ggsave('fxm_m_strength_density_plot.pdf', 
             plot = fxm_m_strength_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
      
      
###############################################################################
##############     7. Female-Female social interaction plots     ##############
###############################################################################      
      
  ### 7.1 Plots of DEGREE (female-female) for FEMALE nodes
    ## a) Density plot of CHR 2015 FEMALE degree
      fxf_f_deg_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'f') %>%
        #ggplot(aes(x = degree.fxf, fill = Year)) +
        ggplot(aes(x = degree.fxf, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.3 ) +
        xlim(0, 22) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('blue4', 'bisque4'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. degree
        geom_vline(xintercept = as.numeric(soc_intx_smmry[5, 6]),
                   color = 'blue4', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[5, 6]) -0.5), 
                 y=0.02, label='sy avg. degree', angle=90) +
        # add vertical line for mean after second year fem. degree
        geom_vline(xintercept = as.numeric(soc_intx_smmry[5, 7]),
                   color = 'bisque4', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[5, 7]) -0.5), 
                 y=0.02, label='asy avg. degree', angle=90) +
        labs(title = 'Density plot of female-female node level degree') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('female degree'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(fxf_f_deg_density_plot)
      
    ## c) Save plot
      ggsave('fxf_f_deg_density_plot.pdf', 
             plot = fxf_f_deg_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)   
      
      
  ### 7.2 Plots of STRENGTH (female-female) for FEMALE nodes
    ## a) Density plot of CHR 2015 FEMALE strength
      fxf_f_strength_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'f') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = strength.fxf, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.15 ) +
        xlim(0, 101) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('blue3', 'bisque3'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. strength
        geom_vline(xintercept = as.numeric(soc_intx_smmry[2, 6]),
                   color = 'blue3', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[2, 6]) -2), 
                 y=0.005, label='sy avg. strength', angle=90) +
        # add vertical line for mean after second year fem. strength
        geom_vline(xintercept = as.numeric(soc_intx_smmry[2, 7]),
                   color = 'bisque3', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[2, 7]) -2.2), 
                 y=0.005, label='asy avg. strength', angle=90) +
        labs(title = 'Density plot of female-female node level strength') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('female strength'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(fxf_f_strength_density_plot)
      
    ## c) Save plot
      ggsave('fxf_f_strength_density_plot.pdf', 
             plot = fxf_f_strength_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)   
      
      
      
###############################################################################
##############       8. Male-Male social interaction plots       ##############
###############################################################################      
      
  ### 8.1 Plots of DEGREE (male-male) for MALE nodes
    ## a) Density plot of CHR 2015 MALE degree
      mxm_m_deg_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'm') %>%
        #ggplot(aes(x = degree.mxm, fill = Year)) +
        ggplot(aes(x = degree.mxm, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.2 ) +
        xlim(0, 22) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('darkslategray4', 'tomato4'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. degree
        geom_vline(xintercept = as.numeric(soc_intx_smmry[5, 8]),
                   color = 'darkslategray4', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[5, 8]) -0.5), 
                 y=0.02, label='sy avg. degree', angle=90) +
        # add vertical line for mean after second year fem. degree
        geom_vline(xintercept = as.numeric(soc_intx_smmry[5, 9]),
                   color = 'tomato4', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[5, 9]) -0.5), 
                 y=0.02, label='asy avg. degree', angle=90) +
        labs(title = 'Density plot of male-male node level degree') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('male degree'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(mxm_m_deg_density_plot)
      
    ## c) Save plot
      ggsave('mxm_m_deg_density_plot.pdf', plot = mxm_m_deg_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)   
      
      
  ### 8.2 Plots of STRENGTH (male-male) for MALE nodes
    ## a) Density plot of CHR 2015 MALE strength
      mxm_m_strength_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'm') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = strength.mxm, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.17 ) +
        xlim(0, 101) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('darkslategray3', 'tomato3'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. strength
        geom_vline(xintercept = as.numeric(soc_intx_smmry[2, 8]),
                   color = 'darkslategray3', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[2, 8]) +1.5), 
                 y=0.005, label='sy avg. strength', angle=90) +
        # add vertical line for mean after second year fem. strength
        geom_vline(xintercept = as.numeric(soc_intx_smmry[2, 9]),
                   color = 'tomato3', size = 2 ) +
        annotate('text', x = (as.numeric(soc_intx_smmry[2, 9]) -2), 
                 y=0.005, label='asy avg. strength', angle=90) +
        labs(title = 'Density plot of male-male node level strength') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('male strength'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(mxm_m_strength_density_plot)
      
    ## c) Save plot
      ggsave('mxm_m_strength_density_plot.pdf', 
             plot = mxm_m_strength_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
   

###############################################################################
##############             9. Female fecundity plots             ##############
###############################################################################           

  ### 9.1 Plots of female total fecundity from nests 1, 2, and 3
    ## a) Frequency distribution of female total fecundity
      f_tot_fecund_freq_dist_plot <- 
        # raw data plot
        ggplot(data = subset(chr15_attrib_df
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
      print(f_tot_fecund_freq_dist_plot)
      
      ## c) Save Plot
      # use ggsave to save the plot
      ggsave('f_tot_fecund_freq_dist_plot.pdf', 
             plot = f_tot_fecund_freq_dist_plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 7,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
    
  ### 9.2 Plots of female total fecundity from nests 2 and 3  
    ## a) Frequency distribution of female fecundity from nests 2 and 3
      f_nest_2_3_fecund_freq_dist_plot <- 
        # raw data plot
        ggplot(data = subset(chr15_attrib_df
                             , Sex == 'f'), aes(x = nest.2.3.fecund)) +
        geom_histogram(binwidth = 1, color = 'black', fill = 'steelblue4') +
        # scale the density plot line
        geom_density(aes(y=after_stat(density)/0.025), color = 'red', 
                     linewidth = 3) +
        # Titles, axes, and legends
        labs(title = 'Frequency distribution of female fecundity, nests 2 and 3') +
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
        xlab(expression(bold('female fecundity (from nest 2 and 3)'))) +
        ylab(expression(bold('count of individuals'))) 
      
    ## b) View plot
      print(f_nest_2_3_fecund_freq_dist_plot)
      
    ## c) Save Plot
      # use ggsave to save the plot
      ggsave('f_nest_2_3_fecund_freq_dist_plot.pdf', 
             plot = f_nest_2_3_fecund_freq_dist_plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 7,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
      
  ### 9.3 Bivariate PLOTS for female total fecundity data by age      
    ## a) Density plot of CHR 2015 FEMALE fecundity for nests 1, 2, and 3
      f_tot_fecund_age_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'f') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = total.fecundity, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.17 ) +
        xlim(0, 20) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('palegreen4', 'steelblue4'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. total fecundity
        geom_vline(xintercept = as.numeric(repro_smmry[2, 3]),
                   color = 'palegreen4', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[2, 3]) + 0.5), 
                 y=0.017, label='sy avg. fecundity', angle=90) +
        # add vertical line for mean after second year fem. total fecundity
        geom_vline(xintercept = as.numeric(repro_smmry[2, 4]),
                   color = 'steelblue4', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[2, 4]) + 0.5), 
                 y=0.017, label='asy avg. fecundity', angle=90) +
        labs(title = 'Density plot of female total fecundity') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('total fecundity'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(f_tot_fecund_age_density_plot)
      
    ## c) Save plot
      ggsave('f_tot_fecund_age_density_plot.pdf', 
             plot = f_tot_fecund_age_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
      
      
    ### 9.4 Bivariate PLOTS for female nest 2-3 fecundity data by age      
      ## a) Density plot of CHR 2015 FEMALE fecundity for nests 2 and 3
      f_nest_2_3_fecund_age_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'f') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = nest.2.3.fecund, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.17 ) +
        xlim(0, 20) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('palegreen4', 'steelblue4'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year fem. nest 2-3 fecundity
        geom_vline(xintercept = as.numeric(repro_smmry[4, 3]),
                   color = 'palegreen4', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[4, 3]) + 0.45), 
                 y=0.025, label='sy avg. fecundity', angle=90) +
        # add vertical line for mean after second year fem. nest 2-3 fecundity
        geom_vline(xintercept = as.numeric(repro_smmry[4, 4]),
                   color = 'steelblue4', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[4, 4]) + 0.45), 
                 y=0.025, label='asy avg. fecundity', angle=90) +
        labs(title = 'Density plot of female fecundity from nest 2 and 3') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('nest 2-3 fecundity'))) +
        ylab(expression(bold('density')))
      
      ## b) Print plot 
      print(f_nest_2_3_fecund_age_density_plot)
      
      ## c) Save plot
      ggsave('f_nest_2_3_fecund_age_density_plot.pdf', 
             plot = f_nest_2_3_fecund_age_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
      
      
      
###############################################################################
##############             10. Male paternity plots              ##############
############################################################################### 
      
  ### 10.1 Plots of male total paternity from nests 2 and 3  
    ## a) Frequency distribution of male total paternity (nests 2 & 3)
      m_nest_2_3_tot_patern_freq_dist_plot <- 
        # raw data plot
        ggplot(data = subset(chr15_attrib_df
                             , Sex == 'm'), aes(x = nest.2.3.tot.pat)) +
        geom_histogram(binwidth = 1, color = 'black', fill = 'wheat3') +
        # scale the density plot line
        geom_density(aes(y=after_stat(density)/0.025), color = 'red', 
                     linewidth = 3) +
        # Titles, axes, and legends
        labs(title = 'Frequency distribution of male total paternity (nests 2-3)') +
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
        xlab(expression(bold('male total paternity (from nests 2 and 3)'))) +
        ylab(expression(bold('count of individuals'))) 
      
    ## b) View plot
      print(m_nest_2_3_tot_patern_freq_dist_plot)
      
    ## c) Save Plot
      # use ggsave to save the plot
      ggsave('m_nest_2_3_tot_patern_freq_dist_plot.pdf', 
             plot = m_nest_2_3_tot_patern_freq_dist_plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 7,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
    
  ### 10.2 Plots of male extra pair paternity (epp) from nests 2 and 3  
    ## a) Frequency distribution of male extra pair paternity (nests 2 & 3)
      m_nest_2_3_epp_freq_dist_plot <- 
        # raw data plot
        ggplot(data = subset(chr15_attrib_df
                             , Sex == 'm'), aes(x = nest.2.3.epp)) +
        geom_histogram(binwidth = 1, color = 'black', fill = 'wheat4') +
        # scale the density plot line
        geom_density(aes(y=after_stat(density)/0.025), color = 'red', 
                     linewidth = 3) +
        # Titles, axes, and legends
        labs(title = 'Frequency distribution of male extra pair paternity (nests 2-3)') +
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
        xlab(expression(bold('male extra pair paternity (epp; from nests 2 and 3)'))) +
        ylab(expression(bold('count of individuals'))) 
      
    ## h) View plot
      print(m_nest_2_3_epp_freq_dist_plot)
      
    ## i) Save Plot
      # use ggsave to save the plot
      ggsave('m_nest_2_3_epp_freq_dist_plot.pdf', 
             plot = m_nest_2_3_epp_freq_dist_plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 7,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
      
  ### 10.3 Plots of male social pair paternity (spp) from nests 2 and 3  
    ## a) Frequency distribution of male socia pair paternity (nests 2 & 3)
      m_nest_2_3_spp_freq_dist_plot <- 
        # raw data plot
        ggplot(data = subset(chr15_attrib_df
                             , Sex == 'm'), aes(x = nest.2.3.spp)) +
        geom_histogram(binwidth = 1, color = 'black', fill = 'wheat2') +
        # scale the density plot line
        geom_density(aes(y=after_stat(density)/0.025), color = 'red', 
                     linewidth = 3) +
        # Titles, axes, and legends
        labs(title = 'Frequency distribution of male social pair paternity (nests 2-3)') +
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
        xlab(expression(bold('male social pair paternity (spp; from nests 2 and 3)'))) +
        ylab(expression(bold('count of individuals'))) 
      
    ## k) View plot
      print(m_nest_2_3_spp_freq_dist_plot)
      
    ## l) Save Plot
      # use ggsave to save the plot
      ggsave('m_nest_2_3_spp_freq_dist_plott.pdf', 
             plot = m_nest_2_3_spp_freq_dist_plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 7,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
      
  ### 10.4 Bivariate PLOTS for male nest 2-3 total paternity data by age      
    ## a) Density plot of CHR 2015 MALE total paternity for nests 2 and 3
      m_nest_2_3_tot_patern_age_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'm') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = nest.2.3.tot.pat, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.07 ) +
        xlim(0, 20) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('lightgoldenrod3', 'indianred3'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year male paternity
        geom_vline(xintercept = as.numeric(repro_smmry[6, 3]),
                   color = 'lightgoldenrod3', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[6, 3]) + 0.6), 
                 y=0.025, label='sy avg. tot. paternity', angle=90) +
        # add vertical line for mean after second year male paternity
        geom_vline(xintercept = as.numeric(repro_smmry[6, 4]),
                   color = 'indianred3', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[6, 4]) + 0.5), 
                 y=0.025, label='asy avg. tot. paternity', angle=90) +
        labs(title = 'Density plot of male total paternity (nest 2-3)') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('nest 2-3 total paternity'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(m_nest_2_3_tot_patern_age_density_plot)
      
    ## c) Save plot
      ggsave('m_nest_2_3_tot_patern_age_density_plot.pdf', 
             plot = m_nest_2_3_tot_patern_age_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
      
  ### 10.5 Bivariate PLOTS for male nest 2-3 extra pair paternity data by age      
    ## a) Density plot of CHR 2015 MALE extra pair paternity for nests 2 and 3
      m_nest_2_3_epp_age_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'm') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = nest.2.3.epp, color = Age.category, label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.7 ) +
        xlim(0, 20) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('lightgoldenrod4', 'indianred4'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year male epp
        geom_vline(xintercept = as.numeric(repro_smmry[8, 3]),
                   color = 'lightgoldenrod4', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[8, 3]) - 0.6), 
                 y=0.02, label='sy avg. epp', angle=90) +
        # add vertical line for mean after second year male paternity
        geom_vline(xintercept = as.numeric(repro_smmry[8, 4]),
                   color = 'indianred4', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[8, 4]) - 0.6), 
                 y=0.02, label='asy avg. epp', angle=90) +
        labs(title = 'Density plot of male node level extra pair paternity (nest 2-3)') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('nest 2-3 extra pair paternity (epp)'))) +
        ylab(expression(bold('density')))
      
    ## h) Print plot 
      print(m_nest_2_3_epp_age_density_plot)
      
    ## i) Save plot
      ggsave('m_nest_2_3_epp_age_density_plot.pdf', 
             plot = m_nest_2_3_epp_age_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
    
      
  ### 10.6 Bivariate PLOTS for male nest 2-3 social pair paternity data by age      
    ## a) Density plot of CHR 2015 MALE social pair paternity for nests 2 and 3
      m_nest_2_3_spp_age_density_plot <- chr15_attrib_df %>%
        filter(Sex == 'm') %>%
        #ggplot(aes(x = strength.fxm, fill = Year)) +
        ggplot(aes(x = nest.2.3.spp, color = Age.category, 
                   label = Age.category)) +
        geom_textdensity(linewidth = 4, size = 8, hjust = 0.5 ) +
        xlim(0, 20) +
        # geom_density()
        # geom_histogram(color='gray50', alpha=0.6, position = 'identity', 
        #                binwidth = 1) +
        scale_color_manual(values=c('lightgoldenrod2', 'indianred2'), 
                           name = 'Age.category',
                           labels = c('second year', 'after second year')) +
        # add vertical line for mean second year male spp
        geom_vline(xintercept = as.numeric(repro_smmry[10, 3]),
                   color = 'lightgoldenrod2', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[10, 3]) - 0.4), 
                 y=0.02, label='sy avg. spp', angle=90) +
        # add vertical line for mean after second year male paternity
        geom_vline(xintercept = as.numeric(repro_smmry[10, 4]),
                   color = 'indianred2', size = 2 ) +
        annotate('text', x = (as.numeric(repro_smmry[10, 4]) + 0.4), 
                 y=0.02, label='asy avg. spp', angle=90) +
        labs(title = 'Density plot of male node level social pair paternity (nests 2-3)') +
        theme(plot.title = element_text(hjust = 0.5, size = 18)) + # center title
        #theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size for axes titles
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
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
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 0)),
              legend.title=element_blank(),
              legend.text=element_text(size=18),
              #legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(bold('nest 2-3 social pair paternity (spp)'))) +
        ylab(expression(bold('density')))
      
    ## b) Print plot 
      print(m_nest_2_3_spp_age_density_plot)
      
    ## c) Save plot
      ggsave('m_nest_2_3_spp_age_density_plot.pdf', 
             plot = m_nest_2_3_spp_age_density_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 6, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
      
      