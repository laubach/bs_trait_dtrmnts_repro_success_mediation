################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############          2. Tidy reproductive success data           #############
#############                                                      #############
#############                  By: Zach Laubach                    #############
#############                created: 20 Aug 2024                  #############
#############             last updated: 24 Feb 2024                #############
################################################################################


  ### PURPOSE: Tidy the reproductive success and data for male and 
             # female barn swallows
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Tidy reproductive success data 
    # 4. Check and tidy missing data 
    # 5. Export data
  


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
      
    # load naniar package
      library ('naniar')
    
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
    ## a) Load reproductive success and nest distance data for CHR 2015
      load(here('data/2_chr_repro_data.RData'))
    
    ## b) Load tidy attribute data
      load(here('data/2_chr_attrib_df.RData'))
      


###############################################################################
##############        3. Tidy reproductive success data          ##############
############################################################################### 
 
  ### 3.1 Tidy chr15_repro_suc and join data to chr15_attrib_df
    ## a) Format data to all lower case
      chr15_repro_suc <- AllCharactersToLower(chr15_repro_suc)
      
    ## b) Replace 'incomplete' with NA
      chr15_repro_suc <- chr15_repro_suc %>% 
        replace_with_na_all(condition = ~.x == 'incomplete')
      
    ## c) Format variable types from character to numeric
      chr15_repro_suc[, 16:21] <- 
        sapply(chr15_repro_suc[, 16:21], as.numeric)
      
      chr15_repro_suc[, 26:32] <- 
        sapply(chr15_repro_suc[, 26:32], as.numeric)
      
      chr15_repro_suc[, 37:43] <- 
        sapply(chr15_repro_suc[, 37:43], as.numeric)

    ## d)  Calculate female total fecundity
      chr15_repro_suc <- chr15_repro_suc  %>%
        rowwise()%>%
        mutate(total.fecundity = case_when(Sex == 'f'
                                           ~ sum(nest.1.fecundity,
                                                 nest.2.fecundity,
                                            nest.3.fecundity)))
      
    ## e) Calculate female fecundity from nests 2 and 3
      chr15_repro_suc <- chr15_repro_suc  %>%
        rowwise()%>%
        mutate(nest.2.3.fecund = case_when(Sex == 'f'
                                           ~ sum(nest.2.fecundity, 
                                                 nest.3.fecundity))) 
  
      
    ## f) Calculate male total paternity from nest 2 and 3
      chr15_repro_suc <- chr15_repro_suc  %>%
        rowwise()%>%
        mutate(nest.2.3.tot.pat = case_when(Sex == 'm'
                                            ~ sum(attmpt.2.tot.pat, 
                                                  attmpt.3.tot.pat))) 
      # NOTE: updated on 6/12/25...exclude egg paternity from initial clutch
      # due to missing data

    ## g) Calculate male extra pair paternity from nest 2 and 3
      chr15_repro_suc <- chr15_repro_suc  %>%
        rowwise()%>%
        mutate(nest.2.3.epp = case_when(Sex == 'm' 
                                        ~ sum(
                                          attmpt.2.epp, 
                                          attmpt.3.epp))) 

    ## h) Calculate male total within nest paternity
      # create a vector of male epp for nest 2 and 3
      chr15_repro_suc <- chr15_repro_suc  %>%
        rowwise()%>%
        mutate(nest.2.3.spp = nest.2.3.tot.pat - nest.2.3.epp)
      
      
    ## i) Add measures of reproductive success chr15_attrib_df
      chr15_attrib_df <- chr15_attrib_df  %>%
        left_join(select(chr15_repro_suc, c(Band.ID, soc.mate.ID.nest.1,
                                            soc.mate.tag.nest.1,
                                            nest.1, 
                                            total.fecundity,
                                            nest.2.3.fecund,
                                            #attmpt.1.tot.pat, attmpt.2.tot.pat,
                                            #total.paternity
                                            nest.2.3.tot.pat,
                                            nest.2.3.epp, nest.2.3.spp)),
                  by = c('Band.ID' = 'Band.ID'), 
                  copy = F) 
      

      
      
      
###############################################################################
##############         4.  Check and tidy missing data           ##############
###############################################################################      
      
  ### 4.1 Create attribute data frames that match the bird IDs in the 
      # social networks df. Remove birds with insufficient pre-manip network data 
      
#*****************************************************************************#
#************************* Manual data cleaning ******************************#
      # NOTE: Remove birds with insufficient pre-manip network data based on
      # checking the earliest vs latest encounter net to ensure networks
      # are based on when all Tags working
      
    ## a) remove bird 49 from both pre and post manip data frames to 
          # match social networks because Tag worked from 6/13-6/14; Tag 108 
          # for female 2640-97169 and tag 69 for male 2591-02869 
          # never worked, so not included in 
          # either pre or post network or attribute file
      chr15_attrib_pre_df <- chr15_attrib_df %>%
        filter(!(Tag == 49 | Tag == 69 | Tag == 108)) 
  

      
      # NOTE: Remove birds with insufficient post manip network data based on
      # checking the earliest vs latest encounter net to ensure networks
      # are based on when all Tags working
      
    ## b) Remove males 55, 61, 67 and females 40, 51, 52 and from post manip 
      # networks. Tag stopped on 6/19
      chr15_attrib_post_df <- chr15_attrib_pre_df %>%
        filter(!(Tag == 55)) %>% # males
        #filter(!(Tag == 57)) %>% already removed for incomplete fecundity
        filter(!(Tag == 61)) %>%
        filter(!(Tag == 67)) %>%
        filter(!(Tag == 40)) %>% # females
        filter(!(Tag == 51)) %>%
        filter(!(Tag == 52)) %>%
        #filter(!(Tag == 116)) %>% already removed for incomplete paternity
        # and birds 56, 80, 82 which have no Tag data by 6/19
        filter(!(Tag == 56)) %>% # male
        filter(!(Tag == 82)) %>%
        filter(!(Tag == 80)) # female
         
#*****************************************************************************#
#*****************************************************************************#  

      
      
###############################################################################
##############                  5. Export data                   ##############
###############################################################################
      
  ### 5.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 pre-manipulation repro. success
      save(file = here('data/3_chr_attrib_data.RData'), 
           list = c('chr15_attrib_pre_df', 
                    'chr15_attrib_post_df'))
      
      


