################################################################################
#############        The role of age and plumage traits as         #############  
#############       determinants of reproductive success and       #############
#############       the mediating role of social interactions      #############
#############                                                      #############
#############          3. Tidy reproductive success data           #############
#############                                                      #############
#############                  By: Zach Laubach                    #############
#############                created: 20 Aug 2024                  #############
#############              last updated: 3 Dec 2024                #############
################################################################################


  ### PURPOSE: Tidy the reproductive success and data for male and 
             # female barn swallows
  
  # Code Blocks
    # 1: Configure work space
    # 2: Import data 
    # 3: Tidy reproductive success data 
    # 4: Check and tidy missing data 
    # 5: Export data files
  


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
      load(here('data/3_chr15_repro_data.RData'))
    
    ## b) Load tidy attribute data
      load(here('data/3_chr15_attribute_df.RData'))
      


###############################################################################
##############        3. Tidy reproductive success data          ##############
############################################################################### 
 
  ### 3.1 Tidy chr15_repro_suc and join data to chr15_attrib_df
    ## a) Format data to all lower case
      chr15_repro_suc <- AllCharactersToLower(chr15_repro_suc)

    ## b)  Calculate female total fecundity
      chr15_repro_suc <- chr15_repro_suc  %>%
        rowwise()%>%
        mutate(total.fecundity = case_when(Sex == 'f'
                                           ~ sum(nest.1.fecundity, 
                                                 nest.2.fecundity, 
                                            nest.3.fecundity)))       
  
    ## c)  Calculate male total paternity
      chr15_repro_suc <- chr15_repro_suc  %>%
        rowwise()%>%
        mutate(total.paternity = case_when(Sex == 'm'
                                           ~ sum(attmpt.1.tot.pat, 
                                                 attmpt.2.tot.pat, 
                                                 attmpt.3.tot.pat)))    

    ## d) Add measures of reproductive success chr15_attrib_df
      chr15_attrib_df <- chr15_attrib_df  %>%
        left_join(select(chr15_repro_suc, c(Band.ID, soc.mate.ID.nest.1,
                                            soc.mate.tag.nest.1,
                                            nest.1, total.fecundity,
                                            attmpt.1.tot.pat, attmpt.2.tot.pat,
                                            total.paternity)),
                  by = c('Band.ID' = 'Band.ID'), 
                  copy = F) 
      
      # NOTEs on reproductive success    
      # Female 
        # total.fecundity = nest.1.fecundity + nest.2.fecundity + nest.3.fecundity
      # Male 
        # orig. clutch paternity = attmpt.1.tot.pat
        # replace clutch paternity = attmpt.2.tot.pat
        #total.paternity = attmpt.1.tot.pat + attmpt.2.tot.pat + attmpt.3.tot.pat
      
      
      
###############################################################################
##############         4.  Check and tidy missing data           ##############
###############################################################################      
      
  ### 4.1 Manual data checks and cleaning
    
#*****************************************************************************#
#************************* Manual data cleaning ******************************#
#*
      # NOTE: Tag ID 57 does not pair until late Aug...fundamentally different
      # than other birds, so remove it here.                
    ## a) Check/remove bird tag 57 from chr15_attrib_pw and chr15_attrib_df
      chr15_attrib_df <- chr15_attrib_df %>%
        filter(!(Tag == 57))   
      
    ## b) Check/remove bird tag remove males who are untagged from 
      # chr15_attrib_df
      chr15_attrib_df <- chr15_attrib_df %>%
        filter(!(Band.ID == '2640-97157' | Band.ID == '2640-97158'))   
      
#*****************************************************************************#
#*****************************************************************************#
  
      
  ### 4.2 Create attribute data frames that match the bird IDs in the 
      # social networks df  
  # Remove birds with insufficient premanip network data or just note and exclude in
      # downstream models 
#*****************************************************************************#
#************************* Manual data cleaning ******************************#
      # NOTE: Remove birds with insufficient premanip network data based on
      # checking the earliest vs latest encounter net to ensure networks
      # are based on when all Tags working
      
    ## a) remove bird 49 from both pre and post manip data frames to mactch
          # social networks 
      # Tag worked from 6/13-6/14
      chr15_attrib_pre_df <- chr15_attrib_df %>%
        filter(!(Tag == 49 )) 
      
      
      # NOTE: Remove birds with insufficient post manip network data based on
      # checking the earliest vs latest encounter net to ensure networks
      # are based on when all Tags working
      
    ## b) Remove birds 67, 40, 61, 51, 52, and 55 from post manip networks 
      # Tag stopped on 6/19
      chr15_attrib_post_df <- chr15_attrib_pre_df %>%
        filter(!(Tag == 67)) %>%
        filter(!(Tag == 40)) %>%
        filter(!(Tag == 61)) %>%
        filter(!(Tag == 51)) %>%
        filter(!(Tag == 52)) %>%
        filter(!(Tag == 55)) %>%
        # and birds 56, 80, 82, 119 which have no Tag data
        filter(!(Tag == 56)) %>%
        filter(!(Tag == 80)) %>%
        filter(!(Tag == 82)) %>%
        filter(!(Tag == 119))
        
      
#*****************************************************************************#
#*****************************************************************************#  

      
      
###############################################################################
##############                5. Export data files               ##############
###############################################################################
      
  ### 5.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 pre-manipulation repro. success
      save(file = here('data/4_chr15_attrib_data.RData'), 
           list = c('chr15_attrib_df', 'chr15_attrib_pre_df', 
                    'chr15_attrib_post_df'))
      
      


