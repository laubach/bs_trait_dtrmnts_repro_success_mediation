################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############               1. Tidy attribute data                 #############
#############                                                      #############
#############                  By: Zach Laubach                    #############
#############                created: 14 Aug 2024                  #############
#############              last updated: 24 Feb 2024               #############
################################################################################


  ### PURPOSE: Tidy the attribute data for male and female barn swallows
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Tidy attribute data frame
    # 4. Export data
  


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
  
    ## a) Load attribute data CHR 2015
      load(here('data/1_chr_attrib_data.RData'))
      


###############################################################################
##############          3. Tidy attribute data frame             ##############
############################################################################### 
  
  ### 3.1 Basic Tidying
    ## a) Format data to all lower case
      chr15_attrib_df <- AllCharactersToLower(chr15_attrib_df)
      
    ## b) create a composite breast brightness variable that includes post
      # treatment breast brightness for manipulated males, and original
      # breast brightness for unmanipulated males.
      chr15_attrib_df <- chr15_attrib_df  %>%
        #rowwise()%>%
        mutate(R.bright.treat.and.orig = ifelse(is.na(Post.R_avg.bright),
                                          R_avg.bright, Post.R_avg.bright))
      
    ## c) clean up df 
      chr15_attrib_df <- chr15_attrib_df  %>%
        select(!c(Mate.band:Mate.tag))
      
    # ## d) clean up df because not using the physiology data
    #   chr15_attrib_df <- chr15_attrib_df  %>%
    #     select(!c(PreCort:PostT))
    ## d) clean up df because not using all the physiology data
      chr15_attrib_df <-  chr15_attrib_df  %>%
        select(!c(Left.wing.mites:Post8OHdG))
      
    ## e) clean up df to remove old reproductive success data
      chr15_attrib_df <- chr15_attrib_df  %>%
        select(!c(Clutch.1.Eggs:Total.offspring..minus.clutch.3.))
      
    ## f) Re-code 'Sex' as *nominal* factor (with ordered levels)
      # Set levels (odering) of variable 
      chr15_attrib_df <- transform(chr15_attrib_df, 
                                   Sex = factor(Sex,
                                                levels = c('f', 'm')))
      
    ## g) Re-code 'Age.category' as *nominal* factor (with ordered levels)
      # Set levels (odering) of variable 
      chr15_attrib_df <- transform(chr15_attrib_df, 
                                   Age.category = factor(Age.category,
                                                    levels = c('sy', 'asy')))
      
    ## h) create a difference in mass (Pre.mass - Post.mass) variable
      # lower values = poorer condition
      chr15_attrib_df <- chr15_attrib_df  %>%
        #rowwise()%>%
        mutate(Diff.mass = Pre.mass - Post.mass)
      
    ## i) create a difference in testosterone (PreT - PostT) variable
      # lower values = poorer condition
      chr15_attrib_df <- chr15_attrib_df  %>%
        #rowwise()%>%
        mutate(Diff.T = PreT - PostT)
      
    ## i) create a difference in Cort (PostCort - PreCort) variable
      # higher values = poorer condition
      chr15_attrib_df <- chr15_attrib_df  %>%
        #rowwise()%>%
        mutate(Diff.Cort = PostCort - PreCort)
    
      

###############################################################################
##############                  4. Export data                   ##############
###############################################################################
      
  ### 4.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 pre-manipulation >= 20 RSSI
      save(file = here('data/2_chr_attrib_df.RData'), 
           list = c('chr15_attrib_df'))
      
      

    
