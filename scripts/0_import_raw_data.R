################################################################################
#############        The role of age and plumage traits as         #############  
#############       determinants of reproductive success and       #############
#############       the mediating role of social interactions      #############
#############                                                      #############
#############                    0. Data Import                    #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 6 Aug 2024                  #############
#############               last updated: 3 Dec 2024               #############
################################################################################


  ### PURPOSE: Load barn swallow trait, social network, and reproductive success
              # data from raw data files
  
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import social intx. data 
    # 3. Import attribute data
    # 4. Import repro. success data
    # 5. Export data files
  



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
 
      # load here packages
        library ('here')

        
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 4.4.2 (2024-10-31)
    # Platform: x86_64-apple-darwin20
    # Running under: macOS Sequoia 15.1.1
    
  
  ### 1.4 Set working directory 
    setwd(here())
    
    
      
###############################################################################
##############            2. Import social intx. data            ##############
###############################################################################    
      
  ### 2.1 Import sample data files
    ## a) Path to subdirectory 
      chr15_intx_path <- here('raw_data/2015_CHR_soc/')
      # a test data set
      #chr15_intx_path_test <- here('raw_data/2015_CHR_soc/test/')
    
    ## a) Wild card to match all text files in subdirectory
      file_pattern <- '*.csv' # regex pattern to match the file name format
    
    ## b) Function to read csv files and create wide data frame
      read_log_files <- function(chr15_intx_path, file_name){
        read_csv(paste0(chr15_intx_path, file_name)) %>% 
          mutate(file_name = file_name) 
      }
      
      # read_log_files_test <- function(chr15_intx_path_test, file_name){
      #   read_csv(paste0(chr15_intx_path_test, file_name)) %>% 
      #     mutate(file_name = file_name) 
      # }
      
    ## c) Use map to apply read_log_files function 
      chr15_intx_df <- 
        list.files(chr15_intx_path, pattern = file_pattern) %>% 
        map_df(~ read_log_files(chr15_intx_path, .))
      
      # chr15_intx_df_test <- 
      #   list.files(chr15_intx_path_test, pattern = file_pattern) %>% 
      #   map_df(~ read_log_files_test(chr15_intx_path_test, .))

      
      
###############################################################################
##############             3. Import attribute data              ##############
###############################################################################    
      
  ### 3.1  Import attribute data
    ## a) Import 2015_chr_attributes.csv as df
      chr15_attrib_df <- read_csv(here('raw_data/2015_CHR_attrib',
                                       '2015_chr_attributes.csv'))
      
    ## b) 2640-97175 is misidentified as a female upon capture. 
      # SO update sex from 'f' to 'm' since it's a male who didn't breed
      chr15_attrib_df[chr15_attrib_df$`Band ID` == '2640-97175', 'Sex'] <- 'm'

  
          
###############################################################################
##############          4. Import repro. success data            ##############
###############################################################################    
      
  ### 4.1 Import reproductive success data
    ## a) Import 2015_chr_egg_paternity.csv as df
      chr15_repro_suc <- read_csv(here('raw_data/2015_CHR_repro',
                                        'chr2015_repro_success.csv'))
      
      
  
###############################################################################
##############                5. Export data files               ##############
###############################################################################
      
  ### 5.1 Export data to an RData file     
    ## a) Save and export raw intx data tables and attribute data
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/1_chr15_intx_data.RData'), 
           list = c('chr15_intx_df', 'chr15_attrib_df'))
      # Save test data set
      # save(file = here('raw_data/2015_CHR_soc/test/1_chr15_raw_intx_test_data.RData'), 
      #      list = c('chr15_intx_df_test'))
      

    ## b) Save and export raw repro. success data 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/3_chr15_repro_data.RData'), 
           list = c('chr15_repro_suc'))
      
      
      

    