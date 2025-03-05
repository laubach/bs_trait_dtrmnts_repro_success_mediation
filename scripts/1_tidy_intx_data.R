################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############              1. Tidy interaction data                #############
#############                                                      #############
#############                  By: Zach Laubach                    #############
#############                created: 7 Aug 2024                   #############
#############             last updated: 23 Feb 2025                #############
################################################################################


  ### PURPOSE: Tidy interaction data and build adjacency matrices for swallow 
          #  pre- and post experimental manipulation at >= 20, 25, and 30 RSSI 
  
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Tidy data
    # 4. Select type of interactions
    # 5. Finalize data selection
    # 6. Build adjacency matrices 
    # 7. Export data 
  


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
    # Running under: macOS Sequoia 15.1.1
    
  
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
  
    ## a) Load raw intx data CHR 2015
      load(here('data/1_chr15_intx_data.RData'))
      # test data
      #load(here('raw_data/2015_CHR/test/1_chr15_raw_intx_test_data.RData'))


      
###############################################################################
##############                   3. Test data                    ##############
###############################################################################       

  #  ### 3.1 Recreate IIL 2015 'C_subset folder by duration.R' script
  #   ## a) Path name
  #      path <- here('raw_data/2015_CHR/test/')
  # 
  #   ## b) Gets a list of all csv files in the folder.
  #     files <- list.files(path=path, pattern="*.csv")
  # 
  #   ## c) Subset the data
  #     for(k in 1:2){
  #         #import data
  #         #data<-read.csv(files[k])
  # 
  #         data <- read.csv(paste0(path, files[k]))
  # 
  #         d.1<-subset(data, data[,9]>=40) #change for subsetting
  # 
  #         t1<-d.1[1,1]
  #         t2<-d.1[1,2]
  # 
  #         #names files
  #         file<-as.character(paste(t1,"b46_2320to40",t2,".csv",sep="")) #dont bother naming meaningful if not final subset
  #         write.csv(d.1,file.path(paste0(path,'filter/',file)),row.names=F)
  #       }
  # 
  #   ### 3.2 Recreate IIL 2015 'D_Counting Interactions and building a matrix.R'
  #     # script
  #     ## a) Make a list of files
  #       files <- list.files(path=paste0(path,'filter/'), pattern="*.csv")
  # 
  #     ## b) Create empty data frame
  #       df1 = data.frame(NA,NA,NA)
  #       names(df1) = c('thisID','enID','Count')
  # 
  #     ## c) set global values
  #       n<-0
  # 
  #     ## d) For loop to count the rows of data / number of intx post filtering
  #       for(k in 1:length(files)){
  # 
  #         #import data
  #         #data<-read.csv(files[k])
  #         data <- read.csv(paste0(path,'filter/', files[k]))
  #         if(is.na(data[1,1])==TRUE){
  #           }
  #         else{
  #           n<-n+1
  #           df1[n,1]<-data[1,1]
  #           df1[n,2]<-data[1,2]
  #           df1[n,3]<-length(data[,1])
  #         }
  #       }
  # 
  #   ## e) view new data frame of intx counts
  #       df1
  # 
  #   ## f) extract ID list from thisID and enID
  #     id_list <- as.character(sort(union(df1$thisID, df1$enID), decreasing = F))
  #     # create an empty matrix (n x n) to store intx counts
  #     m1<-matrix(,nrow=length(id_list),ncol=length(id_list))
  # 
  #     # name matrix rows and columns
  #     rownames(m1)<-id_list
  #     colnames(m1)<-id_list
  # 
  # 
  #   ## g) populate empty matrix from df1
  #     for(i in 1:length(df1[,1])){
  #       print(i)
  #       m1[as.character(df1[i,1]),as.character(df1[i,2])]<-df1[i,3]
  #       m1[as.character(df1[i,2]),as.character(df1[i,1])]<-df1[i,3]
  # 
  #     }
  # 
  #   ## h) Replace NA with 0
  #     m1[is.na(m1)] <- 0
  # 
  # 
  # ### 3.3 Check code update pipeline against prior pipeline
  #   ## a) Subset the data and the select variables of interest
  #     chr15_intx_df_test_filter <- chr15_intx_df_test %>%
  #       filter(duration >=40) %>%
  #       select(thisID, enID)
  # 
  #   ## b) Count the dyadic interactions that pass filtering
  #     chr15_intx_df_intx_test <- chr15_intx_df_test_filter %>%
  #       pmap_dfr(~list(...)[order(c(...))] %>% set_names(names(chr15_intx_df_test_filter))) %>%
  #       group_by_all %>%
  #       count
  # 
  #     # df1
  #     #       thisID enID Count
  #     #   1    104  107    35
  #     #   2    112  104    21
  #     #
  #     #         VS
  #     #
  #     # chr15_intx_df_intx_test
  #     #   A tibble: 2 × 3
  #     # # Groups:   thisID, enID [2]
  #     #       thisID  enID     n
  #     #       <dbl> <dbl> <int>
  #     #   1    104   107    35
  #     #   2    104   112    21

#************************************ NOTE ************************************#
#*************** chr15_intx_df_intx_test matches counts in df1 *****************
#************************************ NOTE ************************************#   
   

###############################################################################
##############                   3. Tidy data                    ##############
###############################################################################    
 
  ### 3.1 Tidy attribute data
    ## a) Format data to all lower case
      chr15_attrib_df <- AllCharactersToLower(chr15_attrib_df)
    
    # ## b) Format variable names to standard R format (sep = .)
    #   chr15_attrib_df <- FormatVarNames(chr15_attrib_df)
    
    ## b) Check variable type
      str(chr15_attrib_df)
    
    ## c) Re-code 'Sex' as *nominal* factor (with ordered levels)
      # Set levels (odering) of variable 
      chr15_attrib_df <- transform(chr15_attrib_df, 
                                 Sex = factor(Sex,
                                              levels = c('f', 'm')))
    
    ## d) Re-code 'Age.category' as *nominal* factor (with ordered levels)
      # Set levels (odering) of variable 
      chr15_attrib_df <- transform(chr15_attrib_df, 
                                 Age.category = factor(Age.category,
                                                    levels = c('sy', 'asy')))
    
      
  ### 3.2 Tidy chr15_intx_df data
      ## a) Rename variables in chr15_intx_df
      chr15_intx_df <- chr15_intx_df %>%
        rename(c(Tag1 = thisID,
                 Tag2 = enID))
      
#*****************************************************************************#
#************************* Manual data cleaning ******************************#
  # # NOTE: Tag ID 116 is not in the attribute data, so remove it here.                
  #   ## b) remove bird tag 116
  #     chr15_intx_df <- chr15_intx_df %>%
  #       filter(!(Tag1 == 116 | Tag2 == 116)) 
  #     
  # # NOTE: Tag ID 57 does not pair until late Aug...fundamentally different
  #     # than other birds, so remove it here.                
  #   ## b) remove bird tag 57 from chr15_intx_df and chr15_attrib_df
  #     chr15_intx_df <- chr15_intx_df %>%
  #       filter(!(Tag1 == 57 | Tag2 == 57))   
  #     
  #     chr15_attrib_df <- chr15_attrib_df %>%
  #       filter(!(Tag == 57))   
  #     
#*****************************************************************************#
#*****************************************************************************#
    
      
      
###############################################################################
##############           4. Select type of interactions          ##############
############################################################################### 


  ### 4.1 Subset pairwise interaction data to include female x male intx only
    ## a) Make a list of female ids
      fem_ids <- chr15_attrib_df[chr15_attrib_df$Sex == 'f', 'Tag']
    
    ## b) Make a list of male ids
      male_ids <- chr15_attrib_df[chr15_attrib_df$Sex == 'm', 'Tag']
    
    # ## c) Make a list of females' social mates
    #   fem_mate_ids <- as.numeric(chr15_attrib_df[chr15_attrib_df$Sex == 'f' , 
    #                                            'Mate.tag'])
    
    ## c) intersect female list with Tag1 chr15_intx_df and where
      # true, populate sex.1 == f, else == m
      chr15_intx_df$sex.1 <- ifelse(chr15_intx_df$Tag1 %in% 
                                             fem_ids, 'f', 'm')
    
    ## d) intersect male list with Tag2 chr15_intx_df and where
      # true, populate sex.2 == m, else f
      chr15_intx_df$sex.2 <- ifelse(chr15_intx_df$Tag2 %in% 
                                             male_ids, 'm', 'f')
    
    ## e) Reorder variables so all females are Tag1 and males are Tag2
      # first, filter males in Tag1 position and switch Tag.# and sex.# order 
      # of variable names
        chr15_intx_df_m <- chr15_intx_df %>%
          filter(sex.1 == 'm') %>%
          rename(Tag1 = Tag2,
                 Tag2 = Tag1,
                 sex.1 = sex.2,
                 sex.2 = sex.1)
    
      # second, reorder columns
        col_order <- c('Tag1', 'Tag2', 'Tstart', 'Tend', 'dyadic', 'RSSImax',
                       'RSSImin', 'RSSImean', 'duration', 'file_name',
                       'sex.1', 'sex.2')
      
        chr15_intx_df_m <- chr15_intx_df_m[, col_order]
    
      # third, filter already in Tag1 position females
        chr15_intx_df <- chr15_intx_df %>%
          filter(sex.1 == 'f') 
    
      # fourth, bind the two reordered data frames together
        chr15_intx_df <- bind_rows(chr15_intx_df, 
                                   chr15_intx_df_m)
  
    ## f) Filter by female x male intx
      chr15_intx_df <- chr15_intx_df %>%
          filter(sex.1 == 'f' & sex.2 == 'm') 
      
    ## g) double check only females in Tag1 sex.1 and males in Tag2 sex.2
      # both should be T, if all T
      all(chr15_intx_df$Tag1 %in% fem_ids)
      all(chr15_intx_df$Tag2 %in% male_ids)
      

  ### 4.2 Determine when all tags active pre-manip
    ## a) Create df of all the earliest and latest TAG signals pre-manip for 
      # females 
      chr15_intx_pre_date_range_f <- chr15_intx_df %>%
        group_by(Tag1) %>%
        filter(day(Tstart) <= 15 & month(Tstart) == 6) %>%
        summarise(start.date.time = min(Tstart),
                  end.date.time = max(Tend)) %>%
        ungroup()
      # and males
      chr15_intx_pre_date_range_m <- chr15_intx_df %>%
        group_by(Tag2) %>%
        filter(day(Tstart) <= 15 & month(Tstart) == 6) %>%
        summarise(start.date.time = min(Tstart),
                  end.date.time = max(Tend)) %>%
        ungroup()
      
      
#*****************************************************************************#
#************************* Manual data cleaning ******************************#
    # NOTE: Check the earliest vs latest encounter net to ensure networks
      # are based on when all Tags working
      
    ## a) remove female 49 from both pre and post manip networks - 
      # Tag worked from 6/13-6/14
      chr15_intx_pre_df <- chr15_intx_df %>%
        filter(!(Tag1 == 49 | Tag2 == 49)) 
      
   # NOTE: Tag 108 for female 2640-97169 never worked, so not included in 
      # either pre or post network or attribute file
      
#*****************************************************************************#
#*****************************************************************************#     
     
       
  ### 4.3 Determine when all tags active post-manip
    ## a) Create df of all the earliest and latest TAG signals post-manip for
      # females
      chr15_intx_post_date_range_f <- chr15_intx_df %>%
        group_by(Tag1) %>%
        filter(day(Tstart) >= 19 & month(Tstart) == 6) %>%
        summarise(start.date.time = min(Tstart),
                  end.date.time = max(Tend)) %>%
        ungroup()
      # and males
      chr15_intx_post_date_range_m <- chr15_intx_df %>%
        group_by(Tag2) %>%
        filter(day(Tstart) >= 19 & month(Tstart) == 6) %>%
        summarise(start.date.time = min(Tstart),
                  end.date.time = max(Tend)) %>%
        ungroup()
      
#*****************************************************************************#
#************************* Manual data cleaning ******************************#
    # NOTE: Check the earliest vs latest encounter net to ensure post manip
      # networks are based on when all Tags working
      
    ## b) Remove males 57, 67, 61, and 55 as well as females 40, 51, 52, 116 
      # from post manip networks...Tags stopped on or before 6/19
      chr15_intx_post_df <- chr15_intx_pre_df %>%
        filter(!(Tag1 == 55 | Tag2 == 55)) %>% # males
        filter(!(Tag1 == 57 | Tag2 == 57)) %>%
        filter(!(Tag1 == 61 | Tag2 == 61)) %>%
        filter(!(Tag1 == 67 | Tag2 == 67)) %>%
        filter(!(Tag1 == 40 | Tag2 == 40)) %>% # females
        filter(!(Tag1 == 51 | Tag2 == 51)) %>%
        filter(!(Tag1 == 52 | Tag2 == 52)) %>%
        filter(!(Tag1 == 116 | Tag2 == 116)) 
      
      
  # NOTE: Tag Tag 56 for male 2600-27559, Tag 82 for male 1741-45779 (left site), 
      # Tag 80 for female 2600-54704, which have no Tag data by 6/19, so
      # they were removed by date filtering
   
#*****************************************************************************#
#*****************************************************************************#      
   
         
      
###############################################################################
##############            5. Finalize data selection             ##############
############################################################################### 
      
  ### 5.1. Subset the pre-manip data by tag signal and select 
      # variables of interest
    ## a) Subset CHR 2015 the data by pre-experiment dates: 6/13, 6/14, 6/15
      chr15_intx_pre_df <- chr15_intx_pre_df %>%
        filter(day(Tstart) <= 15 & month(Tstart) == 6)
      
    ## b)  extract ID list from Tag1 and Tag2 of of IDs to include in 
      # pre-manipulation adjacency matrix.
      # Note this is done before RSSI filtering in case any IDs have no 
      # interactions
      id_list_chr15_intx_pre <- as.character(sort(union(chr15_intx_pre_df$Tag1, 
                                                       chr15_intx_pre_df$Tag2), 
                                                          decreasing = F))
      
    ## c) Subset CHR 2015 pre-manipulation by RSSI mean >= 20
      chr15_intx_pre_df_20rssi <- chr15_intx_pre_df %>%
        filter(RSSImean >= 20) %>%
        select(Tag1, Tag2)  
      
    ## d) create a second data frame that contains the duration of intx at
        # 20RSSI
      chr15_dur_pre_df_20rssi <- chr15_intx_pre_df %>%
        filter(RSSImean >= 20) %>%
        select(Tag1, Tag2, duration)  
      
  
  ### 5.2  Subset the post-manip data by tag signal and select 
      # variables of interest    
    ## a) Subset CHR 2015 the data by post-experiment dates: 6/19 - 6/20
      chr15_intx_post_df <- chr15_intx_post_df %>%
        filter(day(Tstart) == 19 | day(Tstart) == 20 & 
                 month(Tstart) == 6)
      
    ## b)  extract ID list from Tag1 and Tag2 of of IDs to include in 
      # post manipulation adjacency matrix.
      # Note this is done before RSSI filtering in case any IDs have no 
      # interactions
      id_list_chr15_intx_post <- as.character(sort(union(chr15_intx_post_df$Tag1, 
                                                         chr15_intx_post_df$Tag2), 
                                                  decreasing = F))
      
    ## c) Subset CHR 2015 post-manipulation by RSSI mean >= 20
      chr15_intx_post_df_20rssi <- chr15_intx_post_df %>%
        filter(RSSImean >= 20) %>%
        select(Tag1, Tag2)
      
    ## d) create a second data frame that contains the duration of intx at
      # 20RSSI
      chr15_dur_post_df_20rssi <- chr15_intx_post_df %>%
        filter(RSSImean >= 20) %>%
        select(Tag1, Tag2, duration)  
      

      
      
###############################################################################
##############            6. Build adjacency matrices            ##############
###############################################################################          

  ### 6.1 Build CHR 2015 adjacency matrix for >= 20 mean RSSI pre-manipulation
    ## a) Count CHR 2015 dyadic interactions >= 20 mean RSSI
      chr15_intx_pre_df_20rssi <- chr15_intx_pre_df_20rssi %>%
        # pmap_dfr(~list(...)[order(c(...))] %>% 
        #            set_names(names(chr15_intx_pre_df_20rssi))) %>%
        group_by_all %>% 
        count %>%
        as.data.frame # needs to be a data frame and not tibble 
      # counts should sum to 1077
      sum(chr15_intx_pre_df_20rssi$n)
      
    ## b) summarize the total duration each pair spends together
      chr15_dur_pre_df_20rssi <- chr15_dur_pre_df_20rssi %>%
        # pmap_dfr(~list(...)[order(c(...))] %>% 
        #            set_names(names(chr15_intx_pre_df_20rssi))) %>%
        group_by(across(all_of(c('Tag1', 'Tag2')))) %>% 
        summarise(tot.dur = sum(duration)) %>%
        as.data.frame # needs to be a data frame and not tibble 

    ## c) left join tot.dur to the intx data frame  
      chr15_intx_pre_df_20rssi <- chr15_intx_pre_df_20rssi %>%
        left_join(select(chr15_dur_pre_df_20rssi, c(Tag1, Tag2, tot.dur)),
                  by = c('Tag1' = 'Tag1', 'Tag2' = 'Tag2'), 
                  copy = F) 
    
    ## d) create an empty matrix (n x n) to store intx counts using ID list
     # created in 5.1 b)
      chr15_intx_pre_20rssi_mat <- matrix(, nrow = length(id_list_chr15_intx_pre),
                                   ncol = length(id_list_chr15_intx_pre))
      
    ## e) name matrix rows and columns
      rownames(chr15_intx_pre_20rssi_mat)<-id_list_chr15_intx_pre
      colnames(chr15_intx_pre_20rssi_mat)<-id_list_chr15_intx_pre
    
    ## f) populate empty matrix from chr15_intx_pre_df_20rssi
      for(i in 1:length(chr15_intx_pre_df_20rssi[,1])){
        print(i)
        chr15_intx_pre_20rssi_mat[as.character(chr15_intx_pre_df_20rssi[i,1]),
                         as.character(chr15_intx_pre_df_20rssi[i,2])]<-chr15_intx_pre_df_20rssi[i,3]
        chr15_intx_pre_20rssi_mat[as.character(chr15_intx_pre_df_20rssi[i,2]),
                         as.character(chr15_intx_pre_df_20rssi[i,1])]<-chr15_intx_pre_df_20rssi[i,3]
      }
      
    ## g) Replace NA with 0
      chr15_intx_pre_20rssi_mat[is.na(chr15_intx_pre_20rssi_mat)] <- 0
 

  ### 6.2 Build CHR 2015 adjacency matrix for >= 20 mean RSSI post manipulation
    ## a) Count CHR 2015 dyadic interactions >= 20 mean RSSI
      chr15_intx_post_df_20rssi <- chr15_intx_post_df_20rssi %>%
        # pmap_dfr(~list(...)[order(c(...))] %>% 
        #            set_names(names(chr15_intx_post_df_20rssi))) %>%
        group_by_all %>% 
        count %>%
        as.data.frame # needs to be a data frame and not tibble 
      # counts should sum to 405
      sum(chr15_intx_post_df_20rssi$n)
      
    ## b) summarize the total duration each pair spends together
      chr15_dur_post_df_20rssi <- chr15_dur_post_df_20rssi %>%
        # pmap_dfr(~list(...)[order(c(...))] %>% 
        #            set_names(names(chr15_intx_pre_df_20rssi))) %>%
        group_by(across(all_of(c('Tag1', 'Tag2')))) %>% 
        summarise(tot.dur = sum(duration)) %>%
        as.data.frame # needs to be a data frame and not tibble 
      
    ## c) left join tot.dur to the intx data frame
      chr15_intx_post_df_20rssi <- chr15_intx_post_df_20rssi %>%
        left_join(select(chr15_dur_post_df_20rssi, c(Tag1, Tag2, tot.dur)),
                  by = c('Tag1' = 'Tag1', 'Tag2' = 'Tag2'), 
                  copy = F) 
      
    ## d) create an empty matrix (n x n) to store intx counts
      chr15_intx_post_20rssi_mat <- matrix(, nrow = length(id_list_chr15_intx_post),
                                          ncol = length(id_list_chr15_intx_post))
      
    ## e) name matrix rows and columns
      rownames(chr15_intx_post_20rssi_mat)<-id_list_chr15_intx_post
      colnames(chr15_intx_post_20rssi_mat)<-id_list_chr15_intx_post
      
    ## f) populate empty matrix from chr15_intx_post_df_20rssi
      for(i in 1:length(chr15_intx_post_df_20rssi[,1])){
        print(i)
        chr15_intx_post_20rssi_mat[as.character(chr15_intx_post_df_20rssi[i,1]),
                                  as.character(chr15_intx_post_df_20rssi[i,2])]<-chr15_intx_post_df_20rssi[i,3]
        chr15_intx_post_20rssi_mat[as.character(chr15_intx_post_df_20rssi[i,2]),
                                  as.character(chr15_intx_post_df_20rssi[i,1])]<-chr15_intx_post_df_20rssi[i,3]
      }
      
    ## g) Replace NA with 0
      chr15_intx_post_20rssi_mat[is.na(chr15_intx_post_20rssi_mat)] <- 0
      


###############################################################################
##############                  7. Export data                   ##############
###############################################################################
      
  ### 7.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export attribute data 
        # Files are saved in the 'data' folder in the working directory as an
        # RData file.
      save(file = here('data/2_chr15_attribute_data.RData'),
             list = c('chr15_attrib_df'))
      
    ## b) Save and export data for CHR 2015 pre-manipulation >= 20 RSSI
      save(file = here('data/4_chr15_intx_pre_exp_20rssi.RData'), 
           list = c('chr15_intx_pre_df', 'chr15_intx_pre_df_20rssi', 
                    'chr15_intx_pre_20rssi_mat'))
      
    ## c) Save and export data for CHR 2015 post-manipulation >= 20 RSSI
      save(file = here('data/4_chr15_intx_post_exp_20rssi.RData'), 
           list = c('chr15_intx_post_df', 'chr15_intx_post_df_20rssi', 
                    'chr15_intx_post_20rssi_mat'))
      
   
  ### 7.2 Export adjacency matrices to .csv files 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 pre-manipulation >= 20 RSSI
      write.csv(chr15_intx_pre_20rssi_mat, 
                file = here('data/4_chr15_intx_pre_exp_20rssi.csv'), 
                row.names = T)
      
    ## b) Save and export data for CHR 2015 post-manipulation >= 20 RSSI
      write.csv(chr15_intx_post_20rssi_mat, 
                file = here('data/4_chr15_intx_post_exp_20rssi.csv'), 
                row.names = T)
      

      
    