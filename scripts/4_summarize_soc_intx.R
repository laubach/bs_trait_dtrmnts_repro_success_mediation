################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############         4. Summarize social interaction data         #############
#############                                                      #############
#############                  By: Zach Laubach                    #############
#############                created: 21 Aug 2024                  #############
#############             last updated: 25 Feb 2024                #############
################################################################################


  ### PURPOSE: Summarize the social interaction data to create level
            # network measures from pre- and post manipulation social networks
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Tidy data 
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
      
    # load igraph package
      library ('igraph')  
    
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
    ## a) Load pre-manipulation social intx. data for CHR 2015
      load(here('data/4_chr15_intx_pre_exp_20rssi.RData'))
      
    ## b) Load post-manipulation social intx. data for CHR 2015
      load(here('data/4_chr15_intx_post_exp_20rssi.RData'))
      
    ## c) Load adult attribute data for CHR 2015
      load(here('data/4_chr15_attrib_data.RData'))
      
      

###############################################################################
##############                   3. Tidy data                    ##############
############################################################################### 
      
  ### 3.1 Summarize level social intx    
    ## a) Rename social intx. variable pre-manipulation
      chr15_intx_pre_df_20rssi <- chr15_intx_pre_df_20rssi %>%
        rename(cnt.soc.intx.pre = n)
      
    ## b) Summary of each female's social intx pre-manipulation
      female_intx_pre_sum <- chr15_intx_pre_df_20rssi %>%
        group_by(Tag1) %>%
                  # number of unique social partners
        summarise(n.unique.soc.partner.pre = sum(!is.na(cnt.soc.intx.pre)),
                  # max number of social interactions with any one soc. partner
                   max.cnt.soc.intx.pre = round(max(cnt.soc.intx.pre,
                                                na.rm = T), 2),
                  # total number of social interactions with all soc. partners
                   tot.cnt.soc.intx.pre = round (sum(cnt.soc.intx.pre, 
                                                  na.rm = T), 2),
                  # total duration of interactions with all soc. partners
                  tot.dur = round(sum(tot.dur, na.rm = T)/60, 2)) %>%
        rename(Tag = Tag1)
      
    ## c) Summary of each male's social intx pre-manipulation
      male_intx_pre_sum <- chr15_intx_pre_df_20rssi %>%
        group_by(Tag2) %>%
                  # number of unique social partners
        summarise(n.unique.soc.partner.pre = sum(!is.na(cnt.soc.intx.pre)),
                   # max number of social interactions with any one soc. partner
                   max.cnt.soc.intx.pre = round(max(cnt.soc.intx.pre,
                                                na.rm = T), 2),
                   # total number of social interactions with all soc. partners
                   tot.cnt.soc.intx.pre = round (sum(cnt.soc.intx.pre, 
                                                 na.rm = T), 2),
                   # total duration of interactions with all soc. partners
                   tot.dur = round(sum(tot.dur, na.rm = T)/60, 2)) %>%
        rename(Tag = Tag2)
      
    ## d) Combine female and male soc. intx. summaries pre-manipulation
      intx_pre_sum <- bind_rows(female_intx_pre_sum, 
                                 male_intx_pre_sum)
      
    ## e) Left join intx pre-manipulation data to the chr15_attrib_pre_df
      chr15_attrib_pre_df <- chr15_attrib_pre_df  %>%
        left_join(intx_pre_sum,
                  by = c('Tag' = 'Tag'), 
                  copy = F) 
      
    ## f) Rename social intx. variable post-manipulation
      chr15_intx_post_df_20rssi <- chr15_intx_post_df_20rssi %>%
        rename(cnt.soc.intx.post = n)
      
    ## g) Summary of each female's social intx post-manipulation
      female_intx_post_sum <- chr15_intx_post_df_20rssi %>%
        group_by(Tag1) %>%
        # number of unique social partners
        summarise (n.unique.soc.partner.post = sum(!is.na(cnt.soc.intx.post)),
                   # max number of social interactions with any one soc. partner
                   max.cnt.soc.intx.post = round(max(cnt.soc.intx.post,
                                                na.rm = T), 2),
                   # total number of social interactions with all soc. partners
                   tot.cnt.soc.intx.post = round (sum(cnt.soc.intx.post, 
                                                 na.rm = T), 2),
                   # total duration of interactions with all soc. partners
                   tot.dur = round (sum(tot.dur, na.rm = T)/60, 2)) %>%
        rename(Tag = Tag1)
      
    ## h) Summary of each male's social intx post-manipulation
      male_intx_post_sum <- chr15_intx_post_df_20rssi %>%
        group_by(Tag2) %>%
        # number of unique social partners
        summarise (n.unique.soc.partner.post = sum(!is.na(cnt.soc.intx.post)),
                   # max number of social interactions with any one soc. partner
                   max.cnt.soc.intx.post = round(max(cnt.soc.intx.post,
                                                na.rm = T), 2),
                   # total number of social interactions with all soc. partners
                   tot.cnt.soc.intx.post = round (sum(cnt.soc.intx.post, 
                                                 na.rm = T), 2),
                   # total duration of interactions with all soc. partners
                   tot.dur = round (sum(tot.dur, na.rm = T)/60, 2)) %>%
        rename(Tag = Tag2)
      
    ## i) Combine female and male soc. intx. summaries post-manipulation
      intx_post_sum <- bind_rows(female_intx_post_sum, 
                                male_intx_post_sum)
      
    ## j) Left join intx post-manipulation data to the chr15_attrib_post_df
      chr15_attrib_post_df <- chr15_attrib_post_df  %>%
        left_join(intx_post_sum,
                  by = c('Tag' = 'Tag'), 
                  copy = F) 

    ## k) change Tag from double to character
      chr15_attrib_pre_df$Tag <- as.character(chr15_attrib_pre_df$Tag)
      chr15_attrib_post_df$Tag <- as.character(chr15_attrib_post_df$Tag)
 
      
  ### 3.3 Calculate social network measures
    ## a) create an igraph object from pre-manipulation adjacency matrix
      soc_intx_pre_graph <- graph_from_adjacency_matrix(chr15_intx_pre_20rssi_mat, 
                                                    weighted=TRUE,
                                                    mode = 'lower')
      
    ## b) calculate strength (sum of id row and column counts) 
      # pre-manipulation
      soc_intx_pre_graph_strength <- igraph::strength(soc_intx_pre_graph, 
                                          vids = V(soc_intx_pre_graph))
      
    ## c) calculate degree (sum of id row and column cells not zero)
      # pre-manipulation
      soc_intx_pre_graph_degree <- igraph::degree(soc_intx_pre_graph, 
                                          v = V(soc_intx_pre_graph))
    
    ## d) make a tible of the each Tag ID's strength and degree pre-manipulation
      soc_intx_pre_graph_df <- tibble(name=V(soc_intx_pre_graph)$
                                    name, soc_intx_pre_graph_strength, 
                                  soc_intx_pre_graph_degree)
  
    ## e) rename variables in soc_intx_pre_graph_df 
      soc_intx_pre_graph_df <-soc_intx_pre_graph_df %>%
        rename(c('Tag' = 'name',
                 'strength.pre' = 'soc_intx_pre_graph_strength',
                 'degree.pre' = 'soc_intx_pre_graph_degree'))
      
    ## f) create an igraph object from post-manipulation adjacency matrix
      soc_intx_post_graph <- graph_from_adjacency_matrix(chr15_intx_post_20rssi_mat, 
                                                        weighted=TRUE,
                                                        mode = 'lower')
      
    ## g) calculate  strength (sum of id row and column counts) 
      # post-manipulation
      soc_intx_post_graph_strength <- igraph::strength(soc_intx_post_graph, 
                                              vids = V(soc_intx_post_graph))
      
    ## h) calculate  degree (sum of id row and column cells not zero)
      # post-manipulation
      soc_intx_post_graph_degree <- igraph::degree(soc_intx_post_graph, 
                                          v = V(soc_intx_post_graph))
      
    ## i) make a tible of the each Tag ID's strength and degree post-manipulation
      soc_intx_post_graph_df <- tibble(name=V(soc_intx_post_graph)$
                                        name, soc_intx_post_graph_strength, 
                                      soc_intx_post_graph_degree)
      
    ## j) rename variables in soc_intx_post_graph_df 
      soc_intx_post_graph_df <-soc_intx_post_graph_df %>%
        rename(c('Tag' = 'name',
                 'strength.post' = 'soc_intx_post_graph_strength',
                 'degree.post' = 'soc_intx_post_graph_degree'))

    ## k) Left join soc_intx_pre_graph_df to chr15_attrib_pre_df
      chr15_attrib_pre_df <- chr15_attrib_pre_df  %>%
        left_join(soc_intx_pre_graph_df,
                  by = c('Tag' = 'Tag'), 
                  copy = F) 
    
    ## l) replace any NA in pre_df for specified columns
      chr15_attrib_pre_df <- chr15_attrib_pre_df  %>%
        mutate_at(vars('n.unique.soc.partner.pre':'tot.dur'), 
                  ~replace_na(., 0))
      
    ## m) Left join soc_intx_post_graph_df to chr15_attrib_post_df
      chr15_attrib_post_df <- chr15_attrib_post_df  %>%
        left_join(soc_intx_post_graph_df,
                  by = c('Tag' = 'Tag'), 
                  copy = F) 
      
    ## n) replace any NA in post_df for specified columns
      chr15_attrib_post_df <- chr15_attrib_post_df  %>%
        mutate_at(vars('n.unique.soc.partner.post':'tot.dur'), 
                  ~replace_na(., 0))
      

   
###############################################################################
##############                   4. Export data                  ##############
###############################################################################
      
  ### 4.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015  level mediation analysis
      save(file = here('data/5_6_chr15_mediation_data.RData'), 
           list = c('chr15_attrib_pre_df', 
                    'chr15_attrib_post_df', 'soc_intx_pre_graph', 
                    'soc_intx_post_graph', 'chr15_intx_pre_20rssi_mat', 
                    'chr15_intx_post_20rssi_mat'))

      
  ### 4.2 Export data to .csv files 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) data for CHR 2015 pre-manipulation  level
      # mediation analysis to .csv file
      write.csv(chr15_attrib_pre_df, 
                file = here('data/6_chr15_mediation_pre_manip_data.csv'), 
                row.names = T)
      
    ## b) data for CHR 2015 pre-manipulation level
      # mediation analysis to .csv file
      write.csv(chr15_attrib_post_df, 
                file = here('data/6_chr15_mediation_post_manip_data.csv'), 
                row.names = T)
      
      
    