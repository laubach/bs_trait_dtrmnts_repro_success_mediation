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


  ### PURPOSE: Quantify the social interaction data to create node and dyad 
            # level network measures from pre- manipulation social networks
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Quantify and tidy CHR 2015 node intx
    # 4. Join node social intx to attribute 
    # 5. Quantify and tidy female-male dyad intx 
    # 6. Quantify and tidy female-female dyad intx 
    # 7. Quantify and tidy male-male dyad intx 
    # 8. Export data


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
        
    ## b) Network packages    
      # load igraph package
        library ('igraph')  
        
      # load assortnet package
        library ('assortnet')  
    
    ## c) Other packages    
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
    ## a) Load attribute, reproduction, and pre-manipulation social intx. data
      load(here('data/4_chr_intx_summary_data.RData'))
      

      
###############################################################################
##############      3. Quantify and tidy CHR 2015 node intx      ##############
############################################################################### 
      
  ### 3.1 Quantify and tidy female-male node level social intx data CHR 2015
      # - pre-manipulation
    ## a) create an igraph object from adjacency matrix
      chr15_fxm_intx_graph <- graph_from_adjacency_matrix(chr15_fxm_intx_mat,
                                                         weighted=TRUE,
                                                         mode = 'lower')
      
    ## b) calculate strength (sum of id row and column counts) 
      # pre-manipulation
      chr15_fxm_intx_graph_strength <- igraph::strength(chr15_fxm_intx_graph, 
                                          vids = V(chr15_fxm_intx_graph))
      
    ## c) calculate degree (sum of id row and column cells not zero)
      # pre-manipulation
      chr15_fxm_intx_graph_degree <- igraph::degree(chr15_fxm_intx_graph, 
                                          v = V(chr15_fxm_intx_graph))
    
    ## d) make a tibble of the each Tag ID's strength and degree pre-manipulation
      chr15_fxm_intx_graph_df <- tibble(name=V(chr15_fxm_intx_graph)$
                                    name, chr15_fxm_intx_graph_strength, 
                                    chr15_fxm_intx_graph_degree)
  
    ## e) rename variables in soc_fxm_intx_pre_graph_df 
      chr15_fxm_intx_graph_df <-chr15_fxm_intx_graph_df %>%
        rename(c('Tag' = 'name',
                 'strength.fxm' = 'chr15_fxm_intx_graph_strength',
                 'degree.fxm' = 'chr15_fxm_intx_graph_degree'))
      

  ### 3.2 Quantify and tidy female-female node level social intx data CHR 2015 
      # - pre-manipulation
    ## a) create an igraph object from adjacency matrix
      chr15_fxf_intx_graph <- graph_from_adjacency_matrix(chr15_fxf_intx_mat,
                                                          weighted=TRUE,
                                                          mode = 'lower')
      
    ## b) calculate strength (sum of id row and column counts) 
      # pre-manipulation
      chr15_fxf_intx_graph_strength <- igraph::strength(chr15_fxf_intx_graph, 
                                              vids = V(chr15_fxf_intx_graph))
      
    ## c) calculate degree (sum of id row and column cells not zero)
      # pre-manipulation
      chr15_fxf_intx_graph_degree <- igraph::degree(chr15_fxf_intx_graph, 
                                              v = V(chr15_fxf_intx_graph))
      
    ## d) make a tible of the each Tag ID's strength and degree pre-manipulation
      chr15_fxf_intx_graph_df <- tibble(name=V(chr15_fxf_intx_graph)$
                                          name, chr15_fxf_intx_graph_strength, 
                                        chr15_fxf_intx_graph_degree)
      
    ## e) rename variables in soc_fxf_intx_pre_graph_df 
      chr15_fxf_intx_graph_df <-chr15_fxf_intx_graph_df %>%
        rename(c('Tag' = 'name',
                 'strength.fxf' = 'chr15_fxf_intx_graph_strength',
                 'degree.fxf' = 'chr15_fxf_intx_graph_degree'))
      
      
  ### 3.3 Quantify and tidy male-male node level social intx data CHR 2015
      # - pre-manipulation
    ## a) create an igraph object from adjacency matrix
      chr15_mxm_intx_graph <- graph_from_adjacency_matrix(chr15_mxm_intx_mat,
                                                          weighted=TRUE,
                                                          mode = 'lower')
      
    ## b) calculate strength (sum of id row and column counts) 
      # pre-manipulation
      chr15_mxm_intx_graph_strength <- igraph::strength(chr15_mxm_intx_graph, 
                                              vids = V(chr15_mxm_intx_graph))
      
    ## c) calculate degree (sum of id row and column cells not zero)
      # pre-manipulation
      chr15_mxm_intx_graph_degree <- igraph::degree(chr15_mxm_intx_graph, 
                                            v = V(chr15_mxm_intx_graph))
      
    ## d) make a tible of the each Tag ID's strength and degree pre-manipulation
      chr15_mxm_intx_graph_df <- tibble(name=V(chr15_mxm_intx_graph)$
                                          name, chr15_mxm_intx_graph_strength, 
                                        chr15_mxm_intx_graph_degree)
      
    ## e) rename variables in soc_mxm_intx_pre_graph_df 
      chr15_mxm_intx_graph_df <-chr15_mxm_intx_graph_df %>%
        rename(c('Tag' = 'name',
                 'strength.mxm' = 'chr15_mxm_intx_graph_strength',
                 'degree.mxm' = 'chr15_mxm_intx_graph_degree')) 
      
      
  ### 3.1 Quantify and tidy female-male node level social intx data CHR 2015
      # - pre-manipulation
    ## a) create an igraph object from adjacency matrix
      chr15_fxm_intx_post_graph <- 
        graph_from_adjacency_matrix(chr15_fxm_intx_post_mat,
                                    weighted=TRUE,
                                    mode = 'lower')
      
    ## b) calculate strength (sum of id row and column counts) 
      # pre-manipulation
      chr15_fxm_intx_post_graph_strength <- 
        igraph::strength(chr15_fxm_intx_post_graph, 
                         vids = V(chr15_fxm_intx_post_graph))
      
    ## c) calculate degree (sum of id row and column cells not zero)
      # pre-manipulation
      chr15_fxm_intx_post_graph_degree <- 
        igraph::degree(chr15_fxm_intx_post_graph, 
                       v = V(chr15_fxm_intx_post_graph))
      
    ## d) make a tible of the each Tag ID's strength and degree pre-manipulation
      chr15_fxm_intx_graph_post_df <- tibble(name=V
                                      (chr15_fxm_intx_post_graph)$name, 
                                        chr15_fxm_intx_post_graph_strength, 
                                        chr15_fxm_intx_post_graph_degree)
      
    ## e) rename variables in chr15_fxm_intx_graph_post_df 
      chr15_fxm_intx_graph_post_df <-chr15_fxm_intx_graph_post_df %>%
        rename(c('Tag' = 'name',
                 'strength.fxm.post' = 'chr15_fxm_intx_post_graph_strength',
                 'degree.fxm.post' = 'chr15_fxm_intx_post_graph_degree'))
      
      
      
###############################################################################
##############       4. Join node social intx to attribute       ##############
###############################################################################     

  ### 4.1 Join and tidy CHR 2015 data - pre-manipulation
    ## a) Format variable type
      chr15_attrib_df$Tag <- as.character(chr15_attrib_df$Tag)
  
    ## b) Left join chr15_fxm_intx_graph_df to chr15_attrib_df
      chr15_attrib_df <- chr15_attrib_df  %>%
        left_join(chr15_fxm_intx_graph_df,
                  by = c('Tag' = 'Tag'), 
                  copy = F) 
      
    ## c) Left join chr15_fxf_intx_graph_df to chr15_attrib_df
      chr15_attrib_df <- chr15_attrib_df  %>%
        left_join(chr15_fxf_intx_graph_df,
                  by = c('Tag' = 'Tag'), 
                  copy = F)
      
    ## d) Left join chr15_mxm_intx_graph_df to chr15_attrib_df
      chr15_attrib_df <- chr15_attrib_df  %>%
        left_join(chr15_mxm_intx_graph_df,
                  by = c('Tag' = 'Tag'), 
                  copy = F)
      
  
  ### 4.2 Join and tidy CHR 2015 data - post manipulation
    ## a) Format variable type
      chr15_attrib_post_df$Tag <- as.character(chr15_attrib_post_df$Tag)
      
    ## b) Left join chr15_fxm_intx_graph_post_df to chr15_attrib_post_df
      chr15_attrib_post_df <- chr15_attrib_post_df  %>%
        left_join(chr15_fxm_intx_graph_post_df,
                  by = c('Tag' = 'Tag'), 
                  copy = F) 
      

    # ## c) Normalize chr2015 degree and strength based on:
    #   # Metrics for network comparison using egonet feature distributions
    #   # by Carlo Piccardi - 2023
    #   # https://www.nature.com/articles/s41598-023-40938-4
    #   chr15_attrib_df <- chr15_attrib_df  %>%
    #     mutate(norm.strength = round(((strength.fxm - min(strength.fxm, na.rm =T))/
    #                                     (max(strength.fxm, na.rm = T) - 
    #                                        min(strength.fxm, na.rm = T))), 2),
    #            norm.deg = round(((degree.fxm - min(degree.fxm, na.rm = T))/
    #                                (max(degree.fxm, na.rm = T) - 
    #                                   min(degree.fxm, na.rm = T))), 2)) 
      
      
  
  ### 4.3 Manual data cleaning
      
#*****************************************************************************#
#************************* Manual data cleaning ******************************#
  # NOTE: Tag 116, female 2640-97117 is not in the attribute data because she
      # incomplete fecundity (she left site after egg removal) 
    ## a) remove bird tag 116 chr15_attrib_df and chr15_attrib_post_df
      chr15_attrib_df <- chr15_attrib_df %>%
        filter(!(Tag == 116))
      
      chr15_attrib_post_df <- chr15_attrib_post_df %>%
        filter(!(Tag == 116))
     
  # NOTE: Tag 79 for male 2640-97175 does not pair until late Aug...
      # fundamentally different than other birds, so remove it here.
      # He had incomplete paternity - never included in paternity analysis)
    ## b) remove bird tag 79 from chr15_attrib_df and chr15_attrib_post_df
      chr15_attrib_df <- chr15_attrib_df %>%
        filter(!(Tag == 79))
      
      chr15_attrib_post_df <- chr15_attrib_post_df %>%
        filter(!(Tag == 79))

  # NOTE: Tag 57, for male 1921-30295 does not pair until late Aug...
      # fundamentally different than other birds, so remove it here.
      # He had incomplete paternity - never included in paternity analysis)
    ## c) remove bird tag 57 from chr15_attrib_df and chr15_attrib_post_df
      chr15_attrib_df <- chr15_attrib_df %>%
        filter(!(Tag == 57))
      
      chr15_attrib_post_df <- chr15_attrib_post_df %>%
        filter(!(Tag == 57))
      
      
  # Note: males 2640-97153, 2640-97157, and 2640-97158 were never tagged,
    # so not in attribute table and dropped from repro_suc data upon left join
    # to chr15_attrib_df
      
#*****************************************************************************#
#*****************************************************************************#

      
        
###############################################################################
##############                   5. Export data                  ##############
###############################################################################
      
  ### 5.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export node level data
      save(file = here('data/5_6_chr_node_df.RData'), 
           list = c('chr15_attrib_df', 'chr15_attrib_post_df',
                    'chr15_fxm_intx_graph', 'chr15_fxm_intx_mat', 
                    'chr15_fxf_intx_graph', 'chr15_fxf_intx_mat', 
                    'chr15_mxm_intx_graph', 'chr15_mxm_intx_mat',
                    'chr15_fxm_intx_post_graph', 'chr15_fxm_intx_post_mat'))
      
      
  
      


      
      
    