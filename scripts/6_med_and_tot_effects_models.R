################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############        6. Mediation and total effects models         #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 29 Sept 2024                #############
#############               last updated: 6 Dec 2024               #############
################################################################################



### PURPOSE: Run mediation and total effects models for traits and reproductive
           # success
         
  # Code Blocks
    # 1. Configure work space
    # 2. Load RData
    # 3. Trait associations with reproductive success
    # 4. Trait associations with soc. net. metrics 
    # 5. Age association with plumage traits 
    # 6. Soc. net. metrics association with reproduction
    # 7. Age by soc. net. metrics interaction
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
      
      # load broom package
        library ('broom')
      
      # load rempsyc to make publication ready tables
      # see https://cran.r-project.org/web/packages/rempsyc/vignettes/table.html
      library('rempsyc')
      
      pkgs <- c('flextable', 'broom', 'report', 'effectsize')
      install_if_not_installed(pkgs)

   
    ## b) Modeling Packages
      
      # # load MASS (negative binomial model)
      #   library('MASS')

      # load performance
        library('performance')
      
      # load lmPerm
        library('lmPerm')
      
    ## d) Other packages
      # load here packages
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
  

      
###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  
  
  ### 2.1 Load RData
    ## a) Load RData tidy barn swallow data (parental care data)
      load(here('data/5_6_chr_node_df.RData'))
    
      
      
###############################################################################
##############   3. Age associations with reproductive success    #############
###############################################################################

  ### 3.1 Female age association with reproductive success
    ## a) female age association with total fecundity
      f_age_tot_fecund <- lm(total.fecundity ~ 
                                    Age.category,
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                          Sex == 'f' &
                                          !is.na(total.fecundity)))
    
    ## b) model summary
      summary(f_age_tot_fecund)    # model summary
      confint(f_age_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_age_tot_fecund)
      
    ## c) tidy effect estimates
      f_age_tot_fecund_tidy <- 
        tidy(f_age_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_age_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ age')
      
    ## d) female age association with nests 2 and 3fecundity
      f_age_nest_2_3_fecund <- lm(nest.2.3.fecund ~ 
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                       Sex == 'f' &
                                      !is.na(total.fecundity)))
      
    ## e) model summary
      summary(f_age_nest_2_3_fecund)    # model summary
      confint(f_age_nest_2_3_fecund, level = 0.95, method = 'profile')
      #plot(f_age_nest_2_3_fecund)
      
    # ## f) tidy effect estimates
    #   f_age_nest_2_3_fecund_tidy <- 
    #     tidy(f_age_nest_2_3_fecund, effects = 'fixed',
    #          conf.int = T,
    #          conf.method = 'quantile',
    #          conf.level = 0.95) %>%
    #     filter(term != '(Intercept)') %>%
    #     # rename(c(`mix-age` = V2,
    #     #          `asy` = V3)) %>%
    #     select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))
    #   # rownames_to_column(' ')
    #   
    #   # label stratified model
    #   f_age_nest_2_3_fecund_tidy$model <- c('nests 2-3 fecundity')
      
           
  ### 3.2 Male age association with reproductive success
    ## a) male age association with nests 2-3 totat paternity
      m_age_nest_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                      Age.category,
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(nest.2.3.tot.pat)))
     
    ## b) model summary
      summary(m_age_nest_2_3_tot_pat)    # model summary 
      confint(m_age_nest_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_age_nest_2_3_tot_pat)
      
    ## c) tidy effect estimates
      m_age_nest_2_3_tot_pat_tidy <- 
        tidy(m_age_nest_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_age_nest_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ age')
      
    ## d) male age association with nests 2-3 extra pair paternity
      m_age_nest_2_3_epp <- lm(nest.2.3.epp ~ 
                                    Age.category,
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                          Sex == 'm' &
                                          #Color.manipulation. == 'n' &
                                          !is.na(nest.2.3.epp)))
      
    ## e) model summary
      summary(m_age_nest_2_3_epp)    # model summary 
      confint(m_age_nest_2_3_epp, level = 0.95, method = 'profile')
      #plot(m_age_nest_2_3_epp)
      
    ## f) tidy effect estimates
      m_age_nest_2_3_epp_tidy <- 
        tidy(m_age_nest_2_3_epp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_age_nest_2_3_epp_tidy$model <- 
        c('epp ~ age')
      
    ## g) male age association with nests 2-3 social pair paternity
      m_age_nest_2_3_spp <- lm(nest.2.3.spp ~ 
                                    Age.category,
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                           Sex == 'm' &
                                           #Color.manipulation. == 'n' &
                                           !is.na(nest.2.3.spp)))
      
    ## h) model summary
      summary(m_age_nest_2_3_spp)    # model summary 
      confint(m_age_nest_2_3_spp, level = 0.95, method = 'profile')
      #plot(m_age_nest_2_3_spp)
      
    ## i) tidy effect estimates
      m_age_nest_2_3_spp_tidy <- 
        tidy(m_age_nest_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_age_nest_2_3_spp_tidy$model <- 
        c('spp ~ age')
     
      
  ### 3.3 Combine estimates and format for publication 
    ## a) Combine regression estimates into a tidy table
      age_detrmnts_repro_tidy <- bind_rows(f_age_tot_fecund_tidy,
                                           m_age_nest_2_3_tot_pat_tidy,
                                           m_age_nest_2_3_epp_tidy,
                                           m_age_nest_2_3_spp_tidy) 
      
      
    ## b) change order of columns
      age_detrmnts_repro_tidy <- age_detrmnts_repro_tidy %>%
        relocate(model, .before = term)
      
    ## c) format table creating flextable object 
      age_detrmnts_repro_tidy_frmt <- 
        nice_table(age_detrmnts_repro_tidy, separate.header = T)
      
    ## d) preview the table as a word doc    
      #print(age_detrmnts_repro_tidy_frmt, preview = 'docx')
      
    ## e) save age_detrmnts_repro_tidy_frmt as a word doc
      flextable::save_as_docx(age_detrmnts_repro_tidy_frmt, 
                path = here('output/age_detrmnts_repro_tidy_frmt.docx'))   
      
      
      
###############################################################################
############## 4. Plumage associations with reproductive success  #############
###############################################################################    

  ### 4.1 Female plumage trait associations with reproductive success
    ## a) female throat brightness association with 
      # total fecundity from all clutches
      f_T_bright_tot_fecund <- lm(total.fecundity ~ 
                                     scale(T_avg.bright),
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                                 Sex == 'f' &
                                                   !is.na(T_avg.bright) &
                                                   !is.na(total.fecundity)))
      
    ## b) model summary
      summary(f_T_bright_tot_fecund)    # model summary
      confint(f_T_bright_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_T_bright_tot_fecund)
      
    ## c) tidy effect estimates
      f_T_bright_tot_fecund_tidy <- 
        tidy(f_T_bright_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(T_avg.bright)', 
                              'throat brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_T_bright_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ T brightness')
      
    ## d) female breast brightness association with 
      # total fecundity from all clutches
      f_R_bright_tot_fecund <- lm(total.fecundity ~ 
                                    scale(R_avg.bright),
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                          Sex == 'f' &
                                          !is.na(R_avg.bright) &
                                          !is.na(total.fecundity)))
      
    ## e) model summary
      summary(f_R_bright_tot_fecund)    # model summary
      confint(f_R_bright_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_R_bright_tot_fecund)
      
    ## f) tidy effect estimates
      f_R_bright_tot_fecund_tidy <- 
        tidy(f_R_bright_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(R_avg.bright)', 
                              'breast brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_R_bright_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ R brightness')
      
    ## g) female belly brightness association with 
      # total fecundity from all clutches
      f_B_bright_tot_fecund <- lm(total.fecundity ~ 
                                        scale(B_avg.bright),
                                        #family = 'gaussian',
                                        data = subset(chr15_attrib_df,
                                                Sex == 'f' &
                                                !is.na(B_avg.bright) &
                                                !is.na(total.fecundity)))

      summary(f_B_bright_tot_fecund)    # model summary 
      confint(f_B_bright_tot_fecund, level = 0.95, method = 'profile') 
      #plot(f_B_bright_tot_fecund)       # check residuals
      
    ## h) tidy effect estimates
      f_B_bright_tot_fecund_tidy <- 
        tidy(f_B_bright_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_B_bright_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ B brightness')
      
    ## i) female tail streamer length association with 
      # total fecundity from all clutches
      f_TS_tot_fecund <- lm(total.fecundity ~ 
                                    scale(Mean.TS),
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                          Sex == 'f' &
                                          !is.na(Mean.TS) &
                                          !is.na(total.fecundity)))
      
      summary(f_TS_tot_fecund)    # model summary 
      confint(f_TS_tot_fecund, level = 0.95, method = 'profile') 
      #plot(f_TS_tot_fecund)       # check residuals
      
    ## j) tidy effect estimates
      f_TS_tot_fecund_tidy <- 
        tidy(f_TS_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_TS_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ TS')
  
      
  ### 4.2 Male plumage trait associations with total paternity (nest 2-3)
    ## a) male throat brightness association with 
      # total paternity (nests 2-3)
      m_T_bright_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                     scale(T_avg.bright),
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                                 Sex == 'm' &
                                                   !is.na(T_avg.bright) &
                                                   !is.na(nest.2.3.tot.pat)))
      
    ## b) model summary
      summary(m_T_bright_2_3_tot_pat)    # model summary
      confint(m_T_bright_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_T_bright_2_3_tot_pat)
      
    ## c) tidy effect estimates
      m_T_bright_2_3_tot_pat_tidy <- 
        tidy(m_T_bright_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(T_avg.bright)', 
                              'throat brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_T_bright_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ T brightness')
      
    ## d) male breast brightness association with 
      # total paternity (nests 2-3) 
      m_R_bright_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                  scale(R_avg.bright),
                                  # plumage manip. exper.
                                  #scale(R.bright.treat.and.orig), 
                                  #Color.manipulation.,
                                  #family = 'gaussian',
                                  data = subset(chr15_attrib_df,
                                        Sex == 'm' &
                                        !is.na(R_avg.bright) &
                                        # plumage manip. exper.
                                        # !is.na(R.bright.treat.and.orig) &
                                        !is.na(nest.2.3.tot.pat)))
      
    ## e) model summary
      summary(m_R_bright_2_3_tot_pat)    # model summary
      confint(m_R_bright_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_R_bright_2_3_tot_pat)
      
    ## f) tidy effect estimates
      m_R_bright_2_3_tot_pat_tidy <- 
        tidy(m_R_bright_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(R_avg.bright)', 
                              'breast brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_R_bright_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ R brightness')
      
    ## g) male belly brightness association with 
      # total paternity (nests 2-3) 
      m_B_bright_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                     scale(B_avg.bright),
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                                 Sex == 'm' &
                                                   !is.na(B_avg.bright) &
                                                   !is.na(nest.2.3.tot.pat)))
      
      summary(m_B_bright_2_3_tot_pat)    # model summary 
      confint(m_B_bright_2_3_tot_pat, level = 0.95, method = 'profile') 
      #plot(m_B_bright_2_3_tot_pat)       # check residuals
      
    ## h) tidy effect estimates
      m_B_bright_2_3_tot_pat_tidy <- 
        tidy(m_B_bright_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_B_bright_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ B brightness')
      
    ## i) male tail streamer length association with 
      # total paternity (nests 2-3) 
      m_TS_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                               scale(Mean.TS),
                             #family = 'gaussian',
                             data = subset(chr15_attrib_df,
                                           Sex == 'm' &
                                             !is.na(Mean.TS) &
                                             !is.na(nest.2.3.tot.pat)))
      
      summary(m_TS_2_3_tot_pat)    # model summary 
      confint(m_TS_2_3_tot_pat, level = 0.95, method = 'profile') 
      #plot(m_TS_2_3_tot_pat)       # check residuals
      
    ## j) tidy effect estimates
      m_TS_2_3_tot_pat_tidy <- 
        tidy(m_TS_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_TS_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ TS')
      
      
  ### 4.3 Male plumage trait associations with nests 2-3 extra pair paternity
    ## a) male throat brightness association with 
      # nests 2-3 extra pair paternity
      m_T_bright_2_3_epp <- lm(nest.2.3.epp ~ 
                                    scale(T_avg.bright),
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                            Sex == 'm' &
                                            !is.na(T_avg.bright) &
                                            !is.na(nest.2.3.epp)))
      
    ## b) model summary
      summary(m_T_bright_2_3_epp)    # model summary
      confint(m_T_bright_2_3_epp, level = 0.95, method = 'profile')
      #plot(m_T_bright_2_3_epp)
    
    ## c) tidy effect estimates
      m_T_bright_2_3_epp_tidy <- 
        tidy(m_T_bright_2_3_epp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(T_avg.bright)', 
                              'throat brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_T_bright_2_3_epp_tidy$model <- 
        c('epp ~ T brightness')
    
    ## d) male breast brightness association with 
      # nests 2-3 extra pair paternity
      m_R_bright_2_3_epp <- lm(nest.2.3.epp ~ 
                                scale(R_avg.bright),
                                # plumage manip. exper.
                                #scale(R.bright.treat.and.orig), 
                                #Color.manipulation., 
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                        Sex == 'm' &
                                        !is.na(R_avg.bright) &
                                        # plumage manip. exper.
                                        #!is.na(R.bright.treat.and.orig) &
                                        !is.na(nest.2.3.epp)))
    
    ## e) model summary
      summary(m_R_bright_2_3_epp)    # model summary
      confint(m_R_bright_2_3_epp, level = 0.95, method = 'profile')
      #plot(m_R_bright_2_3_epp)
    
    ## f) tidy effect estimates
      m_R_bright_2_3_epp_tidy <- 
      tidy(m_R_bright_2_3_epp, effects = 'fixed',
           conf.int = T,
           conf.method = 'quantile',
           conf.level = 0.95) %>%
      filter(term != '(Intercept)') %>%
      # rename(c(`mix-age` = V2,
      #          `asy` = V3)) %>%
      select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(R_avg.bright)', 
                              'breast brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_R_bright_2_3_epp_tidy$model <- 
        c('epp ~ R brightness')
    
    ## g) male belly brightness association with 
      # nests 2-3 extra pair paternity
      m_B_bright_2_3_epp <- lm(nest.2.3.epp ~ 
                                    scale(B_avg.bright),
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                          Sex == 'm' &
                                          !is.na(B_avg.bright) &
                                          !is.na(nest.2.3.epp)))
      
      summary(m_B_bright_2_3_epp)    # model summary 
      confint(m_B_bright_2_3_epp, level = 0.95, method = 'profile') 
      #plot(m_B_bright_2_3_epp)       # check residuals
    
    ## h) tidy effect estimates
      m_B_bright_2_3_epp_tidy <- 
        tidy(m_B_bright_2_3_epp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_B_bright_2_3_epp_tidy$model <- 
        c('epp ~ B brightness')
      
    ## i) male tail streamer length association with 
      # nests 2-3 extra pair paternity
      m_TS_2_3_epp <- lm(nest.2.3.epp ~ 
                              scale(Mean.TS),
                              #family = 'gaussian',
                              data = subset(chr15_attrib_df,
                                        Sex == 'm' &
                                        !is.na(Mean.TS) &
                                        !is.na(nest.2.3.epp)))
      
      summary(m_TS_2_3_epp)    # model summary 
      confint(m_TS_2_3_epp, level = 0.95, method = 'profile') 
      #plot(m_TS_2_3_tot_pat)       # check residuals
    
    ## j) tidy effect estimates
      m_TS_2_3_epp_tidy <- 
        tidy(m_TS_2_3_epp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_TS_2_3_epp_tidy$model <- 
        c('epp ~ TS')
    
      
  ### 4.4 Male plumage trait associations with nests 2-3 social pair paternity
    ## a) male throat brightness association with 
      # nests 2-3 social pair paternity
      m_T_bright_2_3_spp <- lm(nest.2.3.spp ~ 
                                scale(T_avg.bright),
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                        Sex == 'm' &
                                        !is.na(T_avg.bright) &
                                        !is.na(nest.2.3.spp)))
      
    ## b) model summary
      summary(m_T_bright_2_3_spp)    # model summary
      confint(m_T_bright_2_3_spp, level = 0.95, method = 'profile')
      #plot(m_T_bright_2_3_spp)
      
    ## c) tidy effect estimates
      m_T_bright_2_3_spp_tidy <- 
        tidy(m_T_bright_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(T_avg.bright)', 
                              'throat brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_T_bright_2_3_spp_tidy$model <- 
        c('spp ~ T brightness')
      
    ## d) male breast brightness association with 
      # nests 2-3 social pair paternity
      m_R_bright_2_3_spp <- lm(nest.2.3.spp ~ 
                            scale(R_avg.bright),
                            # plumage manip. exper.
                            #scale(R.bright.treat.and.orig), 
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                        Sex == 'm' &
                                        !is.na(R_avg.bright) &
                                        # plumage manip. exper.
                                        #!is.na(R.bright.treat.and.orig) &
                                        !is.na(nest.2.3.spp)))
      
    ## e) model summary
      summary(m_R_bright_2_3_spp)    # model summary
      confint(m_R_bright_2_3_spp, level = 0.95, method = 'profile')
      #plot(m_R_bright_2_3_spp)
      
    ## f) tidy effect estimates
      m_R_bright_2_3_spp_tidy <- 
        tidy(m_R_bright_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'scale(R_avg.bright)', 
                              'breast brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_R_bright_2_3_spp_tidy$model <- 
        c('spp ~ R brightness')
      
    ## g) male belly brightness association with 
      # nests 2-3 social pair paternity
      m_B_bright_2_3_spp <- lm(nest.2.3.spp ~ 
                            scale(B_avg.bright),
                            #family = 'gaussian',
                            data = subset(chr15_attrib_df,
                                    Sex == 'm' &
                                    !is.na(B_avg.bright) &
                                    !is.na(nest.2.3.spp)))
      
      summary(m_B_bright_2_3_spp)    # model summary 
      confint(m_B_bright_2_3_spp, level = 0.95, method = 'profile') 
      #plot(m_B_bright_2_3_spp)       # check residuals
      
    ## h) tidy effect estimates
      m_B_bright_2_3_spp_tidy <- 
        tidy(m_B_bright_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_B_bright_2_3_spp_tidy$model <- 
        c('spp ~ B brightness')
      
    ## i) male tail streamer length association with 
      # nests 2-3 social pair paternity
      m_TS_2_3_spp <- lm(nest.2.3.spp ~ 
                          scale(Mean.TS),
                          #family = 'gaussian',
                          data = subset(chr15_attrib_df,
                                  Sex == 'm' &
                                  !is.na(Mean.TS) &
                                  !is.na(nest.2.3.spp)))
      
      summary(m_TS_2_3_spp)    # model summary 
      confint(m_TS_2_3_spp, level = 0.95, method = 'profile') 
      #plot(m_TS_2_3_tot_pat)       # check residuals
      
    ## j) tidy effect estimates
      m_TS_2_3_spp_tidy <- 
        tidy(m_TS_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_TS_2_3_spp_tidy$model <- 
        c('spp ~ TS')

      
  ### 4.5 Female and male plumage trait associations with reproductive success,
      # independent of age - treating age as a confounder
      
  # Rationale - if a plumage trait is associated with reproductive success,
  #             then assess if estimate holds after age adjustment
      
    ## a) female tail streamer length association with 
      # total fecundity from all clutches - control for age
      f_TS_tot_fecund_indpn_age <- lm(total.fecundity ~ 
                              scale(Mean.TS)
                             + Age.category,,
                             #family = 'gaussian',
                             data = subset(chr15_attrib_df,
                                      Sex == 'f' &
                                      !is.na(Mean.TS) &
                                      !is.na(total.fecundity)))
      
      summary(f_TS_tot_fecund_indpn_age)    # model summary 
      confint(f_TS_tot_fecund_indpn_age, level = 0.95, method = 'profile') 
      #plot(f_TS_tot_fecund_indpn_age)       # check residuals
      
    ## b) tidy effect estimates
      f_TS_tot_fecund_indpn_age_tidy <- 
        tidy(f_TS_tot_fecund_indpn_age, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
        f_TS_tot_fecund_indpn_age_tidy$model <- 
        c('tot. fecund. ~ TS - age adjusted')
        
    ## c) male tail streamer length association with 
      # total paternity (nests 2-3) - control for age
      m_TS_2_3_tot_pat_indpn_age <- lm(nest.2.3.tot.pat ~ 
                              scale(Mean.TS)
                              + Age.category,
                              #family = 'gaussian',
                              data = subset(chr15_attrib_df,
                                      Sex == 'm' &
                                      !is.na(Mean.TS) &
                                      !is.na(nest.2.3.tot.pat)))
      
      summary(m_TS_2_3_tot_pat_indpn_age)    # model summary 
      confint(m_TS_2_3_tot_pat_indpn_age, level = 0.95, method = 'profile') 
      #plot(m_TS_2_3_tot_pat_indpn_age)       # check residuals
      
    ## d) tidy effect estimates
      m_TS_2_3_tot_pat_indpn_age_tidy <- 
        tidy(m_TS_2_3_tot_pat_indpn_age, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
        m_TS_2_3_tot_pat_indpn_age_tidy$model <- 
        c('tot. pat. ~ TS - age adjusted')
      
    ## e) male tail streamer length association with 
      # epp (nests 2-3) - control for age
      m_TS_2_3_epp_indpn_age <- lm(nest.2.3.epp ~ 
                                    scale(Mean.TS)
                                    + Age.category,
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                            Sex == 'm' &
                                            !is.na(Mean.TS) &
                                            !is.na(nest.2.3.epp)))
      
      summary(m_TS_2_3_epp_indpn_age)    # model summary 
      confint(m_TS_2_3_epp_indpn_age, level = 0.95, method = 'profile') 
      #plot(m_TS_2_3_epp_indpn_age)       # check residuals
      
    ## f) tidy effect estimates
      m_TS_2_3_epp_indpn_age_tidy <- 
        tidy(m_TS_2_3_epp_indpn_age, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_TS_2_3_epp_indpn_age_tidy$model <- 
        c('epp ~ TS - age adjusted')

    ## g) male belly brightness length association with 
      # spp (nests 2-3) - control for age
      m_B_bright_2_3_spp_indpn_age <- lm(nest.2.3.spp ~ 
                                scale(B_avg.bright)
                                + Age.category,
                              #family = 'gaussian',
                              data = subset(chr15_attrib_df,
                                      Sex == 'm' &
                                      !is.na(B_avg.bright) &
                                      !is.na(nest.2.3.epp)))
      
      summary(m_B_bright_2_3_spp_indpn_age)    # model summary 
      confint(m_B_bright_2_3_spp_indpn_age, level = 0.95, method = 'profile') 
      #plot(m_B_bright_2_3_spp_indpn_age)       # check residuals
      
    ## h) tidy effect estimates
      m_B_bright_2_3_spp_indpn_age_tidy <- 
        tidy(m_B_bright_2_3_spp_indpn_age, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_B_bright_2_3_spp_indpn_age_tidy$model <- 
        c('spp ~ B bright. - age adjusted')
      
      
  ### 4.6 Combine estimates and format for publication 
    ## a) Combine regression estimates into a tidy table
      plumage_detrmnts_repro_tidy <- bind_rows(f_T_bright_tot_fecund_tidy,
                                           f_R_bright_tot_fecund_tidy,
                                           f_B_bright_tot_fecund_tidy,
                                           f_TS_tot_fecund_tidy,
                                           f_TS_tot_fecund_indpn_age_tidy,
                                           m_T_bright_2_3_tot_pat_tidy,
                                           m_R_bright_2_3_tot_pat_tidy,
                                           m_B_bright_2_3_tot_pat_tidy,
                                           m_TS_2_3_tot_pat_tidy,
                                           m_TS_2_3_tot_pat_indpn_age_tidy,
                                           m_T_bright_2_3_epp_tidy,
                                           m_R_bright_2_3_epp_tidy,
                                           m_B_bright_2_3_epp_tidy,
                                           m_TS_2_3_epp_tidy,
                                           m_TS_2_3_epp_indpn_age_tidy,
                                           m_T_bright_2_3_spp_tidy,
                                           m_R_bright_2_3_spp_tidy,
                                           m_B_bright_2_3_spp_tidy,
                                           m_B_bright_2_3_spp_indpn_age_tidy,
                                           m_TS_2_3_spp_tidy) 
      

    ## b) change order of columns
      plumage_detrmnts_repro_tidy <- plumage_detrmnts_repro_tidy %>%
        relocate(model, .before = term)
      
    ## c) format table creating flextable object 
      plumage_detrmnts_repro_tidy_frmt <- 
        nice_table(plumage_detrmnts_repro_tidy, separate.header = T)
      
    ## d) preview the table as a word doc    
      #print(plumage_detrmnts_repro_tidy_frmt, preview = 'docx')
      
    ## e) save plumage_detrmnts_repro_tidy_frmt as a word doc
      flextable::save_as_docx(plumage_detrmnts_repro_tidy_frmt, 
                path = here('output/plumage_detrmnts_repro_tidy_frmt.docx'))   

 
    
###############################################################################
##############   5. Trait associations with soc. net. metrics    ##############
###############################################################################  
  
  # Rationale - Among age and plumage measures (traits) that are associated 
  #             with reproductive success, determine their association 
  #             with social network measures, which includes the following
      # 5.1 soc network ~ female age
      # 5.2 soc network ~ male age
      # 5.3 soc network ~ female TS
      # 5.4 soc network ~ male B bright
      # 5.5 soc network ~ male TS
  
  ### 5.1 Female age associations with node-level social network measures
    ## a) female age association with strength - female-male networks
      f_fxm_age_strength <- lm(strength.fxm ~ 
                                  Age.category,
                                  #family = 'gaussian',
                                  data = subset(chr15_attrib_df,
                                        Sex == 'f' &
                                        !is.na(strength.fxm)))
      
    ## b) model summary
      summary(f_fxm_age_strength)    # model summary 
      confint(f_fxm_age_strength, level = 0.95, method = 'profile')  
      #plot(f_fxm_age_strength)       # check residuals
      
    ## c) tidy effect estimates
      f_fxm_age_strength_tidy <- 
        tidy(f_fxm_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_fxm_age_strength_tidy$model <- 
        c('strength ~ female age - fxm')
      
      # use lmPerm to compute permutation p-values to check for issues w/
      # non-independence in social network measures
  ### Note. B. Bolker - https://mac-theobio.github.io/QMEE/lectures/permutation_examples.notes.html
    # lmp() seems to automatically change the contrast settings from the 
    # default treatment contrast to sum-to-zero contrasts, so that the reported 
    # effect size is half what it was (3.75/2), because it is computing the 
    # difference between the (unweighted) average of the two groups and the 
    # first group (field).
      # set.seed (1234)
      # fem.age.strength.pre.perm <- lmp(strength.pre ~ 
      #                               Age.category,
      #                             data = subset(chr15_attrib_df,
      #                                           Sex == 'f' &
      #                                             !is.na(strength.pre) &
      #                                             !is.na(total.fecundity)))
      # summary(fem.age.strength.pre.perm)    # model summary 
   
    ## d) female age association with degree - female-male networks
      f_fxm_age_degree <- lm(degree.fxm ~ 
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                        Sex == 'f' &
                                        !is.na(degree.fxm)))
      
    ## e) model summary
      summary(f_fxm_age_degree)    # model summary 
      confint(f_fxm_age_degree, level = 0.95, method = 'profile')  
      #plot(f_fxm_age_degree)       # check residuals
      
    ## f) tidy effect estimates
      f_fxm_age_degree_tidy <- 
        tidy(f_fxm_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_fxm_age_degree_tidy$model <- 
        c('degree ~ female age - fxm')
      
    ## g) female age association with strength - female-female networks
      f_fxf_age_strength <- lm(strength.fxf ~ 
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                      Sex == 'f' &
                                      !is.na(strength.fxf)))
      
    ## h) model summary
      summary(f_fxf_age_strength)    # model summary 
      confint(f_fxf_age_strength, level = 0.95, method = 'profile')  
      #plot(f_fxf_age_strength)       # check residuals
      
    ## i) tidy effect estimates
      f_fxf_age_strength_tidy <- 
        tidy(f_fxf_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_fxf_age_strength_tidy$model <- 
        c('strength ~ female age - fxf')
      
    ## j) female age association with degree - female-female networks
      f_fxf_age_degree <- lm(degree.fxf ~ 
                              Age.category,
                              #family = 'gaussian',
                              data = subset(chr15_attrib_df,
                                    Sex == 'f' &
                                    !is.na(degree.fxf)))
      
    ## k) model summary
      summary(f_fxf_age_degree)    # model summary 
      confint(f_fxf_age_degree, level = 0.95, method = 'profile')  
      #plot(f_fxf_age_degree)       # check residuals
      
    ## l) tidy effect estimates
      f_fxf_age_degree_tidy <- 
        tidy(f_fxf_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_fxf_age_degree_tidy$model <- 
        c('degree ~ female age - fxf')
      

  ### 5.2 Male age associations with node-level social network measures
    ## a) male age association with strength - female-male networks
      m_fxm_age_strength <- lm(strength.fxm ~ 
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                      Sex == 'm' &
                                      !is.na(strength.fxm)))
      
      ## b) model summary
      summary(m_fxm_age_strength)    # model summary 
      confint(m_fxm_age_strength, level = 0.95, method = 'profile')  
      #plot(m_fxm_age_strength)       # check residuals
      
      ## c) tidy effect estimates
      m_fxm_age_strength_tidy <- 
        tidy(m_fxm_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_fxm_age_strength_tidy$model <- 
        c('strength ~ male age - fxm')
      
    ## d) male age association with degree - female-male networks
      m_fxm_age_degree <- lm(degree.fxm ~ 
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                      Sex == 'm' &
                                      !is.na(degree.fxm)))
      
    ## e) model summary
      summary(m_fxm_age_degree)    # model summary 
      confint(m_fxm_age_degree, level = 0.95, method = 'profile')  
      #plot(m_fxm_age_degree)       # check residuals
      
    ## f) tidy effect estimates
      m_fxm_age_degree_tidy <- 
        tidy(m_fxm_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_fxm_age_degree_tidy$model <- 
        c('degree ~ male age - fxm')
     
    ## g) male age association with strength - male-male networks
      m_mxm_age_strength <- lm(strength.mxm ~ 
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                        Sex == 'm' &
                                        !is.na(strength.mxm)))
      
    ## h) model summary
      summary(m_mxm_age_strength)    # model summary 
      confint(m_mxm_age_strength, level = 0.95, method = 'profile')  
      #plot(m_mxm_age_strength)       # check residuals
      
    ## i) tidy effect estimates
      m_mxm_age_strength_tidy <- 
        tidy(m_mxm_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_mxm_age_strength_tidy$model <- 
        c('strength ~ male age - mxm')
      
    ## j) male age association with degree - male-male networks
      m_mxm_age_degree <- lm(degree.mxm ~ 
                                Age.category,
                              #family = 'gaussian',
                              data = subset(chr15_attrib_df,
                                      Sex == 'm' &
                                      !is.na(degree.mxm)))
      
    ## k) model summary
      summary(m_mxm_age_degree)    # model summary 
      confint(m_mxm_age_degree, level = 0.95, method = 'profile')  
      #plot(m_mxm_age_degree)       # check residuals
      
    ## l) tidy effect estimates
      m_mxm_age_degree_tidy <- 
        tidy(m_mxm_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_mxm_age_degree_tidy$model <- 
        c('degree ~ male age - mxm')
      

  ### 5.3 Female tail streamer length associations with node-level social 
      # network measures independent of age
    ## a) female TS association with strength - female-male networks
      f_fxm_TS_age_strength <- lm(strength.fxm ~ 
                                scale(Mean.TS) +
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                      Sex == 'f' &
                                      !is.na(Mean.TS) &
                                      !is.na(strength.fxm)))
      
      ## b) model summary
      summary(f_fxm_TS_age_strength)    # model summary 
      confint(f_fxm_TS_age_strength, level = 0.95, method = 'profile')  
      #plot(f_fxm_TS_age_strength)       # check residuals
      
    ## c) tidy effect estimates
      f_fxm_TS_age_strength_tidy <- 
        tidy(f_fxm_TS_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_fxm_TS_age_strength_tidy$model <- 
        c('strength ~ female TS - fxm; age adjusted')
      
    ## d) female TS association with degree - female-male networks
        f_fxm_TS_age_degree <- lm(degree.fxm ~ 
                                       scale(Mean.TS) +
                                       Age.category,
                                     #family = 'gaussian',
                                     data = subset(chr15_attrib_df,
                                            Sex == 'f' &
                                            !is.na(Mean.TS) &
                                            !is.na(degree.fxm)))
      
    ## e) model summary
      summary(f_fxm_TS_age_degree)    # model summary 
      confint(f_fxm_TS_age_degree, level = 0.95, method = 'profile')  
      #plot(f_fxm_TS_age_degree)       # check residuals
      
    ## f) tidy effect estimates
      f_fxm_TS_age_degree_tidy <- 
        tidy(f_fxm_TS_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_fxm_TS_age_degree_tidy$model <- 
        c('degree ~ female TS - fxm; age adjusted')
      
      
    ## g) female TS association with strength - female-female networks
      f_fxf_TS_age_strength <- lm(strength.fxf ~ 
                                     scale(Mean.TS) +
                                     Age.category,
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                           Sex == 'f' &
                                           !is.na(Mean.TS) &
                                           !is.na(strength.fxf)))
      
    ## h) model summary
      summary(f_fxf_TS_age_strength)    # model summary 
      confint(f_fxf_TS_age_strength, level = 0.95, method = 'profile')  
      #plot(f_fxf_TS_age_strength)       # check residuals
      
    ## i) tidy effect estimates
      f_fxf_TS_age_strength_tidy <- 
        tidy(f_fxf_TS_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_fxf_TS_age_strength_tidy$model <- 
        c('strength ~ female TS - fxf; age adjusted')
   
      
    ## j) female TS association with degree - female-female networks
      f_fxf_TS_age_degree <- lm(degree.fxf ~ 
                                   scale(Mean.TS) +
                                   Age.category,
                                 #family = 'gaussian',
                                 data = subset(chr15_attrib_df,
                                         Sex == 'f' &
                                         !is.na(Mean.TS) &
                                         !is.na(degree.fxf)))
      
    ## k) model summary
      summary(f_fxf_TS_age_degree)    # model summary 
      confint(f_fxf_TS_age_degree, level = 0.95, method = 'profile')  
      #plot(f_fxf_TS_age_degree)       # check residuals
      
    ## l) tidy effect estimates
      f_fxf_TS_age_degree_tidy <- 
        tidy(f_fxf_TS_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_fxf_TS_age_degree_tidy$model <- 
        c('degree ~ female TS - fxf; age adjusted')
      
  
  ### 5.4 Male belly brightness associations with node-level social 
      # network measures independent of age
    ## a) male B bright association with strength - female-male networks
      m_fxm_B_bright_age_strength <- lm(strength.fxm ~ 
                                     scale(B_avg.bright) +
                                     Age.category,
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                         Sex == 'm' &
                                         !is.na(B_avg.bright) &
                                         !is.na(strength.fxm)))
      
    ## b) model summary
      summary(m_fxm_B_bright_age_strength)    # model summary 
      confint(m_fxm_B_bright_age_strength, level = 0.95, method = 'profile')  
      #plot(m_fxm_B_bright_age_strength)       # check residuals
      
    ## c) tidy effect estimates
      m_fxm_B_bright_age_strength_tidy <- 
        tidy(m_fxm_B_bright_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness')) 
      # rownames_to_column(' ')
      
      # label stratified model
      m_fxm_B_bright_age_strength_tidy$model <- 
        c('strength ~ male B bright. - fxm; age adjusted')
      
    ## d) male B bright association with degree - female-male networks
      m_fxm_B_bright_age_degree <- lm(degree.fxm ~ 
                                   scale(B_avg.bright) +
                                   Age.category,
                                 #family = 'gaussian',
                                 data = subset(chr15_attrib_df,
                                       Sex == 'm' &
                                       !is.na(B_avg.bright) &
                                       !is.na(degree.fxm)))
      
    ## e) model summary
      summary(m_fxm_B_bright_age_degree)    # model summary 
      confint(m_fxm_B_bright_age_degree, level = 0.95, method = 'profile')  
      #plot(m_fxm_B_bright_age_degree)       # check residuals
      
    ## f) tidy effect estimates
      m_fxm_B_bright_age_degree_tidy <- 
        tidy(m_fxm_B_bright_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_fxm_B_bright_age_degree_tidy$model <- 
        c('degree ~ male B bright. - fxm; age adjusted')
      
    ## g) male B bright association with strength - male-male networks
      m_mxm_B_bright_age_strength <- lm(strength.mxm ~ 
                                     scale(B_avg.bright) +
                                     Age.category,
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                         Sex == 'm' &
                                         !is.na(B_avg.bright) &
                                         !is.na(strength.mxm)))
      
    ## h) model summary
      summary(m_mxm_B_bright_age_strength)    # model summary 
      confint(m_mxm_B_bright_age_strength, level = 0.95, method = 'profile')  
      #plot(m_mxm_B_bright_age_strength)       # check residuals
      
    ## i) tidy effect estimates
      m_mxm_B_bright_age_strength_tidy <- 
        tidy(m_mxm_B_bright_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_mxm_B_bright_age_strength_tidy$model <- 
        c('strength ~ male B bright. - mxm; age adjusted')
      
    ## j) male B bright association with degree - male-male networks
      m_mxm_B_bright_age_degree <- lm(degree.mxm ~ 
                                   scale(B_avg.bright) +
                                   Age.category,
                                 #family = 'gaussian',
                                 data = subset(chr15_attrib_df,
                                       Sex == 'm' &
                                       !is.na(B_avg.bright) &
                                       !is.na(degree.mxm)))
      
    ## k) model summary
      summary(m_mxm_B_bright_age_degree)    # model summary 
      confint(m_mxm_B_bright_age_degree, level = 0.95, method = 'profile')  
      #plot(m_mxm_B_bright_age_degree)       # check residuals
      
    ## l) tidy effect estimates
      m_mxm_B_bright_age_degree_tidy <- 
        tidy(m_mxm_B_bright_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(B_avg.bright)', 
                              'belly brightness'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_mxm_B_bright_age_degree_tidy$model <- 
        c('degree ~ male B bright. - mxm; age adjusted')

      
  ### 5.5 Male tail streamer length associations with node-level social 
      # network measures independent of age
   ## a) male TS association with strength - female-male networks
      m_fxm_TS_age_strength <- lm(strength.fxm ~ 
                                     scale(Mean.TS) +
                                     Age.category,
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                         Sex == 'm' &
                                         !is.na(Mean.TS) &
                                         !is.na(strength.fxm)))
      
    ## b) model summary
      summary(m_fxm_TS_age_strength)    # model summary 
      confint(m_fxm_TS_age_strength, level = 0.95, method = 'profile')  
      #plot(m_fxm_TS_age_strength)       # check residuals
      
    ## c) tidy effect estimates
      m_fxm_TS_age_strength_tidy <- 
        tidy(m_fxm_TS_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_fxm_TS_age_strength_tidy$model <- 
        c('strength ~ male TS - fxm; age adjusted')
      
    ## d) male TS association with degree - female-male networks
      m_fxm_TS_age_degree <- lm(degree.fxm ~ 
                                   scale(Mean.TS) +
                                   Age.category,
                                 #family = 'gaussian',
                                 data = subset(chr15_attrib_df,
                                       Sex == 'm' &
                                       !is.na(Mean.TS) &
                                       !is.na(degree.fxm)))
      
    ## e) model summary
      summary(m_fxm_TS_age_degree)    # model summary 
      confint(m_fxm_TS_age_degree, level = 0.95, method = 'profile')  
      #plot(m_fxm_TS_age_degree)       # check residuals
      
    ## f) tidy effect estimates
      m_fxm_TS_age_degree_tidy <- 
        tidy(m_fxm_TS_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_fxm_TS_age_degree_tidy$model <- 
        c('degree ~ male TS - fxm; age adjusted')
      
    ## g) male TS association with strength - male-male networks
      m_mxm_TS_age_strength <- lm(strength.mxm ~ 
                                     scale(Mean.TS) +
                                     Age.category,
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                         Sex == 'm' &
                                         !is.na(Mean.TS) &
                                         !is.na(strength.mxm)))
      
    ## h) model summary
      summary(m_mxm_TS_age_strength)    # model summary 
      confint(m_mxm_TS_age_strength, level = 0.95, method = 'profile')  
      #plot(m_mxm_TS_age_strength)       # check residuals
      
    ## i) tidy effect estimates
      m_mxm_TS_age_strength_tidy <- 
        tidy(m_mxm_TS_age_strength, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_mxm_TS_age_strength_tidy$model <- 
        c('strength ~ male TS - mxm; age adjusted')
      
    ## j) male TS association with degree - male-male networks
      m_mxm_TS_age_degree <- lm(degree.mxm ~ 
                                   scale(Mean.TS) +
                                   Age.category,
                                 #family = 'gaussian',
                                 data = subset(chr15_attrib_df,
                                       Sex == 'm' &
                                       !is.na(Mean.TS) &
                                       !is.na(degree.mxm)))
      
    ## k) model summary
      summary(m_mxm_TS_age_degree)    # model summary 
      confint(m_mxm_TS_age_degree, level = 0.95, method = 'profile')  
      #plot(m_mxm_TS_age_degree)       # check residuals
      
    ## l) tidy effect estimates
      m_mxm_TS_age_degree_tidy <- 
        tidy(m_mxm_TS_age_degree, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'scale(Mean.TS)', 
                              'tail streamer length'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_mxm_TS_age_degree_tidy$model <- 
        c('degree ~ male TS - mxm; age adjusted')
      
      
  ### 5.6 Combine estimates and format for publication 
    ## a) Combine regression estimates into a tidy table
      trait_detrmnts_socia_intx_tidy <- bind_rows(f_fxm_age_strength_tidy,
                                              f_fxm_age_degree_tidy,
                                              f_fxf_age_strength_tidy,
                                              f_fxf_age_degree_tidy,
                                              f_fxm_TS_age_strength_tidy,
                                              f_fxm_TS_age_degree_tidy,
                                              f_fxf_TS_age_strength_tidy,
                                              f_fxf_TS_age_degree_tidy,
                                              m_fxm_age_strength_tidy,
                                              m_fxm_age_degree_tidy,
                                              m_mxm_age_strength_tidy,
                                              m_mxm_age_degree_tidy,
                                              m_fxm_B_bright_age_strength_tidy,
                                              m_fxm_B_bright_age_degree_tidy,
                                              m_mxm_B_bright_age_strength_tidy,
                                              m_mxm_B_bright_age_degree_tidy,
                                              m_fxm_TS_age_strength_tidy,
                                              m_fxm_TS_age_degree_tidy,
                                              m_mxm_TS_age_strength_tidy,
                                              m_mxm_TS_age_degree_tidy) 
      
    ## b) change order of columns
      trait_detrmnts_socia_intx_tidy <- trait_detrmnts_socia_intx_tidy %>%
        relocate(model, .before = term)
      
    ## c) format table creating flextable object 
      trait_detrmnts_socia_intx_frmt <- 
        nice_table(trait_detrmnts_socia_intx_tidy, separate.header = T)
      
    ## d) preview the table as a word doc    
      #print(trait_detrmnts_socia_intx_frmt, preview = 'docx')
      
    ## e) save trait_detrmnts_socia_intx_frmt as a word doc
      flextable::save_as_docx(trait_detrmnts_socia_intx_frmt, 
                    path = here('output/trait_detrmnts_socia_intx_frmt.docx'))  
      
   
      
###############################################################################
############  6. Soc. net. metrics association with reproduction  #############
###############################################################################        
      
  ### 6.1 Female social network measures associated with reproductive success
    ## a) female strength association with total fecundity - fxm network
      f_strength_fxm_tot_fecund <- lm(total.fecundity ~ 
                                scale(strength.fxm),
                              #family = 'gaussian',
                              data = subset(chr15_attrib_df,
                                      Sex == 'f' &
                                      !is.na(strength.fxm) &
                                      !is.na(total.fecundity)))
      
    ## b) model summary
      summary(f_strength_fxm_tot_fecund)    # model summary
      confint(f_strength_fxm_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_strength_fxm_tot_fecund)
      
    ## c) tidy effect estimates
      f_strength_fxm_tot_fecund_tidy <- 
        tidy(f_strength_fxm_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(strength.fxm)', 
                              'strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_strength_fxm_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ strength - fxm')
 
    ## d) female degree association with total fecundity - fxm network
      f_degree_fxm_tot_fecund <- lm(total.fecundity ~ 
                                         scale(degree.fxm),
                                       #family = 'gaussian',
                                       data = subset(chr15_attrib_df,
                                             Sex == 'f' &
                                             !is.na(degree.fxm) &
                                             !is.na(total.fecundity)))
      
    ## e) model summary
      summary(f_degree_fxm_tot_fecund)    # model summary
      confint(f_degree_fxm_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_degree_fxm_tot_fecund)
      
    ## f) tidy effect estimates
      f_degree_fxm_tot_fecund_tidy <- 
        tidy(f_degree_fxm_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(degree.fxm)', 
                              'degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_degree_fxm_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ degree - fxm')

    ## g) female strength association with total fecundity - fxf network
      f_strength_fxf_tot_fecund <- lm(total.fecundity ~ 
                                         scale(strength.fxf),
                                       #family = 'gaussian',
                                       data = subset(chr15_attrib_df,
                                             Sex == 'f' &
                                             !is.na(strength.fxf) &
                                             !is.na(total.fecundity)))
      
    ## h) model summary
      summary(f_strength_fxf_tot_fecund)    # model summary
      confint(f_strength_fxf_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_strength_fxf_tot_fecund)
      
    ## i) tidy effect estimates
      f_strength_fxf_tot_fecund_tidy <- 
        tidy(f_strength_fxf_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(strength.fxf)', 
                              'strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_strength_fxf_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ strength - fxf')
      
    ## j) female degree association with total fecundity - fxf network
      f_degree_fxf_tot_fecund <- lm(total.fecundity ~ 
                                       scale(degree.fxf),
                                     #family = 'gaussian',
                                     data = subset(chr15_attrib_df,
                                           Sex == 'f' &
                                           !is.na(degree.fxf) &
                                           !is.na(total.fecundity)))
      
    ## k) model summary
      summary(f_degree_fxf_tot_fecund)    # model summary
      confint(f_degree_fxf_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_degree_fxf_tot_fecund)
      
    ## l) tidy effect estimates
      f_degree_fxf_tot_fecund_tidy <- 
        tidy(f_degree_fxf_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(degree.fxf)', 
                              'degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_degree_fxf_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ degree - fxf')
      
   
       
  ### 6.2 Male social network measure associations w/ total paternity (nest 2-3)
    ## a) male strength association with nests 2-3 total paternity - fxm network
      m_strength_fxm_nest_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                         scale(strength.fxm),
                                       #family = 'gaussian',
                                       data = subset(chr15_attrib_df,
                                             Sex == 'm' &
                                             !is.na(strength.fxm) &
                                             !is.na(nest.2.3.tot.pat)))
      
    ## b) model summary
      summary(m_strength_fxm_nest_2_3_tot_pat)    # model summary
      confint(m_strength_fxm_nest_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_strength_fxm_nest_2_3_tot_pat)
      
    ## c) tidy effect estimates
      m_strength_fxm_nest_2_3_tot_pat_tidy <- 
        tidy(m_strength_fxm_nest_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(strength.fxm)', 
                              'strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_strength_fxm_nest_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ strength - fxm')
      
    ## d) male degree association with nests 2-3 total paternity - fxm network
      m_degree_fxm_nest_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                       scale(degree.fxm),
                                     #family = 'gaussian',
                                     data = subset(chr15_attrib_df,
                                           Sex == 'm' &
                                           !is.na(degree.fxm) &
                                           !is.na(nest.2.3.tot.pat)))
      
    ## e) model summary
      summary(m_degree_fxm_nest_2_3_tot_pat)    # model summary
      confint(m_degree_fxm_nest_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_degree_fxm_nest_2_3_tot_pat)
      
    ## f) tidy effect estimates
      m_degree_fxm_nest_2_3_tot_pat_tidy <- 
        tidy(m_degree_fxm_nest_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(degree.fxm)', 
                              'degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_degree_fxm_nest_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ degree - fxm')
      
    ## g) male strength association with nests 2-3 total paternity - mxm network
      m_strength_mxm_nest_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                         scale(strength.mxm),
                                       #family = 'gaussian',
                                       data = subset(chr15_attrib_df,
                                             Sex == 'm' &
                                             !is.na(strength.mxm) &
                                             !is.na(nest.2.3.tot.pat)))
      
    ## h) model summary
      summary(m_strength_mxm_nest_2_3_tot_pat)    # model summary
      confint(m_strength_mxm_nest_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_strength_mxm_nest_2_3_tot_pat)
      
    ## i) tidy effect estimates
      m_strength_mxm_nest_2_3_tot_pat_tidy <- 
        tidy(m_strength_mxm_nest_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(strength.mxm)', 
                              'strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_strength_mxm_nest_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ strength - mxm')
      
    ## j) male degree association with nests 2-3 total paternity - mxm network
      m_degree_mxm_nest_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                       scale(degree.mxm),
                                     #family = 'gaussian',
                                     data = subset(chr15_attrib_df,
                                           Sex == 'm' &
                                           !is.na(degree.mxm) &
                                           !is.na(nest.2.3.tot.pat)))
      
      ## k) model summary
      summary(m_degree_mxm_nest_2_3_tot_pat)    # model summary
      confint(m_degree_mxm_nest_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_degree_mxm_nest_2_3_tot_pat)
      
    ## l) tidy effect estimates
      m_degree_mxm_nest_2_3_tot_pat_tidy <- 
        tidy(m_degree_mxm_nest_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(degree.mxm)', 
                              'degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_degree_mxm_nest_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ degree - mxm')
      
      
  ### 6.3 Male social network measure associations w/ extra pair 
      # paternity (nest 2-3)
    ## a) male strength association with nests 2-3 epp - fxm network
      m_strength_fxm_nest_2_3_epp <- lm(nest.2.3.epp ~ 
                                               scale(strength.fxm),
                                             #family = 'gaussian',
                                             data = subset(chr15_attrib_df,
                                                   Sex == 'm' &
                                                   !is.na(strength.fxm) &
                                                   !is.na(nest.2.3.epp)))
      
    ## b) model summary
      summary(m_strength_fxm_nest_2_3_epp)    # model summary
      confint(m_strength_fxm_nest_2_3_epp, level = 0.95, method = 'profile')
      #plot(m_strength_fxm_nest_2_3_epp)
      
    ## c) tidy effect estimates
      m_strength_fxm_nest_2_3_epp_tidy <- 
        tidy(m_strength_fxm_nest_2_3_epp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(strength.fxm)', 
                              'strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_strength_fxm_nest_2_3_epp_tidy$model <- 
        c('epp ~ strength - fxm')
      
    ## d) male degree association with nests 2-3 epp - fxm network
      m_degree_fxm_nest_2_3_epp <- lm(nest.2.3.epp ~ 
                                             scale(degree.fxm),
                                           #family = 'gaussian',
                                           data = subset(chr15_attrib_df,
                                                 Sex == 'm' &
                                                 !is.na(degree.fxm) &
                                                 !is.na(nest.2.3.epp)))
      
    ## e) model summary
      summary(m_degree_fxm_nest_2_3_epp)    # model summary
      confint(m_degree_fxm_nest_2_3_epp, level = 0.95, method = 'profile')
      #plot(m_degree_fxm_nest_2_3_epp)
      
    ## f) tidy effect estimates
      m_degree_fxm_nest_2_3_epp_tidy <- 
        tidy(m_degree_fxm_nest_2_3_epp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(degree.fxm)', 
                              'degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_degree_fxm_nest_2_3_epp_tidy$model <- 
        c('epp ~ degree - fxm')
      
    ## g) male strength association with nests 2-3 epp - mxm network
      m_strength_mxm_nest_2_3_epp <- lm(nest.2.3.epp ~ 
                                               scale(strength.mxm),
                                             #family = 'gaussian',
                                             data = subset(chr15_attrib_df,
                                                   Sex == 'm' &
                                                   !is.na(strength.mxm) &
                                                   !is.na(nest.2.3.epp)))
      
    ## h) model summary
      summary(m_strength_mxm_nest_2_3_epp)    # model summary
      confint(m_strength_mxm_nest_2_3_epp, level = 0.95, method = 'profile')
      #plot(m_strength_mxm_nest_2_3_epp)
      
    ## i) tidy effect estimates
      m_strength_mxm_nest_2_3_epp_tidy <- 
        tidy(m_strength_mxm_nest_2_3_epp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(strength.mxm)', 
                              'strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_strength_mxm_nest_2_3_epp_tidy$model <- 
        c('epp ~ strength - mxm')
      
    ## j) male degree association with nests 2-3 epp - mxm network
      m_degree_mxm_nest_2_3_epp <- lm(nest.2.3.epp ~ 
                                             scale(degree.mxm),
                                           #family = 'gaussian',
                                           data = subset(chr15_attrib_df,
                                                 Sex == 'm' &
                                                 !is.na(degree.mxm) &
                                                 !is.na(nest.2.3.epp)))
      
      ## k) model summary
      summary(m_degree_mxm_nest_2_3_epp)    # model summary
      confint(m_degree_mxm_nest_2_3_epp, level = 0.95, method = 'profile')
      #plot(m_degree_mxm_nest_2_3_epp)
      
    ## l) tidy effect estimates
      m_degree_mxm_nest_2_3_epp_tidy <- 
        tidy(m_degree_mxm_nest_2_3_epp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(degree.mxm)', 
                              'degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_degree_mxm_nest_2_3_epp_tidy$model <- 
        c('epp ~ degree - mxm')
      
     
  ### 6.4 Male social network measure associations w/ social pair 
      # paternity (nest 2-3)
    ## a) male strength association with nests 2-3 spp - fxm network
      m_strength_fxm_nest_2_3_spp <- lm(nest.2.3.spp ~ 
                                           scale(strength.fxm),
                                         #family = 'gaussian',
                                         data = subset(chr15_attrib_df,
                                               Sex == 'm' &
                                               !is.na(strength.fxm) &
                                               !is.na(nest.2.3.spp)))
      
   ## b) model summary
      summary(m_strength_fxm_nest_2_3_spp)    # model summary
      confint(m_strength_fxm_nest_2_3_spp, level = 0.95, method = 'profile')
      #plot(m_strength_fxm_nest_2_3_spp)
      
    ## c) tidy effect estimates
      m_strength_fxm_nest_2_3_spp_tidy <- 
        tidy(m_strength_fxm_nest_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(strength.fxm)', 
                              'strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_strength_fxm_nest_2_3_spp_tidy$model <- 
        c('spp ~ strength - fxm')
      
    ## d) male degree association with nests 2-3 spp - fxm network
      m_degree_fxm_nest_2_3_spp <- lm(nest.2.3.spp ~ 
                                         scale(degree.fxm),
                                       #family = 'gaussian',
                                       data = subset(chr15_attrib_df,
                                               Sex == 'm' &
                                               !is.na(degree.fxm) &
                                               !is.na(nest.2.3.spp)))
      
    ## e) model summary
      summary(m_degree_fxm_nest_2_3_spp)    # model summary
      confint(m_degree_fxm_nest_2_3_spp, level = 0.95, method = 'profile')
      #plot(m_degree_fxm_nest_2_3_spp)
      
    ## f) tidy effect estimates
      m_degree_fxm_nest_2_3_spp_tidy <- 
        tidy(m_degree_fxm_nest_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(degree.fxm)', 
                              'degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_degree_fxm_nest_2_3_spp_tidy$model <- 
        c('spp ~ degree - fxm')
      
    ## g) male strength association with nests 2-3 spp - mxm network
      m_strength_mxm_nest_2_3_spp <- lm(nest.2.3.spp ~ 
                                           scale(strength.mxm),
                                         #family = 'gaussian',
                                         data = subset(chr15_attrib_df,
                                                 Sex == 'm' &
                                                 !is.na(strength.mxm) &
                                                 !is.na(nest.2.3.spp)))
      
    ## h) model summary
      summary(m_strength_mxm_nest_2_3_spp)    # model summary
      confint(m_strength_mxm_nest_2_3_spp, level = 0.95, method = 'profile')
      #plot(m_strength_mxm_nest_2_3_spp)
      
    ## i) tidy effect estimates
      m_strength_mxm_nest_2_3_spp_tidy <- 
        tidy(m_strength_mxm_nest_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(strength.mxm)', 
                              'strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_strength_mxm_nest_2_3_spp_tidy$model <- 
        c('spp ~ strength - mxm')
      
    ## j) male degree association with nests 2-3 spp - mxm network
      m_degree_mxm_nest_2_3_spp <- lm(nest.2.3.spp ~ 
                                         scale(degree.mxm),
                                       #family = 'gaussian',
                                       data = subset(chr15_attrib_df,
                                               Sex == 'm' &
                                               !is.na(degree.mxm) &
                                               !is.na(nest.2.3.spp)))
      
    ## k) model summary
      summary(m_degree_mxm_nest_2_3_spp)    # model summary
      confint(m_degree_mxm_nest_2_3_spp, level = 0.95, method = 'profile')
      #plot(m_degree_mxm_nest_2_3_spp)
      
    ## l) tidy effect estimates
      m_degree_mxm_nest_2_3_spp_tidy <- 
        tidy(m_degree_mxm_nest_2_3_spp, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(degree.mxm)', 
                              'degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_degree_mxm_nest_2_3_spp_tidy$model <- 
        c('spp ~ degree - mxm')
      
      
  ### 6.5 Combine estimates and format for publication 
    ## a) Combine regression estimates into a tidy table
      social_detrmnts_repro_tidy <- bind_rows(f_strength_fxm_tot_fecund_tidy,
                                          f_degree_fxm_tot_fecund_tidy,
                                          f_strength_fxf_tot_fecund_tidy,
                                          f_degree_fxf_tot_fecund_tidy,
                                          m_strength_fxm_nest_2_3_tot_pat_tidy,
                                          m_degree_fxm_nest_2_3_tot_pat_tidy,
                                          m_strength_mxm_nest_2_3_tot_pat_tidy,
                                          m_degree_mxm_nest_2_3_tot_pat_tidy,
                                          m_strength_fxm_nest_2_3_epp_tidy,
                                          m_degree_fxm_nest_2_3_epp_tidy,
                                          m_strength_mxm_nest_2_3_epp_tidy,
                                          m_degree_mxm_nest_2_3_epp_tidy,
                                          m_strength_fxm_nest_2_3_spp_tidy,
                                          m_degree_fxm_nest_2_3_spp_tidy,
                                          m_strength_mxm_nest_2_3_spp_tidy,
                                          m_degree_mxm_nest_2_3_spp_tidy) 
      
    ## b) change order of columns
      social_detrmnts_repro_tidy <- social_detrmnts_repro_tidy %>%
        relocate(model, .before = term)
      
    ## c) format table creating flextable object 
      social_detrmnts_repro_frmt <- 
        nice_table(social_detrmnts_repro_tidy, separate.header = T)
      
    ## d) preview the table as a word doc    
      #print(social_detrmnts_repro_frmt, preview = 'docx')
      
    ## e) save social_detrmnts_repro_frmt as a word doc
      flextable::save_as_docx(social_detrmnts_repro_frmt, 
                    path = here('output/social_detrmnts_repro_frmt.docx'))  
      
    
     
      
###############################################################################
##############      7. Age by soc. net. metrics interaction      ##############
###############################################################################  
      
      
  # Rationale - Among age and plumage measures (traits) that are associated 
  #             with reproductive success in at least one sex, determine
  #             if their are exposure mediation interactions. For mediators,
  #             limit where age or plumage is associated with soc. network 
  #             in at least one sex
      # 7.1 age * strength(fxm)
      # 7.2 TS * strength(fxm)
      # 7.3 TS * degree(fxm)
      
  ### 7.1 Interaction models: Female trait by node level social network 
        # measures with reproductive success
    ## a) Interaction model: female age by strength with total fecundity
      f_age_x_strength_fxm_tot_fecund <- lm(total.fecundity ~ 
                                     Age.category * scale(strength.fxm),
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                          Sex == 'f' &
                                          !is.na(strength.fxm) &
                                          !is.na(total.fecundity)))
      
    ## b) model summary
      summary(f_age_x_strength_fxm_tot_fecund)    # model summary
      confint(f_age_x_strength_fxm_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_age_x_strength_fxm_tot_fecund)
      
    ## c) tidy effect estimates
      f_age_x_strength_fxm_tot_fecund_tidy <- 
        tidy(f_age_x_strength_fxm_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'Age.categoryasy:scale(strength.fxm)', 
                              'age * strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_age_x_strength_fxm_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ age * strength')
      
    ## d) Interaction model: female tail streamer length by strength 
      # with total fecundity
      f_TS_x_strength_fxm_tot_fecund <- lm(total.fecundity ~ 
                                         scale(Mean.TS) * scale(strength.fxm),
                                        #family = 'gaussian',
                                        data = subset(chr15_attrib_df,
                                             Sex == 'f' &
                                             !is.na(Mean.TS) &
                                             !is.na(total.fecundity)))
      
    ## e) model summary
      summary(f_TS_x_strength_fxm_tot_fecund)    # model summary
      confint(f_TS_x_strength_fxm_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_TS_x_strength_fxm_tot_fecund)
      
    ## f) tidy effect estimates
      f_TS_x_strength_fxm_tot_fecund_tidy <- 
        tidy(f_TS_x_strength_fxm_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(Mean.TS):scale(strength.fxm)', 
                              'TS * strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_TS_x_strength_fxm_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ TS * strength')
      
    ## g) Interaction model: female tail streamer length by degree 
      # with total fecundity
      f_TS_x_degree_fxm_tot_fecund <- lm(total.fecundity ~ 
                                      scale(Mean.TS) * scale(degree.fxm),
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                          Sex == 'f' &
                                          !is.na(Mean.TS) &
                                          !is.na(total.fecundity)))
      
    ## h) model summary
      summary(f_TS_x_degree_fxm_tot_fecund)    # model summary
      confint(f_TS_x_degree_fxm_tot_fecund, level = 0.95, method = 'profile')
      #plot(f_TS_x_degree_fxm_tot_fecund)
      
    ## i) tidy effect estimates
      f_TS_x_degree_fxm_tot_fecund_tidy <- 
        tidy(f_TS_x_degree_fxm_tot_fecund, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(Mean.TS):scale(degree.fxm)', 
                              'TS * degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_TS_x_degree_fxm_tot_fecund_tidy$model <- 
        c('tot. fecund. ~ TS * degree')
      

  ### 7.2 Interaction models: Male age by node level social network measures
      # associations with reproductive success
    ## a) Interaction model: male age by strength with total fecundity
      m_age_x_strength_fxm_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                       Age.category * scale(strength.fxm),
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_df,
                                             Sex == 'm' &
                                             !is.na(strength.fxm) &
                                             !is.na(nest.2.3.tot.pat)))
      
    ## b) model summary
      summary(m_age_x_strength_fxm_2_3_tot_pat)    # model summary
      confint(m_age_x_strength_fxm_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_age_x_strength_fxm_2_3_tot_pat)
      
    ## c) tidy effect estimates
      m_age_x_strength_fxm_2_3_tot_pat_tidy <- 
        tidy(m_age_x_strength_fxm_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'Age.categoryasy:scale(strength.fxm)', 
                              'age * strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_age_x_strength_fxm_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ age * strength')
      
    ## d) Interaction model: female tail streamer length by strength 
      # with total fecundity
      m_TS_x_strength_fxm_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                      scale(Mean.TS) * scale(strength.fxm),
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                          Sex == 'm' &
                                          !is.na(Mean.TS) &
                                          !is.na(nest.2.3.tot.pat)))
      
    ## e) model summary
      summary(m_TS_x_strength_fxm_2_3_tot_pat)    # model summary
      confint(m_TS_x_strength_fxm_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_TS_x_strength_fxm_2_3_tot_pat)
      
    ## f) tidy effect estimates
      m_TS_x_strength_fxm_2_3_tot_pat_tidy <- 
        tidy(m_TS_x_strength_fxm_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(Mean.TS):scale(strength.fxm)', 
                              'TS * strength'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_TS_x_strength_fxm_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ TS * strength')
      
    ## g) Interaction model: female tail streamer length by degree 
      # with total fecundity
      m_TS_x_degree_fxm_2_3_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                            scale(Mean.TS) * scale(degree.fxm),
                                          #family = 'gaussian',
                                          data = subset(chr15_attrib_df,
                                                        Sex == 'm' &
                                                          !is.na(Mean.TS) &
                                                          !is.na(nest.2.3.tot.pat)))
      
    ## h) model summary
      summary(m_TS_x_degree_fxm_2_3_tot_pat)    # model summary
      confint(m_TS_x_degree_fxm_2_3_tot_pat, level = 0.95, method = 'profile')
      #plot(m_TS_x_degree_fxm_2_3_tot_pat)
      
    ## i) tidy effect estimates
      m_TS_x_degree_fxm_2_3_tot_pat_tidy <- 
        tidy(m_TS_x_degree_fxm_2_3_tot_pat, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high'))%>%
        mutate(term = replace(term, term == 'scale(Mean.TS):scale(degree.fxm)', 
                              'TS * degree'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_TS_x_degree_fxm_2_3_tot_pat_tidy$model <- 
        c('tot. pat. ~ TS * degree')
      
  
  ### 6.5 Combine estimates and format for publication 
    ## a) Combine regression estimates into a tidy table
      exp_med_intx_models_tidy <- bind_rows(f_age_x_strength_fxm_tot_fecund_tidy,
                                    f_TS_x_strength_fxm_tot_fecund_tidy,
                                    f_TS_x_degree_fxm_tot_fecund_tidy,
                                    m_age_x_strength_fxm_2_3_tot_pat_tidy,
                                    m_TS_x_strength_fxm_2_3_tot_pat_tidy,                                
                                    m_TS_x_degree_fxm_2_3_tot_pat_tidy) 
      
    ## b) change order of columns
      exp_med_intx_models_tidy <- exp_med_intx_models_tidy %>%
        relocate(model, .before = term)
      
    ## c) format table creating flextable object 
      exp_med_intx_models_frmt <- 
        nice_table(exp_med_intx_models_tidy, separate.header = T)
      
    ## d) preview the table as a word doc    
      #print(exp_med_intx_models_frmt, preview = 'docx')
      
    ## e) save exp_med_intx_models_frmt as a word doc
      flextable::save_as_docx(exp_med_intx_models_frmt, 
                        path = here('output/exp_med_intx_models_frmt.docx'))  
      
      
      
 
###############################################################################
##############                  8. Export data                   ##############
###############################################################################  

      
  ### 8.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 mediation analysis
      save(file = here('data/7_chr_effect_decomp_data.RData'), 
           list = c('chr15_attrib_df'))
      
    ## b) Save and export model objects for 2015 mediation analysis 
      # from section 5
      # age --> strength(fxm) 
      # TS --> strength(fxm)  - age adjusted
      # TS --> degree(fxm)    - age adjusted
      save(file = here('data/7_chr_effect_decomp_mods.RData'), 
           list = c('f_fxm_age_strength', 'f_fxm_TS_age_strength',
                    'f_fxm_TS_age_degree', 
                    'm_fxm_age_strength', 'm_fxm_TS_age_strength',
                    'm_fxm_TS_age_degree'))
      