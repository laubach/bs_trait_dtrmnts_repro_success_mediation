################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############    6.5 Age-plumage and Med-Mod interaction models    #############
#############               Supplementary Materials                #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 29 Sept 2024                #############
#############               last updated: 14 Oct 2025              #############
################################################################################



### PURPOSE: Run mediation and total effects models for traits and reproductive
           # success
         
  # Code Blocks
    # 1. Configure work space
    # 2. Load RData
    # 3. Age association with plumage traits 
    # 4. Age and TS by soc. net. metrics interaction

    

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
##############      3. Age associations with plumage traits      ##############
###############################################################################  
    
  ### 3.1 Female age associations with plumage traits    
    ## a) Female age association with throat brightness
      f_age_T_bright <- lm(T_avg.bright ~ Age.category,
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                           Sex == 'f' &
                                           !is.na(Age.category) &
                                           !is.na(T_avg.bright)))
      
      summary(f_age_T_bright)    # model summary 
      confint(f_age_T_bright, level = 0.95, method = 'profile') 
      #plot(f_age_T_bright)       # check residuals
      
    ## b) tidy effect estimates
      f_age_T_bright_tidy <- 
        tidy(f_age_T_bright, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        #filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
    # rownames_to_column(' ')
    
    # label stratified model
      f_age_T_bright_tidy$model <- 
        c('T_bright ~ age - female')
    
    ## c) Female age association with breast brightness
      f_age_R_bright <- lm(R_avg.bright ~ Age.category,
                           #family = 'gaussian',
                           data = subset(chr15_attrib_df,
                                         Sex == 'f' &
                                           !is.na(Age.category) &
                                           !is.na(R_avg.bright)))
      
      summary(f_age_R_bright)    # model summary 
      confint(f_age_R_bright, level = 0.95, method = 'profile') 
      #plot(f_age_R_bright)       # check residuals
      
    ## d) tidy effect estimates
      f_age_R_bright_tidy <- 
        tidy(f_age_R_bright, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        #filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_age_R_bright_tidy$model <- 
        c('R_bright ~ age - female')
      
    ## e) Female age association with belly brightness
      f_age_B_bright <- lm(B_avg.bright ~ Age.category,
                           #family = 'gaussian',
                           data = subset(chr15_attrib_df,
                                         Sex == 'f' &
                                           !is.na(Age.category) &
                                           !is.na(B_avg.bright)))
      
      summary(f_age_B_bright)    # model summary 
      confint(f_age_B_bright, level = 0.95, method = 'profile') 
      #plot(f_age_B_bright)       # check residuals
      
    ## f) tidy effect estimates
      f_age_B_bright_tidy <- 
        tidy(f_age_B_bright, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        #filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_age_B_bright_tidy$model <- 
        c('B_bright ~ age - female')
    
    ## g) female age association with tail streamer length
      f_age_TS <- lm(Mean.TS ~ Age.category,
                     #family = 'gaussian',
                     data = subset(chr15_attrib_df,
                                   Sex == 'f' &
                                     !is.na(Age.category) &
                                     !is.na(Mean.TS)))
      
      summary(f_age_TS)    # model summary 
      confint(f_age_TS, level = 0.95, method = 'profile') 
      #plot(f_age_TS)       # check residuals
      
    ## h) tidy effect estimates
      f_age_TS_tidy <- 
        tidy(f_age_TS, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        #filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      f_age_TS_tidy$model <- 
        c('TS ~ age - female')
    
      
  ### 3.2. Male age associations with plumage traits    
    ## a) Male age association with throat brightness
      m_age_T_bright <- lm(T_avg.bright ~ Age.category,
                           #family = 'gaussian',
                           data = subset(chr15_attrib_df,
                                         Sex == 'm' &
                                           !is.na(Age.category) &
                                           !is.na(T_avg.bright)))
      
      summary(m_age_T_bright)    # model summary 
      confint(m_age_T_bright, level = 0.95, method = 'profile') 
      #plot(m_age_T_bright)       # check residuals
      
    ## b) tidy effect estimates
      m_age_T_bright_tidy <- 
        tidy(m_age_T_bright, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        #filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_age_T_bright_tidy$model <- 
        c('T_bright ~ age - male')
      
    ## c) Male age association with breast brightness
      m_age_R_bright <- lm(R_avg.bright ~ Age.category,
                           #family = 'gaussian',
                           data = subset(chr15_attrib_df,
                                         Sex == 'm' &
                                           !is.na(Age.category) &
                                           !is.na(R_avg.bright)))
      
      summary(m_age_R_bright)    # model summary 
      confint(m_age_R_bright, level = 0.95, method = 'profile') 
      #plot(m_age_R_bright)       # check residuals
      
    ## d) tidy effect estimates
      m_age_R_bright_tidy <- 
        tidy(m_age_R_bright, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        #filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_age_R_bright_tidy$model <- 
        c('R_bright ~ age - male')
      
    ## e) Male age association with belly brightness
      m_age_B_bright <- lm(B_avg.bright ~ Age.category,
                           #family = 'gaussian',
                           data = subset(chr15_attrib_df,
                                         Sex == 'm' &
                                           !is.na(Age.category) &
                                           !is.na(B_avg.bright)))
      
      summary(m_age_B_bright)    # model summary 
      confint(m_age_B_bright, level = 0.95, method = 'profile') 
      #plot(m_age_B_bright)       # check residuals
      
    ## f) tidy effect estimates
      m_age_B_bright_tidy <- 
        tidy(m_age_B_bright, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        #filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_age_B_bright_tidy$model <- 
        c('B_bright ~ age - male')
      
    ## g) Male age association with tail streamer length
      m_age_TS <- lm(Mean.TS ~ Age.category,
                     #family = 'gaussian',
                     data = subset(chr15_attrib_df,
                                   Sex == 'm' &
                                     !is.na(Age.category) &
                                     !is.na(Mean.TS)))
      
      summary(m_age_TS)    # model summary 
      confint(m_age_TS, level = 0.95, method = 'profile') 
      #plot(m_age_TS)       # check residuals
      
    ## h) tidy effect estimates
      m_age_TS_tidy <- 
        tidy(m_age_TS, effects = 'fixed',
             conf.int = T,
             conf.method = 'quantile',
             conf.level = 0.95) %>%
        filter(term != '(Intercept)') %>%
        # rename(c(`mix-age` = V2,
        #          `asy` = V3)) %>%
        select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
        #filter(term != 'Age.categoryasy') %>%
        mutate(term = replace(term, term == 'Age.categoryasy', 
                              'age (ASY)'))
      # rownames_to_column(' ')
      
      # label stratified model
      m_age_TS_tidy$model <- 
        c('TS ~ age - male')  
    
  ### 5.6 Combine estimates and format for publication 
    ## a) Combine regression estimates into a tidy table
      age_trait_assoc_tidy <- bind_rows(f_age_T_bright_tidy,
                                                  f_age_R_bright_tidy,
                                                  f_age_B_bright_tidy,
                                                  f_age_TS_tidy,
                                                  m_age_T_bright_tidy,
                                                  m_age_R_bright_tidy,
                                                  m_age_B_bright_tidy,
                                                  m_age_TS_tidy) 
      
    ## b) change order of columns
      age_trait_assoc_tidy <- age_trait_assoc_tidy %>%
        relocate(model, .before = term)
    
    ## c) format table creating flextable object 
      age_trait_assoc_frmt <- 
        nice_table(age_trait_assoc_tidy, separate.header = T)
    
    ## d) preview the table as a word doc    
      #print(age_trait_assoc_frmt, preview = 'docx')
    
    ## e) save age_trait_assoc_frmt as a word doc
      flextable::save_as_docx(age_trait_assoc_frmt, 
                              path = here('output/age_trait_assoc_frmt.docx'))  
      
    
    
###############################################################################
##############  4. Age and TS by soc. net. metrics interaction   ##############
###############################################################################  
      
      
  # Rationale - Among age and plumage measures (traits) that are associated 
  #             with reproductive success, determine if their are exposure
  #             mediator interactions. For mediators, limit where age or 
  #             plumage is associated with soc. network in at least one sex
    # 4.1 age * strength(fxm)
    # 4.2 TS * strength(fxm)
      
  ### 4.1 Interaction models: Female trait by node level social network 
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
      
  
      
  ### 4.2 Combine estimates and format for publication 
    ## a) Combine regression estimates into a tidy table
      exp_med_intx_models_tidy <- bind_rows(f_age_x_strength_fxm_tot_fecund_tidy,
                                    f_TS_x_strength_fxm_tot_fecund_tidy) 
      
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

      
