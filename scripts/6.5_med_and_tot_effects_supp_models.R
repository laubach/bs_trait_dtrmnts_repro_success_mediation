################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############        6.5 Mediation and total effects models        #############
#############               Supplementary Materials                #############
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
    # 3. Trait associations with soc. net. metrics 
    # 4. Age association with plumage traits 

    

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
##############   3. Trait associations with soc. net. metrics    ##############
###############################################################################  
  
is age associated with traits - focus on traits that are associated with fitness
    
    
  ## a) female age association with tail streamer length
    f_age_TS <- lm(Mean.TS ~ Age.category,,
                                     #family = 'gaussian',
                                     data = subset(chr15_attrib_df,
                                                   Sex == 'f' &
                                                     !is.na(Age.category) &
                                                     !is.na(Mean.TS)))
    
    summary(f_age_TS)    # model summary 
    confint(f_age_TS, level = 0.95, method = 'profile') 
    #plot(f_age_TS)       # check residuals
    
    ## b) tidy effect estimates
    f_age_TS_tidy <- 
      tidy(f_age_TS, effects = 'fixed',
           conf.int = T,
           conf.method = 'quantile',
           conf.level = 0.95) %>%
      filter(term != '(Intercept)') %>%
      # rename(c(`mix-age` = V2,
      #          `asy` = V3)) %>%
      select(c('term', 'estimate', 'p.value', 'conf.low', 'conf.high')) %>%
      filter(term != 'Age.categoryasy') %>%
      mutate(term = replace(term, term == 'Age.categoryasy', 
                            'age (ASY)'))
    # rownames_to_column(' ')
    
    # label stratified model
    f_age_TS_tidy$model <- 
      c('TS ~ age - female')
    
#*****************************************************************************#     

    # NEED TO ADD sensitivity looking at trait associ with post network measures
    # and post network measures with reproductive success
    
    
    ## e) female post manip strength association with 
    # total fecundity
    fem.strength.post.tot.fecund <- lm(total.fecundity ~ 
                                          scale(strength.post)
                                        # sensitivity, is strength associated with
                                        # with fecundity independent of age
                                        #+ Age.category
                                        ,
                                        #family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                                      Sex == 'f' &
                                                        !is.na(strength.post) &
                                                        !is.na(total.fecundity)))
    
    
    ## f) model summary
    summary(fem.strength.post.tot.fecund)    # model summary 
    confint(fem.strength.post.tot.fecund, 
            level = 0.95, method = 'profile')
    #plot(fem.strength.post.tot.fecund)       # check residuals
    # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
    #check_overdispersion(fem.strength.post.tot.fecund) 
    
    ## g) female post manip strength association with 
    # total fecundity
    fem.deg.post.tot.fecund <- lm(total.fecundity ~ 
                                     scale(degree.post)
                                   # sensitivity, is strength associated with
                                   # with fecundity independent of age
                                   #+ Age.category
                                   ,
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_post_df,
                                                 Sex == 'f' &
                                                   !is.na(degree.post) &
                                                   !is.na(total.fecundity)))
    
    ## h) model summary
    summary(fem.deg.post.tot.fecund)    # model summary 
    confint(fem.deg.post.tot.fecund, 
            level = 0.95, method = 'profile')
    #plot(fem.deg.post.tot.fecund)       # check residuals
    # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
    #check_overdispersion(fem.deg.post.tot.fecund) 
    
    
    
    
    
    ## d) male post manip strength association with 
    # clutch 2 total paternity
    m.strength.post.clutch.2.pat <- lm(attmpt.2.tot.pat ~ 
                                          scale(strength.post)
                                        # sensitivity, is strength associated with
                                        # with fecundity independent of age
                                        #+ Age.category
                                        ,
                                        #family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                                      Sex == 'm' &
                                                        !is.na(strength.post) &
                                                        !is.na(nest.2.3.tot.pat)))
    
    
    ## f) model summary
    summary(m.strength.post.clutch.2.pat)    # model summary 
    confint(m.strength.post.clutch.2.pat, 
            level = 0.95, method = 'profile')
    #plot(m.strength.post.clutch.2.pat)       # check residuals
    # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
    #check_overdispersion(m.strength.post.clutch.2.pat) 
    
    ## g) male post manip degree association with 
    # clutch 2 total paternity
    m.deg.post.clutch.2.pat <- lm(attmpt.2.tot.pat ~ 
                                     scale(degree.post)
                                   # sensitivity, is strength associated with
                                   # with fecundity independent of age
                                   #+ Age.category
                                   ,
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_post_df,
                                                 Sex == 'm' &
                                                   !is.na(degree.post) &
                                                   !is.na(nest.2.3.tot.pat)))
    
    ## h) model summary
    summary(m.deg.post.clutch.2.pat)    # model summary 
    confint(m.deg.post.clutch.2.pat, 
            level = 0.95, method = 'profile')
    #plot(m.deg.post.clutch.2.pat)       # check residuals
    # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
    #check_overdispersion(m.deg.post.clutch.2.pat) 
    
    ## i) male pre-manip strength association with 
    # all cluches total paternity
    m.strength.pre.tot.pat <- lm(total.paternity ~ 
                                    scale(strength.pre)
                                  # sensitivity, is strength associated with
                                  # with fecundity independent of age
                                  #+ Age.category
                                  ,
                                  #family = 'gaussian',
                                  data = subset(chr15_attrib_df,
                                                Sex == 'm' &
                                                  !is.na(strength.pre) &
                                                  !is.na(total.paternity)))
    
    ## j) model summary
    summary(m.strength.pre.tot.pat)    # model summary 
    confint(m.strength.pre.tot.pat, 
            level = 0.95, method = 'profile')
    #plot(m.strength.pre.tot.pat)       # check residuals
    # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
    #check_overdispersion(m.strength.pre.tot.pat) 
    
    ## k) male pre-manip degree association with 
    # all clutches total paternity
    m.deg.pre.tot.pat <- lm(total.paternity ~ 
                               scale(degree.pre)
                             # sensitivity, is strength associated with
                             # with fecundity independent of age
                             #+ Age.category
                             ,
                             #family = 'gaussian',
                             data = subset(chr15_attrib_df,
                                           Sex == 'm' &
                                             !is.na(degree.pre) &
                                             !is.na(total.paternity)))
    
    ## l) model summary
    summary(m.deg.pre.tot.pat)    # model summary 
    confint(m.deg.pre.tot.pat, 
            level = 0.95, method = 'profile')
    #plot(m.deg.pre.tot.pat)       # check residuals
    # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
    #check_overdispersion(m.deg.pre.tot.pat) 
    
    
    
#*****************************************************************************#    
    
    
    
  # Rationale - Among age and plumage measures (traits) that are NOT 
  #             associated  with reproductive success, determine
  #             their association with social network measures for Sup Mat
  

  ### 3.1 Female plumage trait associations with node level 
      # social network measures
    ## a) female breast brightness association with 
        # association strength pre-manip
      fem.Bbright.strength.pre <- lm(strength.pre ~ 
                                          scale(R_avg.bright),
                                          #family = 'gaussian',
                                          data = subset(chr15_attrib_df,
                                                Sex == 'f' &
                                                !is.na(strength.pre) &
                                                !is.na(R_avg.bright) &
                                                !is.na(total.fecundity)))
    
    ## b) model summary
      summary(fem.Bbright.strength.pre)    # model summary 
      confint(fem.Bbright.strength.pre, level = 0.95, method = 'profile')
      #plot(fem.Bbright.strength.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.Bbright.strength.pre)  
    
    ## c) female belly brightness association with 
        # association strength pre-manip
      fem.Bbright.strength.pre <- lm(strength.pre ~ 
                                          scale(B_avg.bright),
                                          #family = 'gaussian',
                                          data = subset(chr15_attrib_df,
                                                Sex == 'f' &
                                                !is.na(strength.pre) &
                                                !is.na(B_avg.bright) &
                                                !is.na(total.fecundity)))
    
    ## d) model summary
      summary(fem.Bbright.strength.pre)    # model summary 
      confint(fem.Bbright.strength.pre, level = 0.95, method = 'profile')
      #plot(fem.Bbright.strength.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.Bbright.strength.pre)  
    
    ## e) female breast brightness association 
        # with association degree pre-manip
      fem.Rbright.degree.pre <- lm(degree.pre ~ 
                                      scale(R_avg.bright),
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_df,
                                             Sex == 'f' &
                                             !is.na(degree.pre) &
                                             !is.na(R_avg.bright) &
                                             !is.na(total.fecundity)))
    
    ## f) model summary
      summary(fem.Rbright.degree.pre)    # model summary 
      confint(fem.Rbright.degree.pre, level = 0.95, method = 'profile')
      #plot(fem.Rbright.degree.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.Rbright.degree.pre)  
    
    ## g) female belly brightness association 
        # with association degree pre-manip
      fem.Bbright.degree.pre <- lm(degree.pre ~ 
                                      scale(B_avg.bright),
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_df,
                                            Sex == 'f' &
                                            !is.na(degree.pre) &
                                            !is.na(B_avg.bright) &
                                            !is.na(total.fecundity)))
    
    ## h) model summary
      summary(fem.Bbright.degree.pre)    # model summary 
      confint(fem.Bbright.degree.pre, level = 0.95, method = 'profile')
      #plot(fem.Bbright.degree.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.Bbright.degree.pre)  
      
    ## i) female breast brightness association with 
      # association strength post manip
      fem.Rbright.strength.post <- lm(strength.post ~ 
                                          scale(R_avg.bright),
                                          #family = 'gaussian',
                                          data = subset(chr15_attrib_post_df,
                                                Sex == 'f' &
                                                !is.na(strength.post) &
                                                !is.na(R_avg.bright) &
                                                !is.na(total.fecundity)))
      
    ## j) model summary
      summary(fem.Rbright.strength.post)    # model summary 
      confint(fem.Rbright.strength.post, level = 0.95, method = 'profile')
      #plot(fem.Rbright.strength.post)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.Rbright.strength.post)  
      
    ## k) female belly brightness association with 
      # association strength post manip
      fem.Bbright.strength.post <- lm(strength.post ~ 
                                          scale(B_avg.bright),
                                          #family = 'gaussian',
                                          data = subset(chr15_attrib_post_df,
                                                Sex == 'f' &
                                                !is.na(strength.post) &
                                                !is.na(B_avg.bright) &
                                                !is.na(total.fecundity)))
      
    ## l) model summary
      summary(fem.Bbright.strength.post)    # model summary 
      confint(fem.Bbright.strength.post, level = 0.95, method = 'profile')
      #plot(fem.Bbright.strength.post)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.Bbright.strength.post)  
      
    ## m) female breast brightness association 
      # with association degree post manip
      fem.Rbright.degree.post <- lm(degree.post ~ 
                                        scale(R_avg.bright),
                                        #family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                              Sex == 'f' &
                                              !is.na(degree.post) &
                                              !is.na(R_avg.bright) &
                                              !is.na(total.fecundity)))
      
    ## n) model summary
      summary(fem.Rbright.degree.post)    # model summary 
      confint(fem.Rbright.degree.post, level = 0.95, method = 'profile')
      #plot(fem.Rbright.degree.post)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.Rbright.degree.post)  
      
    ## o) female belly brightness association 
      # with association degree post manip
      fem.Bbright.degree.post <- lm(degree.post ~ 
                                        scale(B_avg.bright),
                                        #family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                              Sex == 'f' &
                                              !is.na(degree.post) &
                                              !is.na(B_avg.bright) &
                                              !is.na(total.fecundity)))
      
    ## p) model summary
      summary(fem.Bbright.degree.post)    # model summary 
      confint(fem.Bbright.degree.post, level = 0.95, method = 'profile')
      #plot(fem.Bbright.degree.post)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.Bbright.degree.post)  

      
  ### 3.2 Male plumage trait associations with node level 
    # social network measures
    ## a) male breast brightness association with 
        # association strength pre-manip
      m.Rbright.strength.pre <- lm(strength.pre ~ 
                                      scale(R_avg.bright),
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_df,
                                             Sex == 'm' &
                                             #Color.manipulation. == 'n' &
                                             !is.na(strength.pre) &
                                             !is.na(R_avg.bright) &
                                             !is.na(nest.2.3.tot.pat)))
      
    ## b) model summary
      summary(m.Rbright.strength.pre)    # model summary 
      confint(m.Rbright.strength.pre, level = 0.95, method = 'profile')
      #plot(m.Rbright.strength.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Rbright.strength.pre)  
    
    ## c) male belly brightness association with 
        # association strength pre-manip
      m.Bbright.strength.pre <- lm(strength.pre ~ 
                                      scale(B_avg.bright),
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_df,
                                             Sex == 'm' &
                                             #Color.manipulation. == 'n' &
                                             !is.na(strength.pre) &
                                             !is.na(B_avg.bright) &
                                             !is.na(nest.2.3.tot.pat)))
    
    ## d) model summary
      summary(m.Bbright.strength.pre)    # model summary 
      confint(m.Bbright.strength.pre, level = 0.95, method = 'profile')
     #plot(m.Bbright.strength.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Bbright.strength.pre)  
    
    ## e) male breast brightness association with 
        # association degree pre-manip
      m.Rbright.degree.pre <- lm(degree.pre ~ 
                                    scale(R_avg.bright),
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                           Sex == 'm' &
                                           #Color.manipulation. == 'n' &
                                           !is.na(degree.pre)&
                                           !is.na(R_avg.bright) &
                                           !is.na(nest.2.3.tot.pat)))
      
    ## f) model summary
      summary(m.Rbright.degree.pre)    # model summary 
      confint(m.Rbright.degree.pre, level = 0.95, method = 'profile')
      #plot(m.Rbright.degree.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Rbright.degree.pre)  
      
    ## g) male belly brightness association with 
        # association degree pre-manip
      m.Bbright.degree.pre <- lm(degree.pre ~ 
                                    scale(B_avg.bright),
                                    #family = 'gaussian',
                                    data = subset(chr15_attrib_df,
                                           Sex == 'm' &
                                           #Color.manipulation. == 'n' &
                                           !is.na(degree.pre) &
                                           !is.na(B_avg.bright) &
                                           !is.na(nest.2.3.tot.pat)))
    
    ## h) model summary
      summary(m.Bbright.degree.pre)    # model summary 
      confint(m.Bbright.degree.pre, level = 0.95, method = 'profile')
      #plot(m.Bbright.degree.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Bbright.degree.pre)  
      
    ## i) male experimental breast brightness association with 
      # association strength post manip
      m.post.Rbright.strength.post <- lm(strength.post ~ 
                                          scale(R.bright.treat.and.orig),
                                          #family = 'gaussian',
                                          data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              !is.na(strength.post) &
                                              !is.na(R.bright.treat.and.orig) &
                                              !is.na(attmpt.2.tot.pat)))
      
    ## j) model summary
      summary(m.post.Rbright.strength.post)    # model summary 
      confint(m.post.Rbright.strength.post, level = 0.95, method = 'profile') 
       #plot(m.post.Rbright.strength.post)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.post.Rbright.strength.post) 
      
    ## k) male belly brightness association with 
      # association strength post manip
      m.Bbright.strength.post <- lm(strength.post ~ 
                                        scale(B_avg.bright),
                                        #family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              #Color.manipulation. == 'n' &
                                              !is.na(strength.post) &
                                              !is.na(B_avg.bright) &
                                              !is.na(attmpt.2.tot.pat)))
      
      ## l) model summary
      summary(m.Bbright.strength.post)    # model summary 
      confint(m.Bbright.strength.post, level = 0.95, method = 'profile')
      #plot(m.Bbright.strength.post)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Bbright.strength.post)  
      
    ## m) male experimental breast brightness association with 
      # association degree post manip
      m.post.Rbright.degree.post <- lm(degree.post ~ 
                                        scale(R.bright.treat.and.orig),
                                        #family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              !is.na(degree.post) &
                                              !is.na(R.bright.treat.and.orig) &
                                              !is.na(attmpt.2.tot.pat)))
      
    ## n) model summary
      summary(m.post.Rbright.degree.post)    # model summary 
      confint(m.post.Rbright.degree.post, level = 0.95, method = 'profile') 
      #plot(m.post.Rbright.degree.post)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.post.Rbright.degree.post) 
      
    ## o) male belly brightness association with 
      # association degree post manip
      m.Bbright.degree.post <- lm(degree.post ~ 
                                      scale(B_avg.bright),
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_post_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(degree.post) &
                                            !is.na(B_avg.bright) &
                                            !is.na(attmpt.2.tot.pat)))
      
    ## p) model summary
      summary(m.Bbright.degree.post)    # model summary 
      confint(m.Bbright.degree.post, level = 0.95, method = 'profile')
      #plot(m.Bbright.degree.post)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Bbright.degree.post) 
      
      
      
###############################################################################
##############      4. Age association with plumage traits       ##############
###############################################################################         
      
  ### 4.1 Female age associations with plumage traits
    ## a) female age association with breast brightness
      fem.age.Rbright <- lm(R_avg.bright ~ 
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                      Sex == 'f' &
                                      !is.na(R_avg.bright) &
                                      !is.na(total.fecundity)))
      
    ## b) model summary
      summary(fem.age.Rbright)    # model summary 
      confint(fem.age.Rbright, level = 0.95, method = 'profile')  
      #plot(fem.age.Rbright)       # check residuals  
      
    ## c) female age association with belly brightness
      fem.age.Bbright <- lm(B_avg.bright ~ 
                                Age.category,
                                #family = 'gaussian',
                                data = subset(chr15_attrib_df,
                                      Sex == 'f' &
                                      !is.na(B_avg.bright) &
                                      !is.na(total.fecundity)))
      
    ## d) model summary
      summary(fem.age.Bbright)    # model summary 
      confint(fem.age.Bbright, level = 0.95, method = 'profile')  
      #plot(fem.age.Bbright)       # check residuals  
      
      
  ### 4.2 Male age associations with plumage traits
    ## a) Male age association with breast brightens
      m.age.Rbright <- lm(R_avg.bright ~ 
                              Age.category,
                              #family = 'gaussian',
                              data = subset(chr15_attrib_df,
                                    Sex == 'm' &
                                    #Color.manipulation. == 'n' &
                                    !is.na(R_avg.bright)  &
                                    !is.na(nest.2.3.tot.pat)))
      
    ## b) model summary
      summary(m.age.Rbright)    # model summary 
      confint(m.age.Rbright, level = 0.95, method = 'profile')  
      #plot(m.age.Rbright)       # check residuals  
      
    ## c) Male age association with belly brightens
      m.age.Bbright <- lm(B_avg.bright ~ 
                              Age.category,
                              #family = 'gaussian',
                              data = subset(chr15_attrib_df,
                                    Sex == 'm' &
                                    #Color.manipulation. == 'n' &
                                    !is.na(B_avg.bright)  &
                                    !is.na(nest.2.3.tot.pat)))
      
    ## d) model summary
      summary(m.age.Bbright)    # model summary 
      confint(m.age.Bbright, level = 0.95, method = 'profile')  
      #plot(m.age.Bbright)       # check residuals  
      

      
