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

      # load here packages
        library ('here')
      
    ## b) Modeling Packages
      
      # load MASS (negative binomial model)
        library('MASS')

      # load performance
        library('performance')
      
      # load lmPerm
        library('lmPerm')

        
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
      load(here('data/5_6_chr15_mediation_data.RData'))
    
      
      
###############################################################################
##############   3.Trait associations with reproductive success   #############
###############################################################################

  ### 3.1 Female age association with reproductive success
    ## a) female age association with total fecundity
      fem.age.tot.fecund <- glm(total.fecundity ~ 
                                    Age.category,
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_pre_df,
                                            Sex == 'f' &
                                            !is.na(total.fecundity)))
    
    ## b) model summary
      summary(fem.age.tot.fecund)    # model summary
      confint(fem.age.tot.fecund, level = 0.95, method = 'profile')
      #plot(fem.age.tot.fecund)
      
      
  ### 3.2 Male age association with reproductive success
    ## a) male age association with first clutch paternity
      m.age.clutch.1.pat <- glm(attmpt.1.tot.pat ~ 
                                      Age.category,
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(attmpt.1.tot.pat)))
     
    ## b) model summary
      summary(m.age.clutch.1.pat)    # model summary 
      confint(m.age.clutch.1.pat, level = 0.95, method = 'profile')
      #plot(m.age.clutch.1.pat)

    ## c) male age association with second clutch paternity
      m.age.clutch.2.pat <- glm(attmpt.2.tot.pat ~ 
                                  Age.category,
                                  family = 'gaussian',
                                  data = subset(chr15_attrib_post_df,
                                         Sex == 'm' &
                                         #Color.manipulation. == 'n' &
                                         !is.na(attmpt.2.tot.pat)))
     
    ## d) model summary
      summary(m.age.clutch.2.pat)    # model summary 
      confint(m.age.clutch.2.pat, level = 0.95, method = 'profile') 
     #plot(m.age.clutch.2.pat)       # check residuals

    ## e) male age association with all clutches, total paternity
      m.age.tot.pat <- glm(total.paternity ~ 
                                    Age.category,
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_pre_df,
                                          Sex == 'm' &
                                          #Color.manipulation. == 'n' &
                                          !is.na(total.paternity)))
      
      
      summary(m.age.tot.pat)    # model summary 
      confint(m.age.tot.pat, level = 0.95, method = 'profile') 
      #plot(m.age.tot.pat)       # check residuals
      
      
  ### 3.3 Female plumage trait associations with reproductive success
    ## a) female breast brightness association with 
      # total fecundity from all clutches
      fem.Rbright.tot.fecund <- glm(total.fecundity ~ 
                                        scale(R_avg.bright),
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_pre_df,
                                                Sex == 'f' &
                                                !is.na(R_avg.bright) &
                                                !is.na(total.fecundity)))

    ## b) model summary
      summary(fem.Rbright.tot.fecund)    # model summary 
      confint(fem.Rbright.tot.fecund, level = 0.95, method = 'profile') 
      #plot(fem.Rbright.tot.fecund)       # check residuals

    ## c) female belly brightness association with 
      # total offspring count
      fem.Bbright.tot.fecund <- glm(total.fecundity ~ 
                                        scale(B_avg.bright),
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_pre_df,
                                                Sex == 'f' &
                                                !is.na(B_avg.bright) &
                                                !is.na(total.fecundity)))

      summary(fem.Bbright.tot.fecund)    # model summary 
      confint(fem.Bbright.tot.fecund, level = 0.95, method = 'profile') 
      #plot(fem.Bbright.tot.fecund)       # check residuals
  
      
  ### 3.4 Male plumage trait associations with reproductive success
    ## a) male breast brightness association with 
      # clutch 1 paternity
      m.Rbright.clutch.1.pat <- glm(attmpt.1.tot.pat ~ 
                                      scale(R_avg.bright),
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(R_avg.bright) &
                                            !is.na(attmpt.1.tot.pat)))
      
    ## b) model summary
      summary(m.Rbright.clutch.1.pat)    # model summary 
      confint(m.Rbright.clutch.1.pat, level = 0.95, method = 'profile') 
      #plot(m.Rbright.clutch.1.pat)       # check residuals
      
    ## c) male belly brightness association with 
      # clutch 1 paternity
      m.Bbright.clutch.1.pat <- glm(attmpt.1.tot.pat ~ 
                                        scale(B_avg.bright),
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_pre_df,
                                              Sex == 'm' &
                                              #Color.manipulation. == 'n' &
                                              !is.na(B_avg.bright) &
                                              !is.na(attmpt.1.tot.pat)))
      
    ## d) model summary
      summary(m.Bbright.clutch.1.pat)    # model summary 
      confint(m.Bbright.clutch.1.pat, level = 0.95, method = 'profile') 
      #plot(m.Bbright.clutch.1.pat)       # check residuals
    
    ## e) male experimental breast brightness association 
      # with clutch 2 paternity
      m.Rbright.clutch.2.pat <- glm(attmpt.2.tot.pat ~ 
                                        scale(R.bright.treat.and.orig),
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              #Color.manipulation. == 'n' &
                                              !is.na(R.bright.treat.and.orig) &
                                              !is.na(attmpt.2.tot.pat)))
    ## f) model summary
      summary(m.Rbright.clutch.2.pat)    # model summary 
      confint(m.Rbright.clutch.2.pat, level = 0.95, method = 'profile')
     #plot(m.Rbright.clutch.2.pat)       # check residuals
      
    ## g) male belly brightness association
      # with clutch 2 paternity
      m.Bbright.clutch.2.pat <- glm(attmpt.2.tot.pat ~ 
                                        scale(B_avg.bright),
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                             Sex == 'm' &
                                             #Color.manipulation. == 'n' &
                                             !is.na(B_avg.bright) &
                                             !is.na(attmpt.2.tot.pat)))
    
    ## h) model summary
      summary(m.Bbright.clutch.2.pat)    # model summary 
      confint(m.Bbright.clutch.2.pat, level = 0.95, method = 'profile')  
      #plot(m.Bbright.clutch.2.pat)       # check residuals
     
    ## i) male breast brightness association with 
      # total paternity from all clutches
      m.Rbright.tot.pat <- lm(total.paternity ~ 
                                        scale(R_avg.bright) 
                                        + Color.manipulation.
                                        ,
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_pre_df,
                                               Sex == 'm' &
                                               #Color.manipulation. == 'n' &
                                               !is.na(R_avg.bright) &
                                               !is.na(attmpt.2.tot.pat)))
      
    ## j) model summary
      summary(m.Rbright.tot.pat)    # model summary 
      confint(m.Rbright.tot.pat, level = 0.95, method = 'profile')
      #plot(m.Rbright.tot.pat)       # check residuals
      
    ## k) male belly brightness association with 
      # total paternity from all clutches
      m.Bbright.tot.pat <- glm(total.paternity ~ 
                                  scale(B_avg.bright) 
                                  # + Color.manipulation.
                                  ,
                                  family = 'gaussian',
                                  data = subset(chr15_attrib_pre_df,
                                        Sex == 'm' &
                                        #Color.manipulation. == 'n' &
                                        !is.na(R_avg.bright) &
                                        !is.na(attmpt.2.tot.pat)))
      
    ## l) model summary
      summary(m.Bbright.tot.pat)    # model summary 
      confint(m.Bbright.tot.pat, level = 0.95, method = 'profile')
      #plot(m.Bbright.tot.pat)       # check residuals

 
    
###############################################################################
##############   4. Trait associations with soc. net. metrics    ##############
###############################################################################    
  
  ### 4.1 Female age associations with node level social network measures
    ## a) female age association with association strength
      # pre-manip
      fem.age.strength.pre <- glm(strength.pre ~ 
                                  Age.category,
                                  family = 'gaussian',
                                  data = subset(chr15_attrib_pre_df,
                                        Sex == 'f' &
                                        !is.na(strength.pre) &
                                        !is.na(total.fecundity)))
      
    ## b) model summary
      summary(fem.age.strength.pre)    # model summary 
      confint(fem.age.strength.pre, level = 0.95, method = 'profile')  
      #plot(fem.age.strength.pre)       # check residuals
      
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
      #                             data = subset(chr15_attrib_pre_df,
      #                                           Sex == 'f' &
      #                                             !is.na(strength.pre) &
      #                                             !is.na(total.fecundity)))
      # summary(fem.age.strength.pre.perm)    # model summary 
   
    ## c) female age association with association degree
      # pre-manip
      fem.age.degree.pre <- glm(degree.pre ~ 
                                  Age.category,
                                  family = 'gaussian',
                                  data = subset(chr15_attrib_pre_df,
                                          Sex == 'f' &
                                          !is.na(degree.pre) &
                                          !is.na(total.fecundity)))
    
    ## d) model summary
      summary(fem.age.degree.pre)    # model summary 
      confint(fem.age.degree.pre, level = 0.95, method = 'profile')
      #plot(fem.age.degree.pre)       # check residuals
      
    ## e) female age association with association strength
      # post manip
      fem.age.strength.post <- glm(strength.post ~ 
                                        Age.category,
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                              Sex == 'f' &
                                              !is.na(strength.post) &
                                              !is.na(total.fecundity)))
      
    ## f) model summary
      summary(fem.age.strength.post)    # model summary 
      confint(fem.age.strength.post, level = 0.95, method = 'profile')  
      #plot(fem.age.strength.post)       # check residuals
      
    ## g) female age association with association degree
      # post manip
      fem.age.degree.post <- glm(degree.post ~ 
                                      Age.category,
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_post_df,
                                            Sex == 'f' &
                                            !is.na(degree.post) &
                                            !is.na(total.fecundity)))
      
    ## h) model summary
      summary(fem.age.degree.post)    # model summary 
      confint(fem.age.degree.post, level = 0.95, method = 'profile')
      #plot(fem.age.degree.post)       # check residuals
      

  ### 4.2 Male age associations with node level social network measures
    ## a) male age association with with 
      # association strength pre-manip
      m.age.strength.pre <- glm(strength.pre ~ 
                                      Age.category,
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
                                              Sex == 'm' &
                                              !is.na(strength.pre) &
                                              !is.na(attmpt.1.tot.pat)))
      
    ## b) model summary
      summary(m.age.strength.pre)    # model summary 
      confint(m.age.strength.pre, level = 0.95, method = 'profile') 
      #plot(m.age.strength.pre)       # check residuals
          
    ## c) male age association with association degree
      # pre-manip
      m.age.degree.pre <- glm(degree.pre ~ 
                                    Age.category,
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_pre_df,
                                           Sex == 'm' &
                                           #Color.manipulation. == 'n' &
                                           !is.na(degree.pre) &
                                           !is.na(attmpt.1.tot.pat)))
        
    ## d) model summary
      summary(m.age.degree.pre)    # model summary 
      confint(m.age.degree.pre, level = 0.95, method = 'profile')
      #plot(m.age.degree.pre)       # check residuals

    ## e) male age association with with 
      # association strength post manip
      m.age.strength.post<- glm(strength.post ~ 
                                      Age.category,
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              #Color.manipulation. == 'n' &
                                              !is.na(strength.post) &
                                              !is.na(attmpt.2.tot.pat)))
        
    ## f) model summary
      summary(m.age.strength.post)    # model summary 
      confint(m.age.strength.post, level = 0.95, method = 'profile') 
      #plot(m.age.strength.post)       # check residuals
        
    ## g) male age association with association degree
      # post manip
      m.age.degree.post <- glm(degree.post ~ 
                                    Age.category,
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_post_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(degree.post) &
                                            !is.na(attmpt.2.tot.pat)))
        
    ## h) model summary
      summary(m.age.degree.post)    # model summary 
      confint(m.age.degree.post, level = 0.95, method = 'profile')
      #plot(m.age.degree.post)       # check residuals
        
     
  ### 4.3 Female plumage trait associations with node level 
      # social network measures
    ## a) female breast brightness association with 
        # association strength pre-manip
      fem.Bbright.strength.pre <- glm(strength.pre ~ 
                                          scale(R_avg.bright),
                                          family = 'gaussian',
                                          data = subset(chr15_attrib_pre_df,
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
      fem.Bbright.strength.pre <- glm(strength.pre ~ 
                                          scale(B_avg.bright),
                                          family = 'gaussian',
                                          data = subset(chr15_attrib_pre_df,
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
      fem.Rbright.degree.pre <- glm(degree.pre ~ 
                                      scale(R_avg.bright),
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
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
      fem.Bbright.degree.pre <- glm(degree.pre ~ 
                                      scale(B_avg.bright),
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
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
      fem.Rbright.strength.post <- glm(strength.post ~ 
                                          scale(R_avg.bright),
                                          family = 'gaussian',
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
      fem.Bbright.strength.post <- glm(strength.post ~ 
                                          scale(B_avg.bright),
                                          family = 'gaussian',
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
      fem.Rbright.degree.post <- glm(degree.post ~ 
                                        scale(R_avg.bright),
                                        family = 'gaussian',
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
      fem.Bbright.degree.post <- glm(degree.post ~ 
                                        scale(B_avg.bright),
                                        family = 'gaussian',
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

      
  ### 4.4 Male plumage trait associations with node level 
    # social network measures
    ## a) male breast brightness association with 
        # association strength pre-manip
      m.Rbright.strength.pre <- glm(strength.pre ~ 
                                      scale(R_avg.bright),
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
                                             Sex == 'm' &
                                             #Color.manipulation. == 'n' &
                                             !is.na(strength.pre) &
                                             !is.na(R_avg.bright) &
                                             !is.na(attmpt.1.tot.pat)))
      
    ## b) model summary
      summary(m.Rbright.strength.pre)    # model summary 
      confint(m.Rbright.strength.pre, level = 0.95, method = 'profile')
      #plot(m.Rbright.strength.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Rbright.strength.pre)  
    
    ## c) male belly brightness association with 
        # association strength pre-manip
      m.Bbright.strength.pre <- glm(strength.pre ~ 
                                      scale(B_avg.bright),
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
                                             Sex == 'm' &
                                             #Color.manipulation. == 'n' &
                                             !is.na(strength.pre) &
                                             !is.na(B_avg.bright) &
                                             !is.na(attmpt.1.tot.pat)))
    
    ## d) model summary
      summary(m.Bbright.strength.pre)    # model summary 
      confint(m.Bbright.strength.pre, level = 0.95, method = 'profile')
     #plot(m.Bbright.strength.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Bbright.strength.pre)  
    
    ## e) male breast brightness association with 
        # association degree pre-manip
      m.Rbright.degree.pre <- glm(degree.pre ~ 
                                    scale(R_avg.bright),
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_pre_df,
                                           Sex == 'm' &
                                           #Color.manipulation. == 'n' &
                                           !is.na(degree.pre)&
                                           !is.na(R_avg.bright) &
                                           !is.na(attmpt.1.tot.pat)))
      
    ## f) model summary
      summary(m.Rbright.degree.pre)    # model summary 
      confint(m.Rbright.degree.pre, level = 0.95, method = 'profile')
      #plot(m.Rbright.degree.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Rbright.degree.pre)  
      
    ## g) male belly brightness association with 
        # association degree pre-manip
      m.Bbright.degree.pre <- glm(degree.pre ~ 
                                    scale(B_avg.bright),
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_pre_df,
                                           Sex == 'm' &
                                           #Color.manipulation. == 'n' &
                                           !is.na(degree.pre) &
                                           !is.na(B_avg.bright) &
                                           !is.na(attmpt.1.tot.pat)))
    
    ## h) model summary
      summary(m.Bbright.degree.pre)    # model summary 
      confint(m.Bbright.degree.pre, level = 0.95, method = 'profile')
      #plot(m.Bbright.degree.pre)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.Bbright.degree.pre)  
      
    ## i) male experimental breast brightness association with 
      # association strength post manip
      m.post.Rbright.strength.post <- glm(strength.post ~ 
                                          scale(R.bright.treat.and.orig),
                                          family = 'gaussian',
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
      m.Bbright.strength.post <- glm(strength.post ~ 
                                        scale(B_avg.bright),
                                        family = 'gaussian',
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
      m.post.Rbright.degree.post <- glm(degree.post ~ 
                                        scale(R.bright.treat.and.orig),
                                        family = 'gaussian',
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
      m.Bbright.degree.post <- glm(degree.post ~ 
                                      scale(B_avg.bright),
                                      family = 'gaussian',
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
##############      5. Age association with plumage traits       ##############
###############################################################################         
      
  ### 5.1 Female age associations with plumage traits
    ## a) female age association with breast brightness
      fem.age.Rbright <- glm(R_avg.bright ~ 
                                Age.category,
                                family = 'gaussian',
                                data = subset(chr15_attrib_pre_df,
                                      Sex == 'f' &
                                      !is.na(R_avg.bright) &
                                      !is.na(total.fecundity)))
      
    ## b) model summary
      summary(fem.age.Rbright)    # model summary 
      confint(fem.age.Rbright, level = 0.95, method = 'profile')  
      #plot(fem.age.Rbright)       # check residuals  
      
    ## c) female age association with belly brightness
      fem.age.Bbright <- glm(B_avg.bright ~ 
                                Age.category,
                                family = 'gaussian',
                                data = subset(chr15_attrib_pre_df,
                                      Sex == 'f' &
                                      !is.na(B_avg.bright) &
                                      !is.na(total.fecundity)))
      
    ## d) model summary
      summary(fem.age.Bbright)    # model summary 
      confint(fem.age.Bbright, level = 0.95, method = 'profile')  
      #plot(fem.age.Bbright)       # check residuals  
      
      
  ### 5.2 Male age associations with plumage traits
    ## a) Male age association with breast brightens
      m.age.Rbright <- glm(R_avg.bright ~ 
                              Age.category,
                              family = 'gaussian',
                              data = subset(chr15_attrib_pre_df,
                                    Sex == 'm' &
                                    #Color.manipulation. == 'n' &
                                    !is.na(R_avg.bright)  &
                                    !is.na(attmpt.1.tot.pat)))
      
    ## b) model summary
      summary(m.age.Rbright)    # model summary 
      confint(m.age.Rbright, level = 0.95, method = 'profile')  
      #plot(m.age.Rbright)       # check residuals  
      
    ## c) Male age association with belly brightens
      m.age.Bbright <- lm(B_avg.bright ~ 
                              Age.category,
                              family = 'gaussian',
                              data = subset(chr15_attrib_pre_df,
                                    Sex == 'm' &
                                    #Color.manipulation. == 'n' &
                                    !is.na(B_avg.bright)  &
                                    !is.na(attmpt.1.tot.pat)))
      
    ## d) model summary
      summary(m.age.Bbright)    # model summary 
      confint(m.age.Bbright, level = 0.95, method = 'profile')  
      #plot(m.age.Bbright)       # check residuals  
      

      
###############################################################################
############  6. Soc. net. metrics association with reproduction  #############
###############################################################################        
      
  ### 6.1 Female social network measures associated with reproductive success
      ## a) female pre-manip strength association with 
          # total fecundity
        fem.strength.pre.tot.fecund <- glm(total.fecundity ~ 
                                        scale(strength.pre),
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_pre_df,
                                              Sex == 'f' &
                                              !is.na(strength.pre) &
                                              !is.na(total.fecundity)))
      
      
    ## b) model summary
      summary(fem.strength.pre.tot.fecund)    # model summary 
      confint(fem.strength.pre.tot.fecund, 
              level = 0.95, method = 'profile')
      #plot(fem.strength.pre.tot.fecund)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.strength.pre.tot.fecund) 
 
    ## c) female pre-manip degree association with 
      # total fecundity
      fem.deg.pre.tot.fecund <- glm(total.fecundity ~ 
                                    scale(degree.pre),
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_pre_df,
                                          Sex == 'f' &
                                          !is.na(degree.pre) &
                                          !is.na(total.fecundity)))
    
    ## d) model summary
      summary(fem.deg.pre.tot.fecund)    # model summary 
      confint(fem.deg.pre.tot.fecund, 
              level = 0.95, method = 'profile')
      #plot(fem.deg.pre.tot.fecund)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(fem.deg.pre.tot.fecund) 
      
    ## e) female post manip strength association with 
      # total fecundity
      fem.strength.post.tot.fecund <- glm(total.fecundity ~ 
                                          scale(strength.post),
                                          family = 'gaussian',
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
      fem.deg.post.tot.fecund <- glm(total.fecundity ~ 
                                        scale(degree.post),
                                        family = 'gaussian',
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
      
      
  ### 6.2 Male social network measures associated with reproductive success
    ## a) male pre-manip strength association with 
      # clutch 1 total paternity
      m.strength.pre.clutch.1.pat <- glm(attmpt.1.tot.pat ~ 
                                          scale(strength.pre),
                                          family = 'gaussian',
                                          data = subset(chr15_attrib_pre_df,
                                                Sex == 'm' &
                                                !is.na(strength.pre) &
                                                !is.na(attmpt.1.tot.pat)))
      
      
    ## b) model summary
      summary(m.strength.pre.clutch.1.pat)    # model summary 
      confint(m.strength.pre.clutch.1.pat, 
              level = 0.95, method = 'profile')
      #plot(m.strength.pre.clutch.1.pat)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.strength.pre.clutch.1.pat) 
      
    ## c) male pre-manip degree association with 
      # clutch 1 total paternity
      m.deg.pre.clutch.1.pat <- glm(attmpt.1.tot.pat ~ 
                                    scale(degree.pre),
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_pre_df,
                                          Sex == 'm' &
                                          !is.na(degree.pre) &
                                          !is.na(attmpt.1.tot.pat)))
      
      
    ## d) model summary
      summary(m.deg.pre.clutch.1.pat)    # model summary 
      confint(m.deg.pre.clutch.1.pat, 
              level = 0.95, method = 'profile')
      #plot(m.deg.pre.clutch.1.pat)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.deg.pre.clutch.1.pat)
      
    ## d) male post manip strength association with 
      # clutch 2 total paternity
      m.strength.post.clutch.2.pat <- glm(attmpt.2.tot.pat ~ 
                                          scale(strength.post),
                                          family = 'gaussian',
                                          data = subset(chr15_attrib_post_df,
                                                Sex == 'm' &
                                                !is.na(strength.post) &
                                                !is.na(attmpt.1.tot.pat)))
      
      
    ## f) model summary
      summary(m.strength.post.clutch.2.pat)    # model summary 
      confint(m.strength.post.clutch.2.pat, 
              level = 0.95, method = 'profile')
      #plot(m.strength.post.clutch.2.pat)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.strength.post.clutch.2.pat) 

    ## g) male post manip degree association with 
      # clutch 2 total paternity
      m.deg.post.clutch.2.pat <- glm(attmpt.2.tot.pat ~ 
                                      scale(degree.post),
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_post_df,
                                             Sex == 'm' &
                                             !is.na(degree.post) &
                                             !is.na(attmpt.1.tot.pat)))

    ## h) model summary
      summary(m.deg.post.clutch.2.pat)    # model summary 
      confint(m.deg.post.clutch.2.pat, 
              level = 0.95, method = 'profile')
      #plot(m.deg.post.clutch.2.pat)       # check residuals
      # check for over/under dispersion; dispers ratio >> 1 over; << 1 under
      #check_overdispersion(m.deg.post.clutch.2.pat) 

    ## i) male pre-manip strength association with 
      # all cluches total paternity
      m.strength.pre.tot.pat <- glm(total.paternity ~ 
                                      scale(strength.pre),
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
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
      m.deg.pre.tot.pat <- glm(total.paternity ~ 
                                    scale(degree.pre),
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_pre_df,
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
      
      

###############################################################################
##############      7. Age by soc. net. metrics interaction      ##############
###############################################################################        
     
  ### 7.1 Interaction models: Female trait by node level social network measures
        # associations with reproductive success
      
    ## a) Interaction model: female age by pre-manip strength association 
      # with total fecundity
      fem.age.x.pre.strength.fecund.intx <- glm(total.fecundity ~ 
                                      Age.category * strength.pre,
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'f' &
                                            !is.na(strength.pre) &
                                            !is.na(total.fecundity)))
        
    ## b) model summary
      summary(fem.age.x.pre.strength.fecund.intx)    # model summary 
      #plot(fem.age.x.pre.strength.fecund.intx)       # check residuals

    ## c) Interaction model: female age by post-manip degree association 
      # with total fecundity
      fem.age.x.post.degree.fecund.intx <- glm(total.fecundity ~ 
                                    Age.category * strength.post,
                                    family = 'gaussian',
                                    data = subset(chr15_attrib_post_df,
                                           Sex == 'f' &
                                           !is.na(strength.post) &
                                           !is.na(total.fecundity)))
      
    ## d) model summary
      summary(fem.age.x.post.degree.fecund.intx)    # model summary 
      #plot(fem.age.x.post.degree.fecund.intx)       # check residuals

      
  ### 7.2 Interaction models: Male age by node level social network measures
      # associations with reproductive success
    ## a) Interaction model: male age by pre-manip strength association 
      # with clutch 1 total paternity
      m.age.x.pre.strength.cluth.1.intx <- glm(attmpt.1.tot.pat ~ 
                                  Age.category * strength.pre,
                                  family = 'gaussian',
                                  data = subset(chr15_attrib_pre_df,
                                        Sex == 'm' &
                                        !is.na(strength.pre) &
                                        !is.na(attmpt.1.tot.pat)))
      
    ## b) model summary
      summary(m.age.x.pre.strength.cluth.1.intx)    # model summary 
      #plot(m.age.x.pre.strength.cluth.1.intx)       # check residuals

    ## c) Interaction model: male age by pre manip strength association 
      # with all clutches total paternity  
      m.age.x.pre.strength.tot.intx <- glm(total.paternity ~ 
                                  Age.category * strength.pre,
                                  family = 'gaussian',
                                  data = subset(chr15_attrib_pre_df,
                                         Sex == 'm' &
                                         !is.na(strength.pre) &
                                         !is.na(total.paternity)))
      
    ## d) model summary
      summary(m.age.x.pre.strength.tot.intx)    # model summary 
      #plot(m.age.x.pre.strength.tot.intx)       # check residuals
      
    ## e) Interaction model: male age by post manip degree association 
      # with clutch 2 total paternity  
      m.age.x.post.degree.cluth.2.intx <- glm(attmpt.2.tot.pat ~ 
                                      Age.category * degree.post,
                                      family = 'gaussian',
                                      data = subset(chr15_attrib_post_df,
                                            Sex == 'm' &
                                            !is.na(degree.post) &
                                            !is.na(attmpt.2.tot.pat)))
      
    ## f) model summary
      summary(m.age.x.post.degree.cluth.2.intx)    # model summary 
      #plot(m.age.x.post.degree.cluth.2.intx)       # check residuals

      
      
 
###############################################################################
##############                  8. Export data                   ##############
###############################################################################  

      
  ### 8.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 mediation analysis
      save(file = here('data/7_chr15_effect_decomp_data.RData'), 
           list = c('chr15_attrib_pre_df', 'chr15_attrib_post_df'))
      
    ## b) Save and export model objects CHR 2015 mediation analysis
      save(file = here('data/7_chr15_effect_decomp_mods.RData'), 
           list = c('fem.age.strength.pre', 'fem.age.degree.post', 
                    'm.age.strength.pre', 'm.age.degree.post'))
      