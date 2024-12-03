################################################################################
#############        The role of age and plumage traits as         #############  
#############       determinants of reproductive success and       #############
#############       the mediating role of social interactions      #############
#############                                                      #############
#############                 6. Mediation analysis                #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 29 Sept 2024                #############
#############               last updated: 3 Dec 2024               #############
################################################################################



### PURPOSE: Run mediation models for node level assortative traits
  ## Model Sets
    # Model context 1: tot repro ~ original color & pre-manip soc intx.
                    #  egg/1st clutch ~ original color & pre-manip soc intx.
          # NOTE: for tot repro, male models should be subset to those NOT manip.

    # Model context 2: replace clutch ~ post-manip color & post-manip soc intx.
          # NOTE: for replace clutch, male models should be subset to those manip.

    # Model context 3: replace clutch ~ post-manip color & pre-manip soc intx.
          # NOTE: for replace clutch, male models should be subset to those manip.

  ## Social network measures, for both pre and post manip., 
    # Social intx metrics 1: all females and males
    
    # Social intx metrics 2: females and males among social mates

    # Social intx metrics 3: females and males excluding social mates

# Should social interactions be restricted to include only manipulated vs 
# unmanipulated males when reproductive success is also limited to males that
# are unmanipulated vs manipulated
          
    
                        
  
  
  # Code Blocks
    # 1. Configure work space
    # 2. Load RData
    # 3. Trait associations with reproductive success
    # 4. Trait associations with soc. net. metrics 
    # 5. Age association with plumage traits 
    # 6. Soc. net. metrics association with reproduction
    # 7. Age by soc. net. metrics interaction
    # 8. Effects decomposition

    

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
     
    ## b) Graph Plotting and Visualization Packages
      # load ggplot2 packages
      library ('ggplot2')
      
      # load ggeffects packages
        library('ggeffects')
      
      # load gridExtra packages
        library('gridExtra')
      
    ## c) Modeling Packages
      
      # load MASS (negative binomial model)
        library('MASS')
      
      # load mediation
        library('mediation')
      
      # load mediation for count based mediation analysis
        library('maczic')
    
      # load performance
        library('performance')
      
      # library('broom')
        library('broom.mixed')
      
      # load dharma
        library('DHARMa')

        
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
  
    # NOTEs on reproductive success    
    # Female total.offspring = Clutch.1.eggs + Clutch.2.Kids + Clutch.3.Kids
    # Male total.offspring = Clutch1.Offspring + Clutch.2.WP.Kids + Clutch.2.EP   
    # + Clutch.2.WP.Kids + Clutch.2.EP 
    
  ### 3.1 Female age association with reproductive success
    ## a) Unadjusted model: female age association with total fecundity
      fem.age.tot.fecund.nb <- glm.nb(total.fecundity ~ 
                                    Age.category,
                                    data = subset(chr15_attrib_pre_df,
                                            Sex == 'f' &
                                            !is.na(total.fecundity)))
    
      fem.age.tot.fecund.nb.theta <- glm(total.fecundity ~ 
                                      Age.category,
                                      family = negative.binomial(theta = 1),
                                      data = subset(chr15_attrib_pre_df,
                                              Sex == 'f' &
                                              !is.na(total.fecundity)))
      
      fem.age.tot.fecund.poiss <- glm(total.fecundity ~ 
                                       Age.category,
                                      family = 'poisson',
                                      data = subset(chr15_attrib_pre_df,
                                                Sex == 'f' &
                                                !is.na(total.fecundity)))
      
      
      fem.age.tot.fecund.lm <- lm(total.fecundity ~ 
                                      Age.category,
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'f' &
                                            !is.na(total.fecundity)))
    
    ## b) model summary
      summary(fem.age.tot.fecund.nb)    # model summary
      summary(fem.age.tot.fecund.nb.theta)
      summary(fem.age.tot.fecund.poiss)
      summary(fem.age.tot.fecund.lm)
      plot(fem.age.tot.fecund.lm)
      
  NEED TO UPDATE
      confint(fem.age.tot.fecund.unadj, level = 0.95, method = 'profile')
      exp(fem.age.tot.fecund.unadj$coefficients)
      exp(confint(fem.age.tot.fecund.unadj, level = 0.95, method = 'profile'))
      
      #plot(fem.age.tot.fecund.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.age.tot.fecund.unadj) 
      check_overdispersion(fem.age.tot.fecund.unadj.poiss) 
      # check boot strap CI
      # set.seed(2019) 
      # confint(fem.age.tot_offsprng.unadj, level = 0.95 , method = 'boot', 
      #         boot.type = 'perc')
      
    # ## c) Unadjusted model: female age association with clutch 1 egg count
    #   fem.age.nest1.eggs.unadj <- glm(Clutch.1.Eggs ~ 
    #                                       Age.category,
    #                                   family = 'poisson',
    #                                   data = subset(chr15_attrib_df,
    #                                         Sex == 'f' &
    #                                         !is.na(Clutch.1.Eggs)))
    # 
    # ## d) model summary
    #   summary(fem.age.nest1.eggs.unadj)    # model summary 
    #   confint(fem.age.nest1.eggs.unadj, level = 0.95, method = 'profile')
    #   exp(fem.age.nest1.eggs.unadj$coefficients)
    #   exp(confint(fem.age.nest1.eggs.unadj, level = 0.95, method = 'profile'))
    #   
    #   #plot(fem.age.nest1.eggs.unadj)       # check residuals
    #   # check for over/under dispersion; dispers ration >> 1 over; << 1 under
    #   check_overdispersion(fem.age.nest1.eggs.unadj)
    #   
    # ## c) Unadjusted model: female age association with clutch 2 egg count
    #   fem.age.nest2.kids.unadj <- glm(Clutch.2.Kids ~ 
    #                                     Age.category,
    #                                   family = 'poisson',
    #                                   data = subset(chr15_attrib_df,
    #                                                 Sex == 'f' &
    #                                                 !is.na(Clutch.2.Kids)))
    #   
    # ## d) model summary
    #   summary(fem.age.nest2.kids.unadj)    # model summary 
    #   confint(fem.age.nest2.kids.unadj, level = 0.95, method = 'profile')
    #   exp(fem.age.nest2.kids.unadj$coefficients)
    #   exp(confint(fem.age.nest2.kids.unadj, level = 0.95, method = 'profile'))
    #   
    #   #plot(fem.age.nest2.kids.unadj)       # check residuals
    #   # check for over/under dispersion; dispers ration >> 1 over; << 1 under
    #   check_overdispersion(fem.age.nest2.kids.unadj)
    # 
    # 
  ### 3.2 Male age association with reproductive success
    ## a) Unadjusted model: male age association with first clutch paternity
      m.age.clutch.1.pat.nb <- glm.nb(attmpt.1.tot.pat ~ 
                                      Age.category,
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(attmpt.1.tot.pat)))
      
      m.age.clutch.1.pat.nb.theta <- glm(attmpt.1.tot.pat ~ 
                                        Age.category,
                                        family = negative.binomial(theta = 1),
                                        data = subset(chr15_attrib_pre_df,
                                                  Sex == 'm' &
                                                  #Color.manipulation. == 'n' &
                                                  !is.na(attmpt.1.tot.pat)))
      
      m.age.clutch.1.pat.poiss <- glm(attmpt.1.tot.pat ~ 
                                      Age.category,
                                      family = 'poisson',
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(attmpt.1.tot.pat)))
      
      
      m.age.clutch.1.pat.lm <- lm(attmpt.1.tot.pat ~ 
                                  Age.category,
                                  data = subset(chr15_attrib_pre_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(attmpt.1.tot.pat)))
      
    ## b) model summary
      summary(m.age.clutch.1.pat.nb)    # model summary 
      summary(m.age.clutch.1.pat.nb.theta)
      summary(m.age.clutch.1.pat.poiss)
      summary(m.age.clutch.1.pat.lm)
      plot(m.age.clutch.1.pat.lm)
      
  NEED TO UPDATE
      confint(m.age.clutch.1.pat.unadj, level = 0.95, method = 'profile')
      exp(m.age.clutch.1.pat.unadj$coefficients)
      exp(confint(m.age.clutch.1.pat.unadj, level = 0.95, method = 'profile'))
      
      #plot(m.age.clutch.1.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.age.clutch.1.pat.unadj)
      check_overdispersion(m.age.clutch.1.pat.unadj.poiss)
      
      
      
#**** NOTE should we use the post manip data here    
    ## c) Unadjusted model: male age association with second clutch paternity
      m.age.clutch.2.pat.unadj <- glm.nb(attmpt.2.tot.pat ~ 
                                  Age.category,
                                  data = subset(chr15_attrib_post_df,
                                         Sex == 'm' &
                                         #Color.manipulation. == 'n' &
                                         !is.na(attmpt.2.tot.pat)))
      
      m.age.clutch.2.pat.unadj.poiss <- glm(attmpt.2.tot.pat ~ 
                                    Age.category,
                                    family = 'poisson',
                                    data = subset(chr15_attrib_post_df,
                                           Sex == 'm' &
                                           #Color.manipulation. == 'n' &
                                           !is.na(attmpt.2.tot.pat)))
      
    ## d) model summary
      summary(m.age.clutch.2.pat.unadj)    # model summary 
      confint(m.age.clutch.2.pat.unadj, level = 0.95, method = 'profile') 
      exp(m.age.clutch.2.pat.unadj$coefficients)
      exp(confint(m.age.clutch.2.pat.unadj, level = 0.95, method = 'profile'))
      
      #plot(m.age.clutch.2.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.age.clutch.2.pat.unadj)
      check_overdispersion(m.age.clutch.2.pat.unadj.poiss)
      
    ## e) Unadjusted model: male age association with all clutches, total paternity
      m.age.tot.pat.unadj <- glm.nb(total.paternity ~ 
                                    Age.category,
                                    data = subset(chr15_attrib_pre_df,
                                          Sex == 'm' &
                                          #Color.manipulation. == 'n' &
                                          !is.na(total.paternity)))
      
      m.age.tot.pat.unadj.poiss <- glm(total.paternity ~ 
                                  Age.category,
                                  family = 'poisson',
                                  data = subset(chr15_attrib_pre_df,
                                          Sex == 'm' &
                                          #Color.manipulation. == 'n' &
                                          !is.na(total.paternity)))
      
    ## f) model summary
      summary(m.age.tot.pat.unadj)    # model summary 
      confint(m.age.tot.pat.unadj, level = 0.95, method = 'profile') 
      exp(m.age.tot.pat.unadj$coefficients)
      exp(confint(m.age.tot.pat.unadj, level = 0.95, method = 'profile'))
      
      #plot(m.age.tot.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.age.tot.pat.unadj)
      check_overdispersion(m.age.tot.pat.unadj.poiss)
      
    #   
    # ## g) Unadjusted model: male age association with clutch 2 paternity
    #   m.age.nest2.pat.unadj <- glm(replace.pat ~ 
    #                                   Age.category,
    #                                 family = 'n',
    #                                 data = subset(chr15_attrib_df,
    #                                             Sex == 'm' &
    #                                             #Color.manipulation. == 'n' &
    #                                             !is.na(replace.pat)))
    #   
    # ## h) model summary
    #   summary(m.age.nest2.pat.unadj)    # model summary 
    #   confint(m.age.nest2.pat.unadj, level = 0.95, method = 'profile') 
    #   exp(m.age.nest2.pat.unadj$coefficients)
    #   exp(confint(m.age.nest2.pat.unadj, level = 0.95, method = 'profile'))
    #   
    #   #plot(m.age.nest2.pat.unadj)       # check residuals
    #   # check for over/under dispersion; dispers ration >> 1 over; << 1 under
    #   check_overdispersion(m.age.nest2.pat.unadj)
    # 
    #   
  ### 3.3 Female plumage trait associations with reproductive success
    ## a) Unadjusted model: female breast brightness association with 
      # total fecundity from all clutches
      fem.Rbright.tot.fecund.unadj <- glm.nb(total.fecundity ~ 
                                        scale(R_avg.bright),
                                        data = subset(chr15_attrib_pre_df,
                                                Sex == 'f' &
                                                !is.na(R_avg.bright) &
                                                !is.na(total.fecundity)))

      # fem.Rbright.tot_offsprng.unadj <- glm(total.offspring ~ 
      #                                     scale(R_avg.bright),
      #                                   family = 'poisson',
      #                                   data = subset(chr15_attrib_df,
      #                                         Sex == 'f' &
      #                                         !is.na(total.offspring) &
      #                                         !is.na(R_avg.bright)))
      
    ## b) model summary
      summary(fem.Rbright.tot.fecund.unadj)    # model summary 
      confint(fem.Rbright.tot.fecund.unadj, level = 0.95, method = 'profile') 
      exp(fem.Rbright.tot.fecund.unadj$coefficients)
      exp(confint(fem.Rbright.tot.fecund.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Rbright.tot.fecund.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Rbright.tot.fecund.unadj)

    ## c) Unadjusted model: female belly brightness association with 
      # total offspring count
      fem.Bbright.tot.fecund.unadj <- glm.nb(total.fecundity ~ 
                                        scale(B_avg.bright),
                                        data = subset(chr15_attrib_pre_df,
                                                Sex == 'f' &
                                                !is.na(B_avg.bright) &
                                                !is.na(total.fecundity)))
      
      # fem.Bbright.tot_offsprng.unadj <- glm(total.offspring ~ 
      #                                         scale(B_avg.bright),
      #                                       family = 'poisson',
      #                                       data = subset(chr15_attrib_df,
      #                                             Sex == 'f' &
      #                                             !is.na(total.offspring) &
      #                                             !is.na(B_avg.bright)))
      
    ## d) model summary
      summary(fem.Bbright.tot.fecund.unadj)    # model summary 
      confint(fem.Bbright.tot.fecund.unadj, level = 0.95, method = 'profile') 
      exp(fem.Bbright.tot.fecund.unadj$coefficients)
      exp(confint(fem.Bbright.tot.fecund.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Bbright.tot.fecund.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Bbright.tot.fecund.unadj)
      
    # ## e) Unadjusted model: female breast brightness association 
    #     # with clutch 1 egg count
    #   fem.Rbright.nest1.eggs.unadj <- glm(Clutch.1.Eggs ~ 
    #                                     scale(R_avg.bright),
    #                                   family = 'poisson',
    #                                   data = subset(chr15_attrib_df,
    #                                         Sex == 'f' &
    #                                         !is.na(Clutch.1.Eggs) &
    #                                         !is.na(R_avg.bright)))
    #   
    # ## f) model summary
    #   summary(fem.Rbright.nest1.eggs.unadj)    # model summary 
    #   confint(fem.Rbright.nest1.eggs.unadj, level = 0.95, method = 'profile')
    #   exp(fem.Rbright.nest1.eggs.unadj$coefficients)
    #   exp(confint(fem.Rbright.nest1.eggs.unadj, level = 0.95, 
    #               method = 'profile'))
    #   
    #   #plot(fem.Rbright.nest1.eggs.unadj)       # check residuals
    #   # check for over/under dispersion; dispers ration >> 1 over; << 1 under
    #   check_overdispersion(fem.Rbright.nest1.eggs.unadj)
    #   
    # ## g) Unadjusted model: female belly brightness association 
    #   # with clutch 1 egg count
    #   fem.Bbright.nest1.eggs.unadj <- glm(Clutch.1.Eggs ~ 
    #                                         scale(B_avg.bright),
    #                                       family = 'poisson',
    #                                       data = subset(chr15_attrib_df,
    #                                             Sex == 'f' &
    #                                             !is.na(Clutch.1.Eggs) &
    #                                             !is.na(B_avg.bright)))
    #   
    # ## h) model summary
    #   summary(fem.Bbright.nest1.eggs.unadj)    # model summary 
    #   confint(fem.Bbright.nest1.eggs.unadj, level = 0.95, method = 'profile')
    #   exp(fem.Bbright.nest1.eggs.unadj$coefficients)
    #   exp(confint(fem.Bbright.nest1.eggs.unadj, level = 0.95, 
    #               method = 'profile'))
    #   
    #   #plot(fem.Bbright.nest1.eggs.unadj)       # check residuals
    #   # check for over/under dispersion; dispers ration >> 1 over; << 1 under
    #   check_overdispersion(fem.Bbright.nest1.eggs.unadj)
    #   
    # ## i) Unadjusted model: female breast brightness association 
    #   # with clutch 2 offspring count
    #   fem.Rbright.nest2.kids.unadj <- glm(Clutch.2.Kids ~ 
    #                                         scale(R_avg.bright),
    #                                       family = 'poisson',
    #                                       data = subset(chr15_attrib_df,
    #                                                 Sex == 'f' &
    #                                                 !is.na(Clutch.2.Kids) &
    #                                                 !is.na(R_avg.bright)))
    #   
    # ## j) model summary
    #   summary(fem.Rbright.nest2.kids.unadj)    # model summary 
    #   confint(fem.Rbright.nest2.kids.unadj, level = 0.95, method = 'profile')
    #   exp(fem.Rbright.nest2.kids.unadj$coefficients)
    #   exp(confint(fem.Rbright.nest2.kids.unadj, level = 0.95, 
    #               method = 'profile'))
    #   
    #   #plot(fem.Rbright.nest2.kids.unadj)       # check residuals
    #   # check for over/under dispersion; dispers ration >> 1 over; << 1 under
    #   check_overdispersion(fem.Rbright.nest2.kids.unadj)
    #   
    # ## k) Unadjusted model: female belly brightness association 
    #   # with clutch 2 offspring count
    #   fem.Bbright.nest2.kids.unadj <- glm(Clutch.2.Kids ~ 
    #                                         scale(B_avg.bright),
    #                                       family = 'poisson',
    #                                       data = subset(chr15_attrib_df,
    #                                                 Sex == 'f' &
    #                                                 !is.na(Clutch.2.Kids) &
    #                                                 !is.na(B_avg.bright)))
    #   
    # ## l) model summary
    #   summary(fem.Bbright.nest2.kids.unadj)    # model summary 
    #   confint(fem.Bbright.nest2.kids.unadj, level = 0.95, method = 'profile')
    #   exp(fem.Bbright.nest2.kids.unadj$coefficients)
    #   exp(confint(fem.Bbright.nest2.kids.unadj, level = 0.95, 
    #               method = 'profile'))
    #   
    #   #plot(fem.Bbright.nest2.kids.unadj)       # check residuals
    #   # check for over/under dispersion; dispers ration >> 1 over; << 1 under
    #   check_overdispersion(fem.Bbright.nest2.kids.unadj)
      
      
  ### 3.4 Male plumage trait associations with reproductive success
    ## a) Unadjusted model: male breast brightness association with 
      # clutch 1 paternity
      m.Rbright.clutch.1.pat.unadj <- glm.nb(attmpt.1.tot.pat ~ 
                                      scale(R_avg.bright),,
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'm' &
                                            #Color.manipulation. == 'n' &
                                            !is.na(R_avg.bright) &
                                            !is.na(attmpt.1.tot.pat)))
      
      # m.Rbright.tot_offsprng.unadj <- glm(total.offspring ~ 
      #                                 scale(R_avg.bright),
      #                                 family = 'poisson',
      #                                 data = subset(chr15_attrib_df,
      #                                       Sex == 'm' &
      #                                       #Color.manipulation. == 'n' &
      #                                       !is.na(total.offspring) &
      #                                       !is.na(R_avg.bright)))
      
    ## b) model summary
      summary(m.Rbright.clutch.1.pat.unadj)    # model summary 
      confint(m.Rbright.clutch.1.pat.unadj, level = 0.95, method = 'profile') 
      exp(m.Rbright.clutch.1.pat.unadj$coefficients)
      exp(confint(m.Rbright.clutch.1.pat.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Rbright.clutch.1.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Rbright.clutch.1.pat.unadj)
      
    ## c) Unadjusted model: male belly brightness association with 
      # clutch 1 paternity
      m.Bbright.clutch.1.pat.unadj <- glm.nb(attmpt.1.tot.pat ~ 
                                        scale(B_avg.bright),,
                                        data = subset(chr15_attrib_pre_df,
                                              Sex == 'm' &
                                              #Color.manipulation. == 'n' &
                                              !is.na(B_avg.bright) &
                                              !is.na(attmpt.1.tot.pat)))
      
      # m.Bbright.tot_offsprng.unadj <- glm(total.offspring ~ 
      #                                       scale(B_avg.bright),
      #                                     family = 'poisson',
      #                                     data = subset(chr15_attrib_df,
      #                                           Sex == 'm' &
      #                                           #Color.manipulation. == 'n' &
      #                                           !is.na(total.offspring) &
      #                                           !is.na(B_avg.bright)))
      
    ## d) model summary
      summary(m.Bbright.clutch.1.pat.unadj)    # model summary 
      confint(m.Bbright.clutch.1.pat.unadj, level = 0.95, method = 'profile') 
      exp(m.Bbright.clutch.1.pat.unadj$coefficients)
      exp(confint(m.Bbright.clutch.1.pat.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Bbright.clutch.1.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Bbright.clutch.1.pat.unadj)
    
    ## e) Unadjusted model: male experimental breast brightness association 
      # with clutch 2 paternity
      m.Rbright.clutch.2.pat.unadj <- glm.nb(attmpt.2.tot.pat ~ 
                                        scale(R.bright.treat.and.orig),
                                        data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              #Color.manipulation. == 'n' &
                                              !is.na(R.bright.treat.and.orig) &
                                              !is.na(attmpt.2.tot.pat)))
    ## f) model summary
      summary(m.Rbright.clutch.2.pat.unadj)    # model summary 
      confint(m.Rbright.clutch.2.pat.unadj, level = 0.95, method = 'profile')
      exp(m.Rbright.clutch.2.pat.unadj$coefficients)
      exp(confint(m.Rbright.clutch.2.pat.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Rbright.clutch.2.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Rbright.clutch.2.pat.unadj)
      
    ## g) Unadjusted model: male belly brightness association
      # with clutch 2 paternity
      m.Bbright.clutch.2.pat.unadj <- glm.nb(attmpt.2.tot.pat ~ 
                                        scale(B_avg.bright),
                                        data = subset(chr15_attrib_post_df,
                                             Sex == 'm' &
                                             #Color.manipulation. == 'n' &
                                             !is.na(B_avg.bright) &
                                             !is.na(attmpt.2.tot.pat)))
    
    ## h) model summary
      summary(m.Bbright.clutch.2.pat.unadj)    # model summary 
      confint(m.Bbright.clutch.2.pat.unadj, level = 0.95, method = 'profile')  
      exp(m.Bbright.clutch.2.pat.unadj$coefficients)
      exp(confint(m.Bbright.clutch.2.pat.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Bbright.clutch.2.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Bbright.clutch.2.pat.unadj)
      
    ## i) Unadjusted model: male breast brightness association with 
      # total paternity from all clutches
      m.Rbright.tot.pat.unadj <- glm.nb(total.paternity ~ 
                                        scale(R_avg.bright) 
                                        + Color.manipulation.
                                        ,
                                        data = subset(chr15_attrib_pre_df,
                                               Sex == 'm' &
                                               #Color.manipulation. == 'n' &
                                               !is.na(R_avg.bright) &
                                               !is.na(attmpt.2.tot.pat)))
      
    ## j) model summary
      summary(m.Rbright.tot.pat.unadj)    # model summary 
      confint(m.Rbright.tot.pat.unadj, level = 0.95, method = 'profile')
      exp(m.Rbright.tot.pat.unadj$coefficients)
      exp(confint(m.Rbright.tot.pat.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Rbright.tot.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Rbright.tot.pat.unadj)
      
    ## k) Unadjusted model: male belly brightness association with 
      # total paternity from all clutches
      m.Bbright.tot.pat.unadj <- glm.nb(total.paternity ~ 
                                  scale(B_avg.bright) 
                                  # + Color.manipulation.
                                  ,
                                  data = subset(chr15_attrib_pre_df,
                                        Sex == 'm' &
                                        #Color.manipulation. == 'n' &
                                        !is.na(R_avg.bright) &
                                        !is.na(attmpt.2.tot.pat)))
      
    ## l) model summary
      summary(m.Bbright.tot.pat.unadj)    # model summary 
      confint(m.Bbright.tot.pat.unadj, level = 0.95, method = 'profile')
      exp(m.Bbright.tot.pat.unadj$coefficients)
      exp(confint(m.Bbright.tot.pat.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Bbright.tot.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Bbright.tot.pat.unadj)
      
 
    
###############################################################################
##############   4. Trait associations with soc. net. metrics    ##############
###############################################################################    
  
      
  ### 4.1 Female age associations with node level social network measures
    ## a) Unadjusted model: female age association with association strength
      # pre-manip
      fem.age.strength.pre.unadj <- glm.nb(strength.pre ~ 
                                      Age.category,
                                      data = subset(chr15_attrib_pre_df,
                                             Sex == 'f' &
                                             !is.na(strength.pre)))
      
    ## b) model summary
      summary(fem.age.strength.pre.unadj)    # model summary 
      confint(fem.age.strength.pre.unadj, level = 0.95, method = 'profile')  
      #plot(fem.age.strength.pre.unadj)       # check residuals
      exp(fem.age.strength.pre.unadj$coefficients)
      exp(confint(fem.age.strength.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.age.strength.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.age.strength.pre.unadj)
     
    ## c) Unadjusted model: female age association with association degree
      # pre-manip
      fem.age.degree.pre.unadj <- glm.nb(degree.pre ~ 
                                      Age.category,
                                      data = subset(chr15_attrib_pre_df,
                                             Sex == 'f' &
                                             !is.na(degree.pre)))
    
    ## d) model summary
      summary(fem.age.degree.pre.unadj)    # model summary 
      confint(fem.age.degree.pre.unadj, level = 0.95, method = 'profile')
      #plot(fem.age.degree.pre.unadj)       # check residuals
      exp(fem.age.degree.pre.unadj$coefficients)
      exp(confint(fem.age.degree.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.age.degree.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.age.degree.pre.unadj)
    
    ## e) Unadjusted model: female age association with association strength
      # post manip
      fem.age.strength.post.unadj <- glm.nb(strength.post ~ 
                                          Age.category,
                                        data = subset(chr15_attrib_post_df,
                                                  Sex == 'f' &
                                                  !is.na(strength.post)))
      
    ## f) model summary
      summary(fem.age.strength.post.unadj)    # model summary 
      confint(fem.age.strength.post.unadj, level = 0.95, method = 'profile')  
      #plot(fem.age.strength.post.unadj)       # check residuals
      exp(fem.age.strength.post.unadj$coefficients)
      exp(confint(fem.age.strength.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.age.strength.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.age.strength.post.unadj)
      
    ## g) Unadjusted model: female age association with association degree
      # post manip
      fem.age.degree.post.unadj <- glm.nb(degree.post ~ 
                                        Age.category,
                                      data = subset(chr15_attrib_post_df,
                                                Sex == 'f' &
                                                !is.na(degree.post)))
      
    ## h) model summary
      summary(fem.age.degree.post.unadj)    # model summary 
      confint(fem.age.degree.post.unadj, level = 0.95, method = 'profile')
      #plot(fem.age.degree.post.unadj)       # check residuals
      exp(fem.age.degree.post.unadj$coefficients)
      exp(confint(fem.age.degree.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.age.degree.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.age.degree.post.unadj)
      
      
    ### 4.2 Male age associations with node level social network measures
      ## a) Unadjusted model: male age association with with 
          # association strength pre-manip
        m.age.strength.pre.unadj <- glm.nb(strength.pre ~ 
                                          Age.category,
                                        data = subset(chr15_attrib_pre_df,
                                               Sex == 'm' &
                                               #Color.manipulation. == 'n' &
                                               !is.na(strength.pre)))
      
      ## b) model summary
        summary(m.age.strength.pre.unadj)    # model summary 
        confint(m.age.strength.pre.unadj, level = 0.95, method = 'profile') 
        #plot(m.age.strength.pre.unadj)       # check residuals
        exp(m.age.strength.pre.unadj$coefficients)
        exp(confint(m.age.strength.pre.unadj, level = 0.95, 
                    method = 'profile'))
        
        #plot(m.age.strength.pre.unadj)       # check residuals
        # check for over/under dispersion; dispers ration >> 1 over; << 1 under
        check_overdispersion(m.age.strength.pre.unadj)
      
      ## c) Unadjusted model: male age association with association degree
        # pre-manip
        m.age.degree.pre.unadj <- glm.nb(degree.pre ~ 
                                      Age.category,
                                      data = subset(chr15_attrib_pre_df,
                                             Sex == 'm' &
                                             #Color.manipulation. == 'n' &
                                             !is.na(degree.pre)))
        
      ## d) model summary
        summary(m.age.degree.pre.unadj)    # model summary 
        confint(m.age.degree.pre.unadj, level = 0.95, method = 'profile')
        #plot(m.age.degree.pre.unadj)       # check residuals
        exp(m.age.degree.pre.unadj$coefficients)
        exp(confint(m.age.degree.pre.unadj, level = 0.95, 
                    method = 'profile'))
        
        #plot(m.age.degree.pre.unadj)       # check residuals
        # check for over/under dispersion; dispers ration >> 1 over; << 1 under
        check_overdispersion(m.age.degree.pre.unadj)
      
      ## e) Unadjusted model: male age association with with 
        # association strength post manip
        m.age.strength.post.unadj <- glm.nb(strength.post ~ 
                                          Age.category,
                                        data = subset(chr15_attrib_post_df,
                                                Sex == 'm' &
                                                #Color.manipulation. == 'n' &
                                                !is.na(strength.post)))
        
      ## f) model summary
        summary(m.age.strength.post.unadj)    # model summary 
        confint(m.age.strength.post.unadj, level = 0.95, method = 'profile') 
        #plot(m.age.strength.post.unadj)       # check residuals
        exp(m.age.strength.post.unadj$coefficients)
        exp(confint(m.age.strength.post.unadj, level = 0.95, 
                    method = 'profile'))
        
        #plot(m.age.strength.post.unadj)       # check residuals
        # check for over/under dispersion; dispers ration >> 1 over; << 1 under
        check_overdispersion(m.age.strength.post.unadj)
        
      ## g) Unadjusted model: male age association with association degree
        # post manip
        m.age.degree.post.unadj <- glm.nb(degree.post ~ 
                                        Age.category,
                                      data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              #Color.manipulation. == 'n' &
                                              !is.na(degree.post)))
        
      ## h) model summary
        summary(m.age.degree.post.unadj)    # model summary 
        confint(m.age.degree.post.unadj, level = 0.95, method = 'profile')
        #plot(m.age.degree.post.unadj)       # check residuals
        exp(m.age.degree.post.unadj$coefficients)
        exp(confint(m.age.degree.post.unadj, level = 0.95, 
                    method = 'profile'))
        
        #plot(m.age.degree.post.unadj)       # check residuals
        # check for over/under dispersion; dispers ration >> 1 over; << 1 under
        check_overdispersion(m.age.degree.post.unadj)  
        
        
  ### 4.3 Female plumage trait associations with node level 
      # social network measures
    ## a) Unadjusted model: female breast brightness association with 
        # association strength pre-manip
      fem.Bbright.strength.pre.unadj <- glm.nb(strength.pre ~ 
                                            scale(R_avg.bright),
                                            data = subset(chr15_attrib_pre_df,
                                                   Sex == 'f' &
                                                   !is.na(strength.pre) &
                                                   !is.na(R_avg.bright)))
    
    ## b) model summary
      summary(fem.Bbright.strength.pre.unadj)    # model summary 
      confint(fem.Bbright.strength.pre.unadj, level = 0.95, method = 'profile')
      exp(fem.Bbright.strength.pre.unadj$coefficients)
      exp(confint(fem.Bbright.strength.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Bbright.strength.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Bbright.strength.pre.unadj)  
    
    ## c) Unadjusted model: female belly brightness association with 
        # association strength pre-manip
      fem.Bbright.strength.pre.unadj <- glm.nb(strength.pre ~ 
                                            scale(B_avg.bright),
                                            data = subset(chr15_attrib_pre_df,
                                                   Sex == 'f' &
                                                   !is.na(strength.pre) &
                                                   !is.na(B_avg.bright)))
    
    ## d) model summary
      summary(fem.Bbright.strength.pre.unadj)    # model summary 
      confint(fem.Bbright.strength.pre.unadj, level = 0.95, method = 'profile')
      exp(fem.Bbright.strength.pre.unadj$coefficients)
      exp(confint(fem.Bbright.strength.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Bbright.strength.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Bbright.strength.pre.unadj)  
    
    ## e) Unadjusted model: female breast brightness association 
        # with association degree pre-manip
      fem.Rbright.degree.pre.unadj <- glm.nb(degree.pre ~ 
                                      scale(R_avg.bright),
                                      data = subset(chr15_attrib_pre_df,
                                             Sex == 'f' &
                                             !is.na(degree.pre) &
                                             !is.na(R_avg.bright)))
    
    ## f) model summary
      summary(fem.Rbright.degree.pre.unadj)    # model summary 
      confint(fem.Rbright.degree.pre.unadj, level = 0.95, method = 'profile')
      exp(fem.Rbright.degree.pre.unadj$coefficients)
      exp(confint(fem.Rbright.degree.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Rbright.degree.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Rbright.degree.pre.unadj)  
    
    ## g) Unadjusted model: female belly brightness association 
        # with association degree pre-manip
      fem.Bbright.degree.pre.unadj <- glm.nb(degree.pre ~ 
                                          scale(B_avg.bright),
                                          data = subset(chr15_attrib_pre_df,
                                                 Sex == 'f' &
                                                 !is.na(degree.pre) &
                                                 !is.na(B_avg.bright)))
    
    ## h) model summary
      summary(fem.Bbright.degree.pre.unadj)    # model summary 
      confint(fem.Bbright.degree.pre.unadj, level = 0.95, method = 'profile')
      exp(fem.Bbright.degree.pre.unadj$coefficients)
      exp(confint(fem.Bbright.degree.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Bbright.degree.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Bbright.degree.pre.unadj)  
      
    ## i) Unadjusted model: female breast brightness association with 
      # association strength post manip
      fem.Bbright.strength.post.unadj <- glm.nb(strength.post ~ 
                                              scale(R_avg.bright),
                                            data = subset(chr15_attrib_post_df,
                                                    Sex == 'f' &
                                                    !is.na(strength.post) &
                                                    !is.na(R_avg.bright)))
      
    ## j) model summary
      summary(fem.Bbright.strength.post.unadj)    # model summary 
      confint(fem.Bbright.strength.post.unadj, level = 0.95, method = 'profile')
      exp(fem.Bbright.strength.post.unadj$coefficients)
      exp(confint(fem.Bbright.strength.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Bbright.strength.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Bbright.strength.post.unadj)  
      
    ## k) Unadjusted model: female belly brightness association with 
      # association strength post manip
      fem.Bbright.strength.post.unadj <- glm.nb(strength.post ~ 
                                              scale(B_avg.bright),
                                            data = subset(chr15_attrib_post_df,
                                                    Sex == 'f' &
                                                    !is.na(strength.post) &
                                                    !is.na(B_avg.bright)))
      
    ## l) model summary
      summary(fem.Bbright.strength.post.unadj)    # model summary 
      confint(fem.Bbright.strength.post.unadj, level = 0.95, method = 'profile')
      exp(fem.Bbright.strength.post.unadj$coefficients)
      exp(confint(fem.Bbright.strength.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Bbright.strength.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Bbright.strength.post.unadj)  
      
    ## m) Unadjusted model: female breast brightness association 
      # with association degree post manip
      fem.Rbright.degree.post.unadj <- glm.nb(degree.post ~ 
                                            scale(R_avg.bright),
                                          data = subset(chr15_attrib_post_df,
                                                  Sex == 'f' &
                                                  !is.na(degree.post) &
                                                  !is.na(R_avg.bright)))
      
    ## n) model summary
      summary(fem.Rbright.degree.post.unadj)    # model summary 
      confint(fem.Rbright.degree.post.unadj, level = 0.95, method = 'profile')
      exp(fem.Rbright.degree.post.unadj$coefficients)
      exp(confint(fem.Rbright.degree.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Rbright.degree.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Rbright.degree.post.unadj)  
      
    ## o) Unadjusted model: female belly brightness association 
      # with association degree post manip
      fem.Bbright.degree.post.unadj <- glm.nb(degree.post ~ 
                                            scale(B_avg.bright),
                                          data = subset(chr15_attrib_post_df,
                                                  Sex == 'f' &
                                                  !is.na(degree.post) &
                                                  !is.na(B_avg.bright)))
      
    ## p) model summary
      summary(fem.Bbright.degree.post.unadj)    # model summary 
      confint(fem.Bbright.degree.post.unadj, level = 0.95, method = 'profile')
      exp(fem.Bbright.degree.post.unadj$coefficients)
      exp(confint(fem.Bbright.degree.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(fem.Bbright.degree.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.Bbright.degree.post.unadj)  

    
  ### 4.4 Male plumage trait associations with node level 
    # social network measures
    ## a) Unadjusted model: male breast brightness association with 
        # association strength pre-manip
      m.Rbright.strength.pre.unadj <- glm.nb(strength.pre ~ 
                                          scale(R_avg.bright),
                                          data = subset(chr15_attrib_pre_df,
                                                 Sex == 'm' &
                                                 #Color.manipulation. == 'n' &
                                                 !is.na(strength.pre) &
                                                 !is.na(R_avg.bright)))
      
    ## b) model summary
      summary(m.Rbright.strength.pre.unadj)    # model summary 
      confint(m.Rbright.strength.pre.unadj, level = 0.95, method = 'profile')
      exp(m.Rbright.strength.pre.unadj$coefficients)
      exp(confint(m.Rbright.strength.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Rbright.strength.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Rbright.strength.pre.unadj)  
    
    ## c) Unadjusted model: male belly brightness association with 
        # association strength pre-manip
      m.Bbright.strength.pre.unadj <- glm.nb(strength.pre ~ 
                                          scale(B_avg.bright),
                                          data = subset(chr15_attrib_pre_df,
                                                 Sex == 'm' &
                                                 #Color.manipulation. == 'n' &
                                                 !is.na(strength.pre) &
                                                 !is.na(B_avg.bright)))
    
    ## d) model summary
      summary(m.Bbright.strength.pre.unadj)    # model summary 
      confint(m.Bbright.strength.pre.unadj, level = 0.95, method = 'profile')
      exp(m.Bbright.strength.pre.unadj$coefficients)
      exp(confint(m.Bbright.strength.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Bbright.strength.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Bbright.strength.pre.unadj)  
    
    ## e) Unadjusted model: male breast brightness association with 
        # association degree pre-manip
      m.Rbright.degree.pre.unadj <- glm.nb(degree.pre ~ 
                                        scale(R_avg.bright),
                                        data = subset(chr15_attrib_pre_df,
                                               Sex == 'm' &
                                               #Color.manipulation. == 'n' &
                                               !is.na(degree.pre)&
                                               !is.na(R_avg.bright)))
      
    ## f) model summary
      summary(m.Rbright.degree.pre.unadj)    # model summary 
      confint(m.Rbright.degree.pre.unadj, level = 0.95, method = 'profile')
      exp(m.Rbright.degree.pre.unadj$coefficients)
      exp(confint(m.Rbright.degree.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Rbright.degree.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Rbright.degree.pre.unadj)  
      
    ## g) Unadjusted model: male belly brightness association with 
        # association degree pre-manip
      m.Bbright.degree.pre.unadj <- glm.nb(degree.pre ~ 
                                        scale(B_avg.bright),
                                        data = subset(chr15_attrib_pre_df,
                                               Sex == 'm' &
                                               #Color.manipulation. == 'n' &
                                               !is.na(degree.pre) &
                                               !is.na(B_avg.bright)))
    
    ## h) model summary
      summary(m.Bbright.degree.pre.unadj)    # model summary 
      confint(m.Bbright.degree.pre.unadj, level = 0.95, method = 'profile')
      exp(m.Bbright.degree.pre.unadj$coefficients)
      exp(confint(m.Bbright.degree.pre.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Bbright.degree.pre.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Bbright.degree.pre.unadj)  
      
    ## i) Unadjusted model: male experimental breast brightness association with 
      # association strength post manip
      m.post.Rbright.strength.post.unadj <- glm.nb(strength.post ~ 
                                          scale(R.bright.treat.and.orig),
                                          data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              !is.na(strength.post) &
                                              !is.na(R.bright.treat.and.orig)))
      
    ## j) model summary
      summary(m.post.Rbright.strength.post.unadj)    # model summary 
      confint(m.post.Rbright.strength.post.unadj, level = 0.95, method = 'profile') 
      exp(m.post.Rbright.strength.post.unadj$coefficients)
      exp(confint(m.post.Rbright.strength.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.post.Rbright.strength.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.post.Rbright.strength.post.unadj) 
      
      
    ## k) Unadjusted model: male belly brightness association with 
      # association strength post manip
      m.Bbright.strength.post.unadj <- glm.nb(strength.post ~ 
                                            scale(B_avg.bright),
                                          data = subset(chr15_attrib_post_df,
                                                  Sex == 'm' &
                                                  #Color.manipulation. == 'n' &
                                                  !is.na(strength.post) &
                                                  !is.na(B_avg.bright)))
      
      ## l) model summary
      summary(m.Bbright.strength.post.unadj)    # model summary 
      confint(m.Bbright.strength.post.unadj, level = 0.95, method = 'profile')
      exp(m.Bbright.strength.post.unadj$coefficients)
      exp(confint(m.Bbright.strength.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Bbright.strength.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Bbright.strength.post.unadj)  
      
    ## m) Unadjusted model: male experimental breast brightness association with 
      # association degree post manip
      m.post.Rbright.degree.post.unadj <- glm.nb(degree.post ~ 
                                            scale(R.bright.treat.and.orig),
                                            data = subset(chr15_attrib_post_df,
                                              Sex == 'm' &
                                              !is.na(degree.post) &
                                              !is.na(R.bright.treat.and.orig)))
      
    ## n) model summary
      summary(m.post.Rbright.degree.post.unadj)    # model summary 
      confint(m.post.Rbright.degree.post.unadj, level = 0.95, method = 'profile') 
      exp(m.post.Rbright.degree.post.unadj$coefficients)
      exp(confint(m.post.Rbright.degree.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.post.Rbright.degree.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.post.Rbright.degree.post.unadj) 
      
      
    ## o) Unadjusted model: male belly brightness association with 
      # association degree post manip
      m.Bbright.degree.post.unadj <- glm.nb(degree.post ~ 
                                          scale(B_avg.bright),
                                        data = subset(chr15_attrib_post_df,
                                                Sex == 'm' &
                                                #Color.manipulation. == 'n' &
                                                !is.na(degree.post) &
                                                !is.na(B_avg.bright)))
      
    ## p) model summary
      summary(m.Bbright.degree.post.unadj)    # model summary 
      confint(m.Bbright.degree.post.unadj, level = 0.95, method = 'profile')
      exp(m.Bbright.degree.post.unadj$coefficients)
      exp(confint(m.Bbright.degree.post.unadj, level = 0.95, 
                  method = 'profile'))
      
      #plot(m.Bbright.degree.post.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.Bbright.degree.post.unadj) 
      

#***** NOTE
  # for the mediation model repro ~ age + strength what would our expectation be 
  # should we also look at the subset of individuals who are treated or 
  # an intx between treatment and breast brightness
      
###############################################################################
##############      5. Age association with plumage traits       ##############
###############################################################################         
      
  ### 5.1 Female age associations with plumage traits
    ## a) Unadjusted model: female age association with breast brightnes
      fem.age.Rbright.unadj <- lm(R_avg.bright ~ 
                                    Age.category,
                                    data = subset(chr15_attrib_pre_df,
                                          Sex == 'f' &
                                          !is.na(R_avg.bright)))
      
    ## b) model summary
      summary(fem.age.Rbright.unadj)    # model summary 
      confint(fem.age.Rbright.unadj, level = 0.95, method = 'profile')  
      #plot(base.gluc.wing.low.care.lmm)       # check residuals  
      
    ## c) Unadjusted model: female age association with breast brightnes
      fem.age.Bbright.unadj <- lm(B_avg.bright ~ 
                                  Age.category,
                                  data = subset(chr15_attrib_pre_df,
                                        Sex == 'f' &
                                        !is.na(B_avg.bright)))
      
    ## d) model summary
      summary(fem.age.Bbright.unadj)    # model summary 
      confint(fem.age.Bbright.unadj, level = 0.95, method = 'profile')  
      #plot(base.gluc.wing.low.care.lmm)       # check residuals  
      
      
  ### 5.2 Male age associations with plumage traits
    ## a) Unadjusted model: Male age association with breast brightnes
      m.age.Rbright.unadj <- lm(R_avg.bright ~ 
                                  Age.category,
                                  data = subset(chr15_attrib_pre_df,
                                        Sex == 'm' &
                                        #Color.manipulation. == 'n' &
                                        !is.na(R_avg.bright)))
      
    ## b) model summary
      summary(m.age.Rbright.unadj)    # model summary 
      confint(m.age.Rbright.unadj, level = 0.95, method = 'profile')  
      #plot(base.gluc.wing.low.care.lmm)       # check residuals  
      
    ## c) Unadjusted model: Male age association with breast brightnes
      m.age.Bbright.unadj <- lm(B_avg.bright ~ 
                                  Age.category,
                                  data = subset(chr15_attrib_pre_df,
                                        Sex == 'm' &
                                        #Color.manipulation. == 'n' &
                                        !is.na(B_avg.bright)))
      
    ## d) model summary
      summary(m.age.Bbright.unadj)    # model summary 
      confint(m.age.Bbright.unadj, level = 0.95, method = 'profile')  
      #plot(base.gluc.wing.low.care.lmm)       # check residuals  
      

      
###############################################################################
############  6. Soc. net. metrics association with reproduction  #############
###############################################################################        
      
  ### 6.1 Female social network measures associated with reproductive success
      ## a) Unadjusted model: female pre-manip strength association with 
          # total fecundity
        fem.strength.pre.tot.fecund.unadj <- glm.nb(total.fecundity ~ 
                                        scale(strength.pre),
                                        data = subset(chr15_attrib_pre_df,
                                              Sex == 'f' &
                                              !is.na(strength.pre) &
                                              !is.na(total.fecundity)))
      
      
    ## b) model summary
      summary(fem.strength.pre.tot.fecund.unadj)    # model summary 
      confint(fem.strength.pre.tot.fecund.unadj, 
              level = 0.95, method = 'profile')
      exp(fem.strength.pre.tot.fecund.unadj$coefficients)
      exp(confint(fem.strength.pre.tot.fecund.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(fem.strength.pre.tot.fecund.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.strength.pre.tot.fecund.unadj) 
      
      
    ## c) Unadjusted model: female pre-manip degree association with 
      # total fecundity
      fem.deg.pre.tot.fecund.unadj <- glm.nb(total.fecundity ~ 
                                    scale(degree.pre),
                                    data = subset(chr15_attrib_pre_df,
                                          Sex == 'f' &
                                          !is.na(degree.pre) &
                                          !is.na(total.fecundity)))
      
      
    ## d) model summary
      summary(fem.deg.pre.tot.fecund.unadj)    # model summary 
      confint(fem.deg.pre.tot.fecund.unadj, 
              level = 0.95, method = 'profile')
      exp(fem.deg.pre.tot.fecund.unadj$coefficients)
      exp(confint(fem.deg.pre.tot.fecund.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(fem.deg.pre.tot.fecund.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.deg.pre.tot.fecund.unadj) 
      
    ## e) Unadjusted model: female post manip strength association with 
      # total fecundity
      fem.strength.post.tot.fecund.unadj <- glm.nb(total.fecundity ~ 
                                          scale(strength.post),
                                          data = subset(chr15_attrib_post_df,
                                                Sex == 'f' &
                                                !is.na(strength.post) &
                                                !is.na(total.fecundity)))
      
      
    ## f) model summary
      summary(fem.strength.post.tot.fecund.unadj)    # model summary 
      confint(fem.strength.post.tot.fecund.unadj, 
              level = 0.95, method = 'profile')
      exp(fem.strength.post.tot.fecund.unadj$coefficients)
      exp(confint(fem.strength.post.tot.fecund.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(fem.strength.post.tot.fecund.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.strength.post.tot.fecund.unadj) 
      
    ## g) Unadjusted model: female post manip strength association with 
      # total fecundity
      fem.deg.post.tot.fecund.unadj <- glm.nb(total.fecundity ~ 
                                          scale(degree.post),
                                          data = subset(chr15_attrib_post_df,
                                                Sex == 'f' &
                                                !is.na(degree.post) &
                                                !is.na(total.fecundity)))
      
      
    ## h) model summary
      summary(fem.deg.post.tot.fecund.unadj)    # model summary 
      confint(fem.deg.post.tot.fecund.unadj, 
              level = 0.95, method = 'profile')
      exp(fem.deg.post.tot.fecund.unadj$coefficients)
      exp(confint(fem.deg.post.tot.fecund.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(fem.deg.post.tot.fecund.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(fem.deg.post.tot.fecund.unadj) 
      
      
  ### 6.2 Male social network measures associated with reproductive success
    ## a) Unadjusted model: male pre-manip strength association with 
      # clutch 1 total paternity
      m.strength.pre.clutch.1.pat.unadj <- glm.nb(attmpt.1.tot.pat ~ 
                                          scale(strength.pre),
                                          data = subset(chr15_attrib_pre_df,
                                                Sex == 'm' &
                                                !is.na(strength.pre) &
                                                !is.na(attmpt.1.tot.pat)))
      
      
    ## b) model summary
      summary(m.strength.pre.clutch.1.pat.unadj)    # model summary 
      confint(m.strength.pre.clutch.1.pat.unadj, 
              level = 0.95, method = 'profile')
      exp(m.strength.pre.clutch.1.pat.unadj$coefficients)
      exp(confint(m.strength.pre.clutch.1.pat.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(m.strength.pre.clutch.1.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.strength.pre.clutch.1.pat.unadj) 
      
      
    ## c) Unadjusted model: male pre-manip degree association with 
      # clutch 1 total paternity
      m.deg.pre.clutch.1.pat.unadj <- glm.nb(attmpt.1.tot.pat ~ 
                                    scale(degree.pre),
                                    data = subset(chr15_attrib_pre_df,
                                          Sex == 'm' &
                                          !is.na(degree.pre) &
                                          !is.na(attmpt.1.tot.pat)))
      
      
      ## d) model summary
      summary(m.deg.pre.clutch.1.pat.unadj)    # model summary 
      confint(m.deg.pre.clutch.1.pat.unadj, 
              level = 0.95, method = 'profile')
      exp(m.deg.pre.clutch.1.pat.unadj$coefficients)
      exp(confint(m.deg.pre.clutch.1.pat.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(m.deg.pre.clutch.1.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.deg.pre.clutch.1.pat.unadj) 
      
      
      
    ## e) Unadjusted model: male post manip strength association with 
      # clutch 2 total paternity
      m.strength.post.clutch.2.pat.unadj <- glm.nb(attmpt.2.tot.pat ~ 
                                          scale(strength.post),
                                          data = subset(chr15_attrib_post_df,
                                                Sex == 'm' &
                                                !is.na(strength.post) &
                                                !is.na(attmpt.1.tot.pat)))
      
      
      ## f) model summary
      summary(m.strength.post.clutch.2.pat.unadj)    # model summary 
      confint(m.strength.post.clutch.2.pat.unadj, 
              level = 0.95, method = 'profile')
      exp(m.strength.post.clutch.2.pat.unadj$coefficients)
      exp(confint(m.strength.post.clutch.2.pat.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(m.strength.post.clutch.2.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.strength.post.clutch.2.pat.unadj) 
      
      
    ## g) Unadjusted model: male post manip strength association with 
      # clutch 2 total paternity
      m.deg.post.clutch.2.pat.unadj <- glm.nb(attmpt.2.tot.pat ~ 
                                      scale(degree.post),
                                      data = subset(chr15_attrib_post_df,
                                             Sex == 'm' &
                                             !is.na(degree.post) &
                                             !is.na(attmpt.1.tot.pat)))
      
      
    ## h) model summary
      summary(m.deg.post.clutch.2.pat.unadj)    # model summary 
      confint(m.deg.post.clutch.2.pat.unadj, 
              level = 0.95, method = 'profile')
      exp(m.deg.post.clutch.2.pat.unadj$coefficients)
      exp(confint(m.deg.post.clutch.2.pat.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(m.deg.post.clutch.2.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.deg.post.clutch.2.pat.unadj) 
      
    ## i) Unadjusted model: male pre-manip strength association with 
      # all cluches total paternity
      m.strength.pre.tot.pat.unadj <- glm.nb(total.paternity ~ 
                                     scale(strength.pre),
                                     data = subset(chr15_attrib_pre_df,
                                          Sex == 'm' &
                                          !is.na(strength.pre) &
                                          !is.na(total.paternity)))
      
      
    ## j) model summary
      summary(m.strength.pre.tot.pat.unadj)    # model summary 
      confint(m.strength.pre.tot.pat.unadj, 
              level = 0.95, method = 'profile')
      exp(m.strength.pre.tot.pat.unadj$coefficients)
      exp(confint(m.strength.pre.tot.pat.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(m.strength.pre.tot.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.strength.pre.tot.pat.unadj) 
      
    ## k) Unadjusted model: male pre-manip degree association with 
      # all clutches total paternity
      m.deg.pre.tot.pat.unadj <- glm.nb(total.paternity ~ 
                                     scale(degree.pre),
                                     data = subset(chr15_attrib_pre_df,
                                          Sex == 'm' &
                                          !is.na(degree.pre) &
                                          !is.na(total.paternity)))
      
      
    ## l) model summary
      summary(m.deg.pre.tot.pat.unadj)    # model summary 
      confint(m.deg.pre.tot.pat.unadj, 
              level = 0.95, method = 'profile')
      exp(m.deg.pre.tot.pat.unadj$coefficients)
      exp(confint(m.deg.pre.tot.pat.unadj, 
                  level = 0.95, method = 'profile'))
      
      #plot(m.deg.pre.tot.pat.unadj)       # check residuals
      # check for over/under dispersion; dispers ration >> 1 over; << 1 under
      check_overdispersion(m.deg.pre.tot.pat.unadj) 
      
      
      
      
    
###############################################################################
##############      7. Age by soc. net. metrics interaction      ##############
###############################################################################        
     
  ### 7.1 Interaction models: Female trait by node level social network measures
        # associations with reproductive success
      
    ## a) Interaction model: female age by pre-manip strength association 
      # with total fecundity
      fem.age.x.pre.strength.fecund.intx <- glm.nb(total.fecundity ~ 
                                      Age.category * scale(strength.pre),
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'f' &
                                            !is.na(strength.pre) &
                                            !is.na(total.fecundity)))
        
    ## b) model summary
      summary(fem.age.x.pre.strength.fecund.intx)    # model summary 
      #plot(fem.age.x.pre.strength.fecund.intx)       # check residuals
      
     
    ## c) Interaction model: female age by post-manip degree association 
      # with total fecundity
      fem.age.x.post.degree.fecund.intx <- glm.nb(total.fecundity ~ 
                                    Age.category * scale(strength.post),
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
      m.age.x.pre.strength.cluth.1.intx <- glm.nb(attmpt.1.tot.pat ~ 
                                  Age.category * scale(strength.pre),
                                  data = subset(chr15_attrib_pre_df,
                                        Sex == 'm' &
                                        !is.na(strength.pre) &
                                        !is.na(attmpt.1.tot.pat)))
      
    ## b) model summary
      summary(m.age.x.pre.strength.cluth.1.intx)    # model summary 
      #plot(m.age.x.pre.strength.cluth.1.intx)       # check residuals
      class(m.age.x.pre.strength.cluth.1.intx)
    
      
      
    ## c) Interaction model: male age by pre manip strength association 
      # with all clutches total paternity  
      m.age.x.pre.strength.tot.intx <- glm.nb(total.paternity ~ 
                                  Age.category * scale(strength.pre),
                                  data = subset(chr15_attrib_pre_df,
                                         Sex == 'm' &
                                         !is.na(strength.pre) &
                                         !is.na(total.paternity)))
      
    ## d) model summary
      summary(m.age.x.pre.strength.tot.intx)    # model summary 
      #plot(m.age.x.pre.strength.tot.intx)       # check residuals
      
    ## e) Interaction model: male age by post manip degree association 
      # with clutch 2 total paternity  
      m.age.x.post.degree.cluth.2.intx <- glm.nb(attmpt.2.tot.pat ~ 
                                      Age.category * scale(degree.post),
                                      data = subset(chr15_attrib_post_df,
                                            Sex == 'm' &
                                            !is.na(degree.post) &
                                            !is.na(attmpt.2.tot.pat)))
      
    ## f) model summary
      summary(m.age.x.post.degree.cluth.2.intx)    # model summary 
      #plot(m.age.x.post.degree.cluth.2.intx)       # check residuals

      
      
 
###############################################################################
##############              7. Effects decomposition             ##############
###############################################################################  

      
  ### 7.1 Meditation models: Female age/trait plus node level social network 
      # measures associations with reproductive success
      
    ## a) Female mediation model: Female age plus strength association 
      # with total fecundity
      fem.age.pre.strength.tot.med <- glm.nb(total.fecundity ~ 
                                      Age.category + scale(strength.pre),
                                      data = subset(chr15_attrib_pre_df,
                                            Sex == 'f' &
                                            !is.na(strength.pre) &
                                            !is.na(total.fecundity)))
      
    ## b) model summary
      summary(fem.age.pre.strength.tot.med)    # model summary 
      confint(fem.age.pre.strength.tot.med, level = 0.95, method = 'profile')  
      #plot(base.gluc.wing.low.care.lmm)       # check residuals
      
      exp(fem.age.pre.strength.tot.med$coefficients)
      exp(confint(fem.age.pre.strength.tot.med, 
                  level = 0.95, method = 'profile'))
      
      
# ********************************************* NOT WORKING
     # test<- mediate_zi(fem.age.strength.pre.unadj, fem.age.pre.strength.tot.med, 
     #             sims = 1000, boot = FALSE, 
     #             treat = "Age.category", mediator = "scale(strength.pre)", 
     #             covariates = NULL, outcome = NULL, control = NULL, 
     #             conf.level = 0.95, control.value = 0, treat.value = 1, 
     #             long = TRUE, dropobs = FALSE, robustSE = FALSE, 
     #             cluster = NULL)
     # 
     # y <- data$total.fecundity
     # z <- data$Age.category
     # m <- data$scale(strength.pre)
     # 
     # test <- mediate_iv(y, z, m, 
     #                    ydist = 'negbin', mtype = 'continuous',
     #                    x.IV, tol = 0.5, 
     #                    n.init = 1, control = list(maxit = 15, ftol = 0.5, 
     #                                               gtol = 0.5, trace = FALSE), 
     #                    sims = 3)
     # 
     # Error in mediate(fem.age.strength.pre.unadj, fem.age.strength.tot.med,  : 
     #                    unsupported glm family
     #  
     #  #average causal mediation effects (ACME)
     #  #average direct effect (ADE)
     #  
     #  med.out <- mediate(med.fit, out.fit, treat = "treat", 
     #                     mediator = "emo", + robustSE = TRUE, sims = 100) 
     #  summary(med.out)     
     #  
     #  
     #  med.out <- mediate(fem.age.strength.pre.unadj, fem.age.strength.tot.med, 
     #                     treat = "Age.category", 
     #                     mediator = "scale(strength.pre)" #, + robustSE = T
     #                     , sims = 100) 
     #  
     #  # using sy and asy as control and treatment, respectively
     #  summary(med.out)     
     #  summary(fem.age.Rbright.unadj)
     #  (fem.age.Rbright.tot.med)
# ********************************************* NOT WORKING      
      
    ## c) Female mediation model: Female age plus post manipulation degree  
      # association with total fecundity
      fem.age.post.degree.tot.med <- glm.nb(total.fecundity ~ 
                                Age.category + scale(degree.post),
                                data = subset(chr15_attrib_post_df,
                                       Sex == 'f' &
                                       !is.na(degree.post) &
                                       !is.na(total.fecundity)))
      
    ## b) model summary
      summary(fem.age.post.degree.tot.med)    # model summary 
      confint(fem.age.post.degree.tot.med, level = 0.95, method = 'profile')  
      #plot(fem.age.post.degree.tot.med)       # check residuals
      
      exp(fem.age.post.degree.tot.med$coefficients)
      exp(confint(fem.age.post.degree.tot.med, 
                  level = 0.95, method = 'profile')) 
      

      
  ### 7.2 Mediation models: Male age plus node level social network 
      # measures associations with reproductive success
      
    ## a) Male mediation model: male age plus pre-manip strength association 
      # with clutch 1 total paternity
      m.age.pre.strength.cluth.1.med <- glm.nb(attmpt.1.tot.pat ~ 
                                     Age.category + scale(strength.pre),
                                     data = subset(chr15_attrib_pre_df,
                                            Sex == 'm' &
                                            !is.na(strength.pre) &
                                            !is.na(attmpt.1.tot.pat)))
      
    ## b) model summary
      summary(m.age.pre.strength.cluth.1.med)    # model summary 
      confint(m.age.pre.strength.cluth.1.med, level = 0.95, method = 'profile')  
      #plot(m.age.pre.strength.cluth.1.med)       # check residuals
      
      exp(m.age.pre.strength.cluth.1.med$coefficients)
      exp(confint(m.age.pre.strength.cluth.1.med, 
                  level = 0.95, method = 'profile'))
      
      
    ## c) Male mediation model: male age plus pre manip strength association 
      # with all clutches total paternity  
      m.age.pre.strength.tot.med <- glm.nb(total.paternity ~ 
                                  Age.category + scale(strength.pre)
                                  # mediator - outcome confounder
                                  + Color.manipulation.
                                  ,
                                  data = subset(chr15_attrib_pre_df,
                                         Sex == 'm' &
                                         !is.na(strength.pre) &
                                         !is.na(total.paternity)))
      
    ## d) model summary
      summary(m.age.pre.strength.tot.med)    # model summary 
      confint(m.age.pre.strength.tot.med, level = 0.95, method = 'profile')  
      #plot(m.age.pre.strength.tot.med)       # check residuals
      
      exp(m.age.pre.strength.tot.med$coefficients)
      exp(confint(m.age.pre.strength.tot.med, 
                  level = 0.95, method = 'profile')) 
      
      
    ## e) Male mediation model: male age plus post manip degree association 
      # with clutch 2 total paternity  
      m.age.post.degree.cluth.2.med <- glm.nb(attmpt.2.tot.pat ~ 
                                   Age.category + scale(degree.post)
                                   # mediator - outcome confounder
                                   + Color.manipulation.
                                   ,
                                   data = subset(chr15_attrib_post_df,
                                          Sex == 'm' &
                                          !is.na(degree.post) &
                                          !is.na(attmpt.2.tot.pat)))
      
    ## f) model summary
      summary(m.age.post.degree.cluth.2.med)    # model summary 
      confint(m.age.post.degree.cluth.2.med, level = 0.95, method = 'profile')  
      #plot(m.age.post.degree.cluth.2.med)       # check residuals
      
      exp(m.age.post.degree.cluth.2.med$coefficients)
      exp(confint(m.age.post.degree.cluth.2.med, 
                  level = 0.95, method = 'profile')) 
      
     
   