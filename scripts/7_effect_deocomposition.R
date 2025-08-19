################################################################################
#############       Testing the mediating role of female-male      #############  
#############    social interactions on the relationship between   #############
#############             age and reproductive success.            #############
#############                                                      #############
#############                7. Effect decomposition               #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 29 Sept 2024                #############
#############               last updated: 6 Dec 2024               #############
################################################################################



### PURPOSE: Decompose the total effect into the direct and indirect effects

  # Code Blocks
    # 1. Configure work space
    # 2. Load RData
    # 3. Female effect decomposition models
    # 4. Male effect decomposition models
    # 5. Tidy effect decomposition
    # 6. Graph effect decomposition

    

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
      
      # load patchwork packages
        library('patchwork')
   
      
    ## c) Modeling Packages
      
      # load mediation
        library('mediation')
      
      # library('broom')
        library('broom.mixed')
      
        
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
    ## a) Load RData mediation data
      load(here('data/7_chr_effect_decomp_data.RData'))
    
    ## b) Load RData mediation models
      load(here('data/7_chr_effect_decomp_mods.RData'))
    
      

###############################################################################
##############       3. Female effect decomposition models       ##############
###############################################################################  
      # age --> strength(fxm) --> repro
      # TS --> strength(fxm) --> repro - age adjusted
      # TS --> degree(fxm)  --> repro - age adjusted
      
  ### 3.1 Female mediation model: female age + pre-manip strength
      # with total fecundity.
      # exposure = age
      # mediators = fxm strength
      # outcomes = total fecundity (females)
      
      # Estimate the average causal mediation effects (ACME) and the 
      # average direct effect (ADE); 'mediation' package
      # mediation model = f_fxm_age_strength; Script 6, section 5.1
      # outcome model = f_fxm_age_strength_fecund
      
    ## a) Female Age - Strength fxm full model
      f_fxm_age_strength_fecund <- lm(total.fecundity ~ 
                                          Age.category + strength.fxm,
                                          # test for intx
                                          #Age.category * strength.fxm,
                                          #family = 'gaussian',
                                          data = subset(chr15_attrib_df,
                                                Sex == 'f' &
                                                !is.na(Age.category) &
                                                !is.na(strength.fxm) &
                                                !is.na(total.fecundity)))
    ## b) full model summary
      summary(f_fxm_age_strength_fecund)
      
    ## c) Estimate average causal mediation effects
      med_out_f_age_strength_fecund <- mediate(f_fxm_age_strength, 
                                        f_fxm_age_strength_fecund, 
                                      treat = "Age.category", 
                                      mediator = "strength.fxm",
                                      robustSE = TRUE, sims = 1000) 
      
    ## d) test for exposure*mediator intx
      # test.TMint(med_out_f_age_strength_fecund,
      #            conf.level = .95)
      
    ## e) Average causal mediation model summary
      summary(med_out_f_age_strength_fecund)
      
    ## f) sensitivity analysis for unmeasured confounder between mediator and 
       # outcome
      sens_f_age_strength_fecund <- medsens(med_out_f_age_strength_fecund,
                          rho.by = 0.1, effect.type = 'indirect',
                          sims = 1000)
      summary(sens_f_age_strength_fecund)
      plot(sens_f_age_strength_fecund)
    
      
##*** MODEL FITTING TEST using 'medflex' package
    # ## g) Expanded data and weights for female pre-manip strength by age
    #   # 'expand the data' by estimating mediator for different levels of the 
    #   # exposure to generate weights; implemented in medflex 
    #   exp_wght_f_age_strength <- neWeight(f_fxm_age_strength)
    #   weights(exp_wght_f_age_strength)
    #   
    # # h) Ne_impute method; implemented in medflex using full model
    #   # f_fxm_age_strength_fecund 
    #   exp_imp_f_age_strength_fecund <- neImpute(f_fxm_age_strength_fecund)
    # 
    # 
    # ## i) Natural effects model female strength by age 
    #   # calculate bootstrapped SE based resampling original data 
    #   # with replacement (1000 iterations); implemented in medflex 
    #   
    #   # Weighting method   
    #   ne_wght_f_age_strength_fecund <- neModel(total.fecundity ~ 
    #                                   Age.category0 + Age.category1,
    #                                   family = gaussian, 
    #                                   expData = exp_wght_f_age_strength)
    #   
    #   summary(ne_wght_f_age_strength_fecund)
    #   
    #   # Imputation method  
    #   ne_imp_f_age_strength_fecund <- neModel(total.fecundity ~ 
    #                                   Age.category0 + Age.category1,
    #                                   family = gaussian, 
    #                                   expData = exp_wght_f_age_strength)
    #   
    #   summary(ne.imp.fem.age.strength.pre.fecund)
    # 
    # ## j) Decomposition of the natural effects based on weights method:
    #   
    #   # Weighting method   
    #   eff.decomp.fem.age.strength.pre.fecund.wght <- 
    #     neEffdecomp(ne.wght.fem.age.strength.pre.fecund)
    #   
    #   summary(eff.decomp.fem.age.strength.pre.fecund.wght)
    #   confint(eff.decomp.fem.age.strength.pre.fecund.wght)
    #   
    #   # Imputation method  
    #   eff.decomp.fem.age.strength.pre.fecund.imp <- 
    #     neEffdecomp(ne.imp.fem.age.strength.pre.fecund)
    #   
    #   summary(eff.decomp.fem.age.strength.pre.fecund.imp)
    #   confint(eff.decomp.fem.age.strength.pre.fecund.imp)

      
  ### 3.2 Female mediation model: female tail streamer length + strength/degree
        # with total fecundity - age addujsted.
        # exposure = age
        # mediators = fxm strength / degree
        # outcomes = total fecundity (females)  
      
        # Estimate the average causal mediation effects (ACME) and the 
        # average direct effect (ADE); 'mediation' package
        # mediation models = f_fxm_TS_age_strength / f_fxm_age_degree; 
                              # Script 6, section 5.3
        # outcome model = f_fxm_TS_age_strength_fecund / 
                      #   f_fxm_TS_age_degree_fecund
        
    ## a) Female Tail Streamer - Strength fxm full model
      # both full and mediator models adjusted for age
      f_fxm_TS_age_strength_fecund <- lm(total.fecundity ~ 
                                      scale(Mean.TS) + strength.fxm 
                                    # adjust for age as confounder
                                      + Age.category ,
                                    # test for intx
                                      #scale(Mean.TS) * strength.fxm,
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_df,
                                            Sex == 'f' &
                                            !is.na(Mean.TS) &
                                            !is.na(strength.fxm) &
                                            !is.na(Age.category) &
                                            !is.na(total.fecundity)))
    ## b) full model summary
      summary(f_fxm_TS_age_strength_fecund)
      
    ## c) Estimate average causal mediation effects
      med_out_f_TS_age_strength_fecund <- mediate(f_fxm_TS_age_strength, 
                                               f_fxm_TS_age_strength_fecund, 
                                               treat = "scale(Mean.TS)", 
                                               mediator = "strength.fxm",
                                               robustSE = TRUE, sims = 1000) 
      
    ## d) test for exposure*mediator intx
      # test.TMint(med_out_f_TS_age_strength_fecund,
      #            conf.level = .95)
      
    ## e) Average causal mediation model summary
      summary(med_out_f_TS_age_strength_fecund)
      
    ## f) sensitivity analysis for unmeasured confounder between mediator and 
      # outcome
      sens_f_TS_age_strength_fecund <- 
        medsens(med_out_f_TS_age_strength_fecund,
                rho.by = 0.1, effect.type = 'indirect',
                sims = 1000)
      
      summary(sens_f_TS_age_strength_fecund)
      plot(sens_f_TS_age_strength_fecund)
      
    ## g) Female Tail Streamer - Degree fxm full model 
      # both full and mediator models adjusted for age
      f_fxm_TS_age_degree_fecund <- lm(total.fecundity ~ 
                                      scale(Mean.TS) + degree.fxm 
                                    # adjust for age as confounder
                                      + Age.category ,
                                    # test for intx
                                      #scale(Mean.TS) * degree.fxm,
                                      #family = 'gaussian',
                                      data = subset(chr15_attrib_df,
                                            Sex == 'f' &
                                            !is.na(Mean.TS) &
                                            !is.na(degree.fxm) &
                                            !is.na(Age.category) &
                                            !is.na(total.fecundity)))
    ## h) full model summary
      summary(f_fxm_TS_age_degree_fecund)
      
    ## i) Estimate average causal mediation effects
      med_out_f_TS_age_degree_fecund <- mediate(f_fxm_TS_age_degree, 
                                            f_fxm_TS_age_degree_fecund, 
                                                treat = "scale(Mean.TS)", 
                                                mediator = "degree.fxm",
                                                robustSE = TRUE, sims = 1000) 
      
    ## j) test for exposure*mediator intx
      # test.TMint(med_out_f_TS_age_degree_fecund,
      #            conf.level = .95)
      
    ## k) Average causal mediation model summary
      summary(med_out_f_TS_age_degree_fecund)
      
    ## l) sensitivity analysis for unmeasured confounder between mediator and 
      # outcome
      sens_f_TS_age_degree_fecund <- 
        medsens(med_out_f_TS_age_degree_fecund,
                rho.by = 0.1, effect.type = 'indirect',
                sims = 1000)
      
      summary(sens_f_TS_age_degree_fecund)
      plot(sens_f_TS_age_degree_fecund)
      
      

          
###############################################################################
##############        4. Male models effects decomposition         ############
###############################################################################
      

  ### 4.1 Male mediation model: male age + pre-manip strength
        # with first clutch total paternity.
        # exposure = age
        # mediators = pre-manip strength
        # outcomes = clutch 1 total paternity 
        
       # Estimate the average causal mediation effects (ACME) and the 
        # average direct effect (ADE); 'mediation' package
        # mediation model = m_fxm_age_strength; Script 6, section 5.2
        # outcome model = m_fxm_age_strength_tot_pat
        
    ## a) Male Age - Strength fxm full model
      m_fxm_age_strength_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                  Age.category + strength.fxm,
                                # test for intx
                                   #Age.category * strength.fxm,
                                   #family = 'gaussian',
                                  data = subset(chr15_attrib_df,
                                        Sex == 'm' &
                                        !is.na(Age.category) &
                                        !is.na(strength.fxm) &
                                        !is.na(nest.2.3.tot.pat)))
      
    ## b) full model summary
      summary(m_fxm_age_strength_tot_pat)
      
    ## c) Estimate average causal mediation effects
      med_out_m_age_strength_tot_pat <- mediate(m_fxm_age_strength, 
                                               m_fxm_age_strength_tot_pat, 
                                               treat = "Age.category", 
                                               mediator = "strength.fxm",
                                               robustSE = TRUE, sims = 1000) 
      
    ## d) test for exposure*mediator intx
      # test.TMint(med_out_m_age_strength_tot_pat,
      #            conf.level = .95)
      
    ## e) Average causal mediation model summary
      summary(med_out_m_age_strength_tot_pat)
      
    ## f) sensitivity analysis for unmeasured confounder between mediator and 
      # outcome
      sens_m_age_strength_tot_pat <- 
        medsens(med_out_m_age_strength_tot_pat,
                rho.by = 0.1, effect.type = 'indirect',
                sims = 1000)
      
      summary(sens_m_age_strength_tot_pat)
      plot(sens_m_age_strength_tot_pat)
    
 
    ### 4.2 Male mediation model: male tail streamer length + strength/degree
      # with total fecundity - age addujsted.
      # exposure = age
      # mediators = fxm strength / degree
      # outcomes = total fecundity (females)  
      
      # Estimate the average causal mediation effects (ACME) and the 
      # average direct effect (ADE); 'mediation' package
      # mediation models = f_fxm_TS_age_strength / f_fxm_age_degree; 
      # Script 6, section 5.5
      # outcome model = f_fxm_TS_age_strength_fecund / 
      #   f_fxm_TS_age_degree_fecund
      
    ## a) Male Tail Streamer - Strength fxm full model
      # both full and mediator models adjusted for age
      m_fxm_TS_age_strength_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                     scale(Mean.TS) + strength.fxm 
                                   # adjust for age as confounder
                                     + Age.category ,
                                   # test for intx
                                     #scale(Mean.TS) * strength.fxm,
                                     #family = 'gaussian',
                                     data = subset(chr15_attrib_df,
                                           Sex == 'm' &
                                           !is.na(Mean.TS) &
                                           !is.na(strength.fxm) &
                                           !is.na(Age.category) &
                                           !is.na(nest.2.3.tot.pat)))
      
    ## b) full model summary
      summary(m_fxm_TS_age_strength_tot_pat)
      
    ## c) Estimate average causal mediation effects
      med_out_m_TS_age_strength_tot_pat <- mediate(m_fxm_TS_age_strength, 
                                            m_fxm_TS_age_strength_tot_pat, 
                                            treat = "scale(Mean.TS)", 
                                            mediator = "strength.fxm",
                                            robustSE = TRUE, sims = 1000) 
      
    ## d) test for exposure*mediator intx
      # test.TMint(med_out_m_TS_age_strength_tot_pat,
      #            conf.level = .95)
      
    ## e) Average causal mediation model summary
      summary(med_out_m_TS_age_strength_tot_pat)
      
    ## f) sensitivity analysis for unmeasured confounder between mediator and 
      # outcome
      sens_m_TS_age_strength_tot_pat <- 
        medsens(med_out_m_TS_age_strength_tot_pat,
                rho.by = 0.1, effect.type = 'indirect',
                sims = 1000)
      
      summary(sens_m_TS_age_strength_tot_pat)
      plot(sens_m_TS_age_strength_tot_pat)
      
    ## g) Male Tail Streamer - Degree fxm full model 
      # both full and mediator models adjusted for age
      m_fxm_TS_age_degree_tot_pat <- lm(nest.2.3.tot.pat ~ 
                                   scale(Mean.TS) + degree.fxm 
                                 # adjust for age as confounder
                                   + Age.category ,
                                 # test for intx
                                   #scale(Mean.TS) * degree.fxm,
                                   #family = 'gaussian',
                                   data = subset(chr15_attrib_df,
                                         Sex == 'm' &
                                         !is.na(Mean.TS) &
                                         !is.na(degree.fxm) &
                                         !is.na(Age.category) &
                                         !is.na(nest.2.3.tot.pat)))
  
    ## h) full model summary
      summary(m_fxm_TS_age_degree_tot_pat)
      
    ## i) Estimate average causal mediation effects
      med_out_m_TS_age_degree_tot_pat <- mediate(m_fxm_TS_age_degree, 
                                            m_fxm_TS_age_degree_tot_pat, 
                                            treat = "scale(Mean.TS)", 
                                            mediator = "degree.fxm",
                                            robustSE = TRUE, sims = 1000) 
      
    ## j) test for exposure*mediator intx
      # test.TMint(med_out_m_TS_age_degree_tot_pat,
      #            conf.level = .95)
      
    ## k) Average causal mediation model summary
      summary(med_out_m_TS_age_degree_tot_pat)
      
    ## l) sensitivity analysis for unmeasured confounder between mediator and 
      # outcome
      sens_m_TS_age_degree_tot_pat <- 
        medsens(med_out_m_TS_age_degree_tot_pat,
                rho.by = 0.1, effect.type = 'indirect',
                sims = 1000)
      
      summary(sens_m_TS_age_degree_tot_pat)
      plot(sens_m_TS_age_degree_tot_pat)
      
      
      
###############################################################################
##############            5. Tidy effects decomposition            ############
###############################################################################

  ### 5.1 Set standard values for the tidy tibbles
    ## a) a list of effect types
      `Effect type` <- c('indirect', 'direct', 'total', 'Prop. Mediated')
     

  ### 5.2 Tidy female effects decomposition: female age + strength
      # with total fecundity.
    ## a) Assign to generic model
      model <- med_out_f_age_strength_fecund
    
    ## b) create a vector of estimates
      `Estimate` <- c(model$d.avg, model$z.avg, model$tau.coef, model$n.avg) 
      
    ## c) create a vector of lower 95% CI
      `95% CI Lower` <- c(model$d.avg.ci['2.5%'], model$z.avg.ci['2.5%'], 
                          model$tau.ci['2.5%'], model$n.avg.ci['2.5%']) 
      
    ## d) create a vector of upper 95% CI
      `95% CI Upper` <- c(model$d.avg.ci['97.5%'], model$z.avg.ci['97.5%'], 
                          model$tau.ci['97.5%'], model$n.avg.ci['97.5%'])
    
    ## e) create a vector of p-values
      `p-value` <- c(model$d.avg.p, model$z.avg.p, model$tau.p, 
                     model$n.avg.p)
    
    ## f) Combine vectors into a tibble
      tidy_f_age_strength_fecund <- tibble(`Effect type`, `Estimate`,
                                                `95% CI Lower`, `95% CI Upper`,
                                                `p-value`)
  
    ## g) Label the estimates in data frame by sex
      tidy_f_age_strength_fecund$sex <- c('female')
      
    ## h) Label the estimates in data frame by sex
      tidy_f_age_strength_fecund$model <- c('age-strength-fecundity')
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy_f_age_strength_fecund <- 
        transform(tidy_f_age_strength_fecund, 
                  `Effect type` = factor(`Effect type`,
                                  levels = c('total','direct','indirect', 
                                             'Prop. Mediated')))

      
  ### 5.3 Tidy female effects decomposition: female tail streamer length +
      # strength with total fecundity - age adjusted
    ## a) Assign to generic model
      model <- med_out_f_TS_age_strength_fecund
      
    ## b) create a vector of estimates
      `Estimate` <- c(model$d.avg, model$z.avg, model$tau.coef, model$n.avg) 
      
    ## c) create a vector of lower 95% CI
      `95% CI Lower` <- c(model$d.avg.ci['2.5%'], model$z.avg.ci['2.5%'], 
                          model$tau.ci['2.5%'], model$n.avg.ci['2.5%']) 
      
    ## d) create a vector of upper 95% CI
      `95% CI Upper` <- c(model$d.avg.ci['97.5%'], model$z.avg.ci['97.5%'], 
                          model$tau.ci['97.5%'], model$n.avg.ci['97.5%'])
      
    ## e) create a vector of p-values
      `p-value` <- c(model$d.avg.p, model$z.avg.p, model$tau.p, 
                     model$n.avg.p)
      
    ## f) Combine vectors into a tibble
      tidy_f_TS_age_strength_fecund <- tibble(`Effect type`, `Estimate`,
                                           `95% CI Lower`, `95% CI Upper`,
                                           `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy_f_TS_age_strength_fecund$sex <- c('female')
      
    ## h) Label the estimates in data frame by sex
      tidy_f_TS_age_strength_fecund$model <- c('TS-strength-fecundity')
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy_f_TS_age_strength_fecund <- 
        transform(tidy_f_TS_age_strength_fecund, 
                  `Effect type` = factor(`Effect type`,
                                     levels = c('total','direct','indirect', 
                                                'Prop. Mediated')))
      
      
  ### 5.4 Tidy female effects decomposition: female tail streamer length +
      # degree with total fecundity - age adjusted
    ## a) Assign to generic model
      model <- med_out_f_TS_age_degree_fecund
      
    ## b) create a vector of estimates
      `Estimate` <- c(model$d.avg, model$z.avg, model$tau.coef, model$n.avg) 
      
    ## c) create a vector of lower 95% CI
      `95% CI Lower` <- c(model$d.avg.ci['2.5%'], model$z.avg.ci['2.5%'], 
                          model$tau.ci['2.5%'], model$n.avg.ci['2.5%']) 
      
    ## d) create a vector of upper 95% CI
      `95% CI Upper` <- c(model$d.avg.ci['97.5%'], model$z.avg.ci['97.5%'], 
                          model$tau.ci['97.5%'], model$n.avg.ci['97.5%'])
      
    ## e) create a vector of p-values
      `p-value` <- c(model$d.avg.p, model$z.avg.p, model$tau.p, 
                     model$n.avg.p)
      
    ## f) Combine vectors into a tibble
      tidy_f_TS_age_degree_fecund <- tibble(`Effect type`, `Estimate`,
                                              `95% CI Lower`, `95% CI Upper`,
                                              `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy_f_TS_age_degree_fecund$sex <- c('female')
      
    ## h) Label the estimates in data frame by sex
      tidy_f_TS_age_degree_fecund$model <- c('TS-degree-fecundity')
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy_f_TS_age_degree_fecund <- 
        transform(tidy_f_TS_age_degree_fecund, 
                  `Effect type` = factor(`Effect type`,
                                     levels = c('total','direct','indirect', 
                                                'Prop. Mediated')))
      
      
  ### 5.5 Tidy male effects decomposition: male age + strength
      # with nests 2-3 total paternity.
    ## a) Assign to generic model
      model <- med_out_m_age_strength_tot_pat
      
    ## b) create a vector of estimates
      `Estimate` <- c(model$d.avg, model$z.avg, model$tau.coef, model$n.avg) 
      
    ## c) create a vector of lower 95% CI
      `95% CI Lower` <- c(model$d.avg.ci['2.5%'], model$z.avg.ci['2.5%'], 
                          model$tau.ci['2.5%'], model$n.avg.ci['2.5%']) 
      
    ## d) create a vector of upper 95% CI
      `95% CI Upper` <- c(model$d.avg.ci['97.5%'], model$z.avg.ci['97.5%'], 
                          model$tau.ci['97.5%'], model$n.avg.ci['97.5%'])
      
    ## e) create a vector of p-values
      `p-value` <- c(model$d.avg.p, model$z.avg.p, model$tau.p, 
                     model$n.avg.p)
      
    ## f) Combine vectors into a tibble
      tidy_m_age_strength_tot_pat <- tibble(`Effect type`, `Estimate`,
                                           `95% CI Lower`, `95% CI Upper`,
                                           `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy_m_age_strength_tot_pat$sex <- c('male')
      
    ## h) Label the estimates in data frame by sex
      tidy_m_age_strength_tot_pat$model <- c('age-strength-tot. pat.')
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy_m_age_strength_tot_pat <- 
        transform(tidy_m_age_strength_tot_pat, 
                  `Effect type` = factor(`Effect type`,
                                     levels = c('total','direct','indirect', 
                                                'Prop. Mediated')))
      
      
  ### 5.6 Tidy male effects decomposition: male tail streamer length +
      # strength with nests 2-3 total paternity - age adjusted
    ## a) Assign to generic model
      model <- med_out_m_TS_age_strength_tot_pat
      
    ## b) create a vector of estimates
      `Estimate` <- c(model$d.avg, model$z.avg, model$tau.coef, model$n.avg) 
      
    ## c) create a vector of lower 95% CI
      `95% CI Lower` <- c(model$d.avg.ci['2.5%'], model$z.avg.ci['2.5%'], 
                          model$tau.ci['2.5%'], model$n.avg.ci['2.5%']) 
      
    ## d) create a vector of upper 95% CI
      `95% CI Upper` <- c(model$d.avg.ci['97.5%'], model$z.avg.ci['97.5%'], 
                          model$tau.ci['97.5%'], model$n.avg.ci['97.5%'])
      
    ## e) create a vector of p-values
      `p-value` <- c(model$d.avg.p, model$z.avg.p, model$tau.p, 
                     model$n.avg.p)
      
    ## f) Combine vectors into a tibble
      tidy_m_TS_age_strength_tot_pat <- tibble(`Effect type`, `Estimate`,
                                              `95% CI Lower`, `95% CI Upper`,
                                              `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy_m_TS_age_strength_tot_pat$sex <- c('male')
      
    ## h) Label the estimates in data frame by sex
      tidy_m_TS_age_strength_tot_pat$model <- c('TS-strength-tot. pat.')
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy_m_TS_age_strength_tot_pat <- 
        transform(tidy_m_TS_age_strength_tot_pat, 
                  `Effect type` = factor(`Effect type`,
                                    levels = c('total','direct','indirect', 
                                               'Prop. Mediated')))
      
      
  ### 5.7 Tidy male effects decomposition: male tail streamer length +
      # degree with nests 2-3 total paternity - age adjusted
    ## a) Assign to generic model
      model <- med_out_m_TS_age_degree_tot_pat
      
    ## b) create a vector of estimates
      `Estimate` <- c(model$d.avg, model$z.avg, model$tau.coef, model$n.avg) 
      
    ## c) create a vector of lower 95% CI
      `95% CI Lower` <- c(model$d.avg.ci['2.5%'], model$z.avg.ci['2.5%'], 
                          model$tau.ci['2.5%'], model$n.avg.ci['2.5%']) 
      
    ## d) create a vector of upper 95% CI
      `95% CI Upper` <- c(model$d.avg.ci['97.5%'], model$z.avg.ci['97.5%'], 
                          model$tau.ci['97.5%'], model$n.avg.ci['97.5%'])
      
    ## e) create a vector of p-values
      `p-value` <- c(model$d.avg.p, model$z.avg.p, model$tau.p, 
                     model$n.avg.p)
      
    ## f) Combine vectors into a tibble
      tidy_m_TS_age_degree_tot_pat <- tibble(`Effect type`, `Estimate`,
                                            `95% CI Lower`, `95% CI Upper`,
                                            `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy_m_TS_age_degree_tot_pat$sex <- c('male')
      
    ## h) Label the estimates in data frame by sex
      tidy_m_TS_age_degree_tot_pat$model <- c('TS-degree-tot. pat.')
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy_m_TS_age_degree_tot_pat <- 
        transform(tidy_m_TS_age_degree_tot_pat, 
                  `Effect type` = factor(`Effect type`,
                                   levels = c('total','direct','indirect', 
                                              'Prop. Mediated')))
      
      
      
###############################################################################
##############           6. Graph effects decomposition            ############
############################################################################### 
      
  ### 6.1 Female effects decomposition graphs
    ## a) Graph results of total, direct and indirect effects from model of 
      # female age, fxm strength and total fecundity
      f_age_strength_fecund_plot <- 
        ggplot(data = (tidy_f_age_strength_fecund %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        ylim(-2, 8) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('steelblue4', 'darkblue', 'steelblue1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'age - strength - fecundity') +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        theme(axis.line = element_line(colour = 'black',
                                       linewidth = 0.5, linetype = 'solid')) +
        # change axes font style, color, size, angle, margin, and legend
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=22, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=22, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10)),
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type'))) +
        ylab(expression
             (bold('Beta estimate and 95% CI'))) 

    ## b) view plot  
      print(f_age_strength_fecund_plot)
      
  ## c) Graph results of total, direct and indirect effects from model of 
      # tail streamer length, fxm strength and total fecundity - age adjusted
      f_TS_age_strength_fecund_plot <- 
        ggplot(data = (tidy_f_age_strength_fecund %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        ylim(-2, 8) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('steelblue4', 'darkblue', 'steelblue1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'TS - strength - fecundity') +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        theme(axis.line = element_line(colour = 'black',
                                       linewidth = 0.5, linetype = 'solid')) +
        # change axes font style, color, size, angle, margin, and legend
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=22, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=22, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10)),
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type'))) +
        ylab(expression
             (bold('Beta estimate and 95% CI'))) 
      
    ## d) view plot  
      print(f_TS_age_strength_fecund_plot)
      
    ## e) Graph results of total, direct and indirect effects from model of 
      # tail streamer length, fxm degree and total fecundity - age adjusted
      f_TS_age_degree_fecund_plot <- 
        ggplot(data = (tidy_f_TS_age_degree_fecund %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        ylim(-2, 8) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('steelblue4', 'darkblue', 'steelblue1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'TS - degree - fecundity') +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        theme(axis.line = element_line(colour = 'black',
                                       linewidth = 0.5, linetype = 'solid')) +
        # change axes font style, color, size, angle, margin, and legend
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=22, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=22, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10)),
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type'))) +
        ylab(expression
             (bold('Beta estimate and 95% CI'))) 
      
    ## f) view plot  
      print(f_TS_age_degree_fecund_plot)
    
    ## g) make panel plot using patchwork
      female_panel_plot <- (f_age_strength_fecund_plot | 
                            f_TS_age_strength_fecund_plot | 
                            f_TS_age_degree_fecund_plot)
      
    ## h) view panel plot
      female_panel_plot
      
    ## i) save panel plot
      # use ggsave to save the plot
      ggsave('female_panel_plot.pdf', 
             plot = female_panel_plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 16,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE)
      
      
  ### 6.2 Male effects decomposition graphss
    ## a) Graph results of total, direct and indirect effects from model of 
      # male age, fxm strength and total fecundity
      m_age_strength_tot_pat_plot <- 
        ggplot(data = (tidy_m_age_strength_tot_pat %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        ylim(-2, 8) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('steelblue4', 'darkblue', 'steelblue1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'age - strength - tot. pat.') +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        theme(axis.line = element_line(colour = 'black',
                                       linewidth = 0.5, linetype = 'solid')) +
        # change axes font style, color, size, angle, margin, and legend
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=22, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=22, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10)),
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type'))) +
        ylab(expression
             (bold('Beta estimate and 95% CI'))) 
      
      ## b) view plot  
      print(m_age_strength_tot_pat_plot)
      
    ## c) Graph results of total, direct and indirect effects from model of 
      # male tail streamer length, fxm strength and nests 2-3 total paternity 
      # - age adjusted
      m_TS_age_strength_tot_pat_plot <- 
        ggplot(data = (tidy_m_age_strength_tot_pat %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        ylim(-2, 8) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('steelblue4', 'darkblue', 'steelblue1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'TS - strength - tot. pat.') +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        theme(axis.line = element_line(colour = 'black',
                                       linewidth = 0.5, linetype = 'solid')) +
        # change axes font style, color, size, angle, margin, and legend
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=22, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=22, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10)),
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type'))) +
        ylab(expression
             (bold('Beta estimate and 95% CI'))) 
      
      ## d) view plot  
      print(m_TS_age_strength_tot_pat_plot)
      
    ## e) Graph results of total, direct and indirect effects from model of 
      # male tail streamer length, fxm degree and nests 2-3 total paternity 
      # - age adjusted
      m_TS_age_degree_tot_pat_plot <- 
        ggplot(data = (tidy_m_TS_age_degree_tot_pat %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        ylim(-2, 8) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('steelblue4', 'darkblue', 'steelblue1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'TS - degree - tot. pat.') +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        theme(axis.line = element_line(colour = 'black',
                                       linewidth = 0.5, linetype = 'solid')) +
        # change axes font style, color, size, angle, margin, and legend
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=22, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=22, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10)),
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type'))) +
        ylab(expression
             (bold('Beta estimate and 95% CI'))) 
      
    ## f) view plot  
      print(m_TS_age_degree_tot_pat_plot)
      
    ## g) make panel plot using patchwork
      male_panel_plot <- (m_age_strength_tot_pat_plot | 
                              m_TS_age_strength_tot_pat_plot | 
                              m_TS_age_degree_tot_pat_plot)
      
    ## h) view panel plot
      male_panel_plot
      
    ## i) save panel plot
      # use ggsave to save the plot
      ggsave('male_panel_plot.pdf', 
             plot = male_panel_plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 16,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE)
      