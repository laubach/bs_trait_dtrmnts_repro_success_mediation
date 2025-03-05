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
    # 3. Female models effects decomposition
    # 4. Male models effects decomposition
    # 5. Tidy effects decomposition
    # 6. Graph effects decomposition

    

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
      load(here('data/7_chr15_effect_decomp_data.RData'))
    
    ## b) Load RData mediation models
      load(here('data/7_chr15_effect_decomp_mods.RData'))
    
      

###############################################################################
##############      3. Female models effects decomposition       ##############
###############################################################################  

      
  ### 3.1 Female mediation model: female age plus pre-manip strength
      # with total fecundity.
      # exposure = age
      # mediators = pre-manip strength
      # outcomes = total fecundity (females)
      
      # Estimate the average causal mediation effects (ACME) and the 
      # average direct effect (ADE); 'mediation' package
      # mediation model = fem.age.strength.pre;see section 4.1
      # outcome model = fem.age.pre.strength.tot.fecund
      
    ## a) fit full model
      fem.age.pre.strength.tot.fecund <- glm(total.fecundity ~ 
                                          Age.category + strength.pre,
                                          # test for intx
                                          #Age.category * strength.pre,
                                          family = 'gaussian',
                                          data = subset(chr15_attrib_pre_df,
                                                Sex == 'f' &
                                                !is.na(strength.pre) &
                                                !is.na(total.fecundity)))
    ## b) full model summary
      summary(fem.age.pre.strength.tot.fecund)
      
    ## c) Estimate average causal mediation effects
      med.out.fem.age.pre.strength.fecund <- mediate(fem.age.strength.pre, 
                                                fem.age.pre.strength.tot.fecund, 
                                                treat = "Age.category", 
                                                mediator = "strength.pre",
                                                boot = TRUE, sims = 1000) 
    ## d) test for exposure*mediator intx
      # test.TMint(med.out.fem.age.pre.strength.fecund, 
      #            conf.level = .95)
      
    ## e) Average causal mediation model summary
      summary(med.out.fem.age.pre.strength.fecund)
      
    # # ## f) sensitivity analysis for unmeasured confounder between mediator and 
    #   # outcome
    #   sens.out <- medsens(med.out.fem.age.pre.strength.fecund,
    #                       rho.by = 0.1, effect.type = 'indirect',
    #                       sims = 1000)
    #   summary(sens.out)
    
      
##*** MODEL FITTING TEST using 'medflex' package
    # ## g) Expanded data and weights for female pre-manip strength by age
    #   # 'expand the data' by estimating mediator for different levels of the 
    #   # exposure to generate weights; implemented in medflex 
    #   exp.wght.fem.age.strength.pre <- neWeight(fem.age.strength.pre)
    #   weights(exp.wght.fem.age.strength.pre)
    #   
    # # h) Ne_impute method; implemented in medflex using full model
    #   # fem.age.pre.strength.tot.out 
    #   exp.imp.fem.age.strength.pre <- neImpute(fem.age.pre.strength.tot.out)
    # 
    # 
    # ## i) Natural effects model female pre-manip strength by age 
    #   # calculate bootstrapped SE based resampling original data 
    #   # with replacement (1000 iterations); implemented in medflex 
    #   
    #   # Weighting method   
    #   ne.wght.fem.age.strength.pre.fecund <- neModel(total.fecundity ~ 
    #                                   Age.category0 + Age.category1,
    #                                   family = gaussian, 
    #                                   expData = exp.wght.fem.age.strength.pre)
    #   
    #   summary(ne.wght.fem.age.strength.pre.fecund)
    #   
    #   # Imputation method  
    #   ne.imp.fem.age.strength.pre.fecund <- neModel(total.fecundity ~ 
    #                                   Age.category0 + Age.category1,
    #                                   family = gaussian, 
    #                                   expData = exp.imp.fem.age.strength.pre)
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

      
  ### 3.2 Female mediation model: female age plus post manip degree
        # with total fecundity.
        # exposure = age
        # mediators = post manip degree
        # outcomes = total fecundity (females)  
      
        # Estimate the average causal mediation effects (ACME) and the 
        # average direct effect (ADE); 'mediation' package
        # mediation model = fem.age.degree.post; see section 4.1
        # outcome model = fem.age.post.deg.tot.fecund
        
    ## a) fit full model
      fem.age.post.deg.tot.fecund <- glm(total.fecundity ~ 
                                        Age.category + degree.post,
                                        # test for intx
                                        # Age.category * degree.post,
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                                  Sex == 'f' &
                                                  !is.na(degree.post) &
                                                  !is.na(total.fecundity)))
    ## b) full model summary
      summary(fem.age.post.deg.tot.fecund)
        
    ## d) Estimate average causal mediation effects
      med.out.fem.age.post.deg.fecund <- mediate(fem.age.degree.post, 
                                                   fem.age.post.deg.tot.fecund, 
                                                   treat = "Age.category", 
                                                   mediator = "degree.post",
                                                   boot = TRUE, sims = 1000) 
        
    ## d) test for exposure*mediator intx
        # test.TMint(med.out.fem.age.post.deg.fecund, 
        #            conf.level = .95)
        
    ## e) Average causal mediation model summary
      summary(med.out.fem.age.post.deg.fecund)
        
    # ## f) sensitivity analysis for unmeasured confounder between mediator and 
        # outcome
        # sens.out <- medsens(med.out.fem.age.post.deg.fecund, 
        #                     rho.by = 0.1, effect.type = "indirect", 
        #                     sims = 1000)
        # summary(sens.out)
        

          
###############################################################################
##############        4. Male models effects decomposition         ############
###############################################################################
      

  ### 4.1 Male mediation model: male age plus pre-manip strength
        # with first clutch total paternity.
        # exposure = age
        # mediators = pre-manip strength
        # outcomes = clutch 1 total paternity 
        
       # Estimate the average causal mediation effects (ACME) and the 
        # average direct effect (ADE); 'mediation' package
        # mediation model = m.age.strength.pre; see section 4.2
        # outcome model = m.age.pre.strength.attmpt.1.tot.pat.out
        
    ## a) fit full model
      m.age.pre.strength.attmpt.1.tot.pat.out <- glm(attmpt.1.tot.pat ~ 
                                        #Age.category + strength.pre,
                                        # test for intx
                                        Age.category * strength.pre,
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_pre_df,
                                                  Sex == 'm' &
                                                  !is.na(strength.pre) &
                                                  !is.na(attmpt.1.tot.pat)))

    ## b) full model summary
      summary(m.age.pre.strength.attmpt.1.tot.pat.out)
        
    ## c) Estimate average causal mediation effects
      med.out.m.age.pre.strength.attmpt.1.tot.pat <- mediate(m.age.strength.pre, 
                                    m.age.pre.strength.attmpt.1.tot.pat.out, 
                                    treat = "Age.category", 
                                    mediator = "strength.pre",
                                    boot = TRUE, sims = 1000) 
        
    # ## d) test for exposure*mediator intx
    #     test.TMint(med.out.m.age.pre.strength.attmpt.1.tot.pat, 
    #                conf.level = .95)
        
    ## e) Average causal mediation model summary
        summary(med.out.m.age.pre.strength.attmpt.1.tot.pat)
        
    # ## f) sensitivity analysis for unmeasure confounder between mediator and 
        # outcome
        # sens.out <- medsens(med.out.m.age.pre.strength.attmpt.1.tot.pat, 
        #                     rho.by = 0.1, effect.type = "indirect", 
        #                     sims = 1000)
        # summary(sens.out)        
    
 
  ### 4.2 Male mediation model: male age plus post manip degree
        # with second clutch total paternity.
        # exposure = age
        # mediators = post manip degree
        # outcomes = clutch 2 total paternity 
        
        # Estimate the average causal mediation effects (ACME) and the 
        # average direct effect (ADE); 'mediation' package
        # mediation model = m.age.degree.post; see section 4.2
        # outcome model = m.age.post.degree.attmpt.2.tot.pat.out
        
    ## a) fit full model
       m.age.post.degree.attmpt.2.tot.pat.out <- glm(attmpt.2.tot.pat ~ 
                                        Age.category + degree.post,
                                        # test for intx
                                        #Age.category * degree.post,
                                        family = 'gaussian',
                                        data = subset(chr15_attrib_post_df,
                                                  Sex == 'm' &
                                                  !is.na(degree.post) &
                                                  !is.na(attmpt.2.tot.pat)))
        
    ## b) full model summary
      summary(m.age.post.degree.attmpt.2.tot.pat.out)
        
    ## c) Estimate average causal mediation effects
      med.out.m.age.post.degree.attmpt.2.tot.pat <- mediate(m.age.degree.post, 
                                      m.age.post.degree.attmpt.2.tot.pat.out, 
                                      treat = 'Age.category', 
                                      mediator = 'degree.post',
                                      boot = TRUE, sims = 1000) 
        
    # ## d) test for exposure*mediator intx
    #   test.TMint(med.out.m.age.post.degree.attmpt.2.tot.pat, 
    #                conf.level = .95)
        
    ## e) Average causal mediation model summary
      summary(med.out.m.age.post.degree.attmpt.2.tot.pat)
        
    # ## f) sensitivity analysis for unmeasured confounder between mediator and 
        # outcome
        # sens.out <- medsens(med.out.m.age.post.degree.attmpt.2.tot.pat, 
        #                     rho.by = 0.1, effect.type = "indirect", 
        #                     sims = 1000)
        # summary(sens.out)    
        
    
  ### 4.3 Male mediation model: male age plus pre-manip strength
        # with all clutches total paternity.
        # exposure = age
        # mediators = pre-manip strength
        # outcomes = all clutches total paternity 
        
        # Estimate the average causal mediation effects (ACME) and the 
        # average direct effect (ADE); 'mediation' package
        # mediation model = m.age.strength.pre; see section 4.2
        # outcome model = m.age.pre.strength.total.paternity.out
      
    ## a) fit full model
      m.age.pre.strength.total.paternity.out <- glm(total.paternity ~ 
                                              Age.category + strength.pre,
                                              # test for intx
                                              #Age.category * strength.pre,
                                              family = 'gaussian',
                                              data = subset(chr15_attrib_pre_df,
                                                     Sex == 'm' &
                                                     !is.na(strength.pre) &
                                                     !is.na(total.paternity)))
      
    ## b) full model summary
      summary(m.age.pre.strength.total.paternity.out)
      
    ## c) Estimate average causal mediation effects
      med.out.m.age.pre.strength.total.paternity <- mediate(m.age.strength.pre, 
                                        m.age.pre.strength.total.paternity.out, 
                                        treat = "Age.category", 
                                        mediator = "strength.pre",
                                        boot = TRUE, sims = 1000) 
      
    # ## d) test for exposure*mediator intx
    #   test.TMint(med.out.m.age.pre.strength.total.paternity, 
    #              conf.level = .95)
      
    ## e) Average causal mediation model summary
      summary(med.out.m.age.pre.strength.total.paternity)
      
    # ## f) sensitivity analysis for unmeasure confounder between mediator and 
      # outcome
      # sens.out <- medsens(med.out.m.age.pre.strength.total.paternity, 
      #                     rho.by = 0.1, effect.type = "indirect", 
      #                     sims = 1000)
      # summary(sens.out)      
        
 
     
   
      
      
###############################################################################
##############            5. Tidy effects decomposition            ############
###############################################################################

  ### 5.1 Set standard values for the tidy tibbles
    ## a) a list of effect types
      `Effect type` <- c('indirect', 'direct', 'total', 'Prop. Mediated')
     

  ### 5.2 Tidy female effects decomposition: female age plus pre-manip strength
      # with total fecundity.
    ## a) Assign to generic model
      model <- med.out.fem.age.pre.strength.fecund
    
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
      tidy.fem.age.pre.strength.fecund <- tibble(`Effect type`, `Estimate`,
                                                `95% CI Lower`, `95% CI Upper`,
                                                `p-value`)
  
    ## g) Label the estimates in data frame by sex
      tidy.fem.age.pre.strength.fecund$sex <- c('female')
      
    ## h) Label the estimates in data frame by sex
      tidy.fem.age.pre.strength.fecund$model <- c('age-strength-fecundity')
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy.fem.age.pre.strength.fecund <- 
        transform(tidy.fem.age.pre.strength.fecund, 
                  `Effect type` = factor(`Effect type`,
                                  levels = c('total','direct','indirect', 
                                                    'Prop. Mediated'
                                         )))
      
      
  ### 5.3 Tidy female effects decomposition: female age plus post manip degree
      # with total fecundity.
    ## a) Assign to generic model
      model <- med.out.fem.age.post.deg.fecund
      
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
      tidy.fem.age.post.deg.fecund <- tibble(`Effect type`, `Estimate`,
                                                 `95% CI Lower`, `95% CI Upper`,
                                                 `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy.fem.age.post.deg.fecund$sex <- c('female')
      
    ## h) Label the estimates in data frame by sex
      tidy.fem.age.post.deg.fecund$model <- c('age-degree-fecundity')  
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy.fem.age.post.deg.fecund <- 
        transform(tidy.fem.age.post.deg.fecund, 
                  `Effect type` = factor(`Effect type`,
                                levels = c('total','direct','indirect', 
                                                    'Prop. Mediated'
                                         )))
      
      
  ### 5.4 Tidy male effects decomposition: male age plus pre-manip strength
      # with first clutch total paternity.
    ## a) Assign to generic model
      model <- med.out.m.age.pre.strength.attmpt.1.tot.pat
      
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
      tidy.m.age.pre.strength.attmpt.1.tot.pat <- tibble(`Effect type`, 
                                            `Estimate`, `95% CI Lower`, 
                                            `95% CI Upper`, `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy.m.age.pre.strength.attmpt.1.tot.pat$sex <- c('male')
      
    ## h) Label the estimates in data frame by sex
      tidy.m.age.pre.strength.attmpt.1.tot.pat$model <- 
        c('age-strength-clutch1-paternity') 
     
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy.m.age.pre.strength.attmpt.1.tot.pat <- 
        transform(tidy.m.age.pre.strength.attmpt.1.tot.pat, 
                  `Effect type` = factor(`Effect type`,
                                  levels = c('total','direct','indirect', 
                                                    'Prop. Mediated'
                                         )))
      
    
  ### 5.5 Tidy male effects decomposition: male age plus post manip degree
      # with second clutch total paternity.
    ## a) Assign to generic model
      model <- med.out.m.age.post.degree.attmpt.2.tot.pat
      
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
      tidy.m.age.post.deg.attmpt.2.tot.pat <- tibble(`Effect type`, 
                                                         `Estimate`, `95% CI Lower`, 
                                                         `95% CI Upper`, `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy.m.age.post.deg.attmpt.2.tot.pat$sex <- c('male')
      
    ## h) Label the estimates in data frame by sex
      tidy.m.age.post.deg.attmpt.2.tot.pat$model <- 
        c('age-deg-clutch2-paternity')  
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy.m.age.post.deg.attmpt.2.tot.pat <- 
        transform(tidy.m.age.post.deg.attmpt.2.tot.pat, 
                  `Effect type` = factor(`Effect type`,
                                  levels = c('total','direct','indirect', 
                                                    'Prop. Mediated'
                                         )))
      
      
  ### 5.6 Tidy male effects decomposition: male age plus pre-manip strength
      # with all clutches total paternity.
    ## a) Assign to generic model
      model <- med.out.m.age.pre.strength.total.paternity
      
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
      tidy.m.age.pre.strength.total.paternity <- tibble(`Effect type`, 
                                                         `Estimate`, `95% CI Lower`, 
                                                         `95% CI Upper`, `p-value`)
      
    ## g) Label the estimates in data frame by sex
      tidy.m.age.pre.strength.total.paternity$sex <- c('male')
      
    ## h) Label the estimates in data frame by sex
      tidy.m.age.pre.strength.total.paternity$model <- 
        c('age-strength-total-paternity')  
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'Effect type' variable 
      tidy.m.age.pre.strength.total.paternity <- 
        transform(tidy.m.age.pre.strength.total.paternity, 
                  `Effect type` = factor(`Effect type`,
                                levels = c('total','direct','indirect', 
                                                    'Prop. Mediated'
                                         )))
      
      
      
      
###############################################################################
##############           6. Graph effects decomposition            ############
############################################################################### 
      
  ### 6.1 Female effects decomposition graphs
    ## a) Graph results of total, direct and indirect effects from model of 
      # female age, pre-manip strength and total fecundity
      fem.age.pre.strength.fecund.plot <- 
        ggplot(data = (tidy.fem.age.pre.strength.fecund %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
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
        labs(title = 'Female age, pre-manipulation strength, 
and total fecundity model') +
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
      print(fem.age.pre.strength.fecund.plot)
      
    ## c) Graph results of total, direct and indirect effects from model of 
      # female age, post manip degree and total fecundity
      fem.age.post.deg.fecund.plot <- 
        ggplot(data = (tidy.fem.age.post.deg.fecund %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
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
        labs(title = 'Female age, post manipulation degree, 
and total fecundity model') +
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
              axis.title.y = element_blank(), # remove y axis title
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type')))
        # +
        # ylab(expression
        #      (bold('Beta estimate and 95% CI'))) 
      
    ## d) view plot
      print(fem.age.post.deg.fecund.plot)
      
    
    ## e) make panel plot using patchwork
      female.panel.plot <- (fem.age.pre.strength.fecund.plot | 
                       fem.age.post.deg.fecund.plot)
      
    ## f) view panel plot
      female.panel.plot
      
    ## h) save panel plot
      # use ggsave to save the plot
      ggsave('female.panel.plot.pdf', 
             plot = female.panel.plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 12,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE)
      
      
  ### 6.2 Male effects decomposition graphs
    ## a) Graph results of total, direct and indirect effects from model of 
      # male age, pre-manip strength and first clutch paternity
      m.age.pre.strength.attmpt.1.tot.pat.plot <- 
        ggplot(data = (tidy.m.age.pre.strength.attmpt.1.tot.pat %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('seagreen3', 'seagreen4', 'seagreen1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'Male age, pre-manipulation strength, 
and clutch 1 paternity model') +
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
      print(m.age.pre.strength.attmpt.1.tot.pat.plot)
      
    ## c) Graph results of total, direct and indirect effects from model of 
      # male age, post manip degree and clutch 2 paternity
      m.age.post.deg.attmpt.2.tot.pat.plot <- 
        ggplot(data = (tidy.m.age.post.deg.attmpt.2.tot.pat %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('seagreen3', 'seagreen4', 'seagreen1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'Male age, post manipulation degree, 
and clutch 2 paternity model') +
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
              #axis.title.y = element_blank(), # remove y axis title
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type')))+
        ylab(expression
             (bold('Beta estimate and 95% CI'))) 
      
    ## d) view plot
      print(m.age.post.deg.attmpt.2.tot.pat.plot)
      
    ## e) Graph results of total, direct and indirect effects from model of 
      # male age, pre-manip degree and all clutches total paternity
      m.age.strength.total.paternity.plot <- 
        ggplot(data = (tidy.m.age.pre.strength.total.paternity %>%
                         filter(`Effect type` != 'Prop. Mediated')), 
               aes(x = `Effect type`, 
                   y = Estimate, 
                   color = `Effect type`)) +
        geom_hline(yintercept = 0, color = 'red',
                   linetype = 2, size = 2) + # line at null behind coefs
        geom_point(size = 8) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), width = .1, size = 2,
                      position = position_dodge(0.5))+
        scale_color_manual(values=c('seagreen3', 'seagreen4', 'seagreen1'
                                    # , 'firebrick4'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'Male age, pre-manipulation strength, 
and all clutches paternity model') +
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
              axis.title.y = element_blank(), # remove y axis title
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Effect Type')))
      # +
      # ylab(expression
      #      (bold('Beta estimate and 95% CI'))) 
      
    ## f) view plot
      print(m.age.strength.total.paternity.plot)
      
    ## g) make panel plot using patchwork
      male.panel.plot <- (m.age.pre.strength.attmpt.1.tot.pat.plot / 
                            (m.age.post.deg.attmpt.2.tot.pat.plot | 
                            m.age.strength.total.paternity.plot))
      
    ## h) view panel plot
      male.panel.plot
      
    ## i) save panel plot
      # use ggsave to save the plot
      ggsave('male.panel.plot.pdf', 
             plot = male.panel.plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 12,
             height = 12,
             units = c('in'), dpi = 300, limitsize = TRUE)