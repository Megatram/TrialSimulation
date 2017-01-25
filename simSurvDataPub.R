simSurvDataPub <- function(N, # For explanation of input variables see below. 
                        HR,
                        cox,
                        haz.rate=0.11,
                        shape.par=.6,
                        weight=0.5,
                        risk.fact=array(c(2.27,1.87,0.62,1.66,0.26,0.44,0.28,0.50),dim=c(4,2)),
                        max.time=2,
                        return.value='p',
                        use.dose.distr=FALSE,
                        two.dose.var=FALSE,
                        dose.var.cor=0.8, 
                        use.dose.clinical=FALSE,
                        adjust.weib= FALSE,
                        dose.distr.dmld=4.473,
                        dose.distr.sd=1.895,
                        dose.distr.phmld,
                        dose.distr.phvx,
                        dose.distr.dold=12,
                        dose.distr.osd=2,
                        dose.distr.doldyD=0.126, # note OR calc for VXX according to 4*y50/d50 for IMPUT WEIBULL
                        dose.distr.dogamma50=0,
                        dose.distr.dodose50=0,
                        dose.distr.dgamma50=1.19,
                        dose.distr.ddose50=34.4,
                        dose.distr.dmldyD=0.138, # note OR calc for MLD according to 4*y50/d50 for IMPUT WEIBULL
                        ) {
  
  # Function for simulating a set of survival data. Input
  # N = number of patients. MUST BE SPECIFIED
  # HR =  hazard ratio between control and experimental arm. MUST BE SPECIFIED
  # cox = choice of cox model for regression; MUST BE SPECIFIED
  #       1 is simple model with only treatment arm, 
  #       2 is model with all risk factors, 
  #       3 is model with only "summary risk factor", 
  #       4 is model with change in MLD as covariate (and not treatment arm)
  #       5 is a simple log-rank test for difference between treatment arms
  #       6 is model with change in NTCP, photon ntcp and individual risk factors all in the CoxPHM. 
  # haz.rate = baseline hazard (i.e. control arm); assumed to be constant.
  #            If shape.par!=1, then 'haz.rate' is the scale parameter in the Weibull distribution
  #            and the hazard rate is no longer constant over time. Note that the definition used here
  #            is different than for the scale paramter 'b' used in the 'rweibull' function in R; 
  #            b=(haz.rate)^(-1/shape.par))
  # shape.par = the shape parameter in the Weibull distribution; 
  #             for shape.par=1, the distribution becomes exponential 
  # weight = percentage of experimental subjects out of all subjects
  # risk.fact = Set of clinical risk factors (hazard ratios and proportion with risk factor),
  #             arranged such that risk.fact[i,1] is a hazard ratio and 
  #             risk.fact[i,2] is the proportion with the corresponding risk factor 
  #             OR for riskfactors taken from Vogelius and Bentzen 2012. 
  #             prevalance of risk factors taken from Appelt et al 2014(?). 
  # max.time = max follow-up time in years
  # return.value = select output from function; if 'p', then output is p-value for effect of trial arm; 
  #                if 'sd', then output is standard deviation from fit (should be used together with 'cox=4'?);
  #                if 'HR', then output is the hazard ratio for the effect estimate from the Cox regression fit
  #                if 'b', then output is the beta for the effect estimate from the Cox regression fit
  #                if 'full' then output is full data set (i.e. S.data)
  # use.dose.distr = specifies whether a non-constant hazard ratio should be used for the experimental arm
  #                  (corresponds to assuming varying change in effectiveness / in dose to OAR (e.g. MLD) with experimental treatment) 
  # use.dose.clinical = specifies if use of clinical measured differences in mld should be used. 
  #                (must be used together with 'use.dose.distr=TRUE')
  # adjust.weib.rate = set to true if the individual patient risk should be modeled in the weibull distribution 
  # dose.distr.sd = specifies the standard deviation to use in case of non-constant hazard ratio.
  #                 Note that since the distribution used is log-normal, the specified s.d. should be for the normal part of the distribution
  # dose.distr.dmld = Specifies the measured difference in MLD between the two treatment arms, fex: the difference between proton and photon MLD. 
  # two.dose.var = if true, generates a second randomly distributed parameter, correlating to the value 
  #                of use.dose.distr for each patient, and fits to this parameter instead
  #                (must be used together with 'use.dose.distr=TRUE')
  # dose.var.cor = correlation (value from 0 to 1) between the second random parameter and the first (see above)
  #                (must be used together with 'two.dose.var=TRUE' and 'use.dose.clinical=FALSE')
  # dose.distr.phmld = distribution of photon mean lung dose from clinical plans
  #                (must be used together with 'use.dose.clinical=TRUE')
  # dose.distr.phvx = distribution of photon Vxx volume for the lung from clinical plans
  #                (must be used together with 'use.dose.clinical=TRUE')  
  # dose.distr.dold = Specifies the measured difference in Vxx between the two treatment arms, fex: the difference between proton and photon Vxx 
  #                (must be used together with 'two.dose.var=TRUE')
  # dose.distr.osd = specifies the standard deviation to use in case of non-constant hazard ratio.
  #                 Note that since the distribution used is log-normal, the specified s.d. should be for the normal part of the distribution
  # dose.distr.doldyD = note OR (or I_risk) calc for VXX according to 4*y50/d50 for INPUT WEIBULL
  # dose.distr.dogamma50= gamma_50 for other dose parameter ie Vxx.
  #                (must be used together with 'two.dose.var=TRUE' and 'use.dose.clinical=TRUE')
  # dose.distr.dodose50 = D_50 for other dose parameter ie Vxx.
  #                (must be used together with 'two.dose.var=TRUE' and 'use.dose.clinical=TRUE')
  # dose.distr.dgamma50=gamma_50 for other dose parameter ie Vxx.
  #                (must be used together with 'two.dose.var=TRUE' and 'use.dose.clinical=TRUE')
  # dose.distr.ddose50=D_50 for other dose parameter ie Vxx.
  #                (must be used together with 'two.dose.var=TRUE' and 'use.dose.clinical=TRUE')
  # dose.distr.dmldyD=0.138, # note OR calc for MLD according to 4*y50/d50 for INPUT WEIBULL
  
  
  # Generate error messages if input is inconsistent
  if(cox==4&&!use.dose.distr) {
    stop('In order to use \'cox=4\' (i.e. fit to non-constant effect variable), \'use.dose.distr\' must be set to TRUE')
  } else if(two.dose.var&&!use.dose.distr) {
    stop('If \'two.dose.var\' is set to TRUE, then \'use.dose.distr\' should be set to TRUE as well')
  } else if(cox==5&&!return.value=='p') {
    stop('If trials arms are to be compared using log-rank test, only p-value for comparison can be returned')
  } else if (!use.dose.distr&&use.dose.clinical) {
    stop('In order to use clinical dose parameters the \'use.dose.distr\' should be set to TRUE')
  } else if ((length(dose.distr.dmld)!=length(dose.distr.dold))&&(length(dose.distr.dmld)>1)&&two.dose.var) {
    stop('In order to use two.dose.var and bootstrap sampling of two parameters they need to be of equal length')
  } else if ((adjust.weib&&length(dose.distr.phmld)==0) || (length(dose.distr.phmld)!=length(dose.distr.dmld))) {
    stop('In order to use adjust weibull dist with photon mean lung dose(dose.distr.phmld) needs to be adjusted and the vector need the same length as dmld')
  } else if (length(dose.distr.phmld)==1 && use.dose.clinical) {
    stop('In order to use clinical DVH parameters and Bootstrap sampling of these, the variable \'use.dose.distr\ should be set to true')
  }
  #else if(return.value==c('sd')&&cox!=4) {
    #stop('In order to output standard deviatation (sd), \'cox\' must be set to 4')
  #}
  
  num.risk.fact <- length(risk.fact[,1])
  S.data <- NULL
  # Stores input parameters in the S.data
  S.data$inputparam.
  # Generate S.data which is the structure that contains all the information for the trial generation and
  # output
  S.data$data <- array(0,dim=c(N,5+num.risk.fact+use.dose.distr*1+two.dose.var*1))
  # Randomly distribute treatment arms with the specified weight
  S.data$data[,3] <- rbinom(N,1,weight)
  # Pick censoring time randomly from available interval
  c.time <- runif(N,min=0,max=max.time)
  
  # Create linear combination of ORs for event time distribution:
  # Randomly pick risk factors, based on probabilities and prevalences given
  S.data$data[,4:(3+num.risk.fact)] <- vapply(risk.fact[,2],rbinom,integer(N),n=N,size=1)
  # Create "summarized risk factor" used in weibull 
  # Here the ORs from input is converted to logOR
  if(num.risk.fact==1) {
    S.data$data[,4+num.risk.fact] <- S.data$data[,4]*log(risk.fact[1,1])
  } else {
    S.data$data[,4+num.risk.fact] <- S.data$data[,4:(3+num.risk.fact)]%*%log(risk.fact[,1])
  }
  # Create distribution of MLD sparing for experimental arm based on specified HR
  # Different ways of creating distrubiton based on input binary variables: 
  # use.dose.clinical =T and dMLD inp is single value -> generates a distribution based on the dMLD and SD of dMLD
  # use.dose.clinical =T and dMLD inp is multi values -> usees bootstrap sampling based on real patient plans
  if(use.dose.distr) {
    if (use.dose.clinical) {
      
        dMLD = dose.distr.dmld
        dMLDsd = dose.distr.sd
        if (two.dose.var)
        {
          dOLD = dose.distr.dold;
          dOLDsd = dose.distr.osd;
        }
      } else { # This part should never be used. 
        dMLD = -log(HR)/0.138 # Calculate change in MLD corresponding to HR # 0.138
        dMLDsd = 2 # Note: Using 2 Gy as Stddev (From Anes previous)
        rand.gen <- mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,dose.var.cor,dose.var.cor,1),2,2))
      } 
    # Randomly pick change in MLD with experimental treatment
    
    if (length(dose.distr.dmld) > 1){ # If dMLD is longer than 1 ie not the mean_dMLD but vector with dMLD force the script to sample only from insertet vector
      #The part below generates the f which containts the following: 
      #dMLD - difference in mean lung dose, 
      #dOLD - difference in other dose metric in case the two.dose.var is used. 
      #dose.distr.phmld - mean lung dose from photon therapy plan
      #dose.distr.phvx - Vxx parameter from the photon therapy plan in case of the two.dose.var used. 
      if (two.dose.var)
      {
        # Selects random samples from both distributions but selects the same row as variables. 
        
        f <-matrix(c(dMLD,dOLD,dose.distr.phmld,dose.distr.phvx),nrow = length(dMLD),ncol = 4)
        f <-f[sample(nrow(f), N,replace = TRUE), ]
        S.data$dataf <- f
        
      }
      else
      {
        f <-matrix(c(dMLD,dMLD,dose.distr.phmld),nrow = length(dMLD),ncol = 3) # Ugly hack. Please change. 
        f <- f[sample(nrow(f), N,replace = TRUE), ]
        S.data$dataf <- f
      }
      
      S.data$data[,6+num.risk.fact] = f[,1]*S.data$data[,3] # Bootstrap sampling from the inserted vector with dMLDs XXX CHANGED TO SECONDARY VARIABLE
      # This part changes so only the dMLD exists for proton plan patients (S.data$data[,3] = 1) otherwise it is zero. 
      # same as trial arm indicator
      # It also generates a sumOR (According to Appelt et al 201X) for each patient in cases where this can be interesting to study
      y50 <- dose.distr.dgamma50
      d50 <- dose.distr.ddose50
      if (num.risk.fact==1) {
        sumOR <- S.data$data[,4]*risk.fact[,1]
        if (sumOR==0) {
          sumOR <- sumOR+1
        } else {
          sumOR[sumOR == 0] <- NA
          sumOR <- apply(sumOR,1,prod,na.rm=TRUE)
          S.data$sumOR <- sumOR
        }
      } else {
        gblOR <- S.data$data[,4:(3+num.risk.fact)]%*%risk.fact[,1]
        S.data$gblOR <- gblOR
        
        
        gblOR[gblOR == 0] <- NA # Create NA values so if no riskfactors the resulting OR will be 1. 
        sumOR <- apply(gblOR,1,prod,na.rm=TRUE)
        S.data$sumOR <- sumOR
      }
      # calculates NTCPs for photon and proton using a logistic model
      phNTCP <- 1/(1+exp(4*y50*(1-f[,3]/d50)))
      prNTCP <- 1/(1+exp(4*y50*(1-(f[,3]-f[,1])/d50)))
      phNTCP.wrf <- 1/(1+sumOR^(-1)*exp(4*y50*(1-f[,3]/d50)))
      prNTCP.wrf <- 1/(1+sumOR^(-1)*exp(4*y50*(1-(f[,3]-f[,1])/d50)))
      
      dNTCP <- (phNTCP-prNTCP)
      dNTCP.wrf <- (phNTCP.wrf-prNTCP.wrf)
      dNTCPOR <- (prNTCP/(1-prNTCP))/(phNTCP/(1-phNTCP))
      #retrive the corrct NTCP data for photons, protons and difference. 
      
      S.data$ntcp.phNTCP <- phNTCP
      S.data$ntcp.prNTCP <- prNTCP
      S.data$ntcp.dNTCP <- dNTCP*S.data$data[,3]
      S.data$ntcp.phNTCP.wrf <- phNTCP.wrf
      S.data$ntcp.prNTCP.wrf <- prNTCP.wrf
      S.data$ntcp.dNTCP.wrf <- dNTCP.wrf*S.data$data[,3]
      
      #Set all patients in photon arm to OR = 1 Cannot multiply here as it is a OR (ie OR=1 no diff)
      dNTCPOR[S.data$data[,3]==0] <- 1
      S.data$ntcp.dNTCPOR <-dNTCPOR
    } 
    else # in case of single input dMLD (ie mean dMLD from all patients) # dont use this!
    { 
        S.data$data[,6+num.risk.fact] <- (rand.gen[,1]*dMLDsd+dMLD)*S.data$data[,3] # generates the first list of variables
    }
    if(two.dose.var) { # If specified, generate also second variable correlating to first 
      if(use.dose.distr)
      {
        if(!use.dose.clinical)
        {
        S.data$data[,7+num.risk.fact] <- (rand.gen[,2]*dMLDsd+dMLD)*S.data$data[,3] # generates the first list of variables
      }
      if (length(dose.distr.dold) > 1)
      { # If dOtherLD is longer than 1 ie not the mean_dOtherLD but vector with dOtherLD force the script to sample only from insertet vector
        S.data$data[,7+num.risk.fact] = f[,2]*S.data$data[,3] # Bootstrap sampling from the inserted vector with dMLDs
      } 
      else 
      {
        S.data$data[,7+num.risk.fact] <- (rand.gen[,2]*dOLDsd+dOLD)*S.data$data[,3] # generates the second list of variables
      }
      } 
    }

  }
  
  
  # Create linear combination of all factors, including either treatment arm or change in MLD, if this is calculated
  if(use.dose.distr) {
    
    if(two.dose.var) {
      # Note here is where the adjustment total OR is generated for adjustment of the weilbull-distribution. 
      S.data$data[,5+num.risk.fact] <- -S.data$data[,7+num.risk.fact]*dose.distr.doldyD + S.data$data[,4+num.risk.fact] 
    }
    else
    {
      S.data$data[,5+num.risk.fact] <- -S.data$data[,6+num.risk.fact]*dose.distr.dmldyD + S.data$data[,4+num.risk.fact] 
    }
  } else {
    S.data$data[,5+num.risk.fact] <- S.data$data[,3]*log(HR) + S.data$data[,4+num.risk.fact]
  }
  
  # Create list of event times and events:
  # Pick event time from Weibull distr based on baseline hazard and treatment arm
  if (adjust.weib){
    # this part uses input y50 and D50 to calculate a 15%probability for RP from the 
    if (two.dose.var){ #This case is used if the incorrect dose parameter (ie. the Vxx parameter should be used
      #in the weibull generation. )
      rpfit <- smooth.spline(1/(1+exp(4*dose.distr.dogamma50*(1-seq(0,100,length=1000)/dose.distr.dodose50))),seq(0,100,length=1000))
      prob15 <- predict(rpfit,0.15)$y
      haz.rate.corrcfact <- (f[,4]-prob15)*4*dose.distr.dogamma50/dose.distr.dodose50
      #Note: Here we convert the individyal Vxx_proton parameter f[,4] to a risk reduction into the haz.rate.corrcfact
      
    } else {
      # To adjust the baseline hazard rate (which is defined by quantec 17 Gy MLD for 15% prob for RP)
      #Hence the y50/D50 is taken from the quantec paper
      y50 <- dose.distr.dgamma50
      d50 <- dose.distr.ddose50
      rpfit <- smooth.spline(1/(1+exp(4*y50*(1-seq(0,100,length=1000)/d50))),seq(0,100,length=1000))
      prob15 <- predict(rpfit,0.15)$y #finds at what mld the probability is 15 % using a logistic model.
      #prob15 <- 17.03049
   haz.rate.corrcfact <- (f[,3]-prob15)*4*y50/d50
   #Note: Here we convert the individyal MLD_photon parameter f[,4] to a risk reduction into the haz.rate.corrcfact
   
   S.data$haz <- haz.rate.corrcfact
   
   QTNTCP15 <- 1/(1+exp(4*y50*(1-prob15/d50)))
   QTNTCPph <- 1/(1+exp(4*y50*(1-(f[,3])/d50)))
   haz.rate.corrcperc <- QTNTCPph - QTNTCP15
   S.data$hazperc <- haz.rate.corrcperc
   S.data$hazpercph <- QTNTCPph
   S.data$hazperc15 <- QTNTCP15
   
   haz.rate.corrcfactOR <- (QTNTCPph/(1-QTNTCPph))/(QTNTCP15/(1-QTNTCP15))
   S.data$hazOR <- haz.rate.corrcfactOR
    }
    #here the eventime is corrected with the patient individual haz.rate.corrcfact based
    #on either the Vxx if misspec model or the MLD for correct model. 
   event.time <- rweibull(N,shape=shape.par,scale=(haz.rate*exp(S.data$data[,5+num.risk.fact]+haz.rate.corrcfact))^(-1/shape.par))
   #S.data$datax <- S.data$data[,5+num.risk.fact]+haz.rate.corrcfact
  } else {
  event.time <- rweibull(N,shape=shape.par,scale=(haz.rate*exp(S.data$data[,5+num.risk.fact]))^(-1/shape.par))
  #S.data$datax <- S.data$data[,5+num.risk.fact]
  }
  S.data$event.time <- event.time
  # Let time variable be minimum of event and censoring time;
  S.data$data[,1] <- pmin(c.time,event.time)
  # Event indicator based on whether event or censoring time is smallest
  S.data$data[,2] <- (event.time<c.time)*1 
  
  # Estimate survival function from simulated data
  S.data$Surv <- Surv(S.data$data[,1],S.data$data[,2])
  S.data$Surv.fit <- survfit(S.data$Surv~1)
  
  # Fix Cox model dependent on input specification
  if(cox==1) {
    # Fit Cox regression model with only treatment arm as covariate 
    S.data$Surv.cox <- coxph(S.data$Surv ~ S.data$data[,3])
  } else if(cox==2) {
    # Fit Cox regression model with all clinical covariates as well as treatment arm
    S.data$Surv.cox <- coxph(S.data$Surv ~ S.data$data[,3] + S.data$data[,4:(3+num.risk.fact)])
  } else if(cox==3) {
    # Fit Cox regression model with "summary risk factor" as well as treatment arm
    S.data$Surv.cox <- coxph(S.data$Surv ~ S.data$data[,3] + S.data$data[,(4+num.risk.fact)])
  } else if(cox==4) {
    if (two.dose.var) {
      # Fit Cox regression model with change in correlated dose parameter as covariate (and not treatment arm)
      S.data$Surv.cox <- coxph(S.data$Surv ~ S.data$data[,7+num.risk.fact])
    } else {
      # Fit Cox regression model with change in MLD as covariate (and not treatment arm)
      S.data$Surv.cox <- coxph(S.data$Surv ~ S.data$data[,6+num.risk.fact])
    }
  } else if(cox==5) {
    #Logrank part on TX-arm. 
      S.data$Surv.cox <- survdiff(S.data$Surv ~ S.data$data[,3])
  } else if(cox == 6){
    #Retrive cox when using dNTCP, dNTCP-photons and all riskfactors as covariates in the CoxPHM
    S.data$Surv.cox <- coxph(S.data$Surv ~ S.data$ntcp.dNTCP+S.data$ntcp.phNTCP+S.data$data[,4:(3+num.risk.fact)])
  }
  
  # Extract statistics from cox fit
  if (cox==5) {
    p.value <- 1 - pchisq(S.data$Surv.cox$chisq, length(S.data$Surv.cox$n) - 1)
  }  else {
    S.data$Surv.cox.summary <- summary(S.data$Surv.cox)
    p.value <- S.data$Surv.cox.summary$coefficients[1,5]
    sd.value <- S.data$Surv.cox.summary$coefficients[1,3]
    HR.value <- S.data$Surv.cox.summary$coefficients[1,2]
    b.value <- S.data$Surv.cox.summary$coefficients[1,1]
  }
    
  # Plot survival fit
  # plot(S.data$Surv.fit)
  
  # p.values <- c(S.data$Surv.cox1.summary$coefficients[1,5],S.data$Surv.cox2.summary$coefficients[1,5],S.data$Surv.cox3.summary$coefficients[1,5])

  if(return.value=='full') {
    return(S.data)
  } else if(return.value=='sd') {
    return(sd.value)
  } else if(return.value=='HR') { 
    return(HR.value)
  } else if (return.value=='weib') {
    return(sum(event.time<2)/N)
  } else if (return.value =='b' || return.value=='logHR') {
    return(b.value)
  } else {
    return(p.value)
  } 
  
}
