library('shiny')
library('pec')
library('survival')
library('ggplot2')

#these are the characteristics of the "reference patient"
#Pre-existing pulmonary co-morbidity","Mid or inferior tumor location","Current smoker","Old age"
Comorb <- 'No'
TumorSite<-'Superior'
Tobacco<-'Current smoker'
Age<-60.0

# ORs
ORComorb <- 2.27
ORTloc <- 1.87
ORTob <- 0.62
ORAge <- 1.66
mldph <- 15
mldpr <- 12

# Model parameters
betadNTCP <- -7.77
betadNTCPsd <- 4.50
refPatient<-data.frame(Comorb,TumorSite,Tobacco,Age,ORComorb,ORTloc,ORTob,ORAge,mldph,mldpr,betadNTCP,betadNTCPsd)  

#ui defines what is showed in the app as input choises
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel( 
      #tags$head(
      #  tags$style(type="text/css", "input { font-size:10px; width:40px; display:inline-block; }"),
      #  tags$style(type="text/css", "#lml, #lmml, { font-size:10px; width:1000px; display:inline-block; }"),
      #  tags$style(type="text/css", "label { font-size:10px; display:inline-block; }")
      #),
      h4("Pre-existing pulmonary co-morbidity:"),
      fluidRow(
        column(6,selectInput('Comorb','-',c('Yes','No'),multiple = FALSE)), 
        column(6,numericInput('ORComorb','Odds-Radio: ',value = 2.27,min = 0))
        ),
      h4("Tumor location:"),
      fluidRow(
        column(6,selectInput('TumorSite','-',c('Superior','Mid','Inferior'),multiple = FALSE)),
        column(6,numericInput('ORTloc','Odds-Radio: ',value = 1.87,min = 0))
      ),
      h4("Smoking history:"),
      fluidRow(
        column(6,selectInput('Tobacco','-',c('Never/previous smoker','Current smoker'),multiple=FALSE,selected='Current smoker')),
        
        column(6, numericInput('ORTob','Odds-Radio: ',value = 0.62,min = 0))
        ),
      h4("Age:"),
      fluidRow(
        column(9,sliderInput('Age','',value=60.0,min=33, max=88)),
        
          column(3,numericInput('ORAge','Odds-Radio: ',value = 1.66,min = 0))
        ),
      h4("DVH data:"),
      fluidRow(
        column(6,numericInput('mldph','Mean lung dose pHotons [Gy]:',20,min = 0,max = 50)),
        column(6,numericInput('mldpr','Mean lung dose pRotons [Gy]:',15,min = 0,max = 50))
        ),
      fluidRow(
        column(6,numericInput('betadNTCP','Beta for dNTCP from CoxPHM',value=betadNTCP,min = -20,max = 20)),
        column(6,numericInput('betadNTCPsd','SD for beta from dNTCP in CoxPHM',value=betadNTCPsd,min = 0,max = 30))
      ),
      h4("Plot details:"),
      checkboxInput('showtext','Show text in plot', value = FALSE),
      fluidRow(
        column(6,checkboxInput('forcexLim','Fixed x-axis', value = FALSE)), 
        column(6,numericInput('xLim','x-lim value [%]:',35,min = 0,max = 100))
      ),
      h4("Math:"),
      withMathJax(),
      helpText('NTCP calculated according to:
           $$NTCP_{photons} = \\frac{1}{1+OR^{-1}_{sum}exp\\{4\\gamma^{0}_{50} (1-\\frac{MLD}{D^0_{50}})\\}}$$'), 
      p('from',a('Appelt et al 2014',href='http://dx.doi.org/10.3109/0284186X.2013.820341'),'.'), 
      p('NTCP for protons calculated according to:'),
      withMathJax(),
      helpText('$$1-NTCP_{protons} = (1-NTCP_{photons})^{I_{risk}}$$'), 
      p('OR for risk factors from ',a('Vogelius and Bentzen 2012',href='http://dx.doi.org/10.3109/0284186X.2012.718093'),'.') 
      
      
      
                
    ),
    mainPanel(
      titlePanel('Estimated risk for radiation pneumonitis'),
                
                #defines what kind of output you want. Should be matched with the render*() function in "server"
         
      plotOutput('test')
                
                )
  )
)

#function that is called from "server" below. It is used for manipulating (in this case making a data frame) the data chosen in the input.
fun<-function(Comorb,TumorSite,Tobacco,Age,ORComorb,ORTloc,ORTob,ORAge,mldph,mldpr,betadNTCP,betadNTCPsd,xLim,forcexLim) {data.frame(Comorb,TumorSite,Tobacco,Age,ORComorb,ORTloc,ORTob,ORAge,mldph,mldpr,betadNTCP,betadNTCPsd,xLim,forcexLim)}  

#Function for plot creation: 
getProtonProbReductionRP <- function(MLD.ph,
                                     MLD.pr,
                                     y50=1.19, # The y50 and d50 is forced
                                     d50=34.4, # and not editable in the online app. 
                                     beta.dNTCP, 
                                     sd.dNTCP,
                                     rfs = c(0,0,0,0), 
                                     rfs.beta = c(1,1,1,1),
                                     return.value = 'p', 
                                     showtext = TRUE, 
                                     xLim = 35, 
                                     forcexLim = FALSE
) {
  # Function for estimating the possible reduction in RP when using protons instead of photons
  # based on the possible variation in the betavalues for dNTCP taken from a cox proportional hazard model (CoxPHM)
  # MLD.ph = mean lung dose in Gy for PHOTON treatment
  # MLD.pr = mean lung dose in Gy for PROTON treatment
  # y50 = gamma50 for the logistic model (Appelt et al 2014 as standard)
  # d50 = dose50 for the logistic model (Appelt et al 2014 as standard)
  # beta.dNTCP = beta for the deltaNTCP (ie each percentage point reduction in dNTCP) from the CoxPHM.
  # sd.dNTCP = standard deviation of the beta.dNTCP from the coxPHM
  # rfs = Activation of riskfactors (included riskfactors are: 
  #       1. Pre-existing pulmonary co-morbidity, 
  #       2. mid or inferior tumor location, 
  #       3. current smoker, 
  #       4. old age (>63 yrs))
  # rfs.beta = beta for the riskfactos, NOTE: NOT THE HAZARDRATIO but the log(HR)
  # ximM = set the max of xlim to
  # return.value = (not specified yet include plot,data vector etc?)
  #
  
  #Calculate the summarized OR based for the risk factors
  # Here the NTCP for photon treatment isincluded as an OR?
  sumOR <- prod(exp(c(rfs.beta*rfs)))
  #calculate the summarized oddsratio for all the included riskfactors. 
  
  
  ntcp.ph <- 1/(1+sumOR^(-1)*exp(4*y50*(1-MLD.ph/d50))) # new probabilities calculated according to Applet et al 2014 (eq. 6)
  ntcp.pr <- 1/(1+sumOR^(-1)*exp(4*y50*(1-MLD.pr/d50))) # NTCP for protons
  
  dntcp <- ntcp.ph-ntcp.pr
  posb <- rnorm(1000,beta.dNTCP,sd.dNTCP) #Generate a distribtion of all possible Betas. 
  
  rpdist <- sort((1-ntcp.ph)^exp(dntcp*posb))
  rpdist <- 1-rpdist
  df <- data.frame(dist = sort(c(rpdist*100)),sort = rep(c("PossibleBeta"),each = 1000),ecdfx = as.vector(knots(ecdf(sort(rpdist*100)))),ecdfy = seq(0,1,length=1000))
  df$selt[df$dist >= ntcp.ph] <- "Photons"
  df$selt[df$dist < ntcp.ph] <- "Protons"
  
  
  p <- qplot(dist,data=df,geom="density",fill=sort,alpha = I(0.5))+
    geom_vline(aes(xintercept=ntcp.ph*100)) + geom_line(data=df,aes(x=ecdfx,y=ecdfy))
    if (showtext) {p <- p + geom_segment(aes(x = rpdist[500]*100,y = 0,xend = rpdist[500]*100,yend = 0.50)) + geom_segment(aes(x = 0,y = 0.5,xend = rpdist[500]*100,yend = 0.50)) }
    #geom_vline(aes(xintercept=ntcp.pr*100))+
  
    p <- p + theme(legend.position="none")
  p <- p + labs(x = "Probability for radiation pneumonitis [%]",y = expression(paste("Cumulative probability for exp{",beta,"*",Delta,"NTCP}")))
  if (forcexLim) {p <- p + xlim(0, xLim)}
  p <- p + ylim(0,1)
  txcor <- max(ggplot_build(p)$panel$ranges[[1]]$x.range)/100
  p <- p + annotate("text", x = ntcp.ph*100+txcor, y = 0.4, label = paste("Model based prediction of photon toxicity: ",spe_dec(ntcp.ph*100,2),"%",sep=""),srt = 270)
  if (showtext) {
  #p <- p + annotate("text", x = ntcp.pr*100+txcor, y = 0.4, label = paste("Protonprob logistic model: ",spe_dec(ntcp.pr*100,2),"%",sep=""),srt = 270)
  p <- p + annotate("text", x = rpdist[500]*100+txcor, y = 0.4, label = paste("Model based prediction of proton toxicity: ",spe_dec(sort(rpdist)[500]*100,2),"%",sep=""),srt = 270)
  }
  
  if (tolower(return.value)=='p'||tolower(return.value) == 'plot') {
    return(p)
  } else if (tolower(return.value)=='d'||tolower(return.value) == 'data') {
    return(df)
  }
  
}
#function for correcting the ecdf to allow for max-value other than 1. 
stat_ecdf2 <- function(mapping = NULL, data = NULL, geom = "step",
                       position = "identity", n = NULL, show.legend = NA,
                       inherit.aes = TRUE, minval=NULL, maxval=NULL,...) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatEcdf2,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    stat_params = list(n = n, minval=minval,maxval=maxval),
    params = list(...)
  )
}


StatEcdf2 <- ggproto("StatEcdf2", StatEcdf,
                     calculate = function(data, scales, n = NULL, minval=NULL, maxval=NULL, ...) {
                       df <- StatEcdf$calculate(data, scales, n, ...)
                       if (!is.null(minval)) { df$x[1] <- minval }
                       if (!is.null(maxval)) { df$x[length(df$x)] <- maxval }
                       df
                     }
)
#function for generation of specificying number of decimals. 
spe_dec <- function(x, k) {format(round(x, k), nsmall=k)}
#genORnames <- function(riskfactors){ c("Pre-existing pulmonary co-morbidity","Mid or inferior tumor location","Current smoker","Old age")[riskfactors==1]}

#plotfun <- function()
#server defines what kind of output you want and how the input should be used
server <- function(input,output){
  output$text1 <- renderPrint(
    {
      
    }
  )
  #the output is called output$test because I named the plotOutput "test" in the ui function above
  #I chose renderPlot because I want a plot as output, and I defined in ui that I want a plotOutput.
  output$test<-renderPlot({
    
    #rbind is used to merge refPatient and the input parameters (by calling function "fun" above) to a data frame containing two patients.
    newData<-fun(input$Comorb,input$TumorSite,input$Tobacco,input$Age,input$ORComorb,input$ORTloc,input$ORTob,input$ORAge,input$mldph,input$mldpr,input$betadNTCP,input$betadNTCPsd,input$xLim,input$forcexLim)
    riskf <- c(0,0,0,0)
    ORs <- c(1,1,1,1)
    #convert input parameters to riskfactor vector and logOR vector based on input. 
    if (newData$Comorb == "Yes") {
      riskf[1] = 1
      ORs[1] = newData$ORComorb
    }
    if (newData$TumorSite != "Superior") {
      riskf[2] = 1
      ORs[2] = newData$ORTloc
    }
    if (newData$Tobacco == "Current smoker") {
      riskf[3] = 1
      ORs[3] = newData$ORTob
    }
    if (newData$Age > 63) {
      riskf[4] = 1
      ORs[4] = newData$ORAge
    }
    logORs <- log(ORs) # note the OR should be fed with logOR as the combined OR is generated in the function getProtonProbReductionRP
    #generate plot with function above:
    p <- getProtonProbReductionRP(MLD.ph=newData$mldph,MLD.pr=newData$mldpr,beta.dNTCP = newData$betadNTCP,sd.dNTCP = newData$betadNTCPsd,rfs = riskf,rfs.beta = logORs,showtext = input$showtext,xLim = input$xLim,forcexLim = input$forcexLim)
  
    #this is the plot that is finally given to the renderPlot function
    p
    
    })
}

#defines what is input and output in my code
shinyApp(ui=ui, server=server)

#upload application:
#library(rsconnect)
#rsconnect::deployApp('/Users/Jonas/Dropbox/Jobb/PostDoc/P 1_1/Trial Webb/')