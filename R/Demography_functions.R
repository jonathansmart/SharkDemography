#' Create empty template for demography life history data
#' @description The demography functions in this package require a large amount of detailed
#'     life history information. This is provided to these functions as a multi-level list
#'     of various life history parameters with the class `Demography.inputs`. This function
#'     creates the template for this input data which can be filled in after it is created.
#' @param maturity.type The type of maturity estimates from the source life history study.
#'     must be one of `logistic - int/slope`, `logistic - a50/a95`, `normal` or `uniform`.
#' @param t0 Logical argument regarding whether the growth models included "t0" or "L0" as
#'     as a parameter. Default is `FALSE`
#' @examples
#' ######-----------
#' # Example code for Silky sharks
#' ######-----------
#'
#' silky_data <- create_data_input("logistic - int/slope", t0 = FALSE)
#' # Add growth data
#' silky_data$`growth`$model.type <- "logistic"
#' silky_data$growth$pars$Linf <- 268
#' silky_data$growth$pars$k <- 0.14
#' silky_data$growth$pars$L0 <- 82.7

#' silky_data$growth$se$Linf.se <- 5.8
#' silky_data$growth$se$k.se <- 0.006
#' silky_data$growth$se$L0.se <- 1.6
#' silky_data$growth$corr.matrix <- matrix(ncol = 3, nrow = 3,
#'                                         dimnames = list(c("Linf", "k", "L0"),c("Linf", "k", "L0")),
#'                                         data = c(1.0000000, -0.907188, 0.6233407,
#'                                                 -0.9071881,  1.0000000, -0.8572509,
#'                                                 0.6233407,-0.857250, 1.0000000))
#' # Add maturity data
#' silky_data$maturity$pars$intercept <- -15.90
#' silky_data$maturity$pars$slope <- 1.14
#' silky_data$maturity$se$intercept.se <- 2.78258
#' silky_data$maturity$se$slope.se <- 0.1971363
#' silky_data$maturity$corr.matrix <- matrix(ncol = 2, nrow = 2,
#'                                           dimnames = list(c("Intercept", "slope")
#'                                           ,c("Intercept", "slope")),
#'                                           data = c(1.0000000, -0.9922574,
#'                                                    -0.9922574, 1.0000000))
#' # max age lower bound
#' silky_data$max.age$min <- 28
#'
#' # Add fecundity info
#' silky_data$litter.size$mean <- 10
#' silky_data$litter.size$se <- 3
#' silky_data$gest.period <- 1
#' silky_data$repro.cycle <- 2
#'
#' # Add TL conversions (if available and required)
#' silky_data$Lt.type <- "TL"
#' silky_data$Lt.to.Wt$model.type <- "PCL"
#' silky_data$Lt.to.Wt$pars$a <- 2.73e-5
#' silky_data$Lt.to.Wt$pars$b <- 2.86
#'
#' silky_data$convert.TL$model.type <- "PCL"
#' silky_data$convert.TL$pars$a <- 2.08
#' silky_data$convert.TL$pars$b <- 1.32
#' @return A multi-level list of the class `Demography.inputs`
#' @export
#'
create_data_input <- function(maturity.type, t0 = FALSE){
  if(!maturity.type %in% c("logistic - int/slope",
                           "logistic - a50/a95",
                           "normal",
                           "uniform"))
    stop("Maturity must be defined as: either 'logistic - int/slope', 'logistic - a50/a95','normal or 'uniform'")

  data_list <- list()
  if(t0 == TRUE){
    data_list[["growth"]]<- list(model.type = NA,
                                 pars = list(Linf = NA, k = NA, t0 = NA ),
                                 se = list(Linf.se = NA, k.se = NA, t0.se = NA ),
                                 corr.matrix = NA)
  } else{
    data_list[["growth"]]<- list(model.type = NA,
                                 pars = list(Linf = NA, k = NA, L0 = NA ),
                                 se = list(Linf.se = NA, k.se = NA, L0.se = NA ),
                                 corr.matrix = NA)
  }

  if(maturity.type == "logistic - int/slope"){
    data_list[["maturity"]]<- list(model.type = "logistic - int/slope",
                                   pars = list(intercept = NA, slope = NA),
                                   se = list(intercept.se = NA, slope.se = NA),
                                   corr.matrix = NA)
  } else if(maturity.type == "logistic - a50/a95"){
    data_list[["maturity"]]<- list(model.type = "logistic - a50/a95",
                                   pars = list(a50 = NA, a95 = NA),
                                   se = list(a50.se = NA, a95.se = NA),
                                   corr.matrix = NA)
  }  else if(maturity.type == "normal"){
    data_list[["maturity"]]<- list(model.type = "normal",
                                   pars = list(mean = NA, se = NA))
  } else if(maturity.type == "uniform"){
    data_list[["maturity"]]<- list(model.type = "uniform",
                                   pars = list(min = NA, max = NA))
  } else {stop("maturity.type is unspecified")}

  data_list[["max.age"]] <- list(min  = NA, max = NA)

  data_list[["Mortality"]] <- list(mean = NA, se = NA)

  data_list[["Lt.type"]] <- "TL"
  data_list[["Lt.to.Wt"]] <- list(model.type = NA,
                                  pars = list(a = NA, b = NA))

  data_list[["convert.TL"]] <- list(model.type = NA,
                                    pars = list(a = NA, b = NA))

  data_list[["litter.size"]] <- list(mean = NA, se = NA)
  data_list[["gest.period"]] <- NA
  data_list[["repro.cycle"]] <- NA
  class(data_list) <- "Demography.inputs"

  return(data_list)
}


#' The base function for calculating demographic traits using a Leslie matrix method
#'
#' @description Using provided life history estimates from various sources, a demographic
#'     analysis is performed using Leslie matrix methods. This function can be run as
#'     stand alond and will given the outputs of a stochastic matrix model given the provided
#'     life history parameters and their associated error. However, the main purpose of this
#'     function is to be used in other simulation functions provided in this package. In
#'     each use this function will calculate: \describe{
#' \item{lambda}{The finite rate of population growth}
#' \item{R0}{The reproductive value of the population}
#' \item{G}{The generation length of the population}
#' \item{elast.fecund}{The mean elasticity of the fecundity elements}
#' \item{elast.juv.survival}{The mean elasticity of the juvenile survivoship elements}
#' \item{juv.ratio}{The ratio of fecundity to juvenile elasticities}
#' \item{adult.ratio}{The ratio of fecundity to adult elasticities}
#' \item{M.estimator}{The natural mortality estimator randomly used in this analysis}
#' }
#'
#' @param data A multi-level list of the class `Demography.inputs` produced from the
#'     `create_data_input` function and then manually completed.
#' @param AAFC Age-at-first-capture which can be specified by the user or automated in the
#'     `Simulate_AAFC` function. Must be an integer age which can include zero to indicate
#'     the availability of the population to capture from birth.
#' @param F. The instantaneous rate of fishing mortality 'F'. This will be applied to all
#'      ages available to capture as defined by either the AALC or AAFC arguments.
#' @param AALC Age-at-last-capture which can be specified by the user or automated in the
#'     `Simulate_AALC` function. Must be an integer age which can include zero.
#' @param M.estimators Any specific natural mortality estimators to be included in the analysis.
#'     Only one will be used in each run which is randomly selected. In the simulation functions
#'     these will be varied among simulations. Must be a single estimator or a vector of estimators.
#'     These can include: "Pet.Wro","Jensen.mat","Chen.Yuan","Then_hoenig","Then_pauly", "
#'     Jensen.mat","Charnov" or "Chen.Want". If none are specified then all applicable estimators
#'     could be chosen.
#'
#' @return A list of parameters including: \describe{
#' \item{lambda}{The finite rate of population growth}
#' \item{R0}{The reproductive value of the population}
#' \item{G}{The generation length of the population}
#' \item{elast.fecund}{The mean elasticity of the fecundity elements}
#' \item{elast.juv.survival}{The mean elasticity of the juvenile survivoship elements}
#' \item{juv.ratio}{The ratio of fecundity to juvenile elasticities}
#' \item{adult.ratio}{The ratio of fecundity to adult elasticities}
#' \item{M.estimator}{The natural mortality estimator randomly used in this analysis}
#' }
#' @examples
#' # load Silky shark data produced by create_data_input()
#' # Type `?create_data_input()` for details
#' data("Silky_data")
#'
#' # run a single Leslie Matrix analysis with random draws from biological data
#' # distributions. Use `Simulate_demography()` to run Monte Carlo Simulations
#' # using this function
#'
#' Calculate_demography(Silky_data)
#'
#' @export
#' @import popbio
#'
Calculate_demography <- function(data, AAFC = NULL, F. = 0, AALC = NULL, M.estimators = NULL){

  ### Error messages --------------


  if(!inherits(data,"Demography.inputs")) stop("Input data is not the correct based on the demography data list")

  if(!data$growth$model.type %in%  c("Von Bertalanffy", "logistic", "Gompertz"))
    stop("Growth model type not correctly specified")

  possible_estimators <- c("Pet.Wro","Jensen.mat","Chen.Yuan","Then_hoenig","Then_pauly", "Jensen","Charnov","Chen.Want")
  if(any(!M.estimators %in% possible_estimators)) stop("One or more M estimators are not correctly specified")
  if(is.null(M.estimators)){
    M.estimators <- possible_estimators
  }

  # remove all mortality estimates that require VB params if VB is not being used
  if(data$growth$model.type != "Von Bertalanffy") M.estimators <- M.estimators[which(!M.estimators %in% c("Chen.Want", "Charnov", "Jensen", "Then_pauly"))]
  if(length(M.estimators) == 0) stop("Specified mortality estimators could not be applied")

  if(is.null(AALC) & is.null(AAFC) & F. > 0){
    warning("No fishing mortality applied as neither AAFC or AALC are specified")
  }

  ### Estimate growth ----------------

  # If L0 is being used
  if(!any(names(data$growth$pars) %in% "t0")){

    Linf.se <- data$growth$se$Linf.se
    k.se <- data$growth$se$k.se
    L0.se <- data$growth$se$L0.se

    # if the parameter correlation is provided used a MVN dist
    if(all(!is.na(data$growth$corr.matrix))){
      e.cov <- data$growth$corr.matrix*as.matrix(c(Linf.se, k.se, L0.se))%*%t(as.matrix(c(Linf.se, k.se, L0.se)))
      growth_pars <- MASS::mvrnorm(1, mu = c(data$growth$pars$Linf, data$growth$pars$k, data$growth$pars$L0), Sigma = e.cov)
      Linf <-as.numeric(growth_pars["Linf"])
      k <- as.numeric(growth_pars["k"])
      L0 <- as.numeric(growth_pars["L0"])
    } else if(all(!is.na(Linf.se),!is.na(k.se),!is.na(L0.se))){
      #If MVN can't be used, then vary parameters based on normal dists
      Linf <- rnorm(1, data$growth$pars$Linf, Linf.se)
      k <- rnorm(1, data$growth$pars$k, k.se)
      L0 <- rnorm(1, data$growth$pars$L0, L0.se)
    } else{
      # if no errors available then just use the parameters
      warning("Growth is deterministic as no growth parameter errors were provided")
      Linf <- data$growth$pars$Linf
      k <- data$growth$pars$k
      L0 <- data$growth$pars$L0
    }
  } else {
    # If t0 is being used the same approach is applied but without L0
    Linf.se <- data$growth$se$Linf.se
    k.se <- data$growth$se$k.se
    t0.se <- data$growth$se$t0.se

    if(all(!is.na(data$growth$corr.matrix))){
      e.cov <- data$growth$corr.matrix*as.matrix(c(Linf.se, k.se, t0.se))%*%t(as.matrix(c(Linf.se, k.se, t0.se)))
      growth_pars <- MASS::mvrnorm(1, mu = c(data$growth$pars$Linf, data$growth$pars$k, data$growth$pars$t0), Sigma = e.cov)
      Linf <-as.numeric(growth_pars["Linf"])
      k <- as.numeric(growth_pars["k"])
      t0 <- as.numeric(growth_pars["t0"])
    } else if(all(!is.na(Linf.se),!is.na(k.se),!is.na(t0.se))){
      Linf <- rnorm(1, data$growth$pars$Linf, Linf.se)
      k <- rnorm(1, data$growth$pars$k, k.se)
      t0 <- rnorm(1, data$growth$pars$t0, t0.se)
    } else{
      warning("Growth is deterministic as no growth parameter errors were provided")
      Linf <- data$growth$pars$Linf
      k <- data$growth$pars$k
      t0 <- data$growth$pars$t0
    }
  }

  # Set max age lower bound
  if(!is.na(data$max.age$min)){
    max.age.lower.bound <- data$max.age$min
  } else {
    stop("No min age specified")
  }



  ### Age at 99% of linf (or 95% of Von Bertalanffy used). Used as upper max age
  if(data$`growth`$model.type == "logistic"){
    tmax <- log(((Linf*0.99)*(Linf-L0))/(L0*(Linf - (Linf*0.99))))/k
  } else if(data$`growth`$model.type == "Von Bertalanffy"){
    if(!any(names(data$growth$pars) %in% "t0")){
      tmax <- log(1-((Linf*0.95)-L0)/(Linf-L0))/-k
    } else {
      tmax <- (log(1-((Linf*0.95)/Linf))/-k)+t0
    }
  } else if(data$`growth`$model.type == "Gompertz"){
    tmax <- log(log(Linf/L0)/(log(Linf/L0)-log((Linf*0.99)/L0)))/k
  } else { stop("model.type incorrectely specific for growth")}

  # Set max age upper bound if pre-specified
  if(!is.na(data$max.age$max)){
    max.age.upper.bound <- data$max.age$max
  } else {
    max.age.upper.bound <- 1000 # To make sure its ignored if unspecified
  }

  # Set max age based on available parameters
  if(max.age.lower.bound > tmax | is.na(tmax)){
    max.age <- max.age.lower.bound
  }else if(max.age.upper.bound < tmax){
    max.age <- round(runif(1,max.age.lower.bound,max.age.upper.bound))
  } else{
    max.age <- round(runif(1,max.age.lower.bound,tmax))
  }

  Age <- 0:max.age

  # Length-at-age estimates
  if(data$`growth`$model.type == "logistic"){
    Lt <- Linf*L0*exp(k*Age)/(Linf+L0*(exp(k*Age)-1))
  } else if(data$`growth`$model.type == "Von Bertalanffy"){
    if(!any(names(data$growth$pars) %in% "t0")){
      Lt <- Linf-(Linf-L0)*(exp(-k*Age))
    } else {
      Lt <- Linf*(1-exp(-k*(Age-t0)))
    }
  } else if(data$`growth`$model.type == "Gompertz"){
    Lt <- L0*(exp(log(Linf/L0)*(1-exp(-k*Age))))
  } else {
    stop("model.type incorrectely specific for growth")
  }

  #----End growth calcs------------------

  if(!is.na(data$gest.period)){
    gest.period <- data$gest.period
  } else {
    stop("No gestation period specified")
  }
  if(!is.na(data$repro.cycle)){
    repro.cycle <- data$repro.cycle
  } else {
    stop("No reproductive cycle (periodicity) specified")
  }


  #Maturity ogive
  if(data$maturity$model.type == "logistic - int/slope"){
    mat.cov <- data$maturity$corr.matrix*as.matrix(c(data$maturity$se$intercept.se, data$maturity$se$slope.se))%*%
      t(as.matrix(c(data$maturity$se$intercept.se, data$maturity$se$slope.se)))

    Mat_pars <- MASS::mvrnorm(1, mu = c(data$maturity$pars$intercept, data$maturity$pars$slope), Sigma = mat.cov)
    Mat_intercept <- as.numeric(Mat_pars[1])
    Mat_slope <-as.numeric(Mat_pars[2])
    Maturity_ogive <- round(1/(1+exp(-(Age*Mat_slope+Mat_intercept))),2)
    age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
    age_at_first_repro <- round(gest.period+age_at_maturity)
    Fecundity_ogive <- round(1/(1+exp(-((Age-gest.period)*Mat_slope+Mat_intercept))),2)
  } else if(data$maturity$model.type == "logistic - a50/a95"){
    Mat_a50 <- rnorm(1,data$maturity$pars$a50, data$maturity$se$a50.se)
    Mat_a95 <- rnorm(1,data$maturity$pars$a95, data$maturity$se$a95.se)
    Maturity_ogive <- round(1*(1+exp(-log(19)*((Age-Mat_a50)/(Mat_a95-Mat_a50))))^-1,2)
    age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
    age_at_first_repro <- round(gest.period+age_at_maturity)
    Fecundity_ogive <- round( 1*(1+exp(-log(19)*(((Age-gest.period)-Mat_a50)/(Mat_a95-Mat_a50))))^-1,2)
  } else if(data$maturity$model.type == "uniform"){
    age_at_maturity <- runif(1,data$maturity$pars$min, data$maturity$pars$max)
    age_at_first_repro <- round(gest.period+age_at_maturity)
  } else if(data$maturity$model.type == "normal"){
    age_at_maturity <- rnorm(1,data$maturity$pars$mean, data$maturity$pars$se)
    age_at_first_repro <- round(gest.period+age_at_maturity)
  } else {
    stop("Maturity model type is not correctly specified.")
  }

  # If the age at first repro is older than the max age in a sim, reassign it.
  if(max.age < age_at_first_repro){
    age_at_first_repro <- max.age - 1
    age_at_maturity <- age_at_first_repro - gest.period
  }



  # litter size which is transformed to number of female pups per female per year.

  if(!any(is.na(data$litter.size$mean), is.na(data$litter.size$se))){
    litter.size <- rnorm(1, data$litter.size$mean, data$litter.size$se)
  } else {
    stop("No litter size mean and/or se specified")
  }
  n.fem.pups <- litter.size/2
  n.pups.year <- n.fem.pups/repro.cycle

  ## estimate mortality


  if(all(!is.na(data$Mortality$mean), !is.na(data$Mortality$se))){
    M.table <- matrix(nrow = length(Age), ncol = 1)
    colnames(M.table)<-c("Pre-specified")
    M.table[,"Pre-specified"] <- rnorm(1,data$Mortality$mean, data$Mortality$se)
    M.select <- 1
    Mortality <- M.table[,M.select]
  } else {
    M.table <- matrix(nrow = length(Age), ncol = length(M.estimators))
    colnames(M.table)<-c(M.estimators)
    if("Pet.Wro" %in% M.estimators){

      # If the length weight model type is the same as the length type then no conversion
      # is needed
      if(data$Lt.to.Wt$model.type == data$Lt.type){
        Length.for.wt.conv <- Lt
      } else{

        # Otherwise perform, conversion as necessary. Currently, this converts TL to
        # PCL or FL only. If PCL or FL is used for length in the growth model then typically
        # They are also used in the length weight relationship and don't need converting.
        if(data$convert.TL$model.type == "PCL" & data$Lt.to.Wt$model.type == "PCL"){
          Length.for.wt.conv<-(Lt-data$convert.TL$pars$a)/data$convert.TL$pars$b
        } else if(data$convert.TL$model.type == "FL" & data$Lt.to.Wt$model.type == "FL"){
          Length.for.wt.conv<-(Lt*data$convert.TL$pars$a)-data$convert.TL$pars$b
        } else if(data$Lt.to.Wt$model.type == "TL"){
          Length.for.wt.conv <- Lt
        } else {
          Length.for.wt.conv <- NA # correct length cannot be calculated so return NA to nullify this method
        }
      }
      Wt <- (data$Lt.to.Wt$pars$a)*Length.for.wt.conv^data$Lt.to.Wt$pars$b

      M.table[,"Pet.Wro"]<-((Wt^-0.25)*1.92)*0.2

    }
    if("Jensen.mat" %in% M.estimators){
      M.table[,"Jensen.mat"] <- 1.65/age_at_maturity}
    if("Chen.Yuan" %in% M.estimators){
      M.table[,"Chen.Yuan"] <-exp(1.46-1.01*log(tmax))}
    if("Then_hoenig" %in% M.estimators){
      M.table[,"Then_hoenig"] <- 4.899*tmax^-0.916}
    if("Then_pauly" %in% M.estimators &
       data$growth$model.type == "Von Bertalanffy"){
      Then_pauly <- 4.118*k^0.73*Linf^-0.33
      M.table <- cbind(M.table, Then_pauly = Then_pauly)
    }
    if("Jensen" %in% M.estimators&
       data$growth$model.type == "Von Bertalanffy"){
      M.table[,"Jensen"] <- 1.5*k
    }
    if("Charnov" %in% M.estimators&
       data$growth$model.type == "Von Bertalanffy"){
      M.table[,"Charnov"] <- k*(Lt/Linf)^-1.5
    }
    if("Chen.Want" %in% M.estimators&
       data$growth$model.type == "Von Bertalanffy"){

      if(!any(names(data$growth$pars) %in% "t0")){
        t0<- (1/k)*log((Linf-L0)/Linf)
      }

      # Chen and Wanatabe method
      tM <- -(1/k)*log(1-exp(k*t0))+t0
      a0 <- 1-exp(-k*(tM-t0))
      a1 <- k*exp(-k*(tM-t0))
      a2 <- -0.5*k^2*exp(-k*(tM-t0))
      CW  <- NULL
      for(j in 0:max.age){
        if(j <=tM){cw <- k/(1-exp(-k*(j-t0)))}
        else{cw <- k/(a0+a1*(j-tM)+a2*(j-tM)^2)}
        CW <- cbind(CW,cw)
      }

      M.table[,"Chen.Want"] <- as.vector(CW)
    }

    if(ncol(M.table) == 1){
      M.select <- sample(x=ncol(M.table),size=1)
      Mortality <- M.table[,M.select]
    } else {
      M.table <- M.table[,colSums(is.na(M.table))<nrow(M.table)]
      M.select <- sample(x=ncol(M.table),size=1)
      Mortality <- M.table[,M.select]
    }
  }

  # Estimate Survivorship from M and F estimates. If an AAFC or AALC is set then F is only applied to ages
  # older than the AAFC or younger than AALC. If both are zero then F is not included.
  # An error (specified earlier) is thrown if both are non-zero

  # if(is.null(AALC) & is.null(AAFC)){
  #   survivorship <- exp(-(Mortality))
  # } else if(is.null(AALC) & !is.null(AAFC)){
  #   survivorship <- exp(-(Mortality))
  #   for(j in seq_along(Age)){
  #     if(Age[j]>=AAFC){survivorship[j] <- exp(-(Mortality[j]+F.))}
  #   }
  # } else if(!is.null(AALC) & is.null(AAFC)){
  #   survivorship <- exp(-(Mortality))
  #   for(j in seq_along(Age)){
  #     if(Age[j]<=AALC){survivorship[j] <- exp(-(Mortality[j]+F.))}
  #   }
  # } else {
  #   stop("AAFC and AALC cannot be used simulatneously")
  # }

  # Code update. An AALC and AAFC can now be specified which is used in the Harvest
  # slot simulation function

  if(is.null(AALC) & is.null(AAFC)){
    survivorship <- exp(-(Mortality))
  } else if(is.null(AALC) & !is.null(AAFC)){
    survivorship <- exp(-(Mortality))
    for(j in seq_along(Age)){
      if(Age[j]>=AAFC){survivorship[j] <- exp(-(Mortality[j]+F.))}
    }
  } else if(!is.null(AALC) & is.null(AAFC)){
    survivorship <- exp(-(Mortality))
    for(j in seq_along(Age)){
      if(Age[j]<=AALC){survivorship[j] <- exp(-(Mortality[j]+F.))}
    }
  } else {
    survivorship <- exp(-(Mortality))
    for(j in seq_along(Age)){
      if(Age[j]<=AALC & Age[j]>=AAFC){survivorship[j] <- exp(-(Mortality[j]+F.))}
    }

    #stop("AAFC and AALC cannot be used simulatneously")
  }


  # Add a zero at the end of the sx vector to identify that no individuals survive
  # the final age class
  # sx <- survivorship
  # sx[max.age+1] <- 0
  #
  # # Calculate fecundity at age based on knife edge or continous dists.
  #
  # if(data$maturity$model.type %in% c("normal", "uniform")){
  #   mx <- rep(0, length(Age))
  #   mx[age_at_first_repro:max.age] <- n.pups.year
  # } else{
  #   mx <- Fecundity_ogive * n.pups.year
  # }
  #
  # # line up vector so sx and mx are offset by 1 so that eventually fx = mx*sx+1...
  # a <- mx[-1]
  #
  # # add zero to the end so both vectors are the same length again.
  # MX <- append(a,0,after=length(a))
  #
  # # Create top row of matrix
  # fx <- MX*sx
  #
  # # create a matrix entirely of zeros
  # Sm <- matrix(rep(0),nrow=length(survivorship),ncol=length(survivorship))
  #
  # #insert survival along the diagonal
  # diag(Sm) <- survivorship
  #
  # # join fecundity vector to survival matrix without row names
  # A <- cbind(rbind(fx,Sm,deparse.level = 0),0)

  sx <- survivorship

  # sx <- append(sx,0,after = Inf)

  # Calculate fecundity at age based on knife edge or continous dists.

  if(data$maturity$model.type %in% c("normal", "uniform")){
    mx <- rep(0, length(Age))
    mx[age_at_first_repro:max.age+1] <- n.pups.year
  } else{
    mx <- Fecundity_ogive * n.pups.year
  }

  # line up vector so sx and mx are offset by 1 so that eventually fx = mx*sx+1...
  a <- mx[-1]

  # add zero to the end so both vectors are the same length again.
  MX <- append(a,a[max.age],after=length(a))

  # Create top row of matrix
  fx <- MX*sx

  # create a matrix entirely of zeros
  Sm <- matrix(rep(0),nrow=length(survivorship),ncol=length(survivorship))

  #insert survival along the diagonal
  diag(Sm) <- survivorship

  # join fecundity vector to survival matrix without row names
  A <- cbind(rbind(fx,Sm,deparse.level = 0),0)


  ############ Matrix projection ##########

  # calculate lambda
  lambda<- eigen.analysis(A)$lambda

  r <- log(lambda)

  #convert sx to lx values for R0 and G calculations
  lx <- NULL
  for(j in 3:length(Age)){
    lx[1] <- 1
    lx[2] <- survivorship[1]
    lx[j] <- lx[j-1]*survivorship[j-1]
  }

  # Net reproductive rate
  R0 <- sum(mx*lx)

  #Generation time
  G <- sum(Age*(exp(-r*Age)*(mx*lx)))

  ### Elasticities-----------------
  # Sometimes the matrix won't invert. When this happens, assign NA's to elasticities
  elast <- eigen.analysis(A)$elasticities
  elast.fecund <- sum(elast[1,], na.rm = TRUE)
  if( elast.fecund> 1) { elast.fecund <- NA}
  elast.survival <- apply(elast[2:max.age,],2, sum, na.rm = TRUE)
  if( any(elast.survival> 1)) { elast.survival <- rep(NA, length(0:max.age))}

  elast.juv.survival <- sum(elast.survival[1:age_at_first_repro-1], na.rm = TRUE)
  elast.adult.survival <- sum(elast.survival[age_at_first_repro:max.age], na.rm = TRUE)

  juv.ratio <- elast.juv.survival/elast.fecund

  adult.ratio <- elast.adult.survival/elast.fecund

  # IF A is singular then no elasticity should be estimated. return NA's when this happens
  if(any(is.na(elast.fecund), is.na(elast.juv.survival), is.na(elast.adult.survival))){
    elast.fecund <- NA
    elast.juv.survival <- NA
    elast.adult.survival <- NA
    juv.ratio <- NA
    adult.ratio <- NA
  } else if(any(c(elast.fecund, elast.juv.survival, elast.adult.survival) == 0)){
    elast.fecund <- NA
    elast.juv.survival <- NA
    elast.adult.survival <- NA
    juv.ratio <- NA
    adult.ratio <- NA
  }

  Results <- list(lambda = lambda,
                  R0 = R0 ,
                  G = G,
                  elast.fecund = elast.fecund,
                  elast.juv.survival = elast.juv.survival,
                  elast.adult.survival = elast.adult.survival,
                  juv.ratio = juv.ratio,
                  adult.ratio = adult.ratio,
                  M.estimator = names(as.data.frame(M.table))[M.select])

  # IF the elasticities don't calculate, return all NA's as the results are probably suspect
  # In some instances the lambda is very large or small which causes issues.
  if(any(is.na(elast.fecund), is.na(elast.juv.survival), is.na(elast.adult.survival))){
    Results$lambda <- NA
    Results$R0 <- NA
    Results$G <- NA
    Results$elast.fecund <- NA
    Results$elast.juv.survival <- NA
    Results$elast.adult.survival <- NA
    Results$juv.ratio <- NA
    Results$adult.ratio <- NA
  }

  return(Results)

}


#' Estimate the left and right eignvectors (Reproductive value: v and Stable age distribution: w)
#'
#'@description This function performs a Monte Carlo simulation analysis which determines
#'     the distributions of the left and right eigenvector. To facilitate this, maximum age
#'     is fixed as the minimum value + 20\% for the species.
#' @param n The number of simulations to be run
#' @param data A multi-level list of the class `Demography.inputs` produced from the
#'     `create_data_input` function and then manually completed.
#' @param M.estimators Any specific natural mortality estimators to be included in the analysis.
#'      Must be a single estimator or a vector of estimators. These can include: "Pet.Wro",
#'     "Jensen.mat","Chen.Yuan","Then_hoenig","Then_pauly", "Jensen.mat","Charnov" or
#'     "Chen.Want". If none are specified then all applicable estimators could be chosen.
#'
#' @return A dataframe with the mean and 95\% quantiles for each eigenvector for each age
#'     class
#' @examples
#' # load Silky shark data produced by create_data_input()
#' # Type `?create_data_input()` for details
#' data("Silky_data")
#'
#' # Run function to get Eigenvectors from  Monte Carlo simulations for
#' # all available natural mortality estimators. Set n = at least 1000 for full
#' # analysis but use n = 100 for testing
#'
#' Estimate_eigenvectors(n = 100, Silky_data)
#' @export
Estimate_eigenvectors <- function(n = 1000, data, M.estimators = NULL){

  ## Error messages ------------------

  if(!inherits(data,"Demography.inputs")) stop("Input data is not the correct based on the demography data list")

  if(!data$growth$model.type %in%  c("Von Bertalanffy", "logistic", "Gompertz"))
    stop("Growth model type not correctly specified")

  possible_estimators <- c("Pet.Wro","Jensen.mat","Chen.Yuan","Then_hoenig","Then_pauly", "Jensen","Charnov","Chen.Want")
  if(any(!M.estimators %in% possible_estimators)) stop("One or more M estimators are not correctly specified")
  if(is.null(M.estimators)){
    M.estimators <- possible_estimators
  }

  # remove all mortality estimates that require VB params if VB is not being used
  if(data$growth$model.type != "Von Bertalanffy") M.estimators <- M.estimators[which(!M.estimators %in% c("Chen.Want", "Charnov", "Jensen", "Then_pauly"))]
  if(length(M.estimators) == 0) stop("Specified mortality estimators could not be applied")

  #--------------------------

  # Empty containers to fill
  # M.table <- matrix(nrow = n, ncol = length(M.estimators))
  # colnames(M.table)<- M.estimators
  # P_W_table <-  matrix(nrow = 100, ncol = n)
  # Charnov.table <- matrix(nrow = 100, ncol = n)
  # Chen.Want.table <- matrix(nrow = 100, ncol = n)

  if(is.na(data$max.age$min)) stop("No min age specified")
  Stable.age.results <- matrix(nrow = round(data$max.age$min*1.2) + 2, ncol = n)
  Repro.val.results <- matrix(nrow = round(data$max.age$min*1.2) + 2, ncol = n)


  pb <- txtProgressBar(title="Demography progress", label="0% done", min=0, max=100, initial=0, style = 3)
  for(i in 1:n){
    info <- sprintf("%d%% done", round((i/n)*100))
    setTxtProgressBar(pb, (i/n)*100, label=info)
    # If L0 is being used
    if(!any(names(data$growth$pars) %in% "t0")){

      Linf.se <- data$growth$se$Linf.se
      k.se <- data$growth$se$k.se
      L0.se <- data$growth$se$L0.se

      # if the parameter correlation is provided used a MVN dist
      if(all(!is.na(data$growth$corr.matrix))){
        e.cov <- data$growth$corr.matrix*as.matrix(c(Linf.se, k.se, L0.se))%*%t(as.matrix(c(Linf.se, k.se, L0.se)))
        growth_pars <- MASS::mvrnorm(1, mu = c(data$growth$pars$Linf, data$growth$pars$k, data$growth$pars$L0), Sigma = e.cov)
        Linf <-as.numeric(growth_pars["Linf"])
        k <- as.numeric(growth_pars["k"])
        L0 <- as.numeric(growth_pars["L0"])
      } else if(all(!is.na(Linf.se),!is.na(k.se),!is.na(L0.se))){
        #If MVN can't be used, then vary parameters based on normal dists
        Linf <- rnorm(1, data$growth$pars$Linf, Linf.se)
        k <- rnorm(1, data$growth$pars$k, k.se)
        L0 <- rnorm(1, data$growth$pars$L0, L0.se)
      } else{
        # if no errors available then just use the parameters
        warning("Growth is deterministic as no growth parameter errors were provided")
        Linf <- data$growth$pars$Linf
        k <- data$growth$pars$k
        L0 <- data$growth$pars$L0
      }
    } else {
      # If t0 is being used the same approach is applied but without L0
      Linf.se <- data$growth$se$Linf.se
      k.se <- data$growth$se$k.se
      t0.se <- data$growth$se$t0.se

      if(all(!is.na(data$growth$corr.matrix))){
        e.cov <- data$growth$corr.matrix*as.matrix(c(Linf.se, k.se, t0.se))%*%t(as.matrix(c(Linf.se, k.se, t0.se)))
        growth_pars <- MASS::mvrnorm(1, mu = c(data$growth$pars$Linf, data$growth$pars$k, data$growth$pars$t0), Sigma = e.cov)
        Linf <-as.numeric(growth_pars["Linf"])
        k <- as.numeric(growth_pars["k"])
        t0 <- as.numeric(growth_pars["t0"])
      } else if(all(!is.na(Linf.se),!is.na(k.se),!is.na(t0.se))){
        Linf <- rnorm(1, data$growth$pars$Linf, Linf.se)
        k <- rnorm(1, data$growth$pars$k, k.se)
        t0 <- rnorm(1, data$growth$pars$t0, t0.se)
      } else{
        warning("Growth is deterministic as no growth parameter errors were provided")
        Linf <- data$growth$pars$Linf
        k <- data$growth$pars$k
        t0 <- data$growth$pars$t0
      }
    }

    # Set max age lower bound
    if(!is.na(data$max.age$min)){
      max.age <- tmax <- round(data$max.age$min*1.2)
    } else {
      stop("No min age specified")
    }

    Age <- 0:max.age

    # Length-at-age estimates
    if(data$`growth`$model.type == "logistic"){
      Lt <- Linf*L0*exp(k*Age)/(Linf+L0*(exp(k*Age)-1))
    } else if(data$`growth`$model.type == "Von Bertalanffy"){
      if(!any(names(data$growth$pars) %in% "t0")){
        Lt <- Linf-(Linf-L0)*(exp(-k*Age))
      } else {
        Lt <- Linf*(1-exp(-k*(Age-t0)))
      }
    } else if(data$`growth`$model.type == "Gompertz"){
      Lt <- L0*(exp(log(Linf/L0)*(1-exp(-k*Age))))
    } else {
      stop("model.type incorrectely specific for growth")
    }

    #----End growth calcs------------------

    if(!is.na(data$gest.period)){
      gest.period <- data$gest.period
    } else {
      stop("No gestation period specified")
    }
    if(!is.na(data$repro.cycle)){
      repro.cycle <- data$repro.cycle
    } else {
      stop("No reproductive cycle (periodicity) specified")
    }


    #Maturity ogive
    if(data$maturity$model.type == "logistic - int/slope"){
      mat.cov <- data$maturity$corr.matrix*as.matrix(c(data$maturity$se$intercept.se, data$maturity$se$slope.se))%*%
        t(as.matrix(c(data$maturity$se$intercept.se, data$maturity$se$slope.se)))

      Mat_pars <- MASS::mvrnorm(1, mu = c(data$maturity$pars$intercept, data$maturity$pars$slope), Sigma = mat.cov)
      Mat_intercept <- as.numeric(Mat_pars[1])
      Mat_slope <-as.numeric(Mat_pars[2])
      Maturity_ogive <- round(1/(1+exp(-(Age*Mat_slope+Mat_intercept))),2)
      age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
      age_at_first_repro <- round(gest.period+age_at_maturity)
      Fecundity_ogive <- round(1/(1+exp(-((Age-gest.period)*Mat_slope+Mat_intercept))),2)
    } else if(data$maturity$model.type == "logistic - a50/a95"){
      Mat_a50 <- rnorm(1,data$maturity$pars$a50, data$maturity$se$a50.se)
      Mat_a95 <- rnorm(1,data$maturity$pars$a95, data$maturity$se$a95.se)
      Maturity_ogive <- round(1*(1+exp(-log(19)*((Age-Mat_a50)/(Mat_a95-Mat_a50))))^-1,2)
      age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
      age_at_first_repro <- round(gest.period+age_at_maturity)
      Fecundity_ogive <- round( 1*(1+exp(-log(19)*(((Age-gest.period)-Mat_a50)/(Mat_a95-Mat_a50))))^-1,2)
    } else if(data$maturity$model.type == "uniform"){
      age_at_maturity <- runif(1,data$maturity$pars$min, data$maturity$pars$max)
      age_at_first_repro <- round(gest.period+age_at_maturity)
    } else if(data$maturity$model.type == "normal"){
      age_at_maturity <- rnorm(1,data$maturity$pars$mean, data$maturity$pars$se)
      age_at_first_repro <- round(gest.period+age_at_maturity)
    } else {
      stop("Maturity model type is not correctly specified.")
    }

    # If the age at first repro is older than the max age in a sim, reassign it.
    if(max.age < age_at_first_repro){
      age_at_first_repro <- max.age - 1
      age_at_maturity <- age_at_first_repro - gest.period
    }



    # litter size which is transformed to number of female pups per female per year.

    if(!any(is.na(data$litter.size$mean), is.na(data$litter.size$se))){
      litter.size <- rnorm(1, data$litter.size$mean, data$litter.size$se)
    } else {
      stop("No litter size mean and/or se specified")
    }
    n.fem.pups <- litter.size/2
    n.pups.year <- n.fem.pups/repro.cycle

    ## estimate mortality


    if(all(!is.na(data$Mortality$mean), !is.na(data$Mortality$se))){
      M.table <- matrix(nrow = length(Age), ncol = 1)
      colnames(M.table)<-c("Pre-specified")
      M.table[,"Pre-specified"] <- rnorm(1,data$Mortality$mean, data$Mortality$se)
      M.select <- 1
      Mortality <- M.table[,M.select]
    } else {
      M.table <- matrix(nrow = length(Age), ncol = length(M.estimators))
      colnames(M.table)<-c(M.estimators)
      if("Pet.Wro" %in% M.estimators){

        # If the length weight model type is the same as the length type then no conversion
        # is needed
        if(data$Lt.to.Wt$model.type == data$Lt.type){
          Length.for.wt.conv <- Lt
        } else{

          # Otherwise perform, conversion as necessary. Currently, this converts TL to
          # PCL or FL only. If PCL or FL is used for length in the growth model then typically
          # They are also used in the length weight relationship and don't need converting.
          if(data$convert.TL$model.type == "PCL" & data$Lt.to.Wt$model.type == "PCL"){
            Length.for.wt.conv<-(Lt-data$convert.TL$pars$a)/data$convert.TL$pars$b
          } else if(data$convert.TL$model.type == "FL" & data$Lt.to.Wt$model.type == "FL"){
            Length.for.wt.conv<-(Lt*data$convert.TL$pars$a)-data$convert.TL$pars$b
          } else if(data$Lt.to.Wt$model.type == "TL"){
            Length.for.wt.conv <- Lt
          } else {
            Length.for.wt.conv <- NA # correct length cannot be calculated so return NA to nullify this method
          }
        }
        Wt <- (data$Lt.to.Wt$pars$a)*Length.for.wt.conv^data$Lt.to.Wt$pars$b

        M.table[,"Pet.Wro"]<-((Wt^-0.25)*1.92)*0.2

      }
      if("Jensen.mat" %in% M.estimators){
        M.table[,"Jensen.mat"] <- 1.65/age_at_maturity}
      if("Chen.Yuan" %in% M.estimators){
        M.table[,"Chen.Yuan"] <-exp(1.46-1.01*log(tmax))}
      if("Then_hoenig" %in% M.estimators){
        M.table[,"Then_hoenig"] <- 4.899*tmax^-0.916}
      if("Then_pauly" %in% M.estimators &
         data$growth$model.type == "Von Bertalanffy"){
        Then_pauly <- 4.118*k^0.73*Linf^-0.33
        M.table <- cbind(M.table, Then_pauly = Then_pauly)
      }
      if("Jensen" %in% M.estimators&
         data$growth$model.type == "Von Bertalanffy"){
        M.table[,"Jensen"] <- 1.5*k
      }
      if("Charnov" %in% M.estimators&
         data$growth$model.type == "Von Bertalanffy"){
        M.table[,"Charnov"] <- k*(Lt/Linf)^-1.5
      }
      if("Chen.Want" %in% M.estimators&
         data$growth$model.type == "Von Bertalanffy"){

        if(!any(names(data$growth$pars) %in% "t0")){
          t0<- (1/k)*log((Linf-L0)/Linf)
        }

        # Chen and Wanatabe method
        tM <- -(1/k)*log(1-exp(k*t0))+t0
        a0 <- 1-exp(-k*(tM-t0))
        a1 <- k*exp(-k*(tM-t0))
        a2 <- -0.5*k^2*exp(-k*(tM-t0))
        CW  <- NULL
        for(j in 0:max.age){
          if(j <=tM){cw <- k/(1-exp(-k*(j-t0)))}
          else{cw <- k/(a0+a1*(j-tM)+a2*(j-tM)^2)}
          CW <- cbind(CW,cw)
        }

        M.table[,"Chen.Want"] <- as.vector(CW)
      }

      if(ncol(M.table) == 1){
        M.select <- sample(x=ncol(M.table),size=1)
        Mortality <- M.table[,M.select]
      } else {
        M.table <- M.table[,colSums(is.na(M.table))<nrow(M.table)]
        M.select <- sample(x=ncol(M.table),size=1)
        Mortality <- M.table[,M.select]
      }
    }

    # Estimate Survivorship from M and F estimates. If an AAFC or AALC is set then F is only applied to ages
    # older than the AAFC or younger than AALC. If both are zero then F is not included.
    # An error (specified earlier) is thrown if both are non-zero
    survivorship <- exp(-(Mortality))



    sx <- survivorship

    # sx <- append(sx,0,after = Inf)

    # Calculate fecundity at age based on knife edge or continous dists.

    if(data$maturity$model.type %in% c("normal", "uniform")){
      mx <- rep(0, length(Age))
      mx[age_at_first_repro:max.age+1] <- n.pups.year
    } else{
      mx <- Fecundity_ogive * n.pups.year
    }

    # line up vector so sx and mx are offset by 1 so that eventually fx = mx*sx+1...
    a <- mx[-1]

    # add zero to the end so both vectors are the same length again.
    MX <- append(a,a[max.age],after=length(a))

    # Create top row of matrix
    fx <- MX*sx

    # create a matrix entirely of zeros
    Sm <- matrix(rep(0),nrow=length(survivorship),ncol=length(survivorship))

    #insert survival along the diagonal
    diag(Sm) <- survivorship

    # join fecundity vector to survival matrix without row names
    A <- cbind(rbind(fx,Sm,deparse.level = 0),0)


    ############ Matrix projection ##########
    if(any(is.na(eigen.analysis(A)$repro.value))) next # skip this sim if NA's occur
    if(any(is.na(eigen.analysis(A)$stable.stage))) next
    if(any(eigen.analysis(A)$repro.value > 100)) next # skip this sim if extreme values occur
    if(any(eigen.analysis(A)$stable.stage > 10)) next


    Stable.age.results[,i] <- eigen.analysis(A)$stable.stage
    Repro.val.results[,i] <- eigen.analysis(A)$repro.value

  }

  Stable.age.results<-data.frame(Type = "Stable Age Distribution",
                                 Age = 0:(nrow(Stable.age.results)-1),
                                 mean = apply(Stable.age.results,1,mean, na.rm = TRUE),
                                 low = apply(Stable.age.results, 1, quantile, probs = c(.025), na.rm = TRUE),
                                 upp =apply(Stable.age.results, 1, quantile, probs = c(.975), na.rm = TRUE))

  Repro.val.results<-data.frame(Type = "Reproductive Value",
                                Age = 0:(nrow(Repro.val.results)-1),
                                mean = apply(Repro.val.results,1,mean, na.rm = TRUE),
                                low = apply(Repro.val.results, 1, quantile, probs = c(.025), na.rm = TRUE),
                                upp =apply(Repro.val.results, 1, quantile, probs = c(.975), na.rm = TRUE))
  results <- rbind(Stable.age.results, Repro.val.results)

  close(pb)
  return(results)
}

#' Estimate the range of mortality estimates produced by the different estimators
#' @description This function performs a Monte Carlo simulation analysis which determines
#'     the distributions of natural mortality (M) produced by each estimator.
#' @param n The number of simulations to be run
#' @param data A multi-level list of the class `Demography.inputs` produced from the
#'     `create_data_input` function and then manually completed.
#' @param M.estimators Any specific natural mortality estimators to be included in the analysis.
#'      Must be a single estimator or a vector of estimators. These can include: "Pet.Wro",
#'     "Jensen.mat","Chen.Yuan","Then_hoenig","Then_pauly", "Jensen.mat","Charnov" or
#'     "Chen.Want". If none are specified then all applicable estimators could be chosen.
#'
#' @return A list which includes a data.frame of the distributions age invariant estimators
#'     summarised as the mean and 95\% quantiles and dataframes of any age dependent estimators
#'     with the mean and 95\% quantiles for each age class
#' @export
#' @examples
#' # load Silky shark data produced by create_data_input()
#' # Type `?create_data_input()` for details
#' data("Silky_data")
#'
#' # Run function to get natural mortality distributions from  Monte Carlo
#' # simulations for all available natural mortality estimators.
#' # Set n = at least 1000 for full analysis but use n = 100 for testing
#'
#' Estimate_mortality_dists(n = 100, Silky_data)
#' @import dplyr tidyr readr

Estimate_mortality_dists <- function(n = 1000, data, M.estimators = NULL){

  ## Error messages ------------------

  if(!inherits(data,"Demography.inputs")) stop("Input data is not the correct based on the demography data list")

  if(!data$growth$model.type %in%  c("Von Bertalanffy", "logistic", "Gompertz"))
    stop("Growth model type not correctly specified")

  possible_estimators <- c("Pet.Wro","Jensen.mat","Chen.Yuan","Then_hoenig","Then_pauly", "Jensen","Charnov","Chen.Want")
  if(any(!M.estimators %in% possible_estimators)) stop("One or more M estimators are not correctly specified")
  if(is.null(M.estimators)){
    M.estimators <- possible_estimators
  }

  # remove all mortality estimates that require VB params if VB is not being used
  if(data$growth$model.type != "Von Bertalanffy") M.estimators <- M.estimators[which(!M.estimators %in% c("Chen.Want", "Charnov", "Jensen", "Then_pauly"))]
  if(length(M.estimators) == 0) stop("Specified mortality estimators could not be applied")

  #--------------------------

  # Empty containers to fill
  M.table <- matrix(nrow = n, ncol = length(M.estimators))
  colnames(M.table)<- M.estimators
  P_W_table <-  matrix(nrow = 100, ncol = n)
  Charnov.table <- matrix(nrow = 100, ncol = n)
  Chen.Want.table <- matrix(nrow = 100, ncol = n)

  pb <- txtProgressBar(title="Demography progress", label="0% done", min=0, max=100, initial=0, style = 3)
  for(i in 1:n){
    info <- sprintf("%d%% done", round((i/n)*100))
    setTxtProgressBar(pb, (i/n)*100, label=info)
    # If L0 is being used
    if(!any(names(data$growth$pars) %in% "t0")){

      Linf.se <- data$growth$se$Linf.se
      k.se <- data$growth$se$k.se
      L0.se <- data$growth$se$L0.se

      # if the parameter correlation is provided used a MVN dist
      if(all(!is.na(data$growth$corr.matrix))){
        e.cov <- data$growth$corr.matrix*as.matrix(c(Linf.se, k.se, L0.se))%*%t(as.matrix(c(Linf.se, k.se, L0.se)))
        growth_pars <- MASS::mvrnorm(1, mu = c(data$growth$pars$Linf, data$growth$pars$k, data$growth$pars$L0), Sigma = e.cov)
        Linf <-as.numeric(growth_pars["Linf"])
        k <- as.numeric(growth_pars["k"])
        L0 <- as.numeric(growth_pars["L0"])
      } else if(all(!is.na(Linf.se),!is.na(k.se),!is.na(L0.se))){
        #If MVN can't be used, then vary parameters based on normal dists
        Linf <- rnorm(1, data$growth$pars$Linf, Linf.se)
        k <- rnorm(1, data$growth$pars$k, k.se)
        L0 <- rnorm(1, data$growth$pars$L0, L0.se)
      } else{
        # if no errors available then just use the parameters
        warning("Growth is deterministic as no growth parameter errors were provided")
        Linf <- data$growth$pars$Linf
        k <- data$growth$pars$k
        L0 <- data$growth$pars$L0
      }
    } else {
      # If t0 is being used the same approach is applied but without L0
      Linf.se <- data$growth$se$Linf.se
      k.se <- data$growth$se$k.se
      t0.se <- data$growth$se$t0.se

      if(all(!is.na(data$growth$corr.matrix))){
        e.cov <- data$growth$corr.matrix*as.matrix(c(Linf.se, k.se, t0.se))%*%t(as.matrix(c(Linf.se, k.se, t0.se)))
        growth_pars <- MASS::mvrnorm(1, mu = c(data$growth$pars$Linf, data$growth$pars$k, data$growth$pars$t0), Sigma = e.cov)
        Linf <-as.numeric(growth_pars["Linf"])
        k <- as.numeric(growth_pars["k"])
        t0 <- as.numeric(growth_pars["t0"])
      } else if(all(!is.na(Linf.se),!is.na(k.se),!is.na(t0.se))){
        Linf <- rnorm(1, data$growth$pars$Linf, Linf.se)
        k <- rnorm(1, data$growth$pars$k, k.se)
        t0 <- rnorm(1, data$growth$pars$t0, t0.se)
      } else{
        warning("Growth is deterministic as no growth parameter errors were provided")
        Linf <- data$growth$pars$Linf
        k <- data$growth$pars$k
        t0 <- data$growth$pars$t0
      }
    }

    # Set max age lower bound
    if(!is.na(data$max.age$min)){
      max.age.lower.bound <- data$max.age$min
    } else {
      stop("No min age specified")
    }

    ### Age at 99% of linf (or 95% of Von Bertalanffy used). Used as upper max age
    if(data$`growth`$model.type == "logistic"){
      tmax <- log(((Linf*0.99)*(Linf-L0))/(L0*(Linf - (Linf*0.99))))/k
    } else if(data$`growth`$model.type == "Von Bertalanffy"){
      if(!any(names(data$growth$pars) %in% "t0")){
        tmax <- log(1-((Linf*0.95)-L0)/(Linf-L0))/-k
      } else {
        tmax <- (log(1-((Linf*0.95)/Linf))/-k)+t0
      }
    } else if(data$`growth`$model.type == "Gompertz"){
      tmax <- log(log(Linf/L0)/(log(Linf/L0)-log((Linf*0.99)/L0)))/k
    } else { stop("model.type incorrectely specific for growth")}

    # Set max age upper bound if pre-specified
    if(!is.na(data$max.age$max)){
      max.age.upper.bound <- data$max.age$max
    } else {
      max.age.upper.bound <- 1000 # To make sure its ignored if unspecified
    }

    # Set max age based on available parameters
    if(max.age.lower.bound > tmax | is.na(tmax)){
      max.age <- max.age.lower.bound
    }else if( max.age.upper.bound < tmax){
      max.age <- round(runif(1,max.age.lower.bound,max.age.upper.bound))
    } else{
      max.age <- round(runif(1,max.age.lower.bound,tmax))
    }
    Age <- 0:max.age

    # Length-at-age estimates
    if(data$`growth`$model.type == "logistic"){
      Lt <- Linf*L0*exp(k*Age)/(Linf+L0*(exp(k*Age)-1))
    } else if(data$`growth`$model.type == "Von Bertalanffy"){
      if(!any(names(data$growth$pars) %in% "t0")){
        Lt <- Linf-(Linf-L0)*(exp(-k*Age))
      } else {
        Lt <- Linf*(1-exp(-k*(Age-t0)))
      }
    } else if(data$`growth`$model.type == "Gompertz"){
      Lt <- L0*(exp(log(Linf/L0)*(1-exp(-k*Age))))
    } else {
      stop("model.type incorrectely specific for growth")
    }

    #----End growth calcs------------------

    if(!is.na(data$gest.period)){
      gest.period <- data$gest.period
    } else {
      stop("No gestation period specified")
    }

    if(!is.na(data$repro.cycle)){
      repro.cycle <- data$repro.cycle
    } else {
      stop("No reproductive cycle (periodicity) specified")
    }

    #Maturity ogive
    if(data$maturity$model.type == "logistic - int/slope"){
      mat.cov <- data$maturity$corr.matrix*as.matrix(c(data$maturity$se$intercept.se, data$maturity$se$slope.se))%*%
        t(as.matrix(c(data$maturity$se$intercept.se, data$maturity$se$slope.se)))
      Mat_pars <- MASS::mvrnorm(1, mu = c(data$maturity$pars$intercept, data$maturity$pars$slope), Sigma = mat.cov)
      Mat_intercept <- as.numeric(Mat_pars[1])
      Mat_slope <-as.numeric(Mat_pars[2])
      Maturity_ogive <- round(1/(1+exp(-(Age*Mat_slope+Mat_intercept))),2)
      age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
      age_at_first_repro <- round(gest.period+age_at_maturity)
      Fecundity_ogive <- round(1/(1+exp(-((Age-gest.period)*Mat_slope+Mat_intercept))),2)
    } else if(data$maturity$model.type == "logistic - a50/a95"){
      Mat_a50 <- rnorm(1,data$maturity$pars$a50, data$maturity$se$a50.se)
      Mat_a95 <- rnorm(1,data$maturity$pars$a95, data$maturity$se$a95.se)
      Maturity_ogive <- round(1*(1+exp(-log(19)*((Age-Mat_a50)/(Mat_a95-Mat_a50))))^-1,2)
      age_at_maturity <- which.max(Maturity_ogive > 0.5) # assign age-at-maturity based on inflection of ogive
      age_at_first_repro <- round(gest.period+age_at_maturity)
      Fecundity_ogive <- round( 1*(1+exp(-log(19)*(((Age-gest.period)-Mat_a50)/(Mat_a95-Mat_a50))))^-1,2)
    } else if(data$maturity$model.type == "uniform"){
      age_at_maturity <- runif(1,data$maturity$pars$min, data$maturity$pars$max)
      age_at_first_repro <- round(gest.period+age_at_maturity)
    } else if(data$maturity$model.type == "normal"){
      age_at_maturity <- rnorm(1,data$maturity$pars$mean, data$maturity$pars$se)
      age_at_first_repro <- round(gest.period+age_at_maturity)
    } else {
      stop("Maturity model type is not correctly specified.")
    }


    ### Mortality calculations -------------------
    if("Pet.Wro" %in% M.estimators){

      # If the length weight model type is the same as the length type then no conversion
      # is needed
      if(data$Lt.to.Wt$model.type == data$Lt.type){
        Length.for.wt.conv <- Lt
      } else{

        # Otherwise perform, conversion as necessary. Currently, this converts TL to
        # PCL or FL only. If PCL or FL is used for length in the growth model then typically
        # They are also used in the length weight relationship and don't need converting.
        if(data$convert.TL$model.type == "PCL" & data$Lt.to.Wt$model.type == "PCL"){
          Length.for.wt.conv<-(Lt-data$convert.TL$pars$a)/data$convert.TL$pars$b
        } else if(data$convert.TL$model.type == "FL" & data$Lt.to.Wt$model.type == "FL"){
          Length.for.wt.conv<-(Lt*data$convert.TL$pars$a)-data$convert.TL$pars$b
        } else if(data$Lt.to.Wt$model.type == "TL"){
          Length.for.wt.conv <- Lt
        } else {
          Length.for.wt.conv <- NA # correct length cannot be calculated so return NA to nullify this method
        }
      }
      Wt <- (data$Lt.to.Wt$pars$a)*Length.for.wt.conv^data$Lt.to.Wt$pars$b

      P_W_table[1:(max.age+1),i]<-((Wt^-0.25)*1.92)*0.2

    }

    if("Jensen.mat" %in% M.estimators){
      M.table[i,"Jensen.mat"] <- 1.65/age_at_maturity}
    if("Chen.Yuan" %in% M.estimators){
      M.table[i,"Chen.Yuan"] <-exp(1.46-1.01*log(tmax))}
    if("Then_hoenig" %in% M.estimators){
      M.table[i,"Then_hoenig"] <- 4.899*tmax^-0.916}
    if("Then_pauly" %in% M.estimators &
       data$growth$model.type == "Von Bertalanffy"){
      M.table[i,"Then_pauly"] <- 4.118*k^0.73*Linf^-0.33
    }
    if("Jensen" %in% M.estimators&
       data$growth$model.type == "Von Bertalanffy"){
      M.table[i,"Jensen" ]  <- 1.5*k
    }
    if("Charnov" %in% M.estimators&
       data$growth$model.type == "Von Bertalanffy"){
      Charnov.table[1:(max.age+1),i] <- k*(Lt/Linf)^-1.5
    }
    if("Chen.Want" %in% M.estimators&
       data$growth$model.type == "Von Bertalanffy"){

      if(!any(names(data$growth$pars) %in% "t0")){
        t0<- (1/k)*log((Linf-L0)/Linf)
      }

      # Chen and Wanatabe method
      tM <- -(1/k)*log(1-exp(k*t0))+t0
      a0 <- 1-exp(-k*(tM-t0))
      a1 <- k*exp(-k*(tM-t0))
      a2 <- -0.5*k^2*exp(-k*(tM-t0))
      CW  <- NULL
      for(j in 0:max.age){
        if(j <=tM){cw <- k/(1-exp(-k*(j-t0)))}
        else{cw <- k/(a0+a1*(j-tM)+a2*(j-tM)^2)}
        CW <- cbind(CW,cw)
      }

      Chen.Want.table[1:(max.age+1),i] <- as.vector(CW)
    }

  }
  results_list <- list()
  M.table <- M.table[,colSums(is.na(M.table))<nrow(M.table)]
  results_list[["Age_invariant_mortality"]] <- as.data.frame(M.table) %>%
    gather(Method, Value) %>%
    group_by(Method) %>%
    summarise(AVG = mean(Value,na.rm = TRUE),
              low = quantile(Value, na.rm = TRUE, probs = c(0.025)),
              high = quantile(Value,na.rm = TRUE,  probs = c(0.975)),.groups = "drop_last")


  if(!all(is.na(P_W_table))){
    P_W_table <- P_W_table[rowSums(is.na(P_W_table))<nrow(P_W_table),]
    results_list[["Peterson_wroblewski"]] <-  P_W_table %>%
      t() %>%
      as.data.frame() %>%
      gather(Age, Value) %>%
      mutate(Age = readr::parse_number(Age)-1) %>%
      group_by(Age) %>%
      summarise(AVG = mean(Value,na.rm = TRUE),
                low = quantile(Value, na.rm = TRUE, probs = c(0.025)),
                high = quantile(Value,na.rm = TRUE,  probs = c(0.975)),.groups = "drop_last")
  }

  if(!all(is.na(Chen.Want.table))){
    Chen.Want.table <- Chen.Want.table[rowSums(is.na(Chen.Want.table))<ncol(Chen.Want.table),]

    results_list[["Chen_Watanabe"]]  <- Chen.Want.table %>%
      t() %>%
      as.data.frame() %>%
      gather(Age, Value) %>%
      mutate(Age = readr::parse_number(Age)-1) %>%
      group_by(Age) %>%
      summarise(AVG = mean(Value,na.rm = TRUE),
                low = quantile(Value, na.rm = TRUE, probs = c(0.025)),
                high = quantile(Value,na.rm = TRUE,  probs = c(0.975)),.groups = "drop_last")
  }
  if(!all(is.na(Charnov.table))){
    Charnov.table <- Charnov.table[rowSums(is.na(Charnov.table))<ncol(Charnov.table),]
    results_list[["Charnov"]]  <- Charnov.table %>%
      t() %>%
      as.data.frame() %>%
      gather(Age, Value) %>%
      mutate(Age = readr::parse_number(Age)-1) %>%
      group_by(Age) %>%
      summarise(AVG = mean(Value,na.rm = TRUE),
                low = quantile(Value, na.rm = TRUE, probs = c(0.025)),
                high = quantile(Value,na.rm = TRUE,  probs = c(0.975)),.groups = "drop_last")
  }
  close(pb)
  return(results_list)
}

#' Monte Carlo simulations of Leslie matrix models
#' @description This is a wrapper function for `Calculate_demography` which runs this function
#'     the specified number of times.
#' @param n The number of specified Monte Carlo simulations to run
#' @param data A multi-level list of the class `Demography.inputs` produced from the
#'     `create_data_input` function and then manually completed.
#' @param AAFC Age-at-first-capture which can be specified by the user.Must be an integer
#'     age which can include zero to indicate the availability of the population to capture from birth.
#' @param F. The instantaneous rate of fishing mortality 'F'. This will be applied to all
#'     ages available to capture as defined by either the AALC or AAFC arguments.
#' @param AALC Age-at-last-capture which can be specified by the user. Must be an integer
#'     age which can include zero.
#' @param M.estimators Any specific natural mortality estimators to be included in the analysis.
#'     Only one will be used in each run which is randomly selected. Must be a single estimator
#'     or a vector of estimators. These can include: "Pet.Wro","Jensen.mat","Chen.Yuan",
#'     "Then_hoenig","Then_pauly", "Jensen.mat","Charnov" or "Chen.Want". If none are specified
#'     then all applicable estimators could be chosen.
#' @param Verbatim Print summary results to screen if TRUE. When FALSE, the progress bar is also disabled.
#' @return A list with two data.frames. The first is the summary of the Monte Carlo simulations
#'     for all parameters calculated by the `Calculate_demography` function with mean and
#'     95\% quantiles. The second is all of the results for each parameter from individual
#'     simulations so that their distributions can be interrogated further.
#' @examples
#' # load Silky shark data produced by create_data_input()
#' # Type `?create_data_input()` for details
#' data("Silky_data")
#'
#' # Run function to get conduct Monte Carlo Simulations using
#' # `Calculate_demography()`  for all available natural mortality estimators.
#' # Set n = at least 1000 for full analysis but use n = 100 for testing
#'
#' Simulate_demography(n = 100, Silky_data)
#' @export
Simulate_demography <- function(n, data, AALC = NULL, AAFC = NULL, F. = 0, M.estimators = NULL, Verbatim = TRUE){

  if(all(!is.na(data$Mortality$mean), !is.na(data$Mortality$se))){
    warning("M.estimators ignored as data contains a pre-specified mortality distribution")
  }

  lambda.results <- NULL
  R0.results <- NULL
  G.results <- NULL
  elast.fecund <- NULL
  elast.juv.survival <- NULL
  elast.adult.survival <- NULL
  juv.ratio <- NULL
  adult.ratio <- NULL
  M.estimator <- NULL

  if(Verbatim == TRUE){
    pb <- txtProgressBar(title="Demography progress", label="0% done", min=0, max=100, initial=0, style = 3)
  }

  for(i in 1:n){

    if(Verbatim == TRUE){
      info <- sprintf("%d%% done", round((i/n)*100))
      setTxtProgressBar(pb, (i/n)*100, label=info)
    }

     tryCatch(
      {
        tmp <- Calculate_demography(data, AAFC, F., AALC, M.estimators);
        lambda.results[i] <- tmp$`lambda`;
        R0.results[i] <- tmp$`R0`;
        G.results[i] <- tmp$`G`;
        elast.fecund[i] <-  tmp$elast.fecund;
        elast.juv.survival[i] <-  tmp$elast.juv.survival;
        elast.adult.survival[i] <- tmp$elast.adult.survival;
        juv.ratio[i] <- tmp$juv.ratio;
        adult.ratio[i] <- tmp$adult.ratio;
        M.estimator[i] <- tmp$M.estimator;
      },
      error=function(e){
        lambda.results[i] <- NA;
        R0.results[i] <- NA;
        G.results[i] <-NA;
        elast.fecund[i] <-  NA;
        elast.juv.survival[i] <-  NA;
        elast.adult.survival[i] <- NA;
        juv.ratio[i] <- NA;
        adult.ratio[i] <- NA;
        M.estimator[i] <- NA;
      },
      silent = TRUE
    )
    # tmp <- Calculate_demography(data, AAFC, F., AALC, M.estimators)
    # lambda.results[i] <- tmp$`lambda`
    # R0.results[i] <- tmp$`R0`
    # G.results[i] <- tmp$`G`
    # elast.fecund[i] <-  tmp$elast.fecund
    # elast.juv.survival[i] <-  tmp$elast.juv.survival
    # elast.adult.survival[i] <- tmp$elast.adult.survival
    # juv.ratio[i] <- tmp$juv.ratio
    # adult.ratio[i] <- tmp$adult.ratio
    # M.estimator[i] <- tmp$M.estimator
  }

  #=========================================================================
  # Analyse Monte Carlo outputs

  lambda.summary <- c(mean(lambda.results, na.rm = TRUE),
                      quantile(lambda.results, c(.025, .975),na.rm = TRUE))

  G.summary <- c(mean(G.results, na.rm = TRUE),
                 quantile(G.results, c(.025, .975),na.rm = TRUE))

  R0.summary <- c(mean(R0.results, na.rm = TRUE),
                  quantile(R0.results, c(.025, .975),na.rm = TRUE))

  elast.fecund.summary <- c(mean(elast.fecund, na.rm = TRUE),
                            quantile(elast.fecund, c(.025, .975),na.rm = TRUE))
  elast.juv.survival.summary <- c(mean(elast.juv.survival, na.rm = TRUE),
                                  quantile(elast.juv.survival, c(.025, .975),na.rm = TRUE))
  elast.adult.survival.summary <- c(mean(elast.adult.survival, na.rm = TRUE),
                                    quantile(elast.adult.survival, c(.025, .975),na.rm = TRUE))
  juv.ratio.summary <- c(mean(juv.ratio, na.rm = TRUE),
                         quantile(juv.ratio, c(.025, .975),na.rm = TRUE))
  adult.ratio.summary <- c(mean(adult.ratio, na.rm = TRUE),
                           quantile(adult.ratio, c(.025, .975),na.rm = TRUE))

  simulations <- bind_cols(Lambda = lambda.results,
                           R0 = R0.results,
                           G = G.results,
                           elast.fecund = elast.fecund,
                           elast.juv.survival = elast.juv.survival,
                           elast.adult.survival = elast.adult.survival,
                           juv.ratio=juv.ratio,
                           adult.ratio = adult.ratio,
                           M.estimator = M.estimator)

  Results <- rbind(lambda.summary,G.summary,R0.summary,
                   elast.fecund.summary,elast.juv.survival.summary,elast.adult.survival.summary,
                   juv.ratio.summary, adult.ratio.summary)
  colnames(Results) <- c("Mean",".025",".975")
  rownames(Results) <- c("Lambda",
                         "Generation Time",
                         "Net Repro Rate",
                         "Fertility elasisticy",
                         "Juvenile elasticity ",
                         "Adult elasticity ",
                         "Juvenile elasticity ratio",
                         "Adult elasticity ratio")


  if(Verbatim == TRUE){
    close(pb)
    print(Results)
  }

  return(list(MonteCarlo_summary = Results, simulations = simulations))
}

#' Monte Carlo simulations of Leslie matrix models under varying levels of F and AAFC
#' @description This is a wrapper function for `Calculate_demography` which runs this function
#'     the specified number of times using a range of F and Age-at-first-capture (AAFC)
#'     values
#' @param n The number of simulations to be run. 1000 is recomended but smaller numbers
#'     should be run when testing to avoid long run times.
#' @param .data A multi-level list of the class `Demography.inputs` produced from the
#'     `create_data_input` function and then manually completed.
#' @param M.estimators Any specific natural mortality estimators to be included in the analysis.
#'     Only one will be used in each run which is randomly selected. Must be a single estimator
#'     or a vector of estimators. These can include: "Pet.Wro","Jensen.mat","Chen.Yuan",
#'     "Then_hoenig","Then_pauly", "Jensen.mat","Charnov" or "Chen.Want". If none are specified
#'     then all applicable estimators could be chosen.
#' @param max.AAFC The maximum Age class to be run in the analyses. This does not need to
#'     be the maximum age for the population and keeping this number reasonable reduces
#'     run-time.
#' @param n_cores The number of cores to be used for parallel processing. It should be 1 core less than the
#'     maximum number available.
#' @return A list with two data.frames. The first provides Fcritical values for each age class.
#'     This is the value of F where the population growth rate is stable (lambda = 1). The
#'     second dataframe is the mean lambda produced for each combination of AAFC and F
#' @export
#' @examples
#' \donttest{
#' # load Silky shark data produced by create_data_input()
#' # Type `?create_data_input()` for details
#' data("Silky_data")
#'
#' # Run function to get conduct an age-at-first-capture (AAFC) analysis using
#' # Monte Carlo Simulations using for all available natural mortality estimators.
#' # Set n = at least 1000 for full analysis but use n = 10 for testing given long run times
#'
#' Simulate_AAFC(n = 10, Silky_data, n_cores = 1)
#' }
#' @import doParallel foreach iterators parallel doFuture
#'

Simulate_AAFC<- function(n = 1000, .data, M.estimators = NULL, max.AAFC = 15, n_cores = 1){
  pb <- txtProgressBar(title="Demography progress", label="0% done", min=0, max=100, initial=0, style = 3)

  AAFC.range <- seq(0,max.AAFC, 1)
  F.range <- seq(0,1,0.01)

  if(n_cores >  parallel::detectCores()-1) {
    n_cores <- 1
    message("Not enough cores available. Reseting to 1 core")
  }

  cl <- parallel::makeCluster(useXDR=FALSE,n_cores)

  doParallel::registerDoParallel(cl)
  doFuture::registerDoFuture()
  mean.lambda.given.F.at.AAFC <- matrix(ncol = 3,nrow = length(F.range) * length(0:max.AAFC))
  colnames(mean.lambda.given.F.at.AAFC) <- c("AAFC", "F.", "mean.lambda")
  I <- 1:length(F.range)
  for(AAFC in 0:max.AAFC){
    setTxtProgressBar(pb, (AAFC/max.AAFC)*100, label=info)

    results <- suppressWarnings(foreach::foreach(F. = F.range, .combine = rbind, .packages = c("popbio","iterators"),
                                        .inorder = TRUE,
                                        .export = c("Calculate_demography", "AAFC", "data", "M.estimators", "n",ls(envir=globalenv()))) %dopar%{
                                          foreach::foreach(icount(n), .packages = c("iterators","popbio"),.combine = rbind, .inorder = FALSE,
                                                  .export = c("Calculate_demography", ls(envir=globalenv()))) %dopar%{

                                                    tryCatch(
                                                      data.frame(AAFC = AAFC, F. = F. ,
                                                                 lambda = Calculate_demography(data = .data, AALC = NULL, AAFC = AAFC,
                                                                                               F. = F., M.estimators = M.estimators)$`lambda`),
                                                      error=function(e){data.frame(AAFC = AAFC, F. = F. , lambda = NA)},
                                                      silent = TRUE
                                                    )

                                                  }})

    lambda_results <- results %>%
      filter(lambda < 2) %>% # in case irregular values occur and throw off distributions
      group_by(AAFC, F.) %>%
      summarise(mean.lambda = mean(lambda, na.rm = TRUE),.groups = "drop_last") %>%
      as.data.frame()

    mean.lambda.given.F.at.AAFC[I,] <- as.matrix(lambda_results)
    I <- I + length(F.range)
  }
  suppressWarnings(parallel::stopCluster(cl))
  AAFC_F_critical <-  as.data.frame(mean.lambda.given.F.at.AAFC) %>%
    filter( mean.lambda <= 1) %>%
    group_by(AAFC) %>% slice(which.max(mean.lambda)) %>%
    as.data.frame()

  close(pb)
  AAFC.Results <- list(AAFC_F_critical, mean.lambda.given.F.at.AAFC)
  return(AAFC.Results)

}

#' Monte Carlo simulations of Leslie matrix models under varying levels of F and AALC
#' @description This is a wrapper function for `Calculate_demography` which runs this function
#'     the specified number of times using a range of F and Age-at-last-capture (AALC)
#'     values
#' @param n The number of simulations to be run. 1000 is recomended but smaller numbers
#'     should be run when testing to avoid long run times.
#' @param data A multi-level list of the class `Demography.inputs` produced from the
#'     `create_data_input` function and then manually completed.
#' @param M.estimators Any specific natural mortality estimators to be included in the analysis.
#'     Only one will be used in each run which is randomly selected. Must be a single estimator
#'     or a vector of estimators. These can include: "Pet.Wro","Jensen.mat","Chen.Yuan",
#'     "Then_hoenig","Then_pauly", "Jensen.mat","Charnov" or "Chen.Want". If none are specified
#'     then all applicable estimators could be chosen.
#' @param min.AALC The last Age class to be run in the analyses. This does not need to
#'     be the maximum age for the population and keeping this number reasonable reduces
#'     run-time.
#' @param n_cores The number of cores to be used for parallel processing. It should be 1 core less than the
#'     maximum number available.
#' @return A list with two data.frames. The first provides Fcritical values for each age class.
#'     This is the value of F where the population growth rate is stable (lambda = 1). The
#'     second dataframe is the mean lambda produced for each combination of AALC and F
#' @export
#' @examples
#' \donttest{
#' # load Silky shark data produced by create_data_input()
#' # Type `?create_data_input()` for details
#' data("Silky_data")
#'
#' # Run function to get conduct an age-at-last-capture (AALC) analysis using
#' # Monte Carlo Simulations using for all available natural mortality estimators.
#' # Set n = at least 1000 for full analysis but use n = 10 for testing given long run times
#'
#' Simulate_AALC(n = 10, Silky_data, n_cores = 1)
#' }
#' @import doParallel foreach iterators parallel doFuture
#'
Simulate_AALC <- function(n = 1000, data, M.estimators = NULL,  min.AALC = 15, n_cores = 1){
  pb <- txtProgressBar(title="Demography progress", label="0% done", min=0, max=100, initial=0, style = 3)

  AALC.range <- seq(0,min.AALC, 1)
  F.range <- seq(0,1,0.01)

  if(n_cores >  parallel::detectCores()-1) {
    n_cores <- 1
    message("Not enough cores available. Reseting to 1 core")
  }

  cl <- parallel::makeCluster(useXDR=FALSE,n_cores)
  doParallel::registerDoParallel(cl)
  doFuture::registerDoFuture()

  mean.lambda.given.F.at.AALC <- matrix(ncol = 3,nrow = length(F.range) * length(0:min.AALC))
  colnames(mean.lambda.given.F.at.AALC) <- c("AALC", "F.", "mean.lambda")
  I <- 1:length(F.range)
  for(AALC in 0:min.AALC){
    setTxtProgressBar(pb, (AALC/min.AALC)*100, label=info)


    results <- suppressWarnings(
      foreach::foreach(F. = F.range, .combine = rbind, .packages = c("popbio","iterators"),
              .inorder = TRUE,
              .export = c("Calculate_demography", "AALC", "data", "M.estimators", "n",ls(envir=globalenv()))) %dopar%{
                foreach::foreach(icount(n), .packages = c("iterators","popbio"),.combine = rbind, .inorder = FALSE,
                        .export = c("Calculate_demography", "AALC", "data", "M.estimators", "n",ls(envir=globalenv()))) %dopar%{
                          tryCatch(
                            data.frame(AALC = AALC, F. = F. ,
                                       lambda = Calculate_demography(data = data, AALC = AALC,
                                                                     F. = F., M.estimators = M.estimators)$`lambda`),
                            error=function(e){data.frame(AALC = AALC, F. = F. , lambda = NA)},
                            silent = TRUE
                          )
                        }})

    lambda_results <- results %>%
      filter(lambda < 2) %>% # in case irregular values occur and throw off distributions
      group_by(AALC, F.) %>%
      summarise(mean.lambda = mean(lambda, na.rm = TRUE),.groups = "drop_last") %>%
      as.data.frame()

    mean.lambda.given.F.at.AALC[I,] <- as.matrix(lambda_results)
    I <- I + length(F.range)

  }
  suppressWarnings(parallel::stopCluster(cl))


  AALC_F_critical <-  as.data.frame(mean.lambda.given.F.at.AALC) %>%
    filter( mean.lambda <= 1) %>%
    group_by(AALC) %>% slice(which.max(mean.lambda)) %>%
    as.data.frame()

  if(all(AALC_F_critical$AALC > 0)){
    AALC_F_critical <- AALC_F_critical %>% bind_rows(data.frame(AALC = 0, F. = max(F.range), mean.lambda = 1)) %>%
      arrange(AALC)
  }

  close(pb)
  AALC.Results <- list(AALC_F_critical, mean.lambda.given.F.at.AALC)
  return(AALC.Results)

}

#' Estimate F critical through simulations
#' @description This is a wrapper function for `Calculate_demography` which runs this function
#'     the specified number of times using a range of F values accross the enture age range of
#'     the population. This determines the rate of population increase at each increment of F.
#' @param n The number of simulations to be run. 1000 is recomended but smaller numbers
#'     should be run when testing to avoid long run times.
#' @param .data A multi-level list of the class `Demography.inputs` produced from the
#'     `create_data_input` function and then manually completed.
#' @param M.estimators Any specific natural mortality estimators to be included in the analysis.
#'     Only one will be used in each run which is randomly selected. Must be a single estimator
#'     or a vector of estimators. These can include: "Pet.Wro","Jensen.mat","Chen.Yuan",
#'     "Then_hoenig","Then_pauly", "Jensen.mat","Charnov" or "Chen.Want". If none are specified
#'     then all applicable estimators could be chosen.
#' @param max.F The maximum value of F for simulations
#' @param n_cores The number of cores to be used for parallel processing. It should be 1 core less than the
#'     maximum number available.
#' @return A list with two data.frames. The first provides the mean F critical with 95\%
#'     confidence intervals. The second dataframe provides the rate of increase for each
#'     increment of F.
#' @export
#' @examples
#' # load Silky shark data produced by create_data_input()
#' # Type `?create_data_input()` for details
#' data("Silky_data")
#'
#' # Run function to get conduct an F critical analysis using
#' # Monte Carlo Simulations using for all available natural mortality estimators.
#' # Set n = at least 1000 for full analysis but use n = 10 for testing given long run times
#'
#' Simulate_F_critical(n = 10, Silky_data, n_cores = 1)
#' @import doParallel foreach iterators parallel doFuture

Simulate_F_critical <- function(n = 1000, .data, M.estimators = NULL, max.F = 0.3, n_cores = 1){

  F.range <- seq(0,max.F,0.01)

  if(n_cores >  parallel::detectCores()-1) {
    n_cores <- 1
    message("Not enough cores available. Reseting to 1 core")
  }

  cl <- parallel::makeCluster(useXDR=FALSE,n_cores)
  doParallel::registerDoParallel(cl)
  doFuture::registerDoFuture()

  results <- suppressWarnings(foreach::foreach(F. = F.range, .combine = rbind, .packages = c("popbio","iterators"),
                                      .inorder = TRUE,
                                      .export = c("Calculate_demography", "AAFC", "data", "M.estimators", "n",ls(envir=globalenv()))) %dopar%{

                                        foreach::foreach(icount(n), .packages = c("iterators","popbio"),.combine = rbind, .inorder = FALSE,
                                                .export = c("Calculate_demography", ls(envir=globalenv()))) %dopar%{
                                                  data.frame(F. = F. ,
                                                             lambda = Calculate_demography(data = .data, AALC = NULL, AAFC = 0,
                                                                                           F. = F., M.estimators = M.estimators)$`lambda`)
                                                }})


  suppressWarnings(parallel::stopCluster(cl))

  lambda_results <- group_by(results, F.) %>% summarise(AVG = mean(lambda, na.rm = TRUE), lwr = quantile(lambda,.025, na.rm = TRUE),
                                                        upr = quantile(lambda,.975, na.rm = TRUE),.groups = "drop_last") %>%
    as.data.frame()


  F_critical <- lambda_results %>% filter(AVG >= 1) %>%
    summarise(F. = .[which(.$AVG == min(.$AVG)),"F."],
              AVG = .[which(.$AVG == min(.$AVG)),"AVG"],
              lwr = .[which(.$AVG == min(.$AVG)),"lwr"],
              upr = .[which(.$AVG == min(.$AVG)),"upr"],.groups = "drop_last")

  return(list(F_critical = F_critical, Simulations = lambda_results))

}



#' Simulate_harvest_slots
#' @description This is a wrapper function for `Calculate_demography` which runs this function
#'     using different values of AAFC, AALC and F to simulate Havest Slot options. The Fcritical is retuned
#'     for each simulation to show the max level of F needed to sustain a stable population.
#' @param n number of iterations for each combination of Age.mid.point, HS.width and F.
#' @param data A multi-level list of the class `Demography.inputs` produced from the
#'     `create_data_input` function and then manually completed.
#' @param M.estimators Any specific natural mortality estimators to be included in the analysis.
#'     Only one will be used in each run which is randomly selected. Must be a single estimator
#'     or a vector of estimators. These can include: "Pet.Wro","Jensen.mat","Chen.Yuan",
#'     "Then_hoenig","Then_pauly", "Jensen.mat","Charnov" or "Chen.Want". If none are specified
#'     then all applicable estimators could be chosen.
#' @param Age.mid.point A vector of ages to be used in the simulation. Each age is used as
#'     a mid point and will have HS.Width subtracted and added to it to determine AAFC and AALC,
#'     respectively.
#' @param HS.width A vector of widths for the Harvest slots in years. Widths are subtracted and added to mid points to determine
#'     the AAFC and AALC in each sim.
#' @param max.F The maximum value of F for simulations
#'
#' @return A data.frame with three columns: MinAge, MaxAge and 'F.'. These represent the age at the
#'     start of a harvest slot, the age at the end of the harvest slot and the F for that harvest slot.
#' @export
#' @examples
#' \donttest{
#' # load Silky shark data produced by create_data_input()
#' # Type `?create_data_input()` for details
#' data("Silky_data")
#'
#' # Run function to get conduct an F critical analysis for different harvest slots using
#' # Monte Carlo Simulations using for all available natural mortality estimators.
#' # Set n = at least 1000 for full analysis but use n = 10 for testing given long run times
#'
#' Simulate_harvest_slots(n = 10, Silky_data,Age.mid.point = 0:28, HS.width = 0:8)
#' }
#' @import interp

Simulate_harvest_slots <- function(n, data,  M.estimators = NULL,
                                   Age.mid.point = NULL, HS.width = NULL, max.F = 1){
  if(is.null(Age.mid.point)) stop("Age.mid.point must be a range")
  if(length(Age.mid.point)==1) stop("Age.mid.point must be a range")
  if(is.null(HS.width)) stop("Age.mid.point must be a range")
  if(length(HS.width)==1) stop("HS.width must be a range")
  if(length(max.F)!=1) stop("max.F must be a single value")

  Results <- expand.grid(Age.mid.point = Age.mid.point, HS.width = HS.width, F. = seq(0,max.F,.05),
                         Lambda = NA,  low = NA, upp = NA)

  max.age <- max(Age.mid.point,na.rm = TRUE)

  pb <- txtProgressBar(title="Havest Slot sim progress", label="0% done", min=0, max=100, initial=0, style = 3)
  for(Mid in unique(Results$Age.mid.point)){
    info <- sprintf("%d%% done", round((Mid/max(Age.mid.point))*100))
    setTxtProgressBar(pb, (Mid/max(Age.mid.point))*100, label=info)
    for(Width in unique(Results$HS.width)){
      for(F. in unique(Results$F.)){
        AALC <- ifelse(Mid + Width <= max.age, Mid + Width, max.age)
        AAFC <- ifelse(Mid - Width > 0, Mid - Width, 0)

        if(AAFC > AALC) {
          next
        }
        tmp <- Simulate_demography(n, data = data, AALC = AALC, AAFC = AAFC, F. = F., M.estimators = M.estimators,Verbatim = FALSE)
        Results[which(Results$Age.mid.point == Mid & Results$HS.width == Width & Results$F. == F.),c("Lambda", "low", "upp")] <- tmp[[1]]["Lambda",c("Mean",".025" ,".975")]
      }
    }
  }

  Results <- Results %>%
    mutate(Lambda = round(Lambda,1),
           MinAge = as.integer(ifelse(Age.mid.point - HS.width >=0, Age.mid.point-HS.width,0 )),
           MaxAge = as.integer(ifelse(Age.mid.point + HS.width < max.age, Age.mid.point+HS.width,max.age )),
           Range = (MaxAge - MinAge)+1
    )%>%
    filter( Lambda == 1) %>%
    group_by( Age.mid.point, Range) %>%
    slice(which.max(F.))


   tryCatch(
    {
      interp.data <- interp::interp(Results$MinAge, Results$MaxAge, Results$F., nx = 300, ny = 300, duplicate = "mean");
      final_results <- interp.data %>% interp::interp2xyz() %>% as.data.frame() %>% na.omit;
      names(final_results) <- c("MinAge", "MaxAge", "F.");
    },
    error=function(e){
      warning("Interpolation failed. Un-interpolated results returned instead");
      final_results <- Results;

    },
    silent = TRUE
  )


  close(pb)
  return(final_results)
}

