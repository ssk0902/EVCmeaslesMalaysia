#' two-variate joint Bernoulli distribution with given correlation 
#' @param x probability of variable X = 1
#' @param y probability of variable Y = 1
#' @param r correlation coefficient of X and Y
#' @param digits the number of decimal places
joint_Bernoulli <- function(x, y, r, digits = 6L){
  # assuming X and Y are two variables from Bernoulli distribution taking values 0: not vaccinated, and 1: vaccinated
  # E(X) = x, D(X) = x(1-x)
  # E(Y) = y, D(Y) = y(1-y)
  # denote correlation parameter r
  # this function calculates the joint distribution of X and Y, i.e. F(X, Y)
  stopifnot(is.numeric(x) & is.numeric(y) & is.numeric(r))
  stopifnot(is.integer(digits))
  stopifnot(x>=0 & x<=1 & y>=0 & y<=1 & r>=-1 & r <=1) # constrain values of input parameters
  
  m <- x*y
  n <- sqrt(x*(1-x)*y*(1-y))
  l <- min(x, y)
  
  #1. P(X=1, Y=1)
  p11 <- r * n +m
  if (p11 < 0){
    p11 <- 0
  }
  if (p11 > l){
    p11 <- l
  }
  if (x + y - p11 > 1){
    p11 <- x + y -1
  }
  r_adjusted <- (p11 - m)/n ## adjusted correlation coefficient if out of bound
  
  #2. P(X = 0, Y = 1)
  p01 <- y - p11
  
  #3. P(X = 1, Y = 0)
  p10 <- x - p11
  
  #4. P(X = 0, Y = 0)
  p00 = 1 - p11 - p01 - p10
  
  list(p11 = round(p11, digits),
       
       p00 = round(p00, digits), 
       
       p01 = round(p01, digits), 
       
       p10 = round(p10, digits),
       
       px = x,
       
       py = y,
       
       cor = r,
       
       cor_adjusted = round(r_adjusted, digits)
  )
  
}

#' @param vax_history vaccination history
#' @param para_eff vaccine efficacy by age and dose
#' @param para_corr correlation between vaccinated and vaccination
#' @param cohort_size cohort size
#' @param sur_rates burden related rates, natual_infection = FoI. use this when efficacy def is reduction in infection
#' @param year_min
#' @param year_max
effective_coverage_history <- function(vax_history, para_eff, para_corr=NULL, cohort_size, sur_rates = NULL, year_min = NULL, year_max = NULL){
  ## step 0: check inputs
  disease <- unique(vax_history$disease)
  if(length(disease) > 1){
    stop("Make sure only one disease is supplied.")
  }
  country <- unique(vax_history$country)
  if(length(country) > 1){
    stop("Make sure only one country is supplied.")
  }
  vax_history$activity_coverage <- as.numeric(vax_history$activity_coverage)
  if(any(is.na(vax_history$activity_coverage))){
    print("coercion: NA coverage to zero.")
    vax_history$activity_coverage <- ifelse(is.na(vax_history$activity_coverage), 0, vax_history$activity_coverage)
  }
  if(any(vax_history$activity_coverage < 0)){
    stop("Negative coverage data not accepted.")
  }
  para_eff <- para_eff[para_eff$disease == disease, ]
  if(nrow(para_eff) == 0L | is.null(para_eff)){
    stop("Efficacy parameter, para_eff, cannot be NULL")
  }
  if(is.null(para_corr)){
    message("NULL correlation parameters provided - use random correlation parameters.")
    corr <- runif(3, -1, 1)
    corr <- c(corr[order(-corr)], rep(-1, 7))
    para_corr <- data.frame(disease = disease, dose = 1:10, corr)
  } else {
    para_corr <- para_corr[para_corr$disease == disease, ]
    if(nrow(para_corr) == 0L){
      stop("unrecognised disease in para_corr")
    }
  }
  cohort_size <- cohort_size[cohort_size$country == country, ]
  if(nrow(cohort_size) == 0L | is.null(cohort_size)){
    stop("Cohort size, cohort_size, cannot be NULL")
  }
  if(is.null(sur_rates)){
    efficacy_def <- "all-or-nothing"
  } else {
    efficacy_def <- "reduction-in-infection"
    sur_rates <- sur_rates[sur_rates$country == country, ]
    sur_rates$cohort <- sur_rates$year - sur_rates$age
    ## todo: do some checks
  }
  if(is.null(year_min)){
    year_min <- min(vax_history$year)
  }
  if(is.null(year_max)){
    year_max <- max(vax_history$year)
  }
  if(min(vax_history$year) > year_max){
    year_max <- min(vax_history$year)
    print("coercion: no coverage upto year_max, use intro year")
  }
  vax_history$disease <- NULL
  vax_history$coverage_set <- NULL
  para_eff$disease <- NULL
  para_corr$disease <- NULL
  
  t <- vax_history[vax_history$year %in% year_min:year_max, ]
  t$year <- NULL
  t$country <- NULL
  t$group <- NULL
  cc <- cohort_size[cohort_size$year %in% year_min:year_max, ]
  
  zz <- ifelse(disease %in% c("Measles", "Rubella"), 2, 1) ## 2 for 2 routine doses, 1 for one routine dose
  
  ## expand efficacy parameters
  tmp <- expand.grid(age = 0:100, stringsAsFactors = FALSE)
  tmp <- merge(para_eff, tmp, all.y = TRUE)
  tmp <- tmp[order(tmp$age), ]
  i <- is.na(tmp$single_dose_eff)
  tmp$single_dose_eff[i] <- max(para_eff$single_dose_eff)
  tmp$multi_dose_eff[i] <- max(para_eff$multi_dose_eff)
  para_eff <- tmp
  
  # for each cohort do the following
  d_out <- list(NULL)
  cohorts <- unique(t$cohort)
  
  for(j in seq_along(cohorts)){
    v <- t[t$cohort == cohorts[j], ]
    v <- v[order(-xtfrm(v$activity_type), v$age), ]
    # now vaccination history is ordered by routine -> campaign and age
    # assuming correlation decreases by this order
    v$dose <- seq_along(v$age)
    kk <- v$activity_type == "campaign"
    if(any(kk)){
      min_sia_dose <- min(v$dose[kk])
      v$dose[kk] <- seq_len(sum(kk)) + zz 
    }
    v <- merge(v, para_corr, all.x = TRUE)
    v <- v[order(v$age), ]
    
    # initiating outcome columns
    v$coverage <- v$activity_coverage  # proportion of people vaccinated
    v$single_dose <- v$activity_coverage # single dose proportion
    v$multi_dose <- v$activity_coverage # multiple dose proportion
    v$multi_dose[1] <- 0
    v$zero_dose <- 1- v$coverage # zero dose proportion
    
    ## step 1: estimate the proportion of vaccinated people with zero dose, single dose and multiple doses
    ## this considers correlation between vaccinated and vaccination activity
    if(nrow(v) > 1){
      for(i in seq_len(nrow(v))[-1]){
        ## tmp is the joint distribution between cohort coverage and the next vaccination
        ## tmp2 is the joint distribution between multiple-dose coverage and the next vaccination
        ## they share a common correlation parameter
        tmp <- joint_Bernoulli(v$coverage[i-1], v$activity_coverage[i], v$corr[i])
        tmp2 <- joint_Bernoulli(v$multi_dose[i-1], v$activity_coverage[i], v$corr[i])
        v$zero_dose[i] <- round(tmp$p00, 6)
        v$multi_dose[i] <- round(v$multi_dose[i-1] + tmp$p11 - tmp2$p11, 6)
        v$coverage[i] <- 1 - v$zero_dose[i]
        
      }
    }
    v$single_dose <- v$coverage - v$multi_dose
    
    ## step 2: estimate effective vaccination coverage by cohort and age
    ## this considers vaccine efficacy
    d <- merge(v, para_eff, all.x = TRUE)
    if(efficacy_def == "all-or-nothing"){
      d$effective_coverage <- d$single_dose * d$single_dose_eff + d$multi_dose * d$multi_dose_eff
    } else {
      d <- vimpact:::merge_by_common_cols(d, sur_rates[c("cohort", "age", "natural_infection")], all.x = TRUE)
      d$effective_coverage <- (1 - (1 - d$single_dose_eff) * d$natural_infection) * d$single_dose + 
        (1 - (1 - d$multi_dose_eff) * d$natural_infection) * d$multi_dose
    }
    d$effective_coverage[is.na(d$effective_coverage)] <- 0
    d$tmp <- c(0, d$effective_coverage[-length(d$effective_coverage)])
    d$tmp2 <- c(0, d$coverage[-length(d$coverage)])
    d_out[[j]] <- d
  }
  d <- do.call(rbind, d_out)
  d$dose <- NULL
  d$corr <- NULL
  d$single_dose_eff <- NULL
  d$multi_dose_eff <- NULL
  d$added_effective_coverage <- d$effective_coverage - d$tmp
  d$added_coverage <- d$coverage - d$tmp2
  d$tmp <- NULL
  d$tmp2 <- NULL
  
  # there are multiple vaccinations for the same cohort-age combination, 
  # we only need the row for which multiple activities have been finished
  s <- aggregate(effective_coverage ~ cohort + age, d, max, na.rm = TRUE)
  d <- merge(s, d, all.x = TRUE)
  
  ## step 3: expand to life-course
  #message("need to consider only females for HPV")
  d <- vimpact:::merge_by_common_cols(d, cc[c("cohort", "age", "cohort_size")], all.y = TRUE)
  i <- is.na(d$effective_coverage)
  d$effective_coverage[i] <- 0
  d$added_effective_coverage[i] <- 0
  d$added_coverage[i] <- 0
  d$activity_type[i] <- "none"
  d$activity_coverage[i] <- 0
  d$single_dose[i] <- 0
  d$multi_dose[i] <- 0
  d$coverage[i] <- 0
  d$zero_dose[i] <- 0
  
  d <- d[order(d$cohort, d$age), ]
  cohorts <- unique(d$cohort)
  d2 <- list()
  for(i in seq_along(cohorts)){
    t <- d[d$cohort == cohorts[i], ]
    for(j in seq_along(t$age)[-1]){
      t$single_dose[j] <- ifelse(t$effective_coverage[j] == 0, t$single_dose[j-1], t$single_dose[j])
      t$multi_dose[j] <- ifelse(t$effective_coverage[j] == 0, t$multi_dose[j-1], t$multi_dose[j])
      t$coverage[j] <- ifelse(t$effective_coverage[j] == 0, t$coverage[j-1], t$coverage[j])
      t$zero_dose[j] <- ifelse(t$effective_coverage[j] == 0, t$zero_dose[j-1], t$zero_dose[j])
      t$effective_coverage[j] <- ifelse(t$effective_coverage[j] == 0, t$effective_coverage[j-1], t$effective_coverage[j])
    }
    d2[[i]] <- t
  }
  d2 <- do.call(rbind, d2)
  d2$effective_zero_dose <- 1 - d2$effective_coverage
  d2$vaccinated <- d2$coverage * d2$cohort_size
  d2$effective_vaccinated <- d2$cohort_size * d2$effective_coverage
  d2$effective_added_fvps <- round(d2$cohort_size * d2$added_effective_coverage)
  d2$fvps <- round(d2$activity_coverage * d2$cohort_size)
  d2$added_fvps <- round(d2$cohort_size * d2$added_coverage)
  d2$disease <- disease
  d2$year <- d2$cohort + d2$age
  d2$country <- country
  
  return(d2)
}
