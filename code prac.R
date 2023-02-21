
prepare_coverage <- function(COUNTRIES = "MYS", year_current = 2021){
  ### this is to clean and transform who/unicef measles coverage data
  ##0.) raw data
  d0 <- read_xlsx("data/Coverage.xlsx") # this is WUENIC routine coverage
  d <- as.data.frame(d0) 
  
  c0 <- read_xls("data/summary_measles_sias.xls", skip = 1) ## this is WHO measles-rubella campaign data
  c <- as.data.frame(c0)
  
  p0 <- read.csv("data/population.csv") # this is VIMC interpolated population estimates
  
  
  ##1.) transform routine coverage
  i <- d$COVERAGE_CATEGORY == "WUENIC" & d$ANTIGEN %in% c("MCV1", "MCV2") & d$CODE %in% "MYS"
  d <- d[i, ]
  d <- d[c("CODE",  "YEAR", "ANTIGEN", "COVERAGE")]
  names(d) <- c("country", "year", "vaccine", "coverage")
  d$coverage <- d$coverage / 1e2
  
  ## add routine age assuming mcv1 and mcv2 given at age 0 and age 2
  d$age <- NA
  i <- d$vaccine == "MCV1"
  j <- d$vaccine == "MCV2"
  d$age[i] <- 0
  d$age[j] <- 7
  d$age[d$vaccine %in% "MCV2" & d$year %in% 2017:2021] <- 1 #changes as per year reported in WUENIC
  d$coverage[is.na(d$coverage) & d$year != 2021] <- 0
  d$cohort <- d$year - d$age
  
  #picking up pop data for routine coverage and merging tables
  p1 <- p0[p0$country == "MYS" & p0$year %in% min(d$year):max(d$year), ]
  p1$cohort <- p1$year - p1$age
  names(p1)[names(p1) == "value"] <- "cohort_size"
  r <- merge(d, p1, by=c("country", "cohort", "age", "year"), all.x = TRUE)
  r$fvps <- r$coverage * r$cohort_size
  
  ### transform campaigns
  c <- c[c("ISO3 code", "Year", "Age group", "Target population", "Reached population", "% Reached")]
  names(c) <- c("country", "year", "age_group", "target", "fvps", "coverage")
  c <- c[c$country %in% COUNTRIES, ]
  c <- cbind(c, classify_age_range(c$age_group))
  c$age_from[grepl("All ages", c$age_group)] <- 1
  c$age_to[grepl("All ages", c$age_group)] <- 4
  c$age_from[grepl("Aborigines", c$age_group)] <- 1
  c$age_to[grepl("Aborigines", c$age_group)] <- 4
  c$age_range_verbatim <- NULL
  c$index <- seq_along(c$country)
  
  ## pick up demographic data for campaign
  c1 <- c
  p <- p0[p0$country == COUNTRIES & p0$year %in% min(c$year):max(c$year), ]
  p$cohort <- p$year - p$age
  names(p)[names(p) == "value"] <- "cohort_size"
  
  output <- NULL
  for (a in seq_along(c1$index)) {
    #print(a)
    t <- c1[c1$index == c1$index[a], ]
    v <- p[p$age %in% t$age_from:t$age_to & p$year == t$year, ]
    tot_pop <- sum(v$cohort_size)
    s <- merge(v, t, by = c("country", "year"), all.x = TRUE)
    s$tot_pop <- tot_pop
    s$cohort_ratio <- s$cohort_size / s$tot_pop
    s$fvps <- s$fvps * s$cohort_ratio
    s$target <- NULL
    s$age_from <- NULL
    s$age_to <- NULL
    s$age_group <- NULL
    s$cohort_ratio <- NULL
    s$tot_pop <- NULL
    s$coverage <- s$fvps / s$cohort_size
    
    output <- rbind(output, s)
  }
  
  # bind routine and campaign data
  output$vaccine <- "Measles"
  output$index <- NULL
  dat <- rbind(r, output)
  dat <- dat[order(dat$cohort, dat$age), ]
  
  ## combine rows for per cohort/age combination
  i <- duplicated(dat[c("cohort", "age")])
  any(i)
  sum(i)
  
  ## assume vaccine type not important
  dat2 <- aggregate(coverage ~ country + vaccine+ cohort + age + year, dat, sum, na.rm = TRUE )
  dat2$coverage[dat2$coverage > 1] <- 1
  i <- duplicated(dat2[c("cohort", "age")])
  dat2$disease <- "Measles"
  names(dat2) <- c("country", "coverage_set", "cohort", "age", "year", "activity_coverage", "disease")
  return(list(coverage = dat2, cohort_size = p1))
}
