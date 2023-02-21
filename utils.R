read_xlsx <- function(...) {
  oo <- options(warnPartialMatchArgs = FALSE)
  if (!is.null(oo$warnPartialMatchArgs)) {
    on.exit(options(oo))
  }
  readxl::read_xlsx(...)
}
read_xls <- function(...) {
  oo <- options(warnPartialMatchArgs = FALSE)
  if (!is.null(oo$warnPartialMatchArgs)) {
    on.exit(options(oo))
  }
  readxl::read_xls(...)
}
library(readxl)
classify_age_range <- function(x) {
  CAMPAIGN_AGE_MIN <- 0
  CAMPAIGN_AGE_MAX <- 100
  
  is_integer_like <- function(x) {
    abs(x - round(x)) < 1e-3
  }
  convert_floor <- function(x, div = 1, part = 1) {
    floor(as.integer(x) / div)
  }
  convert_ceiling <- function(x, div, inclusive, part = 2) {
    xi <- as.integer(x)
    if (any(!inclusive)) {
      if (any(!is_integer_like((xi / div)[!inclusive]))) {
        stop("FIXME")
      }
    }
    ret <- ceiling(xi / div)
    k <- !inclusive | (div == 12)
    ret[k] <- ret[k] - 1L
    ret
  }
  part <- function(x, re, n) {
    sub(re, sprintf("\\%d", n), x)
  }
  
  age_from <- age_to <- rep(NA_real_, length(x))
  
  i <- is.na(x)
  age_from[i] <- CAMPAIGN_AGE_MIN
  age_to[i] <- CAMPAIGN_AGE_MAX
  
  re <- "^([0-9]+)\\s*-\\s*([0-9]+)\\s*([mMyY]?)$"
  i <- is.na(age_from) & grepl(re, x)
  unit <- tolower(part(x[i], re, 3))
  div <- ifelse(unit == "m", 12, 1)
  age_from[i] <- convert_floor(part(x[i], re, 1), div)
  age_to[i] <- convert_ceiling(part(x[i], re, 2), div, TRUE)
  
  
  re <- "^([0-9]+)\\s*([mMyYo])\\s*-\\s*(<)?([0-9]+)\\s*([yY])\\s*([o])?\\s*(rs)?$"
  i <- is.na(age_from) & grepl(re, x)
  div <- ifelse(tolower(part(x[i], re, 2)) == "m", 12, 1)
  incl <- part(x[i], re, 3) != "<"
  age_from[i] <- convert_floor(part(x[i], re, 1), div)
  age_to[i] <- convert_ceiling(part(x[i], re, 4), 1, incl)
  
  re <- "^([0-9]+)\\s*([mMyY])\\s*-\\s*(<)?([0-9]+)\\s*([mM])$"
  i <- is.na(age_from) & grepl(re, x)
  div <- ifelse(tolower(part(x[i], re, 2)) == "m", 12, 1)
  age_from[i] <- convert_floor(part(x[i], re, 1), div)
  age_to[i] <- convert_ceiling(part(x[i], re, 4), div, TRUE)
  
  re <- "^(<)?\\s*([0-9]+)\\s*([mMyY])\\s*([o])?$"
  i <- is.na(age_from) & grepl(re, x)
  div <- ifelse(tolower(part(x[i], re, 3)) == "m", 12, 1)
  incl <- part(x[i], re, 1) != "<" 
  age_from[i] <- ifelse(incl, convert_ceiling(part(x[i], re, 2), div, incl), 0)
  age_to[i] <- convert_ceiling(part(x[i], re, 2), div, incl)
  
  re <- "^(>)\\s*([0-9]+)\\s*([mMyY])$"
  i <- is.na(age_from) & grepl(re, x)
  div <- ifelse(tolower(part(x[i], re, 3)) == "m", 12, 1)
  v <- as.integer(part(x[i], re, 2)) / div
  age_from[i] <- floor(v)
  age_to[i] <- CAMPAIGN_AGE_MAX
  
  re <- "^([0-9]+)\\s*([mM])\\s*\\+$"
  i <- is.na(age_from) & grepl(re, x)
  age_to[i] <- CAMPAIGN_AGE_MAX
  age_from[i] <- floor(as.integer(part(x[i], re, 1)) / 12)
  
  re <- "^([0-9]+)\\s*-<([0-9]+)\\s*Y$"
  i <- is.na(age_from) & grepl(re, x)
  age_from[i] <- as.integer(part(x[i], re, 1))
  age_to[i] <- as.integer(part(x[i], re, 2)) - 1L
  
  i <- is.na(age_from)
  if (any(i)) {
    message(sprintf("\t* Found %d unclassifiable age ranges:\n%s",
                    sum(i),
                    paste(sprintf("\t  - %s", sort(unique(x[i]))),
                          collapse = "\n")))
    age_from[i] <- CAMPAIGN_AGE_MIN
    age_to[i] <- CAMPAIGN_AGE_MAX
  }
  
  data.frame(age_from = age_from,
             age_to = age_to,
             age_range_verbatim = x,
             stringsAsFactors = FALSE)
}