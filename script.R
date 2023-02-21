rm(list=ls())
source("r/utils.R")
source("r/ID_functions.R")
source("code prac.R")
require("tidyr")
dat <- prepare_coverage()
cohort_size <- dat$cohort_size
coverage <- dat$coverage
coverage <- coverage[coverage$activity_coverage > 0, ]
coverage$activity_type <- ifelse(coverage$coverage_set == "Measles", "campaign", "routine")

year_min = 1980 
year_max = 2021
n_size <- 1000
child_age_max <- 17
### for each disease model do - only looking at Measles for this project
### How many doses can a cohort has received?
# coverage$n = 1
# aggregate(n ~ cohort, coverage, sum)
### start estimating effective vaccination coverage 
## todo for Shuren
## 1. to decide efficacy
## 2. generate samples of para_corr for running methodology (use R to generate sample csvs, as corr1.csv, corr2.csv, etc.)
for(run_id in 1:n_size){
  # theoretically p < c < v
  # re-sampling strategy
  
  p <- runif(1, min = 0.72, max = 0.95)
  min_ <- ifelse(p < 0.848, 0.848, p)
  c <- runif(1, min = min_, max = 0.97)
  min_ <- ifelse(p < 0.883, 0.883, c)
  v <- runif(1, min = min_, max = 0.983)

  i <- runif(1, min = -1, max = 1)
  j <- runif(1, min = -1, max = 1)
  k <- runif(1, min = -1, max = 1)
  
  tempeff <- data.frame(disease = rep("Measles", times= 3), age=0:2, single_dose_eff = c(p, c, v), multi_dose_eff = v, stringsAsFactors = FALSE)
  tempdf <- data.frame(disease = rep("Measles", times= 10), dose=1:10, corr=c(1, i, j, k , -1,-1,-1,-1,-1,-1), stringsAsFactors = FALSE)
  write.csv(tempeff, paste0("meta/eff", run_id, ".csv"), row.names = FALSE)
  write.csv(tempdf, paste0("meta/corr", run_id, ".csv"), row.names = FALSE)
}

dat <- list(NULL)
for(run_id in 1:n_size){
  input_corr_name <- paste0("meta/corr", run_id, ".csv")
  input_eff_name <- paste0("meta/eff", run_id, ".csv")
  output_file_name <- paste0("output/output.csv")
  
  para_eff <- read.csv(input_eff_name, stringsAsFactors = FALSE)
  #para_eff <- read.csv("meta/vaccine_efficacy.csv", stringsAsFactors = FALSE)
  para_corr <- read.csv(input_corr_name, stringsAsFactors = FALSE)
  c <- coverage
  ## estimate effective coverage
  o <- NULL
  c$group <- paste(c$disease, c$country, sep = "_")
  group <- unique(c$group)
  for(g in group){
    t <- c[c$group == g, ]
    d_eff_cov <- effective_coverage_history(vax_history = t, para_eff, para_corr, cohort_size, sur_rates = NULL, year_min, year_max)
    d_eff_cov <- d_eff_cov[c("disease", "country", "year", "cohort", "age", "cohort_size", 
                             "effective_coverage", "single_dose", "multi_dose")]
    o <- rbind(o, d_eff_cov)
  }
  
  ### if threshold of age at 18, the cross_pop_coverage in pre-adults is calculated as
  o$fvps <- o$cohort_size * o$effective_coverage
  o$single_dose_fvps <- o$single_dose * o$cohort_size
  o$multi_dose_fvps <- o$multi_dose * o$cohort_size
  d <- aggregate(cbind(single_dose_fvps, multi_dose_fvps, fvps, cohort_size) ~ disease + country + year, o[o$age <= child_age_max, ], sum)
  d$cross_pop_coverage_children <- d$fvps / d$cohort_size
  d$single_dose_children <- d$single_dose_fvps / d$cohort_size
  d$multi_dose_children <- d$multi_dose_fvps / d$cohort_size
  
  o$fvps <- NULL
  d2 <- merge(o, d[c("country", "year", "single_dose_children", "multi_dose_children", "cross_pop_coverage_children")])
  d2$run_id <- run_id
  dat[[run_id]] <- d2
}
dat <- do.call(rbind, dat)
write.csv(dat, output_file_name, row.names = FALSE)
require("dplyr")

v0 <- dat %>%
  group_by(year) %>%
  summarise(q1 = round(quantile(cross_pop_coverage_children, probs = 0.025), 3),
            mean = round(mean(cross_pop_coverage_children), 3), 
            q3 = round(quantile(cross_pop_coverage_children, probs = 0.975), 3))
v <- dat %>%
  group_by(cohort, age) %>%
  summarise(q1 = round(quantile(effective_coverage, probs = 0.025), 3),
            mean = round(mean(effective_coverage), 3), 
            q3 = round(quantile(effective_coverage, probs = 0.975), 3))
v1 <- v %>%
  mutate(year = cohort + age)

v2 <- dat %>%
  filter(age <= child_age_max & year == 2021) %>%
  group_by(year, age) %>%
  summarise(q1 = round(quantile(single_dose, probs = 0.025), 3),
            mean = round(mean(single_dose), 3), 
            q3 = round(quantile(single_dose, probs = 0.975), 3))

v3 <- dat %>%
  filter(age <= child_age_max & year == 2021) %>%
  group_by(year, age) %>%
  summarise(q1 = round(quantile(multi_dose, probs = 0.025), 3),
            mean = round(mean(multi_dose), 3), 
            q3 = round(quantile(multi_dose, probs = 0.975), 3))

v4 <- dat %>%
  mutate(zero_dose = 1 - single_dose - multi_dose) %>%
  filter(age <= child_age_max & year == 2021) %>%
  group_by(year, age) %>%
  summarise(q1 = round(quantile(zero_dose, probs = 0.025), 3),
            mean = round(mean(zero_dose), 3), 
            q3 = round(quantile(zero_dose, probs = 0.975), 3))


v5 <- dat %>%
  group_by(year) %>%
  summarise(q1 = round(quantile(single_dose_children, probs = 0.025), 3),
            mean = round(mean(single_dose_children), 3), 
            q3 = round(quantile(single_dose_children, probs = 0.975), 3))

v6 <- dat %>%
  group_by(year) %>%
  summarise(q1 = round(quantile(multi_dose_children, probs = 0.025), 3),
            mean = round(mean(multi_dose_children), 3), 
            q3 = round(quantile(multi_dose_children, probs = 0.975), 3))

library(ggplot2)

#visualise under 18 EVC
under18evc <- ggplot(v0, aes(year)) + geom_ribbon(aes(ymin = q1, ymax = q3) , fill = "lightblue") + 
  ylab("Under 18 EVC") + xlab("Year") + geom_line(aes(y=mean), size = 1, col = "darkblue") + 
  theme(panel.background = element_blank()) + theme(axis.line = element_line(colour = "black")) + 
  geom_hline(yintercept = 0.95) + 
  annotate("text", 2012, 0.95, hjust = 0, vjust = -1,label = "Herd Immunity Threshold") +
  ylim(0, 1)

ggsave(under18evc,filename="under18evc.png",height=5.5,width=8.8,units="in",dpi=200)

#build heatmap

v1$cohort<- NULL
dfv1<- v1[v1$year %in% 1980:2021 & v1$age %in% 0:40, ]

v1x <- dfv1 %>%
  # convert year to factor
  mutate(year=factor(year)) %>%
  # convert age to factor
  mutate(age=factor(age)) %>%
  # convert value to numeric (also converts '-' to NA, gives a warning)
  mutate(q3=as.numeric(q3)) %>%
  # convert mean to numeric (also converts '-' to NA, gives a warning)
  mutate(mean=as.numeric(mean)) %>%
  # convert mean to numeric (also converts '-' to NA, gives a warning)
  mutate(q1=as.numeric(q1))

v2x <- v1x %>%
  mutate(countfactor=cut(q3,breaks=c(0,0.2,0.5,0.8,0.9,0.95,max(q3,na.rm=T)),
                         labels=c("0-0.2","0.2-0.5","0.5-0.8","0.8-0.9","0.9-0.95","0.95-1")))

# assign text colour
textcol <- "grey40"

# further modified ggplot
q3  <- ggplot(v2x,aes(x=year,y=age,fill=countfactor))+
  geom_tile(colour="white",size=0.2)+
  guides(fill=guide_legend(title="EVC"))+
  labs(x="Year",y="Age",title="EVC since measles vaccination program began")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0),breaks=c("1980","1985","1990","1995","2000","2005","2010","2015", "2020"))+
  scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#ddf1da"),na.value = "grey90")+
  #coord_fixed()+
  theme_grey(base_size=10)+
  theme(legend.position="right",legend.direction="vertical",
        legend.title=element_text(colour=textcol),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"))

#export figure
ggsave(q3,filename="q3heatmap.png",height=5.5,width=8.8,units="in",dpi=200)

v3x <- v1x %>%
  mutate(countfactor=cut(mean,breaks=c(0,0.2,0.5,0.8,0.9,0.95,max(mean,na.rm=T)),
                         labels=c("0-0.2","0.2-0.5","0.5-0.8","0.8-0.9","0.9-0.95","0.95-1")))

# assign text colour
textcol <- "grey40"

# further modified ggplot
mean  <- ggplot(v3x,aes(x=year,y=age,fill=countfactor))+
  geom_tile(colour="white",size=0.2)+
  guides(fill=guide_legend(title="EVC"))+
  labs(x="Year",y="Age",title="EVC since measles vaccination program began")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0),breaks=c("1980","1985","1990","1995","2000","2005","2010","2015", "2020"))+
  scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#ddf1da"),na.value = "grey90")+
  #coord_fixed()+
  theme_grey(base_size=10)+
  theme(legend.position="right",legend.direction="vertical",
        legend.title=element_text(colour=textcol),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"))

#export figure
ggsave(mean,filename="meanheatmap.png",height=5.5,width=8.8,units="in",dpi=200)

v4x <- v1x %>%
  mutate(countfactor=cut(q1,breaks=c(0,0.2,0.5,0.8,0.9,0.95,max(q1,na.rm=T)),
                         labels=c("0-0.2","0.2-0.5","0.5-0.8","0.8-0.9","0.9-0.95","0.95-1")))

# assign text colour
textcol <- "grey40"

# further modified ggplot
q1  <- ggplot(v4x,aes(x=year,y=age,fill=countfactor))+
  geom_tile(colour="white",size=0.2)+
  guides(fill=guide_legend(title="EVC"))+
  labs(x="Year",y="Age",title="EVC since measles vaccination program began in 1983")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0),breaks=c("1980","1985","1990","1995","2000","2005","2010","2015", "2020"))+
  scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#ddf1da"),na.value = "grey90")+
  #coord_fixed()+
  theme_grey(base_size=10)+
  theme(legend.position="right",legend.direction="vertical",
        legend.title=element_text(colour=textcol),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=7,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        axis.text.x=element_text(size=10,colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"))

#export figure
ggsave(q1,filename="q1heatmap.png",height=5.5,width=8.8,units="in",dpi=200)

