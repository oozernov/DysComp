

Packages <- c("dplyr", "reshape", "magrittr", "tidyr", "ggplot2", 
              "lme4", "lmerTest","emmeans", "sjstats", "plotrix","dabestr","lmerTest",
              "lmPerm","gridExtra", "grid","ggpubr")
lapply(Packages, library, character.only = TRUE)


########################## Load and organize data ##########################

## load adult data
adult = read.csv("~/Dropbox (MIT)/Com_Dys_2016_data/adult_analysis/phonwords_data_a_010719.csv", header = T, na.strings = c("", "NA"))

#exclude the 'd' continuum from adult dataset
adult=adult %>%
  filter(grepl("gi", tolower(soundfile)))

#Or combine d and g (lower Ganong effect)
# adult$soundfile = gsub("dask","giss",adult$soundfile)
# adult$soundfile = gsub("dash","gift",adult$soundfile)
# adult$soundfile = gsub("dath","gith",adult$soundfile)
# 
# adult$leftprompt = gsub("D","G",adult$leftprompt)
# adult$leftprompt = gsub("T","K",adult$leftprompt)
# adult$response = gsub("D","G",adult$response)
# adult$response = gsub("T","K",adult$response)

groups_a = read.csv("~/Dropbox (MIT)/Com_Dys_2016_data/adult_analysis/adult_groups_012418.csv", header = T, na.strings = c("", "NA"))
## load child data
setwd("~/Dropbox (MIT)/Com_Dys_2016_data")
child = read.csv("~/Dropbox (MIT)/Com_Dys_2016_data/child_phonwords_data_122718.csv", header = T, na.strings = c("", "NA"))
groups_c = read.csv("~/Dropbox (MIT)/Com_Dys_2016_data/groups_041817.csv", header = T, na.strings = c("", "NA"))

setwd("~/Dropbox (MIT)/Com_Dys_2016_data/")

# add DD columns to the data
child = merge(child, groups_c, by = "PartID")
adult = merge(adult, groups_a, by = "Subject")
names(adult)[names(adult) == 'Subject'] <- 'PartID'

child$group<-ifelse(child$Typical==1, "Typ","Dys")
child$group<-as.factor(child$group)
adult$X <- NULL
adult$X.y <- NULL
adult$X.1 <- NULL
child$Subject <- NULL
child$DD <- NULL
child$Typical <- NULL

##combine the datasets
adult$age<-"Adult"
child$age<-"Child"

d = bind_rows(adult, child) #changed from 'combine' to 'bind_rows'
#d <- rename(d, c("source" = "age"))
#names(d)[names(d) == 'source'] <- 'age'

groups_c$group = ifelse(groups_c$Typical==1, "Typ","Dys")
groups_a$X <- NULL
groups_a$X.1 <- NULL
groups_c$DD <- NULL
groups_c$Typical <- NULL
#groups_a <- rename(groups_a, c("Subject" = "PartID"))
names(groups_a)[names(groups_a) == 'Subject'] <- 'PartID'

groups = bind_rows(groups_a, groups_c)

#d$PartID = as.factor(d$PartID)
# make new columns that will let us separate by continuum
new_columns = colsplit(d$soundfile, "-", names=c("word", "rest"))
d = cbind(d, new_columns)
new_columns = colsplit(d$rest, "", names=c("n_step","w", "x", "y", "z"))
d = cbind(d, new_columns)
d=d%>%filter(!(PartID %in% c('READ_6103')))

# get rid of columns we don't need
d$rest <- NULL
d$w <- NULL
d$x <- NULL
d$y <- NULL
d$z <- NULL

d$com_cond = paste(d$group, "-", d$age)

d_g = d %>% filter(word %in% c('gift', 'gith', 'giss'))
d_g$phoneme_g = ifelse(d_g$response == 'G', 1, 0)
d_conts = d_g %>% dplyr::group_by(com_cond, word, n_step) %>%
  dplyr::summarise(mean_rt = mean(keyRT, na.rm = T))
d_conts$word_pair = ifelse(d_conts$word == "gift", "gift-kift",
                           ifelse(d_conts$word == "gith", "gith-kith", "giss-kiss"))

m2=(lmer(keyRT~word*group*age+(1|PartID), data=d_g, REML=TRUE))
anova(m2)
ggplot(d_conts, aes(colour=word_pair, y=mean_rt, x=n_step)) + 
  geom_point() + geom_line() + 
  #scale_x_continuous(breaks = x_ticks, limits = c(1, 7)) + 
  #scale_y_continuous(breaks = y_ticks, limits = c(0, 100)) + 
  labs(x = "Step", y = "mean RT", title = "Mean RT") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.position = "right", 
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~com_cond, ncol = 2)


##### 6. Create d_export #####
d_export = d_g %>% filter(word == 'gith') %>% 
  dplyr::group_by(PartID, n_step, group, age) %>% 
  dplyr::summarise(mean_rt = mean(keyRT, na.rm = T))
d_export = distinct(data.frame(d_export[, c("PartID")]))
names(d_export)=c("PartID")

##### 7. glm logistic regression to extract slope/inflection #####
calc_g = d_g %>%
  dplyr::group_by(PartID, word, n_step, group, age, com_cond) 
#dplyr::summarise(prop_g = mean(phoneme_g, na.rm = T))

d_glm_fit = data.frame(PartID = rep(0, 10), 
                       word = rep(0, 10), 
                       group = rep(0, 10), 
                       age = rep(0, 10), 
                       com_cond = rep(0, 10),
                       slope_coef = rep(0, 10), 
                       inflection_step = rep(0, 10), 
                       value_at_gith_i = rep(0, 10), 
                       deviance = rep(0, 10))
index = 0
for (sub in d_export$PartID) { 
  index = index + 1
  temp = calc_g %>% filter(PartID == sub)
  temp_gift = temp %>% filter(word == "gift")
  temp_gith = temp %>% filter(word == "gith")
  temp_giss = temp %>% filter(word == "giss")
  
  temp_group = as.character(temp$group[1])
  temp_age = as.character(temp$age[1])
  temp_com_cond = as.character(temp$com_cond[1])
  
  model_gift <- glm(keyRT ~ n_step, data = temp_gift)
  model_gith <- glm(keyRT ~ n_step, data = temp_gith)
  model_giss <- glm(keyRT ~ n_step, data = temp_giss)
  
  b_gith = as.numeric(model_gith$coefficients[1])
  m_gith = as.numeric(model_gith$coefficients[2])
  i_gith = -b_gith/m_gith #-intercept/slope
  dev_fit_gith = as.numeric(deviance(model_gith))
  gith_x = data.frame('n_step' = i_gith)
  value_gith_frame = data.frame(gith_x$n_step, predict(model_gith, list(n_step = gith_x$n_step), type = 'response'))
  value_gith = value_gith_frame[[2]]
  new2 = c(sub, "gith", temp_group, temp_age, temp_com_cond, m_gith, i_gith, value_gith, dev_fit_gith)
  d_glm_fit[index, ] = new2
  index = index + 1
  
  b_gift = as.numeric(model_gift$coefficients[1])
  m_gift = as.numeric(model_gift$coefficients[2])
  i_gift = -b_gift/m_gift
  dev_fit_gift = as.numeric(deviance(model_gift))
  value_gift_frame = data.frame(gith_x$n_step, predict(model_gift, list(n_step = gith_x$n_step), type = 'response'))
  value_gift = value_gift_frame[[2]]
  new1 = c(sub, "gift", temp_group, temp_age, temp_com_cond, m_gift, i_gift, value_gift, dev_fit_gift)
  d_glm_fit[index, ] = new1
  index = index + 1
  
  b_giss = as.numeric(model_giss$coefficients[1])
  m_giss = as.numeric(model_giss$coefficients[2])
  i_giss = -b_giss/m_giss
  dev_fit_giss = as.numeric(deviance(model_giss))
  value_giss_frame = data.frame(gith_x$n_step, predict(model_giss, list(n_step = gith_x$n_step), type = 'response'))
  value_giss = value_giss_frame[[2]]
  new3 = c(sub, "giss", temp_group, temp_age, temp_com_cond, m_giss, i_giss, value_giss, dev_fit_giss)
  d_glm_fit[index, ] = new3
}

d_glm_fit$slope_coef = as.numeric(d_glm_fit$slope_coef)
d_glm_fit$inflection_step = as.numeric(d_glm_fit$inflection_step)
d_glm_fit$deviance = as.numeric(d_glm_fit$deviance)

m1=(lmer(slope_coef~word*group*age+(1|PartID), 
                         weights=deviance, data=d_glm_fit, REML=TRUE))

anova(m1)
lsmeans(m1, list(pairwise ~ word), adjust = "tukey") 
d_g$group<-as.factor(d_g$group)
d_g$age<-as.factor(d_g$age)
d_m = d_g %>%
  dplyr::group_by(PartID, word)%>%
  dplyr::summarise(mean_rt = mean(keyRT, na.rm = T))%>%
  left_join(d_glm_fit, by = c("PartID", "word"))

str(d_m)
m2=(lmer(mean_rt~word*group*age+(1|PartID), data=d_m, REML=TRUE))
anova(m2)
eta_sq(m2)
lsmeans(m2, list(pairwise ~ group), adjust = "tukey")#RT varies in Dys, but not Typ as function of context
lsmeans(m2, list(pairwise ~ age), adjust = "tukey") #children slower
