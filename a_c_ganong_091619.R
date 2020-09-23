###Ganong Analysis
#Ola Ozernov-Palchik oozernov@mit.edu

Packages <- c("dplyr", "reshape", "magrittr", "tidyr", "ggplot2", 
              "lme4", "lmerTest","emmeans", "sjstats", "plotrix","dabestr","lmerTest",
              "lmPerm","gridExtra", "grid","ggpubr",'sjmisc','relaimpo')
lapply(Packages, library, character.only = TRUE)


########################## Load and organize data ##########################
setwd("~/Dropbox (MIT)/Com_Dys_2016_data/final_for_sharing/final_code_JT/data")

## load adult data
adult = read.csv("phonwords_data_a_010719.csv", header = T, na.strings = c("", "NA"))

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

groups_a = read.csv("adult_groups_012418.csv", header = T, na.strings = c("", "NA"))
## load child data
setwd("Com_Dys_2016_data")
child = read.csv("child_phonwords_data_122718.csv", header = T, na.strings = c("", "NA"))
groups_c = read.csv("groups_041817.csv", header = T, na.strings = c("", "NA"))


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


# define some vectors to use for graphs later
x_ticks = c(1,2,3,4,5,6,7)
y_ticks = c(0,10,20,30,40,50,60,70,80,90,100)
y_new_ticks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

####################### Raw percentage /g/ response ########################
# d$word<-NULL #had a duplication error for 'word'
d_g = d %>% filter(word %in% c('gift', 'gith', 'giss'))
d_g$phoneme_g = ifelse(d_g$response == 'G', 1, 0)
d_conts = d_g %>% dplyr::group_by(com_cond, word, n_step) %>%
  dplyr::summarise(percent_g = mean(phoneme_g, na.rm = T))
d_conts$word_pair = ifelse(d_conts$word == "gift", "gift-kift",
                           ifelse(d_conts$word == "gith", "gith-kith", "giss-kiss"))

ggplot(d_conts, aes(colour=word_pair, y=percent_g*100, x=n_step)) + 
  geom_point() + geom_line() + 
  scale_x_continuous(breaks = x_ticks, limits = c(1, 7)) + 
  scale_y_continuous(breaks = y_ticks, limits = c(0, 100)) + 
  labs(x = "Step", y = "% /g/ response", title = "Raw /g/ Response") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.position = "right", 
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~com_cond, ncol = 2)

########################## Get sample sizes by group ##########################
counts = d %>% dplyr::group_by(PartID, age) %>%
  dplyr::summarise(m = mean(n_step))
counts = merge(counts, groups, by = "PartID")
dplyr::count(counts, age, group)


##### 6. Create d_export #####
d_export = d_g %>% filter(word == 'gith') %>% 
  dplyr::group_by(PartID, n_step, group, age) %>% 
  dplyr::summarise(m_g = mean(phoneme_g, na.rm = T))
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
  
  model_gift <- glm(phoneme_g ~ n_step, family = binomial(link = 'logit'), data = temp_gift)
  model_gith <- glm(phoneme_g ~ n_step, family = binomial(link = 'logit'), data = temp_gith)
  model_giss <- glm(phoneme_g ~ n_step, family = binomial(link = 'logit'), data = temp_giss)
  
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
d_glm_fit$word<-as.factor(d_glm_fit$word)
d_glm_fit$group<-as.factor(d_glm_fit$group)
d_glm_fit$age<-as.factor(d_glm_fit$age)

######################Deal with outliers#######################
#bind inflection step to acceptable values

#d_glm_fit$inflection_step = ifelse(d_glm_fit$inflection_step < 1 | d_glm_fit$inflection_step > 7, 0, d_glm_fit$inflection_step)

##bind inflection to be within range for  steps 
sum(d_glm_fit$inflection_step < 1 |d_glm_fit$inflection_step > 7) #31/288 (11%) #what percentage out of bounds?
d_glm_fit<-d_glm_fit %>% filter(inflection_step > 0 & inflection_step < 8)

# Is deviance predicted by group/age?
t.test(deviance ~ group, data = d_glm_fit)
t.test(deviance ~ age, data = d_glm_fit)
#Transform slope
d_glm_fit$slope_coef_t<-as.numeric(abs(d_glm_fit$slope_coef))
d_glm_fit$slope_coef_t<-log10(d_glm_fit$slope_coef_t)

########################## Compare effects for slope and inflection ##########################
# slope across conditions
model_slope<-lm(slope_coef_t~group*age*word, data=d_glm_fit)
anova(model_slope)
eta_sq(model_slope)

lsmeans(model_slope, list(pairwise ~ group), adjust = "tukey")
lsmeans(model_slope, list(pairwise ~ age), adjust = "tukey")

#slope at gith only
d_glm_fit_gith=d_glm_fit%>%dplyr::filter(word=="gith") 
model_gith<-lm(slope_coef_t~group*age, data=d_glm_fit_gith)
anova(model_gith)
eta_sq(model_gith)
#posthoc
lsmeans(model_gith, list(pairwise ~ age), adjust = "tukey")
lsmeans(model_gith, list(pairwise ~ group|age), adjust = "tukey")



########################## Plot slope effect ##########################
d_glm_m_slope = d_glm_fit %>%
  dplyr::group_by(PartID, group,com_cond) %>%
  dplyr::summarise(mean = mean(slope_coef_t, na.rm = T))

multi.group <- 
  d_glm_m_slope %>%
  dabest(com_cond,mean, 
         idx = list(c("Dys - Adult", "Typ - Adult")
                    ,c("Dys - Child","Typ - Child")),
         paired = FALSE
  )
plot(multi.group, color.column = group,rawplot.ylabel = "Slope")


########################## Get lexicality effect ##########################
gift_pos = d_g %>% filter(word == 'gift') %>% 
  dplyr::group_by(PartID, n_step, group, age) %>% 
  dplyr::summarise(m_g = mean(phoneme_g, na.rm = T))
giss_pos = d_g %>% filter(word == 'giss') %>% 
  dplyr::group_by(PartID, n_step, group, age) %>% 
  dplyr::summarise(m_g = mean(phoneme_g, na.rm = T))
gith_pos = d_g %>% filter(word == 'gith') %>% 
  dplyr::group_by(PartID, n_step, group, age)%>% 
  dplyr::summarise(m_g = mean(phoneme_g, na.rm = T))

lexical_g = data.frame(PartID = rep(0, 10), 
                       step = rep(0, 10), 
                       diff = rep(0, 10), 
                       group = rep(0, 10), 
                       age = rep(0, 10)) 
row_index = 0
for (sub in gift_pos$PartID) {
  row_index = row_index + 1
  current_diff = gift_pos$m_g[row_index] - giss_pos$m_g[row_index]
  new = c(sub, gift_pos$n_step[row_index], 
          current_diff, as.character(gift_pos$group[row_index]), 
          as.character(gift_pos$age[row_index]))
  lexical_g[row_index, ] = new
}

lexical_g$diff = as.numeric(lexical_g$diff)
lexical_g_pos = lexical_g %>% 
  dplyr::group_by(PartID, group, age) %>% 
  dplyr::summarise(mean_diff = mean(diff, na.rm = TRUE))


########################## Analyze lexicality effect ##########################

##analyze inflection effect

### run without gith (same results as all three conditions)
d_glm_fit_n_gith=d_glm_fit%>%filter(word=="gift"|word=="giss")
LME_model_i_nogith=lmer(inflection_step~word*group*age+(1|PartID), 
                         weights=deviance, data=d_glm_fit_n_gith, REML=TRUE)
anova(LME_model_i_nogith)
eta_sq(LME_model_i_nogith)

###Post hoc comparisons
lsmeans(LME_model_i_nogith, list(pairwise ~ group|word|age), adjust = "tukey") #run only for gith
lsmeans(LME_model_i_nogith, list(pairwise ~age), adjust = "tukey") #run only for gith

typing.lsm = lsmeans(LME_model_i_nogith, pairwise ~ word|group|age, glhargs=list())
print(typing.lsm, omit=1)
plot(typing.lsm[[2]]) #plot effects for ganong 

###lexicality by step
lexical_model=(lmer(diff~step*group*age+(1|PartID),data=lexical_g, REML=FALSE))
anova(lexical_model)
eta_sq(lexical_model)

lsmeans(lexical_model, list(pairwise ~ step|group|age), adjust = "tukey") #run only for gith

########################## Plot lexicality effect ##########################
multi.two.group.unpaired_a <- 
  lexical_g_pos %>%
  filter(age=="Adult")%>%
  dabest(group,mean_diff,
         idx = list(c("Dys", "Typ")),
         paired = FALSE
  )
a<-plot(multi.two.group.unpaired_a,
        rawplot.ylim = c(-0.3, 1),
        rawplot.ylabel = "Adult Lexical Effect",
        effsize.ylabel = "Effect Size")


multi.two.group.unpaired_c <- 
  lexical_g_pos %>%
  filter(age=="Child")%>%
  dabest(group,mean_diff,
         idx = list(c("Dys", "Typ")),
         paired = FALSE
  )
c<-plot(multi.two.group.unpaired_c,
        rawplot.ylim = c(-0.3, 1),
        rawplot.ylabel = "Child Lexical Effect",
        effsize.ylabel = "Effect Size")

grid.arrange(a,c,ncol = 2)


# plot lexicality effects by step and by group

lexical_g_total <- lexical_g %>%
  group_by(step, group, age) %>%
  dplyr::summarize(mean.diff = mean(diff, na.rm = TRUE),
                   sd.diff = sd(diff, na.rm = TRUE),
                   n.diff = n()) %>%
  ungroup() %>%
  mutate(se.diff = sd.diff / sqrt(n.diff),
         lower.ci.diff = mean.diff - qt(1 - (0.05 / 2), n.diff - 1) * se.diff,
         upper.ci.diff = mean.diff + qt(1 - (0.05 / 2), n.diff - 1) * se.diff) %>%
  mutate_at(vars(step, group), as.factor)

ggplot(data = lexical_g_total,
       aes(x = step, y = mean.diff,  ymin = lower.ci.diff, ymax = upper.ci.diff,
           group = group, fill = group)) +
  geom_ribbon(alpha = 0.5) +
  geom_point() +
  geom_line() +
  labs(x = 'Step', y = "Lexical Effect", title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = 'white'),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap( ~ as.factor(age), nrow = 2)+theme(
    strip.text = element_text(size=14),
    legend.text = element_text(size=16),
    legend.title = element_blank())


##### 11. Group logistic models-Make Graphs #####
d_g$age<-as.factor(d_g$age)
ganong_models = read.csv("~/Dropbox (MIT)/Com_Dys_2016_data/final_for_sharing/final_code_JT/data/ganong_models.csv", header = T)
raw_df_list = list()
predict_df_list = list()

for (i in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) {
  temp_params = ganong_models[i,]
  temp_pos = d_g %>% 
    filter(word == as.character(temp_params[[1]]) & 
             group == as.character(temp_params[[2]]) & 
             age == as.character(temp_params[[3]])) %>% 
    group_by(n_step) %>% 
    summarise(m_g = mean(phoneme_g, na.rm = TRUE))
  temp_model <- glm(m_g ~ n_step, family = binomial(link = 'logit'), 
                    data = temp_pos)
  new_x <- data.frame(n_step = seq(min(temp_pos$n_step), 
                                   max(temp_pos$n_step), len = 100))
  temp_data = data.frame(new_x$n_step, predict(temp_model, list(n_step = new_x$n_step), type = 'response'))
  names(temp_data) = c("n_step", "value")
  temp_data$word = temp_params[[1]]
  temp_data$group = temp_params[[2]]
  temp_data$age = temp_params[[3]]
  names(temp_pos) = c("n_step", "value")
  temp_pos$word = temp_params[[1]]
  temp_pos$group = temp_params[[2]]
  temp_pos$age = temp_params[[3]]
  
  temp_pos$interact = interaction(temp_pos$word, temp_pos$group, temp_pos$age)
  temp_data$interact = interaction(temp_data$word, temp_data$group, temp_data$age)
  
  raw_df_list[[i]] <- temp_pos
  predict_df_list[[i]] <- temp_data
}

raw_df = do.call("rbind", raw_df_list)
predict_df = do.call("rbind", predict_df_list)

# make four plots
plot1 = ggplot() +
  geom_line(data = predict_df[1:300,], aes(x = n_step, y = value, group = word, color = word), size = 1.05) +
  scale_color_manual(values = c('darkorange', 'purple', 'forestgreen')) +
  geom_line(y = 0.5, color = 'black') +
  geom_line(data = data.frame(raw_df_list[1]), aes(x = n_step, y = value), linetype = "dotted", color = 'darkorange', size = 1) +
  geom_line(data = data.frame(raw_df_list[2]), aes(x = n_step, y = value), linetype = "dotted", color = 'forestgreen', size = 1) +
  geom_line(data = data.frame(raw_df_list[3]), aes(x = n_step, y = value), linetype = "dotted", color = 'purple', size = 1) +
  geom_point(data = data.frame(raw_df_list[1]), aes(x = n_step, y = value), color = 'darkorange', size = 1.5, shape = 1) +
  geom_point(data = data.frame(raw_df_list[2]), aes(x = n_step, y = value), color = 'forestgreen', size = 1.5, shape = 1) +
  geom_point(data = data.frame(raw_df_list[3]), aes(x = n_step, y = value), color = 'purple', size = 1.5, shape = 1) +
  scale_y_continuous(limits = c(0, 1.0)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7)) +
  labs(x = 'Step', y = "Proportion /g/ Response", title = "Dys Adult") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.position = "none")

plot2 = ggplot() +
  geom_line(data = predict_df[301:600,], aes(x = n_step, y = value, group = word, color = word), size = 1.05) +
  scale_color_manual(values = c('darkorange', 'purple', 'forestgreen')) +
  geom_line(y = 0.5, color = 'black') +
  geom_line(data = data.frame(raw_df_list[4]), aes(x = n_step, y = value), linetype = "dotted", color = 'darkorange', size = 1) +
  geom_line(data = data.frame(raw_df_list[5]), aes(x = n_step, y = value), linetype = "dotted", color = 'forestgreen', size = 1) +
  geom_line(data = data.frame(raw_df_list[6]), aes(x = n_step, y = value), linetype = "dotted", color = 'purple', size = 1) +
  geom_point(data = data.frame(raw_df_list[4]), aes(x = n_step, y = value), color = 'darkorange', size = 1.5, shape = 1) +
  geom_point(data = data.frame(raw_df_list[5]), aes(x = n_step, y = value), color = 'forestgreen', size = 1.5, shape = 1) +
  geom_point(data = data.frame(raw_df_list[6]), aes(x = n_step, y = value), color = 'purple', size = 1.5, shape = 1) +
  scale_y_continuous(limits = c(0, 1.0)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7)) +
  labs(x = 'Step', y = "Proportion /g/ Response", title = "Dys Child") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.position = "none")

plot3 = ggplot() +
  geom_line(data = predict_df[601:900,], aes(x = n_step, y = value, group = word, color = word), size = 1.05) +
  scale_color_manual(values = c('darkorange', 'purple', 'forestgreen')) +
  geom_line(y = 0.5, color = 'black') +
  geom_line(data = data.frame(raw_df_list[7]), aes(x = n_step, y = value), linetype = "dotted", color = 'darkorange', size = 1) +
  geom_line(data = data.frame(raw_df_list[8]), aes(x = n_step, y = value), linetype = "dotted", color = 'forestgreen', size = 1) +
  geom_line(data = data.frame(raw_df_list[9]), aes(x = n_step, y = value), linetype = "dotted", color = 'purple', size = 1) +
  geom_point(data = data.frame(raw_df_list[7]), aes(x = n_step, y = value), color = '#F8766D', size = 1.5, shape = 1) +
  geom_point(data = data.frame(raw_df_list[8]), aes(x = n_step, y = value), color = '#619CFF', size = 1.5, shape = 1) +
  geom_point(data = data.frame(raw_df_list[9]), aes(x = n_step, y = value), color = '#00BA38', size = 1.5, shape = 1) +
  scale_y_continuous(limits = c(0, 1.0)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7)) +
  labs(x = 'Step', y = "Proportion /g/ Response", title = "Typ Adult") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.position ="none")

plot4 = ggplot() +
  geom_line(data = predict_df[901:1200,], aes(x = n_step, y = value, group = word, color = word), size = 1.05) +
  scale_color_manual(values = c('darkorange', 'purple', 'forestgreen')) +
  geom_line(y = 0.5, color = 'black') +
  geom_line(data = data.frame(raw_df_list[10]), aes(x = n_step, y = value), linetype = "dotted", color = 'darkorange', size = 1) +
  geom_line(data = data.frame(raw_df_list[11]), aes(x = n_step, y = value), linetype = "dotted", color = 'forestgreen', size = 1) +
  geom_line(data = data.frame(raw_df_list[12]), aes(x = n_step, y = value), linetype = "dotted", color = 'purple', size = 1) +
  geom_point(data = data.frame(raw_df_list[10]), aes(x = n_step, y = value), color = 'darkorange', size = 1.5, shape = 1) +
  geom_point(data = data.frame(raw_df_list[11]), aes(x = n_step, y = value), color = 'forestgreen', size = 1.5, shape = 1) +
  geom_point(data = data.frame(raw_df_list[12]), aes(x = n_step, y = value), color = 'purple', size = 1.5, shape = 1) +
  scale_y_continuous(limits = c(0, 1.0)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7)) +
  labs(x = 'Step', y = "Proportion /g/ Response", title = "Typ Child") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.position = 'none')

ggarrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


######### Create individual logistic graph ##########
subID = 5133 # you can change this to whatever you want
subword = "giss" # you can change this to whatever you want

sub_pos = d_g %>% filter(PartID == subID & word == subword) %>% 
  group_by(n_step) %>% 
  summarise(m_g = mean(phoneme_g, na.rm = TRUE))
sub_model <- glm(m_g ~ n_step, family = binomial(link = 'logit'), data = sub_pos)
sub_x <- data.frame(n_step = seq(min(sub_pos$n_step), max(sub_pos$n_step), len = 100))
sub_data = data.frame(sub_x$n_step, predict(sub_model, list(n_step = sub_x$n_step), type = 'response'))
names(sub_data) = c("n_step", "value")
sub_data$word = subword
names(sub_pos) = c("n_step", "value")
sub_pos$word = subword

ggplot() +
  geom_line(data = sub_data, aes(x = n_step, y = value, group = word, color = word)) +
  geom_line(data = sub_pos, aes(x = n_step, y = value), linetype = "dotted", color = '#F8766D', size = 0.7) +
  geom_point(data = sub_pos, aes(x = n_step, y = value), color = '#F8766D', size = 1, shape = 1) +
  scale_y_continuous(limits = c(0, 1.0)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7)) +
  labs(x = 'Step', y = "Proportion /g/ Response", title = paste(subID, subword, sep = ' ')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.position = 'right')


#write.csv(d_export, "adult_dys_comp_ganong_020518.csv", row.names = FALSE) ## change to current date


m <- glmer(phoneme_g ~ trialNum + age + group + word +(1 | PartID), data = d_g, family = binomial, control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 10)

m <- glmer(phoneme_g ~ trialNum + age + group + word +(1 | PartID), data = d_g, family = binomial, control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 10)
# print the mod results without correlations among fixed effects
print(m, corr = FALSE)



############## dprimes ###############

test<-d_g %>% filter(n_step==1 | n_step==7)
test$hit = ifelse(test$n_step == 1 & test$phoneme_g == 1, 1, 0)
test$fa = ifelse(test$n_step == 1 & test$phoneme_g == 0, 1, 0)
test$miss = ifelse(test$n_step == 7 & test$phoneme_g == 1, 1, 0)
test$cr = ifelse(test$n_step == 7 & test$phoneme_g == 0, 1, 0)

library(dplyr)
detach(package:plyr)
detach(package:Rmisc)

### Calculate d-prime
test2<-tbl_df(test)
test2<-test2 %>% filter(word=='gith')

d_prime <- test2 %>% group_by(PartID) %>%
  summarise(hit = sum(hit), 
            fa = sum(fa), 
            miss = sum(miss),
            cr=sum(cr))
library(neuropsychology)
d_prime<-na.omit(d_prime)
indices<-dprime(d_prime$hit, d_prime$miss, d_prime$fa, d_prime$cr)
d_prime$d<-is.numeric(indices$dprime)
omit<-d_prime%>%filter(d<0 | d==0)

##################individual difference######################################

beh = read.csv("~/Dropbox (MIT)/GitHub/DysComp2/Ind_Diff/all_beh.csv")
names(beh)[names(beh)=="ID"] <- "PartID"
beh$a_c<-as.factor(beh$a_c)
beh$Study<-as.factor(beh$Study)

lexical_g_5<-lexical_g %>%filter(step==5) #lexical effect at step 5
slope_inf = d_glm_fit %>%
  dplyr::group_by(PartID, group) %>%
  dplyr::summarise(mean_s = mean(slope_coef_t, na.rm = T),mean_i=mean(inflection_step))


beh_d=merge(beh,lexical_g_5,"PartID")
beh_slope=merge(beh_d,d_glm_fit_gith,"PartID")
beh_slope2=merge(beh_d,slope_inf,"PartID")

###Model lexical effects
#across everyone
m1<-lm(diff~Age+Study+a_c+Vocab+IQ+Elision+Blending+NonWord+LC, data=beh_d) #elision is main
step <- stepAIC(m1, direction = "both",steps = 1000)
step$anova 

#beh_d2<-beh_d%>%filter(Study!='READER')
beh_d_d<-beh_d%>%filter(DD==1)
beh_d_t<-beh_d%>%filter(DD==0)

#in dyslexia group only
m2<-lm(diff~Age+Study+a_c+Vocab+IQ+Elision+Blending+NonWord+RANO+LC+Vocab*Elision, data=beh_d_d) #elision is main
step2 <- stepAIC(m2, direction = "both",steps = 1000)
step2$anova #diff ~ Study + Vocab + Elision + Blending + Vocab:Elision
theme_set(theme_sjplot())
sjPlot::plot_model(m2, type = "int", terms = c("Elision", "Vocab"),mdrt.values = "minmax") #plot interaction

beh_d_d_a<-beh_d_d%>%filter(a_c=='a')
beh_d_d_c<-beh_d_d%>%filter(a_c=='c')

#in adults with dys only
m4<-lm(diff~Age+IQ+Blending+NonWord+RANO+LC+Vocab*Elision, data=beh_d_d_a) #elision is main
step4 <- stepAIC(m4, direction = "both",steps = 1000)
step4$anova #diff ~ Vocab + Elision + Vocab:Elision
sjPlot::plot_model(m4, type = "int", terms = c("Elision","Vocab"),mdrt.values = "minmax")
#sjPlot::plot_model(m4, type = "int", terms = c("Elision","Vocab"),mdrt.values = "meansd")

#in kids with dys only
m5<-lm(diff~Age+IQ+Blending+NonWord+RANO+LC+Vocab*Elision,data=beh_d_d_c) #elision is main
step5 <- stepAIC(m5, direction = "both",steps = 1000)
step5$anova #diff ~ Blending

m3<-lm(diff~Age+Study+a_c+Vocab+IQ+Elision+Blending+NonWord+RANO+LC+Elision*Vocab, data=beh_d_t) #elision is main
step3 <- stepAIC(m3, direction = "both",steps = 1000)
step3$anova #diff ~ Elision + NonWord

sjPlot::plot_model(m3, type = "pred", terms = c("Elision"))
plot(beh_d_t$Elision,beh_d_t$diff)

 # age is the main vairiable so divide by age

beh_d_d_a<-beh_d_d%>%filter(a_c=='a')
beh_d_d_c<-beh_d_d%>%filter(a_c=='c')

anova(lm(diff~Age+Vocab+IQ+DigitsFw+DigitsBck+Elision+Blending+NonWord+LC, data=beh_d_d)) #age is main
anova(lm(diff~Study+Age+Vocab+IQ+Elision+Blending+NonWord+LC, data=beh_d_d_c)) #el


##OLD
summary(lm(slope_coef_t~group.x*age_group+ppvt+KBITss+ctel+ctbw+ran_2, data=beh_slope))
summary(lm(mean~group.x*age_group+ppvt+KBITss+ctel+ctbw+ran_2, data=beh_slope))

beh_d_dys<-beh_d%>%filter(beh_d$group.x=="Dys")
beh_d_typ<-beh_d%>%filter(beh_d$group.x=="Typ")
beh_d_c<-beh_d%>%filter(beh_d$age_group==1)
beh_d_a<-beh_d%>%filter(beh_d$age_group==2)
beh_a_dys<-beh_d%>%filter(beh_d$age_group==2 &beh_d$group.x=="Dys")
beh_a_typ<-beh_d%>%filter(beh_d$age_group==2 &beh_d$group.x=="Typ")
beh_c_dys<-beh_d%>%filter(beh_d$age_group==1 &beh_d$group.x=="Dys")
beh_c_typ<-beh_d%>%filter(beh_d$age_group==1 &beh_d$group.x=="Typ")

cor.test(beh_d$KBITss,beh_d$diff)

cor.test(beh_d_a$ppvt,beh_d_a$diff)
cor.test(beh_d_c$ppvt,beh_d_c$diff)

cor.test(beh_a_dys$ppvt,beh_a_dys$diff)
cor.test(beh_a_typ$ppvt,beh_a_typ$diff)
cor.test(beh_c_dys$ppvt,beh_c_dys$diff)
cor.test(beh_c_typ$ppvt,beh_c_typ$diff)

##export some information for modeling
d_glm_all = d_glm_fit %>%
  dplyr::group_by(word,com_cond) %>%
  dplyr::summarise(slope = mean(slope_coef_t, na.rm = T),inflection=mean(inflection_step,na.rm=T))

