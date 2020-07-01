# phon_rest.R
# Ola Ozernov-Palchik 
# Jul 2019
#2 Jan 2019
#
# MIT Speech Perception study, Phonemic Restoration
#

#### Setup ####
setwd("~/Dropbox (MIT)/Com_Dys_2016_data/final_for_sharing/final_code_JT")
Packages <- c("dplyr", "readr", "magrittr", "tidyr", "ggplot2", 
              "lme4", "lmerTest","emmeans", "sjstats", "plotrix","dabestr","lmerTest",
              "lmPerm","gridExtra", "grid","ggpubr")
lapply(Packages, library, character.only = TRUE)


#### Load and organize the data ####
# child data
d_c <- read_csv("data/child_phonresto_data_122718.csv")
groups_c <- read_csv("data/groups_041817.csv")

d_c %<>%
    # add group columns to the dataframe
    inner_join(groups_c, by = "PartID") %>%
    # make columns to distinguish condition and word set (i.e. bag-tag)
    separate(cond, into = c("pres", "condition"), sep = "-", remove = FALSE) %>%
    separate(soundfile, into = c("word1", "word2", "pres1", "condition1"), sep = "-", remove = FALSE) %>%
    mutate(word_set = paste(word1, word2, sep = "-"),
           item = paste(word_set, condition, sep = "-"),
           # calculate correct when missing
           correct = ifelse(is.na(correct),
                            ifelse(substr(cond, 1, 1) == substr(response,
                                                                nchar(response),
                                                                nchar(response)),
                                   1, 0),
                            correct)) %>%
    select(-pres, -word1, -word2, -pres1, -condition1)

# remove bad items (based on adult mturk data)
bad_items <- c('brown-crown', 'drill-grill', 'pear-bear', 'pole-goal')
d_c %<>%
    filter(!word_set %in% bad_items) %>%
    na.omit()

# child sample size
counts <- d_c %>%
    group_by(PartID) %>%
    summarize(m = mean(keyRT)) %>%
    left_join(groups_c, by = "PartID") %>%
    mutate(group = ifelse(DD == 1, "Dys", "Typ"))
count(counts, group)

# adult data
d_a <- read_csv("data/phonresto_data_a_010719.csv")
groups_a <- read_csv("data/adult_groups_012418.csv")

key <- read_csv("data/phonresto_key.csv") %>%
    select(soundfile, correct)

d_a %<>%
    # add group columns to the dataframe
    inner_join(groups_a, by = "Subject") %>%
    # make columns to distinguish condition and word set (i.e. bag-tag)
    separate(cond, into = c("condition", "pres"), sep = "-", remove = FALSE) %>%
    separate(soundfile, into = c("word", "condition1", "pres1"), sep = "-", remove = FALSE) %>%
    mutate(item = paste(word, condition1, sep = "-")) %>%
    # get rid of columns we don't need
    dplyr::select(-pres, -word, -pres1, -condition1, -X1) %>%
    # add correct column
    left_join(key, by = "soundfile") %>%
    mutate(is_correct = ifelse(response == 'correct', 1, 0))

# adult sample size
counts <- d_a %>%
    group_by(Subject) %>%
    summarize(m = mean(keyRT)) %>%
    left_join(groups_a, by = "Subject")
count(counts, group)

#### Overall Accuracy ####
# Child d'
d_prime_c <- d_c
d_prime_c$hit <- ifelse(d_prime_c$condition == 'c' & d_prime_c$correct == 1, 1, 0)
d_prime_c$fa <- ifelse(d_prime_c$condition == 'i' & d_prime_c$correct == 0, 1, 0)
d_prime_c$miss <- ifelse(d_prime_c$condition == 'c' & d_prime_c$correct == 0, 1, 0)
d_prime_c$cr <- ifelse(d_prime_c$condition == 'i' & d_prime_c$correct == 1, 1, 0)

# Calculate d-prime
# Child
d_prime_c %<>%
    group_by(PartID) %>%
    summarize(hit = sum(hit, na.rm = TRUE),
              fa = sum(fa, na.rm = TRUE),
              miss = sum(miss, na.rm = TRUE),
              cr = sum(cr, na.rm = TRUE)) %>%
    filter_at(vars(hit, fa, miss, cr), any_vars(. != 0)) %>%
    na.omit()

indices_c <- neuropsychology::dprime(d_prime_c$hit, d_prime_c$miss, d_prime_c$fa, d_prime_c$cr)
d_prime_c$d <- indices_c$dprime

d_prime_c_all <- left_join(d_prime_c, select(d_c, PartID, condition), by = "PartID") %>%
    left_join(groups_c, by = "PartID") %>%
    mutate(group = ifelse(DD == 1, "Dys", "Typ")) %>%
    group_by(PartID, group) %>%
    summarize(d = mean(d)) %>%
    ungroup() %>%
    mutate(group = as.factor(group))

#### Child Analysis ####
# test for normal distribution
shapiro.test(d_prime_c_all$d)
t.test(d ~ group, data = d_prime_c_all) #run stats

#### Child Plot ####
# plot
ggplot(d_prime_c_all, aes(x = group, y = d, fill = group)) +
    geom_boxplot(alpha = 0.4) +
    stat_summary(fun.y = mean, geom = "point", shape = 20, size = 10, color = "red", fill = "red") +
    theme(legend.position = "none") +
    labs(title = "Child",
         x = "Group",
         y = "Accuracy (d')") +
    scale_fill_brewer(palette = "Set3")

#### Adult d' ####
d_prime_a <- d_a
d_prime_a$hit <- ifelse(d_prime_a$condition == 'cc' & d_prime_a$is_correct == 1, 1, 0)
d_prime_a$fa <- ifelse(d_prime_a$condition == 'ic' & d_prime_a$is_correct == 0, 1, 0)
d_prime_a$miss <- ifelse(d_prime_a$condition == 'cc' & d_prime_a$is_correct == 0, 1, 0)
d_prime_a$cr <- ifelse(d_prime_a$condition == 'ic' & d_prime_a$is_correct == 1, 1, 0)

d_prime_a %<>%
    group_by(Subject) %>%
    summarize(hit = sum(hit, na.rm = TRUE),
              fa = sum(fa, na.rm = TRUE),
              miss = sum(miss, na.rm = TRUE),
              cr = sum(cr, na.rm = TRUE)) %>%
    filter_at(vars(hit, fa, miss, cr), any_vars(. != 0)) %>%
    na.omit()
indices_a <- neuropsychology::dprime(d_prime_a$hit, d_prime_a$miss, d_prime_a$fa, d_prime_a$cr)
d_prime_a$d <- indices_a$dprime

d_prime_a_all <- left_join(d_prime_a, select(d_a, Subject, condition), by = "Subject") %>%
    left_join(groups_a, by = "Subject") %>%
    group_by(Subject, group) %>%
    summarize(d = mean(d))

counts <- d_prime_a_all %>%
    group_by(Subject) %>%
    summarize(m = mean(d)) %>%
    left_join(groups_a, by = "Subject")
count(counts, group)

#### Adult Analysis ####
shapiro.test(d_prime_a_all$d) # test for normal distribution
t.test(d ~ group, data = d_prime_a_all) #run stats

#### Adult Plot ####
ggplot(d_prime_a_all, aes(x = group, y = d, fill = group)) +
    geom_boxplot(alpha = 0.4) +
    stat_summary(fun.y = mean, geom = "point", shape = 20, size = 10, color = "red", fill = "red") +
    theme(legend.position = "none") +
    labs(title  =  "Adult",
        x  =  "Group",
        y  =  "Accuracy (d')") +
    scale_fill_brewer(palette = "Set3")

#### Combined Plots ####
d_combined <- bind_rows("Child" = d_prime_c_all,
                        "Adult" = rename(d_prime_a_all, PartID = Subject),
                        .id = "age")
d_combined$age <- factor(d_combined$age, levels = c("Child", "Adult"))

##### RTI

###Differences in RTI
#Prepare data
#Child
d_c2<-d_c%>%filter(condition!='n')%>%
    dplyr::group_by(PartID,condition) %>%
    dplyr::summarize(mean.rti = mean(keyRT, na.rm = TRUE))%>%
    ungroup() %>%
    left_join(groups_c, by = c("PartID"))
d_c2$group<-ifelse(d_c2$DD==1,"Dys","Typ")
d_c2$condition <- ifelse(d_c2$condition=="c","Cong","Incong")

#adult
d_a2<-d_a%>%dplyr::group_by(Subject,condition) %>%
    dplyr::summarize(mean.rti = mean(keyRT, na.rm = TRUE))%>%
    ungroup() %>%
    left_join(groups_a, by = c("Subject"))
d_a2$condition <- ifelse(d_a2$condition=="cc","Cong","Incong")
d_a2$group_cond<-paste(d_a2$group, "-", d_a2$condition)

##Adult
m1<-lmerTest::lmer(mean.rti~condition*group+(1|Subject), data=d_a2,REML=TRUE)
anova(m1)
lsmeans(m1, list(pairwise ~ group), adjust = "tukey")
lsmeans(m1, list(pairwise ~ condition), adjust = "tukey") 
eta_sq(m1)
##Child

d_c$DD<-as.factor(d_c$DD)
d_c$condition<-as.factor(d_c$condition)
m2<-lmer(mean.rti~condition*group+(1|PartID), data=d_c2,REML=TRUE)
anova(m2)
lsmeans(m2, list(pairwise ~ DD), adjust = "tukey")
lsmeans(m2, list(pairwise ~ condition), adjust = "tukey") #highest RT in Cong, Lowest in N
eta_sq(m1)
d_c2<-d_c%>%filter(condition!='n')

m3<-lmer(keyRT~condition*DD+(1|PartID), data=d_c,REML=TRUE)
anova(m3)
eta_sq(m3)

lsmeans(m3, list(pairwise ~ condition), adjust = "tukey") #highest RT in Cong, Lowest in N

##Plot  combined effects by rti


d_c2$group_cond<-paste(d_c2$group, "-", d_c2$condition)



p_c <- 
    d_c2 %>%
    dabest(group_cond,mean.rti,
           idx = list(c("Typ - Cong", "Dys - Cong"), 
                      c("Typ - Incong", "Dys - Incong")),
           paired = FALSE
    )
p1<-plot(p_c,
       #  rawplot.ylim = c(3, 7),
         rawplot.ylabel = "Child RT",color.column = group,
         effsize.ylabel = "Effect")



p_a <- 
    d_a2 %>%
    dabest(group_cond,mean.rti,
           idx = list(c("Typ - Cong", "Dys - Cong"), 
                      c("Typ - Incong", "Dys - Incong")),
           paired = FALSE
    )
p2<-plot(p_a,
       #  rawplot.ylim = c(3, 7),
         rawplot.ylabel = "Adult RT",color.column = group,
         effsize.ylabel = "Effect")
        effsize.ylim= c(0,0.4)
grid.arrange(p1,p2,nrow = 2)

#rti

#ggarrange(a,c,nrow=2, labels=c("A",
                           #"B",font.label = list(size = 16)))

#### Accuracy By Condition ####
# Child
d_c_all <- d_c %>%
    group_by(condition) %>%
    mutate(hit = ifelse(keyResp == 'z' & correct == 1, 1, 0),
           fa = ifelse(keyResp == 'm' & correct == 0, 1, 0)) %>%
    group_by(condition, PartID) %>%
    summarize(hit_rate = mean(hit, na.rm = TRUE),
              fa_rate = mean(fa, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(hit_rate = ifelse(hit_rate == 0, 0.00000001, hit_rate),
           fa_rate = ifelse(fa_rate == 0, 0.00000001, fa_rate))

# dprime function
dprime <- function(hit, fa) {
    qnorm(hit) - qnorm(fa)
}

# calculate dprimes (we do this before separating into groups)
d_c_all$dprime <- dprime(d_c_all$hit_rate, d_c_all$fa_rate)

# merge with group columns
d_c_all %<>% left_join(groups_c, by = "PartID") %>%
    # add condition to each hitfa dataframe
    mutate(condition = ifelse(condition == "c", "Congruent",
                            ifelse(condition == "i", "Incongruent",
                                   "Neutral")))

d_c_combined_hitfa <- d_c_all %>%
    mutate(group = ifelse(DD == 1, "Dys", "Typ"))



# create a dataframe with all dprime info to export
comp_c_all <- d_c_combined_hitfa %>%
    select(PartID, condition, dprime) %>%
    spread(condition, dprime) %>%
    rename(dprime_c = Congruent,
           dprime_i = Incongruent,
           dprime_n = Neutral) %>%
    mutate(phonresto_congruent_minus_incongruent = dprime_c - dprime_i)
d_c_final <- left_join(comp_c_all, groups_c, by = "PartID")

#### Outliers? ####
#Source: https://datascienceplus.com/identify-describe-plot-and-removing-the-outliers-from-the-dataset/
#outlier_values <- boxplot.stats(d_c_final$phonresto_congruent_minus_incongruent)$out  # outlier values.
#boxplot(d_c_final$phonresto_congruent_minus_incongruent, main="d'prime", boxwex=0.1)
#mtext(paste("Outliers: ", paste(outlier_values, collapse=", ")), cex=0.6)
#No outliers in the child data.

#### Analysis ####
shapiro.test(d_c_final$phonresto_congruent_minus_incongruent) # test for normal distribution
LME_model1 <- lmer(dprime ~ condition*group + (1|PartID), data = d_c_combined_hitfa, REML = FALSE)
summary(LME_model1)
anova(LME_model1)
lsmeans(LME_model1, list(pairwise ~ condition), adjust = "tukey")

#### Adult ####
d_a_all <- d_a %>%
    group_by(condition) %>%
    mutate(hit = ifelse(keyResp == 'z' & is_correct == 1, 1, 0),
           fa = ifelse(keyResp == 'm' & is_correct == 0, 1, 0)) %>%
    group_by(condition, Subject) %>%
    summarize(hit_rate = mean(hit, na.rm = TRUE),
              fa_rate = mean(fa, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(hit_rate = ifelse(hit_rate == 0, 0.00000001, hit_rate),
           fa_rate = ifelse(fa_rate == 0, 0.00000001, fa_rate))

# dprime function
dprime <- function(hit, fa) {
    qnorm(hit) - qnorm(fa)
}

# calculate dprimes (we do this before separating into groups)
d_a_all$dprime <- dprime(d_a_all$hit_rate, d_a_all$fa_rate)

d_a_all %<>%
    # merge with group columns
    inner_join(groups_a, by = "Subject") %>%
    # add condition
    mutate(condition = ifelse(condition == "cc", "Congruent", "Incongruent"))

d_a_combined_hitfa <- d_a_all


# create a dataframe with all dprime info to export
comp_a_all <- d_a_combined_hitfa %>%
    select(Subject, condition, dprime) %>%
    spread(condition, dprime) %>%
    rename(dprime_c = Congruent,
           dprime_i = Incongruent) %>%
    mutate(phonresto_congruent_minus_incongruent = dprime_c - dprime_i)

# Export dprimes to csv
d_a_final <- inner_join(comp_a_all, groups_a, by = "Subject")

#### Analysis ####
shapiro.test(d_a_final$phonresto_congruent_minus_incongruent)# test for normal distribution
LME_model2 <- lmer(dprime ~ condition*group + (1|Subject), data = d_a_combined_hitfa, REML = FALSE)
summary(LME_model2)
anova(LME_model2)
lsmeans(LME_model2, list(pairwise ~ condition), adjust = "tukey")

#### Combine and plot ####
pos_all <- bind_rows("Child" = pos_c,
                     "Adult" = pos_a,
                     .id = "age")
pos_all$age <- factor(pos_all$age, levels = c("Child", "Adult"))
combined_plot <- ggplot(pos_all, aes(x = condition, y = m_dprime, fill = group)) +
    theme(legend.position = "none") +
    geom_bar(stat = 'identity', position = position_dodge()) +
    scale_fill_brewer(palette = "Set1") +
    geom_errorbar(aes(ymin = lower, ymax = upper, width = .2), position = pd) +
    coord_cartesian(ylim = c(-5, 5)) +
    labs(title = "Phonemic Restoration Accuracy") +
    theme(plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 12)) +
    facet_wrap(~ age, scales = "free_x") # free scales to hide Neutral on adult x-axis
combined_plot

#### By age analysis ####
d_a_final %<>% rename(PartID = Subject)
d_c_final$group <- ifelse(d_c_final$Typical == 1, "Typ", "Dys")
d_c_final$group <- as.factor(d_c_final$group)

# combine the datasets
d <- bind_rows("child" = d_c_final,
               "adult" = d_a_final,
               .id = "age") %>%
    select(PartID, dprime_c, dprime_i, phonresto_congruent_minus_incongruent,
           group, age) %>%
    mutate_at(vars(PartID, group, age), as.factor)
model1 <- lm(phonresto_congruent_minus_incongruent ~ age*group, data = d)
anova(model1)
lsmeans(model1, list(pairwise ~ age), adjust = "tukey")
sjstats::eta_sq(model1)
