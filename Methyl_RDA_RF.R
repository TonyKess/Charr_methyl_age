library(tidyverse)
library(data.table)
library(vegan)
library(ranger)
setwd("/Users/KESST/Desktop/Projects/CharrMethylAge/")


Meta_Otolith_ID <- fread("CharrMethyl_MBB22_23_KL_19_Meta_OtolithAges_sizefilter.tsv")  %>%  mutate(Pop = str_replace(Dataset, "23","")) %>% 
  mutate(Pop = str_replace(Pop, "22","")) %>%
  mutate(Pop = str_replace(Pop, "19",""))

Meta_info <- colnames(Meta_Otolith_ID)

Mmatrix <- fread("CharrMBBKL.Mmatrix.tsv")
colnames(Mmatrix)[1] <- "Site"

Sites <- Mmatrix$Site
IDs <- data.frame(ID = colnames(Mmatrix %>%  select(-Site)) %>%  str_replace(., ".Me", ""))

Mmatrix_only <- Mmatrix %>%  select(-Site)

Mmatrix_flip <- data.frame(t(Mmatrix_only))
colnames(Mmatrix_flip) <- Sites

Mmatrix_Meta <- Mmatrix_flip %>%  
  mutate(ID = rownames(Mmatrix_flip))  %>% 
           mutate(ID = str_replace(ID, ".Me", "")) %>%   
           mutate(ID = str_replace(ID, "\\.", "-")) 

 
Mmatrix_Meta <- inner_join(Mmatrix_Meta, Meta_Otolith_ID) %>%  drop_na()
Mmatrix_Meta <- Mmatrix_Meta

Mmatrix_Meta %>% select(Pop, `Otolith Age`) %>%  group_by(Pop) %>%  summarise(maxage = max(`Otolith Age`))


test_set_MBB <- test_set %>% filter(Pop %in% "MBB")
test_set_KL <- test_set %>% filter(Pop %in% "KL")

test_set_loci <- test_set[,AgeOLsites]
test_set_meta <- test_set[,Meta_info]
test_set_loci_MBB <- test_set_MBB[,AgeOLsites]
test_set_meta_MBB <- test_set_MBB[,Meta_info]
test_set_loci_KL <- test_set_KL[,AgeOLsites]
test_set_meta_KL <- test_set_KL[,Meta_info]




train_set <-Mmatrix_Meta  %>%
  group_by(Pop) %>%
  sample_frac(size = 0.8) %>%  ungroup() 


train_set_MBB <-Mmatrix_Meta  %>%
  filter(Pop %in% "MBB") %>%
  sample_frac(size = 0.8) %>%  ungroup() 

train_set_KL <-Mmatrix_Meta  %>%
  filter(Pop %in% "KL") %>%
  sample_frac(size = 0.8) %>%  ungroup() 



train_set_loci <- train_set[,Sites]
train_set_meta <- train_set[,Meta_info]

train_set_loci_MBB <- train_set_MBB[,Sites]
train_set_meta_MBB <- train_set_MBB[,Meta_info]

train_set_loci_KL <- train_set_KL[,Sites]
train_set_meta_KL <- train_set_KL[,Meta_info]

rda.age.train <- rda(train_set_loci ~ train_set_meta$`Otolith Age` + Condition(train_set$Pop))
rda.age.train_MBB <- rda(train_set_loci_MBB ~ train_set_meta_MBB$`Otolith Age`)
rda.age.train_KL <- rda(train_set_loci_KL ~ train_set_meta_KL$`Otolith Age` )




train.rda.scores <- as.data.frame(unclass(scores(rda.age.train, display = "sites")))
train.rda.scores.meta <- bind_cols(train_set_meta, train.rda.scores)

train.rda.scores_MBB <- as.data.frame(unclass(scores(rda.age.train_MBB, display = "sites")))
train.rda.scores.meta_MBB <- bind_cols(train_set_meta_MBB, train.rda.scores_MBB)

train.rda.scores_KL <- as.data.frame(unclass(scores(rda.age.train_KL, display = "sites")))
train.rda.scores.meta_KL <- bind_cols(train_set_meta_KL, train.rda.scores_KL)


ggplot() + geom_point(data = train.rda.scores.meta, aes(x = RDA1, y = `Otolith Age`, colour = `Otolith Age`)) + 
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()

ggplot() + geom_point(data = train.rda.scores.meta, aes(x = RDA1, y = `Otolith Age`, colour = Pop)) + 
   theme_classic()


ggplot() + geom_point(data = train.rda.scores.meta_MBB, aes(x = RDA1, y = `Otolith Age`, colour = `Otolith Age`)) + 
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()

ggplot() + geom_point(data = train.rda.scores.meta_KL, aes(x = RDA1, y = `Otolith Age`, colour = `Otolith Age`)) + 
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()


Loadings <- data.frame(rda.age.train$CCA$v)
Hi5<- Loadings %>%  slice_max(RDA1, prop = 0.05)
Lo5 <- Loadings %>%  slice_min(RDA1, prop = 0.05)

Age_outliers <- bind_rows(Hi5, Lo5)
AgeOLsites <- rownames(Age_outliers)

Loadings_MBB <- data.frame(rda.age.train_MBB$CCA$v)
Hi5_MBB<- Loadings_MBB %>%  slice_max(RDA1, prop = 0.05)
Lo5_MBB <- Loadings_MBB %>%  slice_min(RDA1, prop = 0.05)

Age_outliers_MBB <- bind_rows(Hi5_MBB, Lo5_MBB)
AgeOLsites_MBB <- rownames(Age_outliers_MBB)


Loadings_KL <- data.frame(rda.age.train_KL$CCA$v)
Hi5_KL<- Loadings_KL %>%  slice_max(RDA1, prop = 0.05)
Lo5_KL <- Loadings_KL %>%  slice_min(RDA1, prop = 0.05)

Age_outliers_KL <- bind_rows(Hi5_KL, Lo5_KL)
AgeOLsites_KL <- rownames(Age_outliers_KL)

OLsites <- AgeOLsites_KL[AgeOLsites_KL %in% AgeOLsites_MBB]

cor.test( train.rda.scores.meta$RDA1, train.rda.scores.meta$`Otolith Age`)
cor.test( train.rda.scores.meta_MBB$RDA1, train.rda.scores.meta_MBB$`Otolith Age`)
cor.test( train.rda.scores.meta_KL$RDA1, train.rda.scores.meta_KL$`Otolith Age`)


# Get the remaining data as a test set
test_set <- anti_join(Mmatrix_Meta, train_set) 
test_set_MBB <- test_set %>% filter(Pop %in% "MBB")
test_set_KL <- test_set %>% filter(Pop %in% "KL")

test_set_loci <- test_set[,AgeOLsites]
test_set_meta <- test_set[,Meta_info]
test_set_loci_MBB <- test_set_MBB[,AgeOLsites_MBB]
test_set_meta_MBB <- test_set_MBB[,Meta_info]
test_set_loci_KL <- test_set_KL[,AgeOLsites_KL]
test_set_meta_KL <- test_set_KL[,Meta_info]



rda.age.test <- rda(test_set_loci ~ test_set_meta$`Otolith Age` + Condition(test_set$Dataset))
rda.age.test_MBB <- rda(test_set_loci_MBB ~ test_set_meta_MBB$`Otolith Age`)
rda.age.test_KL <- rda(test_set_loci_KL ~ test_set_meta_KL$`Otolith Age`)

test.rda.scores <- as.data.frame(unclass(scores(rda.age.test, display = "sites")))
test.rda.scores.meta <- bind_cols(test_set_meta, test.rda.scores)

test.rda.scores_MBB <- as.data.frame(unclass(scores(rda.age.test_MBB, display = "sites")))
test.rda.scores.meta_MBB <- bind_cols(test_set_meta_MBB, test.rda.scores_MBB)
test.rda.scores_KL <- as.data.frame(unclass(scores(rda.age.test_KL, display = "sites")))
test.rda.scores.meta_KL <- bind_cols(test_set_meta_KL, test.rda.scores_KL)



ggplot() + geom_point(data = test.rda.scores.meta, aes(x = RDA1, y = `Otolith Age`, colour = `Otolith Age`)) + 
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()

ggplot() + geom_point(data = test.rda.scores.meta_MBB, aes(x = RDA1, y = `Otolith Age`, colour = `Otolith Age`)) + 
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()

ggplot() + geom_point(data = test.rda.scores.meta_KL, aes(x = RDA1, y = `Otolith Age`, colour = `Otolith Age`)) + 
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()

cor.test( test.rda.scores.meta$RDA1, test.rda.scores.meta$`Otolith Age`)

cor.test( test.rda.scores.meta_MBB$RDA1, test.rda.scores.meta_MBB$`Otolith Age`)
cor.test( test.rda.scores.meta_KL$RDA1, test.rda.scores.meta_KL$`Otolith Age`)

rf_model_train_MBB <- ranger(x = as.matrix(train_set_loci_MBB[,OLsites]), 
                       y = as.matrix(train_set_meta_MBB$`Otolith Age`),
                       mtry = 100, num.trees = 1000)

rf_model_test_MBB <- predict(rf_model_train_MBB, data = as.matrix(test_set_loci_MBB[,OLsites]))


rf_model_train_MBB_predict <-  bind_cols(train_set_meta_MBB, predicted_age = rf_model_train_MBB$predictions)
predict_rf_MBB <- bind_cols(test_set_meta_MBB, predicted_age = rf_model_test_MBB$predictions)

cor.test(rf_model_train_MBB_predict$`Otolith Age`, rf_model_train_MBB_predict$predicted_age)

cor.test(predict_rf_MBB$`Otolith Age`, predict_rf_MBB$predicted_age)

sd(abs(predict_rf_MBB$predicted_age - predict_rf_MBB$`Otolith Age`))
sd(abs(predict_rf_MBB$predicted_age + predict_rf_MBB$`Otolith Age`))

ggplot() + geom_point(data = predict_rf_MBB, aes(x = predicted_age, y = `Otolith Age`, colour = `Otolith Age`)) +
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()
ggplot() + geom_point(data = rf_model_train_MBB_predict, aes(x = predicted_age, y = `Otolith Age`, colour = `Otolith Age`)) +
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()


rf_model_train_KL <- ranger(x = as.matrix(train_set_loci_KL[,OLsites]), 
                             y = as.matrix(train_set_meta_KL$`Otolith Age`),
                             mtry = 100, num.trees = 1000)

rf_model_test_KL <- predict(rf_model_train_KL, data = as.matrix(test_set_loci_KL[,OLsites]))


rf_model_train_KL_predict <-  bind_cols(train_set_meta_KL, predicted_age = rf_model_train_KL$predictions)
predict_rf_KL <- bind_cols(test_set_meta_KL, predicted_age = rf_model_test_KL$predictions)

cor.test(rf_model_train_KL_predict$`Otolith Age`, rf_model_train_KL_predict$predicted_age)

cor.test(predict_rf_KL$`Otolith Age`, predict_rf_KL$predicted_age)

sd(abs(predict_rf_KL$predicted_age - predict_rf_KL$`Otolith Age`))
sd(abs(predict_rf_KL$predicted_age + predict_rf_KL$`Otolith Age`))

ggplot() + geom_point(data = predict_rf_KL, aes(x = predicted_age, y = `Otolith Age`, colour = `Otolith Age`)) +
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()
ggplot() + geom_point(data = rf_model_train_KL_predict, aes(x = predicted_age, y = `Otolith Age`, colour = `Otolith Age`)) +
  scale_colour_gradient(low = "blue", high = "red") + theme_classic()
