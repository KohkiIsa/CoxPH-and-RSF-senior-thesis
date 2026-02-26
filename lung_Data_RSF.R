version
install.packages("randomForestSRC")
install.packages("survival")
install.packages("dplyr")
library("randomForestSRC")
library("survival")
library("dplyr")
?lung
str(lung) #全てnum
nrow(lung) #228
rm(lung)

lung
str(lung)
head(lung,n = 10)
nrow(lung) #228
df <- lung

# sex -1
df$sex <- df$sex - 1
head(df,n = 10)

# status -1
df$status <- df$status - 1
head(df,n = 10)

# institution codeは分析対象には入れない？
table(df$inst)#他の変数も作る
table(df$time)
table(df$status)
table(df$age)
table(df$sex)
table(df$ph.karno)
table(df$pat.karno)
table(df$meal.cal)
table(df$wt.loss)
table(df$completely_ambulatory)
table(df$in_bed_under50)
table(df$in_bed_over50)
table(df$bedbound)
# bedbound 0件で草　変数消そう
df <- data.frame (time = df$time,
                  status = df$status,
                  age = df$age,
                  sex = df$sex,
                  ph.karno = df$ph.karno,
                  pat.karno = df$pat.karno,
                  meal.cal = df$meal.cal,
                  wt.loss = df$wt.loss)

set.seed(234)
#RSF logrank split
# クロスバリデーション
RSF_1_1 <- rfsrc(Surv(time,status) ~ .,
                  data = df,
                  na.action = "na.impute",
                  ntree = 1000,
                  importance = "permute",
                  nodesize = 15,
                  mtry = 2,
                  proximity = TRUE,#近接度行列の計算
                  splitrule = "logrank")
summary(RSF_1_1)
print(RSF_1_1)

cindex_oob1_1 <- 1 - tail(RSF_1_1$err.rate, 1)
print(cindex_oob1_1)
vimp_tbl1_1 <- sort(RSF_1_1$importance, decreasing = TRUE)
print(vimp_tbl1_1)
# conservation splitなかった

# logrankscore
RSF_1_2 <- rfsrc(Surv(time,status) ~ .,
                  data = df,
                  na.action = "na.impute",
                  ntree = 1000,
                  importance = "permute",
                  nodesize = 15,
                  mtry = 2,
                  proximity = TRUE,
                  splitrule = "logrankscore")
summary(RSF_1_2)
print(RSF_1_2)

cindex_oob1_2 <- 1 - tail(RSF_1_2$err.rate, 1)
print(cindex_oob1_2)
vimp_tbl1_2 <- sort(RSF_1_2$importance, decreasing = TRUE)
print(vimp_tbl1_2)

# tuning
set.seed(12)
tuning_RSF_1 <- tune.rfsrc(
  Surv(time,status) ~ .,
  data = df,
  method = "grid",
  ntreeTry = 1000,
  nodesizeTry = c(1,5,10,15,20,25,30,50),
  mtryStart = c(2,4,6),
  proximity = TRUE,
  splitrule = "logrank",
  importance = "permute"
)

tuning_RSF_1$optimal
best_tuned_1 <- tuning_RSF_1$rf
# 結果を改めて示す
summary(best_tuned_1)
print(best_tuned_1)

cindex_oob_t_1 <- 1 - tail(best_tuned_1$err.rate, 1)
print(cindex_oob_t_1)
vimp_tbl_t_1 <- sort(best_tuned_1$importance, decreasing = TRUE)
print(vimp_tbl_t_1)

# logrankscore tuning
set.seed(77)
tuning_RSF_lrc_1 <- tune.rfsrc(
  Surv(time,status) ~ .,
  data = df,
  method = "grid",
  ntreeTry = 1000,
  nodesizeTry = c(1,5,10,15,20,25,30,50),
  mtryStart = c(2,4,6),
  proximity = TRUE,
  splitrule = "logrankscore",
  importance = "permute"
)

tuning_RSF_lrc_1$optimal
best_tuned_lrc_1 <- tuning_RSF_lrc_1$rf
# 結果を改めて示す
summary(best_tuned_lrc_1)
print(best_tuned_lrc_1)

cindex_oob_lrc_t_1 <- 1 - tail(best_tuned_lrc_1$err.rate, 1)
print(cindex_oob_lrc_t_1)
vimp_tbl_lrc_t_1 <- sort(best_tuned_lrc_1$importance, decreasing = TRUE)
print(vimp_tbl_lrc_t_1)