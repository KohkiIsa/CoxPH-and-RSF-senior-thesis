version
install.packages("randomForestSRC")
install.packages("survival")
install.packages("dplyr")
install.packages("mice")
library("randomForestSRC")
library("survival")
library("dplyr")
library("mice")
?lung
str(lung) #全てnum
nrow(lung) #228
rm(lung)

lung
str(lung)
head(lung,n = 10)
nrow(lung) #228
df <- lung
colSums(is.na(df))
# sex - 1　必要あれば
df$sex <- df$sex - 1
head(df,n = 10)

# status - 1　必要あれば
df$status <- df$status - 1
head(df,n = 10)

str(df)
table(df$status)
summary(df)


# institution codeは分析対象には入れない
#欠損値を探す
which(is.na(df))

df|>
  filter(if_any(everything(), ~ is.na(.)))
colSums(is.na(df))
# phecog1個欠損 inst = 21
res <- df[ df$inst == 21 & !is.na(df$inst), ]
res
# Missing Completely At Randomなので消そう
df <- df %>%
  filter(!is.na(ph.ecog))
# 欠損値確認
df|>
  filter(if_any(everything(), ~ is.na(.)))
colSums(is.na(df))
# phkarno1個欠損 inst == 33
res1 <- df[df$inst == 33 & !is.na(df$inst),]
res1
# 欠損値1個のみなので消す（施設の性質上MCARか怪しいが）
df <- df %>%
  filter(!is.na(ph.karno))
# patkarno3個欠損 inst == 16 3 13
df|>
  filter(if_any(everything(), ~ is.na(df$pat.karno)))
res2 <- df[df$inst == 16 & !is.na(df$inst),]
res2
res3 <- df[df$inst == 3 & !is.na(df$inst),]
res3
res4 <- df[df$inst == 13 & !is.na(df$inst),]
res4
# MCARなので消す
df <- df %>%
  filter(!is.na(pat.karno))
colSums(is.na(df))

# meal.calとwt.lossについて埋める
# mice設計情報
ini <- mice(df, maxit = 0,printFlag = FALSE)

# methodベクトル作成
meth <- ini$method
pred <- ini$predictorMatrix

meth[] <- ""

meth["meal.cal"] <- "pmm"
meth["wt.loss"] <- "pmm"
meth

pred["meal.cal", ] <- 1
pred["meal.cal", "meal.cal"] <- 0
pred["wt.loss", ] <- 1
pred["wt.loss", "wt.loss"] <- 0

set.seed(1234)
imp <- mice(df,
            m = 20,          # 完成データを20セット
            maxit = 50,      # 反復回数
            method = meth,
            predictorMatrix = pred,
            printFlag = FALSE)

imp
df <- complete(imp, 1)
colSums(is.na(df))
nrow(df)
# 223
#欠損値残り
df|>
  filter(if_any(everything(), ~ is.na(.)))


# ph.ecogに関してダミー変数を作る
table(df$ph.ecog)
completely_ambulatory <- ifelse(df$ph.ecog == 1,1,0)
in_bed_under50 <- ifelse(df$ph.ecog == 2,1,0)
in_bed_over50 <- ifelse(df$ph.ecog == 3,1,0)
bedbound <- ifelse(df$ph.ecog == 4,1,0)
df <- data.frame (inst = df$inst,
                  time = df$time,
                  status = df$status,
                  age = df$age,
                  sex = df$sex,
                  ph.karno = df$ph.karno,
                  pat.karno = df$pat.karno,
                  meal.cal = df$meal.cal,
                  wt.loss = df$wt.loss,
                  completely_ambulatory = completely_ambulatory,
                  in_bed_under50 = in_bed_under50,
                  in_bed_over50 = in_bed_over50,
                  bedbound = bedbound)

# データの各カラムで何が何件あるか
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
# instもいれない
df <- data.frame (time = df$time,
                  status = df$status,
                  age = df$age,
                  sex = df$sex,
                  ph.karno = df$ph.karno,
                  pat.karno = df$pat.karno,
                  meal.cal = df$meal.cal,
                  wt.loss = df$wt.loss,
                  completely_ambulatory = completely_ambulatory,
                  in_bed_under50 = in_bed_under50,
                  in_bed_over50 = in_bed_over50)
# 欠損値を全て取り除き回しているのもある(パラメータや分割方法については触れられていない( A comparative study of forest methods for time-to-event data: variable selection and predictive performance))
# デフォルトのままでやってるのもある(Neutral Benchmarking of Survival Models in Health Sciences: Comparative Study of Classical and Machine Learning Techniques)
# RSFを行ったと考えられるがそれについての詳細な記載がない( Accelerated and Interpretable Oblique Random Survival Forests)
#欠損値確認
colSums(is.na(df))


#カプランマイヤー曲線
with(data = df, Surv(time,status))
Surv(df$time, df$status == 1)
survfit1 <- survfit(Surv(df$time,df$status)~1,data = df)
summary(survfit1)
plot(survfit1,las =1 ,xlab = "Survival Time(days)", ylab = "Survival",main = 'カプランマイヤー曲線',col=2:6)
colSums(df)
survfit2 <- survfit(Surv(time,status)~pat.karno,data = df)
legend_labels <- paste("pat.karno=" ,sort(unique(df$pat.karno)))
plot(survfit2,las =1 ,xlab = "Survival Time(days)", ylab = "Survival",main = 'pat.karno毎のカプランマイヤー曲線',col=2:6)
legend("topright", legend = legend_labels, lty=1,
       col=2:5)
survfit3 <- survfit(Surv(time,status)~ph.karno,data = df)
legend_labels <- paste("ph.karno =" ,sort(unique(df$ph.karno)))
plot(survfit3,las =1 ,xlab = "Survival Time(days)", ylab = "Survival",col=2:6)
legend("topright", legend = legend_labels, lty=1,
       col=2:5)
survfit4 <- survfit(Surv(time,status)~pat.karno,data = df)
legend_labels <- paste("pat.karno =" ,sort(unique(df$pat.karno)))
plot(survfit4,las =1 ,xlab = "Survival Time(days)", ylab = "Survival",col=2:5)
legend("topright", legend = legend_labels, lty=1,
       col=2:5)

#cox比例ハザード解析 対数尤度による変数選択をしてくれる関数制作
formula_one <- function(terms){
  right_hands <- if(length(terms) == 0) "1" else paste(terms, collapse = "+")
  as.formula(paste0("Surv(time,status) ~ ",right_hands))}
# p値出力関数
Likehood_p <- function(reduced,one_more){
  a <- anova(reduced,one_more,test = "Chisq")
  pcanditates <- c("Pr(>|Chi|)","Pr(>Chi)","P(>|Chi|)","P(>Chi)","Pr","P")
  pcol <- intersect(pcanditates,colnames(a))
  if (length(pcol) == 0 || nrow(a) < 2) return(NA_real_)
  p <- suppressWarnings(as.numeric(a[2, pcol[1]]))
  if (length(p) != 1 || !is.finite(p)) return(NA_real_)
  p
}

#前進選択ステップワイズ法で変数選択
forward_cox_LRT <- function(data,candidates = NULL,alpha = 0.05,ties = "efron"){
  if (is.null(candidates)){
    candidates <- setdiff(names(data),c("time","status"))
  }
  in_model <- character(0)
  improved <- TRUE
  path <- list()
  
  while(improved){
    improved <- FALSE
    best_var <- NULL
    best_p <- Inf
    
  fit0 <- coxph(formula_one(in_model),data = data,ties = ties)
  for (v in setdiff(candidates,in_model)){
    fit1 <- coxph(formula_one(c(in_model,v)), data = data,ties = ties)
    p <- Likehood_p(fit0,fit1)
    if (!is.na(p) && p < best_p){
      best_p <- p;best_var <- v
    }
    
  }
  if (!is.null(best_var)&&best_p < alpha){
    in_model<- c(in_model,best_var)
    improved <- TRUE
    path[[length(path) + 1]] <- list(step = "ADD", var = best_var, p = best_p)
    }
  }
  final_fit <- coxph(formula_one(in_model),data = data, ties = ties)
  list(model = final_fit,selected = in_model, history = path)
}

# 候補
cand <- setdiff(names(df),c("time","status"))
# 前進選択ステップワイズ法
fwd <- forward_cox_LRT(df,candidates = cand,alpha = 0.05)
summary(fwd$model)
print(fwd$selected)
print(fwd$history)

# cox検定＆プロット
fit_final <- coxph(formula_one(fwd$selected), data = df, ties = "efron")
cox_test <- cox.zph(fit_final)
print(cox_test)
plot(cox_test)   # 曲線が水平＆P>0.05が目安

null_fit <- coxph(Surv(time, status) ~ 1, data = df)
anova(null_fit, fwd$model, test = "Chisq")


#後退消去ステップワイズ法で変数選択
backward_cox_LRT <- function(data,candidates = NULL, alpha = 0.05,ties = "efron"){
  if (is.null(candidates)){
    candidates <- setdiff(names(data),c("time","status"))
  }
  in_model <- candidates
  improved <- TRUE
  path <- list()
  
  while(improved){
    improved <- FALSE
    worst_var <- NULL
    worst_p <- -Inf
    
    fullfit <- coxph(formula_one(in_model),data = data,ties = ties)
    for (v in in_model){
      red_terms <- setdiff(in_model, v)
      fit_remo <- coxph(formula_one(red_terms), data = data,ties = ties)
      p <- Likehood_p(fit_remo,fullfit)
      if (!is.na(p) && p > worst_p){
        worst_p <- p
        worst_var <- v
      }
      
    }
    if (!is.null(worst_var)&&worst_p >= alpha){
      in_model<- setdiff(in_model,worst_var)
      improved <- TRUE
      path[[length(path) + 1]] <- list(step = "DROP", var = worst_var, p = worst_p)
    }
  }
  final_fit2 <- coxph(formula_one(in_model),data = data, ties = ties)
  list(model = final_fit2,selected = in_model, history = path)
}

# 候補
cand <- setdiff(names(df),c("time","status"))
# 後退消去ステップワイズ法
bwd <- backward_cox_LRT(df,candidates = cand,alpha = 0.05,ties = "efron")
summary(bwd$model)
print(bwd$selected)
print(bwd$history)

# cox検定＆プロット
fit_final2 <- coxph(formula_one(bwd$selected), data = df, ties = "efron")
cox_test <- cox.zph(fit_final2)
print(cox_test)
plot(cox_test)   # 曲線が水平＆P>0.05が目安

null_fit <- coxph(Surv(time, status) ~ 1, data = df)
anova(null_fit, bwd$model, test = "Chisq")

#両方消去ステップワイズ法で変数選択
#cox比例ハザード解析
both_cox_LRT <- function(data,candidates = NULL, alpha_in = 0.05,alpha_out = 0.05,ties = "efron",start_terms = character(0),max_steps = 1000){
  if (is.null(candidates)){
    candidates <- setdiff(names(data),c("time","status"))
  }
  in_model <- start_terms
  improved <- TRUE
  path <- list()
  step_counter <- 0
  
  repeat{
    improved <- FALSE
    step_counter <- step_counter + 1
    if(step_counter > max_steps) break
    
    # fowardの部分
    fit0 <- coxph(formula_one(in_model),data = data,ties = ties)
    best_var <- NULL
    best_p <- Inf
    for (v in setdiff(candidates,in_model)){
      fit1 <- coxph(formula_one(c(in_model,v)), data = data,ties = ties)
      p <- Likehood_p(fit0,fit1)
      if (!is.na(p) && p < best_p){
        best_p <- p;best_var <- v
      }
      
    }
    if (!is.null(best_var)&&best_p < alpha_in){
      in_model<- c(in_model,best_var)
      improved <- TRUE
      path[[length(path) + 1]] <- list(step = "ADD", var = best_var, p = best_p)
    }
    # backwardの部分
    drop_more <- TRUE
    while(drop_more && length(in_model) > 0){
      fullfit <- coxph(formula_one(in_model),data = data,ties = ties)
      worst_var <- NULL
      worst_p <- -Inf
      for (v in in_model){
        red_terms <- setdiff(in_model, v)
        fit_remo <- if (length(red_terms) == 0){
          coxph(Surv(time,status) ~ 1,data = data,ties = ties)
          }else{
            coxph(formula_one(red_terms),data = data,ties = ties)
            }
        p <- Likehood_p(fit_remo,fullfit)
        if (!is.na(p) && p > worst_p){
          worst_p <- p
          worst_var <- v
        }
        
      }
      if (!is.null(worst_var)&&worst_p >= alpha_out){
        in_model<- setdiff(in_model,worst_var)
        improved <- TRUE
        path[[length(path) + 1]] <- list(step = "DROP", var = worst_var, p = worst_p)
      }else{
        drop_more <- FALSE
      }
    }
    
    if (!improved)break
  }
  final_fit2 <- coxph(formula_one(in_model),data = data, ties = ties)
  list(model = final_fit2,selected = in_model, history = path)
      
    }

# 候補
cand <- setdiff(names(df),c("time","status"))
# 両方選択ステップワイズ法
both <- both_cox_LRT(df,candidates = cand,alpha = 0.05,alpha_out = 0.05,ties = "efron")
summary(both$model)
print(both$selected)
print(both$history)

# cox検定＆プロット
fit_final3 <- coxph(formula_one(both$selected), data = df, ties = "efron")
cox_test <- cox.zph(fit_final3)
print(cox_test)
plot(cox_test)   # 曲線が水平＆P>0.05が目安

null_fit <- coxph(Surv(time, status) ~ 1, data = df)
anova(null_fit, both$model, test = "Chisq")


set.seed(1234)
#RSF logrank split
# クロスバリデーション
RSF1 <- rfsrc(Surv(time,status) ~ .,
              data = df,
              ntree = 1000,
              importance = "permute",
              nodesize = 15,
              mtry = 2,
              proximity = TRUE,
              splitrule = "logrank")
summary(RSF1)
print(RSF1)

cindex_oob1 <- 1 - tail(RSF1$err.rate, 1)
print(cindex_oob1)
vimp_tbl1 <- sort(RSF1$importance, decreasing = TRUE)
print(vimp_tbl1)

# conservation splitなかった

# logrankscore
RSF2 <- rfsrc(Surv(time,status) ~ .,
              data = df,
              ntree = 1000,
              importance = "permute",
              nodesize = 15,
              mtry = 2,
              proximity = TRUE,
              splitrule = "logrankscore")
summary(RSF2)
print(RSF2)
# こっちが数字いいかも

cindex_oob2 <- 1 - tail(RSF2$err.rate, 1)
print(cindex_oob2)
vimp_tbl2 <- sort(RSF2$importance, decreasing = TRUE)
print(vimp_tbl2)

# logrank tuning
set.seed(6)
tuning_RSF <- tune.rfsrc(
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

tuning_RSF$optimal
best_tuned <- tuning_RSF$rf
# 結果を改めて示す
summary(best_tuned)
print(best_tuned)

cindex_oob_t <- 1 - tail(best_tuned$err.rate, 1)
print(cindex_oob_t)
vimp_tbl_t <- sort(best_tuned$importance, decreasing = TRUE)
print(vimp_tbl_t)

# logrankscore tuning
set.seed(77)
tuning_RSF_lrc <- tune.rfsrc(
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

tuning_RSF_lrc$optimal
best_tuned_lrc <- tuning_RSF_lrc$rf
# 結果を改めて示す
summary(best_tuned_lrc)
print(best_tuned_lrc)

cindex_oob_lrc_t <- 1 - tail(best_tuned_lrc$err.rate, 1)
print(cindex_oob_lrc_t)
vimp_tbl_lrc_t <- sort(best_tuned_lrc$importance, decreasing = TRUE)
print(vimp_tbl_lrc_t)

# PH仮定（Schoenfeld残差）検定