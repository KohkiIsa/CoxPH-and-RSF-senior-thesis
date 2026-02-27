# senior thesis

## Overview
Rのsurvivalパッケージに存在する「北中部癌治療グループによる進行肺癌患者の生存データ」を用いてK-M曲線/Cox比例ハザード解析 / Random Survival Forest をまとめたポートフォリオです。
K-M曲線はカテゴリ特徴量について
Cox比例ハザード解析は前進選択ステップワイズ方法、後退消去ステップワイズ法を用いて
Random Survival Forest関しては分割方法と欠損値に対する処理毎に計算した。

## Approach
- **Features**：institute,time,status,age,sex,ph.ecog,ph.karno,pat.karno,meal.cal,wt.loss
- **Models**：Kaplan-Meier estimator,cox proportional hazards model,Random Survival Forest

## Result
