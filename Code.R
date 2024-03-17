
############################################
# SCHILD Nathan & SOMMER Zoé               #
# Master APE                               #
# Courbe de Phillips                       #
############################################


# ---Installation des packages et chargement de la base de données---
library(dplyr)
library(stargazer)
library(urca)
library(zoo)
library(xts)
library(fBasics)
library(lmtest)
library(dynlm)
library(astsa)
library(AER)
library(readxl)
library(ggplot2)
library(RefManageR)

data <- read_excel("data.xlsx")
head(data)


# ---Traitement de la base de données---

dataFR <- data %>% filter(country=="FRA") # Filtrage de la base 
dataUK <- data %>% filter(country=="UK")

# On fait en sorte que la colonne time de notre dataframe soit reconnue comme étant une date
dataFR$time <- as.Date(paste0("01-", 
                              gsub("Q([0-9]+)-([0-9]+)", "\\1-\\2", 
                                   dataFR$time)), 
                       format = "%d-%m-%Y")
dataUK$time <- as.Date(paste0("01-", 
                              gsub("Q([0-9]+)-([0-9]+)", 
                                   "\\1-\\2", 
                                   dataUK$time)), 
                       format = "%d-%m-%Y")

# On met la base de données sous forme de séries temporelles
dataFR_ts <- ts(dataFR[, c("unemp", "inf")], start = start(dataFR$time), frequency = 4)
dataUK_ts <- ts(dataUK[, c("unemp", "inf")], start = start(dataUK$time), frequency = 4)


# ---Statistiques descriptives et analyse graphique---

statsFR <- dataFR %>% select(unemp,inf)
basicStats(statsFR)
statsUK <- dataUK %>% select(unemp,inf)
basicStats(statsUK) # Statistiques descriptives 

par(mfrow=c(1,2))

# Création des graphiques présentant la relation inflation-chômage 
ggplot(dataFR, aes(x = unemp, y = inf)) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "Relation inflation-chômage en France",
       x = "Taux de chômage",
       y = "Inflation") +
  theme_minimal() 


ggplot(dataUK, aes(x = unemp, y = inf)) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "Relation inflation-chômage au Royaume-Uni",
       x = "Taux de chômage",
       y = "Inflation") +
  theme_minimal() 

par(mfrow=c(1,1))


# ---Les estimations économétriques---

# i) Le cas français

MCOfr <- lm(inf~unemp,data=dataFR)
modelfr <- dynlm(inf~unemp+L(inf,1)+L(unemp,1),data=dataFR_ts)
rob_seFR <- list(sqrt(diag(vcovHC(MCOfr, type = "HC1"))),
                 sqrt(diag(vcovHC(modelfr, type = "HC1")))) # Ecarts-types robustes

stargazer(list(MCOfr, modelfr),
          title = "Résultats des estimations pour la France",
          align = T, 
          type = "latex",
          no.space = F,
          se=rob_seFR)

# ii) Le cas britannique

MCOuk <- lm(inf~unemp,data=dataUK)
modeluk <- dynlm(inf~L(inf,1)+unemp+L(unemp,1),data=dataUK_ts)
rob_seUK <- list(sqrt(diag(vcovHC(MCOuk, type = "HC1"))),
                 sqrt(diag(vcovHC(modeluk, type = "HC1")))) # Ecarts-types robustes

stargazer(list(MCOuk, modeluk),
          title = "Résultat des estimations pour le Royaume-Uni",
          align = T, 
          type = "latex",
          no.space = F,
          se=rob_seUK)


# ---Les tests de stationnarité---

# ---Le cas français---

ts_inflation <- ts(dataFR$inf)
ts_unemp <- ts(dataFR$unemp)

result_df_no_trend_drift_inflation <- ur.df(ts_inflation,type="none")
result_df_with_drift_inflation <- ur.df(ts_inflation,type="drift")
result_df_with_trend_inflation <- ur.df(ts_inflation,type="trend") # TEST DF

summary(result_df_no_trend_drift_inflation)
summary(result_df_with_drift_inflation)
summary(result_df_with_trend_inflation)

result_df_no_trend_drift_unemp <- ur.df(ts_unemp,type="none")
result_df_with_drift_unemp <- ur.df(ts_unemp,type="drift")
result_df_with_trend_unemp <- ur.df(ts_unemp,type="trend") # TEST DF

summary(result_df_no_trend_drift_unemp)
summary(result_df_with_drift_unemp)
summary(result_df_with_trend_unemp)

summary(ur.kpss(ts_inflation,type="mu"))
summary(ur.kpss(ts_unemp,type="mu")) # TEST KPSS

result_pp_inflation <- ur.pp(ts_inflation, type = "Z-tau")
result_pp_unemp <- ur.pp(ts_unemp, type = "Z-tau") # TEST PP
summary(result_pp_inflation)
summary(result_pp_unemp)

par(mfrow=c(1,2))

acf(ts_inflation) # Fonction d'autocorrélation
acf(ts_unemp)

par(mfrow=c(1,1))


# ---Le cas anglais---

ts_inflation2 <- ts(dataUK$inf)
ts_unemp2 <- ts(dataUK$unemp)

result_df_no_trend_drift_inflation2 <- ur.df(ts_inflation2,type="none")
result_df_with_drift_inflation2 <- ur.df(ts_inflation2,type="drift")
result_df_with_trend_inflation2 <- ur.df(ts_inflation2,type="trend") # TEST DF

summary(result_df_no_trend_drift_inflation2)
summary(result_df_with_drift_inflation2)
summary(result_df_with_trend_inflation2)


result_df_no_trend_drift_unemp2 <- ur.df(ts_unemp2,type="none")
result_df_with_drift_unemp2 <- ur.df(ts_unemp2,type="drift")
result_df_with_trend_unemp2 <- ur.df(ts_unemp2,type="trend")

summary(result_df_no_trend_drift_unemp2)
summary(result_df_with_drift_unemp2)
summary(result_df_with_trend_unemp2)

summary(ur.kpss(ts_inflation2,type="mu"))
summary(ur.kpss(ts_unemp2,type="mu")) # TEST KPSS

result_pp_inflation2 <- ur.pp(ts_inflation2, type = "Z-tau")
result_pp_unemp2 <- ur.pp(ts_unemp2, type = "Z-tau") # TEST PP
summary(result_pp_inflation2)
summary(result_pp_unemp2)

par(mfrow=c(1,2))

acf(ts_inflation2) # Fonction d'autocorrélation
acf(ts_unemp2)

par(mfrow=c(1,1))


# ---Les tests sur les différences---

# --- Le cas français---

summary(ur.df(diff(ts_unemp),type="trend"))
summary(ur.df(diff(ts_unemp),type="none"))
summary(ur.df(diff(ts_unemp),type="drift"))

summary(ur.df(diff(ts_inflation),type="trend"))
summary(ur.df(diff(ts_inflation),type="none"))
summary(ur.df(diff(ts_inflation),type="drift"))

summary(ur.kpss(diff(ts_inflation),type="mu"))
summary(ur.kpss(diff(ts_unemp),type="mu"))

summary(ur.pp(diff(ts_inflation), type = "Z-tau"))
summary(ur.pp(diff(ts_unemp), type = "Z-tau"))

# --- Le cas anglais---

summary(ur.df(diff(ts_unemp2),type="trend"))
summary(ur.df(diff(ts_unemp2),type="none"))
summary(ur.df(diff(ts_unemp2),type="drift"))

summary(ur.df(diff(ts_inflation2),type="trend"))
summary(ur.df(diff(ts_inflation2),type="none"))
summary(ur.df(diff(ts_inflation2),type="drift"))

summary(ur.kpss(diff(ts_inflation2),type="mu"))
summary(ur.kpss(diff(ts_unemp2),type="mu"))

summary(ur.pp(diff(ts_inflation2), type = "Z-tau"))
summary(ur.pp(diff(ts_unemp2), type = "Z-tau"))


# ---Test de cointégration---

summary(ur.df(residuals(MCOuk),type="none",lags=3))
summary(ur.df(residuals(MCOuk),type="drift",lags=3))
summary(ur.df(residuals(MCOuk),type="trend",lags=3))

summary(ur.kpss(residuals(MCOuk),type="mu"))
summary(ur.pp(residuals(MCOuk),type = "Z-tau"))

summary(ur.df(residuals(MCOfr),type="none",lags=3))
summary(ur.df(residuals(MCOfr),type="drift",lags=3))
summary(ur.df(residuals(MCOfr),type="trend",lags=3))

summary(ur.kpss(residuals(MCOfr),type="mu"))
summary(ur.pp(residuals(MCOfr),type = "Z-tau"))

# ---Modèle à correction d'erreurs--- 

dataFR$resFR<- residuals(MCOfr)
data_tsFR <- ts(dataFR[, c("inf", "unemp","resFR")], frequency = 1)

MCEFR <- dynlm(diff(inf)~ diff(unemp)+L(resFR,1),data=data_tsFR) # Modèle à correction d'erreurs

rob_seFR2 <- list(sqrt(diag(vcovHC(MCEFR, type = "HC1")))) # Ecarts-types robustes

stargazer(MCEFR,title = "Résultats des estimations pour la France",
          align = T, type = "latex",
          no.space = F,se=rob_seFR2)




