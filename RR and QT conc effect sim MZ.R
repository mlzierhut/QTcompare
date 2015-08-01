###############################################################
###############################################################
#                                                             #
#  Simulate drug induced heart rate effect on C-QTc analysis  #
#  Code from Matt Zierhut
#                                                             #
###############################################################
###############################################################

library(lme4)
library(magrittr)
library(plyr)
library(nlme)
library(scales) #for transparent plots (alpha())

n1 <- 60 #no. of subjects
n2 <- 7  #samples per subject

beta <- 0.4 #true beta (correction factor)
b1 <- 0.01  #true QTc_slope
b2 <- 0.55  #true RR_slope
b3 <- 0.0   #true Tx effect on beta

#create dummy data with random QTref and RRref values
df <- data.frame(
    ID=rep(seq(1,n1),each=n2),
    conc = rexp(n1*n2,0.01),
    eta1 = rep(rnorm(n1,0,20),each=n2),
    eta2 = rep(rnorm(n1,0,100),each=n2),
    RRref = rnorm(n1*n2,750,60),
    QTref = rnorm(n1*n2,410,10)
)

#Set baseline records
df$BL <- 0
df$BL[match(unique(df$ID), df$ID)] <- 1
df$conc[df$BL==1] <- 0

#Introduce Tx effect to QT and RR
df1 <- df %>% mutate(QT = (QTref+b1*conc+eta1) *(((RRref+b2*conc+eta2)/1000)^(beta+b3*I(1-BL)))) %>%
  mutate(RR = RRref+b2*conc+eta2)

#Estimate beta for two-stage approach
temp <- df1[df1$BL==1,]
temp$LNQT <- log(temp$QT)
temp$LNRR <- log(temp$RR)

x <- gls(LNQT~1+LNRR, data=temp)
summary(x)
(betaS <- x$coef[2])
betaB <- 0.5
betaF <- 1/3

#Calculate QTcX for S, F and B
df2 <- df1 %>% mutate(QTcS = (QT)/(((RR)/1000)^(betaS))) 
df2 <- df2 %>% mutate(QTcB = (QT)/(((RR)/1000)^(betaB))) 
df2 <- df2 %>% mutate(QTcF = (QT)/(((RR)/1000)^(betaF))) 


#Assess correction factors (QTcX vs RR)
#QT
plot(QT~RR,data=df2[df2$BL==1,],col=alpha(1,0.5),pch=16,
     xlim=range(df2$RR),ylim=range(df2$QT))
points(QT~RR,data=df2[df2$BL==0,],col=alpha(2,0.5),pch=16)

#QTcS
fitB <- gls(QTcS~RR,data=df2[df2$BL==1,])
fitT <- gls(QTcS~RR,data=df2[df2$BL==0,])
plot(QTcS~RR,data=df2[df2$BL==1,],col=alpha(1,0.5),pch=16,
     xlim=range(df2$RR),ylim=range(df2$QTcS))
points(QTcS~RR,data=df2[df2$BL==0,],col=alpha(2,0.5),pch=16)
abline(fitB,col=1)
abline(fitT,col=2)

#QTcF
fitB <- gls(QTcF~RR,data=df2[df2$BL==1,])
fitT <- gls(QTcF~RR,data=df2[df2$BL==0,])
plot(QTcF~RR,data=df2[df2$BL==1,],col=alpha(1,0.5),pch=16,
     xlim=range(df2$RR),ylim=range(df2$QTcF))
points(QTcF~RR,data=df2[df2$BL==0,],col=alpha(2,0.5),pch=16)
abline(fitB,col=1)
abline(fitT,col=2)

#QTcB
fitB <- gls(QTcB~RR,data=df2[df2$BL==1,])
fitT <- gls(QTcB~RR,data=df2[df2$BL==0,])
plot(QTcB~RR,data=df2[df2$BL==1,],col=alpha(1,0.5),pch=16,
     xlim=range(df2$RR),ylim=range(df2$QTcB))
points(QTcB~RR,data=df2[df2$BL==0,],col=alpha(2,0.5),pch=16)
abline(fitB,col=1)
abline(fitT,col=2)



#Fit QTcX (two-stage approach)
#QTcS
fit <- lme(QTcS ~ 1+conc,
           data=df2,
           random = ~1|ID,
           method="ML")
summary(fit)

#QTcF
fit <- lme(QTcF ~ 1+conc,
           data=df2,
           random = ~1|ID,
           method="ML")
summary(fit)

#QTcB
fit <- lme(QTcB ~ 1+conc,
           data=df2,
           random = ~1|ID,
           method="ML")
summary(fit)

# fit <- nlme(QTcS ~ QTc0,
#             data=df2,
#             fixed = QTc0~conc,
#             random = QTc0~1|ID,
#             start=c(400,0.05),
#             method="ML")
# summary(fit)



#One-stage approach
#slope is only on QTc
fit <- nlme(QT ~ (QTc0+conc*slp)*(RR/1000)^cor,
            data=df2,
            fixed = list(QTc0+slp~1,cor~1),
            random = QTc0~1|ID,
            start=c(400,0.01,0.4),
            method="ML")
summary(fit)

#slope is on beta and QTc
fit <- nlme(QT ~ (QTc0+conc*slp)*(RR/1000)^cor,
            data=df2,
            fixed = list(QTc0+slp~1,cor~I(conc/100)),
            random = QTc0~1|ID,
            start=c(400,0.01,0.4,0.1),
            method="ML")
summary(fit)

#slope on QTc, shift on beta
fit <- nlme(QT ~ (QTc0+conc*slp)*(RR/1000)^cor,
            data=df2,
            fixed = list(QTc0+slp~1,cor~I(1-BL)),
            random = QTc0+slp~1|ID,
            start=c(400,0.01,0.4,0.01),
            method="ML")
summary(fit)

#slope is only on beta (cor)
fit <- nlme(QT ~ (QTc0)*(RR/1000)^cor,
            data=df2,
            fixed = list(QTc0~1,cor~I(conc/100)),
            random = QTc0~1|ID,
            start=c(400,0.4,0.01),
            method="ML")
summary(fit)



#Compare beta to if used all data (not only baseline)
#test beta value
df2$LNQT <- log(df2$QT)
df2$LNRR <- log(df2$RR)
x <- gls(LNQT~1+LNRR, data=df2)
summary(x)
betaS

