###############################################################
###############################################################
#                                                             #
#  Simulate drug induced heart rate effect on C-QTc analysis  # 
#     (parameters used for simulation are derived from        #
#             crizotinib two-stage modeling)                  # 
#                                                             #
###############################################################
###############################################################

library(lme4)
library(magrittr)
library(plyr)
library(nlme)
library(scales) #for transparent plots (alpha())

n1 <- 2000#60 #no. of subjects
n2 <- 5#7  #samples per subject

beta <- 0.4   #true beta (correction factor)
b1 <- 0.05#0.05    #true QTc_slope
b2 <- 0.0#0.5    #true RR_slope
b3 <- -0.00#-0.05   #true Tx effect on beta

#create dummy data with random QTref and RRref values
df <- data.frame(
    ID=rep(seq(1,n1),each=n2),
    conc = rexp(n1*n2,0.01),
    eta1 = rep(rnorm(n1,0,20),each=n2),
    eta2 = rep(rnorm(n1,0,100),each=n2),
    RRref = rnorm(n1*n2,750,60),
    QTref = rnorm(n1*n2,420,10)
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

temp <- aggregate(QT ~ ID,data=df2[df2$BL==1,],mean,na.action=na.omit)
names(temp)[2] <- "BQT"
df2 <- merge(df2[,names(df2)!="BQT"],temp,by="ID", all.x=TRUE)

temp <- aggregate(RR ~ ID,data=df2[df2$BL==1,],mean,na.action=na.omit)
names(temp)[2] <- "BRR"
df2 <- merge(df2[,names(df2)!="BRR"],temp,by="ID", all.x=TRUE)

temp <- aggregate(QTcS ~ ID,data=df2[df2$BL==1,],mean,na.action=na.omit)
names(temp)[2] <- "BQTcS"
df2 <- merge(df2[,names(df2)!="BQTcS"],temp,by="ID", all.x=TRUE)

df2$DQTcS <- df2$QTcS - df2$BQTcS

median(df2$QTcS[df$BL==1])
median(df2$QTcF[df$BL==1])
median(df2$QTcB[df$BL==1])

#Assess correction factors (QTcX vs RR)
#QT
plot(QT~RR,data=df2,col=alpha(2-BL,0.5),pch=16,
     xlim=range(df2$RR),ylim=range(df2$QT))

#QTcS
fitB <- gls(QTcS~RR,data=df2[df2$BL==1,])
fitT <- gls(QTcS~RR,data=df2[df2$BL==0,])
plot(QTcS~RR,data=df2,col=alpha(2-BL,0.5),pch=16,
     xlim=range(df2$RR),ylim=range(df2$QTcS))
abline(fitB,col=1)
abline(fitT,col=2)

#QTcF
fitB <- gls(QTcF~RR,data=df2[df2$BL==1,])
fitT <- gls(QTcF~RR,data=df2[df2$BL==0,])
plot(QTcF~RR,data=df2,col=alpha(2-BL,0.5),pch=16,
     xlim=range(df2$RR),ylim=range(df2$QTcS))
abline(fitB,col=1)
abline(fitT,col=2)



#Fit QTcX (two-stage approach)
#QTcF
fit <- lme(QTcF ~ 1+conc,
           data=df2,
           random = ~1|ID,
           method="ML")
summary(fit)

#QTcS
fit <- lme(QTcS ~ 1+conc,
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

#DQTcS
# fit <- lme(DQTcS ~ conc-1,
#            data=df2,
#            random = ~conc-1|ID,
#            method="ML")
# summary(fit)

fit <- gls(DQTcS ~ conc-1,
           data=df2[df2$BL==0,],
           method="ML")
summary(fit)





#One-stage approach
#only slope on QTc (use BQT and BRR)
# fit <- nlme(QT ~ ((BQT/(BRR/1000)^cor)+conc*slp)*(RR/1000)^cor,
#             data=df2[df2$BL==0,],
#             fixed = list(slp~1,cor~1),
#             random = cor~1|ID,
#             start=c(0.01,0.4),
#             method="ML")
# summary(fit)
#No eta
fit <- nls(QT ~ ((BQT/(BRR/1000)^cor)+conc*slp)*(RR/1000)^cor,
            data=df2[df2$BL==0,],
            start=list(slp=0.01,cor=0.4))
summary(fit)


#only slope on QTc (use BQT and BRR), with beta shift
# fit <- nlme(QT ~ ((BQT/(BRR/1000)^cor)+conc*slp)*(RR/1000)^(cor+beff*I(1-BL)),
#             data=df2[df2$BL==0,],
#             fixed = list(slp+cor+beff~1),
#             random = cor~1|ID,
#             start=c(0.01,0.4,0.01),
#             method="ML")
# summary(fit)
#No eta
fit <- nls(QT ~ ((BQT/(BRR/1000)^cor)+conc*slp)*(RR/1000)^(cor+beff*I(1-BL)),
           data=df2[df2$BL==0,],
           start=list(slp=0.01,cor=0.4,beff=0.01))
summary(fit)


#only slope on QTc
fit <- nlme(QT ~ (QTc0+conc*slp)*(RR/1000)^cor,
            data=df2,
            fixed = list(QTc0+slp~1,cor~1),
            random = QTc0~1|ID,
            start=c(400,0.01,0.4),
            method="ML")
summary(fit)

#slope on QTc, shift on beta
fit <- nlme(QT ~ (QTc0+conc*slp)*(RR/1000)^cor,
            data=df2,
            fixed = list(QTc0+slp~1,cor~I(1-BL)),
            random = QTc0~1|ID,
            start=c(400,0.01,0.4,0.01),
            method="ML")
summary(fit)

#only shift on beta (cor)
fit <- nlme(QT ~ (QTc0)*(RR/1000)^cor,
            data=df2,
            fixed = list(QTc0~1,cor~I(1-BL)),
            random = QTc0~1|ID,
            start=c(400,0.4,0.01),
            method="ML")
summary(fit)









#Miscellaneous

QTcB <- 432
QTcS <- 422
QTcF <- 416
beta <- 0.4

RR <- seq(200,1800,length.out=1000)
QTS <- QTcS*(RR/1000)^beta
QTF <- QTcF*(RR/1000)^(1/3)
QTB <- QTcB*(RR/1000)^(1/2)

plot(QTS~RR,type="l")
lines(QTF~RR,col=2)
lines(QTB~RR,col=4)

QTcS2 <- 432
QTS2 <- QTcS2*(RR/1000)^beta
lines(QTS2~RR,col=3)


beta0=0.4
beta=1/3

#RR constant --> magnitude of QT change (QT2) doesn't matter (d/dQT2 = 0)
100*(((260/(400/1000)^beta - 280/(400/1000)^beta)/
        (260/(400/1000)^beta0 - 280/(400/1000)^beta0))-1)

100*(((500/(1500/1000)^beta - 530/(1500/1000)^beta)/
        (500/(1500/1000)^beta0 - 530/(1500/1000)^beta0))-1)

#No bigger difference than 10%



#RR changes --> correction factor is critical!
100*(((300/(500/1000)^beta - 350/(735/1000)^beta)/
        (300/(500/1000)^beta0 - 350/(735/1000)^beta0))-1)

