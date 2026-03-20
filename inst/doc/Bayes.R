## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ham)

## ----vlogo, echo=FALSE, out.width="30%", fig.align="center"-------------------
knitr::include_graphics("logo.png") 

## ----Post0--------------------------------------------------------------------
blos1 <- Bayes(x=losmcmc)

## ----plotdx1, fig.dim = c(6, 4.5)---------------------------------------------
plot(blos1, y="dxa", parameter="muOfY")

## ----plotdx2, fig.dim = c(6, 4.5)---------------------------------------------
plot(blos1, y="dxd", parameter="muOfY")

## ----plotdx3, fig.dim = c(6, 4.5)---------------------------------------------
plot(blos1, y="dxg", parameter="muOfY")

## ----plotdx4, fig.dim = c(6, 4.5)---------------------------------------------
plot(blos1, y="dxt", parameter="muOfY")

## ----Post1--------------------------------------------------------------------
blos1 <- Bayes(x=losmcmc, y="post", parameter="muOfY", newdata=TRUE)
print(blos1$Posterior.Summary)

## ----plotPost101, fig.dim = c(6, 4.5)-----------------------------------------
plot(x=blos1, y="post", parameter="muOfY", compare=4.5, rope=c(4,5), lcol= c("blue","red"),
bcol="goldenrod", HDItext=.3, main= "Summary of average LOS (muOfY)")

## ----plotPost1, fig.dim = c(6, 4.5)-------------------------------------------
plot(x=blos1, y="post", parameter=list("sigmaOfY", "muOfY" ),math="divide",
bcol="cyan", HDItext=.3, main= "Coefficient of Variation")

## ----plotPPC1, fig.dim = c(6, 4.5)--------------------------------------------
plot(x=blos1, y="check", type="n", data=hosprog, dv="los",
parameter=c("muOfY", "sigmaOfY"), breaks=30, cex.axis=1.3, lwd=3, xlab=NULL,
pline=20, vlim=c(-2, 20), xlim=c(-2, 20), add.legend="topright",
main="Length of Stay", cex.main=1.5, xpt=5, pcol="red", lcol="orange",
cex.legend=1, bcol="cyan")

## ----anovaPlot, echo=FALSE, out.width="60%", fig.align="center"---------------
knitr::include_graphics("taov.png") 

## ----Check1-------------------------------------------------------------------
bco2 <- Bayes(x=co2mcmc, y='mcmc', newdata=TRUE )

## ----plotPPC2, fig.dim = c(6, 4.5)--------------------------------------------
plot(x=bco2, y="check", type="ol", data=CO2, dv="uptake", iv="conc",
parameter=c("b0","b1"), add.data="al", cex.axis=1.3, lwd=1.5, pline=50,
vlim=c(50, 1100), xlim=c(0, 1100), ylim=c(0, 50), cex=2, cex.lab=2,
pcol="magenta", cex.main=2,cex.legend=1.2,  add.legend="topleft",
lcol="steelblue")              #vlim lets me extrapolate a little

## ----plotMulti1, fig.dim = c(6, 4.5)------------------------------------------
plot(x=co2multi, y="multi", level=2, aorder=FALSE,
subset= c("Qn2","Qn3","Qc3","Qc2","Mn3","Mn2","Mc2","Mc3"),
lcol="blue", pcol= c("red", "skyblue"), round.c=1, bcol="yellow",
xlim=c(-.1, 1), legend=NULL, add.legend="topright", lwd=3, cex.lab=1.2,
cex= 2, cex.main=1.25, cex.axis=.75, cex.legend=1.5, X.Lab=NULL)

## ----plotMulti2, fig.dim = c(6, 4.5)------------------------------------------
plot(x=co2multi, y="multi", level=3, aorder=FALSE, lcol="blue", pcol= c("green", "pink"),
round.c=1, bcol="lavender", xlim=c(-.1, 1), legend=NULL, add.legend="right", lwd=3,
cex.lab =1.2, cex= 2, cex.main=1.25, cex.axis=.75, cex.legend=1.5, X.Lab=NULL)

## ----Target1------------------------------------------------------------------
btarget1 <- Bayes(x=losmcmc, y="target", type="n", parameter=c("muOfY","sigmaOfY"),
newdata=TRUE, target=list(p=c(.35,.4,.45, .5, .55),  y=c(3,4))) 
print(btarget1$Target)

## ----plotTarget2, fig.dim = c(7.5, 4.5)---------------------------------------
plot(x=btarget1, y="target", type="n", data=hosprog, dv="los", breaks=30,
cex.axis=1.3, lwd=2, pline=20, vlim=c(-1, 12), xlim=c(-1, 10),
parameter=c("muOfY","sigmaOfY"), add.legend="right", main="Length of Stay",
cex.main=1.5, xpt=5, pcol="black", lcol="cyan", tgtcol="blue", bcol="orange",
cex.legend=1.25, cex.text = 1)

## ----plotTarget1, fig.dim = c(6, 4.5)-----------------------------------------
plot(x=btarget1, y="target", type="n", lcol="purple", tgtcol="blue", xlim=c(3.5, 5))

## ----R21----------------------------------------------------------------------
bR2 <- Bayes(x=co2mcmc, y='r2', data=CO2, iv="uptake", parameter=c("b0", "b1", "sigma"))
# R^2
print(bR2$R2.Summary$R2)
# Variance of predicted outcome
print(bR2$R2.Summary$Variance.Pred.Y)      
# Variance of residuals
print(bR2$R2.Summary$Variance.Residuals)
# A few predicted outcome values
print(head(bR2$R2.Summary$yPRED))          

