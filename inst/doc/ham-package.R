## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ham)

## ----vlogo, echo=FALSE, out.width="30%", fig.align="center"-------------------
knitr::include_graphics("logo.png") 

## ----group--------------------------------------------------------------------
gr1 <- group(x="program", y="los", z="month", data=hosprog, dist="t", cluster=TRUE)
print(gr1$Group.CI)

## ----group2-------------------------------------------------------------------
gr2 <- group(x="age", y="los", z="month", data=hosprog, dist="t", increment=3, rolling=6, quarts=TRUE)
print(gr2$Roll.CI$Rolling)

## ----plotGroup1, fig.dim = c(6, 4.5)------------------------------------------
plot(x=gr1, y="group", order="numeric", lwd=4, gcol= "blue", pcol="red", overall=TRUE,
     oband=TRUE, ocol="gray", tcol="green", tgt=4.5, cex=2, cex.axis=1, cex.lab=1.1,
     cex.text=2, cex.main=1.25, adj.alpha=.2)

## ----clustering---------------------------------------------------------------
print(gr1$Clustering)

## ----plotGroup2, fig.dim = c(6, 4.5)------------------------------------------
plot(x=gr1, y="time", lwd=4, gcol=c("red", "blue"), gband=TRUE, overall=TRUE, oband=TRUE,
     ocol="gray", tgt=4.5, tcol="green", tpline=5, tpcol="yellow", name=TRUE, cex.axis=1,
     cex.lab=1, cex.text=2, cex.main=1, adj.alpha=.3)

## ----ols----------------------------------------------------------------------
# Model summary
summary(assess(hp ~ mpg+wt, data=mtcars, regression="ols")$model)

## ----ols_interpret------------------------------------------------------------
# Model summary
interpret(assess(hp ~ mpg+wt, data=mtcars, regression="ols"))$model

## ----ols topcoding propensity-------------------------------------------------
m1 <- assess(formula=cost ~ month * program, data=hosprog, intervention = "program",
regression="ols", topcode=17150, propensity=c("female","age","risk"),
newdata=TRUE)

## ----ols topcoding propensity results-----------------------------------------
summary(m1$model)

## ----summary topcoding propensity---------------------------------------------
summary(m1$newdata[, c( "cost","top.cost", "pscore")])

## ----importance---------------------------------------------------------------
importance(m1$model)

## ----plotImportance, fig.dim = c(6, 4.5)--------------------------------------
#Consider using these graphical parameters
par(mar=c(4.2, 2, 3.5, 3))
par(oma = c(0, 0, 0, 3))
plot(importance(m1$model))

## ----did model 1--------------------------------------------------------------
dm1 <- assess(formula= los ~ ., data=hosprog, intervention = "program",
              int.time="month", treatment= 5, did="two")

## ----did model 1 results------------------------------------------------------
summary(dm1$DID)

## ----interpret did model 1----------------------------------------------------
interpret(dm1)$did

## ----plotDID1, fig.dim = c(6, 4.5)--------------------------------------------
plot(x=dm1, y="DID", add.legend="topleft", xlim=c(-.05, 1.05), ylim=c(2, 9), main="DID: Length of Stay", col=c("cyan","magenta"), lwd=7, cex=2, cex.axis=2, cex.lab=1.5,
     cex.main=3, arrow=TRUE, xshift=c(.045), cex.text=1.5, coefs=TRUE, round.c=2, cfact=T, conf.int=TRUE, adj.alpha=0.2 )

## ----did model 3--------------------------------------------------------------
dm3 <- assess(formula= rdm30 ~ ., data=hosprog, intervention = "program",
              int.time="month", treatment= 5, did="two")

## ----did model 3 results------------------------------------------------------
summary(dm3$DID)

## ----its model 1--------------------------------------------------------------
im11 <- assess(formula=los ~ ., data=hosp1, intervention = "program",
               int.time="month", interrupt= 5, its="one")

## ----its model 1 results------------------------------------------------------
summary(im11$ITS)

## ----its model 1 interpretations----------------------------------------------
interpret(im11)$its

## ----plotITS1, fig.dim = c(6, 4.5)--------------------------------------------
plot(im11, "ITS", add.legend="topright", xlim=c(-1, 14), ylim=c(2, 9), main="ITS: Intervention LOS", col="thistle", lwd=7, cex=3, cex.axis=2, cex.lab=1.5,  cex.main=3,
     arrow=TRUE, xshift=c(.25, .25), cex.text=1.5, coefs=TRUE, round.c=2, cfact=T, conf.int=TRUE, adj.alpha=0.2, pos.text= list("ITS.Time"=3, "Intercept"=4),
     cex.legend=1.25, add.means=TRUE )

## ----its model 2--------------------------------------------------------------
im12 <- assess(formula=los ~ ., data=hosp1, intervention = "program",
               int.time="month", interrupt= c(5, 9), its="one")

## ----its model 2 results------------------------------------------------------
summary(im12$ITS)

## ----its model 3--------------------------------------------------------------
im22 <- assess(formula=los ~ ., data=hosprog, intervention = "program",
               int.time="month", interrupt= c(5, 9), its="two")

## ----its model 3 results------------------------------------------------------
summary(im22$ITS)

## ----its model 3 interpretations----------------------------------------------
interpret(im22)$its

## ----plotAssess, fig.dim = c(6, 4.5)------------------------------------------
plot(im22, "ITS", add.legend="top", xlim=c(-.75, 13.1), ylim=c(2, 9), 
     main="ITS: Length of Stay", col=c("springgreen","thistle"), lwd=7, cex=2, cex.axis=2,
     cex.lab=1.5,  cex.main=3, cex.legend=1.25, arrow=TRUE, xshift=c(0, .5), cex.text=1,
     coefs=TRUE, round.c=1, pos.text= list("txp5"=3, "post9"=4), tcol="dodgerblue",
     conf.int=TRUE, adj.alpha=0.3, add.means=TRUE)

## ----its model 4--------------------------------------------------------------
id22 <- assess(formula=death30 ~ ., data=hosprog, intervention = "program", int.time="month", interrupt= c(5, 9), its="two")

## ----its model 4 results------------------------------------------------------
summary(id22$ITS)

## ----its model 5--------------------------------------------------------------
#Key interruption periods
key_time <- c(5, 14, 17, 29, 42, 59, 69, 73, 80,92)
im10 <- assess(formula=rate ~ ., data=unemployment, intervention = "usa",
      int.time="year", its="one", interrupt= key_time, newdata=TRUE)


## ----plotITS2, fig.dim = c(7.5, 5.25)-----------------------------------------
plot(im10, "ITS", add.means = TRUE, coefs=TRUE, conf.int=TRUE, 
     adj.alpha= .2, lwd=1.75, col="slategray", tcol= "orange", main="US unemployment rate",
     xlab="Years (1929-2024)", ylab= "Proportion of labor market", cex.main=2, 
     cex.axis = 1.25, cex.lab = 1.25, cex=2, cex.text= .75, pos.text=list("ITS.Time"=4,
"post42"=1,"txp42"=3,"txp92"=3), x.axis=unemployment$Year) 
for(i in 1:length(key_time)) {
  text(key_time[i], .22-(.01*i), cex=.85, labels = 
         paste0(unemployment[ key_time[i], "Year"], ": ", unemployment[ key_time[i], "event"]))
}

## ----cronbachs alpha example--------------------------------------------------
alpha(items=c("i1","i2","i3","i4","i5"), data=cas)

## ----cronbachs alpha interpret------------------------------------------------
interpret(alpha(items=c("i1","i2","i3","i4","i5"), data=cas))

