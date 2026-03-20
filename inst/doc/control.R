## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ham)

## ----vlogo, echo=FALSE, out.width="30%", fig.align="center"-------------------
knitr::include_graphics("logo.png") 

## ----xbar1--------------------------------------------------------------------
spc_x <- control(x="los", time="month", data=hosprog, type="x", n.equal=TRUE)
print(spc_x)

## ----chartX1, fig.dim = c(6, 4.5)---------------------------------------------
plot(spc_x)

## ----pchart1------------------------------------------------------------------
spc_p <- control(x="rdm30", time="month", data=hosprog, type="p", n.equal=FALSE)

## ----chartp1, fig.dim = c(6, 4.5)---------------------------------------------
plot(spc_p, tgt=c(0,.25), tgtcol="green", ylim=c(0,0.4), tpline=c(4,8),
tpcol= c("yellow","black"))

## ----uchart1------------------------------------------------------------------
spc_u <- control(x="HAI", y="PatientDays", time="Month", data=infections,
type="u", n.equal=FALSE, intervention=22)

## ----chartU1, fig.dim = c(6, 4.5)---------------------------------------------
plot(spc_u, main="u-Chart: HAI per 1,000 Patient Days Pre/Post Intervention",
col=c("green","dodgerblue"), trend=TRUE, trcol="red", x.axis=c((1:41+12)), round.c=1,
y.axis=seq(min(spc_u$HAI)*1000, max(spc_u$HAI)*1000, length.out=nrow(spc_u)),
xlab="Months (starting at year 2)", icol="gray", lwd=2, cex=2,
cex.axis=1.1, cex.main=1.25, cex.text=1.25)

