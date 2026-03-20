## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ham)

## ----vlogo, echo=FALSE, out.width="30%", fig.align="center"-------------------
knitr::include_graphics("logo.png") 

## ----modClass1----------------------------------------------------------------
car_m1 <- assess(formula=vs ~ hp + am, data=mtcars, regression="logistic")
d1 <- decide(x=car_m1, threshold= -0.767)
print(d1$Model.Summary$Classification)

## ----modClass2----------------------------------------------------------------
print(d1$DCA)

## ----plotClass1, fig.dim = c(6, 4.5)------------------------------------------
plot(x=d1, y= "cl")

## ----plotClass2, fig.dim = c(7, 4.5)------------------------------------------
plot(x=d1, y= "cl", cex.lab=.75, bcol=c("cyan", "magenta"), add.legend="topleft", cex.legend=1.5)

## ----plotNB1, fig.dim = c(6, 4.5)---------------------------------------------
plot(x=d1, y= "nb")

## ----plotNB2, fig.dim = c(6, 4.5)---------------------------------------------
plot(x=d1, y= "nb", add.legend="topright", lwd=3, lcol=c("green", "slategray", "red"))

## ----plotIS1, fig.dim = c(6, 4.5)---------------------------------------------
plot(x=d1, y= "is")

