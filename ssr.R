library(gsDesign)
library(tibble)
library(ggplot2)
library(scales)
library(gt)

### what is exactly the combination test, and how the CHW method was specified in the gsDesign::ssrCP function?

### calculate conditional power based on alpha, stage and stage 2 sample sizes, assumed treatment effects, and interim z

alpha <- .01 
hr <- .6 
n.fix <- nEvents(hr=hr, alpha = alpha) 

timing <- .6 
x <- gsDesign(k = 2, n.fix = n.fix, alpha = alpha, test.type = 1, sfu = sfHSD, sfupar = -20, timing = timing, delta1 = log(hr))
gsBoundSummary(x) %>% gt %>% tab_header(title = "Time-to-event group sequential design") %>% cols_align("left")

# note that this is a group sequential design with the same sample size and power as of the fixed design
# this is possible due to stringent IA bound

hrpostIA = seq(.4, 1, .05) # effect sizes for which you wish to compute conditional power
powr <- condPower(x = x, z1 = 1, n2 = x$n.I[2]-x$n.I[1], theta = log(hrpostIA)/x$delta1*x$theta[2])
qplot(hrpostIA, powr, geom="line",xlab = "HR post IA",ylab="Conditional power", main ="Conditional power as a function of assumed HR")

### ssrCP 

alpha <- .025
beta <- .1

delta0 <- 0
delta1 <- 1

timing <- .5
sfu <- sfHSD # upper spending function
sfupar <- -12

maxinflation <- 2 # maximum sample size inflation
overrun <- 25 # assumed enrollment overrun at IA

z <- seq(0, 4, .025) # interim z-values for plotting
n.fix <- 100 # fixed design sample size

cpadj <- c(.3, .9) # conditional power interval where sample size is to be adjusted
betastar <- beta # targeted Type II error when adapting sample size
z2 <- z2NC # combination test (built-in options are: z2Z, z2NC, z2Fisher)

# use the above parameters to generate a 2-stage group sequential design
x <- gsDesign(k = 2, n.fix = n.fix, timing = timing, sfu = sfu, sfupar = sfupar, alpha = alpha, beta = beta, delta0 = delta0, delta1 = delta1)
gsBoundSummary(x) %>% gt %>% tab_header(title = "Time-to-event group sequential design") %>% cols_align("left")

# extend this to a conditional power design
xx <- ssrCP(x = x, z1 = z, overrun = overrun, beta = betastar, cpadj = cpadj, maxinc = maxinflation, z2 = z2)
plot(xx) # plot the stage 2 sample size
