library(gsDesign)
library(tibble)
library(ggplot2)
library(scales)
library(gt)

###### events and dropout distributions

alpha <- .025; beta <- .1 

hr0 <- 1 ; hr1 <- .75 

median <- 12 # median control time-to-event
lambda <- log2(12) # exponential events rate per unit time: median events time 12 
eta <- .001 # exponential dropout rate per unit time

###### enrollment and trial duration

T <- 36 # study duration
minfup <- 12 # follow-up duration of last patient enrolled
R <- c(1, 2, 3, 4) # enrollment period duration
gamma <- c(1,1.5,2.5,4) # relative enrollment rates during above periods
ratio <- 1 # randomization ratio

###### deriving design with no interim analyses

x <- nSurv(R = R,
           minfup = minfup,
           gamma = gamma, # this becomes the relative enrollment rate
           T = T, # this will fix the accrual duration and solve for the accrual rate
           eta = eta, lambdaC = lambda, hr = hr1, hr0 = hr0,
           beta = beta, alpha = alpha)

y <- nSurv(R = R,
           minfup = minfup,
           gamma = gamma, # this becomes the absolute enrollment rates
           T = NULL, # this will fix the accrual rate and solve for the accrual duration 
           eta = eta, lambdaC = lambda, hr = hr1, hr0 = hr0,
           beta = beta, alpha = alpha)

###### group sequential design

k <- 3 # number of analyses (interim + final)
timing <- c(.25, .75) # proportion of final events at each interim analysis: k-1
sfu <- sfLDOF # efficacy bound spending function, Lan-DeMets spending function approximating O'Brien-Fleming bound
sfupar <- NULL # no parameter required for this spending function
sfl <- sfHSD # futility bound spending function
sflpar <- -7 # futility bound spending parameter specification

x1 <- gsSurv(R = R, gamma = gamma, minfup = minfup, T = T, 
             eta = eta, lambdaC = lambda, hr = hr1, hr0 = hr0 , beta = beta, alpha = alpha,
             k = k, timing = timing, test.type = 1)

x2 <- gsSurv(R = R, gamma = gamma, minfup = minfup, T = T, 
             eta = eta, lambdaC = lambda, hr = hr1, hr0 = hr0 , beta = beta, alpha = alpha,
             k = k, timing = timing, test.type = 2)

x = x1
tibble(Period =rownames(x$gamma), Rate = as.numeric(x$gamma)) %>% gt %>% tab_header(title = "Enrollment Rate")
x = x2
tibble(Period =rownames(x$gamma), Rate = as.numeric(x$gamma)) %>% gt %>% tab_header(title = "Enrollment Rate")

gsBoundSummary(x1) %>% gt %>% tab_header(title = "Time-to-event group sequential design") %>% cols_align("left")
gsBoundSummary(x2) %>% gt %>% tab_header(title = "Time-to-event group sequential design") %>% cols_align("left")
summary(x)

tibble(Period =rownames(x$gamma), Rate = as.numeric(x$gamma)) %>% gt %>% tab_header(title = "Enrollment Rate")

plot(x, plottype = "power", xlab = "HR") + scale_y_continuous(labels = scales::percent)

Z = c(0.25, 2)

gsCP(x, i = 2, zi = Z[2])

plot(x, plottype = "B", base = TRUE, main = "B-value projection", lty = 1, col = 1:2, xlab = "Events")

analysis = 2
N <- c(0, x$n.I[1:analysis]); B <- c(0, Z * sqrt(x$timing[1:analysis]))
points(x = N, y = B, col = 4); lines(x = N, y = B, col = 4)
slope <- B[analysis+1] / N[analysis + 1]
Nvals <- c(N[analysis+1], max(x$n.I))
lines(x = Nvals, y = B[analysis + 1] + c(0, slope * (Nvals[2] - Nvals[1])), col = 4, lty = 2)

################# tests

### parameter "sided" is used to interpret the alpha
x <- gsSurv(R = R, gamma = gamma, minfup = minfup, T = T, 
            eta = eta, lambdaC = log(2) / median, hr = hr, hr0 = hr0 , beta = beta, alpha = alpha, sided = 2,
            k = k, timing = timing, sfu = sfu, sfupar = sfupar, sfl = sfl, sflpar = sflpar)

### binding vs non-binding
z <- gsSurv(R = R, gamma = gamma, minfup = minfup, T = T, 
            eta = eta, lambdaC = log(2) / median, hr = hr, hr0 = hr0 , beta = beta, alpha = alpha, sided = 1,
            k = 2, timing = c(.75), sfu = sfLDOF, sfl = sfHSD, test.type = 3)

z <- gsSurv(R = R, gamma = gamma, minfup = minfup, T = T, 
            eta = eta, lambdaC = log(2) / median, hr = hr, hr0 = hr0 , beta = beta, alpha = alpha, sided = 1,
            k = 2, timing = c(.75), sfu = sfLDOF, sfl = sfHSD, test.type = 4)

z <- gsSurv(R = R, gamma = gamma, minfup = minfup, T = T, 
            eta = eta, lambdaC = log(2) / median, hr = hr, hr0 = hr0 , beta = beta, alpha = alpha, sided = 1,
            k = 2, timing = c(.75), sfu = sfLDOF, sfl = sfHSD, test.type = 5)

z <- gsSurv(R = R, gamma = gamma, minfup = minfup, T = T, 
            eta = eta, lambdaC = log(2) / median, hr = hr, hr0 = hr0 , beta = beta, alpha = alpha, sided = 1,
            k = 2, timing = c(.75), sfu = sfLDOF, sfl = sfHSD, test.type = 6)

summary(z)
gsBoundSummary(z) %>% gt %>% tab_header(title = "Time-to-event group sequential design")

### BILL script on 024 study 

eEvents(lambda=c(.05,.02,.01), eta=.01, gamma=c(5,10,20),
        R=c(2,1,2), S=c(1,1), T=20)

nSurv(lambdaC=log(2)/6, hr=.6, hr0=1, eta = 0, etaE=NULL,
      gamma=1, R=12, S=NULL, T=NULL, minfup = NULL, ratio = 1,
      alpha = 0.025, beta = 0.10, sided = 1, tol = .Machine$double.eps^0.25)
      
nSurv(alpha=0.025, sided=1, beta=0.1, lambdaC=log(2)/4.5, 
       hr=.6, hr0=1, eta=0, etaE=NULL,
       gamma=1, T=24, minfup=12, ratio=2)

nSurv(alpha=0.025, sided=1, beta=0.2, lambdaC=log(2)/4.5, 
            hr=.6, hr0=1, eta=0, etaE=NULL,
            gamma=1, T=24, minfup=12, ratio=2)

gsSurv(k=2, test.type=1, alpha=0.025, sided=1,
             beta=0.1, astar=0, timing=1, sfu=sfLDOF, r=18, 
             lambdaC=log(2)/4.5, 
             hr=.6, hr0=1, eta=0, etaE=NULL,
             gamma=1, T=24, minfup=12, ratio=2)

x<-gsDesign(k=2)
gsBoundCP(x)

# XW script on 024 study 

x1 = - log(0.8) / 6
x2 = - log(0.6) / 6

y1 = gsSurv(hr0 = 1, hr = 0.6, alpha = 0.025, sided = 1, beta = 0.2, ratio = 2, lambdaC = log(2)/4.5,  # test parameter 
            k = 2, timing = 0.4, test.type = 1, sfl = sfLDOF, # boundary
            R = 22, T = 28, minfup = 6, eta = x1, etaE = x2) # accrual and dropouts

y2 = nSurv(hr0 = 1, hr = 0.6, alpha = 0.025, sided = 1, beta = 0.2, ratio = 2, lambdaC = log(2)/4.5,  # test parameter 
           R = 22, T = 28, minfup = 6, eta = x1, etaE = x2) # accrual and dropouts

y2
# what is the median time to drop-out? eta parameters?

