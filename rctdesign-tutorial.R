library(gsDesign) # Keaven Anderson, AVP at Merck
library(RCTdesign) # Scott Emerson, UW

# how to design studies with unequally-spaced interim analyses?
# use sample.size argument

base <- seqDesign(prob.model = "hazard", arms = 2, ratio = c(1,1), null.hypothesis = 1, alt.hypothesis = 0.67, test.type = "less", alpha = 0.025)
update(base, nbr.analyses = 4)
update(base, nbr.analyses = 4, sample.size = c(0.25, 0.5, 0.75, 1))
update(base, nbr.analyses = 4, sample.size = c(0.10, 0.25, 0.75, 1))

# Definition of candidate designs for the CLL trial

base <- seqDesign(prob.model = "hazard", arms = 2, ratio = c(1,1), null.hypothesis = 1, alt.hypothesis = 0.67, test.type = "less", alpha = 0.025)
base <- update(base, sample.size = 263, power = "calculate")

Fixed.Sample <- update(base, nbr.analyses = 1)
SymmOBF.2 <- update(base, nbr.analyses = 2, P = c(1, 1))
SymmOBF.3 <- update(base, nbr.analyses = 3, P = c(1, 1))
SymmOBF.4 <- update(base, nbr.analyses = 4, P = c(1, 1))

SymmOBF.Power <- update(SymmOBF.4, power = 0.901 )

Futility.5 <- update(SymmOBF.4, P = c(1,.5))
Futility.8 <- update(SymmOBF.4, P = c(1,.8))
Futility.9 <- update(SymmOBF.4, P = c(1,.9))

Eff11.Fut8 <- update(SymmOBF.4, P = c(1.1,.8))
Eff11.Fut9 <- update(SymmOBF.4, P = c(1.1,.9))

Fixed.Power <- update(Fixed.Sample, nbr.analyses = 1, power = 0.8853)

seqBoundary(Eff11.Fut8, scale = "X")
seqBoundary(Eff11.Fut8, scale = "Z")
1 - seqBoundary(Eff11.Fut8, scale = "P")
 
# Figure 1: Comparison of stopping boundaries on crude estimate of treatment effect scale
seqPlotBoundary(SymmOBF.4, Eff11.Fut8, Eff11.Fut9, lty = c(1,3,4), col = 1, stagger = 0, fixed = FALSE )

# Figure 2 : Comparison of statistical power curves
seqPlotPower(SymmOBF.4,SymmOBF.3,SymmOBF.2, lty=1:4, col=1, lwd=2)
seqPlotPower(SymmOBF.4,SymmOBF.3,SymmOBF.2, reference=TRUE, lty=1:4, col=1, lwd=2 )

# Table 1 : Computation of power and alternative tables for the Eff11.Fut8 design
seqOC(Eff11.Fut8, power=c(.8,.9,.95,.975))
seqOC(Eff11.Fut8, theta=c(1,.75,.67,.60))

# Figure 3 : Comparison of sample size distributions
seqPlotASN(SymmOBF.4,Futility.9,Futility.8,Futility.5, fixed=FALSE, lty=c(2,1,3,4), col=1, lwd=2) # How the ASN were derived?

# Figure 4 : Depiction of stopping probabilities
seqPlotStopProb(Eff11.Fut8) # where are the numbers in this graph stored?

# Figure 5 : Statistical inference on the boundaries
plot(seqInference(Eff11.Fut8))

# Figure 6a : Patient accrual patterns (early accrual)
Eff11.Fut8Extd.early <- update(Eff11.Fut8, accrualSize = 400, accrualTime = 3, bShapeAccr = 10, eventQuantiles = 16/12, nPtsSim = 10000, seed = 0)
seqPlotPHNSubjects(Eff11.Fut8Extd.early)
Eff11.Fut8Extd.late <- update(Eff11.Fut8, accrualSize = 400, accrualTime = 3, aShapeAccr = 10, eventQuantiles = 16/12, nPtsSim = 10000, seed = 0)
seqPlotPHNSubjects(Eff11.Fut8Extd.late)

# Simulation of CLL data
n <- 200 # number of events
grp1 <- rexp( n, rate = .75 * log(2) ) # note in exponential distribution, median equals log(2) / lambda
grp2 <- rexp( n, rate=(.75*log(2))*.70 )
trueSurv <- c( grp1, grp2 )
entry <- runif( 2*n, 0, 3 ) # 400 accured in 3 years
grp <- rep( 0:1, each=n )

analysisTime <- 1.5 # First analysis at 1.5 years after study start
obsSurv <- ifelse( trueSurv + entry <= analysisTime, trueSurv, analysisTime-entry )
event <- ifelse( obsSurv == trueSurv, 1, 0 )
cllData <- as.data.frame( cbind( grp, entry, obsSurv, event ) )
cllData <- cllData[ cllData$obsSurv > 0, ]
resp <- Surv( cllData$obsSurv, cllData$event )
interim1 <- seqMonitor( Eff11.Fut8, response=resp, treatment=cllData$grp, future.analyses=c(132,198,263) )

analysisTime <- 2.75 # Second analysis at 2.75 years after study start
obsSurv <- ifelse( trueSurv + entry <= analysisTime, trueSurv, analysisTime-entry )
event <- ifelse( obsSurv == trueSurv, 1, 0 )
cllData <- as.data.frame( cbind( grp, entry, obsSurv, event ) )
cllData <- cllData[ cllData$obsSurv > 0, ]
resp <- Surv( cllData$obsSurv, cllData$event )
interim2 <- seqMonitor( interim1, response=resp, treatment=cllData$grp, future.analyses=c(198,263) )

analysisTime <- 3.5 # Third analysis at 3.5 years after study start
obsSurv <- ifelse( trueSurv + entry <= analysisTime, trueSurv, analysisTime-entry )
event <- ifelse( obsSurv == trueSurv, 1, 0 )
cllData <- as.data.frame( cbind( grp, entry, obsSurv, event ) )
cllData <- cllData[ cllData$obsSurv > 0, ]
resp <- Surv( cllData$obsSurv, cllData$event )
interim3 <- seqMonitor( interim2, response=resp, treatment=cllData$grp, future.analyses=c(263) )
interim3.fit = survfit(resp ~ cllData$grp)
plot(interim3.fit, xlim = c(0, 4))

# Figure 8 : Comparison of implemented and original design
plot( interim3, dsnLbls=c("Implemented Design", "Original Design") )

## ----------------------------------------------------------------------------------------------------------------------------------------------------

