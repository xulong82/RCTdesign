library(RCTdesign)

# Binary endpoints

xxx = seqDesign(prob.model = "proportion", arms = 2, null.hypothesis = .3, alt.hypothesis = .23, 
                nbr.analyses = 1, test.type = "less", power = .8, alpha = .025)

update(xxx, display.scale = "P")
xxx = update(xxx, display.scale = "Z")

kkk = update(xxx, test.type = "two.sided")
update(kkk, test.type = "two.sided", alpha = 0.025)
update(kkk, test.type = "two.sided", alpha = 0.05)

xx2 = update(xxx, sample.size = 1700, power = "calculate")

update(xx2, display.scale = "Z")
update(xx2, display.scale = "P")

obf2 = update(xx2, nbr.analyses = 2)
obf3 = update(xx2, nbr.analyses = 3)
obf4 = update(xx2, nbr.analyses = 4)
obf4x = update(obf4, early.stopping = "alternative")
obf4y = update(obf4, early.stopping = "null") # stop only for futility
obf4z = update(obf4, early.stopping = "null", nonbinding.futility = T) # 

seqPlotBoundary(obf2, obf3, obf4)
seqPlotASN(obf2, obf3, obf4)
seqPlotPower(obf2, obf3, obf4)
seqPlotBoundary(obf4, obf4x)

obf2 = update(xxx, nbr.analyses = 2)
obf3 = update(xxx, nbr.analyses = 3)
obf4 = update(xxx, nbr.analyses = 4)

seqPlotASN(obf2, obf3, obf4)

obf4x = update(obf4, early.stopping = "alternative")
obf4y = update(obf4, early.stopping = "null") # stop only for futility

obf4z = update(obf4, P = c(1, 0.8)) # stop only for futility
seqPlotBoundary(obf4, obf4z)

obf4k = update(obf4, P = c(1, Inf)) # stop only for futility
obf4f = update(obf4, early.stopping = "alternative") # stop only for futility

seqPlotPower(obf4, obf4x, reference = T)

update(xxx, sample.size = 1700, power = "calculate")
update(xxx, sample.size = 1700, power = "calculate", display.scale = "Z")

update(xxx, sample.size = 700, power = "calculate")

update(xxx, null.hypothesis = .4, alt.hypothesis = 0.35)

seqInference(xxx)

# time to event

xxx = seqDesign(prob.model = "hazard", arms = 2, null.hypothesis = 1, alt.hypothesis = 0.67, test.type = "less", power = 0.8)
seqPHSubjects(xxx, controlMedian = 9, accrualTime = 36, followupTime = 12)

# Early stopping only for efficacy, for efficacy and futility, and how would it affect the efficacy boundaries?

base <- seqDesign(prob.model = "hazard", arms = 2, ratio = c(1,1), null.hypothesis = 1, alt.hypothesis = 0.5, test.type = "less", alpha = 0.025, power = 0.9)
(d1J4 = update(base, nbr.analyses = 4, display.scale = "P")) # group sequential design w early stopping for both efficacy and futility
(d1J4x = update(d1J4, early.stopping = "alternative"))

plot(d1J4, d1J4x) 
(d1J4.eval = seqEvaluate(d1J4))

# Operation Characteristics w/o early stopping for futility as oppose to early stopping for both efficacy and futility
# 1. power becomes higher with the same sample size, so maximal sample size decreased w the same target power
# 2. false positive rates becomes higher w the same efficacy boundary, so the efficacy boundary increased w the same target alpha

# Early stopping for futility: binding vs non-binding

base <- seqDesign(prob.model = "hazard", arms = 2, ratio = c(1,1), null.hypothesis = 1, alt.hypothesis = 0.5, test.type = "less", alpha = 0.025, power = 0.9)
(d1J4 = update(base, nbr.analyses = 2, display.scale = "P", nonbinding.futility = F)) 
(d1J4x = update(d1J4, early.stopping = "alternative"))

(d1J4y = update(d1J4, nonbinding.futility = T)) 
(d1J4z = update(d1J4, alpha = 0.02471))

# replicate w gsdesign
gsSurv(ratio = 1, hr0 = 1, hr = 0.5, sided = 1, alpha = 0.025, beta = 0.1, k = 4, T = 48, minfup = 0, gamma = 1, test.type = 3)

#######################################################################################

# Different design family

(d1J4y = update(d1J4, design.family = "OBF"))
(d1J4z = update(d1J4, design.family = "Pocock"))

plot(d1J4y, d1J4z)

# What exact are the different boundary shape parameters (the unified family of group sequential stopping rules)

y1 = gsDesign::nSurv(sided = 1, ratio = 1, alpha = 0.025, beta = 0.1, hr0 = 1, hr = 0.5, lambdaC = 1/5, eta = 1/10, gamma = 20, R = 8, T = 20)
y1$d; y1$n; y1$T

y2 = seqDesign(prob.model= "hazard", arms= 2, ratio= c(1,1), null.hypothesis= 1, alt.hypothesis= 0.5, test.type= "less", size= 0.025, power= 0.90)
y2$parameters

y3 = seqDesign(prob.model= "hazard", arms= 2, ratio= c(1,1), null.hypothesis= 1, alt.hypothesis= 0.5, test.type= "less", size= 0.025, power= 0.90,
                accrualTime = 8, studyTime = 20, eventHazard = 1/5, dropoutHazard = 1/10)
y3$parameters

(d1 = seqDesign(prob.model= "mean", arms= 2, ratio= c(1,1), null.hypothesis= 0, alt.hypothesis= -10, sd= c(30,30), test.type= "less", size= 0.025, power= 0.90))
(d1J2 = update(d1, nbr.analyses = 2))
(d1J4 = update(d1, nbr.analyses = 4))

