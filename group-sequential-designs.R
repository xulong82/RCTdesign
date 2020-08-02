library(RCTdesign)

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

