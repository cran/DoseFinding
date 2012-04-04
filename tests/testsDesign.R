library(DoseFinding)

# Some examples from the JASA paper (for validation)
########################################################################
# Emax model p.1228 l. 5
fmods <- list(emax = 25)
doses <- c(0,150)
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4)
fMod$emax[2] <- 0.6666667
doses <- c(0, 18.75, 150)
weights <- c(1)
deswgts1 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="Nelder-Mead")
deswgts2 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="nlminb")

# Paper p. 1228 l. 2
fmods <- list(emax = 25)
doses <- c(0,150)
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4)
doses <- c(0, 18.75, 150)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.2)
deswgts

########################################################################
#### exponential
# Paper p.1229 2nd line
doses <- c(0, 150)
fmods <- list(exponential=c(85))
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4)
doses <- c(0, 50, 104.52, 150)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="Nelder-Mead")
deswgts

# Paper p.1229 1st line
doses <- c(0, 150)
fmods <- list(exponential=c(65))
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4)
fMod$exponential[2] <- 0.08264711
doses <- c(0, 101.57, 150)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.2)
deswgts

# Paper p.1229 2nd line
doses <- c(0, 150)
fmods <- list(exponential=c(85))
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4)
fMod$exponential[2] <- 0.08264711
doses <- c(0, 104.52, 150)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.2)
deswgts

########################################################################
#### Logistic
#### Paper: p.1230 7th line
doses <- c(0, 150)
fmods <- list(logistic=c(50, 10.881))
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4)
doses <- c(0, 37.29, 64.44, 150)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.05)
deswgts

#### Paper p.1230 line 1
doses <- c(0, 150)
fmods <- list(logistic=c(50, 10.881))
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4)
doses <- c(0, 50.22)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.2)
deswgts

########################################################################
#### beta
# Paper p.1230 line 5
fmods <- list(betaMod = c(0.33, 2.31))
doses <- c(0,150)
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4, scal=200)
doses <- c(0, 0.49, 25.2, 108.07, 150)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.1,
         scal=200, control=list(maxit=1000))
deswgts

# Paper p. 1230 line 10
fmods <- list(betaMod = c(1.39, 1.39))
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4, scal=200)
#doses <- c(0, 10, 25, 50, 100, 150)
doses <- c(0, 27, 94.89, 150)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.1,
         scal=200)
deswgts

# Paper p. 1230 line 1
fmods <- list(betaMod = c(0.23, 2.31))
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4, scal=200)
doses <- c(0, 10, 25, 50, 100, 150)
doses <- c(0, 0.35, 150)
weights <- c(1)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.2,
        scal=200)
deswgts

########################################################################
#### mixed Paper p. 1233, l. 2 (note the off and probably also the
#### scal parameter were treated as unknown in this example in the paper, 
#### hence the results need not be consistent with paper)
doses <- c(0, 9.9, 49.5, 115.4, 150)
fmods <- list(linear = NULL, emax = 25, exponential = 85,
               linlog = NULL, logistic = c(50, 10.8811))
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4, off=1)
weights <- rep(1/5, 5)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.2, scal=200, off=1)
deswgts2 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, scal=200, off=1, method = "nlminb")

# Some other examples 
########################################################################
fmods <- list(emax = c(25, 107.14), linear = NULL,
              logistic = c(150, 45.51), betaMod = c(1,1))
doses <- c(0, 62.5, 125, 250, 500)
fMod <- fullMod(fmods, doses, base=60, maxEff=280, scal=1.2*500)
weights <- rep(0.2, length=5)
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=200, scal=1.2*500)

########################################################################
#### using already allocated patients
fmods <- list(betaMod = c(0.33, 2.31))
doses <- c(0,150)
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4, scal=200)
doses <- c(0, 0.49, 25.2, 108.07, 150)
weights <- c(1)
# no previously allocated patients
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.1, scal=200, control=list(maxit=1000))

# now use previously allocated patients
nold <- c(45, 50, 0, 0, 0)
deswgts2 <- calcOptDesign(fMod, weights, doses, clinRel=0.1, n2=30, scal=200,
         off=1, control=list(maxit=1000), nold=nold)
# the overall design
(30*deswgts2$design+nold)/(30+sum(nold))
deswgts$design

########################################################################
#### Dopt Examples
fmods <- list(emax = c(25, 107.14), logistic = c(150, 45.51),
              linear = NULL, betaMod = c(1,1))
doses <- c(0, 62.5, 125, 250, 500)
fMod <- fullMod(fmods, doses, base=60, maxEff=280, scal=500*1.2)
weights <- rep(0.2, 5)
des1 <- calcOptDesign(fMod, weights, doses, clinRel = 200, scal = 500*1.2)
des2 <- calcOptDesign(fMod, weights, doses, clinRel = 200, scal = 500*1.2, type = "Dopt")
des3 <- calcOptDesign(fMod, weights, doses, clinRel = 200, scal = 500*1.2, type = "MED&Dopt")

########################################################################
#### method = "exact" and "solnp"
fmods <- list(emax = c(25, 107.14), logistic = c(150, 45.51),
              linear = NULL, betaMod = c(1,1))
doses <- c(0, 62.5, 125, 250, 500)
fMod <- fullMod(fmods, doses, base=60, maxEff=280, scal=500*1.2)
weights <- rep(0.2, 5)
des41 <- calcOptDesign(fMod, weights, doses, clinRel = 200, n2 = 10,
                      scal = 500*1.2, method = "exact",
                      lowbnd = c(0.3,0,0,0,0))
des42 <- calcOptDesign(fMod, weights, doses, clinRel = 200, 
                      scal = 500*1.2, method = "solnp",
                      lowbnd = c(0.1,0,0,0,0))
des51 <- calcOptDesign(fMod, weights, doses, clinRel = 200, n2 = 10,
                       scal = 500*1.2, type = "Dopt", method = "exact",
                       uppbnd = rep(0.5,5))
des52 <- calcOptDesign(fMod, weights, doses, clinRel = 200, 
                       scal = 500*1.2, type = "Dopt", method = "solnp",
                       uppbnd = rep(0.5,5))
des61 <- calcOptDesign(fMod, weights, doses, clinRel = 200, n2 = 10,
                       method = "exact", scal = 500*1.2,
                       type = "MED&Dopt")
des62 <- calcOptDesign(fMod, weights, doses, clinRel = 200, 
                       method = "solnp", scal = 500*1.2,
                       type = "MED&Dopt")

########################################################################
#### Example from Padmanabhan and Dragalin, Biometrical Journal 52 (2010)
#### p. 836-852
models <- list(sigEmax = c(4, 5))
doses <- 0:8
fm <- fullMod(models, doses, base=0, maxEff=-1.65)
fm$sigEmax <- c(0, -1.70, 4, 5)
## compare to Figure 1, p. 841
desSED <- calcOptDesign(fm, 1, doses, type="Dopt", method = "solnp")
desSEM <- calcOptDesign(fm, 1, doses, clinRel = -1.3, type="MED", method = "solnp")

## designs underlying Table 2, p. 843 (from an e-mail of Vlad)
## I cannot reproduce the displayed efficiencies exactly
## (most probably due to numerical round-off)
##LDoD
## [1,] 0.246 0.141 0.123 0.000 0.000 0.240    0    0 0.250
## [2,] 0.248 0.233 0.061 0.210 0.000 0.000    0    0 0.248
## [3,] 0.246 0.000 0.000 0.223 0.081 0.204    0    0 0.246
## [4,] 0.250 0.247 0.045 0.210 0.000 0.000    0    0 0.248
## [6,] 0.250 0.249 0.192 0.062 0.000 0.000    0    0 0.246
## MEDoD
## [1,] 0.49 0.01 0.00 0.00 0.00 0.00 0.36 0.14    0
## [2,] 0.49 0.02 0.00 0.15 0.35 0.00 0.00 0.00    0
## [3,] 0.23 0.26 0.01 0.00 0.00 0.46 0.04 0.00    0
## [4,] 0.50 0.00 0.49 0.01 0.00 0.00 0.00 0.00    0
## [6,] 0.49 0.01 0.47 0.02 0.00 0.00 0.00 0.00    0
models1 <- list(sigEmax = c(23.07, 1.18))
models2 <- list(sigEmax = c(2, 2.22))
models3 <- list(sigEmax = c(4, 5))
models4 <- list(sigEmax = c(0.79, 1))
models5 <- list(sigEmax = c(0.74, 1.18))
doses <- 0:8
fm <- list()
fm[[1]] <- fullMod(models1, doses, base=0, maxEff=-1.65);fm[[1]]$sigEmax <- c(0, -7.29, 23.07, 1.18)
fm[[2]] <- fullMod(models2, doses, base=0, maxEff=-1.65);fm[[2]]$sigEmax <- c(-0.08, -1.71, 2, 2.22)
fm[[3]] <- fullMod(models3, doses, base=0, maxEff=-1.65);fm[[3]]$sigEmax <- c(0, -1.70, 4, 5)
fm[[4]] <- fullMod(models4, doses, base=0, maxEff=-1.65);fm[[4]]$sigEmax <- c(0, -1.81, 0.79, 1.00)
fm[[5]] <- fullMod(models5, doses, base=0, maxEff=-1.65);fm[[5]]$sigEmax <- c(-0.03, -1.72, 0.74, 1.18)
desD <- desM <- matrix(ncol = 9, nrow = 5)
for(i in 1:5){
  cc1 <- calcOptDesign(fm[[i]], 1, doses, type="MED", method = "solnp",
                      clinRel = -1.3)
  cc2 <- calcOptDesign(fm[[i]], 1, doses, type="Dopt", method = "solnp")
  desM[i,] <- cc1$design
  desD[i,] <- cc2$design
}
round(desD, 3)
round(desM, 2)

################################################################################
#### look at standardized Dopt and MED&Dopt criteria
mods1 <- list(sigEmax = rbind(c(25, 5), c(107.14, 2)))
mods2 <- list(sigEmax = rbind(c(25, 5), c(107.14, 2)), linear = NULL)
doses <- c(0, 62.5, 125, 250, 500)
fMod1 <- fullMod(mods1, doses, base=60, maxEff=280)
fMod2 <- fullMod(mods2, doses, base=60, maxEff=280)
w1 <- rep(0.5, 2)
w2 <- rep(1/3, 3)
## des1 and des2 should be exactly the same
des1 <- calcOptDesign(fMod1, w1, doses, type = "Dopt", standDopt = FALSE)
des2 <- calcOptDesign(fMod1, w1, doses, type = "Dopt", standDopt = TRUE)

## des1 and des2 should be different (as linear and emax have
## different number of parameters)
## des1 <- calcOptDesign(fMod2, w2, doses, type = "Dopt", standDopt = FALSE,
##                       method = "solnp")
## des2 <- calcOptDesign(fMod2, w2, doses, type = "Dopt", standDopt = TRUE,
##                       method = "solnp")

## same with MED&Dopt criterion
## des1 and des2 will differ (due to different scaling
## of Dopt and MED criteria)
## des1 <- calcOptDesign(fMod1, w1, doses, type = "MED&Dopt",
##                       clinRel = 100, standDopt = FALSE,
##                       method = "solnp")
## des2 <- calcOptDesign(fMod1, w1, doses, type = "MED&Dopt",
##                       clinRel = 100, standDopt = TRUE,
##                       method = "solnp")

## des1 and des2 should be different (different no of parameters
## and different scaling of Dopt and MED criteria)
## des1 <- calcOptDesign(fMod2, w2, doses, type = "MED&Dopt",
##                       clinRel = 100, standDopt = FALSE,
##                       method = "solnp")
## des2 <- calcOptDesign(fMod2, w2, doses, type = "MED&Dopt",
##                       clinRel = 100, standDopt = TRUE,
##                       method = "solnp")
