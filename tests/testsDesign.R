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
deswgts3 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, scal=200, off=1, method = "mult")

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
