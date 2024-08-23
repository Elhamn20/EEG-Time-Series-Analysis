# An analysis of the electroencefalogram (eeg) data from Prado, Ferreira and West (2021)
# using AR(p) models.
# By Marco A. R. Ferreira (2020).
# Elham Nasarian-Fall 2022

library(astsa)
library(polynom)

#Task1
eeg <- scan("eeg.dat")

# Plot the data
tsplot(eeg)

#Task2
# Plot the sample ACF
acf(eeg,lag.max=50)

# Plot the sample PACF
pacf(eeg,lag.max=50)

#Task3
# Compute AIC
ord.max = 8
eeg.model <- ar(eeg,order.max=ord.max,aic=TRUE,method="mle")
eeg.model$aic

# Compute BIC
n = length(eeg)
eeg.BIC <- eeg.model$aic - 2*((0:ord.max)+1) + log(n) * ((0:ord.max)+1)
eeg.BIC       # Choose model with smallest BIC
                     # BIC = Bayesian Information Criterion

# Posterior probabilities of the competing AR models
bic.vec <- eeg.BIC
post.prob <- exp(-0.5*bic.vec) / sum(exp(-0.5*bic.vec))
post.prob

#Task4
# Fit an AR(8) model with maximum likelihood
eeg.model8 <- ar(eeg,order.max=8 ,aic=FALSE,method="mle")
eeg.model8

charac.polyn.coeff <- c(1, -eeg.model8$ar)
charac.polyn.coeff

order = 8
MLE.fit = ar.mle(eeg, order=order, aic=FALSE) 
MLE.fit$x.mean
MLE.fit$ar
sqrt(diag(MLE.fit$asy.var.coef))

#Task5
z <- c(1, -MLE.fit$ar) 
a = polyroot(z)
a
# The periods
paste("The periods are as follows:")
(2*pi) / Arg(a)
paste("The reciprocal of the moduli are as follows:") 
1 / Mod(a)

