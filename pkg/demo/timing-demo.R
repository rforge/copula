require(nacopula)

## ==== Results of using timings() =============================================

set.seed(1) # set seed

n <- 10000
taus1 <- c(0.05,(1:3)/10) # for AMH: cannot attain more concordance
taus2 <- c(0.05,(1:7)/10) # for Frank: theta1 corresponding to tau > 0.8 is so large that 1-exp(-theta1) is 1 which leads to an error in rLog()
taus3 <- c(0.05,(1:9)/10,0.95) # "whole" range

## ==== AMH ====

timing(n,"AMH",taus1)

##       0.05   0.1   0.2   0.3
## 0.05 0.002 0.004 0.004 0.005
## 0.1  0.002    NA 0.004 0.004
## 0.2  0.002    NA    NA 0.006
## 0.3  0.002    NA    NA    NA

## conclusion:
## V0:  uniformly fast
## V01: uniformly fast

## ==== Clayton ====

timing(n,"Clayton",taus3)

##       0.05   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9  0.95
## 0.05 0.002 0.037 0.044 0.055 0.062 0.057 0.058 0.052 0.047 0.037 0.034
## 0.1  0.002    NA 0.053 0.060 0.057 0.058 0.062 0.062 0.058 0.054 0.051
## 0.2  0.002    NA    NA 0.041 0.042 0.038 0.043 0.045 0.043 0.042 0.042
## 0.3  0.002    NA    NA    NA 0.029 0.029 0.029 0.030 0.030 0.026 0.030
## 0.4  0.003    NA    NA    NA    NA 0.021 0.022 0.022 0.019 0.022 0.023
## 0.5  0.002    NA    NA    NA    NA    NA 0.017 0.018 0.018 0.018 0.018
## 0.6  0.003    NA    NA    NA    NA    NA    NA 0.014 0.015 0.016 0.015
## 0.7  0.002    NA    NA    NA    NA    NA    NA    NA 0.013 0.013 0.013
## 0.8  0.002    NA    NA    NA    NA    NA    NA    NA    NA 0.012 0.012
## 0.9  0.002    NA    NA    NA    NA    NA    NA    NA    NA    NA 0.013
## 0.95 0.002    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA

## conclusion:
## V0:  uniformly fast
## V01: faster for larger parameter (overall: fast enough)

## ==== Frank (todo: theta0 large) ====

timing(n,"Frank",taus2)

##       0.05   0.1   0.2   0.3   0.4   0.5   0.6    0.7
## 0.05 0.001 0.116 0.152 0.111 0.111 0.111 0.113  0.111
## 0.1  0.001    NA 0.117 0.119 0.120 0.118 0.119  0.119
## 0.2  0.002    NA    NA 0.148 0.149 0.149 0.147  0.147
## 0.3  0.002    NA    NA    NA 0.244 0.239 0.239  0.241
## 0.4  0.002    NA    NA    NA    NA 0.569 0.566  0.564
## 0.5  0.002    NA    NA    NA    NA    NA 2.523  2.425
## 0.6  0.002    NA    NA    NA    NA    NA    NA 20.337
## 0.7  0.002    NA    NA    NA    NA    NA    NA     NA

## conclusion:
## V0:  uniformly fast
## V01: for a given theta0: uniformly fast over theta1; but increasing in theta1

## ==== Gumbel ====

timing(n,"Gumbel",taus3)

##       0.05   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9  0.95
## 0.05 0.014 0.008 0.008 0.009 0.011 0.009 0.009 0.009 0.009 0.009 0.009
## 0.1  0.008    NA 0.009 0.008 0.009 0.009 0.010 0.009 0.009 0.010 0.009
## 0.2  0.008    NA    NA 0.009 0.009 0.009 0.005 0.009 0.010 0.009 0.009
## 0.3  0.008    NA    NA    NA 0.009 0.009 0.009 0.009 0.010 0.009 0.009
## 0.4  0.008    NA    NA    NA    NA 0.009 0.009 0.009 0.007 0.008 0.010
## 0.5  0.005    NA    NA    NA    NA    NA 0.009 0.009 0.010 0.009 0.009
## 0.6  0.007    NA    NA    NA    NA    NA    NA 0.009 0.009 0.009 0.009
## 0.7  0.007    NA    NA    NA    NA    NA    NA    NA 0.009 0.009 0.009
## 0.8  0.007    NA    NA    NA    NA    NA    NA    NA    NA 0.005 0.009
## 0.9  0.008    NA    NA    NA    NA    NA    NA    NA    NA    NA 0.009
## 0.95 0.009    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA

## conclusion:
## V0:  uniformly fast
## V01: uniformly fast

## ==== Joe ====

timing(n,"Joe",taus3)

##       0.05   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9  0.95
## 0.05 0.001 0.368 0.354 0.401 0.375 0.372 0.372 0.378 0.372 0.382 0.380
## 0.1  0.002    NA 0.365 0.369 0.374 0.384 0.386 0.386 0.385 0.381 0.384
## 0.2  0.003    NA    NA 0.380 0.396 0.407 0.423 0.476 0.433 0.444 0.443
## 0.3  0.004    NA    NA    NA 0.419 0.457 0.502 0.537 0.551 0.575 0.570
## 0.4  0.004    NA    NA    NA    NA 0.550 0.677 0.724 0.746 0.792 0.834
## 0.5  0.005    NA    NA    NA    NA    NA 0.905 1.026 1.228 1.372 1.446
## 0.6  0.006    NA    NA    NA    NA    NA    NA 1.237 1.539 1.928 1.997
## 0.7  0.006    NA    NA    NA    NA    NA    NA    NA 1.776 2.443 2.694
## 0.8  0.006    NA    NA    NA    NA    NA    NA    NA    NA 2.123 2.591
## 0.9  0.006    NA    NA    NA    NA    NA    NA    NA    NA    NA 1.541
## 0.95 0.006    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA

## conclusion:
## V0:  (almost) uniformly fast (at least fast enough)
## V01: increasing in theta0 and theta1
## V01 depends on theta0 and theta1 [but is bounded due to approx-parameter]
