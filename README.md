# GenHyp-cpp
Generalized hyperbolic distribution in C++. On Windows, compile with a command such as

```
g++ -IC:\cpp\boost_1_87_0 -LC:\cpp\boost_1_87_0\stage\lib xgenhyp.cpp
```
if Boost in installed in `C:\cpp\boost_1_87_0`. Output:

```
--------------------------------------------
Parameter Set 1:
  lambda = 1
  alpha  = 2
  beta   = 0.5
  delta  = 1
  mu     = 0

Theoretical Results (at x=0 for density):
  Density at x = 0: 0.430768
  Mean     = 0.475739
  Variance = 1.04554

Numerical Integration Results (over [-10, 10]):
  Integral of density (should be ~1): 0.999999
  Numerical Mean  = 0.475733
  Numerical Variance = 1.04548
--------------------------------------------

--------------------------------------------
Parameter Set 2:
  lambda = 0.5
  alpha  = 3
  beta   = 0.8
  delta  = 2
  mu     = 1

Theoretical Results (at x=0 for density):
  Density at x = 0: 0.0779438
  Mean     = 1.64907
  Variance = 0.882601

Numerical Integration Results (over [-9, 11]):
  Integral of density (should be ~1): 1
  Numerical Mean  = 1.64907
  Numerical Variance = 0.8826
--------------------------------------------

--------------------------------------------
Parameter Set 3:
  lambda = 1.5
  alpha  = 2.5
  beta   = 0.3
  delta  = 1.5
  mu     = -0.5

Theoretical Results (at x=0 for density):
  Density at x = 0: 0.412861
  Mean     = -0.210976
  Variance = 0.983133

Numerical Integration Results (over [-10.5, 9.5]):
  Integral of density (should be ~1): 1
  Numerical Mean  = -0.210976
  Numerical Variance = 0.983132
--------------------------------------------

--------------------------------------------
Parameter Set 4:
  lambda = 2
  alpha  = 4
  beta   = 1
  delta  = 2
  mu     = 0

Theoretical Results (at x=0 for density):
  Density at x = 0: 0.361845
  Mean     = 0.697283
  Variance = 0.756659

Numerical Integration Results (over [-10, 10]):
  Integral of density (should be ~1): 1
  Numerical Mean  = 0.697283
  Numerical Variance = 0.756659
--------------------------------------------
```
