# Barale - Shirke Test
**R** implementation of two samples location/scale multidimensional test proposed by Barale and Shirke.

## About the code

The test is performed by the function called ```baraleShirkeTest```. This function is defined in ```baraleShirke.cpp```. That files contains several functions that are used in the ```baraleShirkeTest``` function. One of them is ```mahDepth``` which computes the mahalanobis depth of a vector according to a mean vector and a covariances matrix. You can replace this function if you want to use a different depth meassure.

The test is implemented in **C++** in order to aproximate the p-value with a higer precision at a lower cost. Therefore, it is necesary to install a **C++** compiler, ```Rcpp``` and ```RcppArmadillo``` **R** packages, and ```Armadillo``` **C++** library.

## How to use it

Assuming that you have installed all the dependencies previously mentioned, you should load the following libraries in **R**:

```R
library(Rcpp)
library(RcppArmadillo)
```
Now, assuming that the ```baraleShirke.cpp``` file is in your working directory, you must compile the file as follows:

```R
sourceCpp("baraleShirke.cpp")
```

Finally you're ready to start performing the test. You must know that the ```baraleShirkeTest``` function has three arguments that must be passed in the following order:

* ```X1``` which is the **matrix** that contains the observations of the first sample on its rows.

* ```X2``` which is the **matrix** that contains the observations of the sample sample on its rows.

* ```B```, the number of iterations that will be used to approximate the p-value.

## Files

Here you can fins a list of the files that are contained on this repository:

* ```baraleShirke.cpp```: Source code un **C++**.
* ```baraleShirke.r```: Example code in **R**.
* ```turtles.csv```: Example data (columns are separated by spaces)-
* ```README.md```: This file.


## References
* M. S. Barale & D. T. Shirke (2020): *A test based on data depth for testinglocation-scale of the two multivariate populations*, Journal of Statistical Computation andSimulation, DOI: 10.1080/00949655.2020.1830285
