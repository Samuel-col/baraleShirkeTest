# Barale - Shirke Test
**R** implementation of two samples location/scale multidimensional test proposed by (Barale, 2020).

## About the code

The test is performed by the function called ```baraleShirkeTest```. This function is defined in ```baraleShirke.cpp```. That files contains several functions that are used in the ```baraleShirkeTest``` function. ```baraleShirke.R``` contains the function ```baraleShirke.test``` which is the function that the user shoud call.

## Dependencies

The test is implemented in **C++** in order to aproximate the p-value with a higher precision at a lower cost. Therefore, it is necesary to install a **C++** compiler, ```ddalpha```, ```Rcpp``` and ```RcppArmadillo``` **R** packages, and ```Armadillo``` **C++** library.

## How to use it

Assuming that you have installed all the dependencies previously mentioned, you should load the function:

```R
source("path/to/baraleShirke.R")
```

Now you're ready to start performing the test. You must know that the ```baraleshirke.test``` function has seven arguments that must be passed in the following order:

* ```X1``` which is the **matrix**, **data.frame** or **array** that contains the observations of the first sample on its rows.

* ```X2``` which is the **matrix**, **data.frame** or **array** that contains the observations of the second sample on its rows.

* ```depth``` which is a string that specifies the depth measure that is going to be used. It should be one of the following:

  - halfspace
  - L2
  - projection
  - potential
  - qhpeeling
  - simplicial
  - spatial
  - zonoid

* ```NIter``` which is the number of iterations that will be used to approximate the p-value.

* ```alpha``` which is significance level of the test.

* ```returnDepths``` which is a boolean that indicates if the depths should be returned.

* ```returnSamples``` which is a boolean that indicates if the resampling samples should be returned.

Arguments 3 to 7 are optional.

The function returns a list with the following entries:

* ```Statistic```: The observed Bmax statistic.
* ```NIter```: Number of iterations used to approximate the p-value.
* ```PValue```: Approximated p-value.
* ```n1```: Number of observations in sample 1.
* ```n2```: Number of observations in sample 2.
* ```alpha```: Significance level of the test.
* ```Depth```: Depth measure used in the test.
* ```DepthVals```: Computed depths (if ```returnDepths == TRUE```). Here the first column contains the depths with respect to sample 1 and the second column, with respect to sample 2.
* ```Samples```: Simulated statistics (if ```returnSamples == TRUE```).
* ```Message```: Nice message showing test output.


## Files

Here you can find a list of the files that are contained on this repository:

* ```baraleShirke.cpp```: Source code in **C++**.
* ```baraleShirke.R```: ```baraleshirke.test``` function definition.
* ```test.r```: Example code in **R**.
* ```turtles.csv```: Example data (columns are separated by spaces).
* ```README.md```: This file.


## References
* M. S. Barale & D. T. Shirke (2020): *A test based on data depth for testinglocation-scale of the two multivariate populations*, Journal of Statistical Computation andSimulation, DOI: 10.1080/00949655.2020.1830285
