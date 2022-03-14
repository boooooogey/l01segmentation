# l01segmentation
l01segmentation provides functions to segment an input signal (y) by solving a L0-fused approximation problem, given below: 

![](https://latex.codecogs.com/svg.image?{\color{white}\beta=\underset{\beta&space;\in&space;\mathbb{R}^{n}}{\operatorname{argmax}}&space;\sum_{i=1}^{N}&space;e(\beta_i,&space;y_i\textbf{})-\lambda_2&space;\sum_{i=1}^{N-1}\left|\beta_{i}&space;-&space;\beta_{i&plus;1}\right|_0}#gh-dark-mode-only)

Possible distributions that e can be derived from are:
- Gaussian
- Poison
- Binomial
- Exponential

## Required
- Rcpp
- BH

## Usage
There is one function *fusedsegmentation* that calls different solvers written in C++ depending on the probability distribution used. 
```
fusedsegmentation(y, lambda2, C, N, weight, l, objective, format)
```
- **y:** input vector.
- **lambda2:** hyper parameter that controls the strength of the L0 fused penalty.
- **C:** for now, only used for binomial case, the user needs to pass total number of trials with this argument.
- **N:** instead of using lambda2, the user may give the number of segments (N) they want to have. This is solved much slower than the problem with lambda2. 
- **weight:** user may choose to assign different weights to indices.
- **l:** two options, 0 or 1. 0 points to L0-fused, 1 points to L1-fused. L1-fused (fused lasso) solution is not available for all distributions.
- **objective:** options are "gauss", "poisson", "binomial", "exponential".
- **format:** options are "full", "compressed". For "full", the function returns a vector of the same length as the input. For "compressed", the function returns a list of segment start-ends and the values.

## Example
**Code:**
```
fusedsegmentation(data, lambda2 = 125, objective="poisson", format="full")
```
**Output:**
```
   start   end     value
1      1  1174 1.0732538
2   1175  2744 0.2522293
3   2745  4451 1.0661980
4   4452  5241 2.7037975
5   5242  5645 0.5123762
6   5646  6529 1.9355204
7   6530  7433 4.8772124
8   7434  9027 3.2013802
9   9028 12525 5.1818182
10 12526 12764 0.2426778
11 12765 14081 3.3697798
12 14082 15162 1.8057354
13 15163 15323 0.0000000
14 15324 16908 1.1356467
```
**Figure:**
![github1](https://user-images.githubusercontent.com/15932827/158254225-143b22b7-c427-4808-bfca-d9fefd545d6e.png)

