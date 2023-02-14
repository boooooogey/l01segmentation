# l01segmentation
l01segmentation provides functions to segment an input signal (y) by solving a L0-fused approximation problem, given below: 

![dark](https://user-images.githubusercontent.com/15932827/158258314-21560f85-d37c-4477-a453-0dc5d989ac17.svg#gh-dark-mode-only)
![light](https://user-images.githubusercontent.com/15932827/158258342-d06edc11-11ba-4a73-aaba-9d9e279cee44.svg#gh-light-mode-only)

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

There are also two end-to-end pipeline for reading from a BigWig file, segmenting the signal, and storing the compressed signal in a BedGraph file.

compressBigwig compresses a signal by both binning and using the L0 solver.  
```
compressBigwig(BigWig, BedGraph, bin, lambda2)
```
- **BigWig:** path to the BigWig file.
- **BedGraph:** path to the BedGraph file in which the output will be stored.
- **bin:** bin size of the initial binning.
- **lambda2:** hyperparameter value for L0 solver.

compressBigwig compresses a signal by binning.
```
compressBigwigbyBinning(BigWig, BedGraph, bin)
```

- **BigWig:** path to the BigWig file.
- **BedGraph:** path to the BedGraph file in which the output will be stored.
- **bin:** bin size.

## Example
**Code:**
```
fusedsegmentation(data, lambda2 = 125, objective="poisson", format="full")
```
**Output (format = "compressed"):**
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

**Code:**
```
compressBigwig("pol2.bw", "pol2_100.bedGraph", lambda2 = 100, bin = 20)
compressBigwigbyBinning("pol2.bw", "pol2_bin.bedGraph", bin = 10000)
```
**Figure:**
![github_igv](https://user-images.githubusercontent.com/15932827/189454800-ea22ac5f-8e0c-442c-bed0-81349c8b496d.png)

