# L0_Branch-and-Bound

This repository contains all the code necessary to replicate the findings described in the paper: [Global optimization for sparse solution of least squares problems](https://hal.archives-ouvertes.fr/hal-02066368/document). If you use it in your research, please cite:

> Ramzi Ben Mhenni, Sébastien Bourguignon, Jordan Ninin. Global Optimization for Sparse Solution of Least Squares Problems. 2019. ⟨hal-02066368⟩}

The dedicated branch-and-bound algorithms contained in this repository are for three possible formulations:
cardinality-constrained and cardinality-penalized least-squares, and cardinality minimization under quadratic constraints.

# Structure of the repository
* _./L0_optimization/_ contains the code for the MIP solver and the BaB solver.
* _./armadillo-8.300.2/_
* _./data/_  


# Running the code

## Dependencies

The code was implemented assuming to be run under c++. We have a dependency on:  

* The Cplex solver to solve the QP during homotopy initialization. Cplex can be obtained from [here](https://www.ibm.com/fr-fr/analytics/cplex-optimizer) and academic licenses are available from [here](https://www.ibm.com/support/pages/ibm-ilog-optimization-academic-initiative).
* The Armadillo linear algebra software library for the C++ programming language. Armadillo can be obtained from [here](https://www.ibm.com/fr-fr/analytics/cplex-optimizer) and he has their own dependency, described in their Readme page.

## Installing everything
sudo apt install make  
sudo apt install cmake  
sudo apt install build-essential  
sudo apt install libopenblas-dev  


cd armadillo-8.300.2/  
% delete CmakeCache.txt  
cmake .  

cd L0_optimization/Debug/  
% edit the makefile (path cplex and path armadillo)  
* -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic  
* -L/home/rbenmhenni/Bureau/these/c++/workspace/armadillo-8.300.2  
* -L/opt/ibm/ILOG/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic  

% edit /src/subdir.mk (path cplex and path armadillo)  
* -I/opt/ibm/ILOG/CPLEX_Studio128/cplex/include  
* -I/home/rbenmhenni/Bureau/these/c++/workspace/armadillo-8.300.2/include  
* -I/opt/ibm/ILOG/CPLEX_Studio128/concert/include  


cd L0_optimization/Debug/  
make


## Running the experiments
If you have setup everything according to the previous instructions, you should be able to replicate the experiments of the paper. To do so, follow the following instructions:

% rename the file L0_optimization to L0_optimization_hom and move it to the folder containing the data  
mv L0_optimization ../../data/L0_optimization_homotopy  

% run the code  
taskset -c 0 ./L0_optimization_homotopy ./N100Q100SA_SNR3_rho0.8_K3_instance1/ 1 50 3


