# Monte Carlo Simulations

* Perform Monte Carlo simulations to estimate the distributions of $ h_hat^T * E * h_hat$ and $sum_y h_hat^T * B_y * h_hat$.
* Compute the mutual information between the point and distance.


* This package contains code files to implement the approach described in the following paper.
```
@article{DBLP:journals/corr/abs-1804-02465,
    author    = {Shuai Huang and Ivan Dokmanic},
    title     = {Reconstructing Point Sets from Distance Distributions},
    journal   = {CoRR},
    volume    = {abs/1804.02465},
    year      = {2018},
    url       = {http://arxiv.org/abs/1804.02465},
    archivePrefix = {arXiv},
    eprint    = {1804.02465},
    timestamp = {Mon, 13 Aug 2018 16:47:54 +0200},
    biburl    = {https://dblp.org/rec/bib/journals/corr/abs-1804-02465},
    bibsource = {dblp computer science bibliography, https://dblp.org}
}
```
If you use this package and find it helpful, please cite the above paper. Thanks :smile:

## Summary
```
    ./m_file            -- This folder contains MATLAB files that are used to run the Monte Carlo simulations
    ./output            -- This folder contains the output of the main program
```

## Usage

Detailed comments on the options and experimental settings are in the MATLAB files. Please also check out the paper for more details. 

In the following example, we generate a configuration with N=100 points in the interval [0,1], the minimum pairwise distance is 5e-3, the maximum pairwise distance is 1, the quantization step lambda=1e-3.

Open `MATLAB` and type the following commands into the console:

* Step 1) Perform Gibbs sampling to estimate the distribution of $ h_hat^T * E * h_hat$
```
    >> addpath(genpath('./m_file'))
    >> compute_turnpike_distribution_hEh('6','10','2000','1','./output/turnpike_hEh_10_2000_1');    % This is for the turnpike problem
    >> compute_beltway_distribution_hEh('6','10','2000','1','./output/beltway_hEh_10_2000_1');    % This is for the beltway problem
```

* Step 2) Perform Gibbs sampling to estimate the distribution of $ sum_y h_hat^T * B_y * h_hat$
```
    >> addpath(genpath('./m_file'))
    >> compute_turnpike_distribution_hBh('6','10','2000','1','./output/turnpike_hBh_10_2000_1');    % This is for the turnpike problem
    >> compute_beltway_distribution_hBh('6','10','2000','1','./output/beltway_hBh_10_2000_1');    % This is for the beltway problem
```

* Step 3) Based on the distributions, the empirical convergence radius can be computed as described in the paper.

* Step 4) Estimate the mutual information 
```    
    >> addpath(genpath('./m_file'))
    >> compute_turnpike_mutual_information('6','10','100','./output/turnpike_mutual_info_10_100');    % This is for the turnpike problem
    >> compute_beltway_mutual_information('6','10','100','./output/beltway_mutual_info_10_100');    % This is for the beltway problem
```

