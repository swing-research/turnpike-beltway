# Solving the partial digestion problem

* Perform partial digestion on the E. Coli K12 MG1655 Genome data from GenBank assembly, which is a nucleotide sequence of length=4,641,652.
* We try to reconstruction the locations of two sequences: `CCC|GGG` and `G|GATCC`, where `|` indicates the location of the restriction site in each sequence.
![partial_digestion](partial_digestion.png){width=30%}

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
    ./m_file            -- This folder contains MATLAB files that are used to create the datset, and to extract point locations
    ./src               -- This folder contains C++ files that are used to perform the reconstruction using the proposed approach
    ./opt               -- This folder contains option files for the main program
    ./data              -- This folder contains the dataset
    ./output            -- This folder contains the output of the main program
```

## Usage

You can follow the following steps to run the program.

* Step 1) Create the datasets. Open `MATLAB` and type the following commands into the console:
```
    >> addpath(genpath('./m_file'))
    >> pd_pairwise_distances;    % Find the locations of each DNA sequence, and compute the pairwise distances in each case
```

* Step 2) Compile the C++ program. For the noiseless recovery, we only need to consider the non-zero distance measurements.
```
    g++ -std=c++11 -fopenmp -O3 -o main ./src/main.cpp ./src/DataReader.cpp ./src/uDGP.cpp ./src/SED.cpp ./src/global.cpp -I ../turnpike/eigen-eigen-3.3.7 -pthread
```

* Step 3) Recover the locations of the sequence `CCC|GGG`.
```    
    ./main ./data/enzyme_pw_distance_1 ./opt/main_options_cont_s_1d_495_1_res0 ./output/enzyme_pos_1_est
```
Take the option file `opt/main_options_cont_s_1d_10_1_res1` for example, it contains the following options:
```
    num_thread      --  The number of threads used for computation
    num_smp         --  The number of points to be recovered
    init_type       --  The initialization type: "1" for spectral initialization, "2" for random initialization, "3" for uniform initialization
    max_ite         --  The maximum number of iterations for the projected gradient descent
    max_sg_ite      --  The maximum number of iterations when computing the spectral initializer using the power method
    min_space_unit  --  The quantization step chosen to discretize the 1D domain space, usually set to a tenth of the minimum distance between two differnet points
    sigma           --  The standard deviation of the Gaussian distribution we use to convolve with the distance measurements.
    tau             --  The maximum deviation of a noisy distance from the noiseless distance, also known as the error range. It should be at least 3 times of sigma. Setting it to be 5 times of sigma would contain 99.9999% of the error.
    perturb_std     --  The standandard deviation of a random Gaussian number that is used to perturb the spectral initializer a bit, usually set to some number around 0.01
    sg_tol          --  The convergence threshold when computing the spectral initializer
    step_ori        --  The initialial step size of the projected gradient descent method
    bkt_rate        --  The backtracking rate used while search for a proper step size, usually set to some number between 0.9 and 1
    step_thd        --  The smallest step size of the projected gradient descent
    cvg_thd         --  The convergence threshold of the projected gradient descent
```
<br/><br/>
Recover the locations of the sequence `G|GATCC`.
```
    ./main ./data/enzyme_pw_distance_2 ./opt/main_options_cont_s_1d_512_2_res0 ./output/enzyme_pos_2_est
```
The point locations are saved in the directory `./output`.


