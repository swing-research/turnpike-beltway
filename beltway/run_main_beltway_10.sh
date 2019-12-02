#!/bin/bash

# the program takes three arguments
# $1 the first argument is the distance file
# $2 the second argument is the option file that contains various parameters
# $3 the third argument is the recovered vector x that corresponds to the probability of the point locations
# The true sigma is 0.001, the sigma is tuned here from the set {0.001, 0.003, 0.005, 0.007, 0.009}

for o in 1 3 5 7 9; do
    ./main ./data/d_beltway_10_1_0.01_1_0.001 ./opt/main_options_cont_s_1d_10_1_res${o} ./output/x_beltway_10_1_0.01_1_0.001_${o}
done
