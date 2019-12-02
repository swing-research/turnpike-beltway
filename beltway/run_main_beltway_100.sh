#!/bin/bash

# the program takes three arguments
# $1 the first argument is the distance file
# $2 the second argument is the option file that contains various parameters
# $3 the third argument is the recovered vector x that corresponds to the probability of the point locations
# The true sigma is 1e-05, the sigma is tuned here from the set {1e-05, 3e-05, 5e-05, 7e-05, 9e-05}

for o in 1 3 5 7 9; do
    ./main ./data/d_beltway_100_1_0.0001_1_1e-05 ./opt/main_options_cont_s_1d_100_1_res${o} ./output/x_beltway_100_1_0.0001_1_1e-05_${o}
done
