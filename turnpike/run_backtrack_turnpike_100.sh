#!/bin/bash

# the program takes four arguments
# $1 the first argument is the distance file. Note that the distances should be quantized to integers
# $2 the second argument is the recovered point locations
# $3 the third argument is the choice of the backtracking algorithm. If set to 0, the algorithm adopts a depth-first search. If set to 1, the algorithm adopts the first type breadth-first search. If set to other integers, the algorithm adopts a second type breadth-first search
# $4 the fourth argument is the error range, usually set to five times of the standandard deviation "sigma" of the Gaussian noise. It is set to 5*sigma/dist_unit, where dist_unit is the quantization step. 
# The true sigma is 1e-05, the sigma is tuned here from the set {1e-05, 3e-05, 5e-05, 7e-05, 9e-05}

./backtrack/PDP ./data/d_bt_turnpike_100_1_0.0001_1_1e-05 ./output/x_bt_turnpike_100_1_0.0001_1_1e-05_5 0 5
./backtrack/PDP ./data/d_bt_turnpike_100_1_0.0001_1_1e-05 ./output/x_bt_turnpike_100_1_0.0001_1_1e-05_15 0 15
./backtrack/PDP ./data/d_bt_turnpike_100_1_0.0001_1_1e-05 ./output/x_bt_turnpike_100_1_0.0001_1_1e-05_25 0 25
./backtrack/PDP ./data/d_bt_turnpike_100_1_0.0001_1_1e-05 ./output/x_bt_turnpike_100_1_0.0001_1_1e-05_35 0 35
./backtrack/PDP ./data/d_bt_turnpike_100_1_0.0001_1_1e-05 ./output/x_bt_turnpike_100_1_0.0001_1_1e-05_45 0 45
