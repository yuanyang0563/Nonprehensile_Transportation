# Nonprehensile_Transportation

## Introduction

This repository implements model predictive control algorithms for distributed dual-arm nonprehensile object transportation.

## Features

The proposed control algorithms are computationally efficient.

Tested on desktop computers with Intel Core i7-12700 (12 cores, 20 threads, base frequency 2.1GHz) and running Ubuntu 20.04 LTS with PREEMPT_RT kernel.

The length of the prediction horizon is set to $N=5$. The average computation time required per control cycle is $8.37$ ms (centralized) and $5.40$ ms (distributed).

## Installation

Install Gurobi 12.0.

Install the Robot Operating System Noetic.

Install the Visual Servoing Platform (ViSP). 

Install CoppeliaSim (for simulations only) and Universal_Robots_ROS_Driver (for experiments only).

To install the project, run the following commands in your terminal:

```bash

git clone https://github.com/yuanyang0563/Nonprehensile_Transportation.git

cd Nonprehensile_Transportation/simulations 
mkdir build && cd build
ccmake ..  # type c twice and g once
cmake ..
make

cd Nonprehensile_Transportation/experiments/***
mkdir build && cd build
ccmake ..  # type c twice and g once
cmake ..
make

```
