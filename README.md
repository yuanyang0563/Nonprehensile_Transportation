# Nonprehensile_Transportation

## Introduction

This project implements model predictive control algorithms for a single-arm system and a dual-arm system to perform nonprehensile transportation of a box payload.

## Features

The proposed control algorithms are computationally efficient.

Tested on computers that have 12th Gen Intel Core i7-12700*20 and install Ubuntu 20.04 LTS with PREEMPT_RT kernel. The commuations and control are within ROS Noetic.

The length of the prediction horizon is set to $N=5$.

In sumulations, only $7~ms$ on average are spent in solving the formulated optimization problems per round. Check it!

In experiments, the maximum time spent in computing the conventional controls and the predictive controls are 0.12 seconds and 0.32 seconds per round. 

## Installation

Installed **Gurobi 12.0**. Check ~/.bashrc and export environment variables **GUROBI_HOME**, **PATH** and **LD_LIBRARY_PATH** correctly.

Install the **cppzmq** package for the program to interface with the CoppeliaSim via the ZeroMQ remote API.

To install the project, run the following commands in your terminal:

```bash
# clone the repository
git clone https://github.com/yuanyang0563/Nonprehensile_Transportation.git
# navigate into the project directory
cd Nonprehensile_Transportation
# make a directory and navigate into it
mkdir build && cd build
# build the project
cmake ..
make

```

## Usage
