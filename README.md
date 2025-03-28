# Nonprehensile_Transportation

## Introduction

This project implements model predictive control algorithms for a single-arm system and a dual-arm system to perform nonprehensile transportation of a box payload. The principles of the model predictive control algorithms are explained here:

[ Enjoy it! ] (https://www.overleaf.com/read/svsbmbmrbbzt#efa959)

## Installation

Make sure that you have installed Gurobi 12.0. Check ~/.bashrc and see if you have exported environment variables GUROBI_HOME, PATH and LD_LIBRARY_PATH correctly.

The program interfaces with the CoppeliaSim via ZeroMQ remote API, which requires the jsoncons and cppzmq packages. All files about the ZeroMQ remote API and the jsoncons package have been included in the include directory. You need to install the cppzmq package.
