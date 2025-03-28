# Nonprehensile_Transportation

## Introduction

This project implements model predictive control algorithms for a single-arm system and a dual-arm system to perform nonprehensile transportation of a box payload.

The principles of the model predictive control algorithms are explained here: <a href="https://www.overleaf.com/read/svsbmbmrbbzt#efa959" target="_blank">Enjoy it!</a>.

## Features

The proposed control algorithms are computationally efficient.

The length of the prediction horizon is set to $N=5$.

In sumulations, only $7~ms$ on average are spent in solving the formulated optimization problems per round. Check it!

## Installation

Installed **Gurobi 12.0**. Check ~/.bashrc and export environment variables **GUROBI_HOME**, **PATH** and **LD_LIBRARY_PATH** correctly.

Install the **cppzmq** package for the program to interface with the CoppeliaSim via the ZeroMQ remote API.

## Usage
