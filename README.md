# Hartmann Flow

This problem shows a simple 2D incompressible Hartmann flow using the vector potential form of MHD. The liquid metal flows in the x direction, in a channel with walls at y = 1 and -1. Because the magnetic field in the direction transverse to the channel is zero, we need only the z component of the vector potential A. 

## Problem Description

As is, the problem has four parameters, Ha, Re, Rm, and Pi:

 * Ha is the Hartmann number
 * Re is the Reynolds number
 * Rm is the magnetic Reynolds number
 * Pi is the background pressure gradient in x that drives the channel flow.

The problem is periodic in x and wall-bounded in y, with nx = ny = 64. It starts with the background pressure gradient Pi=0, and all fields equal to zero. The background pressure gradient then ramps up to Pi with a time constant of tau = 1. By default, the script runs to time = 10.

## How to run
0. `cd` to `hartmann_flow/python/` and make the `scratch` directory,
```
$ cd path/to/hartmann_flow/python/
$ mkdir scratch
```
1. Run the script `hartmann.py`. On a single core of my i5-powered mini desktop, this takes about 180 seconds.
```
$ python3 hartmann.py
```
2. Run the script `plot_vx_yprof.py`, which will plot the x-averaged x velocity field as a function of y.
```
$ python3 plot_vx_yprof.py scratch
```
## Suggested Exercises