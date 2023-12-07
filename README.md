# Numerical Techniques for the Wave Equation
## Christian Lentz, Macalester College Department of Mathematics, Statistics and Computer Science
### Mathematical Modeling, Fall 2023

#### Description 

This project is an investigation of numerical techniques for the wave equation, and in particular finite differences. This project was built over November and December of 2023 for the Methematical Modeling capstone at Macalester College, and includes code to approximate solution to the wave equation using four different implementations of finite differences. Also included in this repository is a brief paper describing the finite differences process in generality. In the paper, the discussion and the implementation of finite differences is used to motivate a discussion of the Fourier Spectral method, which is not implemented in this project. 

#### Running this Project

First, open the terminal of your choice and navigate to the root directory of this repository. There are four different modes to run this finite differences program, corresponding to the four different implementations that 
we include. First, run **python FiniteFiff.py**. If this does not work, you may need to run **python3 FiniteFiff.py**. This will print a message in the console, telling you that you must specify a flag in order to run the project. This flag corresponds to the implementation of finte differences that we want to use to run the program. 

**--onedimvec :** this flag specifies that we wish to plot 1D solutions to the wave equation using an efficient vectorized version (thanks, NumPy). 

**--onedinloop :** this flag specifies that we wish to plot 1D solutions to the wave equation using a brute forcefunction that explicitly lopps over both temporal and spatial variables. The vectorized version only loops over the temporal variable. 

**--twodimvec :** this flag specifies that we wish to plot 2D solutions to the wave equation using the efficient vectrized approach (thanks once again, NumPy)

**--twodimstitch :** this flag specifies that we wish to plot 2D solutions to the wave equation by storing two paths, and at each discrete time step of our computation, apply the PDE to these two paths and "stitch" then togethr to form a surface. This is more accurate than our implementation of the 2D vectorized version, and does less work in temrs of applying the PDE. However, we must loop over the entire computational grid in order to stitch together our surface. 

Note that each implementation assumes periodic boundary conditions, and starts with a randomly generated initial condition / function that is peeriodic with the boundries of the computational grid. 

#### Bugs 

There is current a small bug with the two-dimensional vectorized finite differeneces, which causes solutions to diverge after some time, starting with the boundaries. 
