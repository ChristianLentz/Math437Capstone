# Numerical Techniques for the Wave Equation
## Christian Lentz, Macalester College Department of Mathematics, Statistics and Computer Science
### Mathematical Modeling, Fall 2023

#### Description 

This project is an investigation of numerical techniques for the wave equation, and in particular finite differences. This project was built over November and December of 2023 for the Methematical Modeling capstone at Macalester College, and includes code to approximate solution to the wave equation using finite differences. Also included in this repository is a brief paper describing the finite differences process in generality. In the paper, the discussion and the implementation of finite differences is used to motivate a discussion of the Fourier Spectral method, which is not implemented in this project. 

#### Running this Project

Navigate to the root directory of this repository on your local machine, open the terminal and run **python FiniteFiff.py**. If this does not work, you may need to run **python3 FiniteFiff.py**. This will run the finite differences code for the wave equation in two spatial dimensions. If you wish to see a one dimensional model, then run **python FiniteFiff.py --onedim**. This code utilizes numpy to take a vectorized approach, thus we are able to avoid looping over spatial deimensions and only loop over the temporal variable. Further, we assume periodic boundary conditions, and randomly generate initial conditions which are periodic with the boundaries. 

#### Bugs 

There is current a small bug with the two-dimensional finite differeneces, which causes solutions to diverge after some time, starting with the boundaries. 
