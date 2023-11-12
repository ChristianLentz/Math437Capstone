""" 
This file contains code to approximate solutions to the wave equation 
using finite differences. 

Author: Christian Lentz
"""

import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.animation as ani

class FiniteDiffs: 
    
    """ 
    A finite differences PDE solver which we use in FrontEnd.py
    """
    
    def __init__(self, xbounds, ybounds, tbounds, fig, ax):
        
        """ 
        Constructor for the finite differences solver. 
        """
        
        # constants 
        self.c = 0.06
        
        # spatial and temporal bounds
        self.xmin = xbounds[0]
        self.xmax = xbounds[1]
        self.ymin = ybounds[0]
        self.ymax = ybounds[1]
        self.tmin = tbounds[0]
        self.tmax = tbounds[1]
        
        # discretize spatial and temporal variables  
        xN = 50
        self.deltaX = (self.xmax - self.xmin) / xN
        self.xvals = np.linspace(self.xmin, self.xmax, xN) 
        yN = 50
        self.deltaY = (self.ymax - self.ymin) / yN
        self.yvals = np.linspace(self.ymin, self.ymin, yN) 
        tN = 50 
        self.deltaT = (self.tmax - self.tmin) / tN

        # initial conditions
        self.ICx = np.exp(-((self.xvals - self.xmax/2)/15)**2)
        self.ICy = np.exp(-((self.yvals - self.ymax/2)/15)**2)
        
        # array to hold previous and current solutions
        self.uprev = None             # u(x,y,T - dt)
        self.ucurrent = self.ICx      # u(x,y,T)
        
        # counter for time steps
        self.tc = 0
        
        # figure to plot the function 
        self.fig = fig
        self.ax = ax
        
        # plot the initial solution curve 
        self.ax.plot(self.ICx)
        
        # # animate the solution curves 
        # animation = ani.FuncAnimation(fig = self.fig, 
        #                                 func = FD.oneStep(),
        #                                 fargs = (img, ),  
        #                                 frames = 10, 
        #                                 interval = 50)
        
    def get_dt(self):         
        """
        Use the CLT condition to find a value for dt given dx and dy. 
        """
        
        # this may not be included, we'll see
    
    def get_IC(self): 
        
        """
        Randomly generate an initial condition. 
        We use random linear combintions of approximately normal distributions 
        which are centered on random x and y coordinates. 
        """
        
        # this might need to be included in the FrontEnd class? 
        
    def oneStep(self): 
        
        """ 
        Run one step of the PDE solver by computing a solution for a single time.
        Each piece of the difference quotient that we need is written as follows: 
        
            - xl = x - dx = the previous x coordinate
            - xc = x = the current x coordinate
            - xr = x + dx = the next x coordinate 
            
            - yl = y - dy = the previous y coordinate 
            - yc = y = the current y coordinate 
            - yr = y + dy = the next y coordinate
            
            - tp = t - dt = the previous time step
        """

        unew = 0 * self.ucurrent
        for j in range(unew.size):
            
            # get tp if current time step is not the initial 
            if self.tc != 0: 
                tp = self.uprev[j]
            else: 
                tp = 0
                
            # left boundary 
            if j == 0:
                unew[j] = 0
                
            # right boundary 
            elif j == unew.size-1:
                # get components of difference quotient 
                xl = self.ucurrent[j-1]
                xc = self.ucurrent[j]
                # compute difference quotient
                unew[j] = 2*xc - tp
                unew[j] += (((self.c**2)*(self.deltaT**2))/(self.deltaX**2))*(0 + xl - 2*xc)
            
            else: 
                # get components of difference quotient 
                xl = self.ucurrent[j-1]
                xc = self.ucurrent[j]
                xr = self.ucurrent[j+1]
                # compute difference quotient
                unew[j] = 2*xc - tp 
                unew[j] += (((self.c**2)*(self.deltaT**2))/(self.deltaX**2))*(xr + xl - 2*xc) 
        
        # update our solutions 
        self.uprev = self.ucurrent
        self.ucurrent = unew
        self.tc+=1
         
        # plot the new current solution 
        self.ax.plot(self.ucurrent, '-')
        
# ==============================================================

def main():
    
    """
    The main method here is for testing the finite differences code. 
    """ 
    
    fig,ax = plt.subplots()
    FD = FiniteDiffs([-4,50], [-4, 50], [0,25], fig, ax)
    
    for tj in range(25): 
        FD.oneStep()
        
    plt.show()

if __name__ == '__main__':
    main()