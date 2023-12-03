""" 
This file contains code to approximate solutions to the wave equation 
in two dimensions using finite differences. 

Author: Christian Lentz
"""

import argparse
import numpy as np
import math as math
import random as rand
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from matplotlib.animation import PillowWriter
from mpl_toolkits.mplot3d import Axes3D 

class FiniteDiffs: 
    
    """ 
    A finite differences PDE solver which we use in FrontEnd.py
    """
    
    def __init__(self, xbounds, tbounds, fig, ax, ybounds = None):
        
        """ 
        Constructor for the finite differences solver. 
        """
        
        # ======================= constants and variables 
         
        self.c = 20
        if ybounds == None:
            self.is1D = True 
        else: 
            self.is1D = False
        
        # ======================= spatial and temporal bounds
        
        self.xmin = xbounds[0]
        self.xmax = xbounds[1]
        self.tmin = tbounds[0]
        self.tmax = tbounds[1]
        
        # ======================= discretize spatial and temporal variables  
        
        xN = 2 * (self.xmax - self.xmin)
        self.deltaX = (self.xmax - self.xmin) / xN
        self.xvals = np.linspace(self.xmin, self.xmax, xN, endpoint = False) 
        tN = 1000 * (self.tmax - self.tmin) 
        self.deltaT = (self.tmax - self.tmin) / tN
        
        # ======================= add second spatial dimension if necessary
        
        if not self.is1D:
            self.ymin = ybounds[0]
            self.ymax = ybounds[1]
            yN = 2 * (self.ymax - self.xmin)
            self.deltaY = (self.ymax - self.ymin) / yN
            self.yvals = np.linspace(self.ymin, self.ymax, yN, endpoint = False) 
        
        # ======================= initial conditions 
        
        if self.is1D: 
            IC = self.get_IC(self.is1D)
        else:
            ICx, ICy = self.get_IC(self.is1D)
            IC = self.projectToSurface(ICx,ICy)
        
        # =======================
        
        # variables to track current and previous solutions curves 
        # we need both of these to compute the difference quotients 
        self.uprev = IC
        self.ucurr = IC
       
        # counter for time steps
        self.tc = 0
        
        # figure to plot the wave
        self.fig = fig
        self.ax = ax
        
    def runTest(self): 
        
        """
        Called in the main method of this file. Will animate solutions to 
        a randomly generated initial condition. 
        """ 
        
        # pillow = PillowWriter(fps=30)
        
        # plot 1D solutions 
        if self.is1D:
            self.plotSoln1D(self.ucurr)
            wave1D = self.plotSoln1D(self.ucurr)
            animation = ani.FuncAnimation(fig = self.fig, 
                                        func = self.oneStep1DVec,
                                        fargs = (wave1D, ),  
                                        frames = 200, 
                                        interval = 5) 
            animation.save('1DWave.gif', writer=pillow)
        # plot 2D solutions 
        else:
            wave2D = self.plotSoln2D(self.ucurr)
            animation = ani.FuncAnimation(fig = self.fig, 
                                        func = self.oneStep2D,
                                        fargs = (wave2D, ),  
                                        frames = 200, 
                                        interval = 5)
            # animation.save('2DWave.gif', writer=pillow)
        
        plt.show()
                
    def get_dt(self):  
               
        """
        Use the CLT condition to find a value for dt given dx (and dy). 
        """
    
    def get_IC(self, is1D): 
        
        """
        Randomly generate an initial condition. 
        We use random linear combintions of normal distributions. 
        """
        
        # generate random x (and y) path(s) as initial conditions
        if is1D: 
            return self.getRandHillyPath(self.xvals, self.xmax, self.xmin)
        else: 
            ICx = self.getRandHillyPath(self.xvals, self.xmax, self.xmin)
            ICy = self.getRandHillyPath(self.yvals, self.ymax, self.ymin)
            return ICx, ICy
    
    def getRandHillyPath(self, vals, maximum, minimum):
        
        """ 
        Generate a random hilly path. We use this to generate an initial 
        condition for y and x, and then combine these to get a random 
        surface as our initial wave. 
        """
        
        ic = np.zeros(vals.shape)
        for hill in range(rand.randint(8,10)):
            if hill % 2 == 0: 
                ic = ic + self.getRandSinWave(vals, maximum - minimum)
            else: 
                ic = ic + self.getRandCosWave(vals, maximum - minimum)
            
        return ic
    
    def getRandSinWave(self, vals, w):
        
        """
        Generate a random sine wave whose period alings with the boundaries.  
        """
        
        a = self.getRandScalar()
        return a*np.sin((2*np.pi/w)*vals)
    
    def getRandCosWave(self, vals, w): 
        
        """
        Get a random cosine wave whose period aligns with the boundaries. 
        """
        
        b = self.getRandScalar()
        return b*np.cos((2*np.pi/w)*vals)
    
    def getRandScalar(self):
        
        """
        Get a random real number to scale the height of the random periodic function.  
        """
           
        return rand.random()*rand.randint(-5, 5)
    
    def projectToSurface(self, xvals, yvals): 
        
        """
        From two sequences which represent points along the x and y axis, project 
        these into a surface. 
        """
        
        surf = []
        for xval in xvals:
            curr_x_arr = []
            for yval in yvals: 
                curr_x_arr.append(xval+yval)
            surf.append(curr_x_arr)
            
        return np.array(surf)
        
    def oneStep1DLoops(self): 
        
        """ 
        Run one step of the PDE solver by computing a solution for a single time.
        Each piece of the difference quotient that we need is written as follows: 
        
            - ul = u(x - dx, t) 
            - xc = u(x, t) 
            - ur = u(x + dx, t) 
            - tp = u(x, t - dt) 
        
        This version loops over all of the spatial variables at each time step, 
        and is a slow brute-force approach. 
        
        This is deprecated, but is useful for comparison. 
        """ 

        unew = 0 * self.ucurr
        
        for j in range(unew.size):
            
            # get the previous time step: u(x, T - dt)
            tp = self.uprev[j]
                
            # left boundary -- periodic boundary condition
            if j == 0:
                ul = self.ucurr[-1]
                uc = self.ucurr[j]
                ur = self.ucurr[j+1]
                
            # right boundary -- periodic boundary condition
            elif j == unew.size-1:
                ul = self.ucurr[j-1]
                uc = self.ucurr[j]
                ur = self.ucurr[0]
            else: 
                # get components of difference quotient 
                ul = self.ucurr[j-1]
                uc = self.ucurr[j]
                ur = self.ucurr[j+1]
                
            # compute difference quotient
            unew[j] = 2*uc - tp 
            unew[j] += (((self.c**2)*(self.deltaT**2))/(self.deltaX**2))*(ur + ul - 2*uc) 
        
        # update solutions and time step
        self.uprev = self.ucurr
        self.ucurr = unew
        self.tc+=1
        
    def oneStep1DVec(self, frame, wave): 
        
        """
        Run one step of the PDE solver by computing a solution for a single time.
        Each piece of the difference quotient that we need is written as follows: 
        
            - ul = u(x - dx, t) for all x in self.xvals
            - self.ucurr
            - ur = u(x + dx, y) for all x in self.xvals 
            - self.uprev
        
        This version takes a vectorized approach, and only loops over each time step.
        At each time step, we shift the array ucurr left and right by one index and 
        feed these new vectors into the difference quotient.  
        
        This saves a lot on time and space! 
        """
        
        # get the xl vector 
        ul = self.ucurr * 0 
        ul[0:-1] = self.ucurr[1:]
        ul[-1] = self.ucurr[0]
        # get the xr vector 
        ur = self.ucurr * 0
        ur[0] = self.ucurr[-1]
        ur[1:] = self.ucurr[0:-1]
        # compute difference quotient
        unew = 2*self.ucurr - self.uprev
        unew += (((self.c**2)*(self.deltaT**2))/(self.deltaX**2))*(ur + ul - 2*self.ucurr)
        # update variables and animate solns
        plt.cla()
        self.plotSoln1D(unew)
        self.tc+=1
        self.uprev = self.ucurr + 0
        self.ucurr = unew + 0
        
    def plotSoln1D(self, vals): 
        
        """
        Plot a 1D solution to the wave equation. 
        """
        
        self.ax.plot(self.xvals, vals, '-')
        self.ax.set_ylim((-10, 10))
        self.ax.set_xlabel("x axis")
        self.ax.set_ylabel("y axis")
        
    def oneStep2D(self, frame, wave): 
        
        """ 
        Run one step of the PDE solver by computing a solution for a single time.
        Each piece of the difference quotient that we need is written as follows: 
        
            - un = u(x, y + dy, t)
            - us = u(x, y - dy, t)
            - ue = u(x - dx, y, t)
            - uw = u(x + dx, y, t)
            - tp = u(x, y, t - dt)
        
        This takes a vectorized approach similar to oneStep1DVec. At each time step, 
        we get a 2D array for x - dx, x + dx, y - dy, and y + dy by shifting the 2D 
        array self.ucurr either north, south, east or west. 
        """

        # get un - BROKEN RN ?
        un = self.ucurr * 0
        un[0:-1,:] = self.ucurr[1:,:]
        un[-1,:] = self.ucurr[0,:]
        
        # get us - BROKEN RN ?
        us = self.ucurr * 0
        us[1:,:] = self.ucurr[0:-1,:]
        us[0,:] = self.ucurr[-1,:]
        
        # get ue 
        ue = self.ucurr * 0
        ue[:,0:-1] = self.ucurr[:,1:]
        ue[:,-1] = self.ucurr[:,0]
        # get uw
        uw = self.ucurr * 0
        uw[:,0] = self.ucurr[:,-1]
        uw[:,1:] = self.ucurr[0,0:-1]
        
        # compute difference quotient
        unew = 2*self.ucurr - self.uprev 
        unew += ((self.c**2)*(self.deltaT**2))/(self.deltaX**2)*(uw + ue - 2*self.ucurr) 
        
        # this does not work when you comment out the line above, why?
        unew += ((self.c**2)*(self.deltaT**2)/(self.deltaY**2))*(un + us - 2*self.ucurr)
        
        
        # update variables and animate plot
        plt.cla()
        self.plotSoln2D(unew)
        self.tc+=1
        self.uprev = self.ucurr
        self.ucurr = unew
        
    def plotSoln2D(self, vals): 
        
        """
        Plot a 2D solution to the wave equation. 
        """
        
        X, Y = np.meshgrid(self.xvals, self.yvals) 
        self.ax.plot_surface(X, Y, vals)
        self.ax.set_zlim(-10,10)
        self.ax.set_xlabel("x axis")
        self.ax.set_ylabel("y axis")
        self.ax.set_zlabel("z axis")

# ==============================================================

def main():
    
    """
    The main method here is for testing the finite differences code. 
    """ 
    
    # collect command line input
    parser = argparse.ArgumentParser(description="Runs FiniteDiff.py.")
    parser.add_argument('--onedim', action='store_true', required=False)
    args = parser.parse_args()
    
    # set up and run 1D version if specified
    if args.onedim: 
        fig, ax = plt.subplots()
        FD = FiniteDiffs([-50,50], [0,10], fig, ax)
        FD.runTest()
    
    # set up and run 2D version as defualt
    else: 
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        FD = FiniteDiffs([-50,50], [0,10], fig, ax, ybounds = [-50,50])
        FD.runTest()

if __name__ == '__main__':
    main()