""" 
This file contains code to approximate solutions to the wave equation 
in two dimensions using finite differences. 

I credit the followig source with helping me to finsih the 2D code, as 
I was running into problems with the "shifting" method and then found 
this resource: 

https://beltoforion.de/en/recreational_mathematics/2d-wave-equation.php

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
    A finite differences PDE solver with an animated frontend. 
    """
    
    def __init__(self, xbounds, tbounds, fig, ax, ybounds=None, stitch=None):
        
        """ 
        Constructor for the finite differences solver. 
        """
        
        # ======================= constants and variables 
         
        # speed
        self.c = 30
        
        # determines if run in 1D mode
        if ybounds == None:
            self.is1D = True 
        else: 
            self.is1D = False
        
        # determines which 2D mode to run 
        if stitch: 
            self.stitch = True
        else: 
            self.stitch = False
        
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
        self.XC = ((self.c**2)*(self.deltaT**2))/(self.deltaX**2)
        
        # ======================= add second spatial dimension if necessary
        
        if not self.is1D:
            self.ymin = ybounds[0]
            self.ymax = ybounds[1]
            yN = 2 * (self.ymax - self.xmin)
            self.deltaY = (self.ymax - self.ymin) / yN
            self.yvals = np.linspace(self.ymin, self.ymax, yN, endpoint = False) 
            self.YC = (self.c**2)*(self.deltaT**2)/(self.deltaY**2)
        
        # ======================= initial conditions 
        
        if self.is1D: 
            IC = self.get_IC(self.is1D)
        else:
            # in the 2D case, need to determine which mode we use!
            if self.stitch:
                # this is a tuple of two paths we can stitch together
                # in stitch mode we stitch two 1D solutions together at each step!  
                IC = self.get_IC(self.is1D)
            else: 
                # this is a surface projected from two paths
                ICx, ICy = self.get_IC(self.is1D)
                IC = self.projectToSurface(ICx, ICy)
        
        # =======================
        
        # variables to track current and previous solutions curves 
        # we need both of these to compute the difference quotients
        # we initialize uprev to be the same as ucurr so our initial step is smooth 
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
        
        pillow = PillowWriter(fps=40)
        
        # plot 1D solutions 
        if self.is1D:
            self.plotSoln1D(self.ucurr)
            wave1D = self.plotSoln1D(self.ucurr)
            animation = ani.FuncAnimation(fig = self.fig, 
                                        func = self.oneStep1DVec,
                                        fargs = (wave1D, ),  
                                        frames = 1000, 
                                        interval = 5) 
            
            # uncomment this to save a gif to your local! 
            # animation.save('1DWave.gif', writer=pillow)
            
        # plot 2D solutions 
        else:
            
            if self.stitch:
                
                print("all that worked!")
                print(self.uprev[0])
                print(self.uprev[1])
                # uncomment this to save a gif to your local! 
                # animation.save('2DWaveVec.gif', writer=pillow)
                
            else: 
                wave2D = self.plotSoln2D(self.ucurr)
                animation = ani.FuncAnimation(fig = self.fig, 
                                            func = self.oneStep2DShifting,
                                            fargs = (wave2D, ),  
                                            frames = 200, 
                                            interval = 5)
                
                # uncomment this to save a gif to your local! 
                # animation.save('2DWaveVec.gif', writer=pillow)
        
        plt.show()
                
    def get_dt(self):  
               
        """
        Use the CLT condition to find a value for dt given dx (and dy). 
        """
        
        # never implemented :(
    
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
        Run one step of the PDE solver by computing a solution for a single 
        time. Each piece of the difference quotient that we need is written 
        as follows: 
        
            - ul = u(x - dx, t) 
            - xc = u(x, t) 
            - ur = u(x + dx, t) 
            - tp = u(x, t - dt) 
        
        This version loops over all of the spatial variables at each time 
        step, and is a slow brute-force approach. 
        
        This is not currently in use, but is useful for comparison. 
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
            unew[j] += self.XC*(ur + ul - 2*uc) 
        
        # update solutions and time step
        self.uprev = self.ucurr
        self.ucurr = unew
        self.tc+=1
        
    def oneStep1DVec(self, frame, wave): 
        
        """
        Run one step of the PDE solver by computing a solution for a single 
        time. Each piece of the difference quotient that we need is written 
        as follows: 
        
            - ul = u(x - dx, t) for all x in self.xvals
            - self.ucurr
            - ur = u(x + dx, y) for all x in self.xvals 
            - self.uprev
        
        This version takes a vectorized approach, and only loops over each 
        time step. At each time step, we shift the array ucurr left and 
        right by one index and feed these new vectors into the difference 
        quotient.  
        
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
        unew += self.XC*(ur + ul - 2*self.ucurr)
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
        
    def oneStep2DShifting(self, frame, wave): 
        
        """ 
        Run one step of the PDE solver by computing a solution for a single 
        time. Each piece of the difference quotient that we need is written 
        as follows: 
        
            - un = u(x, y + dy, t)
            - us = u(x, y - dy, t)
            - ue = u(x - dx, y, t)
            - uw = u(x + dx, y, t)
            - tp = u(x, y, t - dt)
        
        This takes a vectorized approach similar to oneStep1DVec. At each 
        time step, we get a 2D array for x - dx, x + dx, y - dy, and 
        y + dy by shifting the 2D array self.ucurr either north, south, 
        east or west.
        
        This is not working entirely correctly, and is currently depricated.  
        """

        # get un
        un = self.ucurr * 0
        un[0:-1,:] = self.ucurr[1:,:]
        un[-1,:] = self.ucurr[0,:]
        
        # get us 
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
        unew += self.XC*(uw + ue - 2*self.ucurr) 
        unew += self.YC*(un + us - 2*self.ucurr)
        
        # update variables and animate plot
        plt.cla()
        self.plotSoln2D(unew)
        self.tc+=1
        self.uprev = self.ucurr
        self.ucurr = unew
        
    def oneStep2DStitch(self, frame, wave):
        
        """
        Each step of the PDE solver proceeds by applying oneStep1D vec to two
        paths and then stitching them together using projectToSurface. This 
        is the solution that we plot at each step. 
        
        There is no new math or logic here, just cleverly resuing work that 
        we already did! 
        
        Note that each time we run this, both self.ucurr and self.uprev are 
        tuples of paths which contain all of the infomration needed to apply 
        the PDE! 
        """
        
        # apply oneStep to the x path 
        
        # apply oneStep to the y path 
        
        # project to a surface
        
        # clear the previous solution 
        
        # plot the new solution 
        
        # update variables 
        
        
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
    parser.add_argument('--twodimvec', action='store_true', required=False)
    parser.add_argument('--twodimstitch', action='store_true', required=False)
    args = parser.parse_args()
    
    # set up and run 1D version if specified
    if args.onedim: 
        fig, ax = plt.subplots()
        FD = FiniteDiffs([-50,50], [0,10], fig, ax)
        FD.runTest()
    
    # set up and run 2D vectorized version
    elif args.twodimvec:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        FD = FiniteDiffs([-50,50], [0,10], fig, ax, ybounds=[-50,50], stitch=False)
        FD.runTest()
    
    # setup and run the 2D path stitching version 
    elif args.twodimstitch:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        FD = FiniteDiffs([-50,50], [0,10], fig, ax, ybounds = [-50,50], stitch=True)
        FD.runTest()
    
    # alert user to run in the terminal
    else: 
        print()
        print("Please run this file in the terminal!")
        print("Specify which mode you would like to run using one of the following flags:") 
        print()
        print("     python FiniteDiff.py --onedim")
        print("     python FiniteDiff.py --twodimvec")
        print("     python FiniteDiff.py --twodimstitch")
        print()
        print("Note that you can see a higher fps version of the wave that you create by" + 
              " saving to a gif. Check the function  FiniteDiffs.runTest in the source code" +
              " to see how to do this!")
        print()

if __name__ == '__main__':
    main()