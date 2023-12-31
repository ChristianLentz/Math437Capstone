""" 
This file contains code to approximate solutions to the wave equation 
in one and two dimensions using finite differences. 

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
from Helpers import helper

class FiniteDiffs: 
    
    """ 
    A finite differences PDE solver with an animated frontend. 
    """
    
    def __init__(self, xbounds, tbounds, fig, ax, ybounds=None, stitch=None, loop=None):
        
        """ 
        Constructor for the finite differences solver. 
        """
        
        # ======================= constants and variables 
        
        # helper methods
        self.helper = helper()
         
        # speed
        self.c = 200
        
        # determines if run in 1D mode
        if ybounds == None:
            self.is1D = True 
        else: 
            self.is1D = False
        
        # detemines which 1D mode to run 
        if loop: 
            self.loop = True
        else: 
            self.loop = False
        
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
            IC = self.helper.get_IC(self.is1D, 
                self.xvals, self.xmin, self.xmax, 
                None, None, None) 
        else:
            # in the 2D case, need to determine which mode we use!
            if self.stitch:
                # IC is a tuple of two paths  
                IC = self.helper.get_IC(self.is1D, 
                self.xvals, self.xmin, self.xmax, 
                self.yvals, self.ymin, self.ymax) 
            else: 
                # IC is a surface
                ICx, ICy = self.helper.get_IC(self.is1D, 
                self.xvals, self.xmin, self.xmax, 
                self.yvals, self.ymin, self.ymax)
                IC = self.helper.projectToSurface(ICx, ICy)
        
        # =======================
        
        # variables to track current and previous solutions curves 
        # we need both of these to compute the difference quotients
        # initialize uprev to be the same as ucurr so our initial step is smooth 
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
        
        pillow = PillowWriter(fps=50)
    
        # plot 1D solutions 
        if self.is1D:
            
            if self.loop: 
                print()
                print("Now plotting 1D solutions with the looping function!")
                print()
                wave1D = self.helper.plotSoln1D(self.xvals, self.ucurr, self.ax)
                animation = ani.FuncAnimation(fig = self.fig, 
                                        func = self.oneStep1DLoops,
                                        fargs = (wave1D, ),  
                                        frames = 1000, 
                                        interval = 5)
                
                # uncomment this to save a gif to your local! 
                # animation.save('1DWaveLoop.gif', writer=pillow)
                
            else: 
                print()
                print("Now plotting 1D solutions with vectorized function!")
                print()
                wave1D = self.helper.plotSoln1D(self.xvals, self.ucurr, self.ax)
                animation = ani.FuncAnimation(fig = self.fig, 
                                            func = self.oneStep1DVec,
                                            fargs = (wave1D, ),  
                                            frames = 1000, 
                                            interval = 5) 
                
                # uncomment this to save a gif to your local! 
                # animation.save('1DWaveVec.gif', writer=pillow)
            
        # plot 2D solutions
        else:
            
            if self.stitch:
                print()
                print("Now plotting 2D solutions with stitching function!")
                print()
                initWave = self.helper.projectToSurface(self.ucurr[0], self.ucurr[1])
                wave2D = self.helper.plotSoln2D(self.xvals, self.yvals, initWave, self.ax)
                animation = ani.FuncAnimation(fig = self.fig, 
                                            func = self.oneStep2DStitch,
                                            fargs = (wave2D, ),  
                                            frames = 500, 
                                            interval = 5)
                
                # uncomment this to save a gif to your local! 
                # animation.save('2DWaveStitch.gif', writer=pillow)
                
            else: 
                print()
                print("Now plotting 2D solutions with vectorized function!")
                print()
                wave2D = self.helper.plotSoln2D(self.xvals, self.yvals, self.ucurr, self.ax)
                animation = ani.FuncAnimation(fig = self.fig, 
                                            func = self.oneStep2DShifting,
                                            fargs = (wave2D, ),  
                                            frames = 500, 
                                            interval = 5)
                
                # uncomment this to save a gif to your local! 
                # animation.save('2DWaveVec.gif', writer=pillow)
        
        plt.show()
        
    def oneStep1DLoops(self, frame, wave): 
        
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
        plt.cla()
        self.helper.plotSoln1D(self.xvals, unew, self.ax)
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
        self.helper.plotSoln1D(self.xvals, unew, self.ax)
        self.tc+=1
        self.uprev = self.ucurr + 0
        self.ucurr = unew + 0
        
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
        uw[:,1:] = self.ucurr[:,0:-1]              # this is where the error was: [0,0:-1] !!
        # compute difference quotient
        unew = 2*self.ucurr - self.uprev 
        unew += self.XC*(uw + ue - 2*self.ucurr) 
        unew += self.YC*(un + us - 2*self.ucurr)
        # update variables and animate plot
        plt.cla()
        self.helper.plotSoln2D(self.xvals, self.yvals, unew, self.ax)
        self.tc+=1
        self.uprev = self.ucurr
        self.ucurr = unew
        
    def oneStep2DStitch(self, frame, wave):
        
        """
        Each step of the PDE solver proceeds by applying a vectorized version 
        of 1D finite differences to two paths, and then stitching these paths 
        together to form a surface using project to surface. 
        
        This stitchig together process is linear in the number of points in the 
        computational grid. 
        
        There is no new math or logic here, just cleverly resuing work that 
        we already did! 
        
        The components needed for the computation are defined as:
        
        - xl = u(x - dx, t) for all x in self.xvals
        - xr = u(x + dx, y) for all x in self.xvals 
        - yl = u(x, y - dy t) for all x in self.yvals
        - yr = u(x, y + dy, t) for all x in self.yvals 
        
        Note that each time we run this, both self.ucurr and self.uprev are 
        tuples of paths which contain all of the infomration needed to apply 
        the PDE! 
        """
        
        # get our current and prev soln soln
        curr_xpath = self.ucurr[0]
        prev_xpath = self.uprev[0]
        curr_ypath = self.ucurr[1]
        prev_ypath = self.uprev[1]
        # apply oneStep helper to curr paths
        xl, xr = self.oneStep2DStitchHelper(curr_xpath)
        yl, yr = self.oneStep2DStitchHelper(curr_ypath)
        # apply the PDE to each shifted path
        unewx = 2*curr_xpath - prev_xpath + self.XC*(xr + xl - 2*curr_xpath)
        unewy = 2*curr_ypath - prev_ypath + self.YC*(yr + yl - 2*curr_ypath)
        # project to a surface, plot update variables
        unew2D = self.helper.projectToSurface(unewx, unewy)
        plt.cla()
        self.helper.plotSoln2D(self.xvals, self.yvals, unew2D, self.ax)
        self.tc+=1
        self.uprev = self.ucurr
        self.ucurr = unewx, unewy
        
    def oneStep2DStitchHelper(self, path): 
        
        """
        Call vectorized 1D finite differencces on a single path. We only return the arrays that
        have been shifted over the boundry conditions, and then apply the wave PDE on the 
        object that we have returned.  
        """
            
        # get the "left-shifted" path  
        ul = path * 0 
        ul[0:-1] = path[1:]
        ul[-1] = path[0]
        # get the "right-shifted" path
        ur = path * 0
        ur[0] = path[-1]
        ur[1:] = path[0:-1]
        
        return ul, ur

# ==============================================================

def main():
    
    """
    The main method here is for testing the finite differences code. 
    """ 
    
    # collect command line input
    parser = argparse.ArgumentParser(description="Runs FiniteDiff.py.")
    parser.add_argument('--onedimvec', action='store_true', required=False)
    parser.add_argument('--onedimloop', action='store_true', required=False)
    parser.add_argument('--twodimvec', action='store_true', required=False)
    parser.add_argument('--twodimstitch', action='store_true', required=False)
    args = parser.parse_args()
    
    # set up and run 1D vectoried version if specified
    if args.onedimvec: 
        fig, ax = plt.subplots()
        FD = FiniteDiffs([-50,50], [0,10], fig, ax, loop=False)
        FD.runTest()
    
    # set up and run 1D loops version if specified
    elif args.onedimloop: 
        fig, ax = plt.subplots()
        FD = FiniteDiffs([-50,50], [0,10], fig, ax, loop=True)
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
        print("     python FiniteDiff.py --onedimvec")
        print("     python FiniteDiff.py --onedimloop")
        print("     python FiniteDiff.py --twodimvec")
        print("     python FiniteDiff.py --twodimstitch")
        print()
        print("Note that you can see a higher fps version of the wave that you create by" + 
              " saving to a gif. Check the function FiniteDiffs.runTest in the source code" +
              " to see how to do this!")
        print()

if __name__ == '__main__':
    main()