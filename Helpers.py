""" 
This file contains code helper functions for this project, such 
as plotting solutions and generating initial conditions.  

Author: Christian Lentz
"""

import numpy as np
import random as rand

class helper: 
    
    def __init__(self): 
        
        """
        Construct the helper class. 
        """
        
    def get_dt(self):  
               
        """
        Use the CLT condition to find a value for dt given dx (and dy). 
        """
        
        # never implemented :(
    
    def get_IC(self, is1D, xvals, xmin, xmax, yvals, ymin, ymax): 
        
        """
        Randomly generate an initial condition. 
        We use random linear combintions of normal distributions. 
        """
        
        if is1D: 
            return self.getRandHillyPath(xvals, xmax, xmin)
        else: 
            ICx = self.getRandHillyPath(xvals, xmax, xmin)
            ICy = self.getRandHillyPath(yvals, ymax, ymin)
            return ICx, ICy
    
    def getRandHillyPath(self, vals, maximum, minimum):
        
        """ 
        Generate a random hilly path. We use this to generate an initial 
        condition for y and x, and then combine these to get a random 
        surface as our initial wave. 
        """
        
        ic = np.zeros(vals.shape)
        for hill in range(rand.randint(15,20)):
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
    
    def plotSoln1D(self, xvals, vals, ax): 
        
        """
        Plot a 1D solution to the wave equation. 
        """
        
        ax.plot(xvals, vals, '-')
        ax.set_ylim((-10, 10))
        ax.set_xlabel("x axis")
        ax.set_ylabel("y axis")
        
    def plotSoln2D(self, xvals, yvals, vals, ax): 
        
        """
        Plot a 2D solution to the wave equation. 
        """
        
        X, Y = np.meshgrid(xvals, yvals) 
        ax.plot_surface(X, Y, vals)
        ax.set_zlim(-10,10)
        ax.set_xlabel("x axis")
        ax.set_ylabel("y axis")
        ax.set_zlabel("z axis")