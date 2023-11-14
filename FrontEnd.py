""" 
This file contains a front end GUI using matplotlib which animates solutions 
to the wave equation with either finite differences or spectral methods, and 
adds swimmer agents to navigate the toroidial wave environment. 

Author: Christian Lentz 
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from FiniteDiff import FiniteDiffs

class FrontEnd: 
    
    def __init__(self): 
        
        """
        Initialize the front end class. 
        """
        
        # set the size of the grid 
        self.xSize = 200
        self.ySize = 200
        
        # create the GUI grid 
        
        # randomly generate initial wave
        
        # the finite diff solver
        
        # the spectral solver
        