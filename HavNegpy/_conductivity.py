# -*- coding: utf-8 -*-
"""
Created on Mon May  9 15:56:07 2022

@author: mkolmang
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import mplcursors
import json




class Conductivity:
    """
    A class to analyze real part of complex conductivity
    
    
    """
    def __init__(self):
        """
        

        Returns
        -------
        None.

        """
 
        pass
        
     
    def create_analysis_file(self):
        """
         Creates a file to save the fit results based on the choice of fit function
         
         Provides option to use an existing file and creates a new file if not found

        Returns
        -------
        None.

        """
        res = input(str("Do you want to use an existing file to save fit results? \n eg: existing file to save HN parameters, y or n:"))
        global ana_file
        if res == 'y':
           ex_file = input("Enter the analysis_file_name:")
           try:
             f = open(ex_file)
           except FileNotFoundError as e:  
                print(f"{e}")
           else:
               if os.path.isfile(ex_file):
                  print("file exists")
         
        else:
            ex_file = input("Enter the analysis_file_name:")
            f = open(ex_file,'w')
            f.write(str('File') + '\t' + str('T') + '\t' + str('fc')+ '\t' + str('DC') + '\t'+ str('n') + '\n')
            print(f'{"file created"}',ex_file)
        ana_file = ex_file
        return()
    
    def select_range(self,x,y):
        """
        Selects the region of interest to fit data using mplcursors
        allows two clicks to select the lower and upper bound of the x-axis
        and returns the selected x and y vaues for fitting
        

        Returns
        -------
        x1 : array
            log frequency
        y1 : array
            log real part of complex conductivity

        """
        x = list(x)
        y = list(y)

        plt.figure(1)
        plt.style.use("seaborn-whitegrid")

        plt.scatter(x,y,marker='s',color='b',facecolors='none', s=100,linewidth=2)
        plt.ylabel('log ( $\sigma´$)')
        plt.xlabel('log f')
        #plt.legend()
        plt.style.use("seaborn-whitegrid")
        
        mplcursors.cursor(hover=True)
        
        zoom_ok = False
        plt.title('zoom or pan')
        
        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress()
            plt.title('press space when ready to select points')

   
            
       
        plt.title('only two clicks are allowed, select the range')
        val = plt.ginput(2)
        val.sort()
        x_min,x_max = val[0][0], val[1][0]

        
        tolerance = 0.03
        p1 = round(x_min,3)
        p2 = round(x_max,3)
        
        low_x = p1 - tolerance
        high_x = p2 + tolerance
        
        #print(low_x, high_x)
        indices = []
        indices.clear()
        for i,j in zip(x,y):
            if i<= high_x and i>=low_x :              
                k = x.index(i)              
                indices.append(k)
        
                

        a,b = indices[0], indices[-1]
        x1 = x[a:b+1]
        y1 = y[a:b+1]
        #print(x1)
        #print(y1)
        x2 = np.array(x1)
        y2 = np.array(y1)
        #print(val)
        print("x_lower_limit",x_min, "x_upper_limit",x_max)
        return x2,y2    
         
     

    



    def jonscher_function(self,x,fc,DC,n):
        """
        Jonscher power law function to fit the conductivity data

        Parameters
        ----------
        x : float
            log frequency
        fc : float
            onset frequency
        DC : float
            DC conductivity value
        n : float
            power law exponent

        Returns
        -------
        y : array
            estimated log conductivity based on supplied parameters.

        """
        o = 10**(x)
        y = DC + (np.log10(1+(o/fc)**(n)))
        return y
    

    def dump_parameters(self):
        """
        dumps the initial fit parameters for fit function as a dictionary in a json file
        to load it during curve fitting

        Returns
        -------
         None

        """
        f  = float(input("enter the fc value:"))
        DC  = float(input("enter the DC conductivity value:"))
        n  = float(input("enter the exponent:"))

        fc = 10**f
        
        
        par = {"fc": fc, "DC": DC,  "n": n}
        
        
        with open('cond.json',"w") as outfile:
            json.dump(par,outfile)
            
        with open('cond.json',"r") as openfile:
            loaded_par = json.load(openfile)
            
        print("dumped_parameters",loaded_par)
        return ()
        

    
    def fit(self,x,y):
        """
        Fits the conductivity data with choice of fit function
        The fit parameters are declared as global variables to be saved
        via save_fit function
        
        The initial fit parameters are taken from json file and the final
        fit parameters are dumped in the same json file to be used for next
        iteration.

        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log conductivity.

        Returns
        -------
        None.

        """


        x1 = np.array(x)
        y1= np.array(y)

        global popt1
        
        global  fit_fc,fit_DC,fit_n
        
        plt.figure()
        
        try:
            
          open('cond.json')

        except  FileNotFoundError as e:  
            print(f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:
            with open('cond.json',"r") as openfile:
               loaded_par = json.load(openfile)
            
           
               
    
            p0 = [loaded_par['fc'], loaded_par['DC'],loaded_par['n']]
        
            cond = self.jonscher_function
            popt1, pcov1 = curve_fit(cond, x1, y1, p0, bounds = ((1e-6,-15,0), (1e8,6,1)))
            yfit= cond(x1,*popt1)
            
            plt.scatter(x1,y1,marker='s',color='b',facecolors='none',label='data',s=100,linewidth=2)
            plt.plot(x1,yfit,'m--', label='fit',linewidth=2)
            plt.xlabel('log ( f [Hz])')
            plt.ylabel('log ( $\sigma´$)')
            plt.legend()
            plt.show()
            fit_fc,fit_DC,fit_n = popt1[:]
            
            fit_par = {"fc": fit_fc, "DC": fit_DC,  "n": fit_n}
    
            
            print("fit parameters:\n", popt1)
            print(f' DC cond = {fit_DC:.03f}')
            
            with open('cond.json',"w") as outfile:
                json.dump(fit_par,outfile)
                
            with open('cond.json',"r") as openfile:
                loaded_par = json.load(openfile)
    
            print("fit parameters dumped for next iteration",loaded_par)
            

        return()

    def save_fit(self,T):
        """
        saves the fit parameters of fit function in  a file, 
        the file must be created via create_analysis_file function
       

        Parameters
        ----------
        T : float
            Temperature,or can also be an integer
            that corresponds to a file number during analysis.
        Returns
        -------
        None.

        """
        

        res_file = open(ana_file,'a')
        res_file.write( f'{T}' + '\t' + f'{fit_fc:.03f}' + '\t' + f'{fit_DC:.03f}' + '\t' + f'{fit_n:.03f}'  +"\n")
        return ()  







