# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:46:06 2022

@author: mkolmang
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import mplcursors
import json


class HN():
    """
    A class to analyze dielectric loss data with HN function
    
    Fit functions include single, double, HN with conductivity, and HN with flank
    
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
        res = input(str(
            "Do you want to use an existing file to save fit results? \n eg: existing file to save HN parameters, y or n:"))
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
            func_name = int(input(
                "Choose the fit function\n 1 -- HN, 2 -- HN with cond, 3 -- HN-flank, 4 -- double HN, 5 -- double HN with cond:"))
            ex_file = input("Enter the analysis_file_name:")
            f = open(ex_file, 'w')
            if func_name == 1 or func_name == 2:
                 f.write(f'{"Temperature"}' + '\t' + f'{"beta"}' + '\t' + f'{"gamma"}' + '\t' +
                         f'{"deps"}' + '\t' + f'{"log (fmax)"}' + '\t' + f'{"cond"}' + '\t' + f'{"s"}' + '\n')
            elif func_name == 3:
                 f.write(f'{"Temperature"}' + '\t' + f'{"beta"}' + '\t' + f'{"gamma"}' + '\t' +
                         f'{"deps"}' + '\t' + f'{"log (fmax)"}' + '\t' + f'{"flank"}' + '\t' + f'{"s"}' + '\n')
            elif func_name == 4 or func_name == 5:
                 f.write(f'{"Temperature"}' + '\t' + f'{"beta1"}' + '\t' + f'{"gamma1"}' + '\t' + f'{"deps1"}' + '\t' + f'{"log (fmax1)"}' + '\t' + \
                         f'{"beta2"}' + '\t' + f'{"gamma2"}' + '\t' + f'{"deps2"}' + '\t' + f'{"log (fmax2)"}' + '\t' + f'{"cond"}' + '\t' + f'{"s"}' + '\n')
        print(f'{"file did not exist, created"}', ex_file)
        ana_file = ex_file
        return()

    def select_range(self, x, y):
        """
        
        Selects the region of interest to fit data using mplcursors
        allows two clicks to select the lower and upper bound of the x-axis
        and returns the selected x and y vaues for fitting


        Returns
        -------
        x1 : array
            log frequency
        y1 : array
            log dielectric loss

        """
        x = list(x)
        y = list(y)

        plt.figure(1)
        plt.style.use("seaborn-whitegrid")

        plt.scatter(x, y, marker='s', color='b',
                    facecolors='none', s=100, linewidth=2)
        plt.ylabel('log $\epsilon$"')
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
        x_min, x_max = val[0][0], val[1][0]

        tolerance = 0.03
        p1 = round(x_min, 3)
        p2 = round(x_max, 3)

        low_x = p1 - tolerance
        high_x = p2 + tolerance

        #print(low_x, high_x)
        indices = []
        indices.clear()
        for i, j in zip(x, y):
            if i <= high_x and i >= low_x:
                k = x.index(i)
                indices.append(k)

        a, b = indices[0], indices[-1]
        x1 = x[a:b+1]
        y1 = y[a:b+1]
        #print(x1)
        #print(y1)
        x2 = np.array(x1)
        y2 = np.array(y1)
        #print(val)
        print("x_lower_limit", x_min, "x_upper_limit", x_max)
        return x2, y2

    def hn(self, x, b, g, fm, deps):
        """
        HN fit function to fit single dielectric loss peak


        Parameters
        ----------
        x : float
            frequency.
        b : float
            symmetric fractional parameter.
        g : float
            asymmetric fractional parameter.
        fm : float
            maximum frequency of the peak.
        deps : float
            dielectric strength.

        Returns
        -------
        y : array
            estimated log dielectric loss based on the supplied parameters.

        """
        f = 10**(x)

        ff = f/fm
        ffb = ff**(b)
        n = (np.pi)/2

        xc = np.cos(b*n)
        xs = np.sin(b*n)
        r = 1 + 2*ffb*xc + ffb**(2)
        r_w = (1/r)**(0.5*g)

        p2 = 1/(ffb)
        phi = np.arctan((xs/(p2+xc)))
        e_w = deps*r_w*(np.sin(g*phi))

        y = np.log10(e_w)
        return y

    def hn_cond(self, x, b, g, fm, deps, cond, s):
        """
        HN fit function to fit single dielectric loss peak together with conductivity
        

        Parameters
        ----------
        x : float
            frequency.
        b : float
            symmetric fractional parameter.
        g : float
            asymmetric fractional parameter.
        fm : float
            maximum frequency of the loss peak.
        deps : float
           dielectric strength.
        cond : float
            DC conductivity contribution.
        s : float
            power law exponent.

        Returns
        -------
        y : array
            estimated log dielectric loss based on the supplied parameters.

        """
        f = 10**(x)

        ff = f/fm
        ffb = ff**(b)
        n = (np.pi)/2

        xc = np.cos(b*n)
        xs = np.sin(b*n)
        r = 1 + 2*ffb*xc + ffb**(2)
        r_w = (1/r)**(0.5*g)

        p2 = 1/(ffb)
        phi = np.arctan((xs/(p2+xc)))
        e_w = deps*r_w*(np.sin(g*phi))
        fs = f**(s)
        condb = cond/fs
        y = np.log10(e_w+condb)
        return y

    def hn_flank(self, x, b, g, fm, deps, N, s):
        """
        HN fit function to fit single dielectric loss peak together with a high frequency flank

        Parameters
        ----------
        x : float
            frequency.
        b : float
            symmetric fractional parameter.
        g : float
            asymmetric fractional parameter.
        fm : float
            maximum frequency of the loss peak.
        deps : float
           dielectric strength.
        N : float
            value of the high frequency flank.
        s : float
            power law exponent.

        Returns
        -------
        y : array
            estimated log dielectric loss based on the supplied parameters.

        """
        f = 10**(x)

        ff = f/fm
        ffb = ff**(b)
        n = (np.pi)/2

        xc = np.cos(b*n)
        xs = np.sin(b*n)
        r = 1 + 2*ffb*xc + ffb**(2)
        r_w = (1/r)**(0.5*g)

        p2 = 1/(ffb)
        phi = np.arctan((xs/(p2+xc)))
        e_w = deps*r_w*(np.sin(g*phi))

        fs = f**(s)
        e_w2 = N*fs

        y = np.log10(e_w+e_w2)
        return y

    def flank_s(self, x, N, s):
        """
        Function to estimate the high frequency flank contribution from the total fit
        
        While fitting, the deconvoluted flank is based on this function.
        
        

        Parameters
        ----------
        x : float
            frequency.
        N : float
            value of the high frequency flank.
        s : float
            power law exponent.


        Returns
        -------
        y : array
            estimated log high frequency flank.

        """
        f = 10**(x)
        fs = f**(s)
        #n = 10**N
        e_w2 = N*fs
        y = np.log10(e_w2)
        return y

    def cond_s(self, x, cond, s):
        """
        Function to estimate the conductivity contribution from the total fit
        
        While fitting, the deconvoluted conductivity is based on this function.

        Parameters
        ----------
        x : float
            frequency.
        cond : float
            DC conductivity contribution.
        s : float
            power law exponent.

        Returns
        -------
        y : array
            estimated log conductivity.

        """
        f = 10**(x)
        fs = f**(s)
        condb = cond/fs
        y = np.log10(condb)
        return y


    def double_hn(self, x, b1, g1, fm1, deps1, b2, g2, fm2,deps2):
        """
        HN fit function to fit two dielectric loss peaks

        Parameters
        ----------
        x : float
            frequency.
        b1 : float
            symmetric fractional parameter of the 1st peak.
        g1 : float
            asymmetric fractional parameter of the 1st peak.
        fm1 : float
            maximum frequency of the 1st peak.
        deps1 : float
           dielectric strength of the 1st peak.
        b2 : float
            symmetric fractional parameter of the 2nd peak.
        g2 : float
            asymmetric fractional parameter of the 2nd peak.
        fm2 : float
            maximum frequency of the 2nd peak.
        deps2 : float
           dielectric strength of the 2nd peak.
        Returns
        -------
        y : array
            estimated log dielectric loss based on the supplied parameters.

        """
        f = 10**(x)

        ff1 = f/fm1
        ffb1 = ff1**(b1)
        n = (np.pi)/2

        xc1 = np.cos(b1*n)
        xs1 = np.sin(b1*n)
        r1 = 1 + 2*ffb1*xc1 + ffb1**(2)
        r_w1 = (1/r1)**(0.5*g1)

        p2 = 1/(ffb1)
        phi1 = np.arctan((xs1/(p2+xc1)))
        e_w1 = deps1*r_w1*(np.sin(g1*phi1))

        ff2 = f/fm2
        ffb2 = ff2**(b2)
        n = (np.pi)/2

        xc2 = np.cos(b2*n)
        xs2 = np.sin(b2*n)
        r2 = 1 + 2*ffb2*xc2 + ffb2**(2)
        r_w2 = (1/r2)**(0.5*g2)

        p3 = 1/(ffb2)
        phi2 = np.arctan((xs2/(p3+xc2)))
        e_w2 = deps2*r_w2*(np.sin(g2*phi2))

        y = np.log10(e_w1+e_w2)
        return y


    def double_hn_cond(self, x, b1, g1, fm1, deps1, b2, g2, fm2,deps2,cond,s):
        """
        HN fit function to fit two dielectric loss peaks together with conductivity

        Parameters
        ----------
        x : float
            frequency.
        b1 : float
            symmetric fractional parameter of the 1st peak.
        g1 : float
            asymmetric fractional parameter of the 1st peak.
        fm1 : float
            maximum frequency of the 1st peak.
        deps1 : float
           dielectric strength of the 1st peak.
        b2 : float
            symmetric fractional parameter of the 2nd peak.
        g2 : float
            asymmetric fractional parameter of the 2nd peak.
        fm2 : float
            maximum frequency of the 2nd peak.
        deps2 : float
           dielectric strength of the 2nd peak.
        cond : float
            DC conductivity contribution.
        s : float
            power law exponent.

        Returns
        -------
        y : array
            estimated log dielectric loss based on the supplied parameters.

        """
        f = 10**(x)

        ff1 = f/fm1
        ffb1 = ff1**(b1)
        n = (np.pi)/2

        xc1 = np.cos(b1*n)
        xs1 = np.sin(b1*n)
        r1 = 1 + 2*ffb1*xc1 + ffb1**(2)
        r_w1 = (1/r1)**(0.5*g1)

        p2 = 1/(ffb1)
        phi1 = np.arctan((xs1/(p2+xc1)))
        e_w1 = deps1*r_w1*(np.sin(g1*phi1))

        ff2 = f/fm2
        ffb2 = ff2**(b2)
        n = (np.pi)/2

        xc2 = np.cos(b2*n)
        xs2 = np.sin(b2*n)
        r2 = 1 + 2*ffb2*xc2 + ffb2**(2)
        r_w2 = (1/r2)**(0.5*g2)

        p3 = 1/(ffb2)
        phi2 = np.arctan((xs2/(p3+xc2)))
        e_w2 = deps2*r_w2*(np.sin(g2*phi2))

        fs = f**(s)
        condb = cond/fs
        y = np.log10(e_w1+e_w2+condb)
        return y

    def dump_parameters_hn(self):
        """
        dumps the initial fit parameters for hn function as a dictionary in a json file
        to load it during curve fitting

        Returns
        -------
         None

        """
        b = float(input("enter the beta value:"))
        g = float(input("enter the gamma value:"))
        lf = float(input("enter the fm:"))
        d = float(input("enter the deps:"))
        cond = float(input("enter the cond:"))
        s = float(input("enter the s:"))
        f = 10**lf

        par = {"beta": b, "gamma": g,  "freq": f,
               "deps": d, "cond": cond, "n": s}

        with open('HN.json', "w") as outfile:
            json.dump(par, outfile)

        with open('HN.json', "r") as openfile:
            loaded_par = json.load(openfile)

        print("dumped_parameters", loaded_par)
        return ()

    def dump_parameters_double_hn(self):
        """
        dumps the initial fit parameters for double hn function as a dictionary in a json file
        to load it during curve fitting

        Returns
        -------
         None

        """
        b1 = float(input("enter the beta1 value:"))
        g1 = float(input("enter the gamma1 value:"))
        lf1 = float(input("enter the fmax1:"))
        d1 = float(input("enter the deps1:"))
        b2 = float(input("enter the beta2 value:"))
        g2 = float(input("enter the gamma2 value:"))
        lf2 = float(input("enter the fmax2:"))
        d2 = float(input("enter the deps2:"))
        cond = float(input("enter the cond:"))
        s = float(input("enter the s:"))
        f1 = 10**lf1
        f2 = 10**lf2


        par = {"beta1": b1, "gamma1": g1, "freq1": f1, "deps1": d1, "beta2":b2,"gamma2":g2,"freq2":f2,"deps2":d2, "cond":cond, "n":s}

        with open('double_HN.json', "w") as outfile:
            json.dump(par, outfile)

        with open('double_HN.json', "r") as openfile:
            loaded_par = json.load(openfile)

        print("dumped_parameters", loaded_par)
        return ()

    def dump_parameters_flank(self):
        """
        dumps the initial fit parameters for hn function with flank as a dictionary in a json file
        to load it during curve fitting

        Returns
        -------
         None

        """
        b = float(input("enter the beta value:"))
        g = float(input("enter the gamma value:"))
        lf = float(input("enter the fm:"))
        d = float(input("enter the deps:"))
        N = float(input("enter the second relaxation rate:"))

        s = float(input("enter the slope:"))
        f = 10**lf

        par = {"beta": b, "gamma": g, "freq": f, "deps": d, "flank": N, "n": s}
        with open('HN.json', "w") as outfile:
            json.dump(par, outfile)

        with open('HN.json', "r") as openfile:
            loaded_par = json.load(openfile)

        print("dumped_parameters", loaded_par)
        return ()

    def sel_function(self):
        """
        A function to select the type of fit function during curve fitting

        Returns
        -------
        func_decision : int
            choice of the fit function.

        """
        func_decision = int(input(
            "Choose the fit function\n 1 -- HN, 2 -- HN with cond, 3 -- Hn-flank, 4 -- double HN, 5 -- double HN with cond:"))
        return func_decision

    def initial_view_flank(self, x, y):
        """
        plots the hn function with flank based on the initial parameters given 
        via the dump_parameters method

        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log dielectric loss.

        Returns
        -------
        None.

        """

        try:

            open('HN.json')

        except FileNotFoundError as e:
            print(
                f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:

           with open('HN.json', "r") as openfile:
            loaded_par = json.load(openfile)

            print("loaded parameters \n", loaded_par)

            
            init_fit_par = loaded_par   
            hn_sub_par = loaded_par['beta'], loaded_par['gamma'], loaded_par[
                'freq'], loaded_par['deps'], loaded_par['cond'], loaded_par['n']
            hn_sub = self.hn_flank(x, *hn_sub_par)
            plt.figure()
            plt.scatter(x, y,label='data')
            plt.plot(x, hn_sub, 'b',label='initial guess')
            plt.legend()
            
        return init_fit_par   

    def initial_view_hn(self, x, y):
        """
        plots the hn function based on the initial parameters given 
        via the dump_parameters method

        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log dielectric loss.

        Returns
        -------
        None.

        """

        try:

            open('HN.json')

        except FileNotFoundError as e:
            print(
                f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:

           with open('HN.json', "r") as openfile:
            loaded_par = json.load(openfile)

            print("loaded parameters \n", loaded_par)
            
            hn_par_index = ['beta','gamma','freq','deps']
            init_fit_par = {key:value for key,value in loaded_par.items() if key in hn_par_index}

            hn_sub_par = loaded_par['beta'], loaded_par['gamma'], loaded_par['freq'], loaded_par['deps']

            hn_sub = self.hn(x, *hn_sub_par)

            plt.figure()
            plt.scatter(x, y,label='data')
            plt.plot(x, hn_sub, 'b',label='initial guess')
            plt.legend()
        return init_fit_par    

    def initial_view_hn_cond(self, x, y):
        """
        plots the hn function with conductivity based on the initial parameters given 
        via the dump_parameters method


        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log dielectric loss.

        Returns
        -------
        None.

        """

        try:

            open('HN.json')

        except FileNotFoundError as e:
            print(
                f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:

           with open('HN.json', "r") as openfile:
            loaded_par = json.load(openfile)
            init_fit_par = loaded_par 
            
            print("loaded parameters \n", loaded_par)

            hn_sub_par = loaded_par['beta'], loaded_par['gamma'], loaded_par['freq'], loaded_par['deps']
            cond_sub_par = loaded_par['cond'], loaded_par['n']
            hn_sub = self.hn(x, *hn_sub_par)
            cond_sub = self.cond_s(x, *cond_sub_par)

            plt.figure()
            plt.scatter(x, y,label='data')
            plt.plot(x, hn_sub, 'b',label='loss peak')
            plt.plot(x, cond_sub, 'r',label='conductivity')
            plt.legend()
        return init_fit_par
    
    def initial_view_double_hn(self, x, y):
        """
        plots the double hn function based on the initial parameters given 
        via the dump_parameters method


        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log dielectric loss.
        Returns
        -------
        None.

        """

        try:

            open('double_HN.json')

        except FileNotFoundError as e:
            print(
                f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:

           with open('double_HN.json', "r") as openfile:
            loaded_par = json.load(openfile)

            print("loaded parameters \n", loaded_par)
            
            double_hn_par_index = ['beta1','gamma1','freq1','deps1','beta2','gamma2','freq2','deps2']
            init_fit_par = {key:value for key,value in loaded_par.items() if key in double_hn_par_index}

            hn_sub_par1 = loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1']
            hn_sub_par2 = loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2']

            hn_sub1 = self.hn(x, *hn_sub_par1)
            hn_sub2 = self.hn(x, *hn_sub_par2)

            plt.figure()
            plt.scatter(x, y,label='data')
            plt.plot(x, hn_sub1, 'b',label='initial guess peak1')
            plt.plot(x, hn_sub2, 'g',label='initial guess peak2')
            plt.legend()
        return init_fit_par
    def initial_view_double_hn_cond(self, x, y):
        """
        plots the double hn function with conductivity based on the initial parameters given 
        via the dump_parameters method


        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log dielectric loss.
        Returns
        -------
        None.

        """

        try:

            open('double_HN.json')

        except FileNotFoundError as e:
            print(
                f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:

           with open('double_HN.json', "r") as openfile:
            loaded_par = json.load(openfile)
            
            init_fit_par = loaded_par 
            print("loaded parameters \n", loaded_par)

            hn_sub_par1 = loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1']
            hn_sub_par2 = loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2']
            cond_sub_par = loaded_par['cond'], loaded_par['n']
            hn_sub1 = self.hn(x, *hn_sub_par1)
            hn_sub2 = self.hn(x, *hn_sub_par2)
            cond_sub = self.cond_s(x, *cond_sub_par)

            plt.figure()
            plt.scatter(x, y,label='data')
            plt.plot(x, hn_sub1, 'b',label='initial guess peak1')
            plt.plot(x, hn_sub2, 'g',label='initial guess peak2')
            plt.plot(x, cond_sub, 'r',label='conductivity')
            plt.legend()
        return init_fit_par     

    def fit(self, x, y):
        """    
         Fits the dielectric loss data with choice of fit function
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
             log dielectric loss.

        Returns
        -------
        fit_par : dictionary
            dictionary containing the fit parameters.

        """


        func_number = self.sel_function()
        x1 = np.array(x)
        y1 = np.array(y)

        global popt1

        plt.figure()

        global popt2
        global b, g, fm, deps, cond, s, l_f, b1, g1,fm1,deps1,b2,g2,fm2,deps2,l_f1,l_f2,flank,fit_par

        if func_number == 1:

            try:

              open('HN.json')

            except FileNotFoundError as e:
                print(
                    f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
            else:
               with open('HN.json', "r") as openfile:
                   loaded_par = json.load(openfile)

               hn_p0 = [loaded_par['beta'], loaded_par['gamma'],
                        loaded_par['freq'], loaded_par['deps']]
               popt1, pcov2 = curve_fit(self.hn, x1, y1, hn_p0, bounds =((0, 0, 1e-7, 0), (1, 1, 1e7, np.inf)),absolute_sigma=True)
               yfit2 = self.hn(x1, *popt1)

               plt.scatter(x1, y1, marker='o', color='r',
                           facecolors='none', label='data', s=100, linewidth=2)
               plt.plot(x1, yfit2, 'm--', label='weighted fit', linewidth=2)
               plt.xlabel('log ( f [Hz])')
               plt.ylabel('log $\epsilon$"')
               plt.legend()

               b, g, fm, deps = popt1[0:4]

               cond = 0
               s = 0
               n_s = np.sin(np.pi*b/(2+2*g))**(1/b)
               n_s2 = 1/(np.sin(np.pi*b*g/(2+2*g))**(1/b))
               fmax = fm*n_s*n_s2
               l_f = np.log10(fmax)
               print(*popt1)
               print('log fmax:', l_f)
               fit_par = {"beta": b, "gamma": g,  "freq": fm,
                          "deps": deps, "cond": cond, "n": s}
               with open('HN.json', "w") as outfile:
                   json.dump(fit_par, outfile)

               with open('HN.json', "r") as openfile:
                   loaded_par = json.load(openfile)

               print("fit parameters dumped for next iteration", loaded_par)

        elif func_number == 2:

            try:

              open('HN.json')

            except FileNotFoundError as e:
                print(
                    f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
            else:
               with open('HN.json', "r") as openfile:
                   loaded_par = json.load(openfile)

               p0 = [loaded_par['beta'], loaded_par['gamma'], loaded_par['freq'],
                     loaded_par['deps'], loaded_par['cond'], loaded_par['n']]
               popt2, pcov2 = curve_fit(self.hn_cond, x1, y1, p0, bounds =((0, 0, 1e-7, 0, 0, 0), (1, 1, 1e7,np.inf,np.inf,1)),absolute_sigma=True)
               yfit3 = self.hn_cond(x1, *popt2)

               plt.scatter(x1, y1, marker='o', color='r',
                           facecolors='none', label='data', s=100, linewidth=2)
               plt.plot(x1, yfit3, 'm--', label='weighted fit', linewidth=2)
               plt.xlabel('log ( f [Hz])')
               plt.ylabel('log $\epsilon$"')
               

               hn_sub_par = popt2[0], popt2[1], popt2[2], popt2[3]
               cond_sub_par = popt2[4], popt2[5]
               hn_sub = self.hn(x1, *hn_sub_par)
               cond_sub = self.cond_s(x1, *cond_sub_par)
               plt.plot(x1, hn_sub, 'b',label='loss peak')
               plt.plot(x1, cond_sub, 'g',label='conductivity')
               plt.legend()
               

               print(*popt2)
               b, g, fm, deps, cond, s = popt2[:]

               n_s = np.sin(np.pi*b/(2+2*g))**(1/b)
               n_s2 = 1/(np.sin(np.pi*b*g/(2+2*g))**(1/b))
               fmax = fm*n_s*n_s2
               l_f = np.log10(fmax)
               print('log fmax:', l_f)
               fit_par = {"beta": b, "gamma": g,  "freq": fm,
                          "deps": deps, "cond": cond, "n": s}
               with open('HN.json', "w") as outfile:
                   json.dump(fit_par, outfile)

               with open('HN.json', "r") as openfile:
                   loaded_par = json.load(openfile)

               print("fit parameters dumped for next iteration", loaded_par)

        elif func_number == 3:
            try:

              open('HN_flank.json')

            except FileNotFoundError as e:
                print(
                    f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
            else:
               with open('HN_flank.json', "r") as openfile:
                   loaded_par = json.load(openfile)

               p0 = [loaded_par['beta'], loaded_par['gamma'], loaded_par['freq'],
                     loaded_par['deps'], loaded_par['flank'], loaded_par['n']]
               popt2, pcov2 = curve_fit(self.hn_flank, x1, y1, p0, bounds =((0, 0, 0, 0, 0, 0), (1, 1, 1e7,np.inf,np.inf,1)), absolute_sigma=True)
               yfit3 = self.hn_flank(x1, *popt2)

               plt.scatter(x1, y1, marker='o', color='r',
                           facecolors='none', label='data', s=100, linewidth=2)
               plt.plot(x1, yfit3, 'm--', label='weighted fit', linewidth=2)
               plt.xlabel('log ( f [Hz])')
               plt.ylabel('log $\epsilon$"')
               

               hn_sub_par = popt2[0], popt2[1], popt2[2], popt2[3]
               n_sub_par = popt2[4], popt2[5]
               hn_sub = self.hn(x1, *hn_sub_par)
               n_sub = self.flank_s(x1, *n_sub_par)
               plt.plot(x1, hn_sub, 'b',label='loss peak')
               plt.plot(x1, n_sub, 'g',labe='flank contribution')
               plt.legend()
               
               print(*popt2)

               b, g, fm, deps, cond, s = popt2[:]
               b = popt2[0]
               g = popt2[1]
               fm = popt2[2]
               deps = popt2[3]
               flank = popt2[4]
               s = popt2[5]
               n_s = np.sin(np.pi*b/(2+2*g))**(1/b)
               n_s2 = 1/(np.sin(np.pi*b*g/(2+2*g))**(1/b))
               fmax = fm*n_s*n_s2
               l_f = np.log10(fmax)

               print('log fmax:', l_f)
               fit_par = {"beta": b, "gamma": g,  "freq": fm,
                          "deps": deps, "flank": flank, "n": s}
               with open('HN.json', "w") as outfile:
                   json.dump(fit_par, outfile)

               with open('HN.json', "r") as openfile:
                   loaded_par = json.load(openfile)

               print("fit parameters dumped for next iteration", loaded_par)

        elif func_number == 4:
            try:

              open('double_HN.json')

            except FileNotFoundError as e:
                print(
                    f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
            else:

                with open('double_HN.json', "r") as openfile:
                    loaded_par = json.load(openfile)

                p0 = [loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1'],
                    loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2']]
                print(p0)
                popt2, pcov2 = curve_fit(self.double_hn, x1, y1, p0, bounds =((0, 0, 1e-7, 0, 0, 0, 1e-7, 0), (1,1,1e7,np.inf,1,1,1e7,np.inf)),absolute_sigma=True)
                yfit4 = self.double_hn(x1, *popt2)
                plt.scatter(x1, y1, marker='o', color='r',
                            facecolors='none', label='data', s=100, linewidth=2)
                plt.plot(x1, yfit4, 'm--', label='weighted fit', linewidth=2)
                plt.xlabel('log ( f [Hz])')
                plt.ylabel('log $\epsilon$"')
                

                hn_sub_par1 = popt2[0], popt2[1], popt2[2], popt2[3]
                hn_sub_par2 = popt2[4], popt2[5], popt2[6], popt2[7]
                hn_sub1 = self.hn(x1, *hn_sub_par1)
                hn_sub2 = self.hn(x1, *hn_sub_par2)
                plt.plot(x1, hn_sub1, 'b',label='peak1')
                plt.plot(x1, hn_sub2, 'g', label='peak2')
                plt.legend()

                b1, g1, fm1, deps1, b2, g2, fm2, deps2 = popt2[:]
                cond, s = 0, 1
                nr_1 = np.sin(np.pi*b1/(2+2*g1))**(1/b1)
                nr_2 = 1/(np.sin(np.pi*b1*g1/(2+2*g1))**(1/b1))
                fmax1 = fm1*nr_1*nr_2
                nr_3 = np.sin(np.pi*b2/(2+2*g2))**(1/b2)
                nr_4 = 1/(np.sin(np.pi*b2*g2/(2+2*g2))**(1/b2))
                fmax2 = fm2*nr_3*nr_4
                l_f1 = np.log10(fmax1)
                l_f2 = np.log10(fmax2)

                print('log fmax1:', l_f1, "\n" 'log fmax2:', l_f2)
                fit_par = {"beta1": b1, "gamma1": g1, "freq1": fm1, "deps1": deps1, "beta2":b2,"gamma2":g2,"freq2":fm2,"deps2":deps2, "cond":cond, "n":s}

                with open('double_HN.json', "w") as outfile:
                    json.dump(fit_par, outfile)
                with open('double_HN.json', "r") as openfile:
                    loaded_par = json.load(openfile)
                print("fit parameters dumped for next iteration", loaded_par)

        elif func_number == 5:
            try:

              open('double_HN.json')

            except FileNotFoundError as e:
                print(
                    f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
            else:

                with open('double_HN.json', "r") as openfile:
                    loaded_par = json.load(openfile)

                p0 = [loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1'], loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2'], loaded_par['cond'],loaded_par['n']]
                print(p0)
                popt2, pcov2 = curve_fit(self.double_hn_cond, x1, y1, p0, bounds =((0, 0, 1e-7, 0, 0, 0, 1e-7, 0, 0,0),(1,1,1e7,np.inf,1,1,1e7,np.inf,np.inf,1)),absolute_sigma=True)
                yfit4 = self.double_hn_cond(x1, *popt2)
                plt.scatter(x1, y1, marker='o', color='r',
                            facecolors='none', label='data', s=100, linewidth=2)
                plt.plot(x1, yfit4, 'm--', label='weighted fit', linewidth=2)
                plt.xlabel('log ( f [Hz])')
                plt.ylabel('log $\epsilon$"')
                

                hn_sub_par1 = popt2[0], popt2[1], popt2[2], popt2[3]
                hn_sub_par2 = popt2[4], popt2[5], popt2[6], popt2[7]
                cond_sub_par = popt2[8], popt2[9]
                hn_sub1 = self.hn(x1, *hn_sub_par1)
                hn_sub2 = self.hn(x1, *hn_sub_par2)
                cond_sub = self.cond_s(x1, *cond_sub_par)
                plt.plot(x, hn_sub1, 'b',label='initial guess peak1')
                plt.plot(x, hn_sub2, 'g',label='initial guess peak2')
                plt.plot(x, cond_sub, 'r',label='conductivity')
                plt.legend()

                b1, g1, fm1, deps1, b2, g2, fm2, deps2, cond,s = popt2[:]

                nr_1 = np.sin(np.pi*b1/(2+2*g1))**(1/b1)
                nr_2 = 1/(np.sin(np.pi*b1*g1/(2+2*g1))**(1/b1))
                fmax1 = fm1*nr_1*nr_2
                nr_3 = np.sin(np.pi*b2/(2+2*g2))**(1/b2)
                nr_4 = 1/(np.sin(np.pi*b2*g2/(2+2*g2))**(1/b2))
                fmax2 = fm2*nr_3*nr_4
                l_f1 = np.log10(fmax1)
                l_f2 = np.log10(fmax2)

                print('log fmax1:', l_f1, "\n" 'log fmax2:', l_f2)
                fit_par = {"beta1": b1, "gamma1": g1, "freq1": fm1, "deps1": deps1, "beta2":b2,"gamma2":g2,"freq2":fm2,"deps2":deps2,"cond":cond, "n":s}

                with open('double_HN.json', "w") as outfile:
                    json.dump(fit_par, outfile)
                with open('double_HN.json', "r") as openfile:
                    loaded_par = json.load(openfile)
                print("fit parameters dumped for next iteration", loaded_par)
                
                
        return fit_par     

    def save_fit_hn(self, T):
        """
        saves the fit parameters of hn function in a file
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
        res_file = open(ana_file, 'a')
        res_file.write(f'{T}' + '\t' + f'{b:.3f}' + '\t' + f'{g:.3f}' + '\t' + f'{deps:.3f}' +
                       '\t' + f'{l_f:.3f}' + '\t' + f'{cond:.3f}'+'\t' + f'{s:.3f}' + "\n")
        return()

    def save_fit_hn_cond(self, T):
        """
        saves the fit parameters of hn function with conductivity
        in  a file, the file must be created via create_analysis_file function


        Parameters
        ----------
        T : float
            Temperature,or can also be an integer
            that corresponds to a file number during analysis.
        Returns
        -------
        None.

        """
        res_file = open(ana_file, 'a')
        res_file.write(f'{T}' + '\t' + f'{b:.3f}' + '\t' + f'{g:.3f}' + '\t' + f'{deps:.3f}' +
                       '\t' + f'{l_f:.3f}' + '\t' + f'{cond:.3f}'+'\t' + f'{s:.3f}' + "\n")
        return()

    def save_fit_hn_flank(self, T):
        """
        saves the fit parameters of hn function with flank
        in  a file, the file must be created via create_analysis_file function


        Parameters
        ----------
        T : float
            Temperature,or can also be an integer
            that corresponds to a file number during analysis.
        Returns
        -------
        None.

        """
        res_file = open(ana_file, 'a')
        res_file.write(f'{T}' + '\t' + f'{b:.3f}' + '\t' + f'{g:.3f}' + '\t' + f'{deps:.3f}' +
                       '\t' + f'{l_f:.3f}' + '\t' + f'{flank:.3f}'+'\t' + f'{s:.3f}' + "\n")
        return()

    def save_fit_double_HN(self, T):
        """
        saves the fit parameters of double hn function in  a file, 
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
        res_file = open(ana_file, 'a')
        res_file.write(f'{T}' + '\t' + f'{b1:.3f}' + '\t' + f'{g1:.3f}' + '\t' + f'{deps1:.3f}' + '\t' + f'{l_f1:.3f}' + '\t' +
                       f'{b2:.3f}' + '\t' + f'{g2:.3f}' + '\t' + f'{deps2:.3f}' + '\t' + f'{l_f2:.3f}' + '\t' + f'{cond:.3f}'+'\t' + f'{s:.3f}'+"\n")
        return()

    plt.show()


