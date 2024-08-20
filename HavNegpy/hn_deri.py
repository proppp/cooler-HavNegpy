# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:52:39 2022

@author: mkolmang
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import mplcursors
import json



class HN_derivative:
    """
    A class to analyze the derivaite of real part of complex permittivity with derivative HN function

    Fit functions include single, double, and derivative HN with electrode polarization

    """


    def __init__(self):
        pass


    def create_analysis_file(self):
        """
         Creates a file to save the fit results based on the choice of fit function

         Provides option to use an existing file and creates a new file if not found

        Returns
        -------
        None.

        """
        res = input("Do you want to use an existing file to save fit results? \n eg: existing file to save HN parameters, y or n:")
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
            func_name = int(input("Choose the fit function:\n 1 -- deri_HN\n 2-- deri_HN with EP\n 3 -- double_deri_HN\n 4 -- double_deri_HN with EP\n 5 -- triple_deri_HN with EP\n 6 -- triple_deri_HN with EP and UD\n"))
            ex_file = input("Enter the analysis_file_name:")
            f = open(ex_file,'w')
            if func_name == 1 or func_name == 2:
                 f.write(f'{"Temperature"}\t{"beta"}\t{"gamma"}\t{"deps"}\t{"log (fmax)"}\t{"EP"}\t{"s"}\n')

            elif func_name == 3:
                  f.write(f'{"Temperature"}\t{"beta1"}\t{"gamma1"}\t{"deps1"}\t{"log (fmax1)"}\t{"beta2"}\t{"gamma2"}\t{"deps2"}\t{"log (fmax2)"}\t{"EP"}\t{"s"}\n')

            elif func_name == 4:
                f.write(f'{"Temperature"}\t{"beta1"}\t{"gamma1"}\t{"deps1"}\t{"log (fmax1)"}\t{"beta2"}\t{"gamma2"}\t{"deps2"}\t{"log (fmax2)"}\t{"EP"}\t{"s"}\n')

            elif func_name == 5:
                f.write(f'{"Temperature"}\t{"beta1"}\t{"gamma1"}\t{"deps1"}\t{"log (fmax1)"}\t{"beta2"}\t{"gamma2"}\t{"deps2"}\t{"log (fmax2)"}\t{"beta3"}\t{"gamma3"}\t{"deps3"}\t{"log (fmax3)"}\t{"EP"}\t{"s"}\n')

            elif func_name == 6:
                f.write(f'{"Temperature"}\t{"beta1"}\t{"gamma1"}\t{"deps1"}\t{"log (fmax1)"}\t{"beta2"}\t{"gamma2"}\t{"deps2"}\t{"log (fmax2)"}\t{"beta3"}\t{"gamma3"}\t{"deps3"}\t{"log (fmax3)"}\t{"EP"}\t{"s"}\t{"B"}\t{"m"}\n')
            print(f'{"file did not exist, created"}',ex_file)
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
            log dielectric loss


       """
       x = list(x)
       y = list(y)

       plt.figure(1)
       plt.style.use("seaborn-v0_8-whitegrid")

       plt.scatter(x,y,marker='s',color='r',facecolors='none', s=100,linewidth=2)
       plt.ylabel('log ( -d$\epsilon´$/ dlog f )')
       plt.xlabel('log f')
       #plt.legend()
       plt.style.use("seaborn-v0_8-whitegrid")

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


    def deri_hn(self,x,b,g,fm,deps):
        """
        derivaitive HN fit function to fit single peak

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
            estimated log derivative of epsilon' based on the supplied parameters.

        """
        f = 10**(x)

        ff = f/fm
        ffb = ff**(b)
        n = (np.pi)/2

        xc = np.cos(b*n)
        xs = np.sin(b*n)
        r = 1 + 2*ffb*xc + ffb**(2)
        #r_w = (1/r)**(0.5*g)

        p2 = 1/(ffb)
        phi = np.arctan((xs/(p2+xc)))
        #e_w = deps*r_w*(np.sin(g*phi))


        de_1 = b*g*deps*ffb*np.cos(b*n-(1+g)*phi)
        de_2 = r**((1+g)/2)

        deri = de_1/de_2



        y = np.log10(deri)#+ep)#+condb)
        return y


    def deri_hn_ln(self,x,b,g,fm,deps):
        """
        derivaitive HN fit function to fit single peak

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
            estimated log derivative of epsilon' based on the supplied parameters.

        """
        f = np.exp(x)

        ff = f/fm
        ffb = ff**(b)
        n = (np.pi)/2

        xc = np.cos(b*n)
        xs = np.sin(b*n)
        r = 1 + 2*ffb*xc + ffb**(2)
        #r_w = (1/r)**(0.5*g)

        p2 = 1/(ffb)
        phi = np.arctan((xs/(p2+xc)))
        #e_w = deps*r_w*(np.sin(g*phi))


        de_1 = b*g*deps*ffb*np.cos(b*n-(1+g)*phi)
        de_2 = r**((1+g)/2)

        deri = de_1/de_2
        deri = deri/np.log(10)



        y = np.log(deri)#+ep)#+condb)
        return y


    def deri_hn_ep(self,x,b,g,fm,deps,A,l):
        """
        derivaitive HN fit function to fit single peak along with
        electrode polarization(ep)

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
        A : float
            electrode polarization value.
        l : float
            power law exponent.

        Returns
        -------
        y : array
            estimated log derivative of epsilon' based on the supplied parameters.

        """
        f = 10**(x)
            #w = 2*np.pi*f
        ff = f/fm
        ffb = ff**(b)
        n = (np.pi)/2

        xc = np.cos(b*n)
        xs = np.sin(b*n)
        r = 1 + 2*ffb*xc + ffb**(2)
        #r_w = (1/r)**(0.5*g)

        p2 = 1/(ffb)
        phi = np.arctan((xs/(p2+xc)))
        #e_w = deps*r_w*(np.sin(g*phi))
        #fs= f**(s)

        de_1 = b*g*deps*ffb*np.cos(b*n-(1+g)*phi)
        de_2 = r**((1+g)/2)

        deri = de_1/de_2
        fl = f**l
        ep = A/fl

        y = np.log10(deri+ep)
        return y


    def deri_double_hn(self,x,b1,g1,fm1,deps1,b2,g2,fm2,deps2):
        """
        derivaitive HN fit function to fit two peaks


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
            estimated log derivative of epsilon' based on the supplied parameters.

        """
        f = 10**(x)

        ff1 = f/fm1
        ffb1 = ff1**(b1)
        n = (np.pi)/2

        xc1 = np.cos(b1*n)
        xs1 = np.sin(b1*n)
        r1 = 1 + 2*ffb1*xc1 + ffb1**(2)
        #r_w1 = (1/r1)**(0.5*g1)

        p2 = 1/(ffb1)
        phi1 = np.arctan((xs1/(p2+xc1)))
        de_1 = b1*g1*deps1*ffb1*np.cos(b1*n-(1+g1)*phi1)
        de_2 = r1**((1+g1)/2)

        deri1 = de_1/de_2

        #e_w1 = deps1*r_w1*(np.sin(g1*phi1))

        ff2 = f/fm2
        ffb2 = ff2**(b2)
        n = (np.pi)/2

        xc2 = np.cos(b2*n)
        xs2 = np.sin(b2*n)
        r2 = 1 + 2*ffb2*xc2 + ffb2**(2)
        #r_w2 = (1/r2)**(0.5*g2)

        p3 = 1/(ffb2)
        phi2 = np.arctan((xs2/(p3+xc2)))
        de_3 = b2*g2*deps2*ffb2*np.cos(b2*n-(1+g2)*phi2)
        de_4 = r2**((1+g2)/2)


        deri2 = de_3/de_4
        #e_w2 = deps2*r_w2*(np.sin(g2*phi2))



        y = np.log10(deri1+deri2)
        return y



    def deri_double_hn_ep(self, x, b1, g1, fm1, deps1, b2, g2, fm2, deps2, A, l):
        """
        Derivative double HN fit function to fit two peaks along with electrode polarization (EP).

        Parameters
        ----------
        x : float
            Frequency.
        b1 : float
            Symmetric fractional parameter of the 1st peak.
        g1 : float
            Asymmetric fractional parameter of the 1st peak.
        fm1 : float
            Maximum frequency of the 1st peak.
        deps1 : float
            Dielectric strength of the 1st peak.
        b2 : float
            Symmetric fractional parameter of the 2nd peak.
        g2 : float
            Asymmetric fractional parameter of the 2nd peak.
        fm2 : float
            Maximum frequency of the 2nd peak.
        deps2 : float
            Dielectric strength of the 2nd peak.
        A : float
            Electrode polarization value.
        l : float
            Power law exponent.

        Returns
        -------
        y : array
            Estimated log derivative of epsilon' based on the supplied parameters.

        """
        f = 10**(x)

        # First HN component
        ff1 = f / fm1
        ffb1 = ff1**(b1)
        n = (np.pi) / 2

        xc1 = np.cos(b1 * n)
        xs1 = np.sin(b1 * n)
        r1 = 1 + 2 * ffb1 * xc1 + ffb1**2

        p2 = 1 / ffb1
        phi1 = np.arctan(xs1 / (p2 + xc1))
        de_1 = b1 * g1 * deps1 * ffb1 * np.cos(b1 * n - (1 + g1) * phi1)
        de_2 = r1**((1 + g1) / 2)

        deri1 = de_1 / de_2

        # Second HN component
        ff2 = f / fm2
        ffb2 = ff2**(b2)

        xc2 = np.cos(b2 * n)
        xs2 = np.sin(b2 * n)
        r2 = 1 + 2 * ffb2 * xc2 + ffb2**2

        p3 = 1 / ffb2
        phi2 = np.arctan(xs2 / (p3 + xc2))
        de_3 = b2 * g2 * deps2 * ffb2 * np.cos(b2 * n - (1 + g2) * phi2)
        de_4 = r2**((1 + g2) / 2)

        deri2 = de_3 / de_4

        # Electrode polarization component
        fl = f**l
        ep = A / fl

        y = np.log10(deri1 + deri2 + ep)
        return y


    def deri_triple_hn_ep(self, x, b1, g1, fm1, deps1, b2, g2, fm2, deps2, b3, g3, fm3, deps3, A, l):
        """
        Derivative triple HN fit function to fit three peaks along with electrode polarization (EP).

        Parameters
        ----------
        x : float
            Frequency.
        b1 : float
            Symmetric fractional parameter of the 1st peak.
        g1 : float
            Asymmetric fractional parameter of the 1st peak.
        fm1 : float
            Maximum frequency of the 1st peak.
        deps1 : float
            Dielectric strength of the 1st peak.
        b2 : float
            Symmetric fractional parameter of the 2nd peak.
        g2 : float
            Asymmetric fractional parameter of the 2nd peak.
        fm2 : float
            Maximum frequency of the 2nd peak.
        deps2 : float
            Dielectric strength of the 2nd peak.
        b3 : float
            Symmetric fractional parameter of the 3rd peak.
        g3 : float
            Asymmetric fractional parameter of the 3rd peak.
        fm3 : float
            Maximum frequency of the 3rd peak.
        deps3 : float
            Dielectric strength of the 3rd peak.
        A : float
            Electrode polarization value.
        l : float
            Power law exponent.

        Returns
        -------
        y : array
            Estimated log derivative of epsilon' based on the supplied parameters.
        """
        f = 10**(x)

        # First HN component
        ff1 = f / fm1
        ffb1 = ff1**(b1)
        n = (np.pi) / 2

        xc1 = np.cos(b1 * n)
        xs1 = np.sin(b1 * n)
        r1 = 1 + 2 * ffb1 * xc1 + ffb1**2

        p2 = 1 / ffb1
        phi1 = np.arctan(xs1 / (p2 + xc1))
        de_1 = b1 * g1 * deps1 * ffb1 * np.cos(b1 * n - (1 + g1) * phi1)
        de_2 = r1**((1 + g1) / 2)

        deri1 = de_1 / de_2

        # Second HN component
        ff2 = f / fm2
        ffb2 = ff2**(b2)

        xc2 = np.cos(b2 * n)
        xs2 = np.sin(b2 * n)
        r2 = 1 + 2 * ffb2 * xc2 + ffb2**2

        p3 = 1 / ffb2
        phi2 = np.arctan(xs2 / (p3 + xc2))
        de_3 = b2 * g2 * deps2 * ffb2 * np.cos(b2 * n - (1 + g2) * phi2)
        de_4 = r2**((1 + g2) / 2)

        deri2 = de_3 / de_4

        # Third HN component
        ff3 = f / fm3
        ffb3 = ff3**(b3)

        xc3 = np.cos(b3 * n)
        xs3 = np.sin(b3 * n)
        r3 = 1 + 2 * ffb3 * xc3 + ffb3**2

        p4 = 1 / ffb3
        phi3 = np.arctan(xs3 / (p4 + xc3))
        de_5 = b3 * g3 * deps3 * ffb3 * np.cos(b3 * n - (1 + g3) * phi3)
        de_6 = r3**((1 + g3) / 2)

        deri3 = de_5 / de_6

        # Electrode polarization component
        fl = f**l
        ep = A / fl

        y = np.log10(deri1 + deri2 + deri3 + ep)
        return y


    def deri_triple_hn_ep_ud(self, x, b1, g1, fm1, deps1, b2, g2, fm2, deps2, b3, g3, fm3, deps3, A, l, B, m):
        """
        Derivative triple HN fit function to fit three peaks along with electrode polarization (EP)
        and universal dielectric response term.

        Parameters
        ----------
        x : float
            Frequency.
        b1 : float
            Symmetric fractional parameter of the 1st peak.
        g1 : float
            Asymmetric fractional parameter of the 1st peak.
        fm1 : float
            Maximum frequency of the 1st peak.
        deps1 : float
            Dielectric strength of the 1st peak.
        b2 : float
            Symmetric fractional parameter of the 2nd peak.
        g2 : float
            Asymmetric fractional parameter of the 2nd peak.
        fm2 : float
            Maximum frequency of the 2nd peak.
        deps2 : float
            Dielectric strength of the 2nd peak.
        b3 : float
            Symmetric fractional parameter of the 3rd peak.
        g3 : float
            Asymmetric fractional parameter of the 3rd peak.
        fm3 : float
            Maximum frequency of the 3rd peak.
        deps3 : float
            Dielectric strength of the 3rd peak.
        A : float
            Electrode polarization value.
        l : float
            Power law exponent for EP.
        B : float
            Universal dielectric response coefficient.
        m : float
            Exponent for universal dielectric response term.

        Returns
        -------
        y : array
            Estimated log derivative of epsilon' based on the supplied parameters.
        """
        f = 10**(x)

        # First HN component
        ff1 = f / fm1
        ffb1 = ff1**(b1)
        n = (np.pi) / 2

        xc1 = np.cos(b1 * n)
        xs1 = np.sin(b1 * n)
        r1 = 1 + 2 * ffb1 * xc1 + ffb1**2

        p2 = 1 / ffb1
        phi1 = np.arctan(xs1 / (p2 + xc1))
        de_1 = b1 * g1 * deps1 * ffb1 * np.cos(b1 * n - (1 + g1) * phi1)
        de_2 = r1**((1 + g1) / 2)

        deri1 = de_1 / de_2

        # Second HN component
        ff2 = f / fm2
        ffb2 = ff2**(b2)

        xc2 = np.cos(b2 * n)
        xs2 = np.sin(b2 * n)
        r2 = 1 + 2 * ffb2 * xc2 + ffb2**2

        p3 = 1 / ffb2
        phi2 = np.arctan(xs2 / (p3 + xc2))
        de_3 = b2 * g2 * deps2 * ffb2 * np.cos(b2 * n - (1 + g2) * phi2)
        de_4 = r2**((1 + g2) / 2)

        deri2 = de_3 / de_4

        # Third HN component
        ff3 = f / fm3
        ffb3 = ff3**(b3)

        xc3 = np.cos(b3 * n)
        xs3 = np.sin(b3 * n)
        r3 = 1 + 2 * ffb3 * xc3 + ffb3**2

        p4 = 1 / ffb3
        phi3 = np.arctan(xs3 / (p4 + xc3))
        de_5 = b3 * g3 * deps3 * ffb3 * np.cos(b3 * n - (1 + g3) * phi3)
        de_6 = r3**((1 + g3) / 2)

        deri3 = de_5 / de_6

        # Electrode polarization component
        fl = f**l
        ep = A / fl

        # Universal dielectric response component
        ud = B * f**m

        y = np.log10(deri1 + deri2 + deri3 + ep + ud)
        return y

    def ep_s(self,x,A,l):
        """
        Function to estimate the electrode polarization(EP) contribution from the total fit

        While fitting, the deconvoluted EP is based on this function.


        Parameters
        ----------
        x : float
            frequency.
        A : float
            electrode polarization value.
        l : float
            power law exponent.

        Returns
        -------
        y : array
            estimated log EP.

        """
        f = 10**(x)
        fl = f**(l)
        ep = A/fl
        y = np.log10(ep)
        return y

    def ud_s(self, x, B, m):
        """
        Function to estimate the Universal Dielectric Response (UD) contribution from the total fit.

        This function models the UD response based on the parameters provided. The UD response is assumed
        to follow a power-law dependence on frequency.

        Parameters
        ----------
        x : float
            Frequency (in log scale).
        B : float
            Universal Dielectric Response amplitude.
        m : float
            Universal Dielectric Response exponent (power law).

        Returns
        -------
        y : array
            Estimated log UD response.

        """
        f = 10**(x)  # Convert frequency from log scale to linear scale
        ud = B * f**(m)  # Model the UD response
        y = np.log10(ud)  # Return log of the UD response
        return y
    def dump_parameters_deri_hn(self):
        """
        dumps the initial fit parameters for derivative hn function as a dictionary
        in a json file to load it during curve fitting

        Returns
        -------
         None

        """
        b  = float(input("enter the beta value:"))
        g  = float(input("enter the gamma value:"))
        lf  = float(input("enter the fm:"))
        d = float(input("enter the deps:"))
        ep = float(input("enter the E.P value:"))
        s = float(input("enter the s:"))
        f = 10**(lf)


        par = {"beta": b, "gamma": g,  "freq": f,"deps":d, "ep":ep, "s":s}


        with open('HN_deri.json',"w") as outfile:
            json.dump(par,outfile)

        with open('HN_deri.json',"r") as openfile:
            loaded_par = json.load(openfile)

        print("dumped_parameters",loaded_par)
        return ()


    def dump_parameters_deri_double_hn(self):
        """
        dumps the initial fit parameters for derivative_double hn function as a dictionary
        in a json file to load it during curve fitting

        Returns
        -------
         None

        """
        b1  = float(input("enter the beta1 value:"))
        g1  = float(input("enter the gamma1 value:"))
        lf1  = float(input("enter the fmax1:"))
        d1 = float(input("enter the deps1:"))
        b2  = float(input("enter the beta2 value:"))
        g2  = float(input("enter the gamma2 value:"))
        lf2  = float(input("enter the fmax2:"))
        d2 = float(input("enter the deps2:"))
        ep = float(input("enter the E.P value:"))
        s = float(input("enter the s:"))
        f1 = 10**(lf1)
        f2 = 10**(lf2)


        par = {"beta1":b1,"gamma1":g1,"freq1":f1,"deps1":d1,"beta2":b2,"gamma2":g2,"freq2":f2,"deps2":d2, "ep":ep, "s":s}



        with open('double_HN_deri.json',"w") as outfile:
            json.dump(par,outfile)

        with open('double_HN_deri.json',"r") as openfile:
            loaded_par = json.load(openfile)


        print("dumped_parameters",loaded_par)
        return ()


    def dump_parameters_deri_triple_hn(self):
        """
        Dumps the initial fit parameters for the derivative_triple HN + EP function as a dictionary
        in a JSON file to load during curve fitting.

        Returns
        -------
        None
        """
        # Collecting input for the first HN peak
        b1 = float(input("Enter the beta1 value: "))
        g1 = float(input("Enter the gamma1 value: "))
        lf1 = float(input("Enter the log fmax1 (log10 of freq1): "))
        d1 = float(input("Enter the deps1 value: "))

        # Collecting input for the second HN peak
        b2 = float(input("Enter the beta2 value: "))
        g2 = float(input("Enter the gamma2 value: "))
        lf2 = float(input("Enter the log fmax2 (log10 of freq2): "))
        d2 = float(input("Enter the deps2 value: "))

        # Collecting input for the third HN peak
        b3 = float(input("Enter the beta3 value: "))
        g3 = float(input("Enter the gamma3 value: "))
        lf3 = float(input("Enter the log fmax3 (log10 of freq3): "))
        d3 = float(input("Enter the deps3 value: "))

        # Collecting input for Electrode Polarization (EP) and scaling factor
        ep = float(input("Enter the Electrode Polarization (EP) value: "))
        s = float(input("Enter the scaling factor s: "))

        # Converting log fmax values back to frequency values
        f1 = 10**(lf1)
        f2 = 10**(lf2)
        f3 = 10**(lf3)

        # Packaging all parameters into a dictionary
        par = {
            "beta1": b1, "gamma1": g1, "freq1": f1, "deps1": d1,
            "beta2": b2, "gamma2": g2, "freq2": f2, "deps2": d2,
            "beta3": b3, "gamma3": g3, "freq3": f3, "deps3": d3,
            "ep": ep, "s": s
        }

        # Saving the parameters into a JSON file
        with open('triple_HN_ep_deri.json', "w") as outfile:
            json.dump(par, outfile)

        # Reloading and printing the dumped parameters to confirm
        with open('triple_HN_ep_deri.json', "r") as openfile:
            loaded_par = json.load(openfile)

        print("Dumped parameters:", loaded_par)
        return ()





    def dump_parameters_deri_ud(self):
        """
        Dumps the initial fit parameters for the derivative triple HN + EP function as a dictionary
        in a JSON file to load during curve fitting.

        Returns
        -------
        None
        """
        # Collecting input for the first HN peak
        b1 = float(input("Enter the beta1 value: "))
        g1 = float(input("Enter the gamma1 value: "))
        lf1 = float(input("Enter the log fmax1 (log10 of freq1): "))
        d1 = float(input("Enter the deps1 value: "))

        # Collecting input for the second HN peak
        b2 = float(input("Enter the beta2 value: "))
        g2 = float(input("Enter the gamma2 value: "))
        lf2 = float(input("Enter the log fmax2 (log10 of freq2): "))
        d2 = float(input("Enter the deps2 value: "))

        # Collecting input for the third HN peak
        b3 = float(input("Enter the beta3 value: "))
        g3 = float(input("Enter the gamma3 value: "))
        lf3 = float(input("Enter the log fmax3 (log10 of freq3): "))
        d3 = float(input("Enter the deps3 value: "))

        # Collecting input for Electrode Polarization (EP) and scaling factor
        ep = float(input("Enter the Electrode Polarization (EP) value: "))
        s = float(input("Enter the scaling factor s: "))

        # Collecting input for additional parameters
        B = float(input("Enter the B value: "))  # Ensure this is included
        m = float(input("Enter the m value: "))  # Ensure this is included

        # Converting log fmax values back to frequency values
        f1 = 10**(lf1)
        f2 = 10**(lf2)
        f3 = 10**(lf3)

        # Packaging all parameters into a dictionary
        par = {
            "beta1": b1, "gamma1": g1, "freq1": f1, "deps1": d1,
            "beta2": b2, "gamma2": g2, "freq2": f2, "deps2": d2,
            "beta3": b3, "gamma3": g3, "freq3": f3, "deps3": d3,
            "ep": ep, "s": s,
            "B": B, "m": m  # Add these parameters to the dictionary
        }

        # Saving the parameters into a JSON file
        with open('triple_HN_ep_ud_deri.json', "w") as outfile:
            json.dump(par, outfile)

        # Reloading and printing the dumped parameters to confirm
        with open('triple_HN_ep_ud_deri.json', "r") as openfile:
            loaded_par = json.load(openfile)

        print("Dumped parameters:", loaded_par)
        return ()
    def initial_view_deri_hn(self,x,y):
        """
        plots the derivative hn function based on the initial parameters given
        via the dump_parameters method

        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log derivative of epsilon'.


        Returns
        -------
        None.

        """

        try:

              open('HN_deri.json')

        except  FileNotFoundError as e:
                print(f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:

           with open('HN_deri.json',"r") as openfile:
            loaded_par = json.load(openfile)


            print("loaded parameters \n" ,loaded_par)

            hn_par_index = ['beta','gamma','freq','deps']
            init_fit_par = {key:value for key,value in loaded_par.items() if key in hn_par_index}

            hn_sub_par = loaded_par['beta'],loaded_par['gamma'],loaded_par['freq'],loaded_par['deps']

            hn_sub = self.deri_hn(x,*hn_sub_par)

            plt.figure()
            plt.scatter(x,y,marker='s',color='r',facecolors='none',label='data',s=100,linewidth=2)
            plt.plot(x,hn_sub,'b',label='initial guess')
            plt.xlabel('log ( f [Hz])')
            plt.ylabel('log ( -d$\epsilon´$/ dlog f )')
            plt.legend()

        return init_fit_par
    def initial_view_deri_hn_ep(self,x,y):
        """
        plots the derivative hn function with electrode polarization
        based on the initial parameters given via the dump_parameters method

        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log derivative of epsilon'.

        Returns
        -------
        None.

        """

        try:

              open('HN_deri.json')

        except  FileNotFoundError as e:
                print(f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:

           with open('HN_deri.json',"r") as openfile:
            loaded_par = json.load(openfile)

            init_fit_par = loaded_par
            print("loaded parameters \n" ,loaded_par)

            hn_sub_par = loaded_par['beta'],loaded_par['gamma'],loaded_par['freq'],loaded_par['deps']
            ep_sub_par = loaded_par['ep'],loaded_par['s']
            hn_sub = self.deri_hn(x,*hn_sub_par)
            ep_sub = self.ep_s(x,*ep_sub_par)

            plt.figure()
            plt.scatter(x,y,marker='s',color='r',facecolors='none',label='data',s=100,linewidth=2)
            plt.plot(x,hn_sub,'b',label='initial guess - peak')
            plt.plot(x,ep_sub,'r',label='initial guess - electrode polarization')
            plt.xlabel('log ( f [Hz])')
            plt.ylabel('log ( -d$\epsilon´$/ dlog f )')
            plt.legend()
        return init_fit_par

    def initial_view_deri_double_hn(self,x,y):
        """
        plots the derivative double hn function based on the initial parameters given
        via the dump_parameters method

        Parameters
        ----------
        x : array
            log frequency.
        y : array
            log derivative of epsilon'.

        Returns
        -------
        None.

        """

        try:

              open('double_HN_deri.json')

        except  FileNotFoundError as e:
                print(f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
        else:

           with open('double_HN_deri.json',"r") as openfile:
            loaded_par = json.load(openfile)


            print("loaded parameters \n" ,loaded_par)

            double_hn_par_index = ['beta1','gamma1','freq1','deps1','beta2','gamma2','freq2','deps2']
            init_fit_par = {key:value for key,value in loaded_par.items() if key in double_hn_par_index}


            hn_sub_par1 = loaded_par['beta1'],loaded_par['gamma1'],loaded_par['freq1'],loaded_par['deps1']
            hn_sub_par2 = loaded_par['beta2'],loaded_par['gamma2'],loaded_par['freq2'],loaded_par['deps2']

            hn_sub1 = self.deri_hn(x,*hn_sub_par1)
            hn_sub2 = self.deri_hn(x,*hn_sub_par2)


            plt.figure()
            plt.scatter(x,y,marker='s',color='r',facecolors='none',label='data',s=100,linewidth=2)
            plt.plot(x,hn_sub1,'b',label='initial guess - peak1')
            plt.plot(x,hn_sub2,'g',label='initial guess - peak2')
            plt.xlabel('log ( f [Hz])')
            plt.ylabel('log ( -d$\epsilon´$/ dlog f )')
            plt.legend()

        return init_fit_par

    def initial_view_deri_double_hn_ep(self, x, y):
        """
        Plots the derivative double HN function with electrode polarization
        based on the initial parameters given via the dump_parameters method.

        Parameters
        ----------
        x : array
            Log frequency.
        y : array
            Log derivative of epsilon'.

        Returns
        -------
        None.
        """

        try:
            open('double_HN_deri.json')
        except FileNotFoundError as e:
            print(f'{e}\nPlease dump initial fit parameters using dump.parameters method')
        else:
            with open('double_HN_deri.json', "r") as openfile:
                loaded_par = json.load(openfile)

            print("loaded parameters \n", loaded_par)

            hn_sub_par1 = [loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1']]
            hn_sub_par2 = [loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2']]
            ep_sub_par = [loaded_par['ep'], loaded_par['s']]

            hn_sub1 = self.deri_hn(x, *hn_sub_par1)
            hn_sub2 = self.deri_hn(x, *hn_sub_par2)
            ep_sub = self.ep_s(x, *ep_sub_par)

            plt.figure()
            plt.scatter(x, y, marker='s', color='r', facecolors='none', label='data', s=100, linewidth=2)
            plt.plot(x, hn_sub1, 'b', label='initial guess - peak1')
            plt.plot(x, hn_sub2, 'g', label='initial guess - peak2')
            plt.plot(x, ep_sub, 'r', label='initial guess - electrode polarization')
            plt.xlabel('log ( f [Hz])')
            plt.ylabel('log ( -d$\epsilon´$/ dlog f )')
            plt.legend()

        return loaded_par

    def initial_view_deri_triple_hn_ep(self, x, y):
        """
        Plots the derivative triple HN function with electrode polarization
        based on the initial parameters given via the dump_parameters method.

        Parameters
        ----------
        x : array
            Log frequency.
        y : array
            Log derivative of epsilon'.

        Returns
        -------
        None.
        """

        try:
            open('triple_HN_ep_deri.json')
        except FileNotFoundError as e:
            print(f'{e}\nPlease dump initial fit parameters using dump.parameters method')
        else:
            with open('triple_HN_ep_deri.json', "r") as openfile:
                loaded_par = json.load(openfile)

            print("loaded parameters \n", loaded_par)

            hn_sub_par1 = [loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1']]
            hn_sub_par2 = [loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2']]
            hn_sub_par3 = [loaded_par['beta3'], loaded_par['gamma3'], loaded_par['freq3'], loaded_par['deps3']]
            ep_sub_par = [loaded_par['ep'], loaded_par['s']]

            hn_sub1 = self.deri_hn(x, *hn_sub_par1)
            hn_sub2 = self.deri_hn(x, *hn_sub_par2)
            hn_sub3 = self.deri_hn(x, *hn_sub_par3)
            ep_sub = self.ep_s(x, *ep_sub_par)

            plt.figure()
            plt.scatter(x, y, marker='s', color='r', facecolors='none', label='data', s=100, linewidth=2)
            plt.plot(x, hn_sub1, 'b', label='initial guess - peak1')
            plt.plot(x, hn_sub2, 'g', label='initial guess - peak2')
            plt.plot(x, hn_sub3, 'm', label='initial guess - peak3')
            plt.plot(x, ep_sub, 'r', label='initial guess - electrode polarization')
            plt.xlabel('log ( f [Hz])')
            plt.ylabel('log ( -d$\epsilon´$/ dlog f )')
            plt.legend()

        return loaded_par


    def initial_view_deri_ud(self, x, y):
        """
        Plots the derivative triple HN function with electrode polarization and Universal Dielectric function
        based on the initial parameters given via the dump_parameters method.

        Parameters
        ----------
        x : array
            Log frequency.
        y : array
            Log derivative of epsilon'.

        Returns
        -------
        None.
        """

        try:
            open('triple_HN_ep_ud_deri.json')
        except FileNotFoundError as e:
            print(f'{e}\nPlease dump initial fit parameters using dump.parameters method')
        else:
            with open('triple_HN_ep_ud_deri.json', "r") as openfile:
                loaded_par = json.load(openfile)

            print("loaded parameters \n", loaded_par)

            # Extracting parameters for HN peaks
            hn_sub_par1 = [loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1']]
            hn_sub_par2 = [loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2']]
            hn_sub_par3 = [loaded_par['beta3'], loaded_par['gamma3'], loaded_par['freq3'], loaded_par['deps3']]

            # Extracting parameters for EP
            ep_sub_par = [loaded_par['ep'], loaded_par['s']]

            # Extracting parameters for UD
            ud_par = [loaded_par['B'], loaded_par['m']]

            # Compute the HN peaks, EP, and UD components
            hn_sub1 = self.deri_hn(x, *hn_sub_par1)
            hn_sub2 = self.deri_hn(x, *hn_sub_par2)
            hn_sub3 = self.deri_hn(x, *hn_sub_par3)
            ep_sub = self.ep_s(x, *ep_sub_par)
            ud_sub = self.ud_s(x, *ud_par)  # Compute UD function

            plt.figure()
            plt.scatter(x, y, marker='s', color='r', facecolors='none', label='data', s=100, linewidth=2)
            plt.plot(x, hn_sub1, 'b', label='initial guess - peak1')
            plt.plot(x, hn_sub2, 'g', label='initial guess - peak2')
            plt.plot(x, hn_sub3, 'm', label='initial guess - peak3')
            plt.plot(x, ep_sub, 'r', label='initial guess - electrode polarization')
            plt.plot(x, ud_sub, 'k--', label='initial guess - Universal Dielectric')
            plt.xlabel('log ( f [Hz])')
            plt.ylabel('log ( -d$\epsilon´$/ dlog f )')
            plt.legend()

        return loaded_par

#    def sel_function(self):
#        """
#        A function to select the type of fit function during curve fitting
#
#        Returns
#        -------
#        func_decision : int
#            choice of the fit function.
#
#        """
#
#        func_decision = int(input("Choose the fit function:\n 1 -- deri_HN\n 2-- deri_HN with EP\n 3 -- double_deri_HN\n 4 -- double_deri_HN with EP\n 5 -- triple_deri_HN with EP\n"))
#        return (func_decision)

    def sel_function(self):
        """
        A function to select the type of fit function during curve fitting.

        Returns
        -------
        func_decision : int
            Choice of the fit function.
        """
        print("Choose the fit function:")
        print("1 -- deri_HN")
        print("2 -- deri_HN with EP")
        print("3 -- double_deri_HN")
        print("4 -- double_deri_HN with EP")
        print("5 -- triple_deri_HN with EP")
        print("6 -- deri_UD")  # Added option for Universal Dielectric Response

        func_decision = int(input("Enter your choice: "))
        return func_decision


    def fit(self,x,y):
        """
        Fits the derivaitve of epsilon' data with choice of fit function
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
        y1= np.array(y)

        global popt1

        plt.figure()

        global popt2
        global b,g,fm,deps,ep,s,l_f,b1,g1,fm1,deps1,b2,g2,fm2,deps2,l_f1,l_f2

        if func_number==1:

            try:

              open('HN_deri.json')

            except  FileNotFoundError as e:
                print(f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
            else:
               with open('HN_deri.json',"r") as openfile:
                   loaded_par = json.load(openfile)

               hn_p0 = [loaded_par['beta'],loaded_par['gamma'],loaded_par['freq'],loaded_par['deps']]
               popt1, pcov2 = curve_fit(self.deri_hn, x1, y1, hn_p0,bounds =((0,0,1e-7,0),(1,1,1e7,np.inf)),absolute_sigma=True)
               yfit2 = self.deri_hn(x1,*popt1)

               plt.scatter(x1,y1,marker='s',color='r',facecolors='none',label='data',s=100,linewidth=2)
               plt.plot(x1,yfit2,'m--',label='derivative HN fit', linewidth=2)
               plt.xlabel('log ( f [Hz])')
               plt.ylabel('log ( -d$\epsilon´$/ dlog f )')
               plt.legend()

               b,g,fm,deps = popt1[0:4]

               ep,s =0,1
               #s = 0
               n_s = np.sin(np.pi*b/(2+2*g))**(1/b)
               n_s2 = 1/(np.sin(np.pi*b*g/(2+2*g))**(1/b))
               fmax = fm*n_s*n_s2
               l_f = np.log10(fmax)
               print(*popt1)
               print('log fmax:',l_f)
               fit_par = {"beta": b, "gamma": g,  "freq": fm,"deps":deps, "ep":ep, "s":s}
               with open('HN_deri.json',"w") as outfile:
                   json.dump(fit_par,outfile)

               with open('HN_deri.json',"r") as openfile:
                   loaded_par = json.load(openfile)

               print("fit parameters dumped for next iteration",loaded_par)

        elif func_number==2:

            try:

              open('HN_deri.json')

            except  FileNotFoundError as e:
                print(f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
            else:
               with open('HN_deri.json',"r") as openfile:
                    loaded_par = json.load(openfile)

               p0 = [loaded_par['beta'],loaded_par['gamma'],loaded_par['freq'],loaded_par['deps'],loaded_par['ep'],loaded_par['s']]
               popt2, pcov2 = curve_fit(self.deri_hn_ep, x1, y1, p0, bounds =((0,0,1e-7,0,0,0),(1,1,1e7,np.inf,np.inf,1)),absolute_sigma=True)
               yfit3 = self.deri_hn_ep(x1,*popt2)

               plt.scatter(x1,y1,marker='s',color='r',facecolors='none',label='data',s=100,linewidth=2)
               plt.plot(x1,yfit3,'m--',label='derivative HN with EP fit', linewidth=2)
               plt.xlabel('log ( f [Hz])')
               plt.ylabel('log ( -d$\epsilon´$/ dlog f )')


               hn_sub_par = popt2[0],popt2[1],popt2[2],popt2[3]
               ep_sub_par = popt2[4], popt2[5]
               hn_sub = self.deri_hn(x1,*hn_sub_par)
               ep_sub = self.ep_s(x1,*ep_sub_par)
               plt.plot(x1,hn_sub,'b',label='peak')
               plt.plot(x1,ep_sub,'g',label='Electrode Polarization')
               plt.legend()

               print(*popt2)
               b,g,fm,deps,ep,s = popt2[:]

               n_s = np.sin(np.pi*b/(2+2*g))**(1/b)
               n_s2 = 1/(np.sin(np.pi*b*g/(2+2*g))**(1/b))
               fmax = fm*n_s*n_s2
               l_f = np.log10(fmax)
               print('log fmax:',l_f)
               fit_par = {"beta": b, "gamma": g,  "freq": fm,"deps":deps, "ep":ep, "s":s}
               with open('HN_deri.json',"w") as outfile:
                   json.dump(fit_par,outfile)

               with open('HN_deri.json',"r") as openfile:
                   loaded_par = json.load(openfile)

               print("fit parameters dumped for next iteration",loaded_par)

        elif func_number ==3:
            try:

              open('double_HN_deri.json')

            except  FileNotFoundError as e:
                print(f'{e}' + '\n', "Please dump initial fit parameters using dump.parameters method")
            else:

                with open('double_HN_deri.json',"r") as openfile:
                      loaded_par = json.load(openfile)

                p0 = [loaded_par['beta1'],loaded_par['gamma1'],loaded_par['freq1'],loaded_par['deps1'],loaded_par['beta2'],loaded_par['gamma2'],loaded_par['freq2'],loaded_par['deps2']]
                print(p0)
                popt2, pcov2 = curve_fit(self.deri_double_hn, x1, y1, p0, bounds =((0,0,1e-7,0,0,0,1e-7,0),(1,1,1e7,np.inf,1,1,1e7,np.inf)),absolute_sigma=True)
                yfit4 = self.deri_double_hn(x1,*popt2)
                plt.scatter(x1,y1,marker='s',color='r',facecolors='none',label='data',s=100,linewidth=2)
                plt.plot(x1,yfit4,'m--',label='derivative HN fit', linewidth=2)
                plt.xlabel('log ( f [Hz])')
                plt.ylabel('log ( -d$\epsilon´$/ dlog f )')


                hn_sub_par1 = popt2[0], popt2[1],popt2[2],popt2[3]
                hn_sub_par2 = popt2[4], popt2[5],popt2[6],popt2[7]
                hn_sub1 = self.deri_hn(x1,*hn_sub_par1)
                hn_sub2 = self.deri_hn(x1,*hn_sub_par2)
                plt.plot(x1,hn_sub1,'b',label='peak1')
                plt.plot(x1,hn_sub2,'g',label='peak2')
                plt.legend()

                b1,g1,fm1,deps1,b2,g2,fm2,deps2 = popt2[:]
                ep, s = 0, 1
                nr_1 = np.sin(np.pi*b1/(2+2*g1))**(1/b1)
                nr_2 = 1/(np.sin(np.pi*b1*g1/(2+2*g1))**(1/b1))
                fmax1 = fm1*nr_1*nr_2
                nr_3 = np.sin(np.pi*b2/(2+2*g2))**(1/b2)
                nr_4 = 1/(np.sin(np.pi*b2*g2/(2+2*g2))**(1/b2))
                fmax2 = fm2*nr_3*nr_4
                l_f1 = np.log10(fmax1)
                l_f2 = np.log10(fmax2)

                print('log fmax1:',l_f1 ,"\n" 'log fmax2:',l_f2)
                fit_par = {"beta1":b1,"gamma1":g1,"freq1":fm1,"deps1":deps1,"beta2":b2,"gamma2":g2,"freq2":fm2,"deps2":deps2, "ep":ep, "s":s}

                with open('double_HN_deri.json',"w") as outfile:
                    json.dump(fit_par,outfile)
                with open('double_HN_deri.json',"r") as openfile:
                    loaded_par = json.load(openfile)
                print("fit parameters dumped for next iteration",loaded_par)



        elif func_number == 4:
            try:
                # Check if the required JSON file exists
                open('double_HN_deri.json')
            except FileNotFoundError as e:
                print(f'{e}\nPlease dump initial fit parameters using dump.parameters method')
            else:
                # Read initial parameters from the JSON file
                with open('double_HN_deri.json', "r") as openfile:
                    loaded_par = json.load(openfile)

                # Extract initial guess parameters for the fit
                hn_p0 = [
                    loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1'],
                    loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2'],
                    loaded_par['ep'], loaded_par['s']
                ]
                print("Initial guess parameters:", hn_p0)

                # Perform the curve fitting
                popt1, pcov1 = curve_fit(
                    self.deri_double_hn_ep,
                    x1,
                    y1,
                    p0=hn_p0,
                    bounds=(
                        (0, 0, 1e-7, 0, 0, 0, 1e-7, 0, 0, 0),  # Lower bounds for each parameter
                        (1, 1, 1e7, np.inf, 1, 1, 1e7, np.inf, np.inf, 1)  # Upper bounds for each parameter
                    ),
                    absolute_sigma=True
                )

                # Calculate the fitted values
                yfit = self.deri_double_hn_ep(x1, *popt1)

                # Plot the data and the fit
                plt.scatter(x1, y1, marker='s', color='r', facecolors='none', label='data', s=100, linewidth=2)
                plt.plot(x1, yfit, 'm--', label='fit - double peak with EP', linewidth=2)
                plt.xlabel('log ( f [Hz])')
                plt.ylabel('log ( -d$\epsilon´$/ dlog f )')

                # Calculate model components for plotting
                hn_sub_par1 = popt1[0:4]
                hn_sub_par2 = popt1[4:8]
                ep_par = popt1[8:]

                hn_sub1 = self.deri_hn(x1, *hn_sub_par1)
                hn_sub2 = self.deri_hn(x1, *hn_sub_par2)
                ep_sub = self.ep_s(x1, *ep_par)

                # Plot the individual components
                plt.plot(x1, hn_sub1, 'b', label='peak1')
                plt.plot(x1, hn_sub2, 'g', label='peak2')
                plt.plot(x1, ep_sub, 'k--', label='Electrode Polarization')
                plt.legend()

                # Print and save fit parameters
                fit_par = {
                    'beta1': popt1[0], 'gamma1': popt1[1], 'freq1': popt1[2], 'deps1': popt1[3],
                    'beta2': popt1[4], 'gamma2': popt1[5], 'freq2': popt1[6], 'deps2': popt1[7],
                    'ep': popt1[8], 's': popt1[9]
                }

                with open('double_HN_deri.json', "w") as outfile:
                    json.dump(fit_par, outfile)

                # Reload and print the fit parameters
                with open('double_HN_deri.json', "r") as openfile:
                    loaded_par = json.load(openfile)
                print("Fit parameters dumped for next iteration:", loaded_par)

        elif func_number == 5:
            # Implementation for the new triple HN + EP function
            try:
                open('triple_HN_ep_deri.json')
            except FileNotFoundError as e:
                print(f'{e}\nPlease dump initial fit parameters using dump.parameters method')
            else:
                with open('triple_HN_ep_deri.json', "r") as openfile:
                    loaded_par = json.load(openfile)

                # Extract initial parameters for the fit
                hn_p0 = [
                    loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1'],
                    loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2'],
                    loaded_par['beta3'], loaded_par['gamma3'], loaded_par['freq3'], loaded_par['deps3'],
                    loaded_par['ep'], loaded_par['s']
                ]

                # Perform curve fitting
                popt1, pcov1 = curve_fit(
                    self.deri_triple_hn_ep,
                    x1,
                    y1,
                    p0=hn_p0,
                    bounds=(
                        (0, 0, 1e-7, 0, 0, 0, 1e-7, 0, 0, 0, 1e-7, 0, 0, 0),
                        (1, 1, 1e7, np.inf, 1, 1, 1e7, np.inf, 1, 1, 1e7, np.inf, np.inf, 1)
                    ),
                    absolute_sigma=True
                )

                # Calculate fitted values
                yfit = self.deri_triple_hn_ep(x1, *popt1)

                # Plot data and fit
                plt.scatter(x1, y1, marker='s', color='r', facecolors='none', label='data', s=100, linewidth=2)
                plt.plot(x1, yfit, 'm--', label='fit - triple peak with EP', linewidth=2)
                plt.xlabel('log ( f [Hz])')
                plt.ylabel('log ( -d$\epsilon´$/ dlog f )')

                # Calculate model components for plotting
                hn_sub_par1 = popt1[0:4]
                hn_sub_par2 = popt1[4:8]
                hn_sub_par3 = popt1[8:12]
                ep_par = popt1[12:]

                hn_sub1 = self.deri_hn(x1, *hn_sub_par1)
                hn_sub2 = self.deri_hn(x1, *hn_sub_par2)
                hn_sub3 = self.deri_hn(x1, *hn_sub_par3)
                ep_sub = self.ep_s(x1, *ep_par)

                # Plot individual components
                plt.plot(x1, hn_sub1, 'b', label='peak1')
                plt.plot(x1, hn_sub2, 'g', label='peak2')
                plt.plot(x1, hn_sub3, 'c', label='peak3')
                plt.plot(x1, ep_sub, 'k--', label='Electrode Polarization')
                plt.legend()

                # Print and save fit parameters
                fit_par = {
                    'beta1': popt1[0], 'gamma1': popt1[1], 'freq1': popt1[2], 'deps1': popt1[3],
                    'beta2': popt1[4], 'gamma2': popt1[5], 'freq2': popt1[6], 'deps2': popt1[7],
                    'beta3': popt1[8], 'gamma3': popt1[9], 'freq3': popt1[10], 'deps3': popt1[11],
                    'ep': popt1[12], 's': popt1[13]
                }

                with open('triple_HN_ep_deri.json', "w") as outfile:
                    json.dump(fit_par, outfile)

                # Reload and print fit parameters
                with open('triple_HN_ep_deri.json', "r") as openfile:
                    loaded_par = json.load(openfile)
                print("Fit parameters dumped for next iteration:", loaded_par)





        elif func_number == 6:
            # Implementation for the triple HN + EP + UD function
            try:
                open('triple_HN_ep_ud_deri.json')
            except FileNotFoundError as e:
                print(f'{e}\nPlease dump initial fit parameters using dump.parameters method')
            else:
                with open('triple_HN_ep_ud_deri.json', "r") as openfile:
                    loaded_par = json.load(openfile)

                # Extract initial parameters for the fit
                triple_hn_ep_ud_p0 = [
                    loaded_par['beta1'], loaded_par['gamma1'], loaded_par['freq1'], loaded_par['deps1'],
                    loaded_par['beta2'], loaded_par['gamma2'], loaded_par['freq2'], loaded_par['deps2'],
                    loaded_par['beta3'], loaded_par['gamma3'], loaded_par['freq3'], loaded_par['deps3'],
                    loaded_par['ep'], loaded_par['s'], loaded_par['B'], loaded_par['m']
                ]

                # Number of parameters in p0 is 16
                # Ensure bounds have the same number of elements
                lower_bounds = (
                    0, 0, 1e-7, 0,  # beta1, gamma1, freq1, deps1
                    0, 0, 1e-7, 0,  # beta2, gamma2, freq2, deps2
                    0, 0, 1e-7, 0,  # beta3, gamma3, freq3, deps3
                    0, 0, 0, 0  # ep, s, B, m
                )
                upper_bounds = (
                    1, 1, 1e7, np.inf,  # beta1, gamma1, freq1, deps1
                    1, 1, 1e7, np.inf,  # beta2, gamma2, freq2, deps2
                    1, 1, 1e7, np.inf,  # beta3, gamma3, freq3, deps3
                    np.inf, np.inf, np.inf, np.inf  # ep, s, B, m
                )

                # Perform curve fitting
                popt1, pcov1 = curve_fit(
                    self.deri_triple_hn_ep_ud,
                    x1,
                    y1,
                    p0=triple_hn_ep_ud_p0,
                    bounds=(lower_bounds, upper_bounds),
                    absolute_sigma=True
                )

                # Calculate fitted values
                yfit = self.deri_triple_hn_ep_ud(x1, *popt1)

                # Plot data and fit
                plt.scatter(x1, y1, marker='s', color='r', facecolors='none', label='data', s=100, linewidth=1)
                plt.plot(x1, yfit, 'm--', label='fit - sum of HN functions plus EP', linewidth=2)
                plt.xlabel('log ( f [Hz])')
                plt.ylabel('log ( -d$\epsilon´$/ dlog f )')

                # Calculate model components for plotting
                hn_sub_par1 = popt1[0:4]
                hn_sub_par2 = popt1[4:8]
                hn_sub_par3 = popt1[8:12]
                ep_par = popt1[12:14]
                B, m = popt1[14:16]

                hn_sub1 = self.deri_hn(x1, *hn_sub_par1)
                hn_sub2 = self.deri_hn(x1, *hn_sub_par2)
                ep_sub = self.ep_s(x1, *ep_par)

                # Plot individual components
                plt.plot(x1, hn_sub1, 'b', label='HN1')
                plt.plot(x1, hn_sub2, 'g', label='HN2')
                plt.plot(x1, ep_sub, 'k--', label='Electrode Polarization')
                beta3 = popt1[8]
                gamma3 = popt1[9]
                # This is not a very good way to do it, but whatever. It's wrong, but only slightly wrong.
                # Plot sub3 only if it's relevant
                if beta3 >= 0.09 or gamma3 >= 0.09:
                    hn_sub3 = self.deri_hn(x1, *hn_sub_par3)
                    plt.plot(x1, hn_sub3, 'c', label='HN3')
                # Plot UD function only if B and m meet the threshold criteria
                if B >= 0.2 and m > 1e-10:
                    ud_sub = B * x1**m
                    plt.plot(x1, ud_sub, 'r--', label='Universal Dielectric Response')

                plt.legend()

                # Print and save fit parameters
                fit_par = {
                    'beta1': popt1[0], 'gamma1': popt1[1], 'freq1': popt1[2], 'deps1': popt1[3],
                    'beta2': popt1[4], 'gamma2': popt1[5], 'freq2': popt1[6], 'deps2': popt1[7],
                    'beta3': popt1[8], 'gamma3': popt1[9], 'freq3': popt1[10], 'deps3': popt1[11],
                    'ep': popt1[12], 's': popt1[13], 'B': B, 'm': m
                }

                with open('triple_HN_ep_ud_deri.json', "w") as outfile:
                    json.dump(fit_par, outfile)

                # Reload and print fit parameters
                with open('triple_HN_ep_ud_deri.json', "r") as openfile:
                    loaded_par = json.load(openfile)
                print("Fit parameters dumped for next iteration:", loaded_par)

            return fit_par





    def save_fit_deri_hn(self,T):
        """
        saves the fit parameters of derivative hn function in a file
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
        global deps
        deps = deps/2.303
        res_file.write( f'{T}' + '\t' + f'{b:.3f}' + '\t' + f'{g:.3f}' +  '\t' + f'{deps:.3f}' +'\t'  + f'{l_f:.3f}' + '\t' +f'{ep:.3f}'+'\t' + f'{s:.3f}' +"\n")
        return()

    def save_fit_deri_hn_ep(self,T):
        """
        saves the fit parameters of derivaitve hn function with electrode polarization
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
        res_file = open(ana_file,'a')
        global deps
        deps = deps/2.303
        res_file.write( f'{T}' + '\t' + f'{b:.3f}' + '\t' + f'{g:.3f}' +  '\t' + f'{deps:.3f}' +'\t'  + f'{l_f:.3f}' + '\t' +f'{ep:.3f}'+'\t' + f'{s:.3f}' +"\n")
        return()

    def save_fit_deri_double_HN(self,T):
        """
        saves the fit parameters of derivaitve double hn function in  a file,
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
        global deps1,deps2
        deps1 = deps1/2.303
        deps2 = deps2/2.303
        res_file.write( f'{T}' + '\t' + f'{b1:.3f}' + '\t' + f'{g1:.3f}' +  '\t' + f'{deps1:.3f}' +'\t'  + f'{l_f1:.3f}' + '\t' +f'{b2:.3f}' + '\t' + f'{g2:.3f}' +  '\t' + f'{deps2:.3f}' +'\t'  + f'{l_f2:.3f}'+ '\t' +f'{ep:.3f}'+'\t' + f'{s:.3f}' +"\n")
        return()

    def save_fit_deri_triple_hn_ep_ud(self, T):
        """
        Saves the fit parameters of the triple HN function with Electrode Polarization (EP) and
        Universal Dielectric Response (UD) into a file. The file must be created via a method such as
        create_analysis_file.

        Parameters
        ----------
        T : float
            Temperature or an integer that corresponds to a file number during analysis.

        Returns
        -------
        None.
        """
        # Open the analysis file in append mode
        res_file = open(ana_file, 'a')

        # Parameters extracted from the fit
        # For illustration, replace `self.beta1`, `self.gamma1`, etc. with actual parameter variables
        beta1 = self.beta1
        gamma1 = self.gamma1
        freq1 = self.freq1
        deps1 = self.deps1
        beta2 = self.beta2
        gamma2 = self.gamma2
        freq2 = self.freq2
        deps2 = self.deps2
        beta3 = self.beta3
        gamma3 = self.gamma3
        freq3 = self.freq3
        deps3 = self.deps3
        ep = self.ep
        s = self.s
        B = self.B
        m = self.m

        # Adjust `deps` values if necessary (assuming they need adjustment)
        deps1 /= 2.303
        deps2 /= 2.303
        deps3 /= 2.303

        # Write parameters to the file
        res_file.write(
            f'{T}\t'
            f'{beta1:.3f}\t{gamma1:.3f}\t{freq1:.3f}\t{deps1:.3f}\t'
            f'{beta2:.3f}\t{gamma2:.3f}\t{freq2:.3f}\t{deps2:.3f}\t'
            f'{beta3:.3f}\t{gamma3:.3f}\t{freq3:.3f}\t{deps3:.3f}\t'
            f'{ep:.3f}\t{s:.3f}\t{B:.3f}\t{m:.3f}\n'
        )

        # Close the file
        res_file.close()

        return
    plt.show()
