# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 18:11:58 2022

@author: mkolmang
"""

def test_jonscher():
    import numpy as np
    import pandas as pd
    import HavNegpy as hn
    #import os
    import json
    import matplotlib.pyplot as plt
    cond = hn.Conductivity()
    
    
    
    filename = 'cond_example_data.txt'
    fit_data_filename = 'fitted_jonscher.txt'
    col_names = ['f','s1','s2']
    col_fit = ['f','fit_y']
    df = pd.read_csv(filename, sep='\t',index_col=False,usecols = [0,1,2],names=col_names,header=None,skiprows=10,encoding='unicode_escape',engine='python')
    df2 = pd.read_csv(fit_data_filename, sep='\t',index_col=False,usecols = [0,1],names=col_fit,header=None,skiprows=1,encoding='unicode_escape',engine='python')
    
    x = np.log10(df['f'])
    y = np.log10(df['s1'])
    
    #cond.dump_parameters()
    #x1,y1 = cond.select_range(x, y)
    #cond.fit(x1,y1)
    with open('cond.json',"r") as f:
       loaded_par = json.load(f)
    p0 = [loaded_par['fc'], loaded_par['DC'],loaded_par['n']]
    actual_y = cond.jonscher_function(x,*p0)
    expected_y = np.array(df2['fit_y'].values)
    np.testing.assert_array_almost_equal(actual_y, expected_y,decimal=4)
    plt.plot(x,actual_y)
    plt.plot(x,expected_y)
    
    return None

    '''
    actual_y = cond.jonscher_function(x, 3, 0, 1)
    expected_y = np.array(df['s'].values)
    #print(expected_y)
    #print(actual_y)
    #print(df)
    #np.testing.assert_array_equal(actual_y, expected_y,err_msg = 'equal')
    np.testing.assert_array_almost_equal(actual_y, expected_y)
    plt.plot(x,actual_y)
    plt.plot(x,expected_y)
    plt.show()
    '''
    
    
    