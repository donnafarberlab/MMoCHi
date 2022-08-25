import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable

from .logger import logg
from . import utils

def run_threshold(markname,adata,data_key,thresh):
    '''Lightweight wrapper to run utils.get_data(), then perform pos/neg/? thresholding, given a list of thresholds to apply.
    markname: str, name of the marker
    adata: anndata object
    data_key: str, name of the key in .obsm[] to look for when searching for markname
    thresh: tuple of numerics, max and minimum thresholds to apply for positive and negative populations
    '''
    data = utils.get_data(adata,markname,data_key)
    return _run_threshold(data, thresh)

def _run_threshold(data, thresh):
    '''Using the thresholds provided, perform thresholding
    data: np array of numerics
    thresh: listlike of numerics, the max and min of this list are used as thresholds
    returns an np.array of "?" the length of data, with 
    "pos" if data is greater than or equal to the maximum threshold, or "neg" if less than the minimum threshold'''
    thresholded = np.repeat('?', len(data)).astype(object)
    thresholded[data >= max(thresh)] = 'pos'
    thresholded[data <= min(thresh)] = 'neg'
    return thresholded

def _plot_threshold(data, thresh, msw, markname_full,title_addition=''):
    '''Using the thresholds provided, perform thresholding
    data: np array of numerics
    thresh: listlike of numerics, the max and min of this list are used as thresholds
    msw: tuple of  m - the means of peaks in the gaussian mixture model, 
                   s - the standard deviation of these peaks, 
                   w - the weights (heights) of these peaks
    Plots the histogram of the thresholds, and if msw is provided, plots the gaussian mixture model as well
    '''
    mask = data > 0
    data = data[mask]
    percent_included = round(sum(mask) / len(mask) * 100,2)
    fig = plt.figure(figsize=(8,2),dpi=75)
    hist = plt.hist((data[data>max(thresh)],data[data<min(thresh)],
                     data[(data>min(thresh)) & (data<max(thresh))]),
                     bins=100,stacked=True,density=True)
    if not msw is None:
        (m,s,w) = msw
        x = np.linspace(0,max(data), 10000)
        plt.plot(x, norm.pdf(x, m[0], s[0])*w[0],c='orange')
        if len(m) > 1:
            loc = np.linspace(m[0],m[1], 10000)
            comb = sum([norm.pdf(loc, mm, ss)*ww for mm, ss, ww in zip(m, s, w)])
        plt.plot(x, norm.pdf(x, m[-1], s[-1])*w[-1],c='blue')
        if len(m) > 2:
            plt.plot(x,norm.pdf(x, m[1], s[1])*w[1],c='green')
        if len(m) > 1:
            plt.plot(loc, comb,c='black')
    
    plt.ylim(top = max((hist[0]).flatten())*1.2)
    plt.xlim(0,max(data))
    plt.subplots_adjust(left=0.05, right=0.9, top=0.9, bottom=0.05)
    plt.margins(x=0)
    plt.yticks([])
    plt.xticks(np.around(np.arange(min(data), max(data), step=(max(data)-min(data))/10),2))
    plt.grid(axis = 'x')
    plt.title(f'''{markname_full} {title_addition}\nShowing {percent_included}% of events''')
    plt.show()
    return

def _interactive_threshold(thresh):
    '''Given a threshold suggestion (thresh), attempts to ask the user for thresholds using the input function
    invalid inputs (such as letters or no input) default the threshold to the suggested threshold.
    returns tuple of numerics for the upper and lower thresholds'''
    upper_thresh = input("What should the upper threshold be?")
    try:
        upper_thresh = float(upper_thresh)
    except:
        upper_thresh = max(thresh)
    lower_thresh = input("What should the lower threshold be?")
    try:
        lower_thresh = float(lower_thresh)
    except:
        lower_thresh = min(thresh)
    return (upper_thresh,lower_thresh)

def _fancy_interactive_threshold(thresh,maximum):
    '''Function to generate sliders beneath a donor for upper and lower thresholds.
    thresh: tuple of numerics for the upper and lower thresholds.
    maxiumum: int, the maximum value of the marker, serves as the maximum on the slider.
    Returns these sliders so that the fancy resolver may read them eventually.
    The 590px size was determined to be the correct size to match the plotted output of _plot_threshold'''
    try:
        from ipywidgets import FloatSlider, Layout
        from IPython.display import display
    except ImportError:
        raise ImportError('Please install ipywdigets using pip install ipywdigets')
    
    slider1 = FloatSlider(value= max(thresh),min=0,max=maximum,step=0.1,readout=True,layout=Layout(width='600px'))
    slider2 = FloatSlider(min(thresh),min=0,max=maximum,step=0.1,readout=True,layout=Layout(width='600px'))
    display(slider1)
    display(slider2)
    return (slider1,slider2)

def _fancy_resolver(promise):
    '''Reads in the value from a tuple of slider "promises" 
    Returns a sorted tuple of thresholds'''
    return tuple(sorted((promise[0].value,promise[1].value)))

def threshold(markname,adata,data_key=None,preset_threshold=None,
              include_zeroes=False,n=0,force_model = False,
              plot=True,title_addition='',interactive=True,fancy=False, run=False):
    '''
    Bimodal thresholding for markers.
    Runs utils.get_data() to identify the data.
    Calculates optimal thresholds
    plot, bool: whether to display plots
    fancy, interactive: bools, Allow the user to input their thresholds (via "fancy" float sliders or the "input" function)
    run, whether to run the thresholding to return thresholded data
    '''
    data, thresh, msw, markname_full = _calc_threshold(markname,adata,data_key,n,include_zeroes,preset_threshold=preset_threshold, force_model=force_model)
    if plot:
        _plot_threshold(data, thresh, msw, markname_full,title_addition=title_addition)
    if fancy:
        return _fancy_interactive_threshold(thresh, max(data))
    if interactive:
        if not plot:
            print(markname_full)
        thresh = _interactive_threshold(thresh)
    if run:
        thresholded = _run_threshold(data, thresh)
        return thresholded, thresh
    else:
        return thresh
    
def _calc_threshold(markname,adata,data_key=None,n=0,include_zeroes=False,preset_threshold=None, force_model=False):
    '''
        Internal function to handle threshold making
        gets the data using utils.get_data, determines if 0's should be kept in the data, then creatures a gaussian mixture model of
        either 1 or 2 compoenents. The fit of these models are compared and if only use the 2 component mix if it fits significantly better than 
        the 1 component model. Calculate means, weights, and standard deviations of the gaussians, then calculates thresholds as the minimum and 
        maximum of the 1st peak + 2-4 std, the 2nd peak + 2-4 std (if there're two peaks), and .05 (if there are no peaks)
        markname: str, name of the marker
        adata: anndata object
        data_key: str, name of the key in .obsm[] to look for when searching for markname
        n: int, 0, 1, or 2, to define how many thresholds should be defined
        include_zeroes: whether to include zeroes when calculating the gaussian mixture
        preset_threshold: tuple of numerics, max and minimum thresholds to apply for positive and negative populations, overrides calculation
        force_model: whether to force the creation of a gaussian mixture model. If preset thresholds are defined, this model is only used to create msw
    '''
    data, markname_full = utils.get_data(adata,markname,data_key,return_source=True)

    if "_gex" in markname_full or include_zeroes:
        include_zeros=True
        mask = ~pd.Series(data).isna()
    else:
        mask = data>0
    if force_model or (preset_threshold is None and not "_gex" in markname_full):
        masked_data = np.array(data[mask])
        mix1 = GaussianMixture(n_components=1,random_state=20).fit(masked_data.reshape(-1,1))
        mix2 = GaussianMixture(n_components=2,random_state=20).fit(masked_data.reshape(-1,1))
        if n == 1:
             mix, stds_needed, n = mix1, 4, 1
        elif n == 2:
            mix, stds_needed, n = mix2, 2, 2
        elif mix1.lower_bound_ < (mix2.lower_bound_ - 0.025):
            mix, stds_needed, n = mix2, 2, 2
        else:
            mix, stds_needed, n = mix1, 4, 1

        means, weights = sum(mix.means_.tolist(),[]), mix.weights_
        std = [ np.sqrt(  np.trace(i)/2) for i in mix.covariances_ ]
        m,s,w = zip(*sorted(zip(means,std,weights)))
        
        thresh = []
        thresh.append(m[0]+stds_needed*s[0])
        if n > 1:
            thresh.append(m[1]-stds_needed*s[1])
        elif n < 1:
            thresh.append(0.05)
    elif "_gex" in markname_full:
        thresh = [0,0.5]
    if not preset_threshold is None:
        thresh = list(preset_threshold)
        
    if force_model or (preset_threshold is None and not "_gex" in markname_full):
        return data, (max(thresh),min(thresh)), (m,s,w), markname_full
    else:
        return data, (max(thresh),min(thresh)), None, markname_full