import numpy as np
import pandas as pd
import anndata
import matplotlib.pyplot as plt
import scipy
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable
import ipywidgets

from .logger import logg
from . import utils

def run_threshold(markname: str, adata: anndata.AnnData,
                  data_key: Optional[str], thresh: Tuple[Union[int, float]]) -> np.array:
    '''
    Lightweight wrapper to run utils.get_data(), then perform pos/neg/? thresholding, given a list of thresholds to apply.
    
    Parameters
    ----------
    markname
        Name of the marker
    adata
        AnnData object containing normalized (and possibly batch corrected) events for thresholding
    data_key
        Name of the key in .obsm[] to look for when searching for markname
    thresh
        Upper and lower thresholds to apply for positive and negative populations
    
    Returns
    -------
    np.array
        length of data, with "pos" if data is >= to the upper threshold, or "neg" if < lower threshold, ? if in between upper and lower
    '''
    data = utils.get_data(adata,markname,data_key)
    return _run_threshold(data, thresh)

def _run_threshold(data: np.array, thresh: Sequence[Union[int,float]]) -> np.array:
    '''
    Using the thresholds provided, perform thresholding
    
    Parameters
    ----------
    data
        Data to threshold
    thresh
        The max and min of this list are used as thresholds
    returns
    -------
    np.array
        length of data, with "pos" if data is >= to the upper threshold, or "neg" if < lower threshold, ? if in between upper and lower
    '''
    thresholded = np.repeat('?', len(data)).astype(object)
    thresholded[data >= max(thresh)] = 'pos'
    thresholded[data <= min(thresh)] = 'neg'
    return thresholded

def _plot_threshold(data: np.array, thresh: Sequence[Union[float,int]],
                    msw: Tuple[Tuple[float], Tuple[float],Tuple[float]], markname_full: str,
                    title_addition: str='') -> None:
    '''
    Using the thresholds provided, perform and plot thresholding
    
    Parameters
    ----------
    data
        Data to threshold and display
    thresh
        The max and min of this list are used as thresholds
    msw
        m: tuple of floats
            the means of peaks in the gaussian mixture model
        s: tuple of floats
            the standard deviation of these peaks, 
        w: tuple of floats
            the weights (heights) of these peaks
    markname_full
        Name of the marker, to be used in the title of the graph
    title_addition
        Other information you would like to include in the title after the marker name
        
    Returns
    -------
    None, plots the histogram of the thresholds, and if msw is provided, plots the gaussian mixture model as well
    '''
    mask = data > 0
    data = data[mask]
    percent_included = round(sum(mask) / len(mask) * 100,2)
    fig = plt.figure(figsize=(8,2),dpi=75)
    hist = plt.hist((data[data>max(thresh)],data[data<min(thresh)],
                     data[(data>min(thresh)) & (data<max(thresh))]),
                     bins=100,stacked=True,density=True, color=['#ff7f0e',  '#1f77b4', '#2ca02c']) 
    if not msw is None:
        (m,s,w) = msw
        x = np.linspace(0,max(data), 10000)
        plt.plot(x, norm.pdf(x, m[0], s[0])*w[0],c='blue')
        if len(m) > 1:
            loc = np.linspace(m[0],m[1], 10000)
            comb = sum([norm.pdf(loc, mm, ss)*ww for mm, ss, ww in zip(m, s, w)])
        plt.plot(x, norm.pdf(x, m[-1], s[-1])*w[-1],c='orange')
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

def _interactive_threshold(thresh: Tuple[float, float]) -> Tuple[float, float]:
    '''
    Given a threshold suggestion (thresh), attempts to ask the user for thresholds using the input function.
    Invalid inputs (such as letters or no input) default the threshold to the suggested threshold.
    
    Parameters
    ----------
    thresh
        Default uppper and lower thresholds to use
    
    Returns
    -------
    New upper and lower thresholds
    '''
    if thresh is None or not isinstance(thresh, tuple) or not len(thresh) == 2:
        thresh = tuple(0,0)
    upper_thresh = input(f"Current upper threshold is: {max(thresh)} What should the upper threshold be?")
    try:
        upper_thresh = float(upper_thresh)
    except:
        logg.warn(f"Invalid input, defaulting to threshold of {max(thresh)}")
        upper_thresh = max(thresh)
    lower_thresh = input(f"Current lower threshold is: {min(thresh)} What should the lower threshold be?")
    try:
        lower_thresh = float(lower_thresh)
    except:
        logg.warn(f"Invalid input, defaulting to threshold of {min(thresh)}")
        lower_thresh = min(thresh)
    return (upper_thresh,lower_thresh)

def _fancy_interactive_threshold(thresh: Tuple[Union[int, float]], maximum:Union[int,float]) -> Tuple[ipywidgets.widgets.widget_float.FloatSlider, ipywidgets.widgets.widget_float.FloatSlider]:
    '''
    Function to generate sliders beneath a donor for upper and lower thresholds.
    
    Parameters
    ----------
    thresh: tuple of numerics
        Upper and lower thresholds
    maxiumum: numeric
        The maximum value of the marker, serves as the maximum on the slider
        
    Returns
    -------
    Returns these sliders so that the fancy resolver may read them eventually.
    The 590px size was determined to be the correct size to match the plotted output of _plot_threshold
    '''
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

def _fancy_resolver(promise: Tuple[float]) -> Tuple[float,float]:
    '''
    Reads in the value from a tuple of slider "promises"
    
    Parameters
    ----------
    promise: tuple of numeric
        slider promise values
        
    Returns
    -------
        A sorted tuple of thresholds
    '''
    return tuple(sorted((promise[0].value,promise[1].value)))

def threshold(markname: str, adata: anndata.AnnData,
              data_key: str=utils.DATA_KEY, preset_threshold: Tuple[Union[int, float]]=None,
              include_zeroes: bool=False, n: int=0,
              force_model: bool=False, plot: bool=True,
              title_addition: str='', interactive: bool=True,
              fancy: bool=False,
              run: bool=False) -> Union[Tuple[float,float], Tuple[Tuple[float,float], List[str]]]:
    '''
    Performs thresholding for markers, and optionally returns thresholds as well as thresholded events.
    Thresholded events are returned as a list of "pos" for positive (above the higher threshold), "neg" for 
    negative (below the lower threshold), and "?" for undefined (between the two thresholds).
    
    Runs utils.get_data() to identify the data.
    Automatically calculates threshold and has options for display of these thresholds on a histogram.
    
    Parameters
    ----------
    markname
        Name of the marker being used. Markers can also have "_lo" or "_hi" appended to them to specify
        multiple thresholds on the same marker at the same level.
    adata
        AnnData object for creating thresholds, containing expression in the .X and/or the .obsm[data_key] 
    data_key
        Name of the key in .obsm to look for when searching for markname
    present_threshold
        Max and minimum thresholds to apply for positive and negative populations, overrides calculation
    included_zeroes
        Whether to include zeroes when calculating the gaussian mixture
    n
        Determines the number of gaussians to fit. Can be 1, 2, or 3.
    force_model
        Whether to force the creation of a gaussian mixture model. If preset thresholds are defined, this model is only used to create msw
    plot
        Whether to display plots from _plot_threshold
    title_addition
        Other information you would like to include in the title after the marker name
    interactive
        Whether to prompt the user to enter thresholds (uses defaults if invalid entry)
    fancy
        Whether to create adjustable float sliders to enter thresholds. Returns a list of float sliders. This should only be used by the 
        run_all_thresholds command.
    run
        Whether to run the thresholding to return thresholded data
    Returns
    -------
    thresh Tuple of numerics
        Threshold values (numerics)
    thresholded:
        if run == true returns thresholded data, else only returns thresh
    
    '''
    if markname.endswith("_lo") or markname.endswith("_hi"):
        ending = markname[-3:]
        markname = markname[:-3]
    else:
        ending = ''
    data, thresh, msw, markname_full = _calc_threshold(markname,adata,data_key,n,include_zeroes,preset_threshold=preset_threshold, force_model=force_model)
    if plot:
        _plot_threshold(data, thresh, msw, markname_full+ending,title_addition=title_addition)
    if fancy:
        return _fancy_interactive_threshold(thresh, max(data))
    if interactive:
        if not plot:
            print(markname_full)
        thresh = _interactive_threshold(thresh)
    if run:
        thresholded = _run_threshold(data, thresh)
        return thresh, thresholded
    else:
        return thresh
    
def _calc_threshold(markname: str, adata: anndata.AnnData,
                    data_key: str=utils.DATA_KEY, n: int=0,
                    include_zeroes: bool=False, preset_threshold: Tuple[Union[int,float]]=None,
                    force_model: bool=False) -> Tuple[anndata.AnnData, Tuple[float,float], Tuple[Tuple[float, ...], Tuple[float,...], Tuple[float, ...]], str]:
    '''
        Internal function to handle threshold making.
        
        Gets the data using utils.get_data, determines if 0's should be kept in the data, then creates a gaussian mixture model of either 1 or 2 components. 
        The fit of these models are compared and if only use the 2 component mix if it fits significantly better than the 1 component model. Calculate means, 
        weights, and standard deviations of the gaussians, then calculates thresholds as the minimum and maximum of the 1st peak + 2-4 std, 
        the 2nd peak + 2-4 std (if there're two peaks), and .05 (if there are no peaks)
        
        Parameters
        ----------
        markname
            name of the marker
        adata
        AnnData object for creating thresholds, containing expression in the .X and/or the .obsm[data_key] 
        data_key
           Name of the key in .obsm[] to look for when searching for markname
        n
            Determines the number of gaussians to fit. Can be 1, 2, or 3.
        include_zeroes
            Whether to include zeroes when calculating the gaussian mixture
        preset_threshold
            Max and minimum thresholds to apply for positive and negative populations, overrides calculation
        force_model
            Whether to force the creation of a gaussian mixture model. If preset thresholds are defined, this model is only used to create msw
        
        Returns
        -------
        data: AnnData object
            Data used to find threshold
        thresh: tuple of numerics
            Population differentiating thresholds, length is equal to n
        msw: tuple of 
            m: tuple of floats
                The means of peaks in the gaussian mixture model
            s: tuple of floats
                The standard deviation of these peaks, 
            w: tuple of floats
                The weights (heights) of these peaks
        markname: str
            Full name of the marker being used, found from utils.get_data()
            
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