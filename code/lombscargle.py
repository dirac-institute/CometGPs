import numpy as np
from astropy.stats import LombScargle
from gatspy import periodic




def make_lsp(time, flux, flux_err, p_min=1/24.0, p_max=1.0, oversampling=5, lsp_kwargs={}, nterms=1, use_gatspy=False):
    """
    Compute the Lomb-Scargle periodogram using the astropy LombScargle class.
    
    Parameters
    ----------
    time : numpy.ndarray
        The time stamps of the time series **IN DAYS**
        
    flux : numpy.ndarray
        Flux measurements corresponding to the time stamps
         
    flux_err : numpy.ndarray
        The flux uncertainties corresponding to the data. 

    p_min, p_max : float, float
        The minimum and maximum period to search, **IN DAYS**
        
    oversampling : int
        The oversampling factor; default is 5
        
    lsp_kwargs : dict
        Optional keyword arguments for astropy.stats.LombScargle
        
    Returns
    -------
    freq : numpy.ndarray
        The frequencies for which the Lomb-Scargle periodogram 
        was computed
        
    power : numpy.ndarray
        Normalized Lomb-Scargle power
        
    """
    # number of frequencies to sample
    nfreq = int(oversampling*(p_max - p_min)/p_min) # number of frequencies

    # make a list of frequencies 
    freq = np.linspace(1/p_max, 1/p_min, nfreq)
    
    if use_gatspy:
        model = periodic.LombScargle(fit_period=True, center_data=True, fit_offset=True, Nterms=nterms)
        model.optimizer.period_range = (p_min, p_max)
        
        model.fit(time, flux, flux_err)
        periods = 1/freq
        power = model.periodogram(periods)
        power = power
        
    else:
        # initialize the lomb-scargle periodogram
        ls = LombScargle(time, flux, dy=flux_err, **lsp_kwargs, nterms=nterms)

        # compute the power at each given frequency
        power = ls.power(freq)
    
    return freq, power

