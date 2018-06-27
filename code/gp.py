import numpy as np
import scipy.optimize as op

import george
from george import kernels

import dynesty




## SET UP THE GP
def setup_periodic_gp(time, flux, flux_err, change_term=False):
    """
    
    Parameters
    ----------
    time : 
    
    flux : 
    
    flux_err : 
    
    change_term : bool, default False
        If False, use only an ExpSine2Kernel (simple squared sine); if True,
        include a ExpSquared kernel (to model longer-term variations)
    
    Returns
    -------
    gp : george.GP object
        The GaussianProcess object
    """
    # starting parameters
    gamma = 1.0
    log_period = np.log(5.8/24.0)
    amp = 3.0

    amp2 = 1.0
    metric = 1.0

    kk = amp * kernels.ExpSine2Kernel(gamma=gamma, log_period=log_period)
    if change_term:
        k2 = amp2 * kernels.ExpSquaredKernel(metric=metric) 
        kk *= k2
    
    gp = george.GP(kk, mean=np.mean(flux), fit_mean=True,
               white_noise=np.log(np.mean(flux_err)))
    gp.compute(time)

    return gp

# sample the GP
def sample_gp(time, flux, flux_err, change_term=False, pred_points=1000,
              nlive=500, bound="multi", sample="rwalk"):

    # set up the GP
    gp = setup_periodic_gp(time, flux, flux_err, change_term=change_term)

    t0 = np.min(time)
    tmax = np.max(time)
    
    t_pred = np.linspace(t0, tmax, pred_points)

    ndim = len(gp.get_parameter_vector())

    # Define the objective function (negative log-likelihood in this case).
    def nll(p):
        gp.set_parameter_vector(p)
        ll = gp.log_likelihood(flux, quiet=True)
        return -ll if np.isfinite(ll) else 1e25

    # And the gradient of the objective function.
    def grad_nll(p):
        gp.set_parameter_vector(p)
        return -gp.grad_log_likelihood(flux, quiet=True)

    # You need to compute the GP once before starting the optimization.
    gp.compute(time)

    # Run the optimization routine.
    p0 = gp.get_parameter_vector()
    results = op.minimize(nll, p0, jac=grad_nll, method="L-BFGS-B")


    # Update the kernel and print the final log-likelihood.
    gp.set_parameter_vector(results.x)

    # set up the sampling from the prior
    def from_prior(u):

        # prior transform for the mean value
        min_mag = 12
        max_mag = 20
        pmean = (max_mag - min_mag) * u[0] + min_mag

        # prior transform for the GP log-constant
        min_logamp = -20
        max_logamp = 20
        plogamp = (max_logamp - min_logamp) * u[1] + min_logamp

        # prior transform for the kernel gamma
        min_loggamma = np.log(0.001)
        max_loggamma = np.log(100.0)
        ploggamma = (max_loggamma - min_loggamma) * u[2] + min_loggamma

        # prior transform for the log-period
        min_period = 1./24.0
        max_period = 1.0

        pperiod = (max_period - min_period) * u[3] + min_period

        if change_term:
            min_logamp2 = -20
            max_logamp2 = 20
            plogamp2 = (max_logamp2 - min_logamp2) * u[4] + min_logamp2

            min_logmetric = -10
            max_logmetric = 10
            plogmetric = (max_logmetric - min_logmetric) * u[5] + min_logmetric
    
            return np.array([pmean, plogamp, ploggamma, pperiod, plogamp2, plogmetric])
        else:
            return np.array([pmean, plogamp, ploggamma, pperiod])
        
    # set up the log-lilkelihood
    def loglike(p):
        # this is if the change-term is included:
        if change_term:
            pnew = [p[0], p[1], np.exp(p[2]), np.log(p[3]), p[4], p[5]]
        else:
            pnew = [p[0], p[1], np.exp(p[2]), np.log(p[3])]

        gp.set_parameter_vector(pnew)
        ll = gp.log_likelihood(flux, quiet=True)
        return ll if np.isfinite(ll) else -1e25
    
    # set up the DynamicNestedSampling sampler
    sampler = dynesty.DynamicNestedSampler(loglike, from_prior, ndim, nlive=nlive, 
                                    bound=bound, sample=sample)

    sampler.run_nested(nlive_init=1000)
    
    return sampler
