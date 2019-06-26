import numpy as np
import matplotlib.pyplot as plt
import george
import emcee
import scipy.stats
import pandas as pd

import plotting

class GPFit():
    def __init__(self, time_stamps, flux, flux_error):
        self.time = time_stamps
        self.flux = flux
        self.flux_err = flux_error
        self.data_pts = len(time_stamps)
        self.true_period = None
        self.params = None
        self.walker_params = None
        self.gp = None
        self.sampler = None
        self.lsp_period = None

    def set_params(self):
        """Calculates initial gp parameter values based on data."""
        mean_flux = np.mean(self.flux)
        log_amp = np.log(self.flux.max()-self.flux.min())
        gamma = 1
        log_period = 0

        parameters = {"mean": mean_flux, "log_amp": log_amp, "gamma": gamma,"log_period": log_period}
        self.params = parameters
        return


    def set_walker_param_matrix(self, nwalkers):
        """Creates a matrix of starting parameters for every walker."""

        if self.params is not None:

            p_start = np.array(list(self.params.values()))
            cov_matrix = np.sqrt(np.diag(p_start)**2)
            p0 = np.random.multivariate_normal(mean=p_start, cov=cov_matrix, size=(nwalkers))

            # equally distributed starting period values
            p0[:,3] = np.random.normal(size=nwalkers)*0.5 + np.log(4/24.)

            self.walker_params = p0

        else:
            print("Please set parameter values first.")

        return

    def set_gp_kernel(self):
        """Sets up the Gaussian Process Kernel that is needed for george."""

        kernel = np.exp(self.params["log_amp"]) * george.kernels.ExpSine2Kernel(gamma = self.params["gamma"], log_period = self.params["log_period"])
        gp = george.GP(kernel, fit_mean=True, mean=self.params["mean"])
        gp.compute(self.time, self.flux_err)

        self.gp = gp

        return

    def run_emcee(self, nwalkers, niter, threads, burn_in):
        """Runs emcee's mcmc code."""

        ndim = 4
        sampler = emcee.EnsembleSampler(nwalkers, ndim, post_lnlikelihood, args=[self.gp, self.time, self.flux, self.flux_err], threads=threads)

        #run steps for a burn-in
        state = sampler.run_mcmc(self.walker_params, burn_in)
        sampler.reset()
        print(state[0])
        sampler.run_mcmc(state[0], niter)
        self.sampler = sampler

        return sampler

    def run_lsp(self, filename, true_period, nterms):
        """Determines the Lomb-Scargle Periodogram."""

        from scipy.signal import argrelextrema

        #get l-s best period estimate
        from lombscargle import make_lsp
        from astropy.stats import LombScargle

        freq, power = make_lsp(self.time, self.flux, self.flux_err, p_max=5.0, nterms=nterms)

        # determine the indices of local power maxima
        best_idx = argrelextrema(power, np.greater)

        # sort these indices based on actual power value
        # reverse list so max is read first
        indices = np.argsort(power[best_idx[0]])[::-1]

        # sort our original indices based on the new
        # power-sorted indices
        best_idx = (best_idx[0]).T[indices]
        best_freqs = freq[best_idx].T

        new_freq = best_freqs[0]
        new_period = 1./new_freq
        new_log_period = np.log(1./new_freq)

        self.true_period = true_period
        self.lsp_period = new_period*24.

        # plot all the frequencies
        fig, (ax, bx) = plt.subplots(1,2, figsize=(12,5))
        fig.set_tight_layout('tight')
        ax.plot((1./freq)*24., power, color="black", alpha=0.7)
        ax.set_xlabel('Period (hrs)')
        ax.set_ylabel("Normalized Power")
        ax.vlines(new_period*24., 0, 1, colors='orange', linestyles='--',
                  label = 'Best fit : ' + str(round(new_period*24., 5)))
        ax.vlines(true_period, 0, 1, colors='blue', linestyles='--',
                  label = 'True fit : ' + str(true_period))
        ax.set_xlim([0,24])
        ax.legend()

        bx = plotting.plot_folded_lightcurve(self.time, self.flux, period=new_period, ax=bx)

        namestr=filename + "_plots"
        plt.savefig(namestr + "_lsp.pdf", format="pdf")

        return


def prior(params):

    """
    Calculated the log of the prior values, given parameter values.

    Parameters
    ----------
    params : list
        List of all kernel parameters

    param[0] : float
        mean (between 0 and 2)

    param[1] : float
        log amplitude (between -10 and 10)

    param[2] : float
        gamma (log gamma between 0.1 and 40)

    param[3] : float
        log period (period between 1h and 24hrs)

    Returns
    -------
    sum_log_prior : int
        sum of all log priors (-inf if a parameter is out of range)

    """

    p_mean = scipy.stats.norm(1, 0.5).logpdf(params[0])
    p_log_amp = scipy.stats.norm(np.log(0.15), np.log(2)).logpdf(params[1])
    p_log_gamma = scipy.stats.norm(np.log(10), np.log(2)).logpdf(np.log(params[2]))
    #print(params[2])
    #print("   " + str(p_log_gamma))
    p_log_period = scipy.stats.norm(np.log(4./24.), (12./24.)).logpdf(params[3])
    # log period (period between 0.5hrs and 36hrs)
    #p_log_period = scipy.stats.uniform(np.log(0.5/24), -(np.log(2/3)+np.log(0.5/24))).logpdf((params[3]))

    sum_log_prior =  p_mean + p_log_amp + p_log_gamma + p_log_period

    if np.isnan(sum_log_prior) == True:
        return -np.inf

    return sum_log_prior


def logl(params, gp, tsample, fsample, flux_err):
     # compute lnlikelihood based on given parameters
     gp.set_parameter_vector(params)


     try:
         gp.compute(tsample, flux_err)
         lnlike = gp.lnlikelihood(fsample)
     except np.linalg.LinAlgError:
         lnlike = -1e25

     return lnlike


def post_lnlikelihood(params, gp, tsample, fsample, flux_err):

    """
    Calculates the posterior likelihood from the log prior and
    log likelihood.

    Parameters
    ----------
    params : list
        List of all kernel parameters

    Returns
    -------
    ln_likelihood : float
        The posterior, unless the posterior is infinite, in which case,
        -1e25 will be returned instead.

    """

    # calculate the log_prior
    log_prior = prior(params)

    # return -inf if parameters are outside the priors
    if np.isneginf(log_prior) == True:
        return -np.inf

    try:
        lnlike = logl(params, gp, tsample, fsample, flux_err)
        ln_likelihood = lnlike+log_prior

    except np.linalg.linalg.LinAlgError:
        ln_likelihood = -1e25

    return ln_likelihood if np.isfinite(ln_likelihood) else -1e25
