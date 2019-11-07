# test new kernel addition
import argparse
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import george
import emcee
import scipy.stats
import corner
import h5py

import run_gp
import GP_class
import run_plotting

def prior(params):

    """
    Calculated the log of the prior values, given parameter values.

    Parameters
    ----------
    params : list
        List of all kernel parameters

    'mean:value',

    'kernel:k1:k1:log_constant',

    'kernel:k1:k2:metric:log_M_0_0',

    'kernel:k2:k1:log_constant',

    'kernel:k2:k2:gamma',

    'kernel:k2:k2:log_period')

    Returns
    -------
    sum_log_prior : int
        sum of all log priors (-inf if a parameter is out of range)

    """

    p_mean = scipy.stats.norm(17, 0.5).logpdf(params[0])
    p_log_amp_k1 = scipy.stats.norm(np.log(2), np.log(10)).logpdf(params[1])
    p_log_metric = scipy.stats.norm(np.log(100), np.log(10)).logpdf(np.log(params[2]))

    p_log_amp_k2 = scipy.stats.norm(np.log(2), np.log(2)).logpdf(params[3])
    p_log_gamma = scipy.stats.norm(np.log(10), np.log(2)).logpdf(np.log(params[4]))
    p_log_period = scipy.stats.norm(np.log(4./24.), (12./24.)).logpdf(params[5])


    sum_log_prior =  p_mean + p_log_amp_k1 + p_log_metric + p_log_amp_k2 + p_log_gamma + p_log_period

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

def main():

    #read in the data
    obs_files = filename#list()

    #import glob, os

    #for file in glob.glob("../data/paper_plots/ztf_lightcurves/*k2.txt"):
     #   obs_files.append(file)
        
    #print(obs_files)
        
    for i in np.arange(len(obs_files)):
        print(obs_files)

        time, flux, flux_err = run_gp.read_data(obs_files, whitespace=False)

        # Calculates initial gp parameter values based on data

        mean_flux = np.mean(flux)

        # k1
        log_amp_k1 = np.log(flux.max()-flux.min())
        metric = 5**2

        # k2
        log_amp_k2 = np.log(0.5)
        gamma = 10
        log_period = np.log(6/24.)

        parameters = {"mean": mean_flux, "log_amp_k1": log_amp_k1, "metric": metric, "log_amp_k2": log_amp_k2, "gamma": gamma,"log_period": log_period}
        params = parameters

        # Creates a matrix of starting parameters for every walker.
        p_start = np.array(list(params.values()))
        cov_matrix = np.sqrt(np.diag(p_start)**2)
        p0 = np.random.multivariate_normal(mean=p_start, cov=cov_matrix, size=(nwalkers))

        # equally distributed starting period values for
        p0[:,-1] = np.random.normal(size=nwalkers)*0.5 + np.log(4/24.)

        walker_params = p0

        """Calculates initial gp parameter values based on data."""
        k1 = np.exp(log_amp_k1) * george.kernels.ExpSquaredKernel(metric=metric)
        k2 = np.exp(log_amp_k2) * george.kernels.ExpSine2Kernel(gamma=gamma, log_period=log_period)

        kernel = k1*k2

        gp = george.GP(kernel, mean=np.mean(flux), fit_mean=True)

        gp.compute(time, flux_err)
                       
        print("GP kernel is set.")

        ndim = 6
        threads = 10
        iterations = niter
        burn_in=10000
                     
        print("Burn-in starting.")

        sampler = emcee.EnsembleSampler(nwalkers, ndim, post_lnlikelihood, args=[gp, time, flux, flux_err])

        #run steps for a burn-in
        state = sampler.run_mcmc(walker_params, burn_in)
        sampler.reset()
                       
        print("Burn-in complete.")
        #print(state[0])
        data = sampler.run_mcmc(state[0], iterations)
                       
        
        print("Saving data.")

        with h5py.File(obs_files + "profile_testing.hdf5", "w") as f:
            f.create_dataset("chain", data=sampler.chain)

            f.attrs['walkers'] = nwalkers
            f.attrs['iterations'] = niter
            f.attrs['data_pts'] = len(flux)
            f.attrs['acceptance_fraction'] = sampler.acceptance_fraction
            #f.attrs['burn_in'] = burn_in
            f.create_dataset("time", data= time)
            f.create_dataset("flux", data = flux)
            f.create_dataset("flux_err", data = flux_err)

    return



if __name__ == "__main__":
    ### DEFINE PARSER FOR COMMAND LINE ARGUMENTS
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=" ", #Bayesian QPO searches for burst light curves.",
                                     epilog=textwrap.dedent("""
    NOTE! The first 3 columns of your input file "-f" must correspond to your
    time, flux, and flux error in that order. All columns beyond column 3 will be ignored.

    Examples
    --------

    Print this help message:

            $> python run_gp.py --help

    Run this script from anywhere on your system:

            $> python /absolute/path/to/CometGP/cwl_sandbox/run_gp.py --help


    Run on example data in the data directory:

            $> python /absolute/path/to/CometGP/cwl_sandbox/run_gp.py -f "2001SC170.csv"
                    -d "absolute/path/to/CometGP/data/asteroid_csv"

    Run on example data (from example data directory) with more walkers, steps, etc.

            $> python ../code/run_gp.py -f "2001SC170.csv" -d "./" -c 50 -i 5000 -c 0.5 -t 2


    """))

    ### other arguments
    parser.add_argument('-f', '--filename', action="store", dest="filename", required=True,
                        help="Data file with observed time (in unit days) and flux.")
    parser.add_argument('-d', '--datadir', action="store", dest="datadir", required=False, default="./",
                        help="Directory with the data (default: current directory).")
    parser.add_argument('-w', '--nwalkers', action="store", dest="nwalkers", required=False, type=int, default=100,
                        help="The number of walkers/chains for the MCMC run (default: 100).")
    parser.add_argument('-i', '--niter', action="store", dest="niter", required=False, type=int, default=1000,
                        help="The number of iterations per chain/walker in the MCMC run (default: 1000).")
    parser.add_argument('-t', '--threads', action="store", dest="threads", required=False, type=int, default=1,
                        help="The numer of threads used for computing the posterior (default: 1).")
    parser.add_argument('-b', '--burn_in', action="store", dest="burn_in", required=False, type=int, default=2000,
                        help="The number of iterations to remove from the head of the MCMC chain walkers.")

    clargs = parser.parse_args()

    filename = clargs.filename
    datadir = clargs.datadir
    nwalkers = clargs.nwalkers
    niter = clargs.niter
    threads = clargs.threads
    burn_in = clargs.burn_in

    main()
