
import numpy as np
import numpy.linalg
import scipy.stats
from scipy.signal import argrelextrema
import pandas as pd

import argparse
import textwrap

import emcee
import george

from emcee_utils import walker_params
from plotting import plot_lightcurve, plot_folded_lightcurve, plot_mcmc_sampling_results, plot_steps

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

    p_mean = scipy.stats.uniform(0,20).logpdf(params[0])
    p_log_amp = scipy.stats.uniform(-10,30).logpdf(params[1])
    p_log_gamma = scipy.stats.uniform(np.log(0.1), (np.log(40)-np.log(0.1))).logpdf(np.log(params[2]))
    p_period = scipy.stats.uniform(np.log(0.5/24), -np.log(0.5/24)).logpdf((params[3]))

    sum_log_prior =  p_mean + p_log_amp + p_log_gamma + p_period

    if np.isnan(sum_log_prior) == True:
        return -np.inf

    return sum_log_prior


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

    # compute lnlikelihood based on given parameters
    gp.set_parameter_vector(params)

    try:
        gp.compute(tsample, flux_err)
        ln_likelihood = gp.lnlikelihood(fsample)+log_prior
    except np.linalg.LinAlgError:
        ln_likelihood = -1e25

    return ln_likelihood if np.isfinite(ln_likelihood) else -1e25


def read_data(filename, datadir="./"):
    """
    Read in light curve data from asteroid.
    """

    data  = pd.read_csv(datadir+filename, header=None)

    tsample = data[0]
    fsample = data[1]
    flux_err = data[2]
    data_pts = len(tsample)

    return tsample, fsample, flux_err, data_pts

def run_gp(filename, datadir="./", nchain=100, niter=100, gamma=1, cov_scale=1, threads=1):
    tsample, fsample, flux_err, data_pts = read_data(filename, datadir)

    #get l-s best period estimate
    from lombscargle import make_lsp
    from astropy.stats import LombScargle

    freq, power = make_lsp(tsample, fsample, flux_err, p_max=5.0)

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



    ndim, nwalkers = 4, nchain

    # initialize walker parameters
    best_log_amp = np.log(fsample.max()-fsample.min())

    params = [np.mean(fsample), best_log_amp, gamma, new_log_period]
    p0, gp = walker_params(params, fsample, flux_err, nwalkers, cov_scale=cov_scale)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, post_lnlikelihood, args=[gp, tsample, fsample, flux_err], threads=threads)
    mcmc_sampling = sampler.run_mcmc(p0, niter)

    def save_chain(file_name, sampler):
        header = str(sampler.chain.shape)
        np.savetxt(file_name, sampler.flatchain, header=header)
        return

    save_chain(filename + "_results", sampler)

    ###AVOID HAVING TO USE pd SO YOU DON'T HAVE TO CHANGE THIS TO A NP.ARRAY
    tsample = np.array(tsample)

    plot_mcmc_sampling_results(tsample, fsample, flux_err, gp, sampler, namestr=filename + "_plots")


    return



def main():
    run_gp(filename, datadir, nchain, niter, gamma, cov_scale, threads)

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

            $> python ../code/run_gp.py -f "2001SC170.csv" -d "./" -c 50 -i 5000 -g 5 -c 0.5 -t 2


    """))

    ### other arguments
    parser.add_argument('-f', '--filename', action="store", dest="filename", required=True,
                        help="Data file with observed time (in unit days) and flux.")
    parser.add_argument('-d', '--datadir', action="store", dest="datadir", required=False, default="./",
                        help="Directory with the data (default: current directory).")
    parser.add_argument('-c', '--nchain', action="store", dest="nchain", required=False, type=int, default=100,
                        help="The number of walkers/chains for the MCMC run (default: 100).")
    parser.add_argument('-i', '--niter', action="store", dest="niter", required=False, type=int, default=100,
                        help="The number of iterations per chain/walker in the MCC run (default: 100).")
    parser.add_argument('-g', '--gamma', action="store", dest="gamma", required=False, type=int, default=1,
                        help="The length scale of variations within a single period (default: 1).")
    parser.add_argument('-s', '--cov_scale', action="store", dest="cov_scale", required=False, type=int, default=1,
                        help="Determines the scatter of the multivariate distribution (default: 1).")
    parser.add_argument('-t', '--threads', action="store", dest="threads", required=False, type=int, default=1,
                        help="The numer of threads used for computing the posterior (default: 1).")

    clargs = parser.parse_args()

    filename = clargs.filename
    datadir = clargs.datadir
    nchain = clargs.nchain
    niter = clargs.niter
    gamma = clargs.gamma
    cov_scale = clargs.cov_scale
    threads = clargs.threads

    main()
