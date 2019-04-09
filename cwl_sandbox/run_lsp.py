
#import os
#os.environ["MKL_NUM_THREADS"] = "3"

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import george
import emcee
import scipy.stats
import pandas as pd


import plotting

import argparse
import textwrap

def read_data(filename, datadir="./"):
    """
    Read in light curve data from asteroid.
    """

    data  = pd.read_csv(datadir+filename, header=None, delim_whitespace=False)

    tsample = data[0]
    fsample = data[1]
    flux_err = data[2]

    return tsample, fsample, flux_err


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
            p0[:,3] = np.log(np.linspace(2,12,nwalkers)/24.)

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

    def run_emcee(self, nwalkers, niter, threads=1):
        """Runs emcee's mcmc code."""

        ndim = 4
        sampler = emcee.EnsembleSampler(nwalkers, ndim, post_lnlikelihood, args=[self.gp, self.time, self.flux, self.flux_err], threads=threads)

        mcmc_sampling = sampler.run_mcmc(self.walker_params, niter)
        self.sampler = sampler

        return sampler

    def run_lsp(self, true_period, nterms):
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

        # plot all the frequencies
        fig, (ax, bx) = plt.subplots(1,2, figsize=(12,5))
        fig.set_tight_layout('tight')
        ax.plot((1./freq)*24.,power)
        ax.set_xlabel('Period (hrs)')
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

    def daniela_lsp(self, true_period, nterms):
        """Testing out Daniela's LSP method."""

        ax = plotting.plot_lsp(self.time, self.flux, self.flux_err, p_max=5.0, true_period=true_period,
                                nterms=nterms)

        namestr=filename + "_plots_daniela"
        plt.savefig(namestr + "_lsp.pdf", format="pdf")

        return


def main():

    # read in the data file
    time, flux, flux_err= read_data(filename, datadir)

    asteroid = GPFit(time, flux, flux_err)
    asteroid.daniela_lsp(true_period, nterms)

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
    parser.add_argument('-p', '--period', action="store", dest="period", required=False, type=float, default=None,
                        help="The true period of an asteroid in hours.")
    parser.add_argument('-n', '--nterms', action="store", dest="nterms", required=False, type=int, default=None,
                        help='The number of harmonics to plot apart from the true period.')

    clargs = parser.parse_args()

    filename = clargs.filename
    datadir = clargs.datadir
    true_period = clargs.period
    nterms = clargs.nterms

    main()
