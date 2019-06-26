
import os
os.environ["MKL_NUM_THREADS"] = "3"

import argparse
import textwrap

import numpy as np
import pandas as pd
import h5py

from plotting import plot_mcmc_sampling_results
from GP_class import GPFit

def read_data(filename, whitespace=False, datadir="./"):
    """
    Read in light curve data from asteroid.
    """

    data  = pd.read_csv(datadir+filename, header=None, delim_whitespace=whitespace)

    tsample = data[0]
    fsample = data[1]
    flux_err = data[2]

    return tsample, fsample, flux_err

def write_data(filename, sampler, asteroid, nwalkers, niter):
    """
    Write the sampler results as an HDF5 file,
    with all the other info you might want.
    """

    with h5py.File(filename+".hdf5", "w") as f:
        f.create_dataset("chain", data=sampler.chain)

        if asteroid.true_period == None:
            f.attrs['true_period'] = 0
        else:
            f.attrs['true_period'] = asteroid.true_period

        f.attrs['walkers'] = nwalkers
        f.attrs['iterations'] = niter
        f.attrs['data_pts'] = asteroid.data_pts
        f.attrs['acceptance_faction'] = sampler.acceptance_fraction
        f.create_dataset("time", data= asteroid.time)
        f.create_dataset("flux", data = asteroid.flux)
        f.create_dataset("flux_err", data = asteroid.flux_err)

def main():
    # read in the data file
    time, flux, flux_err= read_data(filename, whitespace, datadir)

    asteroid = GPFit(time, flux, flux_err)
    asteroid.set_params()
    asteroid.set_walker_param_matrix(nwalkers)
    asteroid.set_gp_kernel()

    sampler = asteroid.run_emcee(niter=niter, nwalkers=nwalkers, burn_in=burn_in, threads=threads)

    write_data(filename, sampler, asteroid, nwalkers, niter)

    #plot_mcmc_sampling_results(np.array(asteroid.time), asteroid.flux, asteroid.flux_err,
#                               asteroid.gp, sampler, namestr=filename + "_plots",
#                               true_period=true_period)

    #if lsp:
    #    asteroid.run_lsp(filename, true_period, nterms)



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
    parser.add_argument('-p', '--period', action="store", dest="period", required=False, type=float, default=None,
                        help="The true period of an asteroid in hours.")
    parser.add_argument('-ws', '--whitespace', action="store_true", dest="whitespace", required=False, default=False,
                        help="The delimeter for the input file, assumed to be whitespace.")
    parser.add_argument('-n', '--nterms', action="store", dest="nterms", required=False, type=int, default=1,
                        help='The number of harmonics to plot apart from the true period.')
    parser.add_argument('-lsp', '--lsp', action="store_true", dest="lsp", required=False, default=True,
                        help="Generates a Lomb-Scargle periodogram.")
    parser.add_argument('-b', '--burnin', action="store", dest="burn_in", required=False, type=int, default=2000,
                        help="The number of iterations to remove from the head of the MCMC chain walkers.")

    clargs = parser.parse_args()

    filename = clargs.filename
    datadir = clargs.datadir
    nwalkers = clargs.nwalkers
    niter = clargs.niter
    threads = clargs.threads
    true_period = clargs.period
    whitespace = clargs.whitespace
    nterms = clargs.nterms
    lsp = clargs.lsp
    burn_in = clargs.burn_in

    main()
