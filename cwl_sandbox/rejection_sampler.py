import argparse
import textwrap

import numpy as np
import pandas as pd
import scipy.stats
import george
import numpy.random as rnd

def read_data(filename, whitespace=False, datadir="./"):
    """
    Read in light curve data from asteroid.
    """

    data  = pd.read_csv(datadir+filename, header=None, delim_whitespace=whitespace)

    tsample = data[0]
    fsample = data[1]
    flux_err = data[2]

    return tsample, fsample, flux_err



def main():
    # read in the data file
    time, flux, flux_err= read_data(filename, whitespace, datadir)

    # 1 : set up the prior distriubtion
    prior = scipy.stats.norm(np.log(4./24.), (12./24.))

    # sample the prior pdf with nsample rvs
    J = prior.rvs(nsamples)

    # set up initial parameter values
    mean_flux = np.mean(flux)
    log_amp = np.log(flux.max()-flux.min())
    gamma = 1
    log_period = 0

    # set up the kernel parameters
    kernel = np.exp(log_amp) * george.kernels.ExpSine2Kernel(gamma = gamma, log_period = log_period)
    gp = george.GP(kernel, fit_mean=True, mean=mean_flux)
    gp.compute(time, flux_err)

    # 2: for each nsample, calculate the log likelihood
    L_results = np.ones(nsamples)

    for i in np.arange(nsamples):
        # input the new sampled period value into the kernel
        params = [mean_flux, log_amp, gamma, J[i]]
        gp.set_parameter_vector(params)

        # calculate the log likelihood
        # and check to see if it is a valid number
        try:
            gp.compute(time, flux_err)
            lnlike = gp.log_likelihood(flux)
        except np.linalg.LinAlgError:
            lnlike = -1e25

        # add log likelihood value to results array
        L_results[i] = lnlike

    # 3 : Pick a random number r out of a uniform distribution between 0 and 1
    #     and compare the normalized (between 0 and 1) likelihood value to it
    uu = rnd.uniform(size=len(L_results))

    good_samples_bool = uu < np.exp(L_results-L_results.max())
    good_samples_idx, = np.where(good_samples_bool)

    print(np.exp(J[good_samples_idx])*24.)

    ### figure out what you want to return

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

            $> python rejection_sampler.py --help

    Run this script from anywhere on your system:

            $> python /absolute/path/to/CometGP/cwl_sandbox/rejection_sampler.py --help


    Run on example data in the data directory:

            $> python /absolute/path/to/CometGP/cwl_sandbox/rejection_sampler.py -f "2001SC170.csv"
                    -d "absolute/path/to/CometGP/data/asteroid_csv"

    Run on example data (from example data directory) with more walkers, steps, etc.

            $> python ../code/run_gp.py -f "2001SC170.csv" -d "./" -j 2**24 -i 5000


    """))

    ### other arguments
    parser.add_argument('-f', '--filename', action="store", dest="filename", required=True,
                        help="Data file with observed time (in unit days) and flux.")
    parser.add_argument('-d', '--datadir', action="store", dest="datadir", required=False, default="./",
                        help="Directory with the data (default: current directory).")
    parser.add_argument('-j', '--nsamples', action="store", dest="nsamples", required=False, type=int, default=10**4,
                        help="The number of j samples you want to retreive from the prior pdf.")
    parser.add_argument('-ws', '--whitespace', action="store_true", dest="whitespace", required=False, default=False,
                        help="The delimeter for the input file, assumed to be whitespace.")


    clargs = parser.parse_args()

    filename = clargs.filename
    datadir = clargs.datadir
    nsamples = clargs.nsamples
    whitespace = clargs.whitespace


    main()
