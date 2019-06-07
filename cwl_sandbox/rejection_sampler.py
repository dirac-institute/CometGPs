import argparse
import textwrap

import numpy as np
import pandas as pd
import scipy.stats
import george
import numpy.random as rnd
import multiprocessing as mp
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import time as tm
import corner
import datetime


def read_data(filename, whitespace=False, datadir="./"):
    """
    Read in light curve data from asteroid.
    """

    data  = pd.read_csv(datadir+filename, header=None, delim_whitespace=whitespace)

    tsample = data[0]
    fsample = data[1]
    flux_err = data[2]

    return tsample, fsample, flux_err


class DataManager:
    def __init__(self, flux, time, flux_err):
        self.flux = flux
        kernel = np.exp(0) * george.kernels.ExpSine2Kernel(gamma = 1, log_period = 0)
        self.gp = george.GP(kernel, fit_mean=True, mean=1)
        self.gp.compute(time, flux_err)

    def calculate_likelihood(self, params):
        self.gp.set_parameter_vector(params)

        try:
            #gp.compute(time, flux_err)
            lnlike = self.gp.log_likelihood(self.flux)
        except np.linalg.LinAlgError:
            lnlike = -1e25

        return lnlike


def main():
    # read in the data file
    time, flux, flux_err= read_data(filename, whitespace, datadir)

    # 1 : set up the prior distriubtion
    prior_mean = scipy.stats.norm(1, 0.5)
    prior_log_amp = scipy.stats.norm(np.log(0.15), np.log(2))
    prior_log_gamma = scipy.stats.norm(np.log(10), np.log(2))
    prior_log_period = scipy.stats.norm(np.log(4./24.), (12./24.))

    # sample the prior pdf with nsample rvs
    J_mean = prior_mean.rvs(nsamples)
    J_log_amp = prior_log_amp.rvs(nsamples)
    J_log_gamma = prior_log_gamma.rvs(nsamples)
    J_log_period= prior_log_period.rvs(nsamples)

    # set up initial parameter values
    mean_flux = np.mean(flux)
    log_amp = np.log(flux.max()-flux.min())
    gamma = 1
    log_period = 0

    manager = DataManager(flux, time, flux_err)
    # 2: for each nsample, calculate the log likelihood

    start_time = tm.time()


    pool = mp.Pool(processors)

    L_results = []
    # L_results = pool.map(manager.calculate_likelihood, [i for i in np.arange(nsamples)])
    L_results = pool.map(manager.calculate_likelihood, zip(J_mean, J_log_amp, np.exp(J_log_gamma), J_log_period))

    print("Pool has finished.")

    L_results = np.array(L_results)

    print("Results have been converted into an array!")

    pool.close()

    end_time = tm.time()
    print("\nTotal time taken : %.2f seconds" %(end_time - start_time))

    # 3 : Pick a random number r out of a uniform distribution between 0 and 1
    #     and compare the normalized (between 0 and 1) likelihood value to it




    good_samples_bool = rnd.uniform(size=len(L_results)) < np.exp(L_results-L_results.max())
    #good_samples_idx, = np.where(good_samples_bool)

    #summing the bool = passes
    print("Number of samples passed : %s" %sum(good_samples_bool))

    # print(np.exp(J[good_samples_idx])*24.)

    #plt.hist(np.exp(J_log_period[good_samples_idx])*24)


    ### figure out what you want to return

    now = datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm")

    data = np.array([np.exp(J_log_period[good_samples_bool])*24, np.exp(J_log_gamma[good_samples_bool]), np.exp(J_log_amp[good_samples_bool]), J_mean[good_samples_bool], L_results[good_samples_bool]]).T
    print(data)

    np.savetxt(filename + "results_%s.txt" %(now), data)

    figure = corner.corner(data, labels=["period", 'gamma', 'amp', 'mean'])#, J_log_gamma[good_samples_idx], J_mean[good_samples_idx])
    plt.savefig(filename + "_rej_sampler_posterior_%s.pdf" %now, format="pdf")


    f = open(filename + "_run_reports_%s.txt" %now, "w+")

    f.write(str(now) + "\n" +
            filename + "\n"
            "nsamples: %d \n" % nsamples +

            #filename "\n"
            #"nsamples %d \n" % nsamples
            "processors: %d \n" % processors +
            "time: %.2f sec \n" % float(end_time - start_time) +
            "uu: %.2f \n" %uu +
            "passed: %d \n\n" %sum(good_samples_bool) +
            str(data))
    f.close()


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
    parser.add_argument('-p', '--processors', action="store", dest="processors", required=False, type=int, default=1,
                        help="The number of processors you want to use for multiprocessing.")
    parser.add_argument('-u', '--uu', action="store", dest="uu", required=False, type=float, default=0.1,
                        help="The limit for calculating what constitutes a good likelihood.")

    clargs = parser.parse_args()

    filename = clargs.filename
    datadir = clargs.datadir
    nsamples = clargs.nsamples
    whitespace = clargs.whitespace
    processors = clargs.processors
    uu = clargs.uu


    main()
