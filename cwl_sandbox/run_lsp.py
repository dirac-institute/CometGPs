
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
from GP_class import GPFit

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

def main():

    # read in the data file
    time, flux, flux_err= read_data(filename, datadir)

    asteroid = GPFit(time, flux, flux_err)

    asteroid.run_lsp(filename, true_period, nterms)

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
    parser.add_argument('-n', '--nterms', action="store", dest="nterms", required=False, type=int, default=1,
                        help='The number of harmonics to plot apart from the true period.')

    clargs = parser.parse_args()

    filename = clargs.filename
    datadir = clargs.datadir
    true_period = clargs.period
    nterms = clargs.nterms

    main()
