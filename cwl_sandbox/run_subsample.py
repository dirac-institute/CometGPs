import numpy as np
import pandas as pd

import argparse
import textwrap

from subsample import subsample #daniela's code


def run_subsample(filename, datadir="./", npoints=100, kind="random", days=1, delay=0):

    data = pd.read_csv(datadir+filename, delimiter=' ',
                     header=None, names=['time','flux'], dtype={'time':float, 'flux':float})

    days, delay = days, delay

    # convert days to points
    span = 2880 * days
    start_pt = 2880 * delay

    time = np.array(data.time[start_pt:span+start_pt])
    flux = np.array(data.flux[start_pt:span+start_pt])

    flux_err = np.ones_like(flux) * np.std(flux)/10.0
    tsample, fsample, flux_err = subsample(time, flux, flux_err=flux_err, npoints=npoints, kind=kind)

    ###NEED TO FIND A GOOD FORMAT FOR RETURNING THE SAMPLED DATA
    df = pd.DataFrame()

    df['time'] = tsample
    df['flux'] = fsample
    df['flux_err'] = flux_err


    np.savetxt(filename+"_sampled_"+str(kind)+"_"+str(days)+"days.txt", df, delimiter=",")

    return




def main():
    run_subsample(filename, datadir, npoints, kind, days, delay)

    return

if __name__ == "__main__":
    ### DEFINE PARSER FOR COMMAND LINE ARGUMENTS
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=" ", #Bayesian QPO searches for burst light curves.",
                                     epilog=textwrap.dedent("""
    Warning! Please make sure the first three columns of your input file correspond to the time, flux, and flux error accordingly.

    Examples
    --------

    Print this help message:

            $> python run_subsample.py --help

    Run this script from anywhere on your system:

            $> python /absolute/path/to/CometGP/cwl_sandbox/run_subsample.py --help


    Run on example data in the data directory:

            $> python /absolute/path/to/CometGP/cwl_sandbox/run_subsample.py -f "221_lc_49627_to_49787.txt"
                    -d "absolute/path/to/CometGP/data/asteroid_csv"

    Run on example data (from example data directory) with more walkers, more iterations, more simulations, just MORE!

            $> python ../code/run_subsample.py -f "221_lc_49627_to_49787.txt" -d "./" -n 50 -k "telescope" -days 3 -delay 2


    """))

    ### other arguments
    parser.add_argument('-f', '--filename', action="store", dest="filename", required=True,
                        help="Data file with observed time (in unit days) and flux.")
    parser.add_argument('-d', '--datadir', action="store", dest="datadir", required=False, default="./",
                        help="Directory with the data (default: current directory).")
    parser.add_argument('-n', '--npoints', action="store", dest="npoints", required=False, type=int, default=100,
                        help="Number of data points to subsample (default: 100).")
    parser.add_argument('-k', '--kind', action="store", dest="kind", required=False, type=str, default="random",
                        help="The type of subsampling to perform.")
    parser.add_argument('-days', '--days', action="store", dest="days", required=False, type=int, default=1,
                        help="The number of days over which to sample (defalut: 1).")
    parser.add_argument('-delay', '--delay', action="store", dest="delay", required=False, type=int, default=0,
                        help="The number of days to delay sampling from the start of the data set (default: 0).")


    clargs = parser.parse_args()

    filename = clargs.filename
    datadir = clargs.datadir
    npoints = clargs.npoints
    kind = clargs.kind
    days = clargs.days
    delay = clargs.delay

    main()
