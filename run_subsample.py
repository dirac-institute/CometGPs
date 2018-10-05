import argparse
import textwrap

def main():
    print("It works!")

    return

if __name__ == "__main__":
    ### DEFINE PARSER FOR COMMAND LINE ARGUMENTS
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=" ", #Bayesian QPO searches for burst light curves.",
                                     epilog=textwrap.dedent("""
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



asteroid = '1291'

true_log_p = {'3200':-1.896021, '1291':-1.45813055,
              '221':-0.8321219, '1388':-0.69789175}
true_p = {'3200':3.603957, '1291':5.58410,
              '221':10.443, '1388':11.9432}

txt = '../data/'+str(asteroid)+'_lc_49627_to_49787.txt'

data = pd.read_csv(txt, delimiter=' ',
                 header=None, names=['time','flux'], dtype={'time':float, 'flux':float})

days, delay = 5, 50

# convert days to points
span = 2880 * days
start_pt = 2880 * delay

time = np.array(data.time[start_pt:span+start_pt])
flux = np.array(data.flux[start_pt:span+start_pt])

flux_err = np.ones_like(flux) * np.std(flux)/10.0
tsample, fsample, flux_err = subsample(time, flux, flux_err=flux_err, npoints=100, kind="telescope")
