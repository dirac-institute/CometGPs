

def set_kernel(time, flux, flux_err):
    """
    Set up the Gaussian Process kernel from original light curve data,
    because it couldn't be saved to an hdf5 file.

    Parameters
    ----------
    time : numpy.ndarray
        The time stamps of the periodic light curve

    flux : numpy.ndarray
        Flux measurements corresponding to the time stamps

    flux_err : numpy.ndarray
        The flux uncertainties corresponding to the data.

    Returns
    -------

    gp : GP object
        The basic Gaussian Process object.

    """

    # calculate hyperparameter estimates from light curve data
    mean_flux = np.mean(flux)
    log_amp = np.log(flux.max()-flux.min())
    gamma = 1
    log_period = 0

    params = {"mean": mean_flux, "log_amp": log_amp, "gamma": gamma,"log_period": log_period}

    # set up a non-stationary kernel using the estimated param values
    kernel = np.exp(params["log_amp"]) * george.kernels.ExpSine2Kernel(gamma = params["gamma"], log_period = params["log_period"])

    # assign the kernel to a GP object and fitting a mean model
    gp = george.GP(kernel, fit_mean=True, mean=params["mean"])

    # compute the covariance matrix and factorize it for a set of times and uncertainties
    gp.compute(time, flux_err)

    return gp

def plot_corner(data, gp, true_period=None, colours=None, zoom=False):
    """
    Plot a corner plot showing the projections of a data set in multi-dimesional space,
    with the different dimensions corresponding to the different kernel parameters.

    Parameters
    ----------
    data : numpy.ndarray
        Results pulled from hdf5 file. Assumes the shape to be [nwalkers, iterations, parameters].

    gp :  GP object
        The basic Gaussian Process object.

    true_period : float
        The true period of the asteroid light curves.

    colours : [str, str, str]
        List of (up to) three colours. First colour is used for the data, the second
        colour for the true underlying data, the third for the models.

    zoom : bool
        Toggles whether the corner plot will show a zoomed in version of the histogram,
        focusing on the densest region of the previous binned histogram.

    Returns
    -------

    figure : matplotlib.figure.Figure
        The object with the plot

    """

    if colours == None:
        colours = ["#000000", "#0072B2", "#E69F00", "#009E73", "#F0E442"]

    #get label names from the gp object
    labels = list(gp.get_parameter_names())
    labels[3] = 'period hours'

    flat_data = data.reshape(data.shape[0]*data.shape[1], data.shape[2])

    if zoom:
        prob, edges = calc_prob(data)
        flat_data = data[(data[:,:,3]>edges[0]) & (data[:,:,3]<edges[1])]

    figure = corner.corner(flat_data, labels=labels,
                           show_titles=True, title_kwargs={"fontsize": 8},
                           truths=[None, None, None, true_period],
                           truth_color=colours[1])

    return figure

def plot_corner_5_95(data, walkers, iterations, gp, ax=None):
    """
    Plot a corner plot showing the projections of the 5-95th percentile of a data set in multi-dimesional space,
    with the different dimensions corresponding to the different kernel parameters.

    Parameters
    ----------
    data : numpy.ndarray
        Results pulled from hdf5 file. Assumes the shape to be [nwalkers, iterations, parameters].

    walkers : numpy.int64
        The number of walkers/chains for the MCMC run.

    iterations : numpy.int64
        The number of iterations per chain/walker in the MCMC run.

    gp :  GP object
        The basic Gaussian Process object.

    Returns
    -------

    ax : matplotlib.Axes object
        The object with the plot

    """

    # mask the data
    lower, upper = np.percentile(data[:,:,3], [5,95])
    masked_data = data[(data[:,:,3]>lower) & (data[:,:,3]<upper)]

    #get label names from the gp object
    labels = list(gp.get_parameter_names())
    labels[3] = 'period hours'
    labels_no_outliers = [s + " 5-95th" for s in labels]

    ax = corner.corner(masked_data, labels=labels_no_outliers, , show_titles=True, title_kwargs={"fontsize": 8})
    return ax
