import matplotlib.pyplot as plt
import numpy as np
import george
import scipy.stats

def walker_params(params, fsample, flux_err, nwalkers, cov_scale=1):

    """
    Sets up the initial parameters for all the walkers using optimized
    parameter values as starting values. The function generates a
    scattered multivatiate gaussian distribution of starting parameter
    values.

    Parameters
    ----------
    params : list
        List of all kernel parameters.

    cov_scale : float
        Determines the scatter of the multivariate distribution.

    Returns
    -------
    p0 : numpy.ndarray
        The initial walker parameters [nwalker, ndim]

    gp : george.gp.GP
        GP kernel set with the optimized parameter values.

    """


    gp_mean, log_amp, gamma, log_period = params
    amp = np.exp(log_amp)

    print('amp : ' + str(amp))
    kernel = amp * george.kernels.ExpSine2Kernel(gamma = gamma, log_period = log_period)
    gp = george.GP(kernel, fit_mean=True, mean=gp_mean)
    gp.compute(fsample, flux_err)

    p_start = np.array(params)/100.
    cov_matrix = np.sqrt(np.diag(p_start)**2)*cov_scale
    print("params : " + str(params))
    print("cov matrix : \n" + str(cov_matrix))

    p0 = np.random.multivariate_normal(mean=params, cov=cov_matrix, size=nwalkers)

    return p0, gp

def plot_gpfit(time, fsample, flux_err, gp, ax):

    """
    Plot a gp fit given a gp class and x, y, and yerr data to fit onto.

    """

    t_possible = np.linspace(time[0], time[-1], len(time))
    pred, pred_var = gp.predict(fsample, t_possible, return_var=True)

    temp_color = np.random.rand(3)

    ax.fill_between(t_possible, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var),
                color="red", alpha=0.4)
    ax.plot(t_possible, pred, "red", lw=1.5, alpha=0.7, label = "GP Fit : " + str(round(gp.parameter_vector[-1], 5)))
    ax.legend()
