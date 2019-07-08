import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

import numpy as np
import copy

from dynesty import plotting as dyplot
from dynesty import utils as dyfunc

import corner

from lombscargle import make_lsp

def plot_lightcurve(time, flux, flux_err=None, true_lightcurve=None,
                    models=None, ax=None, colours=None):

    """
    Plot a light curve, potentially including the true underlying
    model that produced the data (in the case of simulations), or model
    light curves from MCMC.

    Parameters
    ----------
    time : numpy.ndarray
        The time stamps of the periodic light curve

    flux : numpy.ndarray
        Flux measurements corresponding to the time stamps

    flux_err : numpy.ndarray
        The flux uncertainties corresponding to the data.

    true_lightcurve : iterable containing (true_time, true_flux)
        In the case of simulated data, this contains the times and flux values from which the
        simulated data was created (could be higher-resolution than the "data"), useful for
        comparison between models created e.g. from MCMC samples and the true underlying process

    models : iterable of shape (model_time, numpy.ndarray of shape (nsamples, len(model_time)))
        First element here contains the time stamps for the models (which may not be the same
        as for the data), the second is an array of shape (nsamples, ndatapoints), where nsamples
        is the number of model light curves, and ndatapoints == len(model_time)

    ax : matplotlib.Axes object
        An Axes object in which to plot the results. If not given, the code will create a
        new figure.

    legend : bool, default True
        If True, include a legend in the plot

    colours : [str, str, str]
        List of (up to) three colours. First colour is used for the data, the second
        colour for the true underlying data, the third for the models.

    Returns
    -------

    ax : matplotlib.Axes object
        The object with the plot

    """

    if colours is None:
        colours = ["black", "#33B3FF", "#FFB733"]

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8,4))

    if flux_err is None:
        ax.scatter(time, flux, s=4, c=colours[0], marker="o")
    else:
        ax.errorbar(time, flux, yerr=flux_err, fmt="o", markersize=4, c=colours[0])

    min_time = np.min(time)
    max_time = np.max(time)

    if true_lightcurve is not None:
        true_time = true_lightcurve[0]
        true_flux = true_lightcurve[1]

        ax.plot(true_time, true_flux, lw=2, alpha=0.5, color=colours[1])

        min_time = np.min([min_time, np.min(true_time)])
        max_time = np.max([max_time, np.max(true_time)])


    if models is not None:
        m_time = models[0]
        m_flux = models[1]

        for m in m_flux:
            ax.plot(m_time, m, color=colours[2], alpha=0.1)

        min_time = np.min([min_time, np.min(m_time)])
        max_time = np.max([max_time, np.max(m_time)])

    ax.set_xlim(min_time, max_time)
    ax.set_xlabel("Time")
    ax.set_ylabel("Flux");

    return ax

def plot_lsp(time, flux, flux_err=None, true_period=None,
             p_min=1/24.0, p_max=1.0, oversampling=5, lsp_kwargs={},
             ax=None, nharmonics=0, legend=True, use_gatspy=False, nterms=1):

    """
    Plot a Lomb-Scargle periodogram of the data in (time, flux) and optionally
    include the "true" period, if known.

    Parameters
    ----------
    time : numpy.ndarray
        The time stamps of the time series **IN DAYS**

    flux : numpy.ndarray
        Flux measurements corresponding to the time stamps

    flux_err : numpy.ndarray
        The flux uncertainties corresponding to the data.

    true_period : float
        For simulated data, we might know the true period from
        which the data came. If so, this will plot a dashed line
        where the true period should be in the Lomb-Scargle periodogram

    p_min, p_max : float, float
        The minimum and maximum period to search, **IN DAYS**

    oversampling : int
        The oversampling factor; default is 5

    lsp_kwargs : dict
        Optional keyword arguments for astropy.stats.LombScargle

    ax : matplotlib.Axes object
        A matplotlib.Axes object into which to plot the LSP. Useful
        for including an LSP in a multi-panel figure. If not given,
        this function will instantiate a new Figure.

    nharmonics : int, default 0
        The number of harmonics to plot apart from the true period.
        Plotting harmonics can be useful if there's significant power
        e.g. at 1/2 the actual period

    legend : bool, default True
        If true, include a legend


    """

    freq, power = make_lsp(time, flux, flux_err, p_min, p_max, oversampling, use_gatspy=use_gatspy,
                           nterms=nterms)

    # array of corresponding periods
    periods = 1./freq

    # invert the arrays so that my plot goes from
    # small periods --> large periods
    periods = periods[::-1]
    power = power[::-1]

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8,4))

    # plot the Lomb-Scargle periodogram
    ax.semilogx(periods, power, lw=2, color="black", linestyle="steps-mid", alpha=0.7, label="data")

    ylim = ax.get_ylim()

    ax.vlines(true_period, ylim[0], ylim[1], color="purple", linestyle="dashed", lw=3,
              label="true period")

    if nharmonics > 0:
        colours = sns.color_palette("viridis", n_colors=nharmonics)
        for n in range(nharmonics):
            true_period /= 2.0
            ax.vlines(true_period, ylim[0], ylim[1], color=colours[n],
                      linestyle="dashed", lw=3, label="%.2f * true period"%(0.5/(n+1)))

    if legend:
        ax.legend()

    ax.set_xlim(np.min(periods), np.max(periods))

    ax.set_xlabel("Period [days]")
    ax.set_ylabel("Normalized Power")

    return ax



def plot_folded_lightcurve(time, flux, period, flux_err=None, models=None, true_lightcurve=None,
                      ax=None, use_radians=False, legend=True, colours=None):
    """
    Plot a folded periodic light curve, potentially including the true underlying
    model that produced the data (in the case of simulations), or model
    light curves from MCMC.

    Parameters
    ----------
    time : numpy.ndarray
        The time stamps of the periodic light curve

    flux : numpy.ndarray
        Flux measurements corresponding to the time stamps

    flux_err : numpy.ndarray
        The flux uncertainties corresponding to the data.

    period : float
        The period on which to fold **in days**

    models : iterable of shape (model_time, numpy.ndarray of shape (nsamples, len(model_time)))
        First element here contains the time stamps for the models (which may not be the same
        as for the data), the second is an array of shape (nsamples, ndatapoints), where nsamples
        is the number of model light curves, and ndatapoints == len(model_time)

    true_lightcurve : iterable containing (true_time, true_flux)
        In the case of simulated data, this contains the times and flux values from which the
        simulated data was created (could be higher-resolution than the "data"), useful for
        comparison between models created e.g. from MCMC samples and the true underlying process


    ax : matplotlib.Axes object
        An Axes object in which to plot the results. If not given, the code will create a
        new figure.

    use_radians : bool, default False
        If True, the phase will be plotted from (0, 2pi) instead of (0,1), which is the default.

    legend : bool, default True
        If True, include a legend in the plot

    colours : [str, str, str]
        List of (up to) three colours. First colour is used for the data, the second
        colour for the true underlying data, the third for the models.

    Returns
    -------

    ax : matplotlib.Axes object
        The object with the plot

    """

    if colours is None:
        colours = ["black", "#33B3FF", "#FFB733"]

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6,4))

    t0 = np.min(time)
    if models:
        t0 = np.min([t0, np.min(models[0])])

    if true_lightcurve:
        t0 = np.min([t0, np.min(true_lightcurve[0])])

    phase = (time-t0)/period - np.floor((time-t0)/period)

    if use_radians:
        phase *= 2.*np.pi

    if flux_err is None:
        ax.scatter(phase, flux, s=5, color=colours[0], label="data")
    else:
        ax.errorbar(phase, flux, yerr=flux_err, fmt="o", c=colours[0], markersize=5, label="data")

    if true_lightcurve:
        true_time = (true_lightcurve[0] - t0)
        true_flux = true_lightcurve[1]
        true_phase = true_time/period - np.floor(true_time/period)

        if use_radians:
            true_phase *= 2.*np.pi

        # compute the difference from one phase bin to the next
        tdiff = np.diff(true_phase)
        # find all differences < 0, which is where the phase wraps around
        idx = np.where(tdiff < 0)[0]

        # loop through indices where phase goes from 1 (or 2pi) to 0
        # plot each phase light curve separately
        istart = 0
        iend = idx[0]+1

        # first phase cycle also contains the label for the legend
        ax.plot(true_phase[istart:iend], true_flux[istart:iend], alpha=0.3,
                c=colours[1], label="true light curve")

        for i, x in enumerate(idx[:-1]):
            ax.plot(true_phase[istart:iend], true_flux[istart:iend], alpha=0.3, c=colours[1])
            istart = x+1
            iend = idx[i+1]+1

        # last plot
        istart = idx[-1]+1
        ax.plot(true_phase[istart:], true_flux[istart:], alpha=0.3, c=colours[1])

    if models:

        m_time = (models[0] - t0)
        m_flux = models[1]

        m_phase = (m_time/period) - np.floor(m_time/period)
        print("mphase " + str(m_phase))
        if use_radians:
            m_phase *= 2.*np.pi

        # compute the difference from one phase bin to the next
        tdiff = np.diff(m_phase)

        print("tdiff " + str(tdiff))
        # find all differences < 0, which is where the phase wraps around
        idx = np.where(tdiff < 0)[0]
        #if idx.size == 0:
        #    idx = np.array(0)

        # loop through the different samples
        for i,m in enumerate(m_flux):
            # loop through indices where phase goes from 1 (or 2pi) to 0
            # plot each phase light curve separately
            istart = 0
            iend = idx[0]+1

            if i == 0:
                # first phase cycle also contains the label for the legend
                ax.plot(m_phase[istart:iend], m[istart:iend], alpha=0.1,
                        c=colours[2], label="model")

            else:
                ax.plot(m_phase[istart:iend], m[istart:iend], alpha=0.1,
                        c=colours[2])

            for j, x in enumerate(idx[:-1]):
                ax.plot(m_phase[istart:iend], m[istart:iend], alpha=0.1,
                        c=colours[2])

                istart = x+1
                iend = idx[j+1]+1

            # last plot
            istart = idx[-1]+1
            ax.plot(m_phase[istart:], m[istart:], alpha=0.1, c=colours[2])

    if legend:
        ax.legend()
    ax.set_xlabel("Rotational Phase")
    ax.set_ylabel("Flux")
    period_hours = period*24.
    ax.set_title(r"period $P = %.3f$"%period_hours)
    if use_radians:
        ax.set_xlim(0, 2*np.pi)
    else:
        ax.set_xlim(0, 1)
    return ax

def plot_steps(sampler, dims=None, p0=None, data_pts=None):
    fig, ax = plt.subplots(2, 2, figsize=(7,6))
    fig.subplots_adjust(wspace=0.25, hspace=0.3)

    if data_pts is not None:
        fig.suptitle("Data points: " + str(data_pts) + "\nMean acceptance fraction:{0:.3f}".format(np.mean(sampler.acceptance_fraction)))

    else:
        fig.suptitle("Mean acceptance fraction: 0:.3f}".format(np.mean(sampler.acceptance_fraction)))


    axs = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]

    x = np.arange(sampler.iterations)

    for i in range(sampler.dim):
        axs[i].set_xlabel('Step Number')
        axs[i].set_ylabel('{}'.format(dims[i]))

        for j in range(len(sampler.chain)):
            param = sampler.chain[j,:,i]
            axs[i].plot(x, param, 'k-', alpha=0.3)
            # fit might guess period is time range of sampling

        flatchain = sampler.flatchain[:,i]


    return axs


def plot_mcmc_sampling_results(tsample, fsample, flux_err, gp, sampler,
                          t_pred=None, true_lightcurve=None,
                          true_period=None, namestr="test",
                          nmodels=10, npred=1000):


    new_samples = sampler.flatchain
    #new_samples = sampler.chain


    ### plot light curve with example models ###

    # first, get the total number of available samples
    nsamples = new_samples.shape[0]

    # get some random samples from the total steps
    #idx = np.random.choice(np.arange(0, nsamples, 1, dtype=int), size=nmodels)
    idx = np.rint(np.linspace(sampler.iterations-1, nsamples-1, nmodels)).astype(int)


    # if the array for the predictions isn't given, make one
    if t_pred is None:
        t_pred = np.linspace(tsample[0], tsample[-1], npred)

    # empty array for output
    m_all = np.zeros((nmodels, t_pred.shape[0]))

    # loop through the indices of samples, for each sample from the GP
    # conditional on the data points
    for i,j in enumerate(idx):
        p = new_samples[j]
        #p = new_samples[j, -1]
        print(p)
        pnew = [p[0], p[1], p[2], p[3]]

        gp.set_parameter_vector(pnew)
        mean_model = gp.sample_conditional(fsample, t_pred)
        m_all[i] = mean_model

    fig, ax = plt.subplots(1, 1, figsize=(6,4))

    plot_lightcurve(tsample, fsample, true_lightcurve=true_lightcurve,
                        models=(t_pred, m_all), ax=ax)

    plt.tight_layout()
    plt.savefig(namestr + "_lc.pdf", format="pdf")



    ### plot histogram of periods ###
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10,4))
    p_hours = np.exp(new_samples[:,-1])*24
    ax.hist(p_hours, bins=100, density=True,
                label="posterior PDF", color="black", alpha=0.5)

    if true_period is not None:
        ylim = ax.get_ylim()
        ax.vlines(true_period, 0, ylim[-1], lw=1, color="red", linestyle="dashed", label="true period : " + str(true_period))

    ax.set_xlabel("Period in hours")
    ax.set_ylabel("Probability")
    ax.set_title("Period Distibution")
    ax.legend()

    lower, upper = np.percentile(p_hours, [5,95])
    bx.hist(p_hours[(p_hours > lower) & (p_hours < upper)], bins=100, density=True,
                label="posterior PDF", color="black", alpha=0.5)
    bx.set_title("5th - 95th Percentile")
    bx.set_xlabel("Period in hours")
    bx.set_ylabel("Probability")

    if true_period is not None:
        ylim = bx.get_ylim()
        bx.vlines(true_period, 0, ylim[-1], lw=1, color="red", linestyle="dashed", label="true period : " + str(true_period))


    plt.tight_layout()
    plt.savefig(namestr + "_period_pdf.pdf", format="pdf")



    ### plot folded light curve ###

    fig, ax = plt.subplots(1, 1, figsize=(6,4))

    if true_period:
        print("true period " + str(true_period/24.))
        ax = plot_folded_lightcurve(tsample, fsample, true_period/24., flux_err=flux_err,
                          models=[t_pred, m_all],
                          true_lightcurve=true_lightcurve, ax=ax, use_radians=False)
    else:
        guess_period = np.quantile(np.exp(sampler.chain[:,-1,-1]), 0.5)
        print("guess period " + str(guess_period))

        ax = plot_folded_lightcurve(tsample, fsample, guess_period, flux_err=flux_err,
                          models=None,#[t_pred, m_all],
                          true_lightcurve=true_lightcurve, ax=ax, use_radians=False)

    plt.tight_layout()
    plt.savefig(namestr + "_folded_lc.pdf", format="pdf")


    # convert period values from log days to hours
    new_samples.T[3] = copy.deepcopy(np.exp(new_samples.T[3])*24.)
     #= x

    labels = list(gp.get_parameter_names())
    labels[3] = 'period hours'



    ### make a corner plot ###

    # percentage of steps to plot
    percent = 0.1

    # generate list of last steps spanning backwards
    x = list(range(-1,-int(sampler.iterations*percent)-1,-1))

    # only plot x last steps
    # reshape into 2d array instead of 3d
    #sampler.chain[walker, step, param]
    figure = corner.corner(sampler.chain[:,x,:].reshape(int(sampler.iterations*percent*sampler.k),4),
                           labels=labels, quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": 10})
    plt.savefig(namestr + "_corner.pdf", format="pdf")



    ### plot trace plot ###

    ax = plot_steps(sampler, dims = ['mean', 'log_amp', 'gamma', 'period'], data_pts=len(fsample))
    plt.savefig(namestr + "_trace.pdf", format="pdf")
