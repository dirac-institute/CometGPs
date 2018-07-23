def plot_sampling_results(time, flux, flux_err, gp, sampler,
                          t_pred=None, true_lightcurve=None,
                          true_period=None, namestr="test", 
                          nmodels=10, npred=1000):
    
    
    # resample from weights
    new_samples = sampler.flatchain
    
    # make a corner plot
    corner.corner(new_samples, labels=gp.get_parameter_names())

    # save to file
    plt.savefig(namestr + "_corner.pdf", format="pdf")

    # plot some light curves with example models

    # first, get the total number of available samples
    nsamples = new_samples.shape[0]

    # get some random samples from the 

    idx = np.random.choice(np.arange(0, nsamples, 1, dtype=int), size=nmodels)

    # if the array for the predictions isn't given, make one
    if t_pred is None:
        t_pred = np.linspace(time.iloc[0], time.iloc[-1], npred)

    # empty array for output
    m_all = np.zeros((nmodels, t_pred.shape[0]))

    # loop through the indices of samples, for each sample from the GP
    # conditional on the data points
    for i,j in enumerate(idx):
        p = new_samples[j]
        pnew = [p[0], p[1], p[2], p[3]]

        gp.set_parameter_vector(pnew)
        mean_model = gp.sample_conditional(fsample, t_pred)
        m_all[i] = mean_model

    fig, ax = plt.subplots(1, 1, figsize=(6,4))
    plot_lightcurve(tsample, fsample, true_lightcurve=true_lightcurve, 
                        models=(t_pred, m_all), ax=ax)

    plt.tight_layout()
    plt.savefig(namestr + "_lc.pdf", format="pdf")

    # plot histogram of periods
    fig, ax = plt.subplots(1, 1, figsize=(5,4))
    ax.hist(np.exp(new_samples[:,-1])*24, bins=100, normed=True, 
                label="posterior PDF", color="black", alpha=0.5)

    if true_period is not None:
        ylim = ax.get_ylim()
        ax.vlines(true_period, 0, ylim[-1], lw=1, color="red", linestyle="dashed", label="true period : " + str(true_period))

    ax.set_xlabel("Period in hours")
    ax.set_ylabel("Probability")
    ax.legend()

    plt.tight_layout()
    plt.savefig(namestr + "_period_pdf.pdf", format="pdf")

    # plot folded light curve

    fig, ax = plt.subplots(1, 1, figsize=(6,4))


    if true_period:
        ax = plot_folded_lightcurve(tsample, fsample, true_period/24, flux_err=0.01, 
                          models=[t_pred, m_all[:2]], 
                          true_lightcurve=true_lightcurve, ax=ax, use_radians=False)
    else:
        ax = plot_folded_lightcurve(tsample, fsample, best_period, flux_err=ferr, 
                          models=[t_pred, m_all[:2]], 
                          true_lightcurve=true_lightcurve, ax=bx, use_radians=False)

    plt.tight_layout()
    plt.savefig(namestr + "_folded_lc.pdf", format="pdf")
