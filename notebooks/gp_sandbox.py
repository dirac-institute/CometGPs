import george
from george import kernels
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import corner 

def create_data(gamma = 10, period = 1, amp = 1):
    """Create a randomized periodic-sinusoidal function.

    Some example inputs are:
    gamma = 10
    period = 1
    amp = 1
    
    Returns all x-values with corresponding y-values and gp 
    function generated with the given gamma, period, and amplitude.
    """

    #generate 10 random numbers for x
    pre_x = 10 * np.sort(np.random.rand(10))
    #determine the y error and make it an array as long as pre_x
    #have the error scale with the amplitude
    pre_yerr = 0.2 * amp * np.ones_like(pre_x)
    x_possible = np.linspace(0,10, 1000)

    #create a sinusoidal plot based on inputs
    #establish the kernel type we are using
    kernel = amp * kernels.ExpSine2Kernel(gamma,period)
    gp = george.GP(kernel)
    #generate a ExpSine2 function using our x-values (0-10) and our y error
    gp.compute(pre_x, pre_yerr)
    #note: pre_x is never used again after this

    #sample the y-values from the function we made
    pre_y = gp.sample(x_possible) #a subset of possible x-values
    #return our simulated data with the original function
    return(x_possible, pre_y, gp)

def sample_data(x_possible, pre_y, yerr_amp, random_n=0, cluster_n=0, cadence_n=0):
    """Samples a random array of points from the sample data generated in
    the create_data function or another type of light curve sample.

    Sample must be presented as an array of possible:
    x-values (x_possible)
    y-values (pre_y)
     """
    #abort if they don't specify how many points/clusters they want to sample
    assert random_n != 0 or cluster_n != 0 or cadence_n != 0, 'must specify a number of sample points'

    if random_n != 0:
        #pick n number of random points from
        idx_choice = random_set(pre_y, random_n)

    if cluster_n != 0:
        #pick m number of random clusters
        idx_choice = cluster_set(pre_y, cluster_n)

    if cadence_n != 0:
        #pick m number of random clusters
        idx_choice = cadence_set(x_possible, cadence_n)

    #filter out the randomly chosen indicies from earlier
    x = x_possible[idx_choice]
    #get the predicted y values corresponding to x values
    y_base = pre_y[idx_choice]
    yerr = yerr_amp * 0.01 * np.ones_like(x) #turn 0.2 into an input (1%)
    #add noise to initial y values
    y = y_base + yerr * np.random.randn(len(x))

    return(x, y, yerr)

def random_set(pre_y, random_n):
    idx = pre_y.size #pre_y is a long as x_possible (500)
    #randomly select n points from 1-500
    idx_choice = np.random.choice(idx, random_n, replace= False)
    idx_choice.sort()
    return(idx_choice)

def cluster_set(pre_y, cluster_n):
    idx = pre_y.size
    start_idx = []
    for j in np.arange(cluster_n):
        #get 5 before and 5 after
        start_idx.append(np.random.choice(np.arange(int((j)*idx/cluster_n),int((j+1)*idx/cluster_n)),1)[0]) 
    idx_sub = np.hstack((start_idx))
    idx_sub.sort()

    idx_sub = []
    for i in np.arange(len(start_idx)):
        idx_sub.append(np.arange(start_idx[i]-5,start_idx[i]+5))

    idx_choice = np.hstack((idx_sub))
    idx_choice.sort()
    return(idx_choice)

def cadence_set(x_possible, days):
    max_idx = x_possible.size
    start_idx = x_possible.index[0]
    init_idx = x_possible.index[0]
    end_idx = start_idx + 48*20 #8 hours later
    idx_choice = []

    while end_idx <= max_idx+init_idx:
        idx_choice.extend(np.arange(start_idx, end_idx, 20, dtype=int))
        start_idx += 2880 #skip forward 24 hours
        end_idx += 2880
    return(idx_choice)


def log_prior(gamma, period, amp):
    if(gamma > 1.0 and gamma < 15.0): p_gamma = 1
    else: p_gamma = 0

    print('gamma prior: ' + str(p_gamma))

    if(period > 1.0 and period < 10.0): p_period = 1
    else: p_period = 0

    print('period prior: ' + str(p_period))

    if(amp > 1.0 and amp < 10.0): p_amp = 1
    else: p_amp = 0

    print('amp prior: ' + str(p_amp))

    p_prior = p_gamma * p_period * p_amp
    ln_prior = np.log(p_prior)
    return ln_prior


def log_post(ln_prior, ln_likelihood):
    ln_post = ln_prior + ln_likelihood
    return ln_post


def fit_data(x, y, yerr, gamma, period, gp):
    gp.compute(x, yerr)
    x_possible = np.linspace(0,10,500)
    pred, pred_var = gp.predict(y, x_possible, return_var=True)
    ln_likelihood = gp.log_likelihood(y)

    #print(ln_likelihood)

    return gp, ln_likelihood


def optimize(y, gp, ln_likelihood, print_results = True):

    from scipy.optimize import minimize

    def neg_ln_like(p):
        gp.set_parameter_vector(p)
        try:
            negloglike =  -gp.log_likelihood(y)
            return negloglike
        except scipy.linalg.LinAlgError:
            return np.inf

    #print(neg_ln_like)

    def grad_neg_ln_like(p):
        gp.set_parameter_vector(p)
        try:
            grad_loglike =  -gp.grad_log_likelihood(y)
            return grad_loglike
        except scipy.linalg.LinAlgError:
            return np.inf
        #return -gp.grad_log_likelihood(y)

    result = minimize(neg_ln_like, gp.get_parameter_vector(), jac=grad_neg_ln_like) #, method='L-BFGS-B')
    gp.set_parameter_vector(result.x)

    ln_likelihood_opt = gp.log_likelihood(y)
    if (print_results == True):
        print(ln_likelihood, result, ln_likelihood_opt)

    return gp, ln_likelihood_opt, result.fun



def plotting(x, y, yerr,pre_y, x_possible, gp):
    """"""
    pred, pred_var = gp.predict(y, x_possible, return_var=True)
    #plt.figure(figsize=(15,10))
    plt.fill_between(x_possible, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var),
                color="red", alpha=0.4)
    plt.plot(x_possible, pred, "red", lw=1.5, alpha=0.7)
    plt.errorbar(x, y, yerr = yerr, fmt=".k", capsize=0)
    plt.plot(x_possible,pre_y)
    plt.xlim(x_possible.min(), x_possible.max())
    plt.ylim(pre_y.min()-0.5, pre_y.max()+0.5)
    plt.xlabel("x")
    plt.ylabel("y");
    plt.show()


def plot_steps(sampler, p0=None, data_pts=None, from_saved=False):
    fig, ax = plt.subplots(2,2)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    
    axs = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]
    dims = ['mean', 'log_amp', 'gamma', 'period']
    
    if from_saved==True:
        x = np.arange(sampler.shape[1])
        
        for i in range(len(dims)):
            axs[i].set_xlabel('Step Number')
            axs[i].set_ylabel('{}'.format(dims[i]))
            
            for j in range(sampler.shape[0]):
                param = sampler[j,:,i]
                if i == 3: 
                    param = np.exp(param)*24
                axs[i].plot(x, param, 'k-', alpha=0.3)
                
            flatchain = sampler[:,:,i]
            if i == 3: 
                pre_mean = flatchain.mean()
                flatchain = np.exp(flatchain)*24
                axs[i].axhline(flatchain.mean(), linestyle='--', 
                               label=(str((round(flatchain.mean(),5)))+'h \n'+str((round(pre_mean,5)))))
            
            else: axs[i].axhline(flatchain.mean(), linestyle='--', 
                                 label=round(flatchain.mean(),5))
            axs[i].legend(loc=1)

    else:        
        x = np.arange(sampler.iterations)

        print(str(p0[0]) + '\nData points: ' + str(data_pts))
        print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

        for i in range(len(dims)):
            axs[i].set_xlabel('Step Number')
            axs[i].set_ylabel('{}'.format(dims[i]))
            
            for j in range(len(sampler.chain)):
                param = sampler.chain[j,:,i]
                if i == 3: 
                    param =  np.exp(param)*24
                axs[i].plot(x, param, 'k-', alpha=0.3)
                # fit might guess period is time range of sampling
                
            flatchain = sampler.flatchain[:,i]
            if i == 3: 
                pre_mean = flatchain.mean()
                flatchain = np.exp(flatchain)*24
                axs[i].axhline(flatchain.mean(), linestyle='--' , label=(str((round(flatchain.mean(),5)))+'h \n'+str((round(pre_mean,5)))))
            
            else: axs[i].axhline(flatchain.mean(), linestyle='--' , label=round(flatchain.mean(),5))
            axs[i].legend(loc=1)
            
def plot_hist(sampler, from_saved=False):
    """
    Working. :)
    """
    fig, ax = plt.subplots(2,2)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    dims = ['mean', 'log_amp', 'gamma', 'period']
    axs = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]
    
    if from_saved == True:
        for i in range(len(dims)):
            axs[i].set_xlabel('{}'.format(dims[i]))
            param = sampler[:,:,i]
            if i == 3: 
                param =  np.exp(param)*24
            axs[i].hist(param, 100,histtype="step")
        
    else:
        for i in range(len(dims)):
            axs[i].set_xlabel('{}'.format(dims[i]))
            param = sampler.flatchain[:,i]
            if i == 3: 
                param =  np.exp(param)*24
            axs[i].hist(param, 100, color='k', histtype="step")
                
def plot_mcmc_sampling_results(time, flux, flux_err, gp, sampler,
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

