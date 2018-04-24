import george
from george import kernels
import numpy as np
import matplotlib.pyplot as plt

def create_data(gamma = 10, period = 1, amp = 1):
    """Create a randomized periodic-sinusoidal function.

    Some example inputs are:
    gamma = 10
    period = 1
    amp = 1
    """

    #generate 10 random numbers for x
    pre_x = 10 * np.sort(np.random.rand(10))
    #determine the y error and make it an array as long as pre_x
    #have the error scale with the amplitude
    pre_yerr = 0.2 * amp * np.ones_like(pre_x)
    x_possible = np.linspace(0,10,500)

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

def sample_data(x_possible, pre_y, yerr_amp, n, m):
    """Samples a random array of points from the sample data generated in
    the create_data function or another type of light curve sample.

    Sample must be presented as an array of possible:
    x-values (x_possible)
    y-values (pre_y)
     """
    #abort if they don't specify how many points/clusters they want to sample
    assert n != 0 or m != 0, 'must specify a number of sample points'

    if n != 0:
        #pick n number of random points from
        idx_choice = random_set(pre_y, n)

    if m != 0:
        #pick m number of random clusters
        idx_choice = cluster_set(pre_y, m)

    #filter out the randomly chosen indicies from earlier
    x = x_possible[idx_choice]
    #get the predicted y values corresponding to x values
    y_base = pre_y[idx_choice]
    yerr = yerr_amp * 0.01 * np.ones_like(x) #turn 0.2 into an input (1%)
    #add noise to initial y values
    y = y_base + yerr * np.random.randn(len(x))

    return(x, y, yerr)

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
        return -gp.log_likelihood(y)

    def grad_neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.grad_log_likelihood(y)

    result = minimize(neg_ln_like, gp.get_parameter_vector(), jac=grad_neg_ln_like)
    gp.set_parameter_vector(result.x)

    ln_likelihood_opt = gp.log_likelihood(y)
    if (print_results == True):
        print(ln_likelihood, result, ln_likelihood_opt)

    return gp, ln_likelihood_opt, result.fun

def plotting( x, y, yerr,pre_y,x_possible, gp):
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
    plt
    plt.show()

def random_set(pre_y, n):
    idx = pre_y.size #pre_y is a long as x_possible (500)
    #randomly select n points from 1-500
    idx_choice = np.random.choice(idx,n, replace= False)
    idx_choice.sort()
    return(idx_choice)

def cluster_set(pre_y, m):
    idx = pre_y.size
    idx_start = []
    for j in np.arange(m):
        idx_start.append(np.random.choice(np.arange(int((j)*idx/m),int((j+1)*idx/m)),1)[0]) #get 5 before and 5 after
    idx_sub = np.hstack((idx_start))
    idx_sub.sort()

    idx_sub = []
    for i in np.arange(len(idx_start)):
        idx_sub.append(np.arange(idx_start[i]-5,idx_start[i]+5))

    idx_choice = np.hstack((idx_sub))
    idx_choice.sort()
    return(idx_choice)
