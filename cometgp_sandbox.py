import george
from george import kernels
import numpy as np
import matplotlib.pyplot as plt

def create_data(gamma = 10, period = 1, amp = 1):

    #gamma = 10
    #period = 1
    #amp = 1

    #generate 10 random numbers for x
    pre_x = 10 * np.sort(np.random.rand(10))
    pre_yerr = 0.2 * amp * np.ones_like(pre_x)
    x_possible = np.linspace(0,10,500)

    #create a sinusoidal plot based on inputs
    kernel = amp * kernels.ExpSine2Kernel(gamma,period)
    gp = george.GP(kernel)
    gp.compute(pre_x, pre_yerr)

    pre_y = gp.sample(x_possible) #a subset of possible x-values
    return(pre_x, pre_y, x_possible, gp)

def sample_data(pre_x, pre_y, x_possible, n, m):

    assert n != 0 or m != 0, 'must specify a number of sample points'

    if n != 0:
        #pick n number of random points from 0-10
        idx_choice = random_set(pre_y, n)

    if m != 0:
        idx_choice = cluster_set(pre_y, m)

    x = x_possible[idx_choice]
    #get the predicted y values corresponding to x values
    y_base = pre_y[idx_choice]
    yerr = 0.2 * np.ones_like(x)
    #add noise to initial y values
    y = y_base + yerr * np.random.randn(len(x))

    return(x, y, yerr)

def fit_data(x, y, yerr, gamma, period, gp, print_results = True):
    gp.compute(x, yerr)
    x_possible = np.linspace(0,10,500)
    pred, pred_var = gp.predict(y, x_possible, return_var=True)
    ln_likelihood_initial = gp.log_likelihood(y)

    from scipy.optimize import minimize

    def neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.log_likelihood(y)

    def grad_neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.grad_log_likelihood(y)

    result = minimize(neg_ln_like, gp.get_parameter_vector(), jac=grad_neg_ln_like)
    gp.set_parameter_vector(result.x)

    ln_likelihood_final = gp.log_likelihood(y)
    if (print_results == True):
        print(ln_likelihood_initial, result, ln_likelihood_final)

    return gp

def plotting( x, y, yerr,pre_y,x_possible, gp):
    """"""
    pred, pred_var = gp.predict(y, x_possible, return_var=True)
    #plt.figure(figsize=(15,10))
    plt.fill_between(x_possible, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var),
                color="red", alpha=0.4)
    plt.plot(x_possible, pred, "red", lw=1.5, alpha=0.7)
    plt.errorbar(x, y, yerr = yerr, fmt=".k", capsize=0)
    plt.plot(x_possible,pre_y)
    plt.xlim(0, 10)
    plt.ylim(pre_y.min()-0.5, pre_y.max()+0.5)
    plt.xlabel("x")
    plt.ylabel("y");
    plt
    plt.show()

def random_set(pre_y, n):
    idx = np.arange(len(pre_y), dtype=int)
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
        idx_sub.append(np.arange(idx_start[i]-50,idx_start[i]+50))

    idx_choice = np.hstack((idx_sub))
    idx_choice.sort()

    return(idx_choice)
