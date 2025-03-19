import jax.numpy as jnp
import numpy as np

from scipy.stats import bernoulli,norm

#import copula functions
from pr_copula.main_copula_density import fit_copula_density,predict_copula_density,predictive_resample_density,check_convergence_pr

    

def run_copula_algo(y):
    
    mean_y = jnp.mean(y)
    std_y = jnp.std(y)
    y = (y-mean_y)/std_y

    #Plot of true pdf
    y_plot = jnp.arange(-4,4,0.05).reshape(-1,1)
    dy = y_plot[1]-y_plot[0]
    true_pdf = (0.8*norm.pdf(std_y*y_plot+mean_y,loc = -2) + 0.2*norm.pdf(std_y*y_plot+mean_y,loc = 2))*std_y

    #Fit copula obj
    copula_density_obj = fit_copula_density(y,seed = 200)
    print('Bandwidth is {}'.format(copula_density_obj.rho_opt))
    print('Preq loglik is {}'.format(copula_density_obj.preq_loglik))

    #small modifications: we prevent the bandiwdith to generate
    #copula_density_obj.rho_opt = jnp.minimum(copula_density_obj.rho_opt, 0.5)
    rho_new = jnp.maximum(copula_density_obj.rho_opt, 0.5)
    copula_density_obj = copula_density_obj._replace(rho_opt=rho_new)
    print('New Bandwidth is {}'.format(copula_density_obj.rho_opt))
    #Predict on yplot
    logcdf_conditionals,logpdf_joints = predict_copula_density(copula_density_obj,y_plot)
    pdf_cop = jnp.exp(logpdf_joints[:,-1])
    cdf_cop = jnp.exp(logcdf_conditionals[:,-1])

    #Predictive resample
    T_fwdsamples = 5000 #T = N-n
    B_postsamples = 1000
    logcdf_conditionals_pr,logpdf_joints_pr= predictive_resample_density(copula_density_obj,y_plot,B_postsamples, T_fwdsamples, seed = 200) #(seed = 200)


    return logcdf_conditionals_pr,logpdf_joints_pr
    


import sys

# Read command-line arguments
#i = sys.argv[1]  # First argument
i=1
n=200

filename = f"data/data_n{n}_{i}.csv"
y = np.loadtxt(filename, delimiter=",", dtype=str)
y=y.astype(np.float32)
y=y.reshape(-1,1)
n = y.shape[0]


logcdf_conditionals_pr,logpdf_joints_pr = run_copula_algo(y)
nameout = f"output_py/out_python_n{n}_{i}"
jnp.save(nameout.format(n),logpdf_joints_pr)
