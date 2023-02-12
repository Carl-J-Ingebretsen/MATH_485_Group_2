#Carl Ingebretsen Code for MATH485 work for COVID modelling

import numpy as np
import scipy as sc
from matplotlib import pyplot as plt
from scipy.integrate import odeint


def main():
    '''This is the main function for the entire program'''

    #define the parameters of the model
    alpha=0.570 #transmission rate S to I
    beta = 0.011 #transmission rate S to D
    gamma = 0.456 #transmission rate S to A
    delta = 0.011 #transmission rate S to R
    epsilon = 0.171 #Prob rate of detection 
    zeta = 0.125 #
    lambda_param = 0.034
    eta = 0.125
    rho = 0.034
    theta = 0.371 #Prob rate of detection
    mu = 0.017
    kappa = 0.017
    nu = 0.027
    xi = 0.017
    sigma = 0.017
    tau = 0.01

    #package the parameters as an tuple
    param_vec =(alpha,beta,gamma,delta,epsilon,zeta,lambda_param,eta,rho,theta,mu,kappa,nu,xi,sigma,tau)

    #Define the extra parameters to simplify the model
    r_1 = epsilon+zeta+lambda_param
    r_2 = eta+rho
    r_3 = theta+mu+kappa
    r_4 = nu+xi
    r_5 = sigma+tau

    #define the basic reporduction number
    Reprod_0 = alpha/r_1+beta*epsilon/(r_1*r_2)+gamma*zeta/(r_1*r_3)+delta*eta*epsilon/(r_1*r_2*r_4)+delta*zeta*theta/(r_1*r_3*r_4)
    print("Basic Reproduction number: ", Reprod_0)
    #Checks out with paper value

    #Define a time coordinate
    t_final = 200 #in days
    dt = 800 #number of timesteps as integer
    t = np.linspace(0,t_final,dt)

    #Define an initial condition vector
    population = 60e6
    I0=200/population
    D0=20/population
    A0=1/population
    R0=2/population
    S0=population-I0-D0-A0-R0
    T0=0
    H0=0
    E0=0
    y0 = np.array([S0,I0,D0,A0,R0,T0,H0,E0])

    #now integrate with scipy
    solution_1=odeint(system_of_equations,y0,t,args=param_vec)
    print("Complete")
    graph_results(t,solution_1,0)

def system_of_equations(y, t, *param_vec):
    '''This is the system of eight differential equations that describe the model.
    Inputs:
        y: is a numpy array of the eight differencial equations under consideration
        t: is a time coordinate to be integrated
        param_vec: is the array of the paramters to be used in the equations
        
        Result:
        dydt: the array of derivatives of each quantity
        '''
    alpha=param_vec[0] #transmission rate S to I
    beta = param_vec[1] #transmission rate S to D
    gamma = param_vec[2] #transmission rate S to A
    delta = param_vec[3] #transmission rate S to R
    epsilon = param_vec[4] #Prob rate of detection 
    zeta = param_vec[5] #
    lambda_param = param_vec[6]
    eta = param_vec[7]
    rho = param_vec[8]
    theta = param_vec[9] #Prob rate of detection
    mu = param_vec[10]
    kappa = param_vec[11]
    nu = param_vec[12]
    xi = param_vec[13]
    sigma = param_vec[14]
    tau = param_vec[15]

    #Compute the system
    dS = -y[0]*(alpha*y[1]+beta*y[2]+gamma*y[3]+delta*y[4])
    dI = y[0]*(alpha*y[1]+beta*y[2]+gamma*y[3]+delta*y[4])-(epsilon+zeta+lambda_param)*y[1]
    dD = epsilon*y[1]-(eta+rho)*y[2]
    dA = zeta*y[1]-(theta+mu+kappa)*y[3]
    dR = eta*y[2]+theta*y[3]-(nu+xi)*y[4]
    dT = mu*y[3]+nu*y[4]-(sigma+tau)*y[5]
    dH = lambda_param*y[1]+rho*y[2]+kappa*y[3]+xi*y[4]+sigma*y[5]
    dE = tau*y[5]

    #repackage as a new vector
    dydt = np.array([dS,dI,dD,dA,dR,dT,dH,dE])
    return dydt

def check_to_change_R():
    '''This is the function to change R in the program'''

def graph_results(t,y,num):
    '''This function graphs the reuslts for the variable num'''

    plt.plot(t,y[:,num])
    plt.show()

main()