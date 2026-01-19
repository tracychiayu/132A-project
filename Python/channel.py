import numpy as np

def gen_rayleigh_channel_matrix(mu, var, Nr, Nt):
    '''
    Generate a Nr x Nt channel matrix, where each entry h_ij follows a complex Gaussian distribution with mean 'mu' and variance 'var'
    Input
        - mu: mean of complex Gaussian distribution
        - var: variance of complex Gaussian distribution
        - Nr: number of receive antennas (num of rows of H)
        - Nt: number of transmit antennas (num of cols of H)
    Output
        Nr x Nt channel matrix
    '''

    # Re{h_ij} and Im{h_ij} follow Gaussian distribution with mean 'mu' and variance 'var/2'
    return np.random.normal(loc=mu, scale=np.sqrt(var/2), size=(Nr, Nt)) + 1j * np.random.normal(loc=mu, scale=np.sqrt(var/2), size=(Nr, Nt))
