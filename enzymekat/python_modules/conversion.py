import numpy as np


def sem_pct_to_lognormal_sigma(sem_pct, mean, n=3):
    sem = sem_pct/100 * mean
    s = sem * np.sqrt(n)
    return (np.sqrt(np.log(1 + (s ** 2)/(mean ** 2))))
    
