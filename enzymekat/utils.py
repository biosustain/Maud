import numpy as np
import os
import pandas as pd
from typing import Iterable


def match_string_to_file(s: str, path: str):
    if os.path.exists(path):
        return open(path, 'r').read() == s
    else:
        return False


def sem_pct_to_lognormal_sigma(sem_pct, mean, n=3):
    sem = sem_pct/100 * mean
    s = sem * np.sqrt(n)
    return (np.sqrt(np.log(1 + (s ** 2)/(mean ** 2))))
    

def codify(l: Iterable[str]):
    return dict(zip(l, range(1, len(l) + 1)))
