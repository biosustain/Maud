import pandas as pd
import os


def match_string_to_file(s: str, path: str):
    if os.path.exists(path):
        return open(path, 'r').read() == s
    else:
        return False
