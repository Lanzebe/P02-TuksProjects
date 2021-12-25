# -*- coding: utf-8 -*-
from numpy import array

def read_input_file(filename):
    """ Construct a dictionary of KEYWORDS and corresponding data

    KEYWORD is defined by *TEXT, the TEXT is then used to created a key in the dictionary
    followed by the data below the KEYWORD.

    Each KEY in the dictionary is associated with a list of lists, including scalars.
    """
    data_dictionary = {}
    KEYWORD = None
    with open(filename) as f:
        for line in f:
            try:
                row = [float(item) for item in line.split()]
                data_dictionary[KEYWORD].append(row)
            except ValueError:
                KEYWORD = line.strip()[1:]
                data_dictionary[KEYWORD] = []

# Convert all lists to numpy arrays for convenience                
#    for key in data_dictionary:
#        data_dictionary[key] = array(data_dictionary[key])
    return data_dictionary