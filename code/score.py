""" score.py -- Compute individual feature scores.
"""

import argparse
import logging
import numpy as np
import os
import scipy.stats as st 
import sys

class ScoreFeatures(object):
    """ Class for scoring individual features.

    Attributes
    ----------
    X: (num_features, num_samples) np.array
        Design matrix (genotype).
    y: (num_samples, ) np.array
        Outcome (phenotype).
    """

    def __init__(self, X, y):
        """
        Parameters
        ----------
        X: (num_features, num_samples) np.array
            Design matrix (genotype).
        y: (num_samples, ) np.array
            Outcome (phenotype).
        """
        self.X = X
        self.y = y


    def compute_skat(self):
        """
        Compute individual linear SKAT scores for each feature.

        Returns
        -------
        scores: (num_features, ) np.array
            Linear SKAT score for each feature.
        """
        scores = (np.dot(self.X, (self.y - np.mean(self.y))))**2
        return scores
        
        
        
