"""evaluation_framework.py -- All that is needed to evaluate feature selection algorithms."""

import numpy as np

from sklearn import cross_validation as skcv 


def consistency_index(sel1, sel2, num_features):
    """ Compute the consistency index between two sets of features.

    Parameters
    ----------
    sel1: set
        First set of indices of selected features
    sel2: set
        Second set of indices of selected features
    num_features: int
        Total number of features

    Returns
    -------
    cidx: float
        Consistency index between the two sets.

    Reference
    ---------
    Kuncheva, L.I. (2007). A Stability Index for Feature Selection.
    AIAC, pp. 390--395.
    """
    observed = float(len(sel1.intersection(sel2)))
    expected = len(sel1) * len(sel2) / float(num_features)
    maxposbl = float(min(len(sel1), len(sel2)))
    cidx = 0.
    if expected != maxposbl:
        cidx = (observed - expected) / (maxposbl - expected)
    return cidx


def consistency_index_k(sel_list, num_features):
    """ Compute the consistency index between more than 2 sets of features.

    This is done by averaging over all pairwise consistency indices.

    Parameters
    ----------
    sel_list: list of lists
        List of k lists of indices of selected features
    num_features: int
        Total number of features

    Returns
    -------
    cidx: float
        Consistency index between the k sets.

    Reference
    ---------
    Kuncheva, L.I. (2007). A Stability Index for Feature Selection.
    AIAC, pp. 390--395.
    """
    cidx = 0.
    for k1, sel1 in enumerate(sel_list):
        for sel2 in sel_list[k1+1, :]:
            cidx += consistency_index(set(sel1), set(sel2), num_features)
    cidx = 2. / (len(sel_list) * (len(sel_list) - 1)) * cidx
    return cidx



class Framework(object):
    """ Setting up evaluation framework.

    Attributes
    ----------
    self.num_samples: int
        Number of samples.
    self.num_folds: int
        Number of cross-validation folds
    self.num_subsamples: int
        Number of subsamples (to evaluate stability)
    self.xp_indices: dictionary
        fold_idx:{
            'trIndices': list of train indices,
            'teIndices': list of test indices,
            'ssIndices': list of list of subsample indices}

    """
    def __init__(self, num_samples, num_folds, num_subsamples):
        """
        Parameters
        ----------
        num_samples: int
            Number of samples.
        num_folds: int
            Number of cross-validation folds
        num_subsamples: int
            Number of subsamples (to evaluate stability)
        """

        
    def compute_indices(self, seed=None):
        """ Compute the cross-validation folds and subsample indices.

        Parameters
        ----------
        seed: {int, None}, optional
            random seed.
            Will always return the same with the same random seed.

        Modified attributes
        -------------------
        xp_indices: dictionary
            fold_idx:{
                'trIndices': list of train indices,
                'teIndices': list of test indices,
                'ssIndices': list of list of subsample indices}
        """
        # use skcv
        # Generate cross-validation indices

        # For each train set, generate self.num_subsamples subsample sets of indices
        

        
    def save_indices(self, data_dir, simu_id):
        """ Save the cross-validation folds and subsample indices to files.

        Parameters
        ----------

        Generated files
        ---------------
        For each fold_idx:
            <data_dir>/<simu_id>.<fold_idx>.trIndices:
                Space-separated list of training indices.
            <data_dir>/<simu_id>.<fold_idx>.teIndices:
                Space-separated list of test indices.
            For each subsample_idx:
                <data_dir>/<simu_id>.<fold_idx>.<ss_idx>.ssIndices
                    Space-separated lists of subsample indices,
                    one line per list / subsample.
        """
        # use np.savetxt 