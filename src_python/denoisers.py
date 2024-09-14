import mne
import scipy
import matlab.engine
import numpy as np 
from sklearn.base import TransformerMixin, BaseEstimator
from functools import reduce
from mne.stats import combine_adjacency

def transpose_indices(indices, array_shape, axes):
    temp = np.zeros(array_shape)
    temp[indices] = 1
    temp = temp.transpose(axes)
    return np.where(temp>0)

class WaveletPacketDenoiser(BaseEstimator, TransformerMixin):
    '''
        Initializes the WaveletPacketDenoiser.

        Parameters:
        - wavelet: str, default 'db4'
            The wavelet function to be used.

        - level: int or None, default None
            The level of wavelet decomposition.
            Should be between 1 and min(4, log2(number of time points per trial)).
            If level=None, the level will be determined from X during fit: 
            level=min(4, int(np.log2(X.shape[2]))).

        - coef_selection_method: str, default 'union'
            Method for combining coefficient selections.
            Possible values: 'ignore', 'union', 'intersect', 'all'.
            'ignore': ignore classes, produce only one selection from all trials.
            'all': do not combine, keep all selections separated.

        - y_sel: list or None, default None
            The classes to be considered when calculating coefficient selectors.
            Used together with coef_selection_method.

        - stat_fun: function or None, default None
            The statistical function to be used for cluster testing.
            Used as stat_fun parameter in mne.stats.spatio_temporal_cluster_1samp_test 
            (https://mne.tools/stable/generated/mne.stats.spatio_temporal_cluster_1samp_test.html)

        - n_perm: int, default 1000
            The number of permutations to be used for cluster testing.

        - ch_adjacency: None or adjacency matrix, default None
            The adjacency matrix for channels.
            See: https://mne.tools/stable/generated/mne.channels.find_ch_adjacency.html

        - test_p_threshold: float, default 0.05
            The threshold for element-wise statistical testing for cluster forming.
            A threshold parameter for the mne.stats.spatio_temporal_cluster_1samp_test method 
            is automatically computed corresponding to a p-value of test_p_threshold based on the number of trials.

        - p_threshold: float, default 0.05
            The threshold for selecting significant clusters.

        - n_jobs: int, default -1
            The number of parallel jobs to run for statistical testing. 
            n_jobs=-1 means all logical cores on the computer will be used.

        - cluster_test_kwargs: dict, default {}
            Additional keyword arguments for cluster testing.
            See: https://mne.tools/stable/generated/mne.stats.spatio_temporal_cluster_1samp_test.html

        - seed: int or None, default = 42
            The seed for controlling reproducability for the statistical test.
            See: https://mne.tools/stable/generated/mne.stats.spatio_temporal_cluster_1samp_test.html

        - matlab_eng: Matlab engine or None, default None
            The Matlab engine object. If None, a new session will be created.
    '''
    
    def __init__(
        self,
        wavelet='db4',
        level=4,
        coef_selection_method='union',
        y_sel=None,
        stat_fun=None,
        n_perm=1000,
        ch_adjacency=None,
        test_p_threshold=0.05,
        p_threshold=0.05,
        n_jobs=-1,
        cluster_test_kwargs={},
        seed=42,
        matlab_eng=None
    ):
        assert coef_selection_method in ['union', 'intersect', 'ignore', 'all'], \
            "coef_selection_method must be one of ('union', 'intersect', 'ignore', 'all')"

        self.wavelet = wavelet
        self.level = level
        self.coef_selection_method = coef_selection_method
        self.y_sel = y_sel
        self.n_perm = n_perm
        self.ch_adjacency = ch_adjacency
        self.test_p_threshold = test_p_threshold
        self.p_threshold = p_threshold
        self.n_jobs = n_jobs
        self.seed = seed
        self.matlab_eng = matlab_eng
        self.stat_fun = stat_fun
        self.cluster_test_kwargs = cluster_test_kwargs
        
    def fit(self, X, y=None):
        '''
        Fits the WaveletPacketDenoiser model.

        Parameters:
        - X: array-like of shape (n_trials, n_channels, n_times)
            The input data.

        - y: array-like of shape (n_trials,), default None
            The target labels.

        Returns:
        - self: WaveletPacketDenoiser
            The fitted WaveletPacketDenoiser model.
        '''
        eng = self.matlab_eng if self.matlab_eng is not None else matlab.engine.start_matlab()

        if y is None:
            y = np.zeros((1,X.shape[0]))

        classes = np.unique(y)
        self.n_times_ = int(X.shape[2])
        self.n_chs_ = int(X.shape[1])
        self.level_ = float(self.level) if not self.level is None else float(min(4, int(np.log2(X.shape[2]))))
        self.n_freqs_ = int(2**self.level_)
        y_sel = self.y_sel if not self.y_sel is None else classes
        self.clusters_ = []
        self.cluster_p_values_ = []
        
        # compute adjacency
        self.adjacency_ = combine_adjacency(self.n_times_, self.n_freqs_, self.ch_adjacency)
        
        _,_,self.LoD_,self.HiD_ = eng.wfilters(self.wavelet, nargout=4)
        
        # obtain wavelet denoising selectors for all classes
        self.sel_all_ = []
        self.sel_all_p_ = []
        self.n_clusters_ = 0
        for c in classes:
            if not c in y_sel:
                continue
            
            X_class = X[y==c] if not self.coef_selection_method == 'ignore' else X
            
            # obtain wavelet packet decomposition for all trials and channels
            X_dwpt = self.decompose(X_class, eng) # trial x channels x frequencies x times
            X_dwpt = np.transpose(X_dwpt, (0,3,2,1)) # trial x times x frequencies x channels for permutation clustering api
            
            threshold = 0.05 if self.test_p_threshold is None else self.test_p_threshold
            t_threshold = scipy.stats.distributions.t.ppf(1 - threshold / 2, df=X_class.shape[0]-1) # two-tailed threshold
            
            _, clusters, cluster_p_values, _ = \
                mne.stats.spatio_temporal_cluster_1samp_test(
                X_dwpt, 
                adjacency=self.adjacency_, 
                n_jobs=self.n_jobs, 
                stat_fun=self.stat_fun,
                threshold=t_threshold, 
                out_type='indices', 
                n_permutations=self.n_perm, 
                verbose=0, 
                seed=self.seed, 
                **self.cluster_test_kwargs
            )

            self.clusters_.append(clusters)
            self.cluster_p_values_.append(cluster_p_values)

            good_clusters = []
            good_clusters_p_values = []
            for indices, p_value in zip(clusters, cluster_p_values):
                if p_value < self.p_threshold:
                    # transpose the indices to fit normal channels x frequencies x times shape
                    good_clusters.append(transpose_indices(indices, (self.n_times_, self.n_freqs_, self.n_chs_), (2,1,0)))
                    good_clusters_p_values.append(p_value)
            
            self.sel_all_.append(good_clusters)
            self.sel_all_p_.append(good_clusters_p_values)
            self.n_clusters_ += 1
            
            if self.coef_selection_method == 'ignore':
                break
        
        # compute the final masks using the specified coef_selection_method
        self.sel_ = []
        sel_cls_flat = [
            reduce(np.union1d, [np.ravel_multi_index(ind, (self.n_chs_, self.n_freqs_, self.n_times_)) for ind in inds_cls])
            for inds_cls in self.sel_all_
        ] # find all unique significant indices for each class (flattened)

        if (self.coef_selection_method == 'union') or (self.coef_selection_method == 'ignore'):
            self.sel_ = [np.unravel_index(reduce(np.union1d, sel_cls_flat), (self.n_chs_, self.n_freqs_, self.n_times_))]
        elif self.coef_selection_method == 'intersect':
            self.sel_ = [np.unravel_index(reduce(np.intersect1d, sel_cls_flat), (self.n_chs_, self.n_freqs_, self.n_times_))]
        elif self.coef_selection_method == 'all':
            self.sel_ = [np.unravel_index(ind_flat, (self.n_chs_, self.n_freqs_, self.n_times_)) for ind_flat in sel_cls_flat]
        else:
            pass
        
        if self.matlab_eng is None:
            eng.quit() # stop auto-created matlab session for saving resources
        
        return self
    
    def transform(self, X):
        '''
        Applies the wavelet packet denoising to the input data.

        Parameters:
        - X: numpy array of shape (n_trials, n_channels, n_times)
            The input data.

        Returns:
        - X_denoised: numpy array of shape (n_trials, n_channels * n_classes, n_times)
            The denoised data. 
            n_classes=1 if coef_selection_method is one of ('union', 'ignore', 'intersect').
            n_classes=numer of unique layers in y during fitting if coef_selection_method='all'.
        '''
        eng = self.matlab_eng if self.matlab_eng is not None else matlab.engine.start_matlab()
        
        # increase index by one and ravel the indices in Fortran-like index order for matlab compatibility
        selectors = [
            matlab.int64(
                (
                    np.ravel_multi_index(indices, (self.n_chs_, self.n_freqs_, self.n_times_), order='F') + 1
                ).tolist()
            )
            for indices in self.sel_
        ]

        X_denoised = eng.modwptdenoise(
            matlab.double(X.tolist()), 
            selectors, 
            self.LoD_, 
            self.HiD_, 
            self.level_
        )
    
        if self.matlab_eng is None:
            eng.quit() # stop auto-created matlab session for saving resources
            
        return np.asarray(X_denoised)
    
    def transform_by_cluster(self, X, clusters='all'):
        '''
        Applies wavelet packet denoising to the input data based on specific clusters.

        Parameters:
        - X: numpy array of shape (n_trials, n_channels, n_times)
            The input data.
        - clusters: str, int, or list
            - 'all': Return denoised X for all significant clusters.
            - int: Return denoised X using the cluster at the given index.
            - list: Return denoised X using the clusters at the given indices.

        Returns:
        - X_denoised: numpy array or list of numpy arrays
            The denoised data.
        '''
        eng = self.matlab_eng if self.matlab_eng is not None else matlab.engine.start_matlab()
        
        if clusters == 'all':
            indices = np.arange(self.get_num_of_all_clusters())
        elif type(clusters) == int:
            indices = [clusters]
        elif isinstance(clusters, (list, tuple, np.ndarray)):
            indices = clusters
        
        X_denoised = []
        i = 0
        for sel in self.sel_all_:
            for clu in sel:
                if i in indices:
                    selectors = [matlab.int64(
                        (
                            np.ravel_multi_index(clu, (self.n_chs_, self.n_freqs_, self.n_times_), order='F') + 1
                        ).tolist()
                    )]
                    X_denoised.append(
                        np.asarray(
                            eng.modwptdenoise(matlab.double(X.tolist()), selectors, self.LoD_, self.HiD_, self.level_)
                        )
                    )
                i += 1
        
        if self.matlab_eng is None:
            eng.quit() # stop auto-created matlab session for saving resources
        
        if type(clusters) == int:
            return X_denoised[0]
        else:
            return X_denoised
        
    def get_num_of_all_clusters(self):
        '''
        Returns the total number of significant clusters across all classes.

        Returns:
        - num_clusters: int
            The total number of significant clusters.
        '''
        i = 0
        for sel in self.sel_all_:
            for clu in sel:
                i += 1
        return i
    
    def decompose(self, X, eng):
        # obtain wavelet packet decomposition
        return eng.modwptdec(matlab.double(X.tolist()), self.LoD_, self.HiD_, self.level_)


class NZTDenoiser(BaseEstimator, TransformerMixin):
    """
    NZT denoiser.

    Parameters:
    - stim: int, default 0
        The sample index where the stimuli occur in each epoch.
    - level: int, default 4
        The level of decomposition. Should be between 1 and min(4, log2(number of time points per trial)).
    - coef_selection_method: str, default 'all'
        The coefficient selection method. Used together with y_sel.
        - 'ignore': The coefficients are calculated without considering class labels.
        - 'all': Each class specified corresponds to one coefficient selection.
    - y_sel: list or None, default None
        The classes to be considered when calculating coefficient selectors.
    - matlab_eng: Matlab engine object or None, default None
        The Matlab engine object to use for MATLAB calculations. If None, a new engine will be created.
    """
    
    def __init__(
            self, 
            stim=0, 
            level=4, 
            coef_selection_method='all', 
            y_sel=None, 
            matlab_eng=None
    ):
        assert coef_selection_method in ['all', 'ignore'], \
            "coef_selection_method must be one of ('ignore', 'all')"
        
        self.level = level
        self.stim = stim
        self.coef_selection_method = coef_selection_method
        self.y_sel = y_sel
        self.matlab_eng = matlab_eng
        
    def fit(self, X, y=None):
        """
        Fit the NZTDenoiser model to the data.

        Parameters:
        - X: numpy array of shape (n_trials, n_channels, n_times)
            The input data.
        - y: numpy array or None, optional (default=None)
            The class labels. If None, all trials are treated together.

        Returns:
        - self: NZTDenoiser
            The fitted NZTDenoiser object.
        """
        eng = self.matlab_eng if self.matlab_eng is not None else matlab.engine.start_matlab()

        if y is None:
            y = np.zeros((1,X.shape[0]))

        classes = np.unique(y)
        self.n_times_ = int(X.shape[2])
        self.n_chs_ = int(X.shape[1])
        y_sel = self.y_sel if not self.y_sel is None else classes
        
        self.den_coeffs_ = []
        for c in classes:
            if not c in y_sel:
                continue
            X_class = X[y==c] if not self.coef_selection_method == 'ignore' else X
            
            den_coeff = eng.NZT_fit(
                matlab.double(X_class.tolist()), matlab.double(int(self.stim)), matlab.double(int(self.level))
            )
            self.den_coeffs_.append(np.array(den_coeff))
            
            if self.coef_selection_method == 'ignore':
                break
        
        if self.matlab_eng is None:
            eng.quit() # stop auto-created matlab session for saving resources
            
        return self
    
    def transform(self, X):
        """
        Apply the NZTDenoiser model to the data.

        Parameters:
        - X: numpy array of shape (n_trials, n_channels, n_times)
            The input data.

        Returns:
        - X_filt: numpy array of shape (n_trials, n_channels * n_classes, n_times)
            The filtered data.
            n_classes=1 if coef_selection_method='ignore'.
            n_classes=numer of unique layers in y during fitting if coef_selection_method='all'.
        """
        eng = self.matlab_eng if self.matlab_eng is not None else matlab.engine.start_matlab()
        
        X_filt = []
        for den_coeff in self.den_coeffs_:
            X_filt.append(
                np.array(eng.NZT_transform(
                    matlab.double(X.tolist()), [matlab.double(d) for d in den_coeff], matlab.double(self.level)
                ))
            )
        
        if self.matlab_eng is None:
            eng.quit() # stop auto-created matlab session for saving resources
            
        return np.concatenate(X_filt, axis=1)
    