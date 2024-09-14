# Supplementary code for data synthesis and denoising

This package contains supplementary code for our work following [1] (details will be available after been published). The package includes all necessary code and data for synthesising EEG trials (with or without error-related potentials - ErrPs) and performing single trial denoising. Both the novel threshold-free wavelet packet transform-based denoising (TF-WPTD) method (proposed) and the so-called NZT method (baseline) [2] are included.

This package also provides a notebook (`Simulation_and_denoising_example.ipynb`), performing denoising steps and illustrating denoising results. 

The EEG data synthesis was based on the real human EEG recorded in our previous work [1]. The data synthesis settings and the TF-WPTD method are briefly described in `Methods.pdf`.

# How to run the sample notebook

To run the sample notebook, ensure that MATLAB and Python are installed on your machine, and ensure the "Wavelet Toolbox" and "Signal Processing Toolbox" are installed for MATLAB. Follow the steps below to install the necessary Python libraries. Also, edit the MATLAB path to include all MATLAB scripts in this package.

## Install

Create and enter a python virtual environment with python>=3.8:

```bash
conda create -n myenv python=3.8
conda activate myenv
```

Install required python libraries using pip (if your environment already supports jupyter notebook, you may delete "notebook" from requirement.txt):

```bash
pip install -r requirements.txt
```

Install matlab engine separately by following MathWorks' official documentation based on your MATLAB version:
https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

In MATLAB, add the "src_matlab" folder and its subfolders to the path. Also add the "data" folder to the path.

## Run sample notebook
Open `Simulation_and_denoising_example.ipynb` with jupyter notebook, jupyter lab, visual studio code, or equivalent. Then run all cells from top to bottom.

# Folder and file structure
```
│   requirements.txt
│   Simulation_and_denoising_example.ipynb
│
├───data
│       estimated_psd.mat
│       gen_noise.m
│       simulate_EEG.m
│       spatial_mapping.mat
│
├───src_matlab
│   │   modwptdec.m
│   │   modwptdenoise.m
│   │   modwptrec.m
│   │   NZT_fit.m
│   │   NZT_transform.m
│   │   rmodwptphaseshift.m
│   │
│   ├───EP_den_Auto
│   └───simulator
│
└───src_python
        denoisers.py
```

> **`./requirements.txt`** <p>

The required python libraries to run the sample notebook (apart from the MATLAB engine to be installed separately).

> **`./Simulation_and_denoising_example.ipynb`** <p>

A sample notebook with all steps required to recreate figure 6 in the paper.

> **`./data/`*** <p>

The code, along with the spatial distributions of the required ERP components and the power spectral density measurements obtained from real EEG data, for data synthesis. Scripts under this folder used functions implemented in "./src_matlab/simulator". Running the notebook will save the synthesized EEG trials in this folder as well.

> **`./src_matlab/modwptdec.m`** <br>
**`./src_matlab/modwptdenoise.m`** <br> 
**`./src_matlab/modwptrec.m`** <br> 
**`./src_matlab/rmodwptphaseshift.m`** <p>

MATLAB functions to support executing the TF-WPTD method.

> **`./src_matlab/NZT_fit.m`** <br>
**`./src_matlab/NZT_transform.m`** <p> 

MATLAB functions to support executing the NZT method. Scripts under this folder used functions implemented in "./src_matlab/EP_den_Auto".

> **`./src_matlab/EP_den_Auto/*`** <p>

The implementation of the NZT denoising method. Cloned from https://github.com/LightingResearchCenter/oddball_threeStimulus/tree/7ba8566ae5b1bcafa69c8ea199f53aa6ab883638/matlab_threeStimulus/EP_den_auto

> **`./src_matlab/simulator/*`** <p>

The implementation used to generate simulated EEG data. Cloned from https://github.com/pchrapka/phasereset/tree/809996cbbc63a6af71528d5b5bfa5c13765947aa

# Acknowledgement

We would like to thank https://github.com/pchrapka for sharing the MATLAB code implementing the EEG data synthesis method as described in [1].

We would also like to thank https://github.com/LightingResearchCenter for sharing the MATLAB code implementing the NZT method as described in [2].

# Refereces

[1] Tang, Y., Zhang, J. J., Corballis, P. M., & Hallum, L. E. (2021, November). Towards the 
classification of error-related potentials using riemannian geometry. In _2021 43rd Annual International 
Conference of the IEEE Engineering in Medicine & Biology Society (EMBC)_ (pp. 5905-5908). IEEE.

[2] Ahmadi, M., & Quiroga, R. Q. (2013). Automatic denoising of single-trial evoked potentials. 
_NeuroImage_, 66, 672-680.
