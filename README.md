MQ-CEST Fit
------------
------------
This code is designed for fitting multi-quantum chemical exchange saturation
transfer data (MQ-CEST) see Karunanithy et al., J Phys. Chem. Lett., 2020
(https://doi.org/10.1021/acs.jpclett.0c01322).
It is currently designed for the characterisation of restricted rotation in the
guanidinium groups of arginine residues and can give information on the rate of
rotation of Nh nuclei and their individual chemical shifts.

Dependencies
------------
  * [Python=3.5](https://www.python.org/downloads/)
  * [SciPy=1.4.1](https://www.scipy.org/install.html)
  * [NumPy=1.18.1](https://www.scipy.org/scipylib/download.html)
  * [Matplotlib=3.0.3](http://matplotlib.org/users/installing.html)
  * [LmFit=1.0.0](https://lmfit.github.io/lmfit-py/)

The MQ-CEST fitting program has been written and tested in Python3 with the
above dependencies. Performance with other module versions has not been tested.

Running MQ-CEST Fit
------------
The main scripts (ABX.py, ABX_reducedBasis.py, red_run_mqfit.py,run_mqfit.py)
should be placed in the system PATH. The program can then be initiated with the
command:

mqcestFit.py runscript

Where runscript refers to the main input file (described below). To run fits of
the example data run:

mqcestFit.py example.in

from the folder containing 'example.in'. The 'example_data' and 'example_infiles'
folder should also be in this directory.

Main Input Script
------------
This is the main input file required to initiate the fits. It contains the
following parameters with input arguments entered followed by a space:

* infol        : this is a folder that contains all the residue specific input
              files (see below for details of these)

* outfol       : this is the directory that will be created for all results from
              the fits. Within this directory one directory will be created
              for the results from each residue.

* resiList     : this is a list of all residues to be fit. For each residue there
              must be a corresponding residue specific input file (with the
              residue name given here and filename extension '.in') in 'infol'

* fitParams    : this is a list of non-relaxation parameters to include in the
              fit. Typically this will be kex, deltaO and w0 but it is possible
              to fix/fit other parameters as desired. Any parameters not listed
                here will be fixed in the optimization process.

* fitParams_rel : A list of relaxation parameters to include in the fits. options
                are r_Cz (longitudinal carbon relaxation), r_Nz (longitudinal
                nitrogen relaxation) and r_Nxy (transverse nitrogen relaxation).
                Note that if data has been collected at multiple B0 fields a
                separate rate for each field is created.

Default System Parameters

Here we setup some default system parameters. These will be overwritten on a per
residue basis where a value is provided in the residue specific input files:

* r_Cz    : carbon longitudinal relaxation rate

* r_Nz    : nitrogen longitudinal relaxation rate

* r_Nxy   : nitrogen transverse relaxation rate

* J       : Cz-Nh coupling constant

* kex     : rotation rate about Cz-Ne bond

* pb      : population of state B (should fix to 0.5 for guanidinium rotation)

* deltaO  : difference in Nh chemical shifts

Default Experimental Parameters

These parameters cannot be fit. The default parameters here can however be
overwritten by values provided in the residue specific input files:

* cest_time : length of CEST saturation pulse in seconds

* inhom_num : number of B1 values to use to characterise B1 inhomogeneity

* phase     : phase of CEST pulse (0.0 = x; 90.0 = y; 180.0 = -x; 240.0 = -y)

Please see example.in as a reference.

Residue Specific Input File
------------
For each residue listed in the main input file there must be a corresponding
residue specific input file located in the 'infol' with the residue name as
listed in 'resiList' and the filename extension '.in'. Parameters set here will
overwrite those set in the main input script.

The residue specific input files have the following input parameters:

* resi      : the residue name

* DataFiles : A list of locations for the relevant CEST data for this residue.
              The order here is important (see below).

* J         : Cz-Nh coupling constant

* pb        : population of state B (should fix to 0.5 for guanidinium rotation)

* kex       : rotation rate about Cz-Ne bond

* w0        : this is the average of the two Nh chemical shifts (in ppm)

* deltaO    : difference between Nh chemical shifts (in ppm)

* cest_time : length of CEST saturation pulse in seconds

* phase     : phase of CEST pulse (0.0 = x; 90.0 = y; 180.0 = -x; 240.0 = -y)

* inhom     : number of B1 values to use to characterise B1 inhomogeneity

* field     : list of 15N B0 fields in MHz. One field is required per dataset.
              The order of these should line up with the DataFiles list.

* v1        : list of cest B1 fields in Hz. One field is required per dataset.
              The order of these should line up with the DataFiles list.

* v1_std   : list of cest B1 field inhomogeneity in Hz. One value is required
             per dataset. The order of these should line up with the DataFiles
             list.

* carrier  : list of 15N carrier frequencies in ppm. The CEST offsets in the
             dataset files should be given in Hz relative to this. The order of
             these should line up with the DataFiles list.

* r_Cz     : list of carbon longitudinal relaxation rates. One rate is required
             per dataset file.

* r_Nz     : list of nitrogen longitudinal relaxation rates. One rate is
             required per dataset file.

* r_Nxy    : list of nitrogen transverse relaxation rates. One rate is required
             per dataset file.

See example_infiles for residue specific input file references.

CEST data files
------------
The CEST data files should have three columns: offset of saturation pulse in Hz,
(this should be relative to the 15N carrier frequency)  normalised CEST
intensity (I/I0 this is the ratio of intensity in the presence of a saturation
pulse at the stated saturation offset divided by the intensity in the
absence of this pulse) and error (the standard deviation of the intensity
ratio). Any comments should be in lines beginning with '#'.

See example_data for CEST data files.

Outputs
------------
For each dataset file used in the fit an eps figure is created showing the
experimental and fitted CEST curves. A corresponding text file (.out) is created
with the values used to make the figure. Finally the fit report from LMFIT is
given showing the values of the fitted parameters and their errors.
