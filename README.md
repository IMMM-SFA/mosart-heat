# MOSART-Heat
A large-scale stream temperature model

## Current Version
This is where I will put your DOI tag once all has been finalized
Note that the version including the representation of water management is associated with the 2013 version of the WM code available at https://github.com/IMMM-SFA/wm and associated with the following DOI https://doi.org/10.5281/zenodo.1225344. 

## Notice
This repository uses the Git Large File Storage (LFS) extension (see https://git-lfs.github.com/ for details).  Please download GitLFS and install then run the following command before cloning if you do not already have Git LFS installed:
`git lfs install`.

## Contact
For questions please contact:

Hongyi Li: hongyili.jadison@gmail.com

## Overview
MOSART was developed as a scalable framework for representing and studying riverine dynamics of water, energy, and biogeochemistry cycles across local, regional and global scales from an integrated human-earth system perspective (Li et al., 2013). In particular, a dam/reservoir module has been embedded within MOSART (denoted as MOSART-WM) to simulate the effects of dam operations as well as surface-water withdrawal on downstream flow and water temperature across a range of spatial scales (Voisin et al. 2013; also available at https://github.com/IMMM-SFA/wm). Here, an energy-balance based stream temperature module was added on top of MOSART only (denoted as MOSART-heat) for the simulation of stream temperature in natural river systems and on top of MOSART-wm (denoted as MOSART-wm-heat) for the simulation of stream temperature in managed river systems. For more details about MOSART-wm-heat, please refer to Li et al. (2015).

## Getting Started with MOSART-Heat
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. The major purpose is for you to compile and run CLM4-MOSART-heat (stream discharge and temperature modeling simulation under natural conditions) and CLM4-MOSART-wm-heat (stream discharge and temperature modeling under water management) for the historical period (1980-2004) and future period (2005-2099) respectively. For the future period, both RCP4.5 and RCP8.5 scenarios are supported. So in total there could be up to two historical simulations, and four future simulations. Here we provide the input data to support simulations over the contiguous United States at a 1/8th-degree resolution and hourly time step.

### Prerequisites
The current versions of MOSART-heat and MOSART-wm-heat have been developed and used as part of the Community Land Model, version 4.0 (CLM4), which is the land component of Community Climate System Model (CCSM). Therefore, they can be compiled and run as part of CLM4 and do not need additional libraries besides those included in the CLM4 package. The users are expected to be sufficiently familiar with CCSM/CLM4.

#### PreStep 1: Configure CLM4 on your local machine
Please follow through the instruction provided at https://github.com/IMMM-SFA/PRIMA_CLM4. Note that in the shell script for setting the historical simulations (e.g., setup_prima_clm4_hist.sh), please modify the begin/end years to 1980/2004 due to the availability of our water demand data. For the future simulation, the begin/end years will be 2005-2095. You may customize the years if having your own water demand data. For each of your simulations, please provide a unique name and run the shell script to generate a corresponding case. After running the shell script, you should be able to see that CLM4/CCSM has been compiled successfully.

#### PreStep 2: Test run CLM4
A successful compiling does not warrant running CLM4 smoothly. Please do a test run to ensure that CLM4 runs smoothly. Taking the PRIMA CLM4 historical simulation for example,
```
cd $CASE_DIR/your_own_case
sbatch your_own_case.runs
```
If there is any error detected at this step, please work with your own technical support or contact Maoyi Huang (maoyi.huang@pnnl.gov) for help. Please DO NOT proceed until your CLM4 test run here can be completed successfully.

### Running MOSART-Heat and MOSART-wm-heat

#### Step 1: Download MOSART-heat and MOSART-wm-heat input data
The `zhang_et_al_2018_data.zip` file in the `data` directory of this repository should be unzipped and contains the `MOSART_wm_para` and `GCAM_waterdemand_nc` directories which house the MOSART and WM parameter files and the water demand data, respectively.  
These files include:

-	`MOSART_h_NLDAS_8th.nc`: MOSART routing parameters
-	`US_reservoir_8th_hist.nc`: water management parameters for historical simulations
-	`US_reservoir_8th_rcp8.5.nc`: water management parameters for RCP8.5 simulations
-	`US_reservoir_8th_rcp4.5.nc`: water management parameters for RCP4.5 simulations
-	`RCP4.5_GCAM_water_demand*.nc`: water demand data for RCP4.5 MOSART-wm-heat simulation (note the files dated as 1980-2004 can be used for the historical MOSART-wm-heat simulation also)
-	`RCP8.5_GCAM_water_demand*.nc`: water demand data for RCP8.5 MOSART-wm-heat simulation

Please note that those reservoir and water demands input parameter files were created and used in Hejazi et al. (2015) and any further use should be properly cited. "We used the MOSART and WM set ups developed and used in Hejazi et al. (2015). Readers can find more details in Li et al. (2013) for MOSART and Voisin et al. (2013) for WM on how to derive those parameters. 

#### Step 2: Obtain and use MOSART-heat and MOSART-wm-heat source code
Please download and use the MOSART-heat source code for your stream temperature simulations under natural conditions, and the MOSART-wm-heat source code for your stream temperature simulations under water management. Once downloaded, copy the source code to your own CLM4 case directory (e.g., $CASE_DIR/clm4_nldas_hist) into the subfolder SourceMod/src.clm (i.e., $CASE_DIR/clm4_nldas_hist/ SourceMods/src.clm).

#### Step 3: Modify user_nl_clm
In your own CLM4 case directory, modify the content of `user_nl_clm` file as following:
```
&clm_inparm
dtime       = 1800
frivinp_rtm		= '/pic/projects/prima/liho745/inputdata/MOSART_h_NLDAS_8th.nc'
rtm_nsteps             = 6
/
```

#### Step 4: Modify RtmMod.F90
This step is only needed for running MOSART-wm-heat (skip it if you are running MOSART-heat).

Modify Lines 1173-1175 in the `rtmMod.F90` file from `$CASE_DIR/clm4_nldas_hist/ SourceMods/src.clm`, to direct to your local paths for the MOSART-wm-heat parameter file(s) and water demand data.
```
WMctl%paraPath = '/pic/projects/prima/liho745/inputdata/MOSART_wm_parameters/'
WMctl%paraFile = 'US_reservoir_8th_hist.nc'
WMctl%demandPath = '/pic/projects/prima/liho745/inputdata/GCAM_waterdemand/RCP_nc/rcp4.5/RCP4.5_GCAM_water_demand_'
```

#### Step 5: Recompile and run MOSART-heat or MOSART-wm-heat
Use the following commands to recompile and run the MOSART-heat or MOSART-w-heat code:
```
cd $CASE_DIR/your_own_case
your_own_case.clean_build
your_own_case.build
sbatch your_own_case.run
```

## Recommended acknowledgement for using the above code package:

The authors would like to acknowledge H. Li at University of Houston for providing MOSART-heat code, M. Huang at Pacific Northwest National Laboratory (PNNL) for sharing the PRIMA CLM4 compsets and driving scripts, Nathalie Voisin at PNNL for providing the water management code , and Mohamad Hejazi at PNNL for providing the water demand data, all supported by the Platform for Regional Integrated Modeling and Analysis (PRIMA) Initiative and the U.S. Department of Energy, Office of Science as part of research in Multi-Sector Dynamics, Earth and Environmental System Modeling Program.

## References
Li, H., M. S. Wigmosta, H. Wu, M. Huang, Y. Ke, A. M. Coleman, and L. R. Leung (2013), A physically based runoff routing model for land surface and earth system models, J. of Hydromet., 14(3):808-828. doi:10.1175/JHM-D-12-015.1

Li, H., L. R. Leung, A Getirana, M Huang, H Wu, Y Xu, J Guo and N Voisin (2015), Evaluating Global Streamflow Simulations by a Physically-based Routing Model Coupled with the Community Land Model, J. of Hydromet., 16(2):948-971, doi: 10.1175/JHM-D-14-0079.1

Li, H., L. Ruby Leung, T. Tesfa, N. Voisin, M. Hejazi, L. Liu, Y. Liu, J. Rice, H. Wu, and X. Yang (2015), Modeling stream temperature in the Anthropocene: An earth system modeling approach, J. Adv. Model. Earth Syst., 7, doi:10.1002/2015MS000471.

Voisin, N., Li, H., Ward, D., Huang, M., Wigmosta, M., and Leung, L. R., 2013: On an improved sub-regional water resources management representation for integration into earth system models, Hydrol. Earth Syst. Sci., 17, 3605-3622, doi:10.5194/hess-17-3605-2013, 2013

Hejazi MI, N Voisin, L Liu, LM Bramer, DC Fortin, JE Hathaway, M Huang, GP Kyle, LYR Leung, H Li, Y Liu, PL Patel, TC Pulsipher, JS Rice, TK Tesfa, CR Vernon, and Y Zhou.  2015.  "21st Century United States Emissions Mitigation Could Increase Water Stress more than the Climate Change it is Mitigating."  Proceedings of the National Academy of Sciences of the United States of America 112(34):10635-10640.  doi:10.1073/pnas.1421675112

WM code (2013). https://doi.org/10.5281/zenodo.1225344
