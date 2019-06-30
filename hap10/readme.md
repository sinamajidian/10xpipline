


# How to use the MATLAB code of hap10
 

## Installation on LINUX
1. Download the MATLAB Runtime R2019a (9.6) for free from [here](http://ssd.mathworks.com/supportfiles/downloads/R2019a/Release/3/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019a_Update_3_glnxa64.zip)

2- Unzip it using `unzip MATLAB_Runtime_R2019a_Update_3_glnxa64.zip`

3- Install it using `./install`. During the installation you are asked to set the destination directory `<MR_directory>`.

4- Download the hap10 package from [here](https://raw.githubusercontent.com/smajidian/10xpipline/master/hap10/hap10_linux) 


5- Run the hap10 using

```
cd hap10_linux
./hap10.sh <MR_directory>/v96 /pathto/fragment_file.txt k
```
in which k is the ploidy level.


## Installation on  MAC
1. Download the MATLAB Runtime R2018b (9.5) for free from [here](http://ssd.mathworks.com/supportfiles/downloads/R2018b/deployment_files/R2018b/installers/maci64/MCR_R2018b_maci64_installer.dmg.zip)

2- Unzip it using `MCR_R2018b_maci64_installer.dmg.zip`

3- Install it using `./install`. During the installation you are asked to set the destination directory `<MR_directory>`.

4- Download the hap10 package from [here](https://github.com/smajidian/10xpipline/tree/master/hap10/hap10_mac) 

5- Run the hap10 using

```
cd hap10_mac
./hap10.sh <MR_directory>/v95 /pathto/fragment_file.txt k
```
in which k is the ploidy level.











## Compilation
This step is not needed for running haplotype. But If you want to compile from source code you need MATLAB and also MATLAB Compiler.  

Run the following on MAC
```
mkdir hap10

mcc -m hap10.m -d hap10 -a AXfun.p frag2mat.m ops.p Atyfun.p printinfo.p Fnorm.p projSDP.p ProjPS.p matvecAAt.p psqmrNEW.p admmplus.m mec_calculator.m refiner.m  mexADM_rescale.mexmaci64 scaling.p admmplus_main_default.p mexFnorm.mexmaci64 sdp_solver.m blkprojSDP.p mexMatvec.mexmaci64 sdpnalplus.m blktrace.p mexNAL_ADMsigma_update.mexmaci64  sdpnalplus_main_default.p competaK1C1.p mexNAL_rescale.mexmaci64 smat.p competaK2C2.p mexeig.mexmaci64 solving_nal.m competaorg.p mexsmat.mexmaci64 svec.p convertdata.p mexsvec.mexmaci64 validate.p
```



The optimization core is from [SDPNAL+](http://www.math.nus.edu.sg/~mattohkc/SDPNALplus.html)
```
L.Q. Yang, D.F. Sun, and K.C. Toh, SDPNAL+: a majorized semismooth Newton-CG augmented Lagrangian method for semidefinite programming with nonnegative constraints, Mathemtical Programming Computation, 7 (2015), pp. 331-366. arXiv:1406.0942.
```

## Copyright
This version of SDPNAL+ is distributed under the Creative Commons Attribution-ShareAlike 4.0 International Public License.
