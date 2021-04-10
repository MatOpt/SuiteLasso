Welcome to SuiteLasso! 
Authors: Xudong Li, Defeng Sun, and Kim-Chuan Toh. 

A MATALB software package for solving regression problems with various
generalized Lasso regularizations based on
semismooth Newton augmented Lagrangian algorithms.

Set up
(1)unpack the software
(2)Run Matlab in the directory SuiteLasso
(3)In the Matlab command window, type: 
>> startup 
By now, SuiteLasso is ready for you to use.
=========================================================================================================================
Solver for classic lasso problems: ClassicLasso_SSNAL

Run files are provided for demonstration purpose: 
(a) test_ClassicLasso_random: for LASSO problems with randomly generated data 
(b) test_ClassicLasso_sparco: for LASSO problems with sparco data 
(c) test_ClassicLasso_UCI: for LASSO problems with UCI data 

(Note that the UCI dataset is not included in the package. To generate the UCI dataset, 
download the original datasets tested in above paper from UCI data repository 
and put them into the folder [\Util\genUCIdatafun\UCIdataorg]. 
Then, run genUCIdata.m. For demonstration purpose, 
''bodyfat_scale'' is left in the folder [\Util\genUCIdatafun\UCIdataorg].)

===========================================================================================================================
===========================================================================================================================
Solver for fused lasso problems: FusedLasso_SSNAL

Run file for demonstration purpose: 
test_Fused_Lasso 
===========================================================================================================================

Although our code was mainly designed for solving lasso and fused lasso in the high-dimensional setting,  
i.e., m (#of samples) << n (#of features), we provide wrapper files (Classic_Lasso_SSNAL_Wrapper and 
Fused_Lasso_SSNAL_Wrapper.m) to hand the case where m >> n 

You may want to try test_ClassicLasso_wrapper.m and test_FusedLasso_wrapper.m



