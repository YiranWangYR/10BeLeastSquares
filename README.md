# 10BeLeastSquares
This project includes codes to estimate surface exposure age from 10Be depth profiles, using least-squares linear inversion combined with a Monte Carlo simulation.

Detailed description of the codes can be found in 
Wang and Oskin, Combined linear regression and Monte Carlo approach to modelling exposure age depth profiles, 202x

To estimate surface age with known erosion/denudation rate, run "Be10_LS_rate.m", compatible with matlab2014

To estimate surface age with known eroded thickness/denudation length, run "Be10_LS_thickness.m", or "Be10_LS_thickness_v2" 
"Be10_LS_thickness.m" is compatible with matlab2014, it runs together with the function "Be10Newton", so the two files should be placed in the same folder.
"Be10_LS_thickness_v2" is compatible with matlab 2020, the function "Be10Newton" is already integrated in the script.

Notice that both approaches assume constant erosion.

If there is no surface erosion, both "Be10_LS_rate.m" and "Be10_LS_thickness.m" works the same.

We do not exclude negative inheritance, in order to aquire the full distribution of the exposure age (See detailed discussion in Wang and Oskin 202x). 

If you encounter any problem, please feel free to contact me (yrwwang@ucdavis.edu)
