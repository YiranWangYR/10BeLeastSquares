# 10BeLeastSquares
This project includes codes to estimate surface exposure age from 10Be depth profiles, using least-squares linear inversion combined with a Monte Carlo simulation.
Detailed description of the codes can be found in Wang and Oskin, Combined linear regression and Monte Carlo approach to modelling exposure age depth profiles, 202x

These codes can be ran in matlab. 

To estimate surface age with known erosion rate, run "Be10_LS_rate.m"
To estimate surface age with known eroded thickness, run "Be10_LS_thickness.m", this code runs together with the function "Be10Newton", so the two files should be placed in the same folder.
Notice that both approaches assume constant erosion.

If there is no surface erosion, both "Be10_LS_rate.m" and "Be10_LS_thickness.m" works the same.

We do not exclude negative inheritance, in order to aquire the full distribution of the exposure age (See detailed discussion in Wang and Oskin 202x). 

These codes were written with an old matlab version (2014), therefore some lines of the codes maybe simplified if using a newer version.
