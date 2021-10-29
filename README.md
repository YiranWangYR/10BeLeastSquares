# 10BeLeastSquares
This project includes codes to estimate surface exposure age from 10Be depth profiles, using least-squares linear inversion combined with a Monte Carlo simulation.
Detailed description of the codes can be found in Wang and Oskin, Combined linear regression and Monte Carlo approach to modelling exposure age depth profiles, 202x

These codes can be ran in matlab. 

Be10_LS_rate: the code to estimate surface age with known erosion rate
Be10_LS_thickness: the code to estimate surface age with known eroded thickness. This code runs together with the function "Be10Newton", so the two files should be placed in the same folder.
Be10Newton: the function that runs Newton's method. Runs together with "Be10_LS_thickness".

We do not exclude negative inheritance, in order to aquire the full distribution of the exposure age (See detailed discussion in Wang and Oskin 202x). 

These codes were written with an old matlab version (2014), therefore some lines of the codes maybe simplified if using a newer version.
