# godunov2d_distributed_SOU
Solve 2d Euler equations using Godunov FVM. 
Utilized second-order upwind scheme with minmod flux limiter (ala MUSCL approach). 
Flux splliting is based on AUSMup or ROE methods. Time marching is possible by explicit first-order Euler or TVD RK3 schemes. 

The test case is  A Mach 3 Wind Tunnel With a Step.

For more than fifty years supersonic Mach 3 flow over forward step is a canonical benchmark to test different computation schemes. 
It was introduced for the first time by Emery, then was represented by Van Leer, and was popularized most widely by Woodward and Colella (1984).
Special interest is general capability of numerical methods to reproduce complex aerodynamic phenomena,  such as transient interactions of shock 
and rarefaction waves as well as Mach disks arising in the process of irregular interactions of waves with each other and with the channel walls. 
It is well known that the corner of the step  is a singularity point, which can produce  significant artificial numerical errors depend on 
applied boundary conditions and numerical schemes used. 

Figure below demonstrates the performance of the algorithm  to reproduce this challenging benchmark on the baseline 80x20 trinagular grid(s).
We compare two solutions obtained on the different grids: one with the sharp corner and another with the smoothed corner (a and b, respectively).
The implemented solver demonstrates the adequate performance to reproduce the the flow phenomena for the case with the smoothed step corner. 
Results obtained for the default configurtaion with the sharped edges revealed disturbed flow field near the Mach stem. 
This is due to numerical instabilities of entropy oscillations, which are generated  by the numerical schemes at the triple point.  

![alt text](https://github.com/dimaZloy/godunov2d_distributed_SOU/blob/main/results.png)
