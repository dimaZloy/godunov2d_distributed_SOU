
println("set numerics ...");


thermo = THERMOPHYSICS(287.0,1.4,1000.0,0.0,0.0);


solver = SOLVER2D(
 	1, # AUSM+ flux
	1, # 1- mimMod , 2 - van Leer
	0, # 1-FOU, 2-SOU 
	2  # 0-RK2; 2 - RK4
	);


solControls = CONTROLS(
	0.1, #CFL
	2.0e-6, # time step, 
	0, # fixed timeStepMethod (1 - adaptive)
	0.0,  # actual physical time to start simulation
	4.0/sqrt(1.4)/100.0,  # actual physical time to stop simulation 
	1, # flag to plot residuals
	0, # flag to constrain density
	0.01, # minDensityConstrained::Float64;
	10.0 # maxDensityConstrained::Float64;	
	);

pControls = plotCONTROLS(
	1, # 1 - filled contours , 0 - contour lines only
	10, # number of contours to plot 
	1.0, #min density
	2.8, #max density 
	0 # flag to production video 
	);

dynControls = DYNAMICCONTROLS(
	0.0, # actual physical time
	0.0, # cpu time;
	0.0, # tau	
	0, # iterator for verbosity (output)
	0, # global iterator for time steppings
	1.0, #max density in domain
	1.0, #min density in domain;
	0.0, #max velocity in domain
	"", # local path to the code
	"", #path to the test 	
	0, # flag to show if the convergence criteria is satisfied
	1  # flag to run or stop Simulation
);

dynControls.globalPath = pwd();
#dynControls.localTestPath = testdir;


output = outputCONTROLS(
	50, #verbosity::Int8;  
	"Time[s]\t Tau[s]\t Resid1\t Resid2\t Resid3\t Resid4\t CPUtime [s]", 
	0, #saveResiduals::Int8;
	0, #saveResults::Int8; 
	"residuals.dat",#fileNameResults::String;
	"solution.dat", #fileNameResiduals::String;
	0, ## save data to VTK
	"zzz"
);


flowTime = solControls.startTime;






