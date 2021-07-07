

using Distributed;
using PyPlot;

const pname = "testTriMesh2d.bson"
const numThreads = 8;


if (numThreads != 1)

	if (nprocs() == 1)
		addprocs(numThreads,lazy=false); 
		display(workers());
	end
		
end


@everywhere using PyPlot;
@everywhere using WriteVTK;
@everywhere using CPUTime;
@everywhere using DelimitedFiles;
@everywhere using Printf
@everywhere using BSON: @load
@everywhere using BSON: @save
@everywhere using SharedArrays;

include("primeObjects.jl");
include("thermo.jl"); #setup thermodynamics
include("utilsIO.jl");
include("RoeFlux2d.jl")
include("AUSMflux2d.jl"); #AUSM+ inviscid flux calculation 
include("utilsFVM2dp.jl"); #FVM utililities
## utilsFVM2dp::cells2nodesSolutionReconstructionWithStencilsImplicitSA
## utilsFVM2dp::cells2nodesSolutionReconstructionWithStencilsSA
## utilsFVM2dp::phs2dcns2dcellsSA

include("boundaryConditions2d.jl"); 

include("initfields2d.jl");
## initfields2d::distibuteCellsInThreadsSA()
## initfields2d::createFields2d_shared()

include("evaluate2d.jl"); 
## propagate2d::updateResidualSA()
## propagate2d::updateVariablesSA()
## propagate2d::updateOutputSA()

include("limiters.jl");
include("computeslope2d.jl");
include("SOUscheme.jl");
## computeslope2d:: computeInterfaceSlope()
## SOUscheme:: SecondOrderUpwindM2()

include("propagate2d.jl");
## propagate:: calcOneStage() expilict Euler first order
## propagate:: doExplicitRK3TVD() expilict RK3-TVD


function saveResults2VTK(filename::String, 
		testMesh::mesh2d,field::Array{Float64,1})

	vtkfile = vtk_grid(filename, testMesh.xNodes, testMesh.yNodes, testMesh.VTKCells);
	vtk_point_data(vtkfile, field, "density");
	outfiles = vtk_save(vtkfile);
	
end






function godunov2dthreads(pname::String, numThreads::Int64, outputfile::String)

	
	@load "testStep2dBaseTriSmooth.bson" testMesh
	#@load "testStep2dBaseTri.bson" testMesh ## load mesh2d
	
	## only for tri meshes!!!
	triangles = zeros(Int64,testMesh.nCells,3);
	for i = 1:testMesh.nCells
	
		## indexes of nodes in PyPLOT are started from Zero!!!!
	
		triangles[i,1] = testMesh.mesh_connectivity[i,4]-1;
		triangles[i,2] = testMesh.mesh_connectivity[i,5]-1;
		triangles[i,3] = testMesh.mesh_connectivity[i,6]-1;
	end
	
	cellsThreads = distibuteCellsInThreadsSA(numThreads, testMesh.nCells); ## partition mesh 
	

	include("setupSolver2d.jl"); #setup FVM and numerical schemes
	
	
	## init primitive variables 
	println("set initial and boundary conditions ...");
	testfields2d = createFields2d_shared(testMesh, thermo);
	
	println("nCells:\t", testMesh.nCells);
	println("nNodes:\t", testMesh.nNodes);
	
	## init conservative variables 
	UconsCellsOldX = SharedArray{Float64}(testMesh.nCells,4);
	UconsCellsNewX = SharedArray{Float64}(testMesh.nCells,4);
	
	UconsCellsNew1X = SharedArray{Float64}(testMesh.nCells,4);
	UconsCellsNew2X = SharedArray{Float64}(testMesh.nCells,4);
	UconsCellsNew3X = SharedArray{Float64}(testMesh.nCells,4);
	
	
	DeltaX = SharedArray{Float64}(testMesh.nCells,4);
	iFLUXX  = SharedArray{Float64}(testMesh.nCells,4);
	#iFLUXdist  = SharedArray{Float64}(testMesh.nCells,4);
	
	
	phs2dcns2dcellsSA(UconsCellsOldX,testfields2d, thermo.Gamma);	
	phs2dcns2dcellsSA(UconsCellsNewX,testfields2d, thermo.Gamma);	
	
	
	
	testMeshDistr = createMesh2dShared(testMesh);
	
		
	println("Start calculations ...");
	println(output.header);
	
	@everywhere trianglesX = $triangles;
	@everywhere testMeshDistrX = $testMeshDistr; 
	
	#@everywhere testMeshDistrX = $testMesh_shared;
	@everywhere testMeshX = $testMesh; 
	@everywhere thermoX   = $thermo;
	@everywhere cellsThreadsX = $cellsThreads;
	@everywhere testfields2dX  = $testfields2d;
		
	@everywhere dynControlsX = $dynControls;
	@everywhere solControlsX = $solControls;
	@everywhere pControlsX = $pControls;
	@everywhere outputX = $output;
	


	timeVector = [];
	residualsVector1 = []; 
	residualsVector2 = []; 
	residualsVector3 = []; 
	residualsVector4 = []; 
	residualsVectorMax = ones(Float64,4);
	convergenceCriteria= [1e-5;1e-5;1e-5;1e-5;];
	
	
	densityF = zeros(Float64,testMesh.nNodes);
	for i = 1:testMesh.nNodes
		densityF[i] = testfields2dX.densityNodes[i];
	end
				
	
		 
	println("Saving  solution to  ", outputfile);
		saveResults2VTK(outputfile, testMesh, densityF);
	println("done ...  ");	
	
	
	dt::Float64 =  solControls.dt;  
	@everywhere dtX = $dt; 
	
	
		#for h = 1:2
		while (dynControlsX.isRunSimulation == 1)
			CPUtic();	
			
			
			# PROPAGATE STAGE: 
			(dynControlsX.velmax,id) = findmax(testfields2dX.VMAXCells);
			# #dynControls.tau = solControls.CFL * testMesh.maxEdgeLength/(max(dynControls.velmax,1.0e-6)); !!!!
			dynControlsX.tau = solControlsX.CFL * testMeshX.maxArea/(max(dynControlsX.velmax,1.0e-6));
		
	
				
			## Explicit Euler first-order	
			## calcOneStage(1.0, dtX, testMeshDistrX , testfields2dX, thermoX, cellsThreadsX,  UconsCellsOldX, iFLUXX,  UconsCellsNewX);
			
			doExplicitRK3TVD(1.0, dtX, testMeshDistrX , testfields2dX, thermoX, cellsThreadsX,  UconsCellsOldX, iFLUXX,  
			UconsCellsNew1X,UconsCellsNew2X,UconsCellsNew3X,UconsCellsNewX);
			
			@sync @distributed for p in workers()	
	
				beginCell::Int64 = cellsThreadsX[p-1,1];
				endCell::Int64 = cellsThreadsX[p-1,2];
				#println("worker: ",p,"\tbegin cell: ",beginCell,"\tend cell: ", endCell);
														
				updateVariablesSA(beginCell, endCell, thermoX.Gamma,  UconsCellsNewX, UconsCellsOldX, DeltaX, testfields2dX);
		
			end
			
			@everywhere finalize(updateVariablesSA);	
								
		
			cells2nodesSolutionReconstructionWithStencilsImplicitSA(testMeshX, testfields2dX); 
			
	
			(dynControlsX.rhoMax,id) = findmax(testfields2dX.densityCells);
			(dynControlsX.rhoMin,id) = findmin(testfields2dX.densityCells);
			

			push!(timeVector, dynControlsX.flowTime); 
			dynControlsX.curIter += 1; 
			dynControlsX.verIter += 1;
			
			
			updateResidualSA(DeltaX, 
				residualsVector1,residualsVector2,residualsVector3,residualsVector4, residualsVectorMax,  
				convergenceCriteria, dynControlsX);
			
			updateOutputSA(timeVector,residualsVector1,residualsVector2,residualsVector3,residualsVector4, residualsVectorMax, 
				testMeshX, trianglesX, testfields2dX, solControlsX, outputX, dynControlsX);
	
			
			# EVALUATE STAGE:
			
			dynControls.flowTime += dt; 
			
			# if (solControlsX.timeStepMethod == 1)
				# dynControlsX.flowTime += dynControlsX.tau;  	
			# else
				# dynControlsX.flowTime += solControlsX.dt;  
			# end
			

	

			if (flowTime>= solControlsX.stopTime || dynControlsX.isSolutionConverged == 1)
				dynControlsX.isRunSimulation = 0;
		
				if (dynControlsX.isSolutionConverged == true)
					println("Solution converged! ");
				else
					println("Simultaion flow time reached the set Time!");
				end
			
				if (outputX.saveResiduals == 1)
					#println("Saving Residuals ... ");
					#cd(dynControlsX.localTestPath);
					#saveResiduals(output.fileNameResiduals, timeVector, residualsVector1, residualsVector2, residualsVector3, residualsVector4);
					#cd(dynControlsX.globalPath);
				end
				if (outputX.saveResults == 1)
					#println("Saving Results ... ");
					#cd(dynControlsX.localTestPath);
					#saveSolution(output.fileNameResults, testMeshX.xNodes, testMeshX.yNodes, UphysNodes);
					#cd(dynControlsX.globalPath);
				end
			
			end

			dynControlsX.cpuTime  += CPUtoq(); 
			
			if (dynControlsX.flowTime >= solControls.stopTime)
				dynControlsX.isRunSimulation = 0;
			end
			
		end ## end while
		 
		 
		for i = 1:testMesh.nNodes
			densityF[i] = testfields2dX.densityNodes[i];
		end
			
		 
		println("Saving  solution to  ", outputfile);
		saveResults2VTK(outputfile, testMesh, densityF);
		println("done ...  ");	
	
end


#@time godunov2dthreads("testMixedMesh2d.bson",4); 
@time godunov2dthreads("testTriMesh2d.bson",numThreads, "souTriEuler"); 
#@time godunov2dthreads("testQuadMesh2d.bson",4); 


#@time godunov2d("testQuadMesh2d.bson"); 
#@time godunov2d("testMixedMesh2d.bson"); 





