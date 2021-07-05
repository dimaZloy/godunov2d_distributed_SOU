

@everywhere  function computeInterfaceSlope(i::Int64, k::Int64, testMesh::mesh2d_shared, testFields::fields2d_shared, thermo::THERMOPHYSICS,  uLeftp::Array{Float64,1} ):: Array{Float64,1}

	nCells = size(testMesh.cell_stiffness,1);
	
	ek::Int64 = testMesh.cell_stiffness[i,k]; ##; %% get right cell 
	
	ek_type::Int64 = testMesh.mesh_connectivity[i,2];
	
	side::Float64 = testMesh.cell_edges_length[i,k];
	nx::Float64   = testMesh.cell_edges_Nx[i,k];
	ny::Float64   = testMesh.cell_edges_Ny[i,k];
				
	
	uUpp = zeros(Float64,4);
	uDownp = zeros(Float64,4);
	uRightp = zeros(Float64,4);

		
	index::Int64 = 0;
	if (k == 1)
		index = 1;
	elseif (k == 2)
		index = 3;
	elseif (k == 3)
		index = 5;
	elseif (k == 4)
		index = 7;	
	end
				
	pDown1::Int64 = 0;
	pDown2::Int64 = 0;
	
	pUp1::Int64 = 0;
	pUp2::Int64 = 0;
	
	
	if (ek >=1 && ek<=nCells)
								   
								   
		if (ek_type == 3) ## tri element 
		
			pDown1 = testMesh.node2cellsL2down[i,index];
			pUp1 = testMesh.node2cellsL2up[i,index];		
					
			uUpp[1] = testFields.densityNodes[pUp1];
			uUpp[2] = testFields.UxNodes[pUp1];
			uUpp[3] = testFields.UyNodes[pUp1];
			uUpp[4] = testFields.pressureNodes[pUp1];
					
			uDownp[1] = testFields.densityNodes[pDown1];
			uDownp[2] = testFields.UxNodes[pDown1];
			uDownp[3] = testFields.UyNodes[pDown1];
			uDownp[4] = testFields.pressureNodes[pDown1];
		
		elseif (ek_type == 2) ## quad element 
		
			pDown1 = testMesh.node2cellsL2down[i,index];
			pDown2 = testMesh.node2cellsL2down[i,index+1];
			
			pUp1 = testMesh.node2cellsL2up[i,index];		
			pUp2 = testMesh.node2cellsL2up[i,index+1];		
		
			uUpp[1] = 0.5*(testFields.densityNodes[pUp1]  + testFields.densityNodes[pUp2]);
			uUpp[2] = 0.5*(testFields.UxNodes[pUp1]       + testFields.UxNodes[pUp2]);
			uUpp[3] = 0.5*(testFields.UyNodes[pUp1]       + testFields.UyNodes[pUp2]);
			uUpp[4] = 0.5*(testFields.pressureNodes[pUp1] + testFields.pressureNodes[pUp2]);
					
			uDownp[1] = 0.5*(testFields.densityNodes[pDown1]  + testFields.densityNodes[pDown2]);
			uDownp[2] = 0.5*(testFields.UxNodes[pDown1]       + testFields.UxNodes[pDown2]);
			uDownp[3] = 0.5*(testFields.UyNodes[pDown1]       + testFields.UyNodes[pDown2]);
			uDownp[4] = 0.5*(testFields.pressureNodes[pDown1] + testFields.pressureNodes[pDown2]);
		
		
		end
		

		uRightp[1] = testFields.densityCells[ek];
		uRightp[2] = testFields.UxCells[ek];
		uRightp[3] = testFields.UyCells[ek];
		uRightp[4] = testFields.pressureCells[ek];					
					
	else
					
					
		uRightp = ComputeUPhysFromBoundaries(i,k, ek, uLeftp, nx,ny);
					
		##uDownp = deepcopy(uLeftp);
		##uUpp = deepcopy(uRightp);
		
		uDownp[1] = uLeftp[1];
		uDownp[2] = uLeftp[2];
		uDownp[3] = uLeftp[3];
		uDownp[4] = uLeftp[4];
		
		uUpp[1] = uRightp[1];
		uUpp[2] = uRightp[2];
		uUpp[3] = uRightp[3];
		uUpp[4] = uRightp[4];
	
	
	end
				
				
	ksi::Float64 = 1.0e-6;
				
	deltaJI = zeros(Float64,4);
	deltaIm = zeros(Float64,4);
	deltaJp = zeros(Float64,4);
	limLeft  = zeros(Float64,4);
	limRight  = zeros(Float64,4);
	UpRight = zeros(Float64,4);
	UpLeft = zeros(Float64,4);

				
	deltaJI[1] =  uRightp[1] - uLeftp[1];
	deltaJI[2] =  uRightp[2] - uLeftp[2];
	deltaJI[3] =  uRightp[3] - uLeftp[3];
	deltaJI[4] =  uRightp[4] - uLeftp[4];
		
	deltaIm[1] =  uLeftp[1]  - uDownp[1];
	deltaIm[2] =  uLeftp[2]  - uDownp[2];
	deltaIm[3] =  uLeftp[3]  - uDownp[3];
	deltaIm[4] =  uLeftp[4]  - uDownp[4];
		
	deltaJp[1] =  uUpp[1]  - uRightp[1];
	deltaJp[2] =  uUpp[2]  - uRightp[2];
	deltaJp[3] =  uUpp[3]  - uRightp[3];
	deltaJp[4] =  uUpp[4]  - uRightp[4];
		
	# deltaIm[1] = 0.0;
	# deltaIm[2] = 0.0;
	# deltaIm[3] = 0.0;
	# deltaIm[4] = 0.0;
	
	# deltaJp[1] = 0.0;
	# deltaJp[2] = 0.0;
	# deltaJp[3] = 0.0;
	# deltaJp[4] = 0.0;
	
					
				
	limLeft[1] = 0.5*Minmod_LimiterA( deltaIm[1], deltaJI[1], ksi);
	limLeft[2] = 0.5*Minmod_LimiterA( deltaIm[2], deltaJI[2], ksi);
	limLeft[3] = 0.5*Minmod_LimiterA( deltaIm[3], deltaJI[3], ksi);
	limLeft[4] = 0.5*Minmod_LimiterA( deltaIm[4], deltaJI[4], ksi);
						
	limRight[1] = 0.5*Minmod_LimiterA( deltaJI[1], deltaJp[1],  ksi);
	limRight[2] = 0.5*Minmod_LimiterA( deltaJI[2], deltaJp[2],  ksi);
	limRight[3] = 0.5*Minmod_LimiterA( deltaJI[3], deltaJp[3],  ksi);
	limRight[4] = 0.5*Minmod_LimiterA( deltaJI[4], deltaJp[4],  ksi);			
					
	
	UpLeft[1]  = uLeftp[1] + limLeft[1];
	UpLeft[2]  = uLeftp[2] + limLeft[2];
	UpLeft[3]  = uLeftp[3] + limLeft[3];
	UpLeft[4]  = uLeftp[4] + limLeft[4];
					
	UpRight[1] = uRightp[1] - limRight[1];
	UpRight[2] = uRightp[2] - limRight[2];	
	UpRight[3] = uRightp[3] - limRight[3];
	UpRight[4] = uRightp[4] - limRight[4];
												
			
	#return RoeFlux2d(UpRight,UpLeft, nx,ny,side,thermo.Gamma);		
	return AUSMplusFlux2d(UpRight,UpLeft, nx,ny,side,thermo.Gamma);		


end