

# @everywhere function doEulerExplicit(
		# betta::Float64, dtX::Float64,
		# testMeshDistrX::mesh2d_shared, testfields2dX::fields2d_shared, thermoX::THERMOPHYSICS, cellsThreadsX::SharedArray{Int64,2}
		# UconsCellsOldX::SharedArray{Float64,2}, iFLUXX::SharedArray{Float64,2},  UconsCellsNewX::SharedArray{Float64,2})

	# @sync @distributed for p in workers()	
	
		# beginCell::Int64 = cellsThreadsX[p-1,1];
		# endCell::Int64 = cellsThreadsX[p-1,2];
		# #println("worker: ",p,"\tbegin cell: ",beginCell,"\tend cell: ", endCell);
										 
		# SecondOrderUpwindM2(beginCell ,endCell, betta, dtX, testMeshDistrX, testfields2dX, thermoX, UconsCellsOldX, iFLUXX, UconsCellsNewX);
					
	# end
			
	
	# @everywhere finalize(SecondOrderUpwindM2);				

# end


@everywhere function calcOneStage(
		betta::Float64, dtX::Float64,
		testMeshDistrX::mesh2d_shared, testfields2dX::fields2d_shared, thermoX::THERMOPHYSICS, cellsThreadsX::SharedArray{Int64,2},
		UconsCellsOldX::SharedArray{Float64,2}, iFLUXX::SharedArray{Float64,2},  UconsCellsNewX::SharedArray{Float64,2})

	@sync @distributed for p in workers()	
	
		beginCell::Int64 = cellsThreadsX[p-1,1];
		endCell::Int64 = cellsThreadsX[p-1,2];
		#println("worker: ",p,"\tbegin cell: ",beginCell,"\tend cell: ", endCell);
										 
		SecondOrderUpwindM2(beginCell ,endCell, betta, dtX, testMeshDistrX, testfields2dX, thermoX, UconsCellsOldX, iFLUXX, UconsCellsNewX);
					
	end
			
	
	@everywhere finalize(SecondOrderUpwindM2);				

end



@everywhere function doExplicitRK3TVD(
		betta::Float64, dtX::Float64,
		testMeshDistrX::mesh2d_shared, testfields2dX::fields2d_shared, thermoX::THERMOPHYSICS, cellsThreadsX::SharedArray{Int64,2},
		UconsCellsOldX::SharedArray{Float64,2}, iFLUXX::SharedArray{Float64,2},  
		UconsCellsNew1X::SharedArray{Float64,2}, UconsCellsNew2X::SharedArray{Float64,2}, UconsCellsNew3X::SharedArray{Float64,2}, UconsCellsNewX::SharedArray{Float64,2})


			
			## based on
			## JOURNAL OF COMPUTATIONAL PHYSICS 83, 32-78 (1989)
			## Efficient Implementation   of Essentially Non-oscillatory   Shock-Capturing   Schemes, II
			## CHI- WANG SHU  AND   STANLEY OSHERS

	
			# @sync @distributed for p in workers()	
	
				# beginCell::Int64 = cellsThreadsX[p-1,1];
				# endCell::Int64 = cellsThreadsX[p-1,2];
				# #println("worker: ",p,"\tbegin cell: ",beginCell,"\tend cell: ", endCell);
										 
				# SecondOrderUpwindM2(beginCell ,endCell, 1.0, dtX, testMeshDistrX, testfields2dX, thermoX, UconsCellsOldX, iFLUXX, UconsCellsNew1X);
					
			# end
			# @everywhere finalize(SecondOrderUpwindM2);				
			
			calcOneStage(1.0, dtX, testMeshDistrX, testfields2dX , thermoX , cellsThreadsX, UconsCellsOldX, iFLUXX,  UconsCellsNew1X); 	
		

			@sync @distributed for p in workers()	
				a1::Float64  = 3.0/4.0;
				a2::Float64  = 1.0/4.0;
				Gamma::Float64 = thermoX.Gamma;
				
				beginCell::Int64 = cellsThreadsX[p-1,1];
				endCell::Int64 = cellsThreadsX[p-1,2];
				
				
				
				for i = beginCell:endCell
					UconsCellsNew2X[i,1] = UconsCellsOldX[i,1].*a1 .+ UconsCellsNew1X[i,1].*a2;
					UconsCellsNew2X[i,2] = UconsCellsOldX[i,2].*a1 .+ UconsCellsNew1X[i,2].*a2;
					UconsCellsNew2X[i,3] = UconsCellsOldX[i,3].*a1 .+ UconsCellsNew1X[i,3].*a2;
					UconsCellsNew2X[i,4] = UconsCellsOldX[i,4].*a1 .+ UconsCellsNew1X[i,4].*a2;

					testfields2dX.densityCells[i] = UconsCellsNew2X[i,1];
					testfields2dX.UxCells[i] 	  = UconsCellsNew2X[i,2]/UconsCellsNew2X[i,1];		
					testfields2dX.UyCells[i] 	  = UconsCellsNew2X[i,3]/UconsCellsNew2X[i,1];
					testfields2dX.pressureCells[i] = (Gamma-1.0)*( UconsCellsNew2X[i,4] - 0.5*( UconsCellsNew2X[i,2]*UconsCellsNew2X[i,2] + UconsCellsNew2X[i,3]*UconsCellsNew2X[i,3] )/UconsCellsNew2X[i,1] );

					testfields2dX.aSoundCells[i] = sqrt( Gamma * testfields2dX.pressureCells[i]/testfields2dX.densityCells[i] );
					testfields2dX.VMAXCells[i]  = sqrt( testfields2dX.UxCells[i]*testfields2dX.UxCells[i] + testfields2dX.UyCells[i]*testfields2dX.UyCells[i] ) + testfields2dX.aSoundCells[i];
					
				end
				
			
			end		
			
		
			
			# @sync @distributed for p in workers()	
	
				# beginCell::Int64 = cellsThreadsX[p-1,1];
				# endCell::Int64 = cellsThreadsX[p-1,2];
				# #println("worker: ",p,"\tbegin cell: ",beginCell,"\tend cell: ", endCell);
										 
				# SecondOrderUpwindM2(beginCell ,endCell, 1.0/4.0, dtX, testMeshDistrX, testfields2dX, thermoX, UconsCellsNew2X, iFLUXX, UconsCellsNew3X);
					
			# end
			# @everywhere finalize(SecondOrderUpwindM2);			


			calcOneStage(1.0/4.0, dtX, testMeshDistrX, testfields2dX , thermoX , cellsThreadsX, UconsCellsNew2X, iFLUXX,  UconsCellsNew3X); 				
			
			
			@sync @distributed for p in workers()	
				b1::Float64  = 1.0/3.0;
				b2::Float64  = 2.0/3.0;
				Gamma::Float64 = thermoX.Gamma;
				
				beginCell::Int64 = cellsThreadsX[p-1,1];
				endCell::Int64 = cellsThreadsX[p-1,2];
				
				for i = beginCell:endCell
				
					UconsCellsNew2X[i,1] = UconsCellsOldX[i,1].*b1 .+ UconsCellsNew3X[i,1].*b2;
					UconsCellsNew2X[i,2] = UconsCellsOldX[i,2].*b1 .+ UconsCellsNew3X[i,2].*b2;
					UconsCellsNew2X[i,3] = UconsCellsOldX[i,3].*b1 .+ UconsCellsNew3X[i,3].*b2;
					UconsCellsNew2X[i,4] = UconsCellsOldX[i,4].*b1 .+ UconsCellsNew3X[i,4].*b2;

					testfields2dX.densityCells[i] = UconsCellsNew2X[i,1];
					testfields2dX.UxCells[i] 	  = UconsCellsNew2X[i,2]/UconsCellsNew2X[i,1];		
					testfields2dX.UyCells[i] 	  = UconsCellsNew2X[i,3]/UconsCellsNew2X[i,1];
					testfields2dX.pressureCells[i] = (Gamma-1.0)*( UconsCellsNew2X[i,4] - 0.5*( UconsCellsNew2X[i,2]*UconsCellsNew2X[i,2] + UconsCellsNew2X[i,3]*UconsCellsNew2X[i,3] )/UconsCellsNew2X[i,1] );

					testfields2dX.aSoundCells[i] = sqrt( Gamma * testfields2dX.pressureCells[i]/testfields2dX.densityCells[i] );
					testfields2dX.VMAXCells[i]  = sqrt( testfields2dX.UxCells[i]*testfields2dX.UxCells[i] + testfields2dX.UyCells[i]*testfields2dX.UyCells[i] ) + testfields2dX.aSoundCells[i];
					
				end
			
				
			
			end		
			

			# @sync @distributed for p in workers()	
	
				# beginCell::Int64 = cellsThreadsX[p-1,1];
				# endCell::Int64 = cellsThreadsX[p-1,2];
				# #println("worker: ",p,"\tbegin cell: ",beginCell,"\tend cell: ", endCell);
										 
				# SecondOrderUpwindM2(beginCell ,endCell, 2.0/3.0, dtX, testMeshDistrX, testfields2dX, thermoX, UconsCellsNew2X, iFLUXX, UconsCellsNewX);
					
			# end
			# @everywhere finalize(SecondOrderUpwindM2);				
			
			
			calcOneStage(2.0/3.0, dtX, testMeshDistrX, testfields2dX , thermoX , cellsThreadsX, UconsCellsNew2X, iFLUXX,  UconsCellsNewX); 	
			

end