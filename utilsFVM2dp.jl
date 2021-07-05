
# utilities for FVM 


@everywhere @inline function phs2dcns2dcellsSA(
	ACons::SharedArray{Float64,2}, testFields::fields2d_shared, gamma::Float64)

	N::Int64 = size(testFields.densityCells,1);
	
	#ACons = zeros(Float64,N,4);

	for i = 1:N
		ACons[i,1] = testFields.densityCells[i];
		ACons[i,2] = testFields.densityCells[i]*testFields.UxCells[i];
		ACons[i,3] = testFields.densityCells[i]*testFields.UyCells[i];
		ACons[i,4] = testFields.pressureCells[i]/(gamma-1.0) + 0.5*testFields.densityCells[i]*(	testFields.UxCells[i]*testFields.UxCells[i] +  testFields.UyCells[i]*testFields.UyCells[i] );

	end #for
	
	#return ACons;
end


@everywhere function cells2nodesSolutionReconstructionWithStencilsImplicitSA(
		testMesh::mesh2d,testFields::fields2d_shared)	

	node_solution = zeros(Float64,testMesh.nNodes,4); 
	
	for J=1:testMesh.nNodes
	
		det::Float64 = 0.0;
		
		for j = 1:testMesh.nNeibCells
		
			neibCell::Int64 = testMesh.cell_clusters[J,j]; 
			
			if (neibCell !=0)
				 wi::Float64 = testMesh.node_stencils[J,j];
				 #node_solution[J,:] += cell_solution[neibCell,:];
				 node_solution[J,1] += testFields.densityCells[neibCell]*wi;
				 node_solution[J,2] += testFields.UxCells[neibCell]*wi;
				 node_solution[J,3] += testFields.UyCells[neibCell]*wi;
				 node_solution[J,4] += testFields.pressureCells[neibCell]*wi;
				 
				 det += wi;
			end
		end
		if (det!=0)
			node_solution[J,:] = node_solution[J,:]/det; 
			
		end
	end

	
	for J=1:testMesh.nNodes
	
		testFields.densityNodes[J] = node_solution[J,1]; 
		testFields.UxNodes[J] 	   = node_solution[J,2]; 
		testFields.UyNodes[J] 	   = node_solution[J,3]; 
		testFields.pressureNodes[J] =  node_solution[J,4]; 
		
	end

end



@everywhere  @inline function cells2nodesSolutionReconstructionWithStencilsSA(
		testMesh::mesh2d,cell_solution::SharedArray{Float64,1} ) ::Array{Float64,1}

node_solution = zeros(Float64,testMesh.nNodes); 

for J=1:testMesh.nNodes
	det::Float64 = 0.0;
	for j = 1:testMesh.nNeibCells
		neibCell::Int64 = testMesh.cell_clusters[J,j]; 
		if (neibCell !=0)
			wi::Float64 = testMesh.node_stencils[J,j];
			node_solution[J] += cell_solution[neibCell]*wi;
			det += wi;
		end
	end
	if (det!=0)
		node_solution[J] = node_solution[J]/det; 
	end
end

return node_solution;	

end
