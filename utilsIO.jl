
function saveResults2VTK(
		filename::String, 
		testMesh::mesh2d, 
		variable::Array{Float64,1}, variableName::String )

	vtkfile = vtk_grid(filename, testMesh.xNodes, testMesh.yNodes, testMesh.VTKCells);
	vtk_point_data(vtkfile, variable, variableName);
	outfiles = vtk_save(vtkfile);
	
end

function readSolution(fileName, nNodes)
	
	io = open(fileName,"r");
	line  = readline(io);
	line  = readline(io);
	
	xNodes = zeros(Float64,nNodes);
	yNodes = zeros(Float64,nNodes);
	nodesSolution = zeros(Float64,nNodes);
	
	
	for i=1:nNodes
		
		line = readline(io);
		z = split(line);
		xNodes[i] = parse(Float64,z[1]);
		yNodes[i] = parse(Float64,z[2]);
		nodesSolution[i] = parse(Float64,z[3]);
		
	end	
	close(io);
	
	return xNodes, yNodes, nodesSolution; 

end


function saveSolution(fileName, xNodes, yNodes, nodesSolution)

io = open(fileName,"w");
for i=1:size(xNodes,1)
	writedlm(io, [xNodes[i]  yNodes[i] nodesSolution[i,1] nodesSolution[i,2] nodesSolution[i,3] nodesSolution[i,4] ], '\t' );
end


close(io);

end

function saveResiduals(fileName, timeVector, residualsVector1, residualsVector2, residualsVector3, residualsVector4)

io = open(fileName,"w");
	
for i=1:size(timeVector,1)
	writedlm(io, [ timeVector[i]  residualsVector1[i] residualsVector2[i] residualsVector3[i] residualsVector4[i]  ], '\t' );
end

close(io);

end
