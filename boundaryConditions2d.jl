
@everywhere function ComputeUPhysFromBoundaries(i,k,neib_cell, cur_cell, nx,ny)

	bnd_cell = zeros(Float64,4);



		if (neib_cell == -3) #inlet

            bnd_cell[1] = 1.4;
            bnd_cell[2] = 300.00;
            bnd_cell[3] = 0.0;
            bnd_cell[4] = 10000.0;

		elseif (neib_cell == -2) #walls

	 	   	#bnd_cell = cur_cell;	
	        #bnd_cell = updateVelocityFromCurvWall(i,k,bnd_cell,nx,ny);
			bnd_cell = updateVelocityFromCurvWall(i,k,cur_cell,nx,ny);

		elseif (neib_cell == -1) # outlet

			bnd_cell = cur_cell;	
					
		end	

			

	return bnd_cell; 
end


@everywhere function updateVelocityFromCurvWall(i::Int64, k::Int64, U, nx::Float64, ny::Float64)

# High-Order Accurate Implementation of Solid Wall Boundary Conditions in Curved Geometries, 
# Lilia Krivodonova and Marsha Berger, Courant Institute of Mathematical Sciences, New York, NY 10012

# a = U[1]*(ny*ny - nx*nx) - 2.0*nx*ny*U[2];
# b = U[2]*(nx*nx - ny*ny) - 2.0*nx*ny*U[1];


	Un = deepcopy(U); 

        Un[2] = U[2]*(ny*ny - nx*nx) - 2.0*nx*ny*U[3];
        Un[3] = U[3]*(nx*nx - ny*ny) - 2.0*nx*ny*U[2];


	return Un;	
end

