@inline @everywhere function RoeFlux2d(
		primL::Array{Float64,1},primR::Array{Float64,1},
		nx::Float64,ny::Float64,side::Float64, gamma::Float64)::Array{Float64,1}

# [1] P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and
# Difference Schemes, Journal of Computational Physics, 43, pp. 357-372.
#
# [2] H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free,
# Boundary-Layer Resolving, Rotated-Hybrid Riemann Solvers,
# Journal of Computational Physics, 227, pp. 2560-2581, 2008.

# -------------------------------------------------------------------------
#  Input:   primL(1:4) =  left state (rhoL, uL, vL, pL)
#           primR(1:4) = right state (rhoR, uR, vR, pR)
#               njk(2) = Face normal (L -> R). Must be a unit vector.
#
# Output:    flux(1:4) = numerical flux
#                  wsn = half the max wave speed (to be used for time step calculations)
# -------------------------------------------------------------------------



#%Tangent vector (Do you like it? Actually, Roe flux can be implemented without any tangent vector. See "I do like CFD, VOL.1" for details.)

  mx::Float64 = -ny;
  my::Float64 =  nx;


#  Left state
    rhoL::Float64 = primL[1];
      uL::Float64 = primL[2];
      vL::Float64 = primL[3];
     unL::Float64 = uL*nx+vL*ny;
     umL::Float64 = uL*mx+vL*my;
      pL::Float64 = primL[4];
      aL::Float64 = sqrt(gamma*pL/rhoL);
      HL::Float64 = aL*aL/(gamma-1.0) + 0.5*(uL*uL+vL*vL);
#  Right state
    rhoR::Float64 = primR[1];
      uR::Float64 = primR[2];
      vR::Float64 = primR[3];
     unR::Float64 = uR*nx+vR*ny;
     umR::Float64 = uR*mx+vR*my;
      pR::Float64 = primR[4];
      aR::Float64 = sqrt(gamma*pR/rhoR);
      HR::Float64 = aR*aR/(gamma-1.0) + 0.5*(uR*uR+vR*vR);

# First compute the Roe Averages
    RT::Float64 = sqrt(rhoR/rhoL);
   rho::Float64 = RT*rhoL;
     u::Float64 = (uL+RT*uR)/(1.0+RT);
     v::Float64 = (vL+RT*vR)/(1.0+RT);
     H::Float64 = (HL+RT* HR)/(1.0+RT);
     a::Float64 = sqrt( (gamma-1.0)*(H-0.5*(u*u+v*v)) );
    un::Float64 = u*nx+v*ny;
    um::Float64 = u*mx+v*my;

# Wave Strengths
   drho::Float64 = rhoR - rhoL; 
     dp::Float64 =   pR - pL;
    dun::Float64 =  unR - unL;
    dum::Float64 =  umR - umL;

  LdU = zeros(Float64,4);
		
  LdU[1] = (dp - rho*a*dun )/(2.0*a*a);
  LdU[2] = rho*dum;
  LdU[3] =  drho - dp/(a*a);
  LdU[4] = (dp + rho*a*dun )/(2.0*a*a);

# Wave Speed

  ws = zeros(Float64,4);
  
  ws[1] = abs(un-a);
  ws[2] = abs(un);
  ws[3] = abs(un);
  ws[4] = abs(un+a);

# Harten's Entropy Fix JCP(1983), 49, pp357-393: 
# only for the nonlinear fields.
  dws = zeros(Float64,4);
  dws[1]=1.0/5.0; 
  if ws[1] < dws[1]
		ws[1] = 0.5*( ws[1]*ws[1]/dws[1] + dws[1] );
  end
  
  dws[4]=1.0/5.0; 
  if ws[4] < dws[4]
		ws[4] = 0.5*( ws[4]*ws[4]/dws[4] + dws[4] ); 
  end

#Right Eigenvectors
  Rv = zeros(Float64,4,4);	


  Rv[1,1] = 1.0;
  Rv[2,1] = u - a*nx;
  Rv[3,1] = v - a*ny;
  Rv[4,1] = H - un*a;

  Rv[1,2] = 0.0;
  Rv[2,2] = mx;
  Rv[3,2] = my;
  Rv[4,2] = um;

  Rv[1,3] = 1.0;
  Rv[2,3] = u;
  Rv[3,3] = v;
  Rv[4,3] = 0.5*(u*u + v*v);

  Rv[1,4] = 1.0;
  Rv[2,4] = u + a*nx;
  Rv[3,4] = v + a*ny;
  Rv[4,4] = H + un*a;


#Dissipation Term
  diss = zeros(Float64,4);
  
  # @simd for i=1:4;
	# @simd for j=1:4;
		# diss[i] = diss[i] + ws[j]*LdU[j]*Rv[i,j];
    # end 
  # end

  for i=1:4;
	for j=1:4;
		diss[i] = diss[i] + ws[j]*LdU[j]*Rv[i,j];
    end 
  end

  
#Compute the flux.
  fL = zeros(Float64,4)
  fL[1] = rhoL*unL;
  fL[2] = rhoL*unL * uL + pL*nx;
  fL[3] = rhoL*unL * vL + pL*ny;
  fL[4] = rhoL*unL * HL;

  fR = zeros(Float64,4)
  fR[1] = rhoR*unR;
  fR[2] = rhoR*unR * uR + pR*nx;
  fR[3] = rhoR*unR * vR + pR*ny;
  fR[4] = rhoR*unR * HR;

  flux = zeros(Float64,4)
  # flux[1] = 0.5 * (fL[1] + fR[1] - diss[1]);
  # flux[2] = 0.5 * (fL[2] + fR[2] - diss[2]);
  # flux[3] = 0.5 * (fL[3] + fR[3] - diss[3]);
  # flux[4] = 0.5 * (fL[4] + fR[4] - diss[4]);
  #wsn = 0.5*(abs(un) + a);  #Normal max wave speed times half
  
  flux[1] = -0.5 * (fL[1] + fR[1] - diss[1]) * side;
  flux[2] = -0.5 * (fL[2] + fR[2] - diss[2]) * side;
  flux[3] = -0.5 * (fL[3] + fR[3] - diss[3]) * side;
  flux[4] = -0.5 * (fL[4] + fR[4] - diss[4]) * side;
  

  return flux;	
  
end