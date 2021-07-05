
println("set thermo-physical properties ...");

@everywhere struct THERMOPHYSICS
	RGAS::Float64;
	Gamma::Float64;
	Cp::Float64;
	mu::Float64;
	kGas::Float64;		
end

#@everywhere thermo 
#@everywhere const thermoX = $thermo;




