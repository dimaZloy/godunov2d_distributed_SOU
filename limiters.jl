



## template <class T> T Minmod_Limiter(T x, T y){return (sign(x)+sign(y))*fabs(x*y)/(fabs(x)+fabs(y)+1.0e-6);}
## template <class T> T sign(T x){return x<0 ? -1: x>0 ? 1:0;}

# @inline @everywhere function signT(x::Float64)::Float64
		
	# if x<0.0
		# return -1.0;
	# elseif x>0.0
		# return 1.0;
	# end
	
	# return 0.0;
# end

# @inline @everywhere function absT(x::Float64)::Float64
		
	# if x<0.0
		# return -1.0*x;
	# end

	# return x;
# end


# @inline @everywhere function Minmod_LimiterT(x::Float64, y::Float64, c::Float64)::Float64
	# return deepcopy( ( signT(x) + signT(y) )*absT(x*y)/( absT(x) + absT(y) + c));
# end

@inline @everywhere function Minmod_Limiter(x::Float64, y::Float64, c::Float64)::Float64
	return ( sign(x) + sign(y) )*abs(x*y)/( abs(x) + abs(y) + c);
end

# @inline @everywhere function Minmod_Limiter(x::Array{Float64,1}, y::Array{Float64,1},c::Float64)::Array{Float64,1}
   # return (sign(x)+sign(y)).*abs(x.*y)./(abs(x).+abs(y)+1.0e-6);
# end

@inline @everywhere function vanLeer_LimiterA(a::Float64, b::Float64, c::Float64)::Float64
	return (a*b + abs(a*b))/(a+b+c);
end

@inline @everywhere function Minmod_LimiterA(x::Float64, y::Float64, c::Float64)::Float64
	return  0.5*( sign(x)+sign(y) )*min( abs(x), abs(y) );
end


@inline @everywhere function vanAlbada_LimiterA(a,b,c)
	return  (a.*(b.*b+c) + b.*(a.*a+c))./(a.*a+b.*b+2.0*c)
end

@inline @everywhere function vanAlbada_LimiterB(a::Float64,b::Float64,c::Float64)::Float64
	return  (a.*(b.*b+c) + b.*(a.*a+c))./(a.*a+b.*b+2.0*c)
end