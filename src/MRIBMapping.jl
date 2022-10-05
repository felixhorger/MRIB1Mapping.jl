
module MRIBMapping
	
	function afi(recons::AbstractArray{<: Complex, N}, TR_ratio::Real, α_nominal::Real) where N
		@assert TR_ratio > 1
		# Get dimensions
		shape = size(recons)
		spatial_shape = shape[1:N-1]
		# Iterate spatial locations
		relB1 = Array{Float64, N-1}(undef, spatial_shape)
		for X in CartesianIndices(spatial_shape)
			signal_ratio = sqrt(abs2(recons[X, 2]) / abs2(recons[X, 1]))
			cosα = (signal_ratio * TR_ratio - 1) / (TR_ratio - signal_ratio)
			if abs2(cosα) < 1
				relB1[X] = acos(cosα) / α_nominal
			else
				relB1[X] = 0
			end
		end
		return relB1
	end
end

