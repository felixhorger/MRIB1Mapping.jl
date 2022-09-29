
module MRIB1Mapping
	
	function afi(recons::AbstractArray{<: Complex, N}, TRs::NTuple{2, Real}) where N
		# Get dimensions
		shape = size(recon)
		spatial_shape = shape[1:N-1]
		num_TR = shape[N]
		# Iterate spatial locations
		n = TRs[2] / TRs[1]
		α = Array{Float64, N}(undef, spatial_shape)
		for X in CartesianIndices(recon)
			r = recon[X, 1] / recon[X, 2]
			α[X] = acos( (r * n - 1) / (n - r) )
		end
		return α
	end
end

