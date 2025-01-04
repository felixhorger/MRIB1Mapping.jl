
module MRIBMapping

	import GLM

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


	# TODO: Lin2020, but include the coil combination from https://github.com/MagneticResonanceImaging/MRIFieldmaps.jl
	# Maybe need to use full image not only phases
	# TODO: instead of computing phase1, phase2, ... and then computing the difference, it is better to
	# compute angle(signal1 * conj(signal2)), see Robinson2011
	# TODO: Channelwise like in Robinson2011
	"""
		Note the units are rad/[TE]
	"""
	function phases2Δω(ϕ::AbstractArray{<: Real, N}, TEs::AbstractVector{<: Real}) where N
		@assert size(ϕ, N) == 3
		@assert all(x -> x > 0.0, diff(TEs))
		shape = size(ϕ)[1:N-1]
		spatial_indices = CartesianIndices(shape)
		# Get window in which no aliasing occurs
		ΔTEs = @views TEs[2:3] .- TEs[1]
		Δω_aliasing = 2π / ΔTEs[1] # TODO: not minimum ...? plus, a factor of two missing? Nyquist!
		# Remove phase from first TE
		ϕ = @views rem2pi.(ϕ[spatial_indices, 2:3] .- ϕ[spatial_indices, 1], RoundNearest)
		# Note: Use remainder, because we don't expect these to lie outside the window Δω_aliasing
		# Project expected phase and correct accordingly
		ϕ_corrected = @views begin
			# Get estimate on frequency offset
			Δω = ϕ[spatial_indices, 1] ./ ΔTEs[1] # Will be reused
			# Estimate expected phase
			ϕ_expected = Δω .* ΔTEs[2]
			# Compute number of full revolutions done by the magnetisation and add to measured phase
			@views @. ϕ[spatial_indices, 2] += 2π * round(Int, (ϕ_expected - ϕ[spatial_indices, 2]) / 2π)
		end
		# Perform linear fit
		ΔTEs = reshape(ΔTEs, 2, 1)
		Δω_err = similar(Δω)
		for X in spatial_indices
			fit = @views GLM.lm(ΔTEs, ϕ[X, :])
			Δω[X] = GLM.coef(fit)[1]
			Δω_err[X] = GLM.stderror(fit)[1]
		end
		return Δω, Δω_err, Δω_aliasing
	end
end

