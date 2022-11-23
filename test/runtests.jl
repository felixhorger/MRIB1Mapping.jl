
using Revise
import MRIBMapping
import PyPlot as plt

# B1 mapping
# TODO

# B0 mapping
TEs = Float64[0.5, 2, 20]
shape = (128, 128)
Δω = @. 1 - 2 * $rand(shape...)
ϕ0 = @. 2π * (1 - 2 * $rand(shape...))
ϕ = Array{Float64, 3}(undef, shape..., 3);
ϕ_unwrapped = Array{Float64, 3}(undef, shape..., 3);
error = @. π/30 * (1 - 2 * $rand(shape..., 3))
for i = 1:3
	@. ϕ_unwrapped[:, :, i] = ϕ0 + Δω * TEs[i] + error[:, :, i]
	@views @. ϕ[:, :, i] = mod2pi(ϕ_unwrapped[:, :, i])
end

Δω_fit, Δω_fit_err, Δω_aliasing = MRIBMapping.phase2Δω(ϕ, TEs)
fig, axs = plt.subplots(2, 2; sharex=true, sharey=true)
axs[1, 1].imshow(Δω)
axs[1, 2].imshow(Δω_fit)
axs[2, 1].imshow(Δω - Δω_fit; vmin=-0.1, vmax=0.1)
axs[2, 2].imshow(Δω_fit_err; vmin=-0.1, vmax=0.1)

#= Debugging B0 mapping
e = MRIBMapping.phase2Δω(ϕ, TEs)
plt.figure()
i = 1
j = 2
a = ϕ_unwrapped[i, j, :] .- ϕ_unwrapped[i, j, 1]
b = [ϕ_unwrapped[i, j, 2] - ϕ_unwrapped[i, j, 1], e[i, j]]
#b = e[i, j, :]
#a[2] += 2π
plt.plot(TEs, a, "+-")
plt.plot(TEs[2:3], b, "+-")

fig, axs = plt.subplots(1, 2; sharex=true, sharey=true)
axs[1].imshow(-e + ϕ[:, :, 3] - ϕ[:, :, 1], vmin=-1, vmax=1)

plt.figure()
plt.imshow(Δω_fit)
=#


