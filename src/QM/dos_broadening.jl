
@doc raw"""
   dos_broadening(es,siteweight=missing;width=0.05)

Broadening of DOS, or PDOS. Convolutes with a test line-shape.

# Arguments
- `es` is an array with the energy spectrum
- `siteweight` is an array with the contribution of a site orbital to each eigenstate
- `width` graduates the smearing
"""
function dos_broadening(es,siteweight=missing;width=0.05)
   Nshape = 20001
   Ndos   = 10001 # Nshape ÷ 2 + 1
   ϵmin = -maximum(abs.(es)) - 1 # some margin for broadening at extrema.
   ϵmax =  maximum(abs.(es)) + 1 # symmetric around 0
   shape = lorentz_like.(1:Nshape,Ndos,inv(Ndos) * (ϵmax - ϵmin),width)
   # Normalize (number of particles must be conserved)
   shape /= sum(shape)

   # Scaling to map spectra ∈ [ϵmin,ϵmax] ↦ [1,Ndos]
   scaling = (Ndos - 1) / (ϵmax - ϵmin)

   # Default is total DOS (weight=1), otherwise partial DOS
   weight = ifelse(ismissing(siteweight),ones(length(es)),siteweight)

   # Convolute
   dos = zeros(Ndos)
   for i in 1:length(es)
      # 0.5 for truncation, 1 for 1-indexing
      ishape = trunc(Int,scaling*(es[i]-ϵmin) + 1.5)
      # (f * g)(t) = ∫ f(τ) g(t-τ) dτ
      for j in 1:Ndos
         dos[j] += weight[i]*shape[j+Ndos-ishape]
      end
   end

   # Finer energy grid of same size as dos, ∈ [ϵmin,ϵmax]
   e_dos = zeros(Ndos)
   for i in 1:Ndos
      # spacing in dos grid to spacing in ϵ
      Δϵ = inv(Ndos) * (ϵmax - ϵmin)
      e_dos[i] = ϵmin + Δϵ*(i-1)
   end
   e_dos,dos
end

function lorentz_like(i,i₀,Δϵ,width)
  dummy = Δϵ*(i-i₀)/(0.5width)
  inv(1 + dummy^2)
end
