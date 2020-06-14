
using Images, TestImages, FFTW, Plots
using LinearAlgebra

# Direct <- IFFT | FFT -> Recip
# ρ ∈ ℝ  -- FFT --> F = |F| exp(i ϕ) ∈ ℂ 
# F ∈ ℂ  -- FFT --> ρ' ∈ ℂ, real(ρ,l)
# |F| is the modulus
# ϕ   is the phase  = tan⁻¹(Fimag/Freal)
# exp(iϕ) = cos(ϕ) + i sin(ϕ) = Freal + i Fimag
# ϕ = tan⁻¹(i Fimag/(i Freal)) = tan⁻¹(Fimag/Freal) = angle(F)

# 0. Get Fabio and Lena
fabio = testimage("fabio_gray_512")
lena  = testimage("lena_gray_512")
save("fabio.png",fabio)
save("lena.png",lena)

fabio_array = channelview(fabio)
lena_array  = channelview(lena)

# 1. FFT of Fabio and Lena
fabio_array_recip     = fftshift(fft(fabio_array))
lena_array_recip      = fftshift(fft(lena_array))
lena_array_recip_mod  = abs.(  lena_array_recip)
fabio_array_recip_mod = abs.(  fabio_array_recip)
fabio_array_recip_pha = angle.(fabio_array_recip)
lena_array_recip_pha  = angle.(lena_array_recip)

# Reconstruct Fabio
fabio2_array      = ifft(fabio_array_recip    )
fabio2      = Gray.(real(fabio2_array     ))
save("fabio2.png",     fabio2     )


# Reconstruct Lena
lena2_array      = ifft(lena_array_recip    )
lena2      = Gray.(real(lena2_array     ))
save("lena2.png",     lena2     )

# ------------------------------------------------------
# Phases are more important than moduli
# ------------------------------------------------------
# Mix moduli and phase of two images, then reconstruct:
# 1. Perform fft of Fabio and Lena
#     Fabio_fft = fft(Fabio)
#     Lena_fft  = fft(Lena)
# 2. Interchange their phases
#     F_w_L = |Fabio_fft| exp(i ϕ_Lena ) 
#     L_w_F = | Lena_fft| exp(i ϕ_Fabio)
# 3. Reconstruct Fabio and Lena with ifft
#     Fabio_w_Lenaϕ = ifft( F_w_L )
#     Lena_w_Fabioϕ = ifft( L_w_F )

fabiomod_lenapha_array_recip = abs.(fabio_array_recip_mod) .* exp.(im*lena_array_recip_pha)
lenamod_fabiopha_array_recip = abs.(lena_array_recip_mod)  .* exp.(im*fabio_array_recip_pha)

fabiomod_lenapha_array = ifft(fabiomod_lenapha_array_recip)
lenamod_fabiopha_array = ifft(lenamod_fabiopha_array_recip)
fabiomod_lenapha = Gray.(real(fabiomod_lenapha_array))
lenamod_fabiopha = Gray.(real(lenamod_fabiopha_array))
save("fabiomod_lenapha.png", fabiomod_lenapha)
save("lenamod_fabiopha.png", lenamod_fabiopha)

#
#_fabio_array_recip = fabio_array_recip_mod .* exp.(im*fabio_array_recip_pha)
#fabio_recip      = Gray.(log10.(abs.(fabio_array_recip    )))
#fabio_recip_mod  = Gray.(log10.(abs.(fabio_array_recip_mod)))
#fabio_recip_pha  = Gray.(            fabio_array_recip_pha)
#save("fabio_recip.png",    fabio_recip    )
#save("fabio_recip_mod.png",fabio_recip_mod)
#save("fabio_recip_pha.png",fabio_recip_pha)
#fabio2_array_fmod = ifft(fabio_array_recip_mod)
#fabio2_array_fpha = ifft(fabio_array_recip_pha)
#fabio2_fmod = Gray.(real(fabio2_array_fmod))
#fabio2_fpha = Gray.(real(fabio2_array_fpha))
#save("fabio2_fmod.png",fabio2_fmod)
#save("fabio2_fpha.png",fabio2_fpha)

#plot(fabio)
#
#lena_array = channelview(lena)
#lena_recip      = Gray.(log10.(abs.(lena_array_recip    )))
#lena_recip_mod  = Gray.(log10.(abs.(lena_array_recip_mod)))
#lena_recip_pha  = Gray.(            lena_array_recip_pha)
#save("lena_recip.png",    lena_recip    )
#save("lena_recip_mod.png",lena_recip_mod)
#save("lena_recip_pha.png",lena_recip_pha)
#lena2_array_fmod = ifft(lena_array_recip_mod)
#lena2_array_fpha = ifft(lena_array_recip_pha)
#lena2_fmod = Gray.(real(lena2_array_fmod))
#lena2_fpha = Gray.(real(lena2_array_fpha))
#save("lena2_fmod.png",lena2_fmod)
#save("lena2_fpha.png",lena2_fpha)
#
#fabiomod_lenapha_recip = Gray.(log10.(abs.(fabiomod_lenapha_array_recip)))
#lenamod_fabiopha_recip = Gray.(log10.(abs.(lenamod_fabiopha_array_recip)))
#save("fabiomod_lenapha_recip.png", fabiomod_lenapha_recip)
#save("lenamod_fabiopha_recip.png", lenamod_fabiopha_recip)
#
# fabio_array = channelview(fabio)
#fabio_array = float.(fabio)
#lena_array  = float.(lena)
