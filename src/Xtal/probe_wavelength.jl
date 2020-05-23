
# Section 2.1 in Modern Condensed Matter Physics by Girvin, ...
m2pm = 1e12 # meter to picometer
eV2J = 1.602176634e-19 # electronvolt to Joule

c = 299792458 # m/s; speed of light in vacuum (exact)
h = 6.626070e-34 # J⋅s; Planck constant

#c_as = c * m2pm # Å/s; speed of light in vacuum (exact)
#h_evs = 4.135667696e-15 # eV⋅s; Planck constant
#hc_angs_eV = h_evs * c_as
#hc_angs_keV = h_evs * c_as * 1e-3
#λ = hc_angs_keV / E_keV # pm

println("The size of atoms is ∼ 100 pm. Matter structure below the wavelength is not resolved.")
# ==================
# Photons, λ = hc/E
# ==================
# A typical soft X-ray source is ≈ 5 keV
# A typical hard X-ray source is > 10 keV
E = eV2J * 5e3
λp = h * c / E * m2pm
println("Photons of energy  $(E/eV2J) eV have a wavelength = $λp pm.")


# ========================
# Particles, λ = h/√(2ME)
# ========================

# Neutrons, 
# ----------
M_n = 1.67492749804e-27 # kg
# Thermal neutrons have energies around 0.025
E = eV2J * 0.025

λn = h / sqrt(2*M_n*E) * m2pm
println("Neutrons of energy  $(E/eV2J) eV have a wavelength = $λn pm.")

# Electrons, 
# ----------
M_e = 9.1093837015e-31 # kg
# LEED: E ∈ [20,200] eV
# "ultrafast electron diffraction, especially using ~MeV- class electron guns, such as the system developed at SLAC"
E = eV2J * 150

λe = h / sqrt(2*M_e*E) * m2pm
println("Electrons of energy $(E/eV2J) eV have a wavelength = $λe pm.")


