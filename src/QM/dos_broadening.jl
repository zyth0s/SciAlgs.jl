
function dos_broadening(es,siteweight=missing)
  # DOS or PDOS broadening. Convolutes with Gaussian
  # with support ∈ [-10,10]
  # es is an array with the eigenenergies
  # siteweight is an array with the contribution of a site orbital to each eigenstate
  dos = zeros(10001)
  shape = zeros(20001)
  width = 0.05
  for i in 1:20001
    dummy = 0.002*(i-10001)/(0.5*width)
    shape[i] = 1.0/(1.0 + dummy^2)
  end
  suma = sum(shape)
  # Normalize
  for i in 1:20001
    shape[i] /= suma
  end

  for i in 1:length(es)
    weight = 1.0 # Default is total DOS
    if !ismissing(siteweight)
       weight = siteweight[i] # Partial DOS
    end
    ishape = 1 + trunc(Int,500.0*(es[i]+10.0) + 0.5)
    for j in 1:10001
      dos[j] += weight*shape[j+10001-ishape]
    end
  end

  e_dos = zeros(10001)
  for i in 1:10001
    e_dos[i] = -10. + 0.002*(i-1)
  end
  e_dos,dos
end

