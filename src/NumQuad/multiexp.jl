
using Formatting: printfmt
using LinearAlgebra

# Adapted from supplementary routines for the book
# "Spectral methods in Chemistry and Physics"

function multiexp(N)
  #nintervals: #of intervals; npts: # of points per interval [rmin,rmax]
  rmin = 0
  rmax = 1
  npts = 1000 - round(Int,500*((20-N)/20)) # minimum effort (ad hoc)
  nintervals = npts
  ntot = nintervals*npts
  @assert ntot/N > 50
  # a suitable quadrature rule on (rmin,rmax)
  pwmd = multidomain_quadrature(nintervals,npts,rmin,rmax)
  # multi-domain matrix of nodes p and weights w
  p = pwmd[:,1]
  w = pwmd[:,2]
  # Discrete Stieltjes algorithm for the calculation of recurrence coefficients
  # αₖ and βₖ for weight function ω(r) = ln²(r)
  # Sections 2.1-2.2 in 10.1137/0903018
  # πₖ₊₁ (r) = (r-αₖ)πₖ(r) - βₖπₖ₋₁(r);   π₋₁(r) = 0,  π₀(r) = 1
  #
  #      ∫ rπₖ²(r) dλ(r)      s2
  # αₖ = ----------------- ≈  --.        k = 0,1,2,...
  #      ∫ πₖ²(r)  dλ(r)      s1
  #
  #      ∫ πₖ²(r)   dλ(r)     h[k]
  # βₖ = ----------------- ≈ --------,   k = 1,2,3,...; β₀ = ∫ dλ(r)
  #      ∫ πₖ₋₁²(r) dλ(r)     h[k-1]
  #
  # Discretization => dλ(r) = ω(r) dr on (rmin,rmax)
  ω = @. log(p)^2
  #
  P0 = ones(ntot) # π₋₁(r) = 0
  h  = zeros(N)
  α  = zeros(N)
  β  = zeros(N)
  # First two integrals; (2.4) in 10.1137/0903018
  s1 = sum(w .*       ω)  # here, π₀(r) = 1
  s2 = sum(w .* (p .* ω)) # " " "
  k = 1
  h[k] = s1
  α[k] = s2/s1
  β[k] = 0 # β₀ = ∫ dλ(r)
  # Norm and α₁ and β₁
  #println("Norm: ", h[k])
  #println(α[k], " ", β[k])
  # Polynomial π₁(r) = (r-α₀)1
  P1 = p .- α[k]
  # Next two integrals
  s1 = sum(w .*       (ω .* (P1.^2)))
  s2 = sum(w .* (p .* (ω .* (P1.^2))))
  if N > 1
    k = 2
    h[k] = s1
    α[k] = s2/s1
    β[k] = h[k]/h[k-1]
  end
  if N > 2
    # Norm and α₂ and β₂
    #println("Norm: ", h[k])
    #println(α[k], " ", β[k])
    for k in 3:N
      Pma = @. p - α[k-1]
      # Recurrence for the next polynomial
      P2 = @. Pma * P1 - β[k-1]*P0
      s1 = sum(w .*       (ω .* (P2.^2)))
      s2 = sum(w .* (p .* (ω .* (P2.^2))))
      α[k] = s2/s1
      h[k] = s1
      β[k] = h[k]/h[k-1]
      P0 = P1
      P1 = P2
    end
  end

  # Calculate the quadrature points and weights with Golub-Welsch
  rtb = sqrt.(β[2:N])
  J = SymTridiagonal(α,rtb) # Jacobi matrix
  λ, f = eigen(J)
  wt = h[1] * f[1,:].^2
  #println("MultiExp quadrature points and weights")
  #for n in 1:N
  #  printfmt(" {:2d} {:26.10f} {:26.10e}\n", n, λ[n],wt[n])
  #end
  λ, wt
end

# Build a composite quadrature; the rule for each interval
# should not include the endpoints to avoid duplication.
function multidomain_quadrature(nintervals,npts,rmin,rmax)
  ntot = nintervals * npts
  dr = (rmax-rmin)/ntot
  pwmd = Array{Float64,2}(undef,npts,2)
  for i in 1:nintervals
    a = rmin + (i-1)*npts*dr
    b = a + npts*dr
    pw = fejer2(a,b,npts)
    if i == 1
      pwmd = pw
    else
      pwmd = vcat(pwmd,pw)
    end
  end
  sum(pwmd[:,2]) ≈ (rmax-rmin) || error("Weights do not sum up to $(rmax-rmin)")
  pwmd
end

# Fejer quadrature rule
# Nodes of Tₙ (1st kind Chebyshev polynomials)
function fejer2(a,b,N)
  n = N:-1:1
  m = 1:floor(Int,N/2)
  θ = @. (n-0.5)*π / N
  p = cos.(vec(θ))
  w = zeros(N)
  for k in N:-1:1
    s = sum(@. cos(2m*θ[k])/(4(m^2)-1))
    w[k] = 2(1-2s)/N
  end
  # p ∈ [-1,1] -> ps ∈ [a,b]
  r1 = (b-a)/2
  r2 = (a+b)/2
  ps = @. r1*p + r2
  ws = r1*w
  sum(ws) ≈ (b-a) || error("Weights do not sum up to $(b-a)")
  [ps vec(ws)]
end


function test_multiexp()
   # Taken from Table 1 of 10.1002/jcc.10211
   # Also https://rsc.anu.edu.au/~pgill/multiexp.php?n=5
   println("Testing nquad=1")
   r,w = multiexp(1)
   @assert isapprox(r[1],0.1250000000)
   println("Testing nquad=2")
   r,w = multiexp(2)
   @assert isapprox(r[1],0.0598509925)
   @assert isapprox(r[2],0.4536625210)
   println("Testing nquad=3")
   r,w = multiexp(3)
   @assert isapprox(r[1],0.0362633111)
   @assert isapprox(r[2],0.2731486024)
   @assert isapprox(r[3],0.6537110896)
   println("Testing nquad=4")
   r,w = multiexp(4)
   @assert isapprox(r[1],0.0246451318)
   @assert isapprox(r[2],0.1831933310)
   @assert isapprox(r[3],0.4610171077)
   @assert isapprox(r[4],0.7655906466)
   println("Testing nquad=5")
   r,w = multiexp(5)
   @assert isapprox(r[1],0.0179624485)
   @assert isapprox(r[2],0.1317184306)
   @assert isapprox(r[3],0.3395971926)
   @assert isapprox(r[4],0.5945982935)
   @assert isapprox(r[5],0.8320575996)
   println("Testing nquad=6")
   r,w = multiexp(6)
   @assert isapprox(r[1],0.0137303290)
   @assert isapprox(r[2],0.0994431475)
   @assert isapprox(r[3],0.2596678762)
   @assert isapprox(r[4],0.4685229897)
   @assert isapprox(r[5],0.6874245835)
   @assert isapprox(r[6],0.8742037763)
   println("Testing nquad=8")
   r,w = multiexp(8)
   @assert isapprox(r[1],0.0088308098)
   @assert isapprox(r[2],0.0626470137)
   @assert isapprox(r[3],0.1653937470)
   @assert isapprox(r[4],0.3076475309)
   @assert isapprox(r[5],0.4738811643)
   @assert isapprox(r[6],0.6449256028)
   @assert isapprox(r[7],0.8005576159)
   @assert isapprox(r[8],0.9222020825)
   println("Testing nquad=10")
   r,w = multiexp(10)
   @assert isapprox(r[1] ,0.0061869147)
   @assert isapprox(r[2] ,0.0431849645)
   @assert isapprox(r[3] ,0.1143932978)
   @assert isapprox(r[4] ,0.2157443263)
   @assert isapprox(r[5] ,0.3400163758)
   @assert isapprox(r[6] ,0.4777530668)
   @assert isapprox(r[7] ,0.6181540046)
   @assert isapprox(r[8] ,0.7500277506)
   @assert isapprox(r[9] ,0.8627655156)
   @assert isapprox(r[10],0.9472975116)
   println("Testing nquad=15")
   r,w = multiexp(15)
   @assert isapprox(r[1] ,0.0031568454)
   @assert isapprox(r[2] ,0.0214217428)
   @assert isapprox(r[3] ,0.0567152146)
   @assert isapprox(r[4] ,0.1082790024)
   @assert isapprox(r[5] ,0.1744898369)
   @assert isapprox(r[6] ,0.2530501060)
   @assert isapprox(r[7] ,0.3411169468)
   @assert isapprox(r[8] ,0.4354309056)
   @assert isapprox(r[9] ,0.5324530702)
   @assert isapprox(r[10],0.6285097879)
   @assert isapprox(r[11],0.7199410768)
   @assert isapprox(r[12],0.8032477673)
   @assert isapprox(r[13],0.8752323216)
   @assert isapprox(r[14],0.9331298792)
   @assert isapprox(r[15],0.9747402975)
   println("Testing nquad=20")
   r,w = multiexp(20)
   @assert isapprox(r[1] ,0.0019241239)
   @assert isapprox(r[2] ,0.0128189043)
   @assert isapprox(r[3] ,0.0338360585)
   @assert isapprox(r[4] ,0.0647886177)
   @assert isapprox(r[5] ,0.1051527594)
   @assert isapprox(r[6] ,0.1541448603)
   @assert isapprox(r[7] ,0.2107591088)
   @assert isapprox(r[8] ,0.2737989508)
   @assert isapprox(r[9] ,0.3419087328)
   @assert isapprox(r[10],0.4136071132)
   @assert isapprox(r[11],0.4873224017)
   @assert isapprox(r[12],0.5614294486)
   @assert isapprox(r[13],0.6342874604)
   @assert isapprox(r[14],0.7042779985)
   @assert isapprox(r[15],0.7698423653)
   @assert isapprox(r[16],0.8295175821)
   @assert isapprox(r[17],0.8819702175)
   @assert isapprox(r[18],0.9260275322)
   @assert isapprox(r[19],0.9607063554)
   @assert isapprox(r[20],0.9852482390)
end
##          weights
#1         2.0000000000
#2         1.6691361082
#          0.3308638918
#3         1.3638303836
#          0.5658154596
#          0.0703541567
#4         1.1330156422
#          0.6612166786
#          0.1857929500
#          0.0199747293
#5         0.9588537970
#          0.6830020585
#          0.2815660272
#          0.0695856412
#          0.0069924762
#6         0.8247373524
#          0.6701662035
#          0.3457959549
#          0.1269936018
#          0.0294555188
#          0.0028513686
#8         0.6343476124
#          0.6078905783
#          0.4036485580
#          0.2184654657
#          0.0959319490
#          0.0320637571
#          0.0069997943
#          0.0006522851
#10        0.5075100632
#          0.5377370883
#          0.4101581499
#          0.2686821571
#          0.1544536152
#          0.0768964762
#          0.0319269841
#          0.0102578456
#          0.0021782820
#          0.0001993384
#15        0.3253137565
#          0.3961243562
#          0.3593522229
#          0.2926536989
#          0.2219588508
#          0.1582926857
#          0.1061466261
#          0.0665501730
#          0.0385929320
#          0.0203475527
#          0.0094938596
#          0.0037511449
#          0.0011607732
#          0.0002398218
#          0.0000215456
#20        0.2308490189
#          0.3027886725
#          0.2986810329
#          0.2678766357
#          0.2274158666
#          0.1852655647
#          0.1455612505
#          0.1104305753
#          0.0808087449
#          0.0568766276
#          0.0383323430
#          0.0245782243
#          0.0148581882
#          0.0083614715
#          0.0043003770
#          0.0019659640
#          0.0007640520
#          0.0002333787
#          0.0000477502
#          0.0000042614
