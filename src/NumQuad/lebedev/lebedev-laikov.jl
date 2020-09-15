
#using Libdl
#push!(Libdl.DL_LOAD_PATH, @__DIR__)

const LIBLEBEDEV    = joinpath(@__DIR__, "lebedev-laikov.so")
LEBEDEVSOURCE = joinpath(@__DIR__, "lebedev-laikov.c")
LEBEDEVHEADER = joinpath(@__DIR__, "lebedev-laikov.h")

# Compile Lebedev-Laikov routine as a shared library
if !isfile(LIBLEBEDEV) ||
   mtime(LEBEDEVSOURCE) > mtime(LIBLEBEDEV) ||
   mtime(LEBEDEVHEADER) > mtime(LIBLEBEDEV)

   run(`gcc -shared -fPIC -o $LIBLEBEDEV $LEBEDEVSOURCE`)
end

function lebedev_laikov_wrapper(nang_pts::Int)

   valid_orders = [6,    14,   26,   38,   50,   74,   86,   110,  146,  170,  194,
                   230,  266,  302,  350,  434,  590,  770,  974,  1202, 1454, 1730,
                   2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]

   # Interpolation to closest valid order; if unique choose smallest
   valid_nang_pts = valid_orders[argmin(abs.(nang_pts .- valid_orders))]
   if valid_nang_pts != nang_pts 
      @debug "Number of angular grid points changed: $nang_pts -> $valid_nang_pts"
   end

   @assert valid_nang_pts ∈ valid_orders

   x       = zeros(valid_nang_pts)
   y       = zeros(valid_nang_pts)
   z       = zeros(valid_nang_pts)
   weights = zeros(valid_nang_pts)

   @static if VERSION >= v"1.5.0"

      @ccall LIBLEBEDEV.ld_by_order(
         valid_nang_pts :: Int,
         x              :: Ref{Float64},
         y              :: Ref{Float64},
         z              :: Ref{Float64},
         weights        :: Ref{Float64}
      )::Cint
   else

      ccall((:ld_by_order,"lebedev-laikov.so"),
         Int,
         (Int, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
         valid_nang_pts, x, y, z, weights)
   end

   @assert sum(weights) ≈ 1 # as a test?

   [x y z], weights
end

using SciAlgs: cart2spherical

function lebedev_laikov_spherical(n)
   points, weights = lebedev_laikov_wrapper(n)
   θ, ϕ = cart2spherical(points)
   θ, ϕ, weights
end
