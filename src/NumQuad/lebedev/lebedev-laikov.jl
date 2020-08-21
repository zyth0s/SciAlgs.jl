
#using Libdl
#push!(Libdl.DL_LOAD_PATH,pwd())

# TODO: test it against: https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html

# Compile Lebedev-Laikov routine as a shared library
if !isfile("lebedev-laikov.so") ||
   mtime("lebedev-laikov.c") > mtime("lebedev-laikov.so") ||
   mtime("lebedev-laikov.h") > mtime("lebedev-laikov.so")

   run(`cc -shared -o lebedev-laikov.so lebedev-laikov.c`)
end

const LIBLEBEDEV = "lebedev-laikov.so"

function lebedev_laikov_wrapper(nang_pts::Int)

   valid_orders = [6,    14,   26,   38,   50,   74,   86,   110,  146,  170,  194,
                   230,  266,  302,  350,  434,  590,  770,  974,  1202, 1454, 1730,
                   2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]

   # Interpolation to closest valid order; if unique choose smallest
   valid_nang_pts = valid_orders[argmin(abs.(nang_pts .- valid_orders))]
   if valid_nang_pts != nang_pts 
      @debug "Number of angular grid points changed: $nang_pts -> $valid_nang_pts"
   end
   check_nang_pts = 0

   @assert valid_nang_pts ∈ valid_orders

   x       = zeros(valid_nang_pts)
   y       = zeros(valid_nang_pts)
   z       = zeros(valid_nang_pts)
   weights = zeros(valid_nang_pts)

   @static if VERSION >= v"1.5.0"

      check_nang_pts =
               @ccall LIBLEBEDEV.leb_order_routing(
                  valid_nang_pts :: Int,
                  x              :: Ref{Float64},
                  y              :: Ref{Float64},
                  z              :: Ref{Float64},
                  weights        :: Ref{Float64}
               )::Int
   else

      check_nang_pts =
         ccall((:leb_order_routing,"lebedev-laikov.so"),
            Int,
            (Int, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
            valid_nang_pts, x, y, z, weights)
   end
   @assert check_nang_pts == valid_nang_pts

   @assert sum(weights) ≈ 1 # as a test?

   x, y, z, weights
end

# Test
n = 11 # not a valid order -> 6th order
x,y,z,w = lebedev_laikov_wrapper(n)
@info "$(length(w)) Lebedev-Laikov grid points" x' y' z' w'
println()
