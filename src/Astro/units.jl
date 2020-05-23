
# Arclengths
# ===========
function deg2min(degree)
   degree * 60
end

function sec2deg(seconds)
   seconds/3600
end

# Lengths
# ==========
au = 149597870700 # atronomical unit (au) [m]

# Parsec is the distance at which 1 au subtends 1 arcsec.
#         /| 1 arcsec
#       /  |
#     /    |  x au = 1 pc
#   /      |
# /________|
#     1 au

pc = au/tand(sec2deg(1)) # parsec (pc) [m]


function distance_from_parallax_au(p)
   # p is the parallax in arcsec
   # the returned distance is in au units
   1/p
end


# Examples
# ========
# 1) 61 Cygni, p = 0.3 arcsec
#distance_from_parallax(0.3)
# 2) Proxima Centauri, p = 0.762 arcsec
isapprox(distance_from_parallax_au(0.762), 1.31,atol=1e-2) || error("Wrong distance estimation")

