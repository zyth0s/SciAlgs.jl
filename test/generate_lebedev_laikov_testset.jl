
# Download Lebedev-Laikov dataset from
#    https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
# and generate a Julia @testset

import HTTP
using Formatting: fmt

precision_table = collect(3:2:131)

# No dataset available for some of them
filter!( !in([33,37,39,43,45,49,51,55,57,61,63,67,69,73,
              75,79,81,85,87,91,93,97,99,103,105,109,111,
              115,117,121,123,127,129]),
        precision_table)

open("test_lebedev_laikov.jl", "w") do f

   write(f, "\n")
   write(f, "# Datasets: https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html\n")
   write(f, "\n")
   write(f, "using SciAlgs.NumQuad: lebedev_laikov_spherical \n")
   write(f, "\n")
   write(f, "@testset \"Lebedev-Laikov quadrature\" begin\n") 
   write(f, "\n")

   for i in precision_table

      tolerance = max(exp10(-i), 1e-10)

      write(f, "   # Dataset $i\n")
      write(f, "   dt = [\n")
      r = HTTP.get("https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/lebedev_" * fmt("03d", i) * ".txt")
      #r.status == 200 || @warn "Could not download dataset $i"

      lines = split(String(r.body),"\n")
      text = join(map( x -> x, lines[1:end-1]), ";\n")
      text = text * "]\n"
      write(f, text)
      write(f, "   θ_ref = dt[:,1]\n")
      write(f, "   ϕ_ref = dt[:,2]\n")
      write(f, "   w_ref = dt[:,3]\n")
      write(f, "   n = length(w_ref)\n")
      write(f, "   θ, ϕ, w = lebedev_laikov_spherical(n)\n")
      write(f, "   θ, ϕ = rad2deg.(θ), rad2deg.(ϕ)\n")
      write(f, "   @test isapprox(θ, θ_ref, atol=$tolerance)\n")
      write(f, "   @test isapprox(ϕ, ϕ_ref, atol=$tolerance)\n")
      write(f, "   @test isapprox(w, w_ref, atol=$tolerance)\n")
      write(f, "\n")
   end
   write(f, "end\n")
end
