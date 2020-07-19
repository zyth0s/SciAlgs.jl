#!/usr/bin/env julia -q

# Rebuild PyCall using another Python installation

which_python = Dict(
                    "base"  => "/Users/daniel/miniconda3/bin/python",           # conda base
                    "p4env" => "/Users/daniel/miniconda3/envs/p4env/bin/python" # conda p4env (psi4)
                   )

if length(ARGS) < 1
   println("A Python installation must be given. Choose: base | p4env")
   exit()
end

if !haskey(which_python,ARGS[1])
   println("Unknown conda Python installation. Select base | p4env")
   exit()
end

println("Activating $(ARGS[1]) ...")
ENV["PYTHON"] = which_python[ARGS[1]]

using Pkg

Pkg.build("PyCall")

