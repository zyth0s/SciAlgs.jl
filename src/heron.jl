#!/usr/bin/env julia
#

# Analysis of Heron's algorithm to
# compute √x. 
# Theory: "Numerical Analysis" by Ridley Scott

using Formatting: printfmt

function heron(x,y)
  x = .5*(x + y/x)
end

function errheron(x,y)
  println("Value       Error      Rel. Error")
  println("-------  -----------  ------------")
  for i=1:5
    x = heron(x,y)
    err = x - sqrt(y)
    relerr = (x/sqrt(y)) - 1.
    printfmt("{:7.4f}  {:+7.4e}  {:+7.4e}\n",x, err, relerr)
  end
end

x = 1. # Starting point
y = 2.

errheron(x,y)
