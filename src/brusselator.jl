
# TODO: Look at 10.1021/acs.jchemed.7b00703

using DifferentialEquations

#    A      -->  X    ; k₁
#   2X + Y  --> 3X    ; k₂
#    B + X  -->  Y + D; k₃
#    X      -->  E    ; k₄
#    d [A]/dt = d[B]/dt = d[D]/dt = d[E]/dt = 0
#    d [X]
#    ----- = k₁[A] + k₂[X]²[Y] - k₃[B][X] - k₄[X] --> dx[3]
#    dt
#    d [Y]
#    ----- = - k₂[X]²[Y] + k₃[B][X] --> dx[4]
#    dt
function brusselator!(dx,x,k,t)
  k₁, k₂, k₃, k₄ = k
  A, B, X, Y, D, E = x
  # A B D E CONSTANT
  dx[1] = 0.0;
  dx[2] = 0.0;
  dx[3] = k₁*A + k₂*X^2*Y - k₃*B*X - k₄*X;
  dx[4] = -k₂*X^2*Y + k₃*B*X;
  dx[5] = 0.0;
  dx[6] = 0.0;
end


#using ParametrizedFunctions
#brusselator = @ode_def begin 
#  dx[1] = 0.0;
#  dx[2] = 0.0;
#  dx[3] = k[1] * x[1] + k[2] * x[3]^2 * x[4] - k[3] * x[2] *x[3] -k[4] * x[3];
#  dx[4] = -k[2] * x[3]^2 * x[4] + k[3] * x[2] * x[3];
#  dx[5] = 0.0;
#  dx[6] = 0.0;
#end k₁ k₂ k₃ k₄

k = [1.0,1.0,1.0,1.0]
x0 = zeros(Float64,6)
x0[1] = 1.0; # mol / L
x0[2] = 5.5; # mol / L
x0[3] = 1.0; # mol / L
x0[4] = 1.0; # mol / L
x0[5] = 0.0; # mol / L
x0[6] = 0.0; # mol / L
tspan = (0.0,100.0)

prob = ODEProblem(brusselator!,x0,tspan,k)
x = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

using Plots
plot(x, title="k=[$(k[1]) $(k[2]) $(k[3]) $(k[4])],[ ]=[1,5.5,1,1,0,0]", xaxis="[X]",yaxis="[Y]")


#grid;
#print ("@x0(2) = 5.5.png","-dpng");
#pause

