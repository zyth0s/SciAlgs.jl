
using DifferentialEquations

function brusselator!(dx,x,k,t)
  dx[1] = 0.0;
  dx[2] = 0.0;
  dx[3] = k[1] * x[1] + k[2] * x[3]^2 * x[4] - k[3] * x[2] *x[3] -k[4] * x[3];
  dx[4] = -k[2] * x[3]^2 * x[4] + k[3] * x[2] * x[3];
  dx[5] = 0.0;
  dx[6] = 0.0;
end

k = [1.0,1.0,1.0,1.0]
x0 = zeros(Float64,6)
x0[1] = 1.0;
x0[2] = 5.5;
x0[3] = 1.0;
x0[4] = 1.0;
x0[5] = 0.0;
x0[6] = 0.0;
tspan = (0.0,100.0)

prob = ODEProblem(brusselator!,x0,tspan,k)
x = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

using Plots
plot(x, title="k=[1,1,1,1],[ ]=[1,5.5,1,1,0,0]", xaxis="[X]",yaxis="[Y]")


#grid;
#print ("@x0(2) = 5.5.png","-dpng");
#pause

