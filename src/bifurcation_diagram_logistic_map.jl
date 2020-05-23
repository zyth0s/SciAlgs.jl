
using Plots
using Random

Random.seed!(0)
function logistic_map(x,r)
   r*x*(1-x)
end

rrange = range(0.7,stop=4,length=10000)
xarray = []
for r in rrange
   x = rand()
   for i in 1:2050
      x = logistic_map(x,r)
   end
   push!(xarray,x)
end
scatter(rrange, xarray, 
        legend=false,
        xlabel = "Growth factor, r",
        ylabel = "Far future population",
        markershape = :circle, 
        markercolor = :black, 
        ms = 0.2)
