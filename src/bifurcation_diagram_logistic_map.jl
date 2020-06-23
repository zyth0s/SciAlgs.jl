
using Random
import PyPlot
pyplt = PyPlot
mpl = PyPlot.matplotlib
pyplt.matplotlib.style.reload_library()
pyplt.matplotlib.style.use("5in_color")
mpl.use(backend="Qt5Agg")

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

fig = pyplt.figure()
#fig.set_size_inches([15/2.54,10/2.54],forward=true)
fig.set_dpi(260)
ax = fig.add_subplot(111)
# Plot data
ax.scatter(rrange,xarray, label=nothing, s=1,c=:black)
ax.set_ylim(0.0, 1.0)
ax.set_xlim(1.0, 4.1)
# Edit the major and minor tick locations
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
# Set axis labels
ax.set_xlabel("Growth factor, r")
ax.set_ylabel("Far future population")
fig.tight_layout(pad=0.1)
pyplt.plot()
pyplt.savefig("figures/logistic_map.pdf")
