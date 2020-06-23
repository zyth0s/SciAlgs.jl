
import PyPlot
pyplt = PyPlot
mpl = PyPlot.matplotlib

pyplt.matplotlib.style.reload_library()
pyplt.matplotlib.style.use("sci")

fig, axs = pyplt.subplots(ncols=2,nrows=1,sharex=false,sharey=true,
                          figsize=(7,5))

# Create some data
x = range(0, stop=π, length=200)
y1 = sin.(x)
y2 = 1.5sin.(x)
y3 = 2sin.(x)
for ax in axs
   # Plot data
   ax.plot(x, y1, label="A = 1")
   ax.plot(x, y2, label="A = 1.5")
   ax.plot(x, y3, label="A = 2")
   # Set axis labels
   ax.set_xlabel("Wavelength (nm)", labelpad=10)
   #ax.set_xlabel(r"$\mathregular{\lambda}$ (nm)", labelpad=10)
   ax.set_ylabel("Absorbance (O.D.)", labelpad=10)
   # Edit the major and minor tick locations
   ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
   ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
   ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
   ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
   # Create new axes object by cloning the y-axis of the first plot
   ax2 = ax.twiny()
   # Edit the tick parameters of the new x-axis
   ax2.xaxis.set_tick_params(which="major", size=10, width=2, direction="in")
   ax2.xaxis.set_tick_params(which="minor", size=7, width=2, direction="in")
   # Function to convert energy (eV) to wavelength (nm)
   function E_to_WL(E)
       return 1240 ./ E
    end
   # Add ticks manually to energy axis
   ax2.xaxis.set_major_locator(mpl.ticker.FixedLocator(E_to_WL(range(1.5, stop=3.0, length=4))))
   ax2.xaxis.set_minor_locator(mpl.ticker.FixedLocator(E_to_WL(range(1.4, stop=3.2, length=19))))
   # Add tick labels manually to energy axis
   ax2.set_xticklabels(["1.5", "2.0", "2.5", "3.0"])
   # Add energy axis label
   ax2.set_xlabel("Energy (eV)", labelpad=10)
   # Set energy axis limits
   ax2.set_xlim(370, 930)
   # Add legend - loc is a tuple specifying the bottom left corner
   # Add legend to plot
   ax.legend(bbox_to_anchor=(0.8, 0.35), loc=1, frameon=false, fontsize=16)
end

# Save plot
fig.tight_layout(pad=0.1)
#pyplt.show()
pyplt.savefig("pyplot_subplots.pdf")
