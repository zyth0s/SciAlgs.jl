
using Plots

# Inspired by doi:10.1021/acs.jchemed.9b00542
function linear_sweep_voltammogram()
  ###################################
  ## Section 1: Simulation inputs ##
  ###################################
  ### Constants ###
  F = 96485 # Faraday’s constant [C/mol]
  R = 8.314 # Gas constant [C V/(mol K)]
  T = 298 # Temperature [K]
  ### Electrochemical parameters ###
  # Potential waveform
  E_start = 0 # Start potential [V] {default: 0}
  E_end = -0.8 # End potential [V] {default: -0.8}
  E₀ = -0.4 # Standard potential [V] {default: -0.4}
  ν = 10 # Scan rate [mV/s] {default: 10}
  # Redox mediator
  D = 1E-6 # Diffusion coefficient [cm^2/s] {default: 1E-6}
  c_ox = 1.0 # Bulk concentration of oxidized species [M] {default: 1}
  c_red = 1E-3 # Bulk concentration of reduced species [M] {default: 1E-3}
  n = 1 # Number of electrons transferred during reaction {default: 1}
  α = 0.5 # Transfer coefficient (unitless) {default: 0.5}
  # Electrode
  A = 1 # Electrode area [cm²] {default: 1}
  k₀ = 1E-1  # Heterogeneous rate constant [m/s] {default: 1E-1}
  ### Finite difference parameters ###
  npts_x = 100 # Number of mesh points {default: 100}
  npts_t = 500 # Number of time steps {default: 500}
  #####################################################
  # Section 2: Time 0 / Initial Setup and Time Step ##
  #####################################################
  # Discretize x
  total_x = 0.1 # Maximum x value to simulate [cm]
  Δx = total_x/npts_x # x interval
  x = range(0,stop=total_x,length=npts_x)
  # Discretize t
  total_t = abs(E_end - E_start)/(ν/1000) # Maximum t value to simulate (1 LSV)
  Δt = total_t/npts_t # Time interval
  t = range(0, stop = total_t,length=npts_t)
  # Set uniform concentration along discretized x
  cox_x     = fill(c_ox, npts_x,npts_t)
  cred_x    = fill(c_red, npts_x,npts_t)
  cred_temp = zeros(npts_x,npts_t)
  cox_temp = zeros(npts_x,npts_t)
  # Setup empty matrices for time-dependent quantities that will be
  # filled in as the simulation progresses
  E_t    = zeros(npts_t)
  kred_t = zeros(npts_t)
  iₜ    = zeros(npts_t)
  # Calculate initial potential
  E_curr = E_start # E_curr = Current potential (single value)
  E_t[1] = E_curr # E_t = Potential waveform (matrix of values)
  # Calculate formal potential according to the Nernst equation
  E_eq = E₀ - ((R*T)/(n*F) * log(c_ox/c_red))
  # Calculate initial rate constant according to the Butler-Volmer equation
  η = E_curr - E_eq # Overpotential (V)
  kred_curr = k₀*exp((-α*n*F*η)/(R*T))
  # Concentration in box 0 updated due to chemical reaction
  cred_x[1,1] = cred_x[1] + (Δt*(kred_curr * cox_x[1,1])) # Equation 17
  cox_x[1,1] = cox_x[1] - (Δt*(kred_curr * cox_x[1,1])) # Equation 20
  # Current is calculated based on change in concentration
  iₜ[1] = -n*F*(Δx*(cred_x[1,1] - c_red) / Δt)
  ################################################
  ## Section 3: Time > 0: Simulation Continues ##
  ################################################
  for j = 2:npts_t # Iterate over number of time points
    ### Pull concentration profiles from end of previous time step ###
    cred_x[:,j] = cred_x[:,j-1]
    cox_x[:,j] = cox_x[:,j-1]
    ### Update concentration due to diffusion: Equation 26 ###
    for i = 1:npts_x # Iterate over number of boxes
      # Note: a temporary variable is used here to make sure that the
      # concentration for this time step does not change as it is being
      # used in the calculation
      # In the first box (i=1), we use a modified diffusion equation since
      # there is no box with index (i-1)
      if i == 1
        cred_temp[i,j] = cred_x[i,j] + (D*Δt/(Δx^2))* ((cred_x[i+1,j] - cred_x[i,j]))
        cox_temp[i,j] = cox_x[i,j] + (D*Δt/(Δx^2))* ((cox_x[i+1,j] - cox_x[i,j]))
        # In the last box (i=max), we apply our boundary condition of
        # concentration = bulk concentration
      elseif i == npts_x
        cred_temp[i,j] = c_red
        cox_temp[i,j] = c_ox
        # In all other boxes, we use the diffusion expression
        # (Equation 26) as normal
      else
        cred_temp[i,j] = cred_x[i,j] + (D*Δt/Δx^2)*((cred_x[i+1,j] - 2*cred_x[i,j] + cred_x[i-1,j]))
        cox_temp[i,j] = cox_x[i,j] + (D*Δt/Δx^2)*((cox_x[i+1,j] - 2*cox_x[i,j] + cox_x[i-1,j]))
      end
    end
    cred_x[:,j] = cred_temp[:,j]
    cox_x[:,j] = cox_temp[:,j]
    ### Update potential: Equation 12 ###
    E_curr = E_curr - Δt*((ν/1000))
    E_t[j] = E_curr
    ### Update rate constant: Equation 14 ###
    η = E_curr - E_eq # Overpotential (V)
    kred_curr = k₀* exp((-α*n*F*η)/(R*T))
    kred_t[j] = kred_curr
    ### Chemical reaction in box 0: Equations 17 & 20 ###
    cred_x[1,j] = cred_x[1,j-1] + (Δt*(kred_curr * cox_x[1,j-1]))
    cox_x[1,j] = cox_x[1,j-1] - (Δt*(kred_curr * cox_x[1,j-1]))
    ### Calculate current: Equation 22 ###
    iₜ[j] = -n*F*(Δx*(cred_x[1,j] - cred_x[1,j-1]) / Δt)
  end

  ###########################
  ## Section 4: Figures ##
  ###########################
  # Settings applying to the appearance of all graphs
  FontSize = 18
  DataLineWidth = 2
  PlotLineWidth = 2
  # Indices of times to plot concentration profiles
  index_time2 = argmin(iₜ) # Locate time index when current is largest
  index_time1 = trunc(Int,0.9*index_time2) # Select a time point prior to peak
  index_time3 = trunc(Int,1.1*index_time2) # Select a time point after peak
  ## Construct the megafigure
  #figure(’units’,’normalized’,’outerposition’,[0 0 1 0.9])
  #subplot(2,4,[1,2,5,6]) # Voltammogram
  p1 = plot(E_t,iₜ,
    xlabel = "Potential (V)",
    ylabel = "Current (A)",
    title = "Voltammogram"
    ) # ’linewidth’,2)
  #axis([E_end E_start 1.25*min(iₜ) -0.1*min(iₜ)])
  #xL = get(gca,’XLim’)
  #yL = get(gca,’YLim’)
  #ytickformat(’%,.1f’)
  #line([E_t(index_time1) E_t(index_time1)],yL,’linewidth’,DataLineWidth-0.5,’color’,[0.8500
  #0.3250 0.0980])
  #line([E_t(index_time2) E_t(index_time2)],yL,’linewidth’,DataLineWidth-0.5,’color’,[0.9290
  #0.6940 0.1250])
  #line([E_t(index_time3) E_t(index_time3)],yL,’linewidth’,DataLineWidth-0.5,’color’,[0.4940
  #0.1840 0.5560])
  #legend(’i(E)’,’Point 1’,’Point 2’,’Point 3’,’location’,’best’)
  #set(gca,’FontSize’,FontSize,’LineWidth’,PlotLineWidth)
  #subplot(2,4,3) # Potential waveform
  p2 = plot(t,E_t,
    xlabel = "Time (s)",
    ylabel = "Potential (V)",
    title = "Potential Waveform",
    ) #’linewidth’,DataLineWidth)
  #xL = get(gca,’XLim’)
  #line(xL,[E₀, E₀],’linestyle’,’--’,’linewidth’,DataLineWidth,’color’,’k’)
  #line(xL,[E_eq, E_eq],’linestyle’,’:’,’linewidth’,DataLineWidth,’color’,’k’)
  #ytickformat(’%,.1f’)
  #legend(’E(t)’,’E_{0}’,’E_{eq}’,’location’,’best’)
  #set(gca,’FontSize’,FontSize,’LineWidth’,PlotLineWidth)
  #axis([-Inf Inf -Inf Inf])
  #subplot(2,4,4) # Rate constant
  p3 = plot(t,kred_t,
    xlabel = "Time (s)",
    ylabel = "k_{red} (cm/s)",
    title = "Rate constant"
    ) # ’linewidth’,DataLineWidth)
  #set(gca,’FontSize’,FontSize,’LineWidth’,PlotLineWidth)
  #subplot(2,4,7) # Normalized concentration profile of c_red
  p4 = plot(1000 .* x,cred_x[:,index_time1]./c_ox)
  plot!(p4, 1000 .* x,cred_x[:,index_time2]./c_ox)
  plot!(p4, 1000 .* x,cred_x[:,index_time3]./c_ox, lw = 2)
  #set(h(1),’color’,[0.8500 0.3250 0.0980])
  #set(h(2),’color’,[0.9290 0.6940 0.1250])
  #set(h(3),’color’,[0.4940 0.1840 0.5560])
  #xL = get(gca,’XLim’)
  #line(xL,[c_ox c_ox],’linestyle’,’--’,’color’,’k’)
  #xlabel(’x (mm)’)
  #ylabel(’c_R/c_{O (bulk)}’)
  #set(gca,’FontSize’,FontSize,’LineWidth’,PlotLineWidth)
  #axis([0 1000*0.5*max(x) 0.5*c_red 1.2*max(max(cred_x/c_ox))])
  #legend(’Point 1’,’Point 2’,’Point 3’,’c_{ox (bulk)}’,’location’,’southeast’)
  #title(’Concentration Profile’)
  #subplot(2,4,8) # Normalized concentration profile of c_ox
  p5 = plot(1000 .*x,cox_x[:,index_time1]./c_ox)
  plot!(p5, 1000 .*x,cox_x[:,index_time2]./c_ox)
  plot!(p5, 1000 .*x,cox_x[:,index_time3]./c_ox,lw=2)
  #set(h(1),’color’,[0.8500 0.3250 0.0980])
  #set(h(2),’color’,[0.9290 0.6940 0.1250])
  #set(h(3),’color’,[0.4940 0.1840 0.5560])
  #xL = get(gca,’XLim’)
  #line(xL,[c_ox c_ox],’linestyle’,’--’,’color’,’k’)
  #xlabel(’x (mm)’)
  #ylabel(’c_O/c_{O (bulk)}’)
  #ytickformat(’%,.1f’)
  #set(gca,’FontSize’,FontSize,’LineWidth’,PlotLineWidth)
  #axis([0 1000*0.5*max(x) 0 1.2*c_ox])
  #legend(’Point 1’,’Point 2’,’Point 3’,’c_{ox (bulk)}’,’location’,’southeast’)
  #title(’Concentration Profile’)
  l = @layout [ a [b c; d e]]
  display(plot(p1,p2,p3,p4,p5,layout = l))
end

# Simulation of a linear-sweep voltammetric curve with mass transport and ohmic effects
#   Julia implementation of doi:10.1021/ed077p100
#   by Benedetto P. Bozzini
#   Dipartimento di Scienza dei Materiali - Universit_ di Lecce
#   v. Arnesano, I-73100 Lecce (Italy)
#   benedetto. bozzini@polimi.it & bozzini@axpmat.unile.it
#   version 05/10/1998
# Constants and algorithm of this version are an example of:
# codeposition of Co and H₂ with mass transport limitations and ohmic drop
# the program can be tailored to your needs by substituting the suitable electrochemical constants
function linear_sweep_voltammogram_lifelike()
  # input
  # 1. define the array "etac" (mV) of voltages within which the true kinetic overvoltage is to be
  # searched
  # e.g. etac=[-1500:.1:-10]
  # 2. define the value etapd (mV) the voltage imposed in the potentiodynamic scan
  # output
  # the program outputs the following value:
  # 1. MARK: the true kinetic overvoltage
  # 2. INMARK: the c.d. of Co deposition
  # 3. IHMARK: the c.d. of H₂ deposition initialisation
  # (clears arrays whose lenghts can be altered in successive runs)
  # model constants
  R = 8.31
  # R [J/K/mol]: gas constant
  T = 300
  # T [K]: absolute temperature
  F = 96500
  # F [C/eq]: Faraday's constant
  i₀n = 2.0e-3
  # i₀n [mA/cm2]: exchange current density for the metal reaction
  # or, in general, for the first electrochemical reaction of interest
  i₀h = 6.3e-3
  # i₀h [mA/cm2]: exchange current density for the hydrogen reaction
  # or, in general, for the second electrochemical reaction of interest
  iln = -85
  # iln [mA/cm2]: limiting current density for the metal reaction
  # or, in general, for the first electrochemical reaction of interest
  ilh = -0.140
  # ilh [mA/cm2]: limiting current density for the hydrogen reaction
  # or, in general, for the second electrochemical reaction of interest
  an = 0.65
  # an [1]: transfer coefficient for the metal reaction
  # or, in general, for the first electrochemical reaction of interest
  ah = 0.43
  # ah [1]: transfer coefficient for the hydrogen reaction
  # or, in general, for the second electrochemical reaction of interest
  etamix = -212.2321
  # etamix [mV]: mixed potential value of starting potential
  etaeqn = -250
  # etaeqn [mV]: equilibrium potential value of the metal semicell
  # or, in general, for the first electrochemical reaction of interest
  etaeqh = -148.7
  # etaeqh [mV]: equilibrium potential value of the hydrogen semicell
  # or, in general, for the second electrochemical reaction of interest
  res = 1.000
  # res [ohm=mV/mA]: resistance
  area = 4
  # area [cm2]: electrode surface area
  f = F/R/T/1000
  # f: an auxialiary quantity
  # definition of the array to be searched for the true kinetic overvoltage
  etac = range(-910,stop=0,step=1) #[-910:1:0]
  # definition of the array of scanned potentiodynamic voltages
  escan = range(-910,stop=-10,step=1)

  amark   = zeros(length(escan))
  ainmark = zeros(length(escan))
  aihmark = zeros(length(escan))
  aitmark = zeros(length(escan))
  # most external loop (= linear-sweep voltammetric loop) starts here
  # scanning potentiodynamic voltage
  for ibonji = 1:length(escan)
    etapd = escan[ibonji]
    # a running check of runtime: measures the percentage of run currently performed
    runamount = ibonji/length(escan)*100
    #
    etan = etamix .+ etac .- etaeqn
    etah = etamix .+ etac .- etaeqh
    # computes actual overvoltages for each electrodic reaction
    # 1. etac is scanned in the LSV from the instrumental 0 to the instumental final value
    # 2. the actual electrode voltage during scanning is etamix+etac
    # 3. the actual overvoltage of react. A is (actual electrode voltage)-etaeqA
    #
    _inmark = i₀n*(exp.((1-an)*2*f*etan) - exp.(-an*2*f*etan))
    _ihmark = i₀h*(exp.((1-ah)*f*etah) - exp.(-ah*f*etah))
    # pure Butler-Volmer current densities
    # check if pure Butler-Volmer current densitie does not exceed limiting c.d. for
    # each electrodic reaction
    for i = 1:length(etac)
      if _inmark[i] < iln
        ɛ = 0.0000001
        # the added small quantity "ϵ" prevent the check from stopping the loop immediately
        _inmark[i] = iln + ɛ
      end
    end
    for i = 1:length(etac)
      if _ihmark[i] < ilh
        ɛ = 0.0000001
        _ihmark[i] = ilh + ɛ
      end
    end
    # check if the electrodic current density at the current etac value is cathodic
    # and if it is, evaluate the mass transport limitation overvoltage
    y = zeros(length(etac))
    for i = 1:length(etac)
      if _inmark[i] < 0
        etaln = 2/f*log(1-_inmark[i]/iln)
      else
        etaln = 0
      end
      if _ihmark[i] < 0
        etalh = 1/f*log(1-_ihmark[i]/ilh)
      else
        etalh = 0
      end
      # overall voltage balance (must be nought for the true etac)
      y[i] = etapd - etac[i] - etaln - etalh - res*area*(_inmark[i]+_ihmark[i])
    end
    # selects values corresponding to y(i)=0
    mark   = 0.0
    inmark = 0.0
    ihmark = 0.0
    for i = 1:length(y)
      if y[i] >= 0
        # mark: true etac
        mark = etac[i]
        # inmark: c.d. for metal (or first electrochemical reaction, in general) at true etac
        inmark = _inmark[i]
        # ihmark: c.d. for hydrogen (or second electrochemical reaction, in general) at true etac
        ihmark = _ihmark[i]
      end
    end
    amark[ibonji] = mark
    ainmark[ibonji] = inmark
    aihmark[ibonji] = ihmark
    aitmark[ibonji] = inmark + ihmark
  end
  # most external loop (=LSV loop) stops here
  # auxilary quantities for conventional semilog plot of potentiodynamic curve
  nteo = log.(abs.(ainmark/1000))/log(10)
  hteo = log.(abs.(aihmark/1000))/log(10)
  xteo = log.(abs.(aitmark/1000))/log(10)
  yteo = escan/1000
  xteo, yteo, nteo, hteo
end

# First model
linear_sweep_voltammogram()


# Second model
@time xteo, yteo, nteo, hteo = linear_sweep_voltammogram_lifelike()

# output a plot on the screen
plot( xteo,yteo,ls=:solid, lc=:black)
plot!(nteo,yteo,ls=:dot,   lc=:red)
plot!(hteo,yteo,ls=:dash,  lc=:blue,
    xlabel = "log (current density, j / (A cm^-2))",
    ylabel = "Overpotential, eta / V"
    )


