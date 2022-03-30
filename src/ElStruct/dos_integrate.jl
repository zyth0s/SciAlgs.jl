
# Integrate a density of states (DOS) in some energy
# range from a set of data points of the DOS

using REPL.TerminalMenus

function read_dos_file()
   # Enter the filename with the DOS data
   println("* Please, enter the DOS file name:")
   filename = readline(stdin)
   isfile(filename) || error("$filename is not a valid filename")

   # Specify the number of spin channels
   options = ["Single spin channel (α=β)",
              "Double spin channel (α,β)"]
   menu_nr_spin_channels = RadioMenu(options)
   nr_spin_channels = request("* Select the number of spin channels:",
                              menu_nr_spin_channels)

   # Energy and DOS ticks
   energy_array = Float64[]
   dos_α_array  = Float64[]
   dos_β_array  = Float64[]

   open(filename) do file

      # Read line by line until the EOF
      while !eof(file)

         line = readline(file)
         # Skip comments
         startswith(line, "#") && continue
         # Skip empty lines (including whitespaces)
         isempty(strip(line)) && continue

         energy = split(line)[1]
         push!(energy_array, parse(Float64, energy))
         dos_α = split(line)[2]
         push!(dos_α_array, parse(Float64, dos_α))
         dos_β = Missing
         if nr_spin_channels == 2
            dos_β = split(line)[3]
            push!(dos_β_array, parse(Float64, dos_β))
         end
      end

   end
   energy_array, dos_α_array, dos_β_array
end

# Trapezoidal integration rule (should be better than Simpson for rought functions ← DOS is not smooth)
function trapezoidal(x_array, y_array, x_min, x_max)
   ∑ = 0.0
   for i in 1:length(x_array)-1
      if x_min ≤ x_array[i] ≤ x_max
         P = 0.5(y_array[i] + y_array[i+1])
         ∑ += P*(x_array[i+1] - x_array[i]) # = ∑ₐᵇ ½(P(xᵢ₊₁) + P(xᵢ)) Δxᵢ
      end
   end
   ∑
end


# Run example

# → ../data/LiAlSi.dos
energy_array, dos_α_array, dos_β_array = read_dos_file()

# Enter the energy window of integrations
# → -11 -7
println("* Please, enter the energy window [e_min, e_max]: (e.g. e_min e_max)")
e_window_str = split(readline(stdin))
e_min = parse(Float64, e_window_str[1])
e_max = parse(Float64, e_window_str[2])

∑dos_α = trapezoidal(energy_array, dos_α_array, e_min, e_max)
println("→ There are $∑dos_α states (electrons) in the [$e_min, $e_max] energy_range")
# → 2 elec
