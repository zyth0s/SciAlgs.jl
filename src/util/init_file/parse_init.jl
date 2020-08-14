
import Pkg.TOML

# ----------------------------------------
# From a string
input = """
[section1]
var = 1
foo = [1,2,3]

[section2]
alg = "best"

"""

config = TOML.parse(input)
@assert config["section1"]["var"] == 1

# ----------------------------------------
# From a file
config = TOML.parsefile("init.toml")

@info "Reproduce the input file used for the calculation:" 
TOML.print(config)

# Insecure
etol = 0.0
etol = config["section_hf"]["etol"]

# Secure
if haskey(config, "section_hf")

   section_hf = config["section_hf"]

   dtol = get(section_hf, "dtol", 1e-8) # with default value
   symm = get(section_hf, "symm", "C1") # with default value
end

@info "Parsed parameters:" dtol  etol symm
