
#@static if VERSION >= v"1.5.0" # nicer syntax

   #const FORTLIB = "fortranlib.so"
   #five = @ccall FORTLIB.__m_MOD_five()::Int # ??
   #println("Five = ",five)
#else
   println("Using function")
   five = ccall( (:__m_MOD_five, "fortranlib.so"), Int, () )
   println("Five = ",five)
#end
