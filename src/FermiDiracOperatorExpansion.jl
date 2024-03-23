module FermiDiracOperatorExpansion

abstract type Solver end
struct CG <: Solver end
struct NewtonSchulz <: Solver end

end
