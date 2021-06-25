module PyThermo

using PyCall 
using Printf
using Unitful

import Base: convert, ==, isequal, hash, getindex, setindex!, haskey, keys, show

export Thermo, ShockTube
export Species, Mixture
export isentropic_exponent, temperature, pressure, density, molar_density, R_specific, soundspeed

include("init.jl")
const Thermo = PY_THERMO

_unit(x, u) = Float64(ustrip(uconvert(u, x)))

function _SI_TP(kwargs)
    d = Dict{Any, Any}(kwargs)
    T = get!(d, :T, 298.15)
    P = get!(d, :P, 101325.0)
    d[:T] = (T isa Unitful.Temperature) ? _unit(T, u"K") : Float64(T)
    d[:P] = (P isa Unitful.Pressure)    ? _unit(P, u"Pa") : Float64(P)
    return (d...,)
end

abstract type Chemical end

"""
    Species(ID::String; T=298.15, P=101325)

Creates a `Species` object which contains basic information such as
molecular weight and the structure of the species, as well as thermodynamic
and transport properties as a function of temperature and pressure.
If `T` and `P` are given as non-`Unitful` numbers, they must have units of `K` and `Pa`.

Parameters
----------
ID : One of the following [-]:

* Name, in IUPAC form or common form or a synonym registered in PubChem \n
* InChI name, prefixed by "InChI=1S/" or "InChI=1/" \n
* InChI key, prefixed by "InChIKey=" \n
* PubChem CID, prefixed by "PubChem=" \n
* SMILES (prefix with "SMILES=" to ensure smiles parsing) \n
* CAS number
    
T : temperature of the chemical (default 298.15 K) \n
P : pressure of the chemical (default 101325 Pa)

Examples
--------
```jldoctest
julia> He = Species("He")
Species(He, 298.1 K, 1.013e+05 Pa)

julia> density(He)
0.16352343013746262 kg m^-3

julia> using Unitful

julia> He.T = 30u"K"
30 K

julia> density(He)
1.6235030074934973 kg m^-3
```

A wide variety of unexported properties can be accessed from the underlying Python object:
```
julia> SF6 = Species("SF6", P=30u"psi")
Species(SF6, 298.1 K, 2.068e+05 Pa)

julia> SF6.MW
146.055419

help?> SF6.<tab>
A                                 SGs                                __init_subclass__
API                               STEL                               __le__
Am                                STEL_source                        __lt__
Bond                              STEL_sources                       __module__
Bvirial                           S_dep_Tb_P_ref_g                   __ne__
CAS                               S_dep_Tb_Pb_g                      __new__
Capillary                         S_dep_Tb_Pb_l                      __reduce__
Carcinogen                        S_dep_ref_g                        __reduce_ex__
Carcinogen_source                 S_int_Tb_to_T_ref_g                __repr__
Carcinogen_sources                S_int_l_Tm_to_Tb                   __setattr__
[...]
```
"""
struct Species <: Chemical
    o::PyObject
end
Species(chemname::String; kwargs...) = Species(PY_CHEM.Chemical(chemname; _SI_TP(kwargs)...))

Base.show(io::IO, species::Species) = @printf(io, "Species(%s, %0.1f K, %0.3e Pa)", species.ID, species.T, species.P)


"""
    Mixture(chemnames::Vector{String}; kwargs...)
    Mixture(chemnames::Vector{Pair{String, Float64}}, kwargs...)

Creates a `Mixture` object which contains basic information such as
molecular weight and the structure of the species, as well as thermodynamic
and transport properties as a function of temperature and pressure.

The components of the mixture are specified by the names of
the chemicals; the composition can be specified by providing any one of the
following parameters as a keyword argument:

* Mass fractions `ws`
* Mole fractions `zs`
* Liquid volume fractions (based on pure component densities) `Vfls`
* Gas volume fractions (based on pure component densities) `Vfgs`

The composition can also be specified by providing a vector of `"ID" => molefrac` pairs.

Examples
--------
```jldoctest
julia> air = Mixture(["N2" => 0.78, "O2" => 0.21, "Ar" => 0.01])
Mixture({N2: 0.78, O2: 0.21, Ar: 0.01}, 298.1 K, 1.013e+05 Pa)

julia> soundspeed(air)
346.1466532754559 m s^-1
```
"""
struct Mixture <: Chemical
    o::PyObject
end
Mixture(chemnames::Vector{String}; kwargs...) = Mixture(PY_CHEM.Mixture(chemnames; _SI_TP(kwargs)...))

function Mixture(chems::Vector{Pair{String, Float64}}; kwargs...)
    Mixture(first.(chems); kwargs..., zs = last.(chems))
end

function Mixture(chems::Vector{Pair{T1, T2}}; kwargs...) where {T1, T2}
    Mixture(string.(first.(chems)) .=> Float64.(last.(chems)); kwargs...)
end

# Mixture(args...; kwargs...) = Mixture([args...], kwargs...)

function composition_string(mix)
    species_str = try
        mix.IDs
    catch
        mix.names
    end
    s = "{"
    for (species, χ) in zip(species_str, mix.zs)
        s *= @sprintf("%s: %0.3g, ", species, χ)
    end
    s[1:end-2] * "}"
end
Base.show(io::IO, mix::Mixture) = @printf(io, "Mixture(%s, %0.1f K, %0.3e Pa)", composition_string(mix), mix.T, mix.P)

PyObject(c::Chemical) = getfield(c, :o)
convert(::Type{T}, o::PyCall.PyObject) where {T <: Chemical} = T(o)
==(c1::Chemical, c2::Chemical) = PyObject(c1) == PyObject(c2)
hash(c::Chemical) = hash(PyObject(c))
Base.copy(c::T) where {T <: Chemical} = T(PY_COPY.copy(PyObject(c)))

Base.Docs.doc(c::Chemical) = Base.Docs.doc(PyObject(c))
Base.Docs.Binding(c::Chemical, s::Symbol) = Base.Docs.Binding(PyObject(c), s)

Base.getproperty(c::Chemical, s::Symbol) = getproperty(PyObject(c), s)
Base.setproperty!(c::Chemical, s::Symbol, x) = setproperty!(PyObject(c), s, x)
Base.hasproperty(c::Chemical, s::Symbol) = hasproperty(PyObject(c), s)
Base.propertynames(c::Chemical) = propertynames(PyObject(c))
haskey(c::Chemical, x) = haskey(PyObject(c), x)

# Strip units for unitful `setproperty`
Base.setproperty!(c::Chemical, s::Symbol, T::Unitful.Temperature) = setproperty!(PyObject(c), s, _unit(T, u"K"))
Base.setproperty!(c::Chemical, s::Symbol, T::Unitful.Pressure) = setproperty!(PyObject(c), s, _unit(T, u"Pa"))


# Thermodynamic property accessors
temperature(c::Chemical)  = c.T * u"K"
pressure(c::Chemical)  = c.P * u"Pa"
density(c::Chemical)  = c.rho * u"kg/m^3"
molar_density(c::Chemical) = c.rhom * u"mol/m^3"
isentropic_exponent(c::Chemical)  = c.isentropic_exponent
R_specific(c::Chemical)  = c.R_specific * u"J/(kg*K)"
soundspeed(c::Chemical) = sqrt(isentropic_exponent(c) * R_specific(c) * temperature(c)) |> u"m/s"  

include("ShockTube.jl")
end
