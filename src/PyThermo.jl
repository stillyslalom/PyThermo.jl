module PyThermo

# If the user has set JULIA_CONDAPKG_OFFLINE, use that value.  Otherwise,
# if .CondaPkg exists in the current environment, use "true" to skip
# Conda pkg resolution.  Otherwise, use "false" to allow Conda to
# resolve the pkg environment.
# using Pkg
# get!(ENV, "JULIA_CONDAPKG_OFFLINE") do
#     if ispath(joinpath(dirname(Pkg.project().path), ".CondaPkg"))
#         "yes"
#     else
#         "no"
#     end
# end

using CondaPkg
using PythonCall
using Printf
using Unitful

import Base: convert, ==, isequal, hash, getindex, setindex!, haskey, keys, show

export Thermo, ShockTube
export Species, Mixture
export temperature, pressure, isentropic_exponent, R_specific, soundspeed
export density, molar_density, molar_volume, compressibility
export heat_capacity, molar_heat_capacity
export enthalpy, molar_enthalpy, entropy, molar_entropy, internal_energy, molar_internal_energy
export viscosity, kinematic_viscosity, thermal_conductivity, thermal_diffusivity, prandtl, surface_tension
export isobaric_expansion, joule_thomson
export molecular_weight, CAS, formula, phase
export mole_fractions, mass_fractions, components
export T_critical, P_critical, acentric_factor, T_boiling, enthalpy_vaporization, P_saturation
export setstate!

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

const _DOCTEST = """jldoctest pythermo; setup = :(using PyThermo), filter = r"(\\d\\.\\d\\d)\\d*" => s"\\1" """

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
```$_DOCTEST
julia> He = Species("He")
Species(He, 298.1 K, 1.013e+05 Pa)

julia> density(He)
0.16360253235815483 kg m^-3

julia> using Unitful

julia> He.T = 30u"K"
30 K

julia> density(He)
1.623503007493497 kg m^-3
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
    o::Py
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
```$_DOCTEST
julia> air = Mixture(["N2" => 0.78, "O2" => 0.21, "Ar" => 0.01])
Mixture(78% N2, 21% O2, 1% Ar, 298.1 K, 1.013e+05 Pa)

julia> soundspeed(air)
346.1345044609487 m s^-1
```
"""
struct Mixture <: Chemical
    o::Py
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
    # s = ""
    # for (species, χ) in zip(species_str, mix.zs)
    #     s *= @sprintf("%s: %0.3g, ", species, χ)
    # end
    # s[1:end-2] * "}"
    join([@sprintf("%0.3g%s %s", 100*pyconvert(Float64, χ), '%', species) for (species, χ) in zip(species_str, mix.zs)], ", ")
end
composition_string(s::Species) = s.name

Base.show(io::IO, mix::Mixture) = @printf(io, "Mixture(%s, %0.1f K, %0.3e Pa)", composition_string(mix), mix.T, mix.P)

Py(c::Chemical) = getfield(c, :o)
convert(::Type{T}, o::PythonCall.Py) where {T <: Chemical} = T(o)
==(c1::Chemical, c2::Chemical) = pyconvert(Bool, Py(c1) == Py(c2))
Base.isequal(c1::Chemical, c2::Chemical) = isequal(Py(c1), Py(c2))
Base.hash(c::Chemical, h::UInt) = hash(Py(c), h)
Base.copy(c::T) where {T <: Chemical} = T(PY_COPY.copy(Py(c)))

Base.Docs.getdoc(c::Chemical, sig) = Base.Docs.getdoc(Py(c), sig)
Base.Docs.doc(c::Chemical, sig::Type=Union{}) = Base.Docs.getdoc(c, sig)
Base.Docs.Binding(c::Chemical, s::Symbol) = pygetattr(Py(c), String(s))

function Base.getproperty(c::Chemical, s::Symbol)
    if (s === :T) || (s === :P)
        pyconvert(Float64, getproperty(Py(c), s))::Float64
    elseif s === :ID
        pyconvert(String, getproperty(Py(c), s))::String
    elseif s === :IDs
        pyconvert(Vector{String}, getproperty(Py(c), s))::Vector{String}
    else
        pyconvert(Any, getproperty(Py(c), s))
    end
end
Base.setproperty!(c::Chemical, s::Symbol, x) = setproperty!(Py(c), s, x)
Base.hasproperty(c::Chemical, s::Symbol) = pyhasattr(Py(c), s)
Base.propertynames(c::Chemical) = propertynames(Py(c))
haskey(c::Chemical, x) = haskey(Py(c), x)

# Strip units for unitful `setproperty`
Base.setproperty!(c::Chemical, s::Symbol, T::Unitful.Temperature) = setproperty!(Py(c), s, _unit(T, u"K"))
Base.setproperty!(c::Chemical, s::Symbol, T::Unitful.Pressure) = setproperty!(Py(c), s, _unit(T, u"Pa"))

include("properties.jl")

# Thermodynamic property accessors. The bulk of these are macro-generated in
# `properties.jl`; the four below are kept hand-written because they don't
# need pyconvert plumbing (`T`/`P` already come back as `Float64`) or because
# their natural form is the bare dimensionless number returned by Python.

"""
    temperature(c::Chemical) -> Quantity{K}

Temperature at the chemical's current state.
"""
temperature(c::Chemical) = c.T * u"K"

"""
    pressure(c::Chemical) -> Quantity{Pa}

Pressure at the chemical's current state.
"""
pressure(c::Chemical)    = c.P * u"Pa"

"""
    isentropic_exponent(c::Chemical) -> Float64

Isentropic exponent γ = Cp / Cv at the chemical's current state.
Dimensionless.
"""
isentropic_exponent(c::Chemical) = c.isentropic_exponent

"""
    R_specific(c::Chemical) -> Quantity{J/(kg·K)}

Mass-basis specific gas constant R / MW.
"""
R_specific(c::Chemical)  = c.R_specific * u"J/(kg*K)"

"""
    compressibility(c::Chemical) -> Float64
    compressibility(c::Chemical, phase::Symbol) -> Float64

Compressibility factor `Z = PV/(nRT)` at the chemical's current state.
Dimensionless.

For a `Mixture`, forwards to thermo's EOS-aware `Z`.

For a `Species`, the no-argument call uses the attached cubic EOS
(Peng-Robinson by default) to return the real-gas gas-phase factor
(`eos.Z_g`), refreshing the EOS at the current T/P first. For non-gas phases,
or for species whose EOS does not expose `Z_g` (e.g. helium in the
supercritical-at-STP regime), it falls back to thermo's curve-based `Z`, which
is pinned near 1.0.

The explicit phase-argument form (`compressibility(c, :gas/:liquid/:solid)`)
returns thermo's **curve-based** per-phase value (`Zg`/`Zl`/`Zs`); this is not
EOS-derived and may sit near 1.0 for a `Species`. It therefore need not equal
the no-argument result, which is EOS-derived for a gas-phase `Species`.
"""
compressibility(c::Mixture) = pyconvert(Float64, Py(c).Z)

function compressibility(c::Species)
    if pyconvert(String, Py(c).phase) == "g"
        # Refresh the cached EOS at the current T/P before reading `Z_g` — see
        # `soundspeed` below for why direct `c.T = …` / `c.calculate(...)`
        # mutations do not rebuild the EOS object.
        Py(c).set_eos(T=Py(c).T, P=Py(c).P)
        eos = Py(c).eos
        pyhasattr(eos, "Z_g") && return pyconvert(Float64, eos.Z_g)
    end
    # Fallback: thermo's curve-based Z (pinned near 1.0 for a Species).
    pyconvert(Float64, Py(c).Z)
end

function compressibility(c::Chemical, phase::Symbol)
    attr = _phase_attr(:compressibility, phase, (gas=:Zg, liquid=:Zl, solid=:Zs))
    _wrap_strict(pygetattr(Py(c), String(attr)), NoUnits)
end

"""
    soundspeed(c::Chemical) -> Quantity{m/s}

Real-gas sound speed at the chemical's current state.

For a `Mixture`, forwards to thermo's `speed_of_sound` directly.

For a `Species`, uses the attached cubic EOS (Peng-Robinson by default) to
evaluate

    a² = γ_real · (∂P/∂ρ_mass)_T

with `γ_real = (Cpgm + Cp_dep_g) / (Cvgm + Cv_dep_g)` and
`(∂P/∂ρ_mass)_T = eos.dP_drho_g / MW`. The EOS-derived path is gas-phase
only — for non-gas phases or for species whose EOS does not expose
departure derivatives (e.g. helium in the supercritical-at-STP regime), the
ideal-gas formula `a = √(γ R T)` is used instead.
"""
soundspeed(c::Mixture) = pyconvert(Float64, Py(c).speed_of_sound) * u"m/s"

function soundspeed(c::Species)
    phase = pyconvert(String, Py(c).phase)
    if phase == "g"
        # `Chemical.set_eos` rebuilds the attached EOS at the chemical's
        # current T/P. This is necessary because direct `c.T = …` /
        # `c.calculate(...)` mutations do not refresh the cached EOS object,
        # so the departure derivatives read below would otherwise be stale.
        Py(c).set_eos(T=Py(c).T, P=Py(c).P)
        eos = Py(c).eos
        if pyhasattr(eos, "Cp_dep_g") && pyhasattr(eos, "Cv_dep_g") && pyhasattr(eos, "dP_drho_g")
            Cp_real = pyconvert(Float64, Py(c).Cpgm) + pyconvert(Float64, eos.Cp_dep_g)
            Cv_real = pyconvert(Float64, Py(c).Cvgm) + pyconvert(Float64, eos.Cv_dep_g)
            dP_drho_molar = pyconvert(Float64, eos.dP_drho_g)
            MW_kg = pyconvert(Float64, Py(c).MW) * 1e-3
            return sqrt((Cp_real / Cv_real) * dP_drho_molar / MW_kg) * u"m/s"
        end
    end
    # Fallback: ideal-gas formula. Used for non-gas phases and for species
    # whose attached EOS does not expose departure derivatives (e.g. helium
    # in the supercritical-at-STP regime).
    sqrt(isentropic_exponent(c) * R_specific(c) * temperature(c)) |> u"m/s"
end

"""
    setstate!(c::Chemical; T=nothing, P=nothing) -> Chemical

Set temperature and/or pressure and re-flash the chemical's state. Unlike
assignment to `c.T` / `c.P`, this calls the underlying `calculate` method,
which correctly handles phase changes.

Either or both of `T` and `P` may be provided. `Unitful` quantities are
converted to K / Pa; bare `Real` values are taken to already be in K / Pa.

# Examples
```$_DOCTEST
julia> using Unitful

julia> SF6 = Species("SF6")
Species(SF6, 298.1 K, 1.013e+05 Pa)

julia> setstate!(SF6, T=20u"K");

julia> SF6.phase
"s"
```
"""
function setstate!(c::Chemical; T=nothing, P=nothing)
    kwargs = Dict{Symbol, Float64}()
    T !== nothing && (kwargs[:T] = T isa Unitful.Temperature ? _unit(T, u"K") : Float64(T))
    P !== nothing && (kwargs[:P] = P isa Unitful.Pressure    ? _unit(P, u"Pa") : Float64(P))
    Py(c).calculate(; kwargs...)
    return c
end

include("ShockTube.jl")

# if ccall(:jl_generating_output, Cint, ()) == 1
    # __init__()
    # ShockTube.shockcalc(Species("N2"), Mixture(["N2" => 0.78, "O2" => 0.21, "Ar" => 0.01]), 1.8)
# end
end # module
