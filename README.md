# PyThermo

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stillyslalom.github.io/PyThermo.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stillyslalom.github.io/PyThermo.jl/dev)
[![Build Status](https://github.com/stillyslalom/PyThermo.jl/workflows/CI/badge.svg)](https://github.com/stillyslalom/PyThermo.jl/actions)

PyThermo.jl offers access to [Thermo](https://github.com/CalebBell/thermo), a Python library for thermophysical properties with data available for a large number of chemical species. Accessor functions are available for basic properties (density, pressure, gas constant, ...), with returned values in SI units. Many additional properties can be accessed as properties of the underlying `PyObject`, although these internal properties are non-`Unitful`. Those properties are detailed for [`Species`](https://thermo.readthedocs.io/thermo.chemical.html) and [`Mixture`](https://thermo.readthedocs.io/thermo.mixture.html) in Thermo's external documentation.
```julia
julia> using PyThermo

julia> argon = Species("Ar", P = 101325, T = 300)
Species(Ar, 300.0 K, 1.013e+05 Pa)

julia> density(argon)
1.6227671732556135 kg m^-3

julia> argon.molecular_diameter # Angstrom
3.40744
```
If left unspecified, the pressure and temperature are set to STP. if given as `Float64` numbers, `P` and `T` must have units of Pascal and Kelvin, respectively. The state variables can also be set as `Unitful` quantities.
```julia
julia> using Unitful

julia> air = Mixture(["N2" => 0.78, "O2" => 0.21, "Ar" => 0.01], P = 1u"atm")
Mixture(78% N2, 21% O2, 1% Ar, 298.1 K, 1.013e+05 Pa)

julia> air.Cp
1004.3091399287284
```

### Combining mixtures and species

A `Mixture` can be assembled from existing `Species`, other `Mixture`s, or
chemical-name strings, given as `constituent => amount` pairs. Amounts are
relative mole numbers by default (`basis=:mole`; pass `basis=:mass` for
relative masses), they need not sum to one, and a species appearing in more
than one constituent is merged by CAS number.
```julia
julia> air = Mixture(["N2" => 0.79, "O2" => 0.21]);

julia> Mixture([air => 0.8, Species("acetone") => 0.2])
Mixture(63.2% N2, 16.8% O2, 20% acetone, 298.1 K, 1.013e+05 Pa)

julia> Mixture([air => 0.8, Species("acetone") => 0.2]; basis=:mass)
Mixture(70.3% N2, 18.7% O2, 11% acetone, 298.1 K, 1.013e+05 Pa)
```
By default the result is built at STP (override with `T=`/`P=`), using only the
constituents' composition. Passing `adiabatic=true` instead conserves enthalpy
at constant pressure, solving for the equilibrium temperature and accounting
for heat capacity and latent heat — so mixing room-temperature *liquid* acetone
into air evaporatively cools the result and leaves a two-phase state:
```julia
julia> mix = Mixture([air => 0.85, Species("acetone") => 0.15]; adiabatic=true)
Mixture(67.2% N2, 17.8% O2, 15% acetone, 263.0 K, 1.013e+05 Pa)

julia> phase(mix)
:two_phase
```

State-dependent properties are updated for a given phase when the `T` or `P` fields are set, but in case of phase change, the state variables must be updated via the bound `calculate` method.
```julia
julia> SF6 = Species("SF6", P=30u"psi")
Species(SF6, 298.1 K, 2.068e+05 Pa)

julia> SF6.GWP # global warming potential
16300.0

julia> SF6.phase # one of "g", "l", "s"
"g"

julia> SF6.T = 20u"K"
20 K

julia> density(SF6)
181.6743991134152 kg m^-3

julia> SF6.phase # that's not right!
"g"

julia> SF6.calculate(T = 20)

julia> SF6.phase # much better!
"s"
```

### Installation
PyThermo.jl is registered in Julia's general package repository and can be installed with the package manager. [CondaPkg.jl](https://github.com/cjdoris/CondaPkg.jl) is used to automatically install all Python dependencies to a PyThermo-specific Conda environment unless otherwise specified (see CondaPkg docs for details).
```
(v1.8) pkg> add PyThermo
   Resolving package versions...
   Installed PyThermo ─ v0.2.0
```

### Curated accessors

A curated set of Unitful Julia accessors is now available for the bulk of
the commonly-used `thermo` properties — see the
[Property accessors](https://stillyslalom.github.io/PyThermo.jl/stable/properties/)
page of the docs for the full list and a `thermo`-name → PyThermo-name
cheatsheet. State variables and state-dependent properties (T, P, ρ, Cp,
μ, k, …) are strict; per-substance constants with a real likelihood of
being absent from the database (Tc, Pc, ω, Tb, Hvap, Psat) return
`missing` rather than throwing.

`soundspeed` evaluates a real-gas formula via the chemical's cubic EOS
(Peng-Robinson by default) for `Species`, and forwards to thermo's
`speed_of_sound` for `Mixture`. The previous ideal-gas √(γRT) expression is
recovered for near-ideal species but yields noticeably different values for
dense or strongly non-ideal gases (e.g. SF6 above ~1 atm).

A `setstate!(c; T, P)` helper avoids the phase-change footgun shown above
by routing temperature/pressure updates through thermo's `calculate` method.

### Easter eggs
PyThermo includes a sub-module for 1D gas dynamics capable of calculating shock properties for gas species and mixtures.
```julia
julia> using PyThermo, PyThermo.ShockTube, Unitful

julia> driver = Species("He")
Species(He, 298.1 K, 1.013e+05 Pa)

julia> driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
Mixture(95% He, 5% acetone, 291.1 K, 8.500e+04 Pa)

julia> shockjump(driven, 2.2) # find shock jump conditions
(Mixture(95% He, 5% acetone, 623.6 K, 4.818e+05 Pa), 1023.9538914888177 m s^-1)

julia> shockcalc(driver, driven, 2.2) # calculate shock states for M=2.2
  Region    P [MPa] T [K] ρ [kg/m³] cₛ [m/s]
  ––––––––– ––––––– ––––– ––––––––– ––––––––
  Driver     3.732  298.1   6.025     1016
  Driven     0.085  291.1  0.2355    748.1
  Shocked   0.4818  623.6  0.6232     1095
  Reflected  1.755  1003    1.412     1388

  Incident wave: 1646 m/s (Mach 2.2)
  Reflected wave: 1965 m/s (Mach 1.795)
  Post-shock velocity: 1024 m/s

julia> density(ans.shocked) / density(ans.driven)
2.6465904286501147
```