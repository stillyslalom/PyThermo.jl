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
1.6237570273514512 kg m^-3

julia> argon.molecular_diameter # Angstrom
3.40744
```
If left unspecified, the pressure and temperature are set to STP. if given as `Float64` numbers, `P` and `T` must have units of Pascal and Kelvin, respectively. The state variables can also be set as `Unitful` quantities.
```julia
julia> using Unitful

julia> air = Mixture(["N2" => 0.78, "O2" => 0.21, "Ar" => 0.01], P = 1u"atm")
Mixture(78% nitrogen, 21% oxygen, 1% argon, 298.1 K, 1.013e+05 Pa)

julia> air.Cp
1004.1326426200408
```

State-dependent properties are updated for a given phase when the `T` or `P` fields are set, but in case of phase change, the state variables must be updated via the bound `calculate` method.
```julia
julia> SF6 = Species("SF6", P=30u"psi")
Species(SF6, 298.1 K, 2.068e+05 Pa)

julia> SF6.GWP # global warming potential
22800.0

julia> SF6.phase # one of "g", "l", "s"
"g"

julia> SF6.T = 20u"K"
20 K

julia> density(SF6)
181.67446044245906 kg m^-3

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

### Future development

- [ ] Add `Unitful` accessors for more properties
- [ ] Use `missing` instead of `nothing` for missing properties

### Easter eggs
PyThermo includes a sub-module for 1D gas dynamics capable of calculating shock properties for gas species and mixtures.
```julia
julia> using PyThermo, PyThermo.ShockTube, Unitful

julia> driver = Species("He")
Species(He, 298.1 K, 1.013e+05 Pa)

julia> driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
Mixture(95% helium, 5% acetone, 291.1 K, 8.500e+04 Pa)

julia> shockjump(driven, 2.2) # find shock jump conditions
(Mixture(95% helium, 5% acetone, 623.7 K, 4.819e+05 Pa), 1023.9438673559401 m s^-1)

julia> shockcalc(driver, driven, 2.2) # calculate shock states for M=2.2
  Region    P [MPa] T [K] ρ [kg/m³] cₛ [m/s]
  ––––––––– ––––––– ––––– ––––––––– ––––––––
  Driver     3.732  298.1   5.92      1016
  Driven     0.085  291.1  0.2355    748.1
  Shocked   0.4819  623.7  0.6231     1066  
  Reflected  1.754  1003    1.408     1332

Driver gas: helium
Driven gas: 95% helium, 5% acetone
Post-shock velocity: 1024 m/s

julia> density(ans.shocked) / density(ans.driven)
2.64075735975203
```