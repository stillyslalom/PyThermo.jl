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
Mixture({N2: 0.78, O2: 0.21, Ar: 0.01}, 298.1 K, 1.013e+05 Pa)

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
181.67446044245906 kg m^-

julia> SF6.phase # that's not right!
"g"

julia> SF6.calculate(T = 20)

julia> SF6.phase # much better!
"s"
```

### Future development

- [ ] Add `Unitful` accessors for more properties
- [ ] Use `missing` instead of `nothing` for missing properties