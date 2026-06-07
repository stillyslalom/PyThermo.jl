```@meta
CurrentModule = PyThermo
```
# Mixtures

A [`Mixture`](@ref) holds a multi-component fluid and exposes the same
thermodynamic and transport accessors as a [`Species`](@ref). There are two
ways to build one:

1. Flat construction, where you name the leaf components and give their
   composition directly.
2. Combining, where you assemble a mixture from existing `Species`, other
   `Mixture`s, or chemical names, given as `constituent => amount` pairs.

The state variables `T` and `P` follow the same conventions as `Species`.
They default to STP and accept either bare `K`/`Pa` numbers or `Unitful`
quantities, unless an adiabatic mix solves for the temperature instead.

## Flat construction

List the component names and supply the composition through one of thermo's
composition keywords: mole fractions `zs`, mass fractions `ws`, or liquid and
gas volume fractions `Vfls` and `Vfgs`.

```julia
using PyThermo, Unitful

Mixture(["N2", "O2", "Ar"]; zs = [0.78, 0.21, 0.01])
Mixture(["water", "ethanol"]; ws = [0.6, 0.4], T = 60u"°C")
```

As a shorthand for the mole-fraction case, pass a vector of
`"name" => mole_fraction` pairs:

```julia
julia> air = Mixture(["N2" => 0.79, "O2" => 0.21])
Mixture(79% N2, 21% O2, 298.1 K, 1.013e+05 Pa)
```

Fractions need not sum to one, since they are normalized.

## Combining species and mixtures

A `Mixture` can also be assembled from a list of `constituent => amount`
pairs, where each constituent is a [`Species`](@ref), another `Mixture`, or a
chemical-name `String`, and `amount` is a relative quantity (the amounts need
not sum to one). Every constituent is flattened to its components, and like
species are merged by CAS number, so a component appearing in more than one
constituent is combined into a single entry rather than duplicated:

```julia
julia> air = Mixture(["N2" => 0.79, "O2" => 0.21]);

julia> Mixture([air => 0.8, Species("acetone") => 0.2])
Mixture(63.2% N2, 16.8% O2, 20% acetone, 298.1 K, 1.013e+05 Pa)
```

```julia
julia> # air already contains N2; adding more merges instead of duplicating

julia> Mixture([air => 0.5, Species("N2") => 0.5])
Mixture(89.5% N2, 10.5% O2, 298.1 K, 1.013e+05 Pa)
```

A bare string is treated as `Species(name)`, so constituents can be mixed and
matched freely, and a string-only list reproduces the flat
`"name" => mole_fraction` behaviour:

```julia
Mixture([air => 0.8, "acetone" => 0.2])   # string is Species("acetone")
```

By default the combined mixture is built at STP (override with the usual `T=`
and `P=` keywords). Only the constituents' composition is used; their
individual states are ignored unless you opt into [adiabatic mixing](@ref
"Adiabatic mixing").

## Mole vs mass basis

The `basis` keyword sets how each constituent's `amount` is interpreted.
`:mole`, the default, treats amounts as relative moles, while `:mass` treats
them as relative masses. The composition is always displayed as mole
fractions, so the two bases give different results for the same numbers:

```julia
julia> Mixture([air => 0.8, Species("acetone") => 0.2])
Mixture(63.2% N2, 16.8% O2, 20% acetone, 298.1 K, 1.013e+05 Pa)

julia> Mixture([air => 0.8, Species("acetone") => 0.2]; basis = :mass)
Mixture(70.3% N2, 18.7% O2, 11% acetone, 298.1 K, 1.013e+05 Pa)
```

This is distinct from the flat constructor's `zs` and `ws` keywords. Those
give the final leaf-component fractions directly, whereas `basis` governs how
the per-constituent amounts are split, and each constituent may expand into
several components. Use whichever matches how the quantity is known: `:mole`
for a fluid metered by molar flow, `:mass` for one weighed out.

## Adiabatic mixing

Passing `adiabatic = true` performs a constant-pressure, enthalpy-conserving
mix. Instead of placing the result at a chosen temperature, the constituents'
temperatures, phases, heat capacities and latent heats are accounted for and
the equilibrium temperature is solved for, so `T` may not be given.

The headline case is evaporative cooling. Mixing room-temperature liquid
acetone into air lets the acetone evaporate, drawing its heat of vaporization
from the gas and dropping the temperature well below 298 K, leaving a
two-phase state:

```julia
julia> air = Mixture(["N2" => 0.79, "O2" => 0.21]);

julia> mix = Mixture([air => 0.85, Species("acetone") => 0.15]; adiabatic = true)
Mixture(67.2% N2, 17.8% O2, 15% acetone, 263.0 K, 1.013e+05 Pa)

julia> phase(mix)
:two_phase
```

When no phase change occurs the result is a plain sensible-heat balance. For
example, equal moles of argon at 300 K and nitrogen at 500 K equilibrate near
417 K, weighted toward the higher heat capacity of the diatomic:

```julia
julia> Ar = Species("Ar", T = 300u"K"); N2 = Species("N2", T = 500u"K");

julia> Mixture([Ar => 0.5, N2 => 0.5]; adiabatic = true)
Mixture(50% Ar, 50% N2, 417.2 K, 1.013e+05 Pa)
```

The `basis` keyword applies here too. `:mass` conserves mass-specific enthalpy
through thermo's `H` flash rather than molar `Hm`. Because both bases describe
the same physical system, equivalent amounts reach the same equilibrium
temperature.

### Pressure and convergence

If `P` is omitted, it is inherited from the constituents when they all agree,
and otherwise defaults to 1 atm; pass `P=` explicitly to pin it. The
equilibrium is found by thermo's enthalpy flash, which lacks usable data for
some systems, notably strongly supercritical species such as helium. When the
flash fails to converge the constructor raises an informative `ArgumentError`
rather than returning a malformed mixture. Build the mixture at an explicit
`T` in that case.

## Reference

```@docs
Mixture
```
