```@meta
CurrentModule = PyThermo
```

# Property accessors

PyThermo provides curated Julia accessor functions for the most common
thermodynamic, transport, and identity properties of `Species` and `Mixture`.
Each accessor wraps the underlying `thermo` Python attribute, attaches a
`Unitful` unit, and (for constants that may be absent from the database)
returns `missing` rather than `nothing`.

For properties not covered by a curated accessor, the `getproperty`
fallthrough still works: `c.molecular_diameter`, `c.GWP`, `c.Hf`, etc. See
the upstream
[`thermo.Chemical`](https://thermo.readthedocs.io/thermo.chemical.html) and
[`thermo.Mixture`](https://thermo.readthedocs.io/thermo.mixture.html)
references for the full attribute list.

## Two-tier missing-value handling

State-dependent properties (T, P, density, Cp, viscosity, …) are **strict**:
their return type is `Quantity`, and accessing a property that thermo can
not evaluate at the current state throws.

Per-substance constants whose database entry may genuinely be absent (`Tc`,
`Pc`, `omega`, `Tb`, `Hvap`, `Psat`) are **optional**: they return `missing`
when thermo returns `None`, and otherwise a `Quantity`. The one dimensionless
optional, `acentric_factor`, follows the same `NoUnits` convention as the
strict dimensionless accessors (`prandtl`, the curve-based path of
`compressibility`) and so returns `Union{Float64, Missing}` rather than
`Union{Quantity, Missing}`.

## Phase variants

State-dependent properties whose Python attribute has phase-suffixed
variants (e.g. `rhog`, `rhol`, `rhos` for gas/liquid/solid density) accept
an optional positional `Symbol`:

```julia
density(c)            # current phase, equivalent to c.rho
density(c, :gas)      # c.rhog
density(c, :liquid)   # c.rhol
density(c, :solid)    # c.rhos   (where defined)
```

Phases that thermo does not expose for a given property throw
`ArgumentError` — for example `viscosity(c, :solid)` (no `mus` attribute).
The supported phases for each accessor are listed in the cheatsheet below
and in the per-accessor docstring.

## State

```@docs
temperature
pressure
density
molar_density
molar_volume
compressibility
phase
```

`phase` returns `:gas`, `:liquid`, `:solid`, or `:two_phase`. The last covers a
coexisting-phase state, which a `Mixture` reaches when a flash lands inside a
phase envelope (see [adiabatic mixing](@ref "Adiabatic mixing")).

## Caloric

```@docs
heat_capacity
molar_heat_capacity
isochoric_heat_capacity
molar_isochoric_heat_capacity
enthalpy
molar_enthalpy
entropy
molar_entropy
internal_energy
molar_internal_energy
```

## Transport

```@docs
viscosity
kinematic_viscosity
thermal_conductivity
thermal_diffusivity
prandtl
surface_tension
```

## Derived

```@docs
isentropic_exponent
R_specific
soundspeed
isobaric_expansion
joule_thomson
```

## Identity

```@docs
molecular_weight
CAS
formula
```

## Mixture composition

```@docs
mole_fractions
mass_fractions
components
```

## Optional constants

```@docs
T_critical
P_critical
acentric_factor
T_boiling
enthalpy_vaporization
P_saturation
```

## State updates

```@docs
setstate!
```

## thermo → PyThermo cheatsheet

The table below maps `thermo`'s Python attribute names to the corresponding
PyThermo accessors. Phase-variant columns indicate the supported `phase`
argument for the no-suffix accessor.

### State

| thermo attribute(s)                       | PyThermo accessor                  | Unit       |
|-------------------------------------------|------------------------------------|------------|
| `T`                                       | `temperature(c)`                   | K          |
| `P`                                       | `pressure(c)`                      | Pa         |
| `rho` / `rhog` / `rhol` / `rhos`          | `density(c[, :gas/:liquid/:solid])`| kg/m³      |
| `rhom` / `rhogm` / `rholm` / `rhosm`      | `molar_density(c[, phase])`        | mol/m³     |
| `Vm` / `Vmg` / `Vml` / `Vms`              | `molar_volume(c[, phase])`         | m³/mol     |
| `Z` / `Zg` / `Zl` / `Zs` (see [`compressibility`](@ref)) | `compressibility(c[, phase])` | —          |
| `phase`                                   | `phase(c)` *(`:gas`/`:liquid`/`:solid`/`:two_phase`)* | — |

### Caloric

| thermo attribute(s)                       | PyThermo accessor                  | Unit       |
|-------------------------------------------|------------------------------------|------------|
| `Cp` / `Cpg` / `Cpl` / `Cps`              | `heat_capacity(c[, phase])`        | J/(kg·K)   |
| `Cpm` / `Cpgm` / `Cplm` / `Cpsm`          | `molar_heat_capacity(c[, phase])`  | J/(mol·K)  |
| `Cvg` *(gas only)*                        | `isochoric_heat_capacity(c)`       | J/(kg·K)   |
| `Cvgm` *(gas only)*                       | `molar_isochoric_heat_capacity(c)` | J/(mol·K)  |
| `H`                                       | `enthalpy(c)`                      | J/kg       |
| `Hm`                                      | `molar_enthalpy(c)`                | J/mol      |
| `S`                                       | `entropy(c)`                       | J/(kg·K)   |
| `Sm`                                      | `molar_entropy(c)`                 | J/(mol·K)  |
| `U`                                       | `internal_energy(c)`               | J/kg       |
| `Um`                                      | `molar_internal_energy(c)`         | J/mol      |

### Transport

| thermo attribute(s)                       | PyThermo accessor                  | Unit       |
|-------------------------------------------|------------------------------------|------------|
| `mu` / `mug` / `mul`                      | `viscosity(c[, :gas/:liquid])`     | Pa·s       |
| `nu` / `nug` / `nul`                      | `kinematic_viscosity(c[, phase])`  | m²/s       |
| `k` / `kg` / `kl`                         | `thermal_conductivity(c[, phase])` | W/(m·K)    |
| `alpha` / `alphag` / `alphal`             | `thermal_diffusivity(c[, phase])`  | m²/s       |
| `Pr` / `Prg` / `Prl`                      | `prandtl(c[, phase])`              | —          |
| `sigma`                                   | `surface_tension(c)`               | N/m        |

### Derived

| thermo attribute(s)                                | PyThermo accessor              | Unit       |
|----------------------------------------------------|--------------------------------|------------|
| `isentropic_exponent`                              | `isentropic_exponent(c)`       | —          |
| `R_specific`                                       | `R_specific(c)`                | J/(kg·K)   |
| EOS-derived (see [`soundspeed`](@ref))             | `soundspeed(c)`                | m/s        |
| `isobaric_expansion` / `_g` / `_l`                 | `isobaric_expansion(c[, phase])` | 1/K      |
| `JT` / `JTg` / `JTl`                               | `joule_thomson(c[, phase])`    | K/Pa       |

### Identity & composition

| thermo attribute(s)         | PyThermo accessor              | Return type                       |
|-----------------------------|--------------------------------|-----------------------------------|
| `MW`                        | `molecular_weight(c)`          | `Quantity{g/mol}`                 |
| `CAS`                       | `CAS(s)` *(Species only)*      | `String`                          |
| `formula`                   | `formula(s)` *(Species only)*  | `String`                          |
| `zs`                        | `mole_fractions(m)`            | `Vector{Float64}`                 |
| `ws`                        | `mass_fractions(m)`            | `Vector{Float64}`                 |
| `IDs` + `zs`                | `components(m)`                | `Vector{Pair{String, Float64}}`   |

### Optional constants

| thermo attribute | PyThermo accessor              | Unit  |
|------------------|--------------------------------|-------|
| `Tc`             | `T_critical(c)`                | K     |
| `Pc`             | `P_critical(c)`                | Pa    |
| `omega`          | `acentric_factor(c)`           | —     |
| `Tb`             | `T_boiling(c)`                 | K     |
| `Hvap`           | `enthalpy_vaporization(c)`     | J/kg  |
| `Psat`           | `P_saturation(c)`              | Pa    |

## Internals

The bulk of the curated accessors are generated by a small macro. These
internals are unexported and documented for reference only: external packages
*can* use them to define additional wrappers, but `PropSpec`'s field layout
and the macro's argument signature are implementation details and may change
in any release without a breaking-version bump.

```@docs
PyThermo.PropSpec
PyThermo.@thermo_property
```
