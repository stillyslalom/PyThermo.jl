```@meta
CurrentModule = PyThermo
```

# PyThermo

PyThermo.jl provides Julia access to
[Thermo](https://github.com/CalebBell/thermo), a Python library of
thermophysical property data and correlations covering a large number of
chemical species. It wraps thermo's `Chemical` and `Mixture` objects as the
Julia types [`Species`](@ref) and [`Mixture`](@ref), and adds a curated set of
property accessors that return [Unitful](https://github.com/PainterQubits/Unitful.jl)
quantities in SI units. Properties without a curated accessor remain reachable
as fields of the underlying Python object.

## Installation

PyThermo is registered in Julia's General registry:

```julia
using Pkg
Pkg.add("PyThermo")
```

The Python dependencies are installed automatically through CondaPkg.jl the
first time the package is loaded, as described under
[Interaction with Conda](@ref).

## Quickstart

A `Species` is created by name. Temperature and pressure default to STP and may
be given either as bare numbers in K and Pa or as `Unitful` quantities.

```julia
julia> using PyThermo, Unitful

julia> argon = Species("Ar", T = 300u"K")
Species(Ar, 300.0 K, 1.013e+05 Pa)

julia> density(argon)
1.6227671732556135 kg m^-3
```

Properties without a curated accessor are read straight from the underlying
thermo object:

```julia
julia> argon.molecular_diameter   # Ångström
3.40744
```

A `Mixture` is created from its components and a composition. The same
accessors work on mixtures:

```julia
julia> air = Mixture(["N2" => 0.78, "O2" => 0.21, "Ar" => 0.01])
Mixture(78% N2, 21% O2, 1% Ar, 298.1 K, 1.013e+05 Pa)

julia> soundspeed(air)
346.13450470160916 m s^-1
```

## Guide

- [Mixtures](@ref) covers building mixtures, including combining existing
  species and mixtures and the adiabatic, enthalpy-conserving mode.
- [Property accessors](@ref) lists the curated accessors with their units,
  the strict and optional return conventions, and a thermo-to-PyThermo
  cheatsheet.
- [ShockTube module](@ref) provides one-dimensional gas dynamics for shock
  jumps, driver pressures, and gas-gas Riemann problems.

## Species

```@docs
Species
```

## Interaction with Conda

PyThermo's Python dependencies are managed by CondaPkg.jl, which registers a
project's dependencies in CondaPkg.toml (similar to Julia's Project.toml).
These dependencies are installed automatically in a shared Conda environment
located at `~/.julia/conda_environments/Thermo` when PyThermo is first loaded.
To use a different Conda environment, set the corresponding preference as
described in the
[CondaPkg.jl documentation](https://github.com/JuliaPy/CondaPkg.jl?tab=readme-ov-file#preferences).

## Index

```@index
```
