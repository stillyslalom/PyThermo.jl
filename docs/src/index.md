```@meta
CurrentModule = PyThermo
```

# PyThermo

Documentation for [PyThermo](https://github.com/stillyslalom/PyThermo.jl).

## Modules

- **[ShockTube](@ref)** - Shock tube analysis and gas dynamics calculations

```@index
```

```@docs
Species
Mixture
```

## Interaction with Conda
PyThermo's Python dependencies are managed by CondaPkg.jl, which registers
a project's dependencies in CondaPkg.toml (similar to Julia's Project.toml).
These dependencies are installed automatically in a shared Conda environment
located at ~/.julia/conda_environments/Thermo when PyThermo is first loaded.
If you'd like to use a different Conda environment, you can set the corresponding preference
as described in the [CondaPkg.jl documentation](https://github.com/JuliaPy/CondaPkg.jl?tab=readme-ov-file#preferences).

