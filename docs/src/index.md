```@meta
CurrentModule = PyThermo
```

# PyThermo

Documentation for [PyThermo](https://github.com/stillyslalom/PyThermo.jl).

```@index
```

```@docs
Species
Mixture
```

## Interaction with Conda
PyThermo's Python dependencies are managed by CondaPkg.jl, which registers
a project's dependencies in CondaPkg.toml (similar to Julia's Project.toml).
These dependencies are installed automatically when PyThermo is first loaded.
To avoid Conda management overhead during subsequent initialization of PyThermo,
the `JULIA_CONDAPKG_OFFLINE` environment is set to `"true"` by default. This can
be overriden by setting `JULIA_CONDAPKG_OFFLINE` to `"false"` before loading PyThermo.
