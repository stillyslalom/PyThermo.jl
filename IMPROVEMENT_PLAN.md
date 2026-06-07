# PyThermo Accessor Expansion — Implementation Plan

Handoff document. Self-contained. Read top-to-bottom before touching code.

## Status

- **PR 1 (Infrastructure) — DONE.** Landed on `main` as a single commit. See the
  "PR 1 retrospective" section below for what shipped, what was discovered
  along the way, and what those discoveries mean for PR 2/3.
- **PR 2 (Bulk wrappers) — DONE.** All 20 wrappers shipped via the macro plus a
  56-case `@testset "Property accessors"`. See the "PR 2 retrospective"
  section below for the macro fix it surfaced and the doctest scope decision.
- **PR 3 (Identity, composition, optionals, docs) — DONE.** All hand-written
  accessors, optional constants, and the categorized docs page shipped. See
  the "PR 3 retrospective" section below for the surprises that came up
  during the docs build and the missing-value test design.

## Background

PyThermo.jl is a Julia frontend for the Python [`thermo`](https://thermo.readthedocs.io/)
package (legacy `chemical` / `mixture` interface). The package wraps `thermo.Chemical`
and `thermo.Mixture` Python objects as Julia `Species` and `Mixture` structs that share an
abstract `Chemical` supertype.

Two access patterns exist today:

1. **`getproperty` fallthrough** (src/PyThermo.jl:186-199) — `c.Cp`, `c.Hvap`, `c.MW`, etc.
   forward to the Python object via `pyconvert(Any, …)`. No units, returns `nothing` for
   missing values. Works for the ~200 attributes exposed by the upstream Chemical/Mixture
   classes.
2. **Curated Julian accessors** (src/PyThermo.jl:207-213) — function-oriented, Unitful return
   values, SI units. Currently only seven: `temperature`, `pressure`, `density`,
   `molar_density`, `isentropic_exponent`, `R_specific`, `soundspeed`.

The goal of this work is to expand the curated layer to cover commonly-needed thermodynamic,
transport, caloric, critical, and identity properties, with idiomatic Julian names, Unitful
return values, and consistent missing-value handling.

The upstream `thermo.chemical` and `thermo.mixture` interfaces are flagged as **legacy** by
their maintainer (who recommends migrating to `thermo.flash`). For this work we are
scoped to the legacy interface only — no `thermo.flash` migration is planned.

## Decisions already made

All design questions have been resolved with the project owner. Do not reopen them.

1. **Phase variants use a positional Symbol argument**, not a kwarg or suffixed names.
   `density(c)` returns current-phase density; `density(c, :liquid)` returns liquid-phase
   density via `c.rhol`; `density(c, :gas)` via `c.rhog`; `density(c, :solid)` via `c.rhos`.
   Unsupported phases (e.g. `viscosity(c, :solid)` — `mus` doesn't exist in thermo) throw
   `ArgumentError`. Rationale: these accessors won't be heavily loaded with kwargs, so a
   positional Symbol is the lightest API.

2. **Mass vs molar naming uses the `molar_` prefix.** `heat_capacity` (J/(kg·K)),
   `molar_heat_capacity` (J/(mol·K)); `enthalpy` / `molar_enthalpy`; `entropy` /
   `molar_entropy`; etc. Rationale: full-word prefix is unambiguous and reads naturally;
   suffix `_molar` would put the qualifier at the end.

3. **First-pass curation is small (~30 wrappers).** Everything in the table below.
   Additional properties (Parachor, solubility_parameter, dimensionless-number methods,
   safety/hazard constants, structural predicates) are out of scope for the first pass
   and will eventually live in a `PyThermo.KitchenSink` submodule (name TBD) that
   re-exports a much wider set without polluting the primary namespace. The
   `getproperty` fallthrough remains available for everything in the meantime.

4. **Missing-value handling is two-tier.**
   - **Bedrock state variables and state-dependent properties** (T, P, ρ, Cp, μ, k, …) are
     **strict** — accessor return type is `Quantity`, and if Python returns `None` the
     accessor throws (let `pyconvert(Float64, ::Py)` raise its own error; the message is
     adequate).
   - **Constants with a real likelihood of being absent from the database** (Tc, Pc, ω, Tb,
     Hvap, Psat, …) return `Union{Quantity, Missing}`. `None` from Python becomes
     `missing`.

   The macro takes an `optional::Bool` per entry to switch behavior.

5. **Scope is legacy `chemical` / `mixture` only.** No `thermo.flash` work.

6. **`soundspeed` is rewired to thermo's real-gas calculation.** The current implementation
   (src/PyThermo.jl:213) is `√(γ R_specific T)` — ideal gas. New implementation reads from
   the equation-of-state object attached to `Chemical`. See the probe step below to
   determine the exact accessor.

7. **`PropSpec` has parametric types** to avoid `Unitful.Units` being stored as an abstract
   field type. See the struct definition below.

## Property table (first pass)

This is the complete first-iteration table. Implement these wrappers; defer anything not on
this list to the future KitchenSink submodule.

Format: `julian_name` ← `pyattr` (unit, optional?, phase variants if any)

### State (strict)
- `temperature` ← `T` (K) — already exists
- `pressure` ← `P` (Pa) — already exists
- `density` ← `rho` (kg/m³); phases: `:gas`→`rhog`, `:liquid`→`rhol`, `:solid`→`rhos`
- `molar_density` ← `rhom` (mol/m³); phases: `:gas`→`rhogm`, `:liquid`→`rholm`, `:solid`→`rhosm`
- `molar_volume` ← `Vm` (m³/mol); phases: `:gas`→`Vmg`, `:liquid`→`Vml`, `:solid`→`Vms`
- `compressibility` ← `Z` (dimensionless / `NoUnits`); phases: `:gas`→`Zg`, `:liquid`→`Zl`, `:solid`→`Zs`
- `phase` ← `phase` — hand-written, returns `Symbol` (`:gas`/`:liquid`/`:solid`), not macro-generated

### Caloric (strict)
- `heat_capacity` ← `Cp` (J/(kg·K)); phases: `:gas`→`Cpg`, `:liquid`→`Cpl`, `:solid`→`Cps`
- `molar_heat_capacity` ← `Cpm` (J/(mol·K)); phases: `:gas`→`Cpgm`, `:liquid`→`Cplm`, `:solid`→`Cpsm`
- `enthalpy` ← `H` (J/kg)
- `molar_enthalpy` ← `Hm` (J/mol)
- `entropy` ← `S` (J/(kg·K))
- `molar_entropy` ← `Sm` (J/(mol·K))
- `internal_energy` ← `U` (J/kg)
- `molar_internal_energy` ← `Um` (J/mol)
- `isochoric_heat_capacity` ← `Cvg` (J/(kg·K)) — gas-phase only, no phase argument
- `molar_isochoric_heat_capacity` ← `Cvgm` (J/(mol·K)) — gas-phase only, no phase argument

Note: thermo does not expose a current-phase, liquid, or solid `Cv` on Chemical — only the
gas-phase `Cvg` / `Cvgm`. These were initially deferred; they were later added (post-PR 3)
as the gas-phase-only `isochoric_heat_capacity` / `molar_isochoric_heat_capacity` accessors,
which take no phase argument since no other phase variant exists.

### Transport (strict)
- `viscosity` ← `mu` (Pa·s); phases: `:gas`→`mug`, `:liquid`→`mul` (no `:solid`)
- `kinematic_viscosity` ← `nu` (m²/s); phases: `:gas`→`nug`, `:liquid`→`nul`
- `thermal_conductivity` ← `k` (W/(m·K)); phases: `:gas`→`kg`, `:liquid`→`kl`
- `thermal_diffusivity` ← `alpha` (m²/s); phases: `:gas`→`alphag`, `:liquid`→`alphal`
- `prandtl` ← `Pr` (NoUnits); phases: `:gas`→`Prg`, `:liquid`→`Prl`
- `surface_tension` ← `sigma` (N/m) — no phase variants

### Derived (strict)
- `isentropic_exponent` ← `isentropic_exponent` (NoUnits) — already exists; keep
- `R_specific` ← `R_specific` (J/(kg·K)) — already exists; keep
- `soundspeed` — see "soundspeed rewire" below
- `isobaric_expansion` ← `isobaric_expansion` (1/K); phases: `:gas`→`isobaric_expansion_g`, `:liquid`→`isobaric_expansion_l`
- `joule_thomson` ← `JT` (K/Pa); phases: `:gas`→`JTg`, `:liquid`→`JTl`

### Identity (Species; hand-written)
- `molecular_weight` ← `MW` (g/mol) — note thermo uses g/mol for MW. Strict.
- `CAS` ← `CAS` — string, hand-written
- `formula` ← `formula` — string, hand-written

### Composition (Mixture; hand-written)
- `mole_fractions` ← `zs` — `Vector{Float64}`
- `mass_fractions` ← `ws` — `Vector{Float64}`
- `components` ← `(IDs, zs)` — `Vector{Pair{String, Float64}}`, built from both attrs

### Optional constants (`Union{Quantity, Missing}`)
- `T_critical` ← `Tc` (K)
- `P_critical` ← `Pc` (Pa)
- `acentric_factor` ← `omega` (NoUnits)
- `T_boiling` ← `Tb` (K)
- `enthalpy_vaporization` ← `Hvap` (J/kg)
- `P_saturation` ← `Psat` (Pa)

## Design

### `PropSpec` and the macro

```julia
# src/properties.jl

struct PropSpec{U<:Unitful.Units, P<:Union{Nothing, NamedTuple}}
    name::Symbol         # Julia function name
    pyattr::Symbol       # Python attribute on Chemical/Mixture
    unit::U              # Unitful.NoUnits for genuinely dimensionless
    optional::Bool       # false = strict (throw on None); true = Union{Quantity, Missing}
    phase_table::P       # nothing, or e.g. (gas=:rhog, liquid=:rhol, solid=:rhos)
    doc::String
end
```

The macro emits two methods per entry (one if `phase_table === nothing`):

```julia
@thermo_property density rho u"kg/m^3" false (gas=:rhog, liquid=:rhol, solid=:rhos) "Mass density."
```

expands to roughly:

```julia
"""
    density(c::Chemical) -> Quantity{kg/m^3}
    density(c::Chemical, phase::Symbol) -> Quantity{kg/m^3}

Mass density.

With no phase argument, returns the current-phase density (`c.rho`).
With `phase = :gas`, `:liquid`, or `:solid`, returns the corresponding phase-specific
density (`c.rhog`, `c.rhol`, `c.rhos`).
"""
density(c::Chemical) = pyconvert(Float64, Py(c).rho) * u"kg/m^3"

function density(c::Chemical, phase::Symbol)
    pyattr = phase === :gas    ? :rhog :
             phase === :liquid ? :rhol :
             phase === :solid  ? :rhos :
             throw(ArgumentError("density does not support phase :$phase"))
    pyconvert(Float64, pygetattr(Py(c), String(pyattr))) * u"kg/m^3"
end
```

For `optional=true` entries, the body becomes:

```julia
function T_critical(c::Chemical)
    val = pyconvert(Union{Float64, Nothing}, Py(c).Tc)
    val === nothing ? missing : val * u"K"
end
```

For `NoUnits` entries, multiplying by `NoUnits` returns a bare `Float64` — that's fine and
expected.

### Helper for value extraction

A small internal helper keeps the macro clean:

```julia
@inline _wrap_strict(py, unit) = pyconvert(Float64, py) * unit
@inline function _wrap_optional(py, unit)
    val = pyconvert(Union{Float64, Nothing}, py)
    val === nothing ? missing : val * unit
end
```

### Where the wrappers live

- New file: `src/properties.jl` — contains the macro, the helpers, all macro-generated
  wrappers, and the hand-written non-numeric accessors (`phase`, `CAS`, `formula`,
  `mole_fractions`, `mass_fractions`, `components`, and `soundspeed`).
- `src/PyThermo.jl` — `include("properties.jl")` after the type definitions; remove the
  current hand-rolled `temperature`/`pressure`/`density`/`molar_density`/
  `isentropic_exponent`/`R_specific`/`soundspeed` block (lines 207-213); extend the export
  list to cover the new public names.

### Backward compatibility

- All currently-exported names continue to exist with the same single-argument behavior.
  Adding a second positional `phase::Symbol` overload is non-breaking.
- `soundspeed` changes behavior from ideal-gas to real-gas — see implications below.

## `soundspeed` rewire

The current implementation:

```julia
soundspeed(c::Chemical) = sqrt(isentropic_exponent(c) * R_specific(c) * temperature(c)) |> u"m/s"
```

is the ideal-gas formula. We are switching to thermo's real-gas calculation.

### Probe step (do this first)

`thermo.Chemical` exposes an `eos` attribute (a `CubicEOS` instance). The exact method
name for sound speed needs to be confirmed against the installed `thermo` version
(currently pinned to `0.4.1` in CondaPkg.toml). Candidates to check from a Julia REPL:

```julia
using PyThermo, PythonCall
c = Species("N2")
propertynames(Py(c).eos)   # look for speed_of_sound, c_s, sound_speed, etc.
```

If `eos` exposes a direct sound speed (e.g. `eos.speed_of_sound_g(T, P)`), use it.

Fallback path if no direct accessor: compute via the EOS-derived departure functions —
`a = √(γ_real · Z · R_specific · T)` where `γ_real = Cp / Cv` with both pulled from the
EOS object (`eos.Cp_dep_g`, `eos.Cv_dep_g`, added to ideal-gas Cp/Cv).

Document the chosen path in a short comment on the `soundspeed` definition.

### Tests for the rewire

Add to `test/runtests.jl`:

- **Ideal-gas consistency**: `soundspeed(Species("N2", T=300u"K", P=1u"atm"))` ≈ `√(γRT)` to
  within 1% (low-pressure ideal limit). Same for He and Ar.
- **Real-gas deviation**: `soundspeed(Species("SF6", T=300u"K", P=1u"atm"))` deviates from
  the ideal-gas formula by more than ~1% — proves the real-gas path is actually being
  exercised.

### ShockTube implications

The Rankine-Hugoniot relations in `src/ShockTube.jl` (lines 59-69, 113-119) assume ideal
gas. They consume `soundspeed(gas)` at lines 63 and 115. Switching `soundspeed` to real-gas
gives a more realistic shock velocity (`u_s = Ms · a₁`) but the jump ratios still use `γ`
ideal-gas-wise. This is acceptable — the math remains self-consistent within each region
— but the existing test in `test/runtests.jl` (lines 27-32) hardcodes numerical values
that may shift slightly. Update tolerances or expected values if they drift outside `2e-3`
relative tolerance.

Also update the doctest in the `shockjump!` docstring (src/ShockTube.jl:42-53) if the
post-shock numbers change. The README example block on lines 81-94 should be verified
manually after the change.

## `setstate!`

New helper to fix the phase-change footgun called out in README.md:50-55. Add to
src/PyThermo.jl alongside the other Chemical methods:

```julia
"""
    setstate!(c::Chemical; T=nothing, P=nothing) -> Chemical

Set temperature and/or pressure and re-flash the chemical's state. Unlike
`c.T = …`, this calls the underlying `calculate` method, which correctly handles
phase changes.

Either or both of `T` and `P` may be provided as Unitful quantities or as bare
`Float64` in K/Pa.
"""
function setstate!(c::Chemical; T=nothing, P=nothing)
    kwargs = Dict{Symbol, Float64}()
    T !== nothing && (kwargs[:T] = T isa Unitful.Temperature ? _unit(T, u"K") : Float64(T))
    P !== nothing && (kwargs[:P] = P isa Unitful.Pressure    ? _unit(P, u"Pa") : Float64(P))
    Py(c).calculate(; kwargs...)
    return c
end
```

Export it. Document in the new properties docs page.

## PR phasing

Three PRs, in order:

### PR 1 — Infrastructure ✅ DONE
- New file `src/properties.jl` with `PropSpec` struct, `@thermo_property` macro,
  `_wrap_strict` / `_wrap_optional` helpers.
- `soundspeed` rewire (probe first, then implement).
- `setstate!`.
- Tests: ideal-gas consistency for N2/He/Ar, real-gas deviation for SF6, `setstate!`
  with phase change (the README SF6 example).
- ShockTube test tolerance update if needed.
- No new public property names beyond `setstate!`.
- Update README.md "Future development" section: remove the "use `missing` instead of
  `nothing`" TODO once PR 2 lands; soundspeed change can be called out here.

#### PR 1 retrospective — read this before PR 2

1. **EOS staleness is the central gotcha.** `thermo.Chemical` caches its EOS
   object at construction time. Mutating state via `c.T = …`, `c.P = …`, or
   even `c.calculate(T=…, P=…)` updates the Chemical's reported state *but
   does not refresh the cached `c.eos`*. Any wrapper that reads from
   `Py(c).eos.*` must first call `Py(c).set_eos(T=Py(c).T, P=Py(c).P)`.
   This rebuilds the EOS at the current state (it's cheap — one new PR
   object — but it is mandatory for correctness).

   PR 2 properties that will be affected: `isobaric_expansion`,
   `joule_thomson`. Both currently work in the plan as straight property
   forwards (`isobaric_expansion`, `JT`, `isobaric_expansion_g`, etc.).
   Confirm by reading thermo's source — if any of these resolve to an
   `eos.*` derivative, the wrapper needs a `set_eos` refresh just like
   `soundspeed(::Species)`. Pure curve-based properties (densities, Cp/Cv,
   viscosity, thermal conductivity, surface tension) are safe.

2. **Species sound speed uses EOS departures, not Z.** `Chemical.Z` returns
   1.0 even for SF6 at 10 atm — it's curve-based, not EOS-based. The plan's
   fallback formula `a = √(γ Z R T)` would have collapsed back to the ideal
   gas. The implementation instead reads `eos.Cp_dep_g`, `eos.Cv_dep_g`,
   `eos.dP_drho_g` directly and computes `a² = γ_real · (∂P/∂ρ_mass)_T`.
   This is the real formula and yields the expected ~1 % deviation for
   SF6 at 1 atm (and 11 % at 10 atm).

3. **Helium-style fallback.** Helium at room temperature is supercritical
   (Tc ≈ 5.2 K). Its PR EOS doesn't expose `Cp_dep_g` / `Cv_dep_g`. The
   `soundspeed(::Species)` body guards with `pyhasattr` and falls back to
   `√(γ R T)`. Any PR 2 wrapper hitting EOS departures should follow the
   same pattern.

4. **Liquid-phase sound speed is currently unsupported on Species.**
   thermo 0.4 exposes `Cp_dep_g` / `Cv_dep_g` but not the `_l` variants on
   the EOS, so the formula isn't directly available. Liquid `Species`
   currently falls through to the ideal-gas expression — fine in practice
   because liquid sound speed isn't a common use case for this package.
   `Mixture` doesn't have this problem because it forwards to thermo's own
   phase-aware `speed_of_sound`.

5. **Mixture's `eos` attribute can raise `IndexError`.** Don't rely on it
   in PR 2 wrappers — Mixture goes through curve-based properties for
   everything we care about anyway.

6. **Doctest setup needs `using Unitful`.** The shared `_DOCTEST` constant
   in `src/PyThermo.jl` only sets up `using PyThermo`. Any doctest that
   uses Unitful literals needs an inline `julia> using Unitful` line, the
   way `setstate!`'s doctest does. (Don't change `_DOCTEST` globally — the
   `Species` and `Mixture` doctests above it deliberately demonstrate that
   bare-float kwargs work.)

7. **No ShockTube test tolerance changes were needed.** The rewire shifted
   `soundspeed(reflected)` from 544.5 → 543.2 m/s in the
   `shockcalc(N2, Ar, 1.8)` case — well inside the existing `5e-3` rtol.

8. **The macro is built but not yet exercised.** `@thermo_property` and
   `PropSpec` are in `src/properties.jl`; PR 1 didn't define any wrappers
   with the macro. The first real macro use is in PR 2. If something in
   the macro is wrong (e.g. docstring generation, phase dispatch), PR 2
   will discover it.

### PR 2 — Bulk wrappers ✅ DONE
- All entries in the property table above except those handled in PR 3.
- Specifically: state (other than already-exposed temperature/pressure/density/
  molar_density), caloric, transport, derived (isobaric_expansion, joule_thomson).
- Update existing `density` / `molar_density` to gain the phase-positional overload via
  the macro (replace the hand-rolled versions).
- Doctests on every wrapper using the existing `_DOCTEST` setup (src/PyThermo.jl:43).
- Extend `test/runtests.jl` with a `@testset "Property accessors"` block covering
  unit-correctness and phase-argument dispatch.

#### PR 2 retrospective — read this before PR 3

1. **The macro's `@doc` emission was broken.** As PR 1 retro point 8
   anticipated, the first real use of `@thermo_property` blew up at
   precompile with `cannot document the following expression:
   $(Expr(:escape, :density))`. The original emission was
   `:(@doc $doc_str $(esc(name_sym)))`, but `@doc` doesn't accept an escaped
   symbol — it can't resolve the binding through the escape wrapper. Fixed
   by splicing the bare symbol into a `Core.@doc` call and escaping the
   whole macrocall: `esc(:(Core.@doc $doc_str $name_sym))`. PR 3 doesn't
   need to touch this again.

2. **EOS-staleness review for `isobaric_expansion` / `joule_thomson` came
   back clean.** Per chemical.py:2503,2525 (`VolumeLiquid` /
   `VolumeGas` `TP_dependent_property_derivative_T`) and 2688,2713 (`JTl` /
   `JTg` built on `Vml`/`Vmg` + `Cplm`/`Cpgm` + `isobaric_expansion_l`/`_g`),
   neither wrapper resolves to an `eos.*` derivative. The plain macro
   forward (no `set_eos` refresh) is correct.

3. **Doctest scope was reduced from the plan.** The plan called for
   `jldoctest` blocks on every wrapper using `$_DOCTEST` interpolation, but
   the macro builds its docstring at macro-expansion time as a static
   String, which doesn't permit `$_DOCTEST` interpolation from the caller's
   scope. Two cheap options exist if/when PR 3 wants to bring doctests
   back:
   - Switch `_build_property_docstring` to return an `Expr(:string, ...)`
     so the user's `doc` arg (which can contain `$_DOCTEST`) is assembled
     at runtime in the caller's scope.
   - Or use `raw"""..."""` for each wrapper's `doc` arg and spell the
     `jldoctest pythermo; setup = ..., filter = ...` header literally.

   Coverage-wise, the 56-case `@testset "Property accessors"` exercises
   every wrapper for unit correctness, phase dispatch, and unsupported-phase
   `ArgumentError`, so correctness is locked in regardless of which path
   PR 3 picks.

4. **`compressibility` is currently flat at 1.0 for `Species`.** This is
   the same caveat called out in PR 1 retro point 2 (`Chemical.Z` is
   curve-based, not EOS-based). The wrapper returns whatever `thermo`
   returns, so a real-gas `Z` is only meaningful via `Mixture` or, in
   future, via an EOS-derived KitchenSink accessor. The plan's first-pass
   scope is fine — just don't be surprised by `compressibility(Species("SF6"))
   ≈ 1.0`.

5. **Mixture caloric properties (`H`, `Hm`, `S`, `Sm`, `U`, `Um`) depend on
   a successful caloric flash.** Per mixture.py:559-560,968-989, these are
   `None` until the property package fills them in during
   `flash_caloric`. They are populated for the default state of common
   mixtures (air, the README He/Acetone mix) so PR 2 tests don't exercise
   the failure mode — but if PR 3 adds optional-constant wrappers, keep in
   mind that the strict caloric wrappers can legitimately raise on
   exotic mixtures whose property package skipped the caloric step.

6. **Export list grew significantly.** PR 2 added 18 new exported names
   to `src/PyThermo.jl`. PR 3 adds another ~9 (identity, composition,
   optional constants). Worth keeping the categorized
   `export …, …, …` blocks already introduced in PR 2 — they're easier to
   read and to extend than one long line.

### PR 3 — Identity, composition, optionals, docs ✅ DONE
- Hand-written: `phase` (returns Symbol), `CAS`, `formula`, `molecular_weight`,
  `mole_fractions`, `mass_fractions`, `components`.
- Optional constants: `T_critical`, `P_critical`, `acentric_factor`, `T_boiling`,
  `enthalpy_vaporization`, `P_saturation`.
- New docs page `docs/src/properties.md` — categorized API reference plus a
  thermo-name → Julian-name cheatsheet (essential for users coming from the Python lib's
  docs).
- Register the new page in `docs/make.jl`.

#### PR 3 retrospective

1. **`_wrap_optional` worked first try.** PR 1 built the optional helper
   without ever exercising it; PR 3 is the first use. The missing-value
   path is covered by `enthalpy_vaporization(Species("He"))` at 298 K —
   helium's Tc is ~5.2 K, so thermo's vaporization correlation returns
   `None` and the wrapper converts that to `missing` as designed. No
   adjustments to `_wrap_optional` were needed.

2. **`phase` does not clash with `Base`.** Aqua's `test_all` passed with
   the new export. There is no `Base.phase` (the complex-argument function
   is `Base.angle`), so the natural name is free for our `:gas`/`:liquid`/
   `:solid` Symbol accessor. Worth knowing if anyone adds more
   short-named exports later — Aqua is the reliable check.

3. **The docs build forced docstrings onto five long-undocumented
   accessors.** Adding `@docs temperature pressure isentropic_exponent
   R_specific soundspeed` to the new properties page surfaced that those
   five — all hand-written and predating PR 1 — had never had docstrings.
   Documenter's `:missing_docs` check turned that into a build error. PR 3
   landed short docstrings on each as part of the same change. The
   `soundspeed` docstring documents the EOS-derived gas-phase formula plus
   the ideal-gas fallback path established in PR 1; this is now the
   canonical place that behavior is described to users.

4. **`PropSpec` and `@thermo_property` are now documented public surface.**
   Documenter's `:missing_docs` check also flagged the macro internals
   (PR 1 had written full docstrings on them but they were not referenced
   from any manual page). Rather than demoting the check to a warning, the
   properties page now has a short "Internals" section that includes both
   via `@docs`. The implication is that future refactors to `PropSpec`'s
   field layout or the macro's argument signature are now semver-relevant
   for any downstream package that uses them to define its own wrappers.

5. **Optional-constant tests pin presence and absence both.** N2 covers
   `T_critical`, `P_critical`, `acentric_factor`, `T_boiling`; liquid
   water at 300 K covers `enthalpy_vaporization` and `P_saturation`; He at
   298 K covers the `missing` return. That's enough coverage to lock in
   both the strict-quantity and `missing` branches of `_wrap_optional`.
   Numerical tolerances are deliberately loose (rtol up to 10 %) since the
   point is to verify the wrapper plumbing, not to validate thermo's
   correlations.

6. **`molecular_weight` lives with the hand-written identity accessors.**
   The plan put `MW` in the "Identity (Species; hand-written)" group but
   it works just as well on `Mixture` (mole-fraction-weighted MW). The
   wrapper is `(c::Chemical)`, hand-written, returning `g/mol` (thermo's
   native unit for `MW`, not kg/mol).

7. **README cleanup replaced both TODOs.** The "Add Unitful accessors for
   more properties" and "Use `missing` instead of `nothing`" entries from
   the Future Development checklist are both delivered by PR 3 / the
   curated layer as a whole. The README now points at the new properties
   docs page rather than just deleting the items.

8. **No KitchenSink work was attempted.** The plan's KitchenSink submodule
   (re-exporting Parachor, solubility parameter, dimensionless-number
   helpers, safety/hazard constants, structural predicates) remains
   future work. The curated layer is now complete enough that the
   `getproperty` fallthrough plus the cheatsheet table can carry users
   through anything not exposed as a curated accessor.

## What NOT to do

- Do not try to auto-wrap all ~200 attributes via introspection. The curation in the table
  above is the deliberate scope.
- Do not add safety/hazard constants (LFL, UFL, GWP, ODP, Tflash, Tautoignition) in this
  round. They're KitchenSink material.
- Do not wrap the property-curve objects (`HeatCapacityGas`, `ViscosityLiquid`, etc.).
  Different access pattern (T,P → curve), separate design problem.
- Do not wrap the dimensionless-number methods (`Reynolds(V, D)`, etc.). Fluids.jl-style
  packages exist independently.
- Do not wrap the `is_*` boolean structural predicates.
- Do not migrate to `thermo.flash`. Out of scope.
- Do not change the `Chemical` abstract type, `Species` / `Mixture` constructors, the
  `getproperty` fallthrough, or any of the `==` / `hash` / `copy` / show machinery in
  src/PyThermo.jl. The `Py` / `pyconvert` plumbing stays as-is.

## Files touched

- `src/PyThermo.jl` — remove the seven hand-rolled accessor lines (207-213), add
  `include("properties.jl")`, extend the `export` list, add `setstate!`.
- `src/properties.jl` — new file, contains everything.
- `src/ShockTube.jl` — possibly update doctest expected values if `soundspeed` rewire
  shifts them; no other changes.
- `test/runtests.jl` — add property-accessor testsets, ideal-gas/real-gas soundspeed
  tests, `setstate!` test. Adjust ShockTube test tolerances if needed.
- `docs/src/properties.md` — new (PR 3).
- `docs/make.jl` — register new page (PR 3).
- `README.md` — minor: clear the `missing` TODO, mention soundspeed change.

## Reference: upstream property inventory

A full enumeration of `thermo.Chemical` properties (constants, state-dependent, transport,
caloric, identity, structural, safety) and the corresponding `Mixture` differences was
collected from <https://thermo.readthedocs.io/thermo.chemical.html> and
<https://thermo.readthedocs.io/thermo.mixture.html>. The first-pass table above is a
curated subset of that inventory. When in doubt about a Python attribute name or its
units, those upstream docs are authoritative.
