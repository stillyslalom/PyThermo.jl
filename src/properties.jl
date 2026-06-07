# Property accessor infrastructure for PyThermo.
#
# This file defines `PropSpec`, the `@thermo_property` macro, and small
# helpers used by the curated Julian accessors that wrap properties on
# `thermo.Chemical` / `thermo.Mixture` Python objects.

"""
    PropSpec{U, P}

Describes one curated Julian property accessor.

* `name`        — the Julia function name (e.g. `:density`).
* `pyattr`      — the underlying Python attribute name on `Chemical`/`Mixture`
                  (e.g. `:rho`).
* `unit`        — a `Unitful.Units` value (use `Unitful.NoUnits` for genuinely
                  dimensionless quantities).
* `optional`    — `false` for strict accessors that throw on `None`;
                  `true` for accessors that return `missing` when the underlying
                  Python value is `None`.
* `phase_table` — `nothing`, or a `NamedTuple` mapping phase symbols to
                  per-phase Python attribute names, e.g.
                  `(gas = :rhog, liquid = :rhol, solid = :rhos)`.
* `doc`         — short docstring fragment used as the first line of the
                  generated docstring.

!!! warning "Not part of the public API"
    `PropSpec` is documented so downstream packages can read the field
    layout, but it is an implementation detail. Its fields and type
    parameters may change in any release without a breaking-version bump.
"""
struct PropSpec{U<:Unitful.Units, P<:Union{Nothing, NamedTuple}}
    name::Symbol
    pyattr::Symbol
    unit::U
    optional::Bool
    phase_table::P
    doc::String
end

@inline _wrap_strict(py, unit) = pyconvert(Float64, py) * unit

@inline function _wrap_optional(py, unit)
    val = pyconvert(Union{Float64, Nothing}, py)
    val === nothing ? missing : val * unit
end

# Lookup the per-phase Python attribute name; throw a clear ArgumentError
# if the requested phase isn't supported for this property.
@inline function _phase_attr(name::Symbol, phase::Symbol, table::NamedTuple)
    haskey(table, phase) || throw(ArgumentError(
        "$name does not support phase :$phase (supported: $(keys(table)))"))
    return getfield(table, phase)
end

"""
    @thermo_property name pyattr unit optional phase_table doc

Define a curated Julian accessor that forwards to a Python attribute on a
`Chemical` (Species or Mixture).

* `unit` is a `Unitful.Units` value (use `NoUnits` for dimensionless results).
* `optional` is `false` for strict accessors (throw on `None`) and `true` for
  optional accessors (return `missing` when the Python attribute is `None`).
* `phase_table` is either the literal `nothing` or a `NamedTuple` expression
  like `(gas=:rhog, liquid=:rhol, solid=:rhos)`. When provided, a second
  method dispatching on `(::Chemical, ::Symbol)` is also emitted.
* `doc` is a short docstring fragment used as the function summary.

!!! warning "Not part of the public API"
    This macro is documented for reference, but its argument signature is an
    implementation detail and may change in any release without a
    breaking-version bump. Downstream packages that build their own wrappers
    on it do so at their own risk.
"""
macro thermo_property(name, pyattr, unit, optional, phase_table, doc)
    name_sym   = name isa QuoteNode ? name.value : name
    pyattr_sym = pyattr isa QuoteNode ? pyattr.value : pyattr
    wrap = optional === true || optional === :true ? :_wrap_optional : :_wrap_strict

    base_method = quote
        function $(esc(name_sym))(c::Chemical)
            $wrap(pygetattr(Py(c), $(String(pyattr_sym))), $(esc(unit)))
        end
    end

    if phase_table === :nothing || phase_table === nothing
        phase_method = nothing
    else
        phase_method = quote
            function $(esc(name_sym))(c::Chemical, phase::Symbol)
                attr = _phase_attr($(QuoteNode(name_sym)), phase, $(esc(phase_table)))
                $wrap(pygetattr(Py(c), String(attr)), $(esc(unit)))
            end
        end
    end

    doc_str = _build_property_docstring(name_sym, pyattr_sym, optional, phase_table, doc)

    # `@doc "string" $(esc(name))` fails because the escape wrapper hides the
    # binding from the documentation system. Splice the symbol literal into a
    # `@doc` call and esc the whole macrocall so the name resolves against the
    # caller's module.
    doc_call = esc(:(Core.@doc $doc_str $name_sym))

    out = Expr(:block)
    push!(out.args, base_method)
    phase_method === nothing || push!(out.args, phase_method)
    push!(out.args, doc_call)
    return out
end

function _build_property_docstring(name, pyattr, optional, phase_table, doc)
    opt_note = (optional === true || optional === :true) ?
        "Returns `missing` if the underlying value is unavailable." :
        "Throws if the underlying value is unavailable."

    if phase_table === :nothing || phase_table === nothing
        return string("    ", name, "(c::Chemical)\n\n",
                      doc, "\n\n",
                      "Forwards to the Python attribute `", pyattr, "`. ", opt_note, "\n")
    end

    # phase_table here is the original Expr (NamedTuple constructor like
    # (gas=:rhog, liquid=:rhol, solid=:rhos)). Pull out the keyword pairs
    # so the docstring can name the per-phase attributes.
    phase_lines = String[]
    if phase_table isa Expr && phase_table.head === :tuple
        for arg in phase_table.args
            if arg isa Expr && arg.head === :(=)
                pname = arg.args[1]
                pattr = arg.args[2]
                pattr_sym = pattr isa QuoteNode ? pattr.value : pattr
                push!(phase_lines, "* `:$pname` → `$pattr_sym`")
            end
        end
    end
    phases = isempty(phase_lines) ? "" : ("\n\nPhase-specific Python attributes:\n" * join(phase_lines, "\n"))

    return string("    ", name, "(c::Chemical)\n",
                  "    ", name, "(c::Chemical, phase::Symbol)\n\n",
                  doc, "\n\n",
                  "Without a phase argument, forwards to `", pyattr, "`. ", opt_note,
                  phases, "\n")
end

# ---------------------------------------------------------------------------
# Curated property wrappers
#
# Each `@thermo_property` invocation defines a Julian accessor that forwards
# to a `thermo` Python attribute, attaches a Unitful unit, and (when a phase
# table is supplied) also generates a `(c::Chemical, phase::Symbol)` method
# that dispatches to the matching phase-suffixed attribute.
# ---------------------------------------------------------------------------

# --- State (strict) ---

@thermo_property density rho u"kg/m^3" false (gas=:rhog, liquid=:rhol, solid=:rhos) """
Mass density at the chemical's current state.
"""

@thermo_property molar_density rhom u"mol/m^3" false (gas=:rhogm, liquid=:rholm, solid=:rhosm) """
Molar density at the chemical's current state.
"""

@thermo_property molar_volume Vm u"m^3/mol" false (gas=:Vmg, liquid=:Vml, solid=:Vms) """
Molar volume at the chemical's current state.
"""

# `compressibility` is hand-written in PyThermo.jl (next to `soundspeed`)
# rather than macro-generated: the no-argument `Species` call reads the
# real-gas `Z` from the attached cubic EOS, which the plain attribute forward
# (`Chemical.Z`, curve-based and pinned near 1.0) cannot provide.

# --- Caloric (strict) ---

@thermo_property heat_capacity Cp u"J/(kg*K)" false (gas=:Cpg, liquid=:Cpl, solid=:Cps) """
Isobaric specific heat capacity (mass basis) at the chemical's current state.
"""

@thermo_property molar_heat_capacity Cpm u"J/(mol*K)" false (gas=:Cpgm, liquid=:Cplm, solid=:Cpsm) """
Isobaric molar heat capacity at the chemical's current state.
"""

@thermo_property enthalpy H u"J/kg" false nothing """
Specific enthalpy (mass basis) at the chemical's current state, relative to thermo's
internal reference state.
"""

@thermo_property molar_enthalpy Hm u"J/mol" false nothing """
Molar enthalpy at the chemical's current state, relative to thermo's internal reference state.
"""

@thermo_property entropy S u"J/(kg*K)" false nothing """
Specific entropy (mass basis) at the chemical's current state, relative to thermo's
internal reference state.
"""

@thermo_property molar_entropy Sm u"J/(mol*K)" false nothing """
Molar entropy at the chemical's current state, relative to thermo's internal reference state.
"""

@thermo_property internal_energy U u"J/kg" false nothing """
Specific internal energy (mass basis) at the chemical's current state, relative to thermo's
internal reference state.
"""

@thermo_property molar_internal_energy Um u"J/mol" false nothing """
Molar internal energy at the chemical's current state, relative to thermo's internal reference state.
"""

# --- Transport (strict) ---

@thermo_property viscosity mu u"Pa*s" false (gas=:mug, liquid=:mul) """
Dynamic viscosity at the chemical's current state.
"""

@thermo_property kinematic_viscosity nu u"m^2/s" false (gas=:nug, liquid=:nul) """
Kinematic viscosity (μ / ρ) at the chemical's current state.
"""

@thermo_property thermal_conductivity k u"W/(m*K)" false (gas=:kg, liquid=:kl) """
Thermal conductivity at the chemical's current state.
"""

@thermo_property thermal_diffusivity alpha u"m^2/s" false (gas=:alphag, liquid=:alphal) """
Thermal diffusivity (k / (ρ Cp)) at the chemical's current state.
"""

@thermo_property prandtl Pr NoUnits false (gas=:Prg, liquid=:Prl) """
Prandtl number (Cp μ / k) at the chemical's current state. Dimensionless.
"""

@thermo_property surface_tension sigma u"N/m" false nothing """
Liquid-vapor surface tension at the chemical's current state.
"""

# --- Derived (strict) ---

@thermo_property isobaric_expansion isobaric_expansion u"K^-1" false (gas=:isobaric_expansion_g, liquid=:isobaric_expansion_l) """
Isobaric expansion coefficient β = (1/V)(∂V/∂T)_P at the chemical's current state.
"""

@thermo_property joule_thomson JT u"K/Pa" false (gas=:JTg, liquid=:JTl) """
Joule-Thomson coefficient μ_JT = (∂T/∂P)_H at the chemical's current state.
"""

# --- Optional constants (Union{Quantity, Missing}) ---
#
# These are per-substance constants that may genuinely be absent from thermo's
# database for unusual chemicals. The wrappers return `missing` when the
# Python attribute is `None`. For mixtures, where a single critical-point or
# vaporization value may not be defined, `missing` is also the expected
# return.

@thermo_property T_critical Tc u"K" true nothing """
Critical-point temperature. Returns `missing` if not available.
"""

@thermo_property P_critical Pc u"Pa" true nothing """
Critical-point pressure. Returns `missing` if not available.
"""

@thermo_property acentric_factor omega NoUnits true nothing """
Pitzer acentric factor ω. Dimensionless. Returns `missing` if not available.
"""

@thermo_property T_boiling Tb u"K" true nothing """
Normal boiling point (1 atm). Returns `missing` if not available.
"""

@thermo_property enthalpy_vaporization Hvap u"J/kg" true nothing """
Specific enthalpy of vaporization at the chemical's current temperature.
Returns `missing` if the chemical is above its critical point or has no
vaporization data.
"""

@thermo_property P_saturation Psat u"Pa" true nothing """
Vapor pressure at the chemical's current temperature. Returns `missing` if
the chemical is above its critical point or has no Psat correlation.
"""

# --- Identity (hand-written) ---
#
# These return strings or `Symbol`s rather than `Quantity`s, so the macro's
# `_wrap_strict` plumbing doesn't fit. `molecular_weight` could go through
# the macro, but lives here next to the other identity accessors.

"""
    molecular_weight(c::Chemical) -> Quantity{g/mol}

Molecular weight. For mixtures, this is the mole-fraction-weighted average
of the component molecular weights.

Note: thermo's `MW` attribute is in g/mol, not kg/mol; the returned
`Quantity` preserves that unit.
"""
molecular_weight(c::Chemical) = pyconvert(Float64, Py(c).MW) * u"g/mol"

"""
    CAS(s::Species) -> String

CAS registry number. Defined on `Species` only; for the per-component CAS
numbers of a `Mixture`, use the property fallthrough (`m.CASs`).
"""
CAS(s::Species) = pyconvert(String, Py(s).CAS)

"""
    formula(s::Species) -> String

Empirical chemical formula. Defined on `Species` only.
"""
formula(s::Species) = pyconvert(String, Py(s).formula)

"""
    phase(c::Chemical) -> Symbol

Aggregate phase at the chemical's current state — one of `:gas`, `:liquid`,
`:solid`. Throws `ErrorException` if thermo returns an unrecognized phase
string.
"""
function phase(c::Chemical)
    p = pyconvert(Union{String, Nothing}, Py(c).phase)
    p === nothing && error("thermo reported no phase (got `nothing`)")
    p == "g" ? :gas    :
    p == "l" ? :liquid :
    p == "s" ? :solid  :
    error("unknown phase $(repr(p)) — expected \"g\", \"l\", or \"s\"")
end

# --- Composition (Mixture-only, hand-written) ---

"""
    mole_fractions(m::Mixture) -> Vector{Float64}

Mole fractions of the mixture's components, in the same order as `m.IDs`.
"""
mole_fractions(m::Mixture) = pyconvert(Vector{Float64}, Py(m).zs)

"""
    mass_fractions(m::Mixture) -> Vector{Float64}

Mass fractions of the mixture's components, in the same order as `m.IDs`.
"""
mass_fractions(m::Mixture) = pyconvert(Vector{Float64}, Py(m).ws)

"""
    components(m::Mixture) -> Vector{Pair{String, Float64}}

Mixture composition as a vector of `"ID" => mole_fraction` pairs. The order
matches `m.IDs`.
"""
function components(m::Mixture)
    ids = pyconvert(Vector{String}, Py(m).IDs)
    zs  = pyconvert(Vector{Float64}, Py(m).zs)
    return [id => z for (id, z) in zip(ids, zs)]
end
