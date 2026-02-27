module ShockTube

using ..PyThermo
using PyThermo: Chemical
using Unitful
using Markdown
using Printf

export shockjump!, shockjump, driverpressure!, driverpressure, shockcalc!, shockcalc, ShockCalc
export riemann_interface, RiemannSolution

# Include the exact Riemann solver at module load time
include("riemann_simple.jl")

γ(f::PyThermo.Chemical)  = isentropic_exponent(f)

"""
    shockjump!(gas::Chemical, Mach::Real) -> (Chemical, Unitful.Velocity)
    shockjump(gas::Chemical, Mach::Real) -> (Chemical, Unitful.Velocity)

Calculate shock jump conditions across a normal shock wave using the Rankine-Hugoniot relations.

The mutating version `shockjump!` modifies the input `gas` object in place, while `shockjump` 
creates a copy and leaves the original unchanged.

# Parameters
- `gas`: A `Species` or `Mixture` object representing the pre-shock gas state
- `Mach`: Shock Mach number (ratio of shock velocity to sound speed in pre-shock gas)

# Returns
A tuple containing:
- Modified `gas` object with post-shock temperature and pressure
- Post-shock gas velocity relative to the laboratory frame (with units `m/s`)

# Physics
The function calculates:
- Pressure ratio (PR) across the shock
- Temperature ratio (TR) across the shock  
- Post-shock gas velocity using the Rankine-Hugoniot jump conditions

# Examples
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
Mixture(95% He, 5% acetone, 291.1 K, 8.500e+04 Pa)

julia> shocked, velocity = shockjump(driven, 2.2);

julia> round(Int, ustrip(u"K", temperature(shocked)))
624

julia> round(Int, ustrip(u"m/s", velocity))
1024
```

# See Also
- [`shockcalc!`](@ref): Comprehensive shock tube calculation
- [`driverpressure!`](@ref): Calculate required driver pressure
"""
function shockjump!(gas, Mach)
    α1 = (γ(gas) + 1) / (γ(gas) - 1)
    PR = (Mach^2 * (1 + α1) - 1) / α1
    TR = PR * (PR + α1) / (1 + α1 * PR)
    u2 = (soundspeed(gas) * (α1 - 1) * (PR - 1) /
          #-----------------------------
           √((1 + α1) * (1 + α1 * PR)))
    gas.P = gas.P * PR
    gas.T = gas.T * TR
    return gas, uconvert(u"m/s", u2)
end

shockjump(gas, Mach) = shockjump!(copy(gas), Mach)

"""
    driverpressure!(driver::Chemical, driven::Chemical, Ms::Real) -> Chemical
    driverpressure(driver::Chemical, driven::Chemical, Ms::Real) -> Chemical

Calculate the required driver gas pressure to achieve a specified shock Mach number in a shock tube.

The mutating version `driverpressure!` modifies the input `driver` object in place, while 
`driverpressure` creates a copy and leaves the original unchanged.

# Parameters
- `driver`: A `Species` or `Mixture` object representing the driver gas (high-pressure section)
- `driven`: A `Species` or `Mixture` object representing the driven gas (low-pressure test section)
- `Ms`: Desired shock Mach number in the driven section

# Returns
The modified `driver` object with updated pressure that will produce the specified shock Mach number.

# Physics
Uses the shock tube theory relating driver and driven gas properties through the contact 
surface and expansion wave. The calculation accounts for:
- Isentropic exponents (γ) of both gases
- Sound speeds in both sections
- The relationship between shock strength and pressure ratio across the contact surface

# Examples
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> driver = Species("He");

julia> driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa");

julia> driver_calc = driverpressure(driver, driven, 2.2);

julia> round(typeof(1.0u"MPa"), pressure(driver_calc), digits=2)
3.73 MPa
```

# See Also
- [`shockjump!`](@ref): Calculate shock jump conditions
- [`shockcalc!`](@ref): Comprehensive shock tube calculation
"""
function driverpressure!(driver, driven, Ms)
    γ1, γ4 = γ(driven), γ(driver)
    a1, a4 = soundspeed(driven), soundspeed(driver)
    driver.P = driven.P * (1 + 2γ1 / (γ1 + 1) * (Ms^2 - 1)) *
                    (1 + a1 / a4 * (γ4 - 1) / (γ1 + 1) * (1/Ms - Ms)) ^ (2γ4 / (1 - γ4))
    return driver
end
driverpressure(driver, driven, Ms) = driverpressure!(copy(driver), driven, Ms)

# Internal function for reflected shock calculations
function shockreflect!(shocked, Ms)
    γ2 = γ(shocked)
    lhs = Ms/(Ms^2 - 1) * sqrt(1 + 2(γ2 - 1)/(γ2 + 1)^2 * (Ms^2 - 1) * (γ2 + 1/Ms^2))
    Mr = (1 + sqrt(1 +4*lhs^2))/(2lhs)
    first(shockjump!(shocked, Mr)), Mr
end
shockreflect(shocked, Ms) = shockreflect!(copy(shocked), Ms)

"""
    ShockCalc

Container struct for comprehensive shock tube calculation results.

This struct stores the complete state information for all regions in a shock tube experiment,
including the driver section, driven section, post-shock (region 2), and reflected shock regions.

# Fields
- `driver::Chemical`: Driver gas state (high-pressure section, region 4)
- `driven::Chemical`: Initial driven gas state (low-pressure test section, region 1)
- `shocked::Chemical`: Post-incident-shock gas state (region 2)
- `reflected::Chemical`: Post-reflected-shock gas state (region 5)
- `Ms::Float64`: Incident shock Mach number
- `Mr::Float64`: Reflected shock Mach number
- `u2::Unitful.Velocity`: Post-shock gas velocity (region 2)

# Display
When displayed, `ShockCalc` presents results in a formatted table showing pressure, temperature,
density, and sound speed for each region, along with incident and reflected wave velocities.

# Examples
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> driver = Species("He");

julia> driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa");

julia> result = shockcalc(driver, driven, 2.2);

julia> result.Ms
2.2

julia> round(typeof(1.0u"kPa"), pressure(result.shocked), digits=1)
481.8 kPa
```

# See Also
- [`shockcalc!`](@ref): Function that creates `ShockCalc` objects
"""
struct ShockCalc{DR<:PyThermo.Chemical,DV<:PyThermo.Chemical,U<:Unitful.Velocity}
    driver::DR
    driven::DV
    shocked::DV
    reflected::DV
    Ms::Float64
    Mr::Float64
    u2::U
end

function showm(io::IO, m::MIME, sc::ShockCalc)
    ush(u, val) = @sprintf("%0.4g", ustrip(u, val))
    printstate(s) = join((ush(u"MPa", pressure(s)),
                          ush(u"K", temperature(s)),
                          ush(u"kg/m^3", density(s)),
                          ush(u"m/s", soundspeed(s))), " | ")
    show(io, m, Markdown.parse("""
    | Region      |  P [MPa]  | T [K] | ρ [kg/m³] | cₛ [m/s] |
    |:----------- | :------------: | :-------------: | :-------------: | :---------------: |
    | Driver      | $(printstate(sc.driver)) |
    | Driven      | $(printstate(sc.driven)) |
    | Shocked     | $(printstate(sc.shocked)) |
    | Reflected   | $(printstate(sc.reflected)) |

    Incident wave: $(ush(u"m/s", sc.Ms*soundspeed(sc.driven))) m/s (Mach $(@sprintf "%0.4g" sc.Ms))\\
    Reflected wave: $(ush(u"m/s", sc.Mr*soundspeed(sc.shocked))) m/s (Mach $(@sprintf "%0.4g" sc.Mr))\\
    Post-shock velocity: $(ush(u"m/s", sc.u2)) m/s
    """))
end

Base.show(io::IO, sc::ShockCalc) = showm(io, MIME("text/plain"), sc)
Base.show(io::IO, m::MIME"text/html", sc::ShockCalc) = showm(io, m, sc)

"""
    shockcalc!(driver::Chemical, driven::Chemical, Ms::Real) -> ShockCalc
    shockcalc(driver::Chemical, driven::Chemical, Ms::Real) -> ShockCalc

Perform a comprehensive shock tube calculation for a given shock Mach number.

This function calculates all relevant states in a shock tube experiment:
- Required driver pressure to achieve the specified Mach number
- Post-shock conditions in the driven section (region 2)
- Reflected shock conditions at the end wall (region 5)
- Wave velocities and gas velocities

The mutating version `shockcalc!` modifies the input `driver` object in place, while 
`shockcalc` creates a copy and leaves the original unchanged.

# Parameters
- `driver`: A `Species` or `Mixture` object representing the driver gas (high-pressure section)
- `driven`: A `Species` or `Mixture` object representing the driven gas (low-pressure test section)
- `Ms`: Desired incident shock Mach number

# Returns
A [`ShockCalc`](@ref) struct containing:
- All four gas states (driver, driven, shocked, reflected)
- Incident and reflected shock Mach numbers
- Post-shock gas velocity

The returned object displays as a formatted table with pressure, temperature, density,
and sound speed for each region.

# Examples
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> driver = Species("He")
Species(He, 298.1 K, 1.013e+05 Pa)

julia> driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
Mixture(95% He, 5% acetone, 291.1 K, 8.500e+04 Pa)

julia> result = shockcalc(driver, driven, 2.2);

julia> round(density(result.shocked) / density(result.driven), digits=2)
2.65
```

# Physics
The calculation follows standard shock tube theory:
1. Determines driver pressure needed for specified shock Mach number
2. Calculates incident shock jump conditions using Rankine-Hugoniot relations
3. Computes reflected shock conditions at the end wall
4. Returns complete state information for analysis

# See Also
- [`shockjump!`](@ref): Calculate shock jump conditions only
- [`driverpressure!`](@ref): Calculate required driver pressure only
- [`ShockCalc`](@ref): Result container struct
"""
function shockcalc!(driver, driven, Ms)
    shocked, u2 = shockjump(driven, Ms)
    driverpressure!(driver, driven, Ms)
    reflected, Mr = shockreflect(shocked, Ms)
    # return (driver = driver, driven = driven, shocked=shocked, reflected=reflected, u2=u2)
    return ShockCalc(driver, driven, shocked, reflected, Ms, Mr, u2)
end
shockcalc(driver, driven, Ms) = shockcalc!(copy(driver), driven, Ms)

"""
    RiemannSolution

Container struct for exact Riemann solver results at a gas interface.

This struct stores the complete solution of the Riemann problem for two gases in contact,
including the star region properties (interface conditions), star densities for both gases,
and all relevant wave speeds. All properties can be accessed using dot notation.

# Fields
- `p_star::Unitful.Pressure`: Pressure at the interface (star region)
- `u_star::Unitful.Velocity`: Velocity at the interface (star region)
- `rho_star_L::Unitful.Density`: Density of left gas in star region
- `rho_star_R::Unitful.Density`: Density of right gas in star region
- `T_star_L::Unitful.Temperature`: Temperature of left gas in star region
- `T_star_R::Unitful.Temperature`: Temperature of right gas in star region
- `S_L::Unitful.Velocity`: Left wave speed (shock or rarefaction head)
- `S_R::Unitful.Velocity`: Right wave speed (shock or rarefaction head)
- `S_contact::Unitful.Velocity`: Contact discontinuity speed
- `gamma_L::Float64`: Isentropic exponent of left gas
- `gamma_R::Float64`: Isentropic exponent of right gas

# Examples
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> left = Species("N2", P=500u"kPa", T=600u"K");

julia> right = Species("SF6", P=100u"kPa", T=300u"K");

julia> sol = riemann_interface(left, right);

julia> round(u"kPa", sol.p_star, digits=2)
319.59 kPa

julia> round(u"m/s", sol.u_star, digits=1)
155.8 m s^-1
```

# See Also
- [`shockcalc!`](@ref): Function that creates `ShockCalc` objects
"""
struct RiemannSolution{P<:Unitful.Pressure, V<:Unitful.Velocity, D<:Unitful.Density, T<:Unitful.Temperature}
    p_star::P
    u_star::V
    rho_star_L::D
    rho_star_R::D
    T_star_L::T
    T_star_R::T
    S_L::V
    S_R::V
    S_contact::V
    gamma_L::Float64
    gamma_R::Float64
end

function Base.show(io::IO, sol::RiemannSolution)
    @printf(io, "RiemannSolution:\n")
    @printf(io, "  Interface: p* = %0.3g %s, u* = %0.3g %s\n",
            ustrip(u"kPa", sol.p_star), "kPa",
            ustrip(u"m/s", sol.u_star), "m/s")
    @printf(io, "  Left  state: ρ = %0.3g kg/m³, T = %0.3g K, S_L = %0.3g m/s\n",
            ustrip(u"kg/m^3", sol.rho_star_L),
            ustrip(u"K", sol.T_star_L),
            ustrip(u"m/s", sol.S_L))
    @printf(io, "  Right state: ρ = %0.3g kg/m³, T = %0.3g K, S_R = %0.3g m/s",
            ustrip(u"kg/m^3", sol.rho_star_R),
            ustrip(u"K", sol.T_star_R),
            ustrip(u"m/s", sol.S_R))
end

"""
    riemann_interface(left::Chemical, right::Chemical, u_left=0.0u"m/s", u_right=0.0u"m/s") -> RiemannSolution
    riemann_interface(left::Chemical, right::Chemical, Ms::Real) -> RiemannSolution

Solve the exact Riemann problem for two gases in contact at an interface.

This function uses the exact Riemann solver to determine the post-shock state at a gas
interface, including pressure, velocity, densities, temperatures, and wave speeds for
both gases. The solution accounts for different isentropic exponents (γ) on each side
of the interface.

# Parameters

## Method 1: Direct specification of states and velocities
- `left::Chemical`: Left gas state (Species or Mixture)
- `right::Chemical`: Right gas state (Species or Mixture)
- `u_left::Unitful.Velocity`: Velocity of left gas (default 0.0u"m/s")
- `u_right::Unitful.Velocity`: Velocity of right gas (default 0.0u"m/s")

## Method 2: Shock tube convenience method
- `left::Chemical`: Unshocked driven gas (will be passed through shock)
- `right::Chemical`: Right gas state (Species or Mixture, assumed at rest)
- `Ms::Real`: Incident shock Mach number

When using Method 2, the function automatically:
1. Applies shock jump conditions to the left gas at the given Mach number
2. Computes the post-shock velocity
3. Solves the Riemann problem with the shocked left state
4. Assumes right gas is at rest (u_right = 0)

# Returns
A [`RiemannSolution`](@ref) struct containing:
- Interface pressure and velocity
- Star region densities and temperatures for both gases
- Wave speeds (left wave, right wave, contact discontinuity)
- Isentropic exponents for reference

# Physics
The exact Riemann solver:
1. Iteratively solves for the star region pressure using Newton-Raphson
2. Computes the star region velocity from momentum conservation
3. Calculates star region densities from Rankine-Hugoniot relations (shock) or isentropic relations (rarefaction)
4. Determines temperatures from the ideal gas law
5. Computes all wave speeds

# Examples

Direct specification of states:
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> left = Species("N2", P=500u"kPa", T=600u"K");

julia> right = Species("SF6", P=100u"kPa", T=300u"K");

julia> sol = riemann_interface(left, right);

julia> round(ustrip(u"kPa", sol.p_star), digits=1)
319.6

julia> round(ustrip(u"m/s", sol.u_star), digits=1)
155.8
```

With velocities:
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> left = Species("N2", P=200u"kPa", T=400u"K");

julia> right = Species("Ar", P=100u"kPa", T=300u"K");

julia> sol = riemann_interface(left, right, 100.0u"m/s", -50.0u"m/s");

julia> sol.u_star > 0.0u"m/s"
true
```

Shock tube convenience method:
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> driven = Species("N2");

julia> test_gas = Species("SF6");

julia> sol = riemann_interface(driven, test_gas, 1.5);

julia> round(ustrip(u"m/s", sol.S_contact), digits=1)
159.2
```

# See Also
- [`RiemannSolution`](@ref): Result container struct
- [`shockjump!`](@ref): Calculate shock jump conditions
"""
function riemann_interface(left::Chemical, right::Chemical, 
                          u_left::Unitful.Velocity=0.0u"m/s", 
                          u_right::Unitful.Velocity=0.0u"m/s")
    # Extract primitive variables
    rho_L = ustrip(u"kg/m^3", density(left))
    p_L = ustrip(u"Pa", pressure(left))
    rho_R = ustrip(u"kg/m^3", density(right))
    p_R = ustrip(u"Pa", pressure(right))
    
    # Convert velocities to m/s
    u_L = ustrip(u"m/s", u_left)
    u_R = ustrip(u"m/s", u_right)
    
    # Get isentropic exponents
    gamma_L = γ(left)
    gamma_R = γ(right)
    
    # Call exact Riemann solver at interface (x/t = 0)
    WL = [rho_L, u_L, p_L]
    WR = [rho_R, u_R, p_R]
    
    _, p_star, u_star, rho_star_L, rho_star_R, S_L, S_R, S_contact, _ = 
        exact_riemann_solver(WL, WR, gamma_L, gamma_R, 0.0)
    
    # Calculate temperatures using ideal gas law: P = ρ R_specific T
    R_L = ustrip(u"J/(kg*K)", R_specific(left))
    R_R = ustrip(u"J/(kg*K)", R_specific(right))
    T_star_L = p_star / (rho_star_L * R_L)
    T_star_R = p_star / (rho_star_R * R_R)
    
    # Return with units
    return RiemannSolution(
        p_star * u"Pa",
        u_star * u"m/s",
        rho_star_L * u"kg/m^3",
        rho_star_R * u"kg/m^3",
        T_star_L * u"K",
        T_star_R * u"K",
        S_L * u"m/s",
        S_R * u"m/s",
        S_contact * u"m/s",
        gamma_L,
        gamma_R
    )
end

"""
    riemann_interface(left::Chemical, right::Chemical, Ms::Real) -> RiemannSolution

Convenience method for shock tube problems: automatically applies shock jump conditions
to the left gas and solves the Riemann problem.

This method assumes:
- `left` is the unshocked driven gas
- `right` is at rest (u_right = 0)
- A shock of Mach number `Ms` passes through the left gas

The function computes the shocked state and velocity of the left gas, then solves
the exact Riemann problem at the interface.

# Parameters
- `left::Chemical`: Unshocked driven gas (will be passed through shock at Ms)
- `right::Chemical`: Test gas at rest
- `Ms::Real`: Incident shock Mach number

# Returns
A [`RiemannSolution`](@ref) with the interface solution

# Examples
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> driven = Species("N2");

julia> test_gas = Species("SF6");

julia> sol = riemann_interface(driven, test_gas, 1.5);

julia> round(ustrip(u"m/s", sol.S_contact), digits=1)
159.2
```
"""
function riemann_interface(left::Chemical, right::Chemical, Ms::Real)
    # Apply shock jump to left gas
    shocked_left, u_left = shockjump(left, Ms)
    
    # Solve Riemann problem with shocked state moving at u_left, right gas at rest
    return riemann_interface(shocked_left, right, u_left, 0.0u"m/s")
end


end # module
