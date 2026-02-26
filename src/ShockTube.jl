module ShockTube

using ..PyThermo
using PyThermo: Chemical
using Unitful
using Markdown
using Printf

export shockjump!, shockjump, driverpressure!, driverpressure, shockcalc!, shockcalc, ShockCalc
export riemann_interface, RiemannSolution
export interface_velocity, interface_pressure, left_wave_speed, right_wave_speed

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
Mixture(95.0% helium, 5.00% acetone, 291.1 K, 8.500e+04 Pa)

julia> shocked, velocity = shockjump(driven, 2.2);

julia> round(Int, ustrip(u"K", temperature(shocked)))
624

julia> round(Int, ustrip(u"m/s", velocity))
1024
```

# See Also
- [`shockcalc`](@ref): Comprehensive shock tube calculation
- [`driverpressure`](@ref): Calculate required driver pressure
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

julia> round(pressure(driver_calc), digits=2)
3.73 MPa
```

# See Also
- [`shockjump`](@ref): Calculate shock jump conditions
- [`shockcalc`](@ref): Comprehensive shock tube calculation
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

julia> round(pressure(result.shocked), digits=3)
481.9 kPa
```

# See Also
- [`shockcalc`](@ref): Function that creates `ShockCalc` objects
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
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful), filter = r"\\d\\.\\d{3,}\\d*" => s"\\d\\.\\d{2}"
julia> driver = Species("He")
Species(He, 298.1 K, 1.013e+05 Pa)

julia> driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
Mixture(95.0% helium, 5.00% acetone, 291.1 K, 8.500e+04 Pa)

julia> result = shockcalc(driver, driven, 2.2);

julia> density(result.shocked) / density(result.driven)
2.6407573597520297
```

# Physics
The calculation follows standard shock tube theory:
1. Determines driver pressure needed for specified shock Mach number
2. Calculates incident shock jump conditions using Rankine-Hugoniot relations
3. Computes reflected shock conditions at the end wall
4. Returns complete state information for analysis

# See Also
- [`shockjump`](@ref): Calculate shock jump conditions only
- [`driverpressure`](@ref): Calculate required driver pressure only
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

# Helper functions for Riemann solver

"""
    shock_pressure_velocity(gas::Chemical, p_star::Real) -> Float64

Calculate velocity change across a shock wave for given pressure ratio.
Returns velocity change magnitude in m/s (always positive).
"""
function shock_pressure_velocity(gas::Chemical, p_star::Real)
    p0 = gas.P
    γ = isentropic_exponent(gas)
    a0 = ustrip(u"m/s", soundspeed(gas))
    
    # Pressure ratio
    pr = p_star / p0
    
    # Velocity change across shock (Rankine-Hugoniot)
    # Formula from Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics"
    α = (γ + 1) / (2γ)
    β = (γ - 1) / (2γ)
    du = a0 * (pr - 1) * √(1 / (γ * (β + α * pr)))
    
    return abs(du)
end

"""
    expansion_pressure_velocity(gas::Chemical, p_star::Real) -> Float64

Calculate velocity change across an isentropic expansion wave.
Returns velocity change magnitude in m/s (always positive).
"""
function expansion_pressure_velocity(gas::Chemical, p_star::Real)
    p0 = gas.P
    γ = isentropic_exponent(gas)
    a0 = ustrip(u"m/s", soundspeed(gas))
    
    # Pressure ratio
    pr = p_star / p0
    
    # Velocity change across expansion (isentropic)
    du = 2 * a0 / (γ - 1) * (1 - pr^((γ - 1) / (2γ)))
    
    return abs(du)
end

"""
    wave_velocity(gas::Chemical, p_star::Real, wave_type::Symbol) -> Float64

Calculate wave speed (shock or head of expansion fan) in m/s.
"""
function wave_velocity(gas::Chemical, p_star::Real, wave_type::Symbol)
    p0 = gas.P
    γ = isentropic_exponent(gas)
    a0 = ustrip(u"m/s", soundspeed(gas))
    
    if wave_type == :shock
        # Shock speed
        pr = p_star / p0
        α = (γ + 1) / (γ - 1)
        return a0 * √(1 + α * (pr - 1))
    else  # :expansion
        # Head of expansion fan
        pr = p_star / p0
        return a0 * (1 + (γ - 1) / 2 * (1 - pr^((γ - 1) / (2γ))))
    end
end

"""
    RiemannSolution{L<:PyThermo.Chemical, R<:PyThermo.Chemical}

Container struct for Riemann problem solution at material interfaces.

This struct stores the solution to the Riemann problem when a shock or other
discontinuity encounters an interface between two different gases. It contains
the left and right post-wave states, interface conditions, and wave speeds.

# Fields
- `left_state::L`: Post-reflected-wave gas state (left side)
- `right_state::R`: Post-transmitted-wave gas state (right side)
- `interface_velocity::Float64`: Contact surface velocity [m/s]
- `interface_pressure::Float64`: Pressure at contact surface [Pa]
- `left_wave_speed::Float64`: Reflected wave speed [m/s] (negative if leftward)
- `right_wave_speed::Float64`: Transmitted wave speed [m/s] (positive if rightward)
- `left_wave_type::Symbol`: Type of reflected wave (`:shock` or `:expansion`)
- `right_wave_type::Symbol`: Type of transmitted wave (`:shock` or `:expansion`)

# Accessor Functions
Convenient accessor functions return Unitful quantities:
- `interface_velocity(rs)`: Returns velocity with units
- `interface_pressure(rs)`: Returns pressure with units
- `left_wave_speed(rs)`: Returns wave speed with units
- `right_wave_speed(rs)`: Returns wave speed with units

# Examples
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> air = Mixture(["N2" => 0.78, "O2" => 0.21, "Ar" => 0.01]);

julia> shocked_air, _ = shockjump(air, 2.0);

julia> sf6 = Species("SF6");

julia> sol = riemann_interface(shocked_air, sf6);

julia> sol.left_wave_type
:expansion

julia> sol.right_wave_type
:shock
```

# See Also
- [`riemann_interface`](@ref): Function that creates `RiemannSolution` objects
"""
struct RiemannSolution{L<:PyThermo.Chemical, R<:PyThermo.Chemical}
    left_state::L
    right_state::R
    interface_velocity::Float64      # m/s
    interface_pressure::Float64      # Pa
    left_wave_speed::Float64         # m/s (negative if leftward)
    right_wave_speed::Float64        # m/s (positive if rightward)
    left_wave_type::Symbol           # :shock or :expansion
    right_wave_type::Symbol          # :shock or :expansion
end

# Accessor functions for Unitful quantities
"""
    interface_velocity(rs::RiemannSolution) -> Unitful.Velocity

Get the interface velocity with units from a RiemannSolution.
"""
interface_velocity(rs::RiemannSolution) = rs.interface_velocity * u"m/s"

"""
    interface_pressure(rs::RiemannSolution) -> Unitful.Pressure

Get the interface pressure with units from a RiemannSolution.
"""
interface_pressure(rs::RiemannSolution) = rs.interface_pressure * u"Pa"

"""
    left_wave_speed(rs::RiemannSolution) -> Unitful.Velocity

Get the left wave speed with units from a RiemannSolution.
"""
left_wave_speed(rs::RiemannSolution) = rs.left_wave_speed * u"m/s"

"""
    right_wave_speed(rs::RiemannSolution) -> Unitful.Velocity

Get the right wave speed with units from a RiemannSolution.
"""
right_wave_speed(rs::RiemannSolution) = rs.right_wave_speed * u"m/s"

function Base.show(io::IO, rs::RiemannSolution)
    ush(u, val) = @sprintf("%0.4g", ustrip(u, val))
    printstate(s) = join((ush(u"kPa", pressure(s)),
                          ush(u"K", temperature(s)),
                          ush(u"kg/m^3", density(s)),
                          ush(u"m/s", soundspeed(s))), " | ")
    
    left_dir = rs.left_wave_speed < 0 ? "←" : "→"
    right_dir = rs.right_wave_speed > 0 ? "→" : "←"
    
    println(io, "RiemannSolution:")
    println(io, "  Left wave:  $(rs.left_wave_type) $(left_dir) $(@sprintf("%0.4g", abs(rs.left_wave_speed))) m/s")
    println(io, "  Right wave: $(rs.right_wave_type) $(right_dir) $(@sprintf("%0.4g", abs(rs.right_wave_speed))) m/s")
    println(io, "  Interface:  u* = $(@sprintf("%0.4g", rs.interface_velocity)) m/s, p* = $(@sprintf("%0.4g", rs.interface_pressure/1e3)) kPa")
    println(io, "")
    println(io, "  State       |  P [kPa] | T [K] | ρ [kg/m³] | cₛ [m/s]")
    println(io, "  ──────────────────────────────────────────────────────────")
    println(io, "  Left        | $(printstate(rs.left_state))")
    println(io, "  Right       | $(printstate(rs.right_state))")
end

"""
    riemann_interface(left_gas::Chemical, right_gas::Chemical; 
                     p_star_guess=nothing, max_iter=50, tol=1e-6) -> RiemannSolution

Solve the Riemann problem at a material interface between two gases.

This function solves for the waves (shock or expansion) that arise when two gases at
different states meet at an interface. It's particularly useful for analyzing shock
transmission through material interfaces in shock tubes.

# Parameters
- `left_gas`: Initial state of the left gas (e.g., post-shock gas)
- `right_gas`: Initial state of the right gas (e.g., unshocked test gas)
- `p_star_guess`: Optional initial guess for interface pressure (can be Unitful or Float64 in Pa)
- `max_iter`: Maximum iterations for pressure convergence (default: 50)
- `tol`: Convergence tolerance for velocity matching (default: 1e-6 m/s)

# Returns
A [`RiemannSolution`](@ref) struct containing:
- Post-wave states for both gases
- Interface velocity and pressure
- Wave speeds and types (shock or expansion)

# Physics
The solver iterates to find the interface pressure where:
1. Velocities match at the contact surface (u_left = u_right = u*)
2. Pressures match at the contact surface (p_left = p_right = p*)

Wave types are automatically determined:
- If p* > p0: shock wave
- If p* < p0: expansion wave

# Examples
```jldoctest; setup = :(using PyThermo, PyThermo.ShockTube, Unitful)
julia> air = Mixture(["N2" => 0.78, "O2" => 0.21, "Ar" => 0.01]);

julia> shocked_air, _ = shockjump(air, 2.0);

julia> sf6 = Species("SF6");

julia> sol = riemann_interface(shocked_air, sf6);

julia> round(sol.interface_velocity, digits=1)
304.6

julia> sol.right_wave_type
:shock
```

# See Also
- [`RiemannSolution`](@ref): Result container struct
- [`shockjump`](@ref): Calculate simple shock jump conditions
"""
function riemann_interface(left_gas::Chemical, right_gas::Chemical; 
                          p_star_guess=nothing, max_iter=50, tol=1e-6)
    # Extract initial conditions in MKS
    pL = left_gas.P
    pR = right_gas.P
    uL = 0.0  # Assume stationary in lab frame initially
    uR = 0.0
    
    # Initial guess for interface pressure
    p_star = if p_star_guess isa Unitful.Pressure
        Float64(ustrip(u"Pa", p_star_guess))
    elseif p_star_guess isa Number
        Float64(p_star_guess)
    else
        # Default: geometric mean of pressures
        sqrt(pL * pR)
    end
    
    # Iteratively solve for interface pressure using secant method
    p_old = p_star
    f_old = 0.0
    
    for iter in 1:max_iter
        # Determine wave types
        left_type = p_star > pL ? :shock : :expansion
        right_type = p_star > pR ? :shock : :expansion
        
        # Calculate velocity changes (always positive magnitudes)
        duL = if left_type == :shock
            shock_pressure_velocity(left_gas, p_star)
        else
            expansion_pressure_velocity(left_gas, p_star)
        end
        
        duR = if right_type == :shock
            shock_pressure_velocity(right_gas, p_star)
        else
            expansion_pressure_velocity(right_gas, p_star)
        end
        
        # Velocity at interface from left and right
        # Both gases should converge to the same velocity at the interface
        # For gases initially at rest (uL = uR = 0):
        # - Left: expansion accelerates gas rightward → u_star = uL + duL
        # - Right: shock accelerates gas rightward → u_star = uR + duR
        # Both should be positive and equal at convergence
        u_left = uL + duL
        u_right = uR + duR
        
        # Residual: difference in velocities
        f = u_left - u_right
        
        # Check convergence
        if abs(f) < tol
            # Found solution! Create final states
            left_final = copy(left_gas)
            left_final.P = p_star
            if left_type == :shock
                # Use shock relations
                γL = isentropic_exponent(left_gas)
                pr = p_star / pL
                α = (γL + 1) / (γL - 1)
                tr = pr * (pr + α) / (1 + α * pr)
                left_final.T = left_gas.T * tr
            else
                # Use isentropic relations
                γL = isentropic_exponent(left_gas)
                left_final.T = left_gas.T * (p_star / pL)^((γL - 1) / γL)
            end
            
            right_final = copy(right_gas)
            right_final.P = p_star
            if right_type == :shock
                γR = isentropic_exponent(right_gas)
                pr = p_star / pR
                α = (γR + 1) / (γR - 1)
                tr = pr * (pr + α) / (1 + α * pr)
                right_final.T = right_gas.T * tr
            else
                γR = isentropic_exponent(right_gas)
                right_final.T = right_gas.T * (p_star / pR)^((γR - 1) / γR)
            end
            
            u_star = u_left  # Could also use u_right, they're equal
            
            # Calculate wave speeds in lab frame
            left_speed = if left_type == :shock
                uL - wave_velocity(left_gas, p_star, :shock)
            else
                uL - wave_velocity(left_gas, p_star, :expansion)
            end
            
            right_speed = if right_type == :shock
                uR + wave_velocity(right_gas, p_star, :shock)
            else
                uR + wave_velocity(right_gas, p_star, :expansion)
            end
            
            return RiemannSolution(left_final, right_final, u_star, p_star,
                                  left_speed, right_speed, left_type, right_type)
        end
        
        # Secant method update
        if iter > 1
            dp = -f * (p_star - p_old) / (f - f_old)
            p_old = p_star
            f_old = f
            p_star = p_star + dp
        else
            # First iteration: simple perturbation
            f_old = f
            p_old = p_star
            p_star = p_star * 1.01
        end
        
        # Ensure pressure stays positive and reasonable
        p_star = clamp(p_star, min(pL, pR) * 0.01, max(pL, pR) * 100)
    end
    
    error("Riemann solver failed to converge after $max_iter iterations")
end

end # module
