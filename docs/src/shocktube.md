# ShockTube module

The `ShockTube` module provides functions for analyzing shock tube flows and gas-gas interfaces using the Rankine-Hugoniot relations and exact Riemann solvers.

## Overview

Shock tubes are fundamental devices in experimental fluid dynamics, used to generate controlled shock waves and study high-speed gas dynamics. This module provides tools to:

1. Calculate shock jump conditions across moving shock waves
2. Determine required driver section pressures for desired shock strengths
3. Solve complete shock tube problems including reflected shocks
4. Analyze gas-gas interface dynamics using exact Riemann solvers

## Basic shock tube functions

### Shock jump conditions

Calculate the post-shock state for a gas traversed by a shock wave:

```julia
using PyThermo
using PyThermo.ShockTube

# Define the initial state
driven = Species("N2", P=100u"kPa", T=300u"K")

# Calculate shock jump for Mach 2.0 shock
shocked, velocity = shockjump!(driven, 2.0)

println("Post-shock pressure: ", pressure(shocked))
println("Post-shock temperature: ", temperature(shocked))
println("Flow velocity: ", velocity)
```

The `shockjump!` variant modifies the input state in-place:

```julia
driven = Species("N2")
shocked, velocity = shockjump!(driven, 2.0)
# driven is now modified to the shocked state
@assert shocked === driven
```

### Driver pressure calculation

Calculate the required driver section pressure to achieve a desired shock strength:

```julia
driver = Species("He")
driven = Mixture(["He" => 0.95, "acetone" => 0.05], T=18u"°C", P=85u"kPa")
Ms = 2.2

# Calculate required driver pressure
driver_calc = driverpressure!(driver, driven, Ms)
println("Required driver pressure: ", pressure(driver_calc))
# Output: ~3.73 MPa
```

The `driverpressure!` variant modifies the driver state in-place:

```julia
driver = Species("He")
driver_calc = driverpressure!(driver, driven, Ms)
@assert driver_calc === driver
```

### Complete shock tube calculation

Perform a complete shock tube analysis including reflected shock:

```julia
driver = Species("He")
driven = Mixture(["He" => 0.95, "acetone" => 0.05], T=18u"°C", P=85u"kPa")
Ms = 2.2

result = shockcalc!(driver, driven, Ms)

# Access all states
println("Initial driven state: ", result.driven)
println("Shocked state: ", result.shocked)
println("Reflected shock state: ", result.reflected)
println("Driver state: ", result.driver)

# Access flow properties
println("Incident shock Mach: ", result.Ms)
println("Reflected shock Mach: ", result.Mr)
println("Post-shock velocity: ", result.u2)

# Density ratio (useful for test time estimation)
println("Density ratio: ", density(result.shocked) / density(result.driven))
# Output: ~2.64
```

## Riemann interface solver

The exact Riemann solver analyzes the interaction between two gases initially separated by a diaphragm. When the diaphragm is removed (or when a shock reaches the interface), waves propagate into both gases and an interface (contact discontinuity) forms.

### Basic usage: Direct specification

Solve for the interface state between two gases at rest:

```julia
using PyThermo.ShockTube: riemann_interface

# Define left and right states
left = Species("N2", P=500u"kPa", T=600u"K")
right = Species("SF6", P=100u"kPa", T=300u"K")

# Solve the Riemann problem
sol = riemann_interface(left, right)

# Interface properties
println("Interface pressure: ", sol.p_star)
println("Interface velocity: ", sol.u_star)
println("Contact speed: ", sol.S_contact)

# Post-wave states
println("Left star density: ", sol.rho_star_L)
println("Right star density: ", sol.rho_star_R)
println("Left star temperature: ", sol.T_star_L)
println("Right star temperature: ", sol.T_star_R)

# Wave speeds
println("Left wave speed: ", sol.S_L)
println("Right wave speed: ", sol.S_R)

# Gamma values
println("Left gamma: ", sol.gamma_L)
println("Right gamma: ", sol.gamma_R)
```

### With initial velocities

Specify non-zero initial velocities for each gas:

```julia
left = Species("N2", P=200u"kPa", T=400u"K")
right = Species("Ar", P=100u"kPa", T=300u"K")

# Left gas moving right at 100 m/s, right gas moving left at 50 m/s
sol = riemann_interface(left, right, 100.0u"m/s", -50.0u"m/s")

# Or use plain numbers (interpreted as m/s)
sol = riemann_interface(left, right, 100.0, -50.0)
```

### Shock tube convenience method

Automatically apply shock jump and solve interface problem:

```julia
# Define initial states
driven = Species("N2", P=100u"kPa", T=300u"K")
test_gas = Species("SF6", P=100u"kPa", T=300u"K")

# Shock Mach number
Ms = 2.0

# Automatically calculates shock jump in driven gas,
# then solves interface with test gas
sol = riemann_interface(driven, test_gas, Ms)

println("Interface pressure: ", sol.p_star)
println("Interface velocity: ", sol.u_star)
```

This is equivalent to:

```julia
# Manual approach
shocked, u_shocked = shockjump!(driven, Ms)
sol = riemann_interface(shocked, test_gas, u_shocked, 0.0)
```

## Understanding the results

### RiemannSolution structure

The `RiemannSolution` struct contains all information about the solved interface:

**Interface Properties** (same on both sides):
- `p_star`: Interface pressure
- `u_star`: Interface velocity
- `S_contact`: Contact discontinuity speed (equals `u_star`)

**Left Star Region** (between left wave and contact):
- `rho_star_L`: Density
- `T_star_L`: Temperature
- `S_L`: Left wave speed (shock or rarefaction tail)

**Right Star Region** (between contact and right wave):
- `rho_star_R`: Density
- `T_star_R`: Temperature
- `S_R`: Right wave speed (shock or rarefaction tail)

**Original State Information**:
- `gamma_L`: Specific heat ratio of left gas
- `gamma_R`: Specific heat ratio of right gas

### Wave structure

The Riemann solution consists of three waves:

1. **Left Wave** (speed `S_L`): Propagates into left gas
   - If `p_star > p_left`: Shock wave (compression)
   - If `p_star < p_left`: Rarefaction wave (expansion)

2. **Contact Discontinuity** (speed `S_contact = u_star`): Material interface
   - Pressure and velocity continuous across contact
   - Density and temperature discontinuous (different gases)

3. **Right Wave** (speed `S_R`): Propagates into right gas
   - If `p_star > p_right`: Shock wave (compression)
   - If `p_star < p_right`: Rarefaction wave (expansion)

### Physical interpretation

```julia
left = Species("N2", P=400u"kPa", T=500u"K")
right = Species("Ar", P=100u"kPa", T=300u"K")

sol = riemann_interface(left, right)

# High pressure on left causes expansion wave into left gas
# and compression (shock) into right gas
if sol.p_star < pressure(left)
    println("Left wave is rarefaction")
end

if sol.p_star > pressure(right)
    println("Right wave is shock")
end

# Interface moves toward lower pressure region
if sol.u_star > 0.0u"m/s"
    println("Interface moving right")
else
    println("Interface moving left")
end

# Wave speeds bracket the contact
@assert sol.S_L < sol.S_contact < sol.S_R
```

## Display and printing

The `RiemannSolution` struct has a custom display format:

```julia
sol = riemann_interface(left, right)
println(sol)
```

Output:
```
RiemannSolution:
  Interface: p* = 234.5 kPa, u* = 123.4 m/s
  Left:  ρ* = 2.345 kg/m³, T* = 456.7 K, S = -234.5 m/s
  Right: ρ* = 3.456 kg/m³, T* = 345.6 K, S = 345.6 m/s
```

## Examples

### Example 1: Helium-driven shock tube

Complete analysis of a helium-driven shock tube with acetone seeding:

```julia
using PyThermo
using PyThermo.ShockTube

# Define gases
driver = Species("He")
driven = Mixture(["He" => 0.95, "acetone" => 0.05], T=18u"°C", P=85u"kPa")

# Target shock Mach number
Ms = 2.2

# Complete calculation
result = shockcalc!(driver, driven, Ms)

println("=== Shock Tube Results ===")
println("Incident shock Mach: ", result.Ms)
println("Post-shock velocity: ", result.u2)
println("Post-shock pressure: ", pressure(result.shocked))
println("Post-shock temperature: ", temperature(result.shocked))
println("Density ratio: ", density(result.shocked) / density(result.driven))
println("\nReflected shock:")
println("  Mach number: ", result.Mr)
println("  Pressure: ", pressure(result.reflected))
println("  Temperature: ", temperature(result.reflected))
println("\nRequired driver pressure: ", pressure(result.driver))
```

### Example 2: Gas-gas interface analysis

Analyze the interface between shocked nitrogen and sulfur hexafluoride:

```julia
using PyThermo
using PyThermo.ShockTube: riemann_interface

# Initial states
driven = Species("N2", P=100u"kPa", T=300u"K")
test_gas = Species("SF6", P=100u"kPa", T=300u"K")

# Shock strength
Ms = 1.8

# Solve interface problem
sol = riemann_interface(driven, test_gas, Ms)

println("=== Interface Analysis ===")
println("Interface pressure: ", sol.p_star)
println("Interface velocity: ", sol.u_star)
println("\nLeft (shocked N2) properties:")
println("  Density: ", sol.rho_star_L)
println("  Temperature: ", sol.T_star_L)
println("  Wave speed: ", sol.S_L)
println("\nRight (SF6) properties:")
println("  Density: ", sol.rho_star_R)
println("  Temperature: ", sol.T_star_R)
println("  Wave speed: ", sol.S_R)
println("\nContact discontinuity speed: ", sol.S_contact)
```

### Example 3: Different gas combinations

Compare interface behavior for different gas pairs:

```julia
using PyThermo.ShockTube: riemann_interface

# Common left state
left = Species("He", P=300u"kPa", T=500u"K")

# Different right gases
gases = ["Ar", "N2", "SF6", "CO2"]

println("Gas pair comparisons (left = He at 300 kPa, 500 K):")
for gas in gases
    right = Species(gas, P=100u"kPa", T=300u"K")
    sol = riemann_interface(left, right)
    
    println("\nHe / $gas:")
    println("  p* = ", sol.p_star)
    println("  u* = ", sol.u_star)
    println("  ρ*_L = ", sol.rho_star_L)
    println("  ρ*_R = ", sol.rho_star_R)
    println("  Density jump = ", sol.rho_star_L / sol.rho_star_R)
end
```

## API reference

```@docs
shockjump!
driverpressure!
shockcalc!
ShockCalc
riemann_interface
RiemannSolution
```

## Theory background

### Rankine-Hugoniot relations

The shock jump conditions relate pre-shock and post-shock states through conservation of mass, momentum, and energy:

- **Mass**: ρ₁u₁ = ρ₂u₂
- **Momentum**: p₁ + ρ₁u₁² = p₂ + ρ₂u₂²
- **Energy**: h₁ + u₁²/2 = h₂ + u₂²/2

Where subscript 1 denotes pre-shock and subscript 2 denotes post-shock states.

### Riemann problem

The Riemann problem considers the evolution of two initially separated gas regions. The solution consists of three waves:

1. A wave propagating into the left state
2. A contact discontinuity (material interface)
3. A wave propagating into the right state

The interface pressure p* and velocity u* are found by solving:

f_L(p*, u*) = f_R(p*, u*)

where f represents the relations across the left and right waves. The solver iteratively finds p* and u* that satisfy continuity of pressure and velocity at the interface.

## See also

- [PyThermo main documentation](index.md)
- Gas properties: `Species`, `Mixture`
- Property accessors: `pressure`, `temperature`, `density`, `soundspeed`