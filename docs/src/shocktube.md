# ShockTube Module

The `ShockTube` module provides tools for analyzing shock tube experiments and solving related gas dynamics problems.

## Overview

Shock tubes are experimental devices used to study high-speed gas dynamics phenomena. The `ShockTube` module provides functionality for:

- Calculating shock jump conditions across normal shocks
- Determining required driver gas pressures for desired shock strengths
- Computing complete shock tube flow solutions
- Solving Riemann problems at material interfaces

## Shock Tube Calculations

```@docs
ShockCalc
shockcalc
shockcalc!
```

## Jump Conditions

```@docs
shockjump
shockjump!
```

## Driver Pressure Calculations

```@docs
driverpressure
driverpressure!
```

## Riemann Interface Solver

```@docs
RiemannSolution
riemann_interface
interface_velocity
interface_pressure
left_wave_speed
right_wave_speed