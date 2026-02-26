using PyThermo
using Test
using PyThermo.ShockTube: shockcalc
using Unitful
using Aqua

@testset "PyThermo.jl" begin

    @testset "Species" begin
        SF6 = Species("SF6")
        @test isapprox(ustrip(density(SF6)), 5.9699, rtol=2e-3)
        SF6.calculate(T = 500)
        @test isapprox(ustrip(density(SF6)), 3.5657, rtol=2e-3)
    end
    @testset "Mixture" begin
        HeAce = Mixture(["Helium" => 0.95, "Acetone" => 0.05])
        ρ_HeAce = HeAce.rho
        @test !isnothing(ρ_HeAce)
        @test isapprox(ρ_HeAce, 0.2747138795604815, rtol=3e-3)

        # test mixture display for the case where mole fractions do not sum to 1
        @test occursin("Mixture", repr(Mixture(["N2" => 1.0, "Acetone" => 0.1])))
    end
end

@testset "ShockTube.jl" begin
    st = shockcalc(Species("N2"), Species("Ar"), 1.8)
    @test isapprox(ustrip(pressure(st.driver)), 1.426894182009389e6, rtol=2e-3)
    @test isapprox(ustrip(st.u2), 300.0825207636464, rtol=2e-3)
    @test isapprox(ustrip(density(st.shocked)), 3.3867027987446874, rtol=2e-3)
    @test isapprox(ustrip(soundspeed(st.reflected)), 544.5317089952606, rtol=5e-3)
end

@testset "ShockTube Integration Tests" begin
    using PyThermo.ShockTube: shockjump, shockjump!, driverpressure, driverpressure!
    using PyThermo.ShockTube: riemann_interface, RiemannSolution
    using PyThermo.ShockTube: interface_velocity, interface_pressure, left_wave_speed, right_wave_speed
    
    @testset "Complete Shock Tube Workflow" begin
        # Test the full pipeline from README example
        driver = Species("He")
        driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
        Ms = 2.2
        
        # Full calculation
        result = shockcalc(driver, driven, Ms)
        
        # Verify all regions are present
        @test result.driver isa Species
        @test result.driven isa Mixture
        @test result.shocked isa Mixture
        @test result.reflected isa Mixture
        
        # Verify Mach numbers
        @test result.Ms == Ms
        @test result.Mr > 1.0  # Reflected shock should also be supersonic
        
        # Verify pressure increases through shock
        @test pressure(result.shocked) > pressure(result.driven)
        @test pressure(result.reflected) > pressure(result.shocked)
        @test pressure(result.driver) > pressure(result.shocked)
        
        # Verify temperature increases through shock
        @test temperature(result.shocked) > temperature(result.driven)
        @test temperature(result.reflected) > temperature(result.shocked)
        
        # Verify density ratios (from README)
        @test isapprox(density(result.shocked) / density(result.driven), 2.6407, rtol=5e-3)
        
        # Verify post-shock velocity is positive
        @test ustrip(u"m/s", result.u2) > 0
        @test isapprox(ustrip(u"m/s", result.u2), 1024, rtol=5e-3)
    end
    
    @testset "shockjump Function" begin
        # Test jump conditions calculation
        driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
        Ms = 2.2
        
        # Non-mutating version
        shocked, velocity = shockjump(driven, Ms)
        
        # Verify original is unchanged
        @test isapprox(ustrip(u"K", temperature(driven)), 291.1, rtol=1e-2)
        @test isapprox(ustrip(u"kPa", pressure(driven)), 85.0, rtol=1e-2)
        
        # Verify shocked state (from README)
        @test isapprox(ustrip(u"K", temperature(shocked)), 624, rtol=5e-2)
        @test isapprox(ustrip(u"kPa", pressure(shocked)), 482, rtol=5e-2)
        @test isapprox(ustrip(u"m/s", velocity), 1024, rtol=5e-3)
        
        # Test mutating version
        driven_copy = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
        shocked2, velocity2 = shockjump!(driven_copy, Ms)
        
        # Verify mutation occurred
        @test shocked2 === driven_copy
        @test isapprox(ustrip(u"K", temperature(driven_copy)), 624, rtol=5e-2)
        
        # Verify results match non-mutating version
        @test isapprox(ustrip(u"K", temperature(shocked2)), ustrip(u"K", temperature(shocked)), rtol=1e-6)
        @test isapprox(ustrip(u"m/s", velocity2), ustrip(u"m/s", velocity), rtol=1e-6)
    end
    
    @testset "driverpressure Function" begin
        driver = Species("He")
        driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
        Ms = 2.2
        
        # Non-mutating version
        driver_calc = driverpressure(driver, driven, Ms)
        
        # Verify original is unchanged
        @test isapprox(ustrip(u"Pa", pressure(driver)), 101325, rtol=1e-3)
        
        # Verify calculated pressure (from docstring example)
        @test isapprox(ustrip(u"MPa", pressure(driver_calc)), 3.73, rtol=5e-2)
        
        # Test mutating version
        driver_copy = Species("He")
        driver_calc2 = driverpressure!(driver_copy, driven, Ms)
        
        # Verify mutation occurred
        @test driver_calc2 === driver_copy
        @test isapprox(ustrip(u"MPa", pressure(driver_copy)), 3.73, rtol=5e-2)
    end
    
    @testset "Riemann Interface Solver - Simple Case" begin
        # Test with same species on both sides (should match shock jump conditions)
        # Left: post-shock argon, Right: pre-shock argon
        argon_driven = Species("Ar")
        Ms = 2.0
        argon_shocked, u_shock = shockjump(argon_driven, Ms)
        
        # Create right side at same initial conditions
        argon_right = Species("Ar")
        
        # Solve interface problem
        sol = riemann_interface(argon_shocked, argon_right)
        
        # Verify structure
        @test sol isa RiemannSolution
        @test sol.left_state isa Species
        @test sol.right_state isa Species
        
        # Verify wave types are identified
        @test sol.left_wave_type in (:shock, :expansion)
        @test sol.right_wave_type in (:shock, :expansion)
        
        # For same species: expect expansion on left (reflected), shock on right (transmitted)
        @test sol.left_wave_type == :expansion
        @test sol.right_wave_type == :shock
        
        # Verify interface velocity is positive (gas moving to the right)
        @test sol.interface_velocity > 0
        
        # Verify wave speeds have correct signs
        @test sol.left_wave_speed < 0  # Reflected wave moves left
        @test sol.right_wave_speed > 0  # Transmitted wave moves right
        
        # Verify pressure and velocity continuity at interface (within tolerance)
        @test isapprox(ustrip(u"Pa", pressure(sol.left_state)), sol.interface_pressure, rtol=1e-6)
        @test isapprox(ustrip(u"Pa", pressure(sol.right_state)), sol.interface_pressure, rtol=1e-6)
        
        # Test accessor functions return Unitful quantities
        @test interface_velocity(sol) isa Unitful.Velocity
        @test interface_pressure(sol) isa Unitful.Pressure
        @test left_wave_speed(sol) isa Unitful.Velocity
        @test right_wave_speed(sol) isa Unitful.Velocity
        
        # Verify accessor values match struct fields
        @test ustrip(u"m/s", interface_velocity(sol)) == sol.interface_velocity
        @test ustrip(u"Pa", interface_pressure(sol)) == sol.interface_pressure
        @test ustrip(u"m/s", left_wave_speed(sol)) == sol.left_wave_speed
        @test ustrip(u"m/s", right_wave_speed(sol)) == sol.right_wave_speed
    end
    
    @testset "Riemann Solver with Pressure Guess" begin
        # Test with same species and explicit pressure guess
        argon_driven = Species("Ar")
        argon_shocked, _ = shockjump(argon_driven, 2.0)
        argon_right = Species("Ar")
        
        # Solve with Unitful pressure guess
        p_guess = (pressure(argon_shocked) + pressure(argon_right)) / 2
        sol = riemann_interface(argon_shocked, argon_right, p_star_guess=p_guess)
        
        @test sol isa RiemannSolution
        @test sol.interface_pressure > 0
        
        # Solve with Float64 pressure guess (in Pa)
        sol2 = riemann_interface(argon_shocked, argon_right, p_star_guess=ustrip(u"Pa", p_guess))
        
        @test sol2 isa RiemannSolution
        @test sol2.interface_pressure > 0
    end
    
    @testset "Physical Consistency - Rankine-Hugoniot" begin
        # Test conservation across shock
        driven = Species("Ar")
        Ms = 2.0
        
        shocked, u2 = shockjump(driven, Ms)
        
        # Get properties
        ρ1 = ustrip(u"kg/m^3", density(driven))
        ρ2 = ustrip(u"kg/m^3", density(shocked))
        p1 = ustrip(u"Pa", pressure(driven))
        p2 = ustrip(u"Pa", pressure(shocked))
        a1 = ustrip(u"m/s", soundspeed(driven))
        u2_val = ustrip(u"m/s", u2)
        
        # Shock speed
        us = Ms * a1
        
        # Mass conservation: ρ1 * us = ρ2 * (us - u2)
        @test isapprox(ρ1 * us, ρ2 * (us - u2_val), rtol=1e-3)
        
        # Momentum conservation: p1 + ρ1 * us^2 = p2 + ρ2 * (us - u2)^2
        @test isapprox(p1 + ρ1 * us^2, p2 + ρ2 * (us - u2_val)^2, rtol=1e-3)
    end
    
    @testset "Edge Cases" begin
        # Test weak shock (M ≈ 1)
        driven = Species("N2")
        Ms_weak = 1.1
        
        shocked_weak, _ = shockjump(driven, Ms_weak)
        
        # Properties should be close to initial
        @test isapprox(pressure(shocked_weak) / pressure(driven), 1.0, rtol=0.3)
        @test isapprox(temperature(shocked_weak) / temperature(driven), 1.0, rtol=0.2)
        
        # Test strong shock (M >> 1)
        Ms_strong = 5.0
        shocked_strong, _ = shockjump(driven, Ms_strong)
        
        # Properties should be significantly different
        @test pressure(shocked_strong) / pressure(driven) > 10.0
        @test temperature(shocked_strong) / temperature(driven) > 5.0
    end
    
    @testset "Different Gas Combinations" begin
        # Monatomic gas
        he = Species("He")
        shocked_he, _ = shockjump(he, 2.0)
        @test temperature(shocked_he) > temperature(he)
        
        # Diatomic gas
        n2 = Species("N2")
        shocked_n2, _ = shockjump(n2, 2.0)
        @test temperature(shocked_n2) > temperature(n2)
        
        # Polyatomic mixture
        air_acetone = Mixture(["N2" => 0.78, "O2" => 0.21, "acetone" => 0.01])
        shocked_mix, _ = shockjump(air_acetone, 2.0)
        @test temperature(shocked_mix) > temperature(air_acetone)
        
        # Different gases should have different density ratios for same Mach
        ρ_ratio_he = density(shocked_he) / density(he)
        ρ_ratio_n2 = density(shocked_n2) / density(n2)
        @test abs(ρ_ratio_he - ρ_ratio_n2) > 0.1  # Should be noticeably different
    end
    
    @testset "README Examples Validation" begin
        # Verify the main README example works exactly as documented
        driver = Species("He")
        driven = Mixture(["He" => 0.95, "acetone" => 0.05], T = 18u"°C", P = 85u"kPa")
        
        # shockjump example
        shocked, velocity = shockjump(driven, 2.2)
        @test isapprox(ustrip(u"K", temperature(shocked)), 624, rtol=5e-2)
        @test isapprox(ustrip(u"kPa", pressure(shocked)), 482, rtol=5e-2)
        @test isapprox(ustrip(u"m/s", velocity), 1024, rtol=5e-3)
        
        # shockcalc example
        result = shockcalc(driver, driven, 2.2)
        @test isapprox(density(result.shocked) / density(result.driven), 2.64, rtol=5e-2)
    end
end

@testset "Aqua" begin
    Aqua.test_all(PyThermo; ambiguities=false, persistent_tasks=false)
end
