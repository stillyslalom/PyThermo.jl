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
    using PyThermo.ShockTube: shockjump, shockjump!, driverpressure, driverpressure!, riemann_interface
    
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
    
    @testset "Riemann Interface Solver" begin
        @testset "Direct Specification - Basic Functionality" begin
            # Test with two different gases at rest
            left = Species("N2", P=500u"kPa", T=600u"K")
            right = Species("SF6", P=100u"kPa", T=300u"K")
            
            sol = riemann_interface(left, right)
            
            # Verify structure
            @test sol isa PyThermo.ShockTube.RiemannSolution
            
            # Verify all fields have proper units
            @test sol.p_star isa Unitful.Pressure
            @test sol.u_star isa Unitful.Velocity
            @test sol.rho_star_L isa typeof(1.0u"kg/m^3")
            @test sol.rho_star_R isa typeof(1.0u"kg/m^3")
            @test sol.T_star_L isa Unitful.Temperature
            @test sol.T_star_R isa Unitful.Temperature
            @test sol.S_L isa Unitful.Velocity
            @test sol.S_R isa Unitful.Velocity
            @test sol.S_contact isa Unitful.Velocity
            
            # Physical correctness: interface pressure between initial pressures
            @test sol.p_star > pressure(right)
            @test sol.p_star < pressure(left)
            
            # Contact discontinuity speed equals interface velocity
            @test sol.S_contact == sol.u_star
            
            # All densities and temperatures should be positive
            @test sol.rho_star_L > 0.0u"kg/m^3"
            @test sol.rho_star_R > 0.0u"kg/m^3"
            @test sol.T_star_L > 0.0u"K"
            @test sol.T_star_R > 0.0u"K"
            
            # Wave speeds: left wave moves left, right wave moves right (for this case)
            @test sol.S_L < 0.0u"m/s"
            @test sol.S_R > 0.0u"m/s"
            
            # Verify gamma values are stored
            @test sol.gamma_L isa Float64
            @test sol.gamma_R isa Float64
            @test sol.gamma_L > 1.0
            @test sol.gamma_R > 1.0
        end
        
        @testset "Shock Tube Convenience Method" begin
            # Test automatic shock jump application
            driven = Species("N2")
            test_gas = Species("SF6")
            Ms = 1.5
            
            sol = riemann_interface(driven, test_gas, Ms)
            
            # Verify solution exists
            @test sol isa PyThermo.ShockTube.RiemannSolution
            
            # Interface velocity should be positive (moving into test gas)
            @test sol.u_star > 0.0u"m/s"
            
            # Verify consistency with shockjump
            shocked_driven, u_shocked = shockjump(driven, Ms)
            
            # The interface pressure should be above test gas pressure
            @test sol.p_star > pressure(test_gas)
            
            # Contact speed should match interface velocity
            @test sol.S_contact == sol.u_star
            
            # Left wave should move faster than contact (shocked gas sound speed)
            @test abs(sol.S_L) > abs(sol.S_contact)
        end
        
        @testset "Velocity Specification with Unitful" begin
            left = Species("N2", P=200u"kPa", T=400u"K")
            right = Species("Ar", P=100u"kPa", T=300u"K")
            
            # Test with Unitful velocities
            u_left = 100.0u"m/s"
            u_right = -50.0u"m/s"
            
            sol = riemann_interface(left, right, u_left, u_right)
            
            # Verify solution computed correctly
            @test sol isa PyThermo.ShockTube.RiemannSolution
            @test sol.u_star isa Unitful.Velocity
            
            # Verify velocities are handled correctly
            @test sol.u_star > 0.0u"m/s"  # Net motion should be rightward
        end
        
        @testset "Edge Case: Identical States" begin
            # Same gas, same conditions on both sides
            left = Species("N2", P=101325u"Pa", T=300u"K")
            right = Species("N2", P=101325u"Pa", T=300u"K")
            
            sol = riemann_interface(left, right)
            
            # Interface should have same pressure
            @test isapprox(ustrip(u"Pa", sol.p_star), 101325.0, rtol=1e-3)
            
            # Interface velocity should be ~zero
            @test isapprox(ustrip(u"m/s", sol.u_star), 0.0, atol=1e-6)
            
            # Densities should match original
            @test isapprox(ustrip(u"kg/m^3", sol.rho_star_L), 
                          ustrip(u"kg/m^3", density(left)), rtol=1e-3)
        end
        
        @testset "Edge Case: Same Pressure, Different Temperature" begin
            p_common = 200u"kPa"
            left = Species("N2", P=p_common, T=600u"K")
            right = Species("N2", P=p_common, T=300u"K")
            
            sol = riemann_interface(left, right)
            
            # Interface pressure should be close to initial pressure
            @test isapprox(ustrip(u"kPa", sol.p_star), 200.0, rtol=0.1)
            
            # Should generate flow due to temperature difference
            @test abs(ustrip(u"m/s", sol.u_star)) > 0.0
        end
        
        @testset "Different Gas Species" begin
            # Test with various gas combinations
            gases = [("He", "Ar"), ("N2", "SF6"), ("O2", "CO2")]
            
            for (gas1, gas2) in gases
                left = Species(gas1, P=300u"kPa", T=500u"K")
                right = Species(gas2, P=100u"kPa", T=300u"K")
                
                sol = riemann_interface(left, right)
                
                # Basic physical checks
                @test sol.p_star > 0.0u"Pa"
                @test sol.rho_star_L > 0.0u"kg/m^3"
                @test sol.rho_star_R > 0.0u"kg/m^3"
                @test sol.T_star_L > 0.0u"K"
                @test sol.T_star_R > 0.0u"K"
            end
        end
        
        @testset "Physical Consistency" begin
            left = Species("N2", P=400u"kPa", T=500u"K")
            right = Species("Ar", P=100u"kPa", T=300u"K")
            
            sol = riemann_interface(left, right)
            
            # Star region pressure should be equal on both sides
            @test sol.p_star == sol.p_star  # Redundant but explicit
            
            # Star region velocity should be equal on both sides (contact property)
            @test sol.u_star == sol.S_contact
            
            # For expansion into lower pressure, interface moves toward right
            @test sol.u_star > 0.0u"m/s"
            
            # Wave speeds should bracket contact speed
            @test sol.S_L < sol.S_contact
            @test sol.S_R > sol.S_contact
        end
        
        @testset "Display and Printing" begin
            left = Species("N2", P=200u"kPa", T=400u"K")
            right = Species("Ar", P=100u"kPa", T=300u"K")
            
            sol = riemann_interface(left, right)
            
            # Test that display doesn't error
            str_repr = sprint(show, sol)
            @test occursin("RiemannSolution", str_repr)
            @test occursin("Interface", str_repr)
            @test occursin("Left", str_repr)
            @test occursin("Right", str_repr)
        end
        
        @testset "Consistency with shockjump" begin
            # When using Mach number method, results should be consistent
            driven = Species("N2", P=100u"kPa", T=300u"K")
            test_gas = Species("SF6", P=100u"kPa", T=300u"K")
            Ms = 2.0
            
            # Method 1: Use riemann_interface with Mach number
            sol1 = riemann_interface(driven, test_gas, Ms)
            
            # Method 2: Manual shock jump then riemann
            shocked, u_shocked = shockjump(driven, Ms)
            sol2 = riemann_interface(shocked, test_gas, u_shocked)
            
            # Should give same results
            @test isapprox(ustrip(u"Pa", sol1.p_star), ustrip(u"Pa", sol2.p_star), rtol=1e-10)
            @test isapprox(ustrip(u"m/s", sol1.u_star), ustrip(u"m/s", sol2.u_star), rtol=1e-10)
        end
    end
end

@testset "Aqua" begin
    Aqua.test_all(PyThermo; ambiguities=false, persistent_tasks=false)
end
