using PyThermo
using Test
using PyThermo.ShockTube: shockcalc
using Unitful
using Aqua

@testset "PyThermo.jl" begin

    @testset "Species" begin
        SF6 = Species("SF6")
        @test isapprox(ustrip(density(SF6)), 6.0383, rtol=2e-3)
        SF6.calculate(T = 500)
        @test isapprox(ustrip(density(SF6)), 3.5657, rtol=2e-3)
    end
    @testset "Mixture" begin
        HeAce = Mixture(["Helium" => 0.95, "Acetone" => 0.05])
        ρ_HeAce = HeAce.rho
        @test !isnothing(ρ_HeAce)
        @test isapprox(ρ_HeAce, 0.2747138795604815, rtol=3e-3)
    end
end

@testset "ShockTube.jl" begin
    st = shockcalc(Species("N2"), Species("Ar"), 1.8)
    @test isapprox(ustrip(pressure(st.driver)), 1.426894182009389e6, rtol=2e-3)
    @test isapprox(ustrip(st.u2), 300.0825207636464, rtol=2e-3)
    @test isapprox(ustrip(density(st.shocked)), 3.3867027987446874, rtol=2e-3)
    @test isapprox(ustrip(soundspeed(st.reflected)), 544.5317089952606, rtol=5e-3)
end

@testset "Aqua" begin
    Aqua.test_all(PyThermo; ambiguities=false)
end
