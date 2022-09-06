module ShockTube

using ..PyThermo
using Unitful
using Markdown
using Printf

export shockjump!, shockjump, driverpressure!, driverpressure, shockcalc!, shockcalc, ShockCalc

γ(f::PyThermo.Chemical)  = isentropic_exponent(f)

# Shock jump conditions
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

function driverpressure!(driver, driven, Ms)
    γ1, γ4 = γ(driven), γ(driver)
    a1, a4 = soundspeed(driven), soundspeed(driver)
    driver.P = driven.P * (1 + 2γ1 / (γ1 + 1) * (Ms^2 - 1)) *
                    (1 + a1 / a4 * (γ4 - 1) / (γ1 + 1) * (1/Ms - Ms)) ^ (2γ4 / (1 - γ4))
    return driver
end
driverpressure(driver, driven, Ms) = driverpressure!(copy(driver), driven, Ms)

function shockreflect!(shocked, Ms)
    γ2 = γ(shocked)
    lhs = Ms/(Ms^2 - 1) * sqrt(1 + 2(γ2 - 1)/(γ2 + 1)^2 * (Ms^2 - 1) * (γ2 + 1/Ms^2))
    Mr = (1 + sqrt(1 +4*lhs^2))/(2lhs)
    first(shockjump!(shocked, Mr)), Mr
end
shockreflect(shocked, Ms) = shockreflect!(copy(shocked), Ms)

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

function shockcalc!(driver, driven, Ms)
    shocked, u2 = shockjump(driven, Ms)
    driverpressure!(driver, driven, Ms)
    reflected, Mr = shockreflect(shocked, Ms)
    # return (driver = driver, driven = driven, shocked=shocked, reflected=reflected, u2=u2)
    return ShockCalc(driver, driven, shocked, reflected, Ms, Mr, u2)
end
shockcalc(driver, driven, Ms) = shockcalc!(copy(driver), driven, Ms)

end # module