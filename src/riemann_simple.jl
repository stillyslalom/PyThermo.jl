# ============================================================================
# EXACT RIEMANN SOLVER FOR TWO IDEAL GASES WITH DIFFERENT GAMMAS
# ============================================================================

function exact_riemann_solver(WL, WR, gamma_L, gamma_R, x_over_t)
    """
    WL = [rho_L, u_L, p_L]     # Left state
    WR = [rho_R, u_R, p_R]     # Right state
    gamma_L, gamma_R           # Specific heat ratios
    x_over_t                   # Sampling position (x/t for self-similar solution)
    
    Returns: (W_sample, p_star, u_star, rho_star_L, rho_star_R, S_L, S_R, S_contact, species)
    """
    
    # Extract primitive variables
    rho_L, u_L, p_L = WL
    rho_R, u_R, p_R = WR
    
    # Calculate sound speeds
    a_L = sqrt(gamma_L * p_L / rho_L)
    a_R = sqrt(gamma_R * p_R / rho_R)
    
    # Calculate gamma-dependent constants for left state
    G1_L = (gamma_L - 1.0) / (2.0 * gamma_L)
    G2_L = (gamma_L + 1.0) / (2.0 * gamma_L)
    G3_L = 2.0 * gamma_L / (gamma_L - 1.0)
    G4_L = 2.0 / (gamma_L - 1.0)
    G5_L = 2.0 / (gamma_L + 1.0)
    G6_L = (gamma_L - 1.0) / (gamma_L + 1.0)
    G7_L = (gamma_L - 1.0) / 2.0
    
    # Calculate gamma-dependent constants for right state
    G1_R = (gamma_R - 1.0) / (2.0 * gamma_R)
    G2_R = (gamma_R + 1.0) / (2.0 * gamma_R)
    G3_R = 2.0 * gamma_R / (gamma_R - 1.0)
    G4_R = 2.0 / (gamma_R - 1.0)
    G5_R = 2.0 / (gamma_R + 1.0)
    G6_R = (gamma_R - 1.0) / (gamma_R + 1.0)
    G7_R = (gamma_R - 1.0) / 2.0
    
    # ---------------------------------------------------------------------------
    # STEP 1: COMPUTE STAR REGION PRESSURE (p_star) AND VELOCITY (u_star)
    # ---------------------------------------------------------------------------
    
    # Initial guess for p_star (Two-Rarefaction approximation)
    p_star = ((a_L + a_R - G7_L * (u_R - u_L)) / 
              (a_L / p_L^G1_L + a_R / p_R^G1_R))^(1.0/G1_L)
    
    # Ensure positive pressure
    TOL = 1e-10
    p_star = max(TOL, p_star)
    
    # Newton-Raphson iteration
    max_iter = 100
    tolerance = 1e-6
    local f_L, f_R
    
    for iter in 1:max_iter
        # Compute pressure functions and derivatives for left state
        if p_star > p_L  # Shock wave
            A_L = G5_L / rho_L
            B_L = G6_L * p_L
            f_L = (p_star - p_L) * sqrt(A_L / (p_star + B_L))
            df_L = sqrt(A_L / (p_star + B_L)) * 
                   (1.0 - (p_star - p_L) / (2.0 * (B_L + p_star)))
        else  # Rarefaction wave
            f_L = G4_L * a_L * ((p_star / p_L)^G1_L - 1.0)
            df_L = (1.0 / (rho_L * a_L)) * (p_star / p_L)^(-G2_L)
        end
        
        # Compute pressure functions and derivatives for right state
        if p_star > p_R  # Shock wave
            A_R = G5_R / rho_R
            B_R = G6_R * p_R
            f_R = (p_star - p_R) * sqrt(A_R / (p_star + B_R))
            df_R = sqrt(A_R / (p_star + B_R)) * 
                   (1.0 - (p_star - p_R) / (2.0 * (B_R + p_star)))
        else  # Rarefaction wave
            f_R = G4_R * a_R * ((p_star / p_R)^G1_R - 1.0)
            df_R = (1.0 / (rho_R * a_R)) * (p_star / p_R)^(-G2_R)
        end
        
        # Total pressure function and derivative
        f = f_L + f_R + (u_R - u_L)
        df = df_L + df_R
        
        # Newton-Raphson update
        p_star_new = p_star - f / df
        p_star_new = max(TOL, p_star_new)
        
        # Check convergence
        change = 2.0 * abs(p_star_new - p_star) / (p_star_new + p_star)
        if change < tolerance
            break
        end
        
        p_star = p_star_new
    end
    
    # Compute star region velocity
    u_star = 0.5 * (u_L + u_R) + 0.5 * (f_R - f_L)
    
    # ---------------------------------------------------------------------------
    # STEP 2: COMPUTE STAR REGION DENSITIES
    # ---------------------------------------------------------------------------
    
    # Left star region density
    if p_star > p_L  # Shock
        rho_star_L = rho_L * ((p_star / p_L + G6_L) / (G6_L * p_star / p_L + 1.0))
    else  # Rarefaction
        rho_star_L = rho_L * (p_star / p_L)^(1.0 / gamma_L)
    end
    
    # Right star region density
    if p_star > p_R  # Shock
        rho_star_R = rho_R * ((p_star / p_R + G6_R) / (G6_R * p_star / p_R + 1.0))
    else  # Rarefaction
        rho_star_R = rho_R * (p_star / p_R)^(1.0 / gamma_R)
    end
    
    # ---------------------------------------------------------------------------
    # STEP 3: COMPUTE WAVE SPEEDS
    # ---------------------------------------------------------------------------
    
    # Left wave speed
    if p_star > p_L  # Shock
        S_L = u_L - a_L * sqrt(G2_L * p_star / p_L + G1_L)
    else  # Rarefaction (head)
        S_L = u_L - a_L
    end
    
    # Right wave speed
    if p_star > p_R  # Shock
        S_R = u_R + a_R * sqrt(G2_R * p_star / p_R + G1_R)
    else  # Rarefaction (tail)
        S_R = u_R + a_R
    end
    
    # Contact discontinuity speed
    S_contact = u_star
    
    # ---------------------------------------------------------------------------
    # STEP 4: SAMPLE THE SOLUTION AT x/t
    # ---------------------------------------------------------------------------
    
    if x_over_t < S_contact
        # ===== SAMPLE LEFT OF CONTACT =====
        species = "LEFT"
        
        if p_star > p_L
            # ===== LEFT SHOCK =====
            S_L_shock = u_L - a_L * sqrt(G2_L * p_star / p_L + G1_L)
            
            if x_over_t <= S_L_shock
                # Left of shock: original left state
                W_sample = [rho_L, u_L, p_L]
            else
                # Right of shock: star region
                W_sample = [rho_star_L, u_star, p_star]
            end
        else
            # ===== LEFT RAREFACTION =====
            a_star_L = sqrt(gamma_L * p_star / rho_star_L)
            S_head = u_L - a_L
            S_tail = u_star - a_star_L
            
            if x_over_t <= S_head
                # Left of rarefaction fan: original left state
                W_sample = [rho_L, u_L, p_L]
            elseif x_over_t > S_tail
                # Right of rarefaction fan: star region
                W_sample = [rho_star_L, u_star, p_star]
            else
                # Inside rarefaction fan
                rho = rho_L * (G5_L + G6_L / a_L * (u_L - x_over_t))^G4_L
                u = G5_L * (a_L + G7_L * u_L + x_over_t)
                p = p_L * (G5_L + G6_L / a_L * (u_L - x_over_t))^G3_L
                W_sample = [rho, u, p]
            end
        end
        
    else
        # ===== SAMPLE RIGHT OF CONTACT =====
        species = "RIGHT"
        
        if p_star > p_R
            # ===== RIGHT SHOCK =====
            S_R_shock = u_R + a_R * sqrt(G2_R * p_star / p_R + G1_R)
            
            if x_over_t >= S_R_shock
                # Right of shock: original right state
                W_sample = [rho_R, u_R, p_R]
            else
                # Left of shock: star region
                W_sample = [rho_star_R, u_star, p_star]
            end
        else
            # ===== RIGHT RAREFACTION =====
            a_star_R = sqrt(gamma_R * p_star / rho_star_R)
            S_head = u_R + a_R
            S_tail = u_star + a_star_R
            
            if x_over_t >= S_head
                # Right of rarefaction fan: original right state
                W_sample = [rho_R, u_R, p_R]
            elseif x_over_t < S_tail
                # Left of rarefaction fan: star region
                W_sample = [rho_star_R, u_star, p_star]
            else
                # Inside rarefaction fan
                rho = rho_R * (G5_R - G6_R / a_R * (u_R - x_over_t))^G4_R
                u = G5_R * (-a_R + G7_R * u_R + x_over_t)
                p = p_R * (G5_R - G6_R / a_R * (u_R - x_over_t))^G3_R
                W_sample = [rho, u, p]
            end
        end
    end
    
    return W_sample, p_star, u_star, rho_star_L, rho_star_R, S_L, S_R, S_contact, species
end


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

# Initial conditions
# WL = [1.0, 0.0, 1.0]      # Left state: [rho, u, p]
# WR = [0.125, 0.0, 0.1]    # Right state: [rho, u, p]
# gamma_L = 1.4             # Air
# gamma_R = 5.0/3.0         # Monatomic gas

# # Solve at contact discontinuity
# x_over_t = 0.0
# W_sample, p_star, u_star, rho_star_L, rho_star_R, S_L, S_R, S_contact, species = 
#     exact_riemann_solver(WL, WR, gamma_L, gamma_R, x_over_t)

# println("Star region pressure: ", p_star)
# println("Star region velocity: ", u_star)
# println("Left star density: ", rho_star_L)
# println("Right star density: ", rho_star_R)
# println("Left wave speed: ", S_L)
# println("Right wave speed: ", S_R)
# println("Contact speed: ", S_contact)
# println("Sampled state: ", W_sample)
# println("Species at x/t=", x_over_t, ": ", species)

# # ============================================================================
# EXAMPLE: SAMPLE ENTIRE SOLUTION
# ============================================================================

# println("\n--- Sampling entire solution ---")
# x_over_t_values = range(-1.0, 1.0, length=11)

# for xt in x_over_t_values
#     W, _, _, _, _, _, _, _, sp = exact_riemann_solver(WL, WR, gamma_L, gamma_R, xt)
#     println("x/t = ", round(xt, digits=2), " => rho = ", round(W[1], digits=4), 
#             ", u = ", round(W[2], digits=4), ", p = ", round(W[3], digits=4), 
#             " (", sp, ")")
# end

## Test on shock tube problem
# using PyThermo
# using PyThermo.ShockTube
# # Create initial conditions
# N2 = Species("N2")
# Ms = 1.3
# N2_shocked, u2 = shockjump(N2, Ms)
# SF6 = Species("SF6")
# # Sample solution at x/t = 0 (contact discontinuity)
# W_sample, p_star, u_star, rho_star_L, rho_star_R, S_L, S_R, S_contact, species = 
#     exact_riemann_solver([ustrip(u"kg/m^3", density(N2_shocked)), ustrip(u"m/s", u2), ustrip(u"Pa", pressure(N2_shocked))], 
#                          [ustrip(u"kg/m^3", density(SF6)), 0.0, ustrip(u"Pa", pressure(SF6))], 
#                          isentropic_exponent(N2_shocked), isentropic_exponent(SF6), 
#                          0.0)
# # Sample solution across x/t range
# x_over_t_values = range(-500.0, 500.0, length=11)
# for xt in x_over_t_values
#     W, _, _, _, _, _, _, _, sp = exact_riemann_solver([ustrip(u"kg/m^3", density(N2_shocked)), ustrip(u"m/s", u2), ustrip(u"Pa", pressure(N2_shocked))], 
#                                                      [ustrip(u"kg/m^3", density(SF6)), 0.0, ustrip(u"Pa", pressure(SF6))], 
#                                                      isentropic_exponent(N2_shocked), isentropic_exponent(SF6), 
#                                                      xt)
#     println("x/t = ", round(xt, digits=1), " => rho = ", round(W[1], digits=4), 
#             ", u = ", round(W[2], digits=4), ", p = ", round(W[3]/1e3, digits=4), " kPa", 
#             " (", sp, ")")
# end