/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,

    c_Mom,
    c_Ham_abs_sum,
    c_rho_phi, //<---- c_rho
    c_rho_g2,
    c_rho_g3,
    c_rho_GB,
    c_rho_exc,

    c_sqrt_gam, // sqrt(gamma)=pow(chi,-3/2) volume factor of spatial metric
    c_rho_scaled,
    c_S_scaled,
    c_K_scaled,
    c_A2,
    c_rho_contrast,
    c_sqrt_gam_exc,

    c_weak_coupling_condition_g2,
    c_weak_coupling_condition_g3,
    c_weak_coupling_condition_GB,

    c_friction,
    c_gradient,
    c_GB,
    c_potential,
    c_g2,

    c_Veff,

    c_Discriminant,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom",
    "Ham_abs_sum",

    "rho_phi", "rho_g2", "rho_g3", "rho_GB",
    
    
    "rho_exc", 
    "sqrt_gam", "rho_scaled", "S_scaled", "K_scaled", "A2", 
    "rho_contrast", "sqrt_gam_exc",
    
    "weak_coupling_condition_g2",
    "weak_coupling_condition_g3",
    "weak_coupling_condition_GB",

    "friction", "gradient", "GB", "potential", "g2",

    "Veff", "Discriminant"
    
    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */
