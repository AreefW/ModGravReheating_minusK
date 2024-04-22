/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MODIFIEDGRAVITYCONSTRAINTS_HPP_)
#error                                                                         \
    "This file should only be included through ModifiedGravityConstraints.hpp"
#endif

#ifndef MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_
#define MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class theory_t>
ModifiedGravityConstraints<theory_t>::ModifiedGravityConstraints(
    const theory_t a_theory, double dx,
    const std::array<double, CH_SPACEDIM> a_center, double G_Newton, double rho_mean /*added for rho contrast*/,
    int a_c_Ham, const Interval &a_c_Moms, int a_c_Ham_abs_terms /* defaulted*/,
    const Interval &a_c_Moms_abs_terms /*defaulted*/)
    : Constraints(dx, a_c_Ham, a_c_Moms, a_c_Ham_abs_terms, a_c_Moms_abs_terms,
                  0.0 /*No cosmological constant*/),
      my_theory(a_theory), m_center(a_center), m_G_Newton(G_Newton), m_rho_mean(rho_mean) /*added for rho contrast*/
{
}

template <class theory_t>
template <class data_t>
void ModifiedGravityConstraints<theory_t>::compute(
    Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<BSSNTheoryVars>();
    const auto d1 = m_deriv.template diff1<BSSNTheoryVars>(current_cell);
    const auto d2 = m_deriv.template diff2<BSSNTheoryVars>(current_cell);

    //const auto Tensor<1, data_t> shift = {current_cell.load_vars(c_shift1),current_cell.load_vars(c_shift2),current_cell.load_vars(c_shift3)};
    //----
   // const auto theory_vars = current_cell.template load_vars<Vars>();
    const auto advec =
        this->m_deriv.template advection<BSSNTheoryVars>(current_cell, vars.shift);

    // Inverse metric and Christoffel symbol
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    // Coordinates
    Coordinates<data_t> coords{current_cell, m_deriv.m_dx, m_center};

    // Get the GR terms for the constraints
    Vars<data_t> out = constraint_equations(vars, d1, d2, h_UU, chris);

    // Energy Momentum Tensor
    RhoAndSi<data_t> rho_and_Si =
        my_theory.compute_rho_and_Si(vars, d1, d2, coords);
    SijTFAndS<data_t> Sij_TF_and_S =
        my_theory.compute_Sij_TF_and_S(vars, d1, d2, advec, coords);
    //AllRhos<data_t> all_rhos = my_theory.compute_all_rhos(vars, d1, d2, coords);
    // Hamiltonian constraint
    if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0)
    {
        out.Ham += -16. * M_PI * m_G_Newton * rho_and_Si.rho;
        out.Ham_abs_terms += 16. * M_PI * m_G_Newton * abs(rho_and_Si.rho);
    }

    // Momentum constraints
    if (m_c_Moms.size() > 0 || m_c_Moms_abs_terms.size() > 0)
    {
        FOR(i)
        {
            out.Mom[i] += -8. * M_PI * m_G_Newton * rho_and_Si.Si[i];
            out.Mom_abs_terms[i] +=
                8. * M_PI * m_G_Newton * abs(rho_and_Si.Si[i]);
        }
    }
    out.rho = rho_and_Si.rho;
    //out.rho = all_rhos.phi+all_rhos.g2+all_rhos.g3+all_rhos.GB;
    //out.S = Sij_TF_and_S.S;
    //out.rho_scaled = (all_rhos.phi+all_rhos.g2+all_rhos.g3+all_rhos.GB)/ pow(vars.chi, 3. / 2.);
    out.rho_scaled = rho_and_Si.rho / pow(vars.chi, 3. / 2.);
    current_cell.store_vars(out.rho_scaled, c_rho_scaled);
    //out.S_scaled = Sij_TF_and_S.S / pow(vars.chi, 3. / 2.);
    out.rho_contrast = ((rho_and_Si.rho/m_rho_mean)-1);
    current_cell.store_vars(out.rho_contrast, c_rho_contrast);
    // Write the constraints into the output FArrayBox
    store_vars(out, current_cell);
}

#endif /* MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_ */
