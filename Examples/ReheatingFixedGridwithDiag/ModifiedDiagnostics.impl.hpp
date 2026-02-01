/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MODIFIEDDIAGNOSTICS_HPP_)
#error "This file should only be included through ModifiedDiagnostics.hpp"
#endif

#ifndef MODIFIEDDIAGNOSTICS_IMPL_HPP_
#define MODIFIEDDIAGNOSTICS_IMPL_HPP_

template <class theory_t, class gauge_t, class deriv_t>
ModifiedDiagnostics<theory_t, gauge_t, deriv_t>::ModifiedDiagnostics(
    theory_t a_theory, modified_params_t a_params, gauge_t a_gauge, double a_dx,
    double a_sigma, double a_K_mean, double a_rho_mean, 
    const std::array<double, CH_SPACEDIM> a_center,
    double a_G_Newton)
    : ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>(
          a_theory, a_params, a_gauge, a_dx, a_sigma, a_K_mean, a_rho_mean, a_center, a_G_Newton)
{
}

template <class theory_t, class gauge_t, class deriv_t>
template <class data_t>
void ModifiedDiagnostics<theory_t, gauge_t, deriv_t>::compute(
    Cell<data_t> current_cell) const
{
    // copy data from chombo gridpoint into local variables
    const auto theory_vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        this->m_deriv.template advection<Vars>(current_cell, theory_vars.shift);

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(theory_vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Call CCZ4 RHS - work out GR RHS, no dissipation
    Vars<data_t> theory_rhs;
    this->rhs_equation(theory_rhs, theory_vars, d1, d2, advec);

    Coordinates<data_t> coords{current_cell, this->m_deriv.m_dx,
                               this->m_center};

    // add functions a(x) and b(x) of the modified gauge
    this->template add_a_and_b_rhs(theory_rhs, theory_vars, d1, d2, advec,
                                   coords);

    // add RHS theory terms from EM Tensor
    this->template add_emtensor_rhs(theory_rhs, theory_vars, d1, d2, advec,
                                    coords);

    // add evolution of theory fields themselves
    this->my_theory.add_theory_rhs(theory_rhs, theory_vars, d1, d2, advec,
                                   coords);

    // solve linear system for the theory fields that require it (e.g. 4dST)
    this->my_theory.solve_lhs(theory_rhs, theory_vars, d1, d2, advec, coords);

    data_t RGB = compute_RGB(theory_rhs, theory_vars, d1, d2, advec, coords);

    data_t dfdphi, d2fdphi2, g2, dg2dphi, V_of_phi, dVdphi;

    this->my_theory.my_coupling_and_potential.compute_coupling_and_potential(dfdphi, d2fdphi2,
                                        g2, dg2dphi, V_of_phi, dVdphi, theory_vars, coords);

    data_t my_Veff = V_of_phi-RGB*(this->my_theory.my_coupling_and_potential.get_coupling(theory_vars));
    current_cell.store_vars(my_Veff, c_Veff);

    //Calculate different components of d^2phi/dt^2

    data_t friction_term = theory_vars.lapse * theory_vars.K * theory_vars.Pi;

    data_t gradient_terms = 0.;

    FOR(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        gradient_terms += -h_UU[i][j] * (0.5 * d1.chi[j] * theory_vars.lapse * d1.phi[i] -
                                 theory_vars.chi * theory_vars.lapse * d2.phi[i][j] -
                                 theory_vars.chi * d1.lapse[i] * d1.phi[j]);
        FOR(k)
        {
            gradient_terms += -theory_vars.chi * theory_vars.lapse * h_UU[i][j] * chris.ULL[k][i][j] *
                      d1.phi[k];
        }
    }

    data_t GB_term = dfdphi * RGB * theory_vars.lapse;

    data_t potential_term = -theory_vars.lapse * dVdphi;

    data_t g2_term = theory_rhs.Pi - friction_term - gradient_terms - GB_term - potential_term;

    current_cell.store_vars(friction_term, c_friction);
    current_cell.store_vars(gradient_terms, c_gradient);
    current_cell.store_vars(GB_term, c_GB);
    current_cell.store_vars(potential_term, c_potential);
    current_cell.store_vars(g2_term, c_g2);

}

template <class theory_t, class gauge_t, class deriv_t>
template <class data_t>
data_t ModifiedDiagnostics<theory_t, gauge_t, deriv_t>::
compute_RGB(const Vars<data_t> &theory_rhs,
            const Vars<data_t> &vars,
            const Vars<Tensor<1, data_t>> &d1,
            const Diff2Vars<Tensor<2, data_t>> &d2,
            const Vars<data_t> &advec,
            const Coordinates<data_t> &coords) const
{
    // set the coupling and potential values
    data_t dfdphi = 0.;
    data_t d2fdphi2 = 0.;
    data_t g2 = 0.;
    data_t dg2dphi = 0.;
    data_t V_of_phi = 0.;
    data_t dVdphi = 0.;

    // compute coupling and potential
    this->my_theory.my_coupling_and_potential.compute_coupling_and_potential(
        dfdphi, d2fdphi2, g2, dg2dphi, V_of_phi, dVdphi, vars, coords);

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Compute useful quantities for the Gauss-Bonnet sector

    const data_t chi_regularised = simd_max(1e-6, vars.chi);
    const data_t lapse_regularised = simd_max(1e-6, vars.lapse);

    ScalarVectorTensor<data_t> SVT = this->my_theory.compute_M_Ni_and_Mij(vars, d1, d2);
    data_t M = SVT.scalar;
    Tensor<1, data_t> Ni = SVT.vector;
    Tensor<2, data_t> Mij = SVT.tensor;

    data_t divshift = compute_trace(d1.shift);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse_times_chi;
    FOR(k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR(m) covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m];
        covd2lapse_times_chi[k][l] =
            vars.chi * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                   vars.h[k][l] * dlapse_dot_dchi);
    }
    data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR1(i)
    {
        tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
        FOR1(j)
        {
            tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                           d1.lapse[i] * d1.chi[j]);
        }
    }

    Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);

    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2 = compute_trace(vars.A, A_UU);

    // F_{ij} = \chi{\mathcal L}_nAphys_{ij} + \chi D_iD_j\alpha/\alpha
    //+ A_{ik}A^k_{~j})

    Tensor<2, data_t> Fij;
    FOR(i, j)
    {
        Fij[i][j] =
            (theory_rhs.A[i][j] - advec.A[i][j] + covd2lapse_times_chi[i][j]) /
                lapse_regularised -
            2. / 3. * vars.A[i][j] * (vars.K - divshift / lapse_regularised);
        FOR(k)
        {
            Fij[i][j] -= (vars.A[k][i] * d1.shift[k][j] +
                          vars.A[k][j] * d1.shift[k][i]) /
                         lapse_regularised;
            FOR(l)
            Fij[i][j] += h_UU[k][l] * vars.A[i][k] * vars.A[l][j];
        }
    }

    // F = {\mathcal L}_nK + D^i_D_i\alpha/\alpha - K_{ij}K^{ij}
    data_t F = (theory_rhs.K - advec.K + tr_covd2lapse) / lapse_regularised -
               (tr_A2 + vars.K * vars.K / 3.);

    // other useful quantities
    Tensor<3, data_t> covdtilde_A;
    Tensor<3, data_t> covd_Aphys_times_chi;
    FOR(i, j, k)
    {
        covdtilde_A[j][k][i] = d1.A[j][k][i];
        FOR(l)
        {
            covdtilde_A[j][k][i] += -chris.ULL[l][i][j] * vars.A[l][k] -
                                    chris.ULL[l][i][k] * vars.A[l][j];
        }
        covd_Aphys_times_chi[j][k][i] =
            covdtilde_A[j][k][i] +
            (vars.A[i][k] * d1.chi[j] + vars.A[i][j] * d1.chi[k]) /
                (2. * chi_regularised);
        FOR(l, m)
        {
            covd_Aphys_times_chi[j][k][i] -=
                h_UU[l][m] * d1.chi[m] / (2. * chi_regularised) *
                (vars.h[i][j] * vars.A[k][l] + vars.h[i][k] * vars.A[j][l]);
        }
    }

    Tensor<2, data_t> Mij_TF = Mij;
    make_trace_free(Mij_TF, vars.h, h_UU);
    Tensor<2, data_t> Mij_TF_UU_over_chi =
        raise_all(Mij_TF, h_UU); // raise all indexs
    FOR(i, j) Mij_TF_UU_over_chi[i][j] *= vars.chi;

    // rhs of the Gauss-Bonnet curvature (multiplied by the lapse)
    data_t RGB = -4. / 3. * M * F;
    FOR(i, j)
    {
        RGB += 8. * Mij_TF_UU_over_chi[i][j] * Fij[i][j] +
               16. / 3. * vars.chi * h_UU[i][j] * d1.K[i] *
                   (Ni[j] + 1. / 3. * d1.K[j]) +
               8. * vars.chi * h_UU[i][j] * Ni[i] * Ni[j];
        FOR(k, l, m, n)
        RGB -= 8. * vars.chi * h_UU[i][l] * h_UU[j][m] * h_UU[k][n] *
               covd_Aphys_times_chi[m][n][l] *
               (covd_Aphys_times_chi[j][k][i] - covd_Aphys_times_chi[i][j][k]);
    }

    return RGB;
}

#endif /* MODIFIEDDIAGNOSTICS_IMPL_HPP_ */