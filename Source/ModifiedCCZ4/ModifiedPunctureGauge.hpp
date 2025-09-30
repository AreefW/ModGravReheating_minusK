/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODIFIEDPUNCTUREGAUGE_HPP_
#define MODIFIEDPUNCTUREGAUGE_HPP_

// #include "Coordinates"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * This class implements a slightly more generic version of the moving puncture
 * gauge. In particular it uses a Bona-Masso slicing condition of the form
 * f(lapse) = -c*lapse^(p-2)
 * and a Gamma-driver shift condition
 **/

// template <class theory_t>
class ModifiedPunctureGauge
{
  public:
    struct params_t
    {
        // lapse params:
        double lapse_advec_coeff = 0.; //!< Switches advection terms in
                                       //! the lapse condition on/off
        double lapse_power = 1.; //!< The power p in \f$\partial_t \alpha = - c
                                 //!\alpha^p(K-2\Theta)\f$
        double lapse_coeff = 2.; //!< The coefficient c in \f$\partial_t \alpha
                                 //!= -c \alpha^p(K-2\Theta)\f$
        // shift params:
        double shift_Gamma_coeff = 0.75; //!< Gives the F in \f$\partial_t
                                         //!  \beta^i =  F B^i\f$
        double shift_advec_coeff = 0.;   //!< Switches advection terms in the
                                         //! shift condition on/off
        double eta = 1.; //!< The eta in \f$\partial_t B^i = \partial_t \tilde
                         //!\Gamma - \eta B^i\f$
        double a0 = 0.;  //!< constant value of a(x)
        double b0 = 0.;  //!< constant value of b(x)
    };


    // template <class data_t> struct RhoAndSi
    // {
    //   Tensor<1, data_t> Si; //!< S_i = T_ia_n^a
    //   data_t rho;           //!< rho = T_ab n^a n^b
    // };

  protected:
    params_t m_params;

  public:
    ModifiedPunctureGauge(const params_t &a_params) : m_params(a_params) {}

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    inline void rhs_gauge(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                          const vars_t<Tensor<1, data_t>> &d1,
                          const diff2_vars_t<Tensor<2, data_t>> &d2,
                          const vars_t<data_t> &advec,
                          double m_K_mean
                          ) const //---- added m_K_mean
    {
        // RhoAndSi<data_t> rho_and_Si = my_theory.compute_rho_and_Si(theory_vars, d1, d2, coords);
        // auto rho_and_Si = my_theory.compute_rho_and_Si(theory_vars, d1, d2, coords);

        rhs.lapse = m_params.lapse_advec_coeff * advec.lapse;
        // rhs.lapse = m_params.lapse_advec_coeff * advec.lapse -
        //             m_params.lapse_coeff *
        //                 pow(vars.lapse, m_params.lapse_power) *
        //                 (vars.K - pow(24.0 * M_PI * rho_and_Si.rho, 0.5) - 2 * vars.Theta); //--- added -m_K_mean 
        FOR(i)
        {
            rhs.shift[i] = m_params.shift_advec_coeff * advec.shift[i] +
                           m_params.shift_Gamma_coeff * vars.Gamma[i] -
                           m_params.eta * vars.shift[i];
            rhs.B[i] = 0.;
        }
    }

    template <class data_t, template <typename> class coords_t>
    inline void compute_a_and_b(data_t &a_of_x, data_t &b_of_x,
                                const coords_t<data_t> &coords) const
    {
        a_of_x = m_params.a0;
        b_of_x = m_params.b0;
    }

  // Class members
  // theory_t my_theory; //!< The theory object, e.g. 4dST

};

#endif /* MODIFIEDPUNCTUREGAUGE_HPP_ */
