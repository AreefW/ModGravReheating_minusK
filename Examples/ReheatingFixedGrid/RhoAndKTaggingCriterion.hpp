/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RHOANDKTAGGINGCRITERION_HPP_
#define RHOANDKTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class RhoAndKTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const int m_level;
    const double m_threshold_rho;

  public:
    RhoAndKTaggingCriterion(double dx, const int a_level, double threshold_rho)
        : m_dx(dx), m_deriv(dx), m_level(a_level),
          m_threshold_rho(threshold_rho){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Tensor<1, data_t> d1_rho;
        FOR(idir) m_deriv.diff1(d1_rho, current_cell, idir, c_rho);
        data_t mod_d1_rho = 0;
        FOR(idir) { mod_d1_rho += d1_rho[idir] * d1_rho[idir]; }
        data_t criterion = m_dx * sqrt(mod_d1_rho) / m_threshold_rho;

        // auto rho = current_cell.load_vars(c_deltarho);
        // data_t criterion = rho / (m_threshold_rho * pow(100.0, m_level));

        // pout() << criterion << endl;
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* RHOANDKTAGGINGCRITERION_HPP_ */
