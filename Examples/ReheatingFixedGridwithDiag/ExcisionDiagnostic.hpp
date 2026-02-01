/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONDIAGNOSTICS_HPP_
#define EXCISIONDIAGNOSTICS_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"

//! Does excision for outer part of object (collect only part within object)

template <class theory_t> class Excision95Density
// class Excision95Density
{
  protected:
    const double m_dx;                                  //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_obj_center; //!< The object center
    const FourthOrderDerivatives m_deriv;
    const double m_obj_rho_max;      //!< Max density
    const double m_obj_bound_cutoff; //!< Define boundary of object as fraction
                                     //!< of max density
    const double m_rho_mean;         //!< mean value of energy density
    const double m_obj_rad;          //!< radius of object

  public:
    Excision95Density(const double a_dx,
                      const std::array<double, CH_SPACEDIM> a_obj_center,
                      const double a_obj_rho_max,
                      const double a_obj_bound_cutoff, const double a_rho_mean,
                      const double a_obj_rad)
        : m_dx(a_dx), m_deriv(m_dx), m_obj_center(a_obj_center),
          m_obj_rho_max(a_obj_rho_max), m_obj_bound_cutoff(a_obj_bound_cutoff),
          m_rho_mean(a_rho_mean), m_obj_rad(a_obj_rad)
    {
    }

    template <class data_t> void compute(const Cell<data_t> current_cell) const
    {
        // const Coordinates<double> coords(current_cell, m_dx, m_center);
        Coordinates<double> coords{current_cell, m_dx, m_obj_center};
        const data_t r = coords.get_radius();
        auto outside_obj = simd_compare_gt(r, m_obj_rad);
        auto rho_now = current_cell.load_vars(c_rho_total);
        auto sqrt_gam_now = current_cell.load_vars(c_sqrt_gam);
        double rho_cutoff = m_obj_bound_cutoff * m_obj_rho_max;

        //  Compute Desity contrast
        auto rho_contrast = (rho_now / m_rho_mean) - 1.0;
        current_cell.store_vars(rho_contrast, c_rho_contrast);
        // rho_now < rho_cutoff ||
        if (outside_obj)
        {
            current_cell.store_vars(0., c_rho_exc);
            current_cell.store_vars(0., c_sqrt_gam_exc);
        }
        else
        {
            current_cell.store_vars(rho_now, c_rho_exc);
            current_cell.store_vars(sqrt_gam_now, c_sqrt_gam_exc);
        }
    }
};

#endif /* EXCISIONDIAGNOSTICS_HPP_ */