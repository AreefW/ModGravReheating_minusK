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
// template <class theory_t> class ExcisionDiagnostics
// {
//   protected:
//     const double m_dx;                              //!< The grid spacing
//     const std::array<double, CH_SPACEDIM> m_center; //!< The object center
//     const FourthOrderDerivatives m_deriv;
//     const double m_obj_r; //!< Radius of object

//   public:
//     ExcisionDiagnostics(const double a_dx,
//                         const std::array<double, CH_SPACEDIM> a_center,
//                         const double a_obj_r)
//         : m_dx(a_dx), m_deriv(m_dx), m_center(a_center),
//           m_obj_r(a_obj_r)
//     {
//     }

//     void compute(const Cell<double> current_cell) const
//     {
//         const Coordinates<double> coords(current_cell, m_dx, m_center);
        
//         if (coords.get_radius() < m_obj_r)
//         {
//           auto rho_in = current_cell.load_vars(c_rho_scaled);
//           current_cell.store_vars(rho_in, c_rho_exc);

//           auto chi_in = current_cell.load_vars(c_chi);
//           current_cell.store_vars(chi_in, c_chi_exc);
//             // current_cell.store_vars(0.0, c_rho_GB);
//             // current_cell.store_vars(0.0, c_rho_phi);
//             // current_cell.store_vars(0.0, c_rho_scaled);
//         } // else do nothing
//         else {
//           current_cell.store_vars(0., c_rho_exc);
//           current_cell.store_vars(0., c_chi_exc);
//         }

//       // data_t current_rho;

//       // current_rho = current_cell.load_vars(c_rho_scaled);
//       // if (current_rho > max_rho) {
//       //     current_cell.store_vars(0.0, c_rho_scaled);
//       // }
//     }
// };

template <class theory_t> class Excision95Density
{
  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center >> calculate r from origin
    const FourthOrderDerivatives m_deriv;
    const double m_obj_rho_max; //!< Max density
    const double m_obj_bound_cutoff; //!< Define boundary of object as fraction of max density
    const double m_rho_mean; //!< mean value of energy density

  public:
    Excision95Density(const double a_dx, const std::array<double, CH_SPACEDIM> a_center,
                        const double a_obj_rho_max, const double a_obj_bound_cutoff, const double a_rho_mean)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center), m_obj_rho_max(a_obj_rho_max), m_obj_bound_cutoff(a_obj_bound_cutoff), m_rho_mean(a_rho_mean)
    {
    }


    void compute(const Cell<double> current_cell) const
    {
        //const Coordinates<double> coords(current_cell, m_dx, m_center);
        Coordinates<double> coords{current_cell, m_dx, m_center};
        auto rho_now = current_cell.load_vars(c_rho_scaled);
        //auto chi_now = current_cell.load_vars(c_chi);
        auto sqrt_gam_now = current_cell.load_vars(c_sqrt_gam);
        //int obj_count = 0;
        //pout() << "Rho_now: " << rho_now << std::endl;
        //pout() << "Position of rho_now: " << coords << std::endl;
        //pout() << "Count: " << obj_count << std::endl;
                //  Compute Desity contrast
        auto rho_contrast = (rho_now/m_rho_mean) - 1.0 ;
        current_cell.store_vars(rho_contrast, c_rho_contrast);
        
        double rho_cutoff = m_obj_bound_cutoff * m_obj_rho_max;
        //pout() << "rho95 : " << rho95 << endl;
        //pout() << "rho_now : " << rho_now << endl;
        if (rho_now < rho_cutoff){
          current_cell.store_vars(0., c_rho_exc);
          current_cell.store_vars(0., c_sqrt_gam_exc);
        }
        else {
          current_cell.store_vars(rho_now, c_rho_exc);
          current_cell.store_vars(sqrt_gam_now, c_sqrt_gam_exc);

        //  double diff = abs(rho_now - m_obj_rho_max);

        //   if (diff < 1e-6){
        //   pout() << "Position of rho_max: " << coords << std::endl;
        //   std::array<double, CH_SPACEDIM> pos_rho_max = {coords.x, coords.y, coords.z};
        //   double rad_rho_max = coords.get_radius();
        //   pout() << "Stored position: " << 
        //         "(" << pos_rho_max[0] << "," << pos_rho_max[1] << "," << pos_rho_max[2] << ")" 
        //         << " r: " << rad_rho_max << std::endl;
        //   Coordinates<double> coords_obj{current_cell, m_dx, pos_rho_max};
        // }

        //   double obj_rad = rad_rho_max  
        //   if (){


        //   }
          //pout() << coords << std::endl;

        }
        //obj_count ++;
    }
};


#endif /* EXCISIONDIAGNOSTICS_HPP_ */