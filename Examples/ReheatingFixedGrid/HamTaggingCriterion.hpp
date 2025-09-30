/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef HAMTAGGINGCRITERION_HPP_
#define HAMTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "Coordinates.hpp"

class HamTaggingCriterion
{
  protected:
    const double m_dx;
    std::array<double, CH_SPACEDIM> m_obj_center;
    double m_rad;

  public:
    HamTaggingCriterion(double dx, std::array<double, CH_SPACEDIM> obj_center,
                        double rad)
        : m_dx(dx), m_obj_center(obj_center), m_rad(rad){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto Ham_abs_sum = current_cell.load_vars(c_Ham_abs_sum);
        auto sqrt_gam = current_cell.load_vars(c_sqrt_gam);

        const Coordinates<data_t> coords(current_cell, m_dx, m_obj_center);
        const data_t r = coords.get_radius();

        data_t criterion = Ham_abs_sum * sqrt_gam * m_dx;
        // pout() << "Criterion = " << criterion << endl;
        auto regrid = simd_compare_gt(r, m_rad);

        // data_t criterion = 0.0;
        // auto regrid = simd_compare_lt(r, m_rad);

        criterion = simd_conditional(regrid, 0.0, criterion);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* HAMTAGGINGCRITERION_HPP_ */
