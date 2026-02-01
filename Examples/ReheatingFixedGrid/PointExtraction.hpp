/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POINTEXTRACTION_HPP_
#define POINTEXTRACTION_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SimulationParametersBase.hpp"
#include "SmallDataIO.hpp"
#include "SphericalHarmonics.hpp"
#include "UserVariables.hpp"
#include <fstream>
#include <iostream>

//!  The class allows extraction of the values of components at
//!  specified points
class PointExtraction
{
  private:
    //! Params for extraction
    const int m_comp1;
    const int m_comp2;
    const int m_num_points;
    const double m_L;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_dt;
    const double m_time;

  public:
    //! The constructor
    PointExtraction(int a_comp1, int a_comp2, int a_num_points, double a_L,
                    std::array<double, CH_SPACEDIM> a_center, double a_dt,
                    double a_time)
        : m_comp1(a_comp1), m_comp2(a_comp2), m_num_points(a_num_points),
          m_center(a_center), m_L(a_L), m_dt(a_dt), m_time(a_time)
    {
    }

    //! Destructor
    ~PointExtraction() {}

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator,
                       std::string a_file_prefix) const
    {
        CH_TIME("CustomExtraction::execute_query");
        if (a_interpolator == nullptr)
        {
            MayDay::Error("Interpolator has not been initialised.");
        }
        // std::vector<double> interp_var_data(m_num_points);
        std::vector<double> interp_var1(m_num_points);
        std::vector<double> interp_var2(m_num_points);
        std::vector<double> interp_x(m_num_points);
        std::vector<double> interp_y(m_num_points);
        std::vector<double> interp_z(m_num_points);

        // Work out the coordinates
        // go out radially in diagonal dircetion to L/2
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            interp_x[idx] =
                m_center[0] + (double(idx) / double(m_num_points) * 0.5 * m_L);
            interp_y[idx] =
                m_center[1] + (double(idx) / double(m_num_points) * 0.5 * m_L);
            interp_z[idx] =
                m_center[2] + (double(idx) / double(m_num_points) * 0.5 * m_L);
        }

        // set up the query
        InterpolationQuery query1(m_num_points);
        query1.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data())
            .addComp(m_comp1, interp_var1.data(), Derivative::LOCAL,
                     VariableType::diagnostic); // evolution/diagnostic

        InterpolationQuery query2(m_num_points);
        query2.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data())
            .addComp(m_comp2, interp_var2.data(), Derivative::LOCAL,
                     VariableType::diagnostic); // evolution/diagnostic

        // submit the query
        a_interpolator->interp(query1);
        a_interpolator->interp(query2);

        // now write out
        bool first_step = (m_time == 0.0);
        double restart_time = 0.0;
        SmallDataIO output_file(a_file_prefix, m_dt, m_time, restart_time,
                                SmallDataIO::APPEND, first_step);

        std::vector<std::string> header_line(2);
        if (first_step)
        {
            header_line[0] = "rho_c";
            header_line[1] = "contrast_c";
            output_file.write_header_line(header_line);
        }

        for (int idx = 0; idx < m_num_points; ++idx)
        {
            std::vector<double> data(2);
            data[0] = interp_var1[idx];
            data[1] = interp_var2[idx];
            output_file.write_time_data_line(data);
        }
    }
};

#endif /* POINTEXTRACTION_HPP_ */
