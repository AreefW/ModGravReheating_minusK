/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "KerrBH4dSTLevel.hpp"
#include "BoxLoops.hpp"
#include "ChiTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedGravityConstraints.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "RhoDiagnostics.hpp"
#include "RHSDiagnostics.hpp" // WCC diagnostics
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"
#include "AMRReductions.hpp"

//Tagging
#include "HamTaggingCriterion.hpp"

// write output
#include "SmallDataIO.hpp"
// Custom extraction
#include "CustomExtraction.hpp"
#include "ExcisionDiagnostic.hpp"

// Initial data
#include "GammaCalculator.hpp"
#include "KerrBH.hpp"

void KerrBH4dSTLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

void KerrBH4dSTLevel::initialData()
{
    CH_TIME("KerrBH4dSTLevel::initialData");
    if (m_verbosity)
        pout() << "KerrBH4dSTLevel::initialData " << m_level << endl;

    // First set everything to zero then calculate initial data  Get the Kerr
    // solution in the variables, then calculate the \tilde\Gamma^i numerically
    // as these are non zero and not calculated in the Kerr ICs
    BoxLoops::loop(
        make_compute_pack(SetValue(0.),
                          InitialScalarData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
    
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);

}

// Things to do when restarting from a checkpoint file <-- Set new K as it depends on couplings
void KerrBH4dSTLevel::postRestart()
{
    
    // only want to do this on the first restart and also every restart
    if (m_time == 0.0)
    {   
        // data hierarchy should be set up but fill ghosts just to be sure  
        fillAllGhosts();
        CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
        FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
	    ModifiedGravityConstraints<FourDerivScalarTensorWithCouplingAndPotential>
        constraints(fdst, m_dx, m_p.center, m_p.G_Newton, m_bh_amr.m_rho_mean/*added for rho contrast*/, c_Ham,
                    Interval(c_Mom, c_Mom),c_Ham_abs_sum);
        RhoDiagnostics<FourDerivScalarTensorWithCouplingAndPotential>
        rho_diagnostics(fdst, m_dx, m_p.center);
        auto compute_pack = make_compute_pack(constraints, rho_diagnostics);

        BoxLoops::loop(compute_pack, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);

        pout() << "Setting K mean on restart at t = " << m_time << " on level "
               << m_level << endl;

        AMRReductions<VariableType::diagnostic> amr_reductions_diag(m_bh_amr);
            double phys_vol = amr_reductions_diag.sum(c_sqrt_gam);
            //m_bh_amr.m_rho_mean = amr_reductions_diag.sum(c_rho_phi) / phys_vol;
            //m_bh_amr.m_rho_mean = amr_reductions_diag.sum(c_rho_scaled); // / phys_vol;
            double rho_mean_all = (amr_reductions_diag.sum(c_rho_phi) + amr_reductions_diag.sum(c_rho_g2) + amr_reductions_diag.sum(c_rho_g3) + amr_reductions_diag.sum(c_rho_GB)) / phys_vol;
            m_bh_amr.m_K_mean = - sqrt(3.0 * rho_mean_all);

        // pout() << "Calculated K mean as " << m_bh_amr.m_K_mean
        //        << " at t = " << m_time << " on restart at level " << m_level
        //        << endl;
        // pout() << "PostRestart c_rho_contrast = " << amr_reductions_diag.sum(c_rho_contrast) << endl;
         pout() << "rho phi = " << amr_reductions_diag.sum(c_rho_phi) / phys_vol << endl;
         pout() << "rho g2 = " << amr_reductions_diag.sum(c_rho_g2) / phys_vol << endl;
         pout() << "rho g3 = " << amr_reductions_diag.sum(c_rho_g3) / phys_vol << endl;
         pout() << "rho GB = " << amr_reductions_diag.sum(c_rho_GB) / phys_vol << endl;
         pout() << "rho mean all = " << rho_mean_all << endl;
         pout() << "phys_vol = " << phys_vol << endl;
         pout() << "K mean = " << m_bh_amr.m_K_mean << endl;
    }
}

#ifdef CH_USE_HDF5
void KerrBH4dSTLevel::prePlotLevel()
{
    bool last_step = (m_time == m_p.stop_time);
    /*CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);*/
    fillAllGhosts();
    //---
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);

    //---
    ModifiedGravityConstraints<FourDerivScalarTensorWithCouplingAndPotential>
        constraints(fdst, m_dx, m_p.center, m_p.G_Newton, m_bh_amr.m_rho_mean/*added for rho contrast*/, c_Ham,
                    Interval(c_Mom, c_Mom),c_Ham_abs_sum);
    RhoDiagnostics<FourDerivScalarTensorWithCouplingAndPotential>
        rho_diagnostics(fdst, m_dx, m_p.center);
    auto compute_pack = make_compute_pack(constraints, rho_diagnostics);

    BoxLoops::loop(compute_pack, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
    int min_level = 0;
    bool calculate_diagnostics = at_level_timestep_multiple(min_level);
    if (m_level == min_level && calculate_diagnostics)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions_diag(m_bh_amr);
            //pout() << "PrePlotLevel c_rho_contrast = " << amr_reductions_diag.sum(c_rho_contrast)<< endl;
            //double rho_max = amr_reductions_diag.max(c_rho_scaled);
        double rho_max = amr_reductions_diag.max(c_rho_scaled);
                BoxLoops::loop(
            Excision95Density<FourDerivScalarTensorWithCouplingAndPotential>(m_dx, m_p.center, rho_max),
             m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
             disable_simd());
        /*if (last_step){
            pout() << "prePlot last step called" << endl;
            double rho_max = amr_reductions_diag.max(c_rho_scaled);
                BoxLoops::loop(
            Excision95Density<FourDerivScalarTensorWithCouplingAndPotential>(m_dx, m_p.center, rho_max),
             m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
             disable_simd());
            }
        */
        }

    // BoxLoops::loop(
    //         ExcisionDiagnostics<FourDerivScalarTensorWithCouplingAndPotential>(
    //             m_dx, m_p.center_obj, m_p.obj_r),
    //         m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
    //         disable_simd());
}
#endif /* CH_USE_HDF5 */

void KerrBH4dSTLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                      const double a_time)
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate ModifiedCCZ4 right hand side with theory_t =
    // FourDerivScalarTensor
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
    RhoDiagnostics<FourDerivScalarTensorWithCouplingAndPotential>
        rho_diagnostics(fdst, m_dx, m_p.center);

    BoxLoops::loop(rho_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
    ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);
    if (m_p.max_spatial_derivative_order == 4)
    {
        ModifiedCCZ4RHS<FourDerivScalarTensorWithCouplingAndPotential,
                        ModifiedPunctureGauge, FourthOrderDerivatives>
            my_modified_ccz4(fdst, m_p.modified_ccz4_params,
                             modified_puncture_gauge, m_dx, m_p.sigma,
                             m_bh_amr.m_K_mean, m_bh_amr.m_rho_mean,
                             m_p.center, m_p.G_Newton);
        BoxLoops::loop(my_modified_ccz4, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        ModifiedCCZ4RHS<FourDerivScalarTensorWithCouplingAndPotential,
                        ModifiedPunctureGauge, SixthOrderDerivatives>
            my_modified_ccz4(fdst, m_p.modified_ccz4_params,
                             modified_puncture_gauge, m_dx, m_p.sigma,
                             m_bh_amr.m_K_mean, m_bh_amr.m_rho_mean,
                             m_p.center, m_p.G_Newton);
        BoxLoops::loop(my_modified_ccz4, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

void KerrBH4dSTLevel::specificUpdateODE(GRLevelData &a_soln,
                                        const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void KerrBH4dSTLevel::preTagCells()
{
    // Pre tagging - fill ghost cells and calculate Ham terms
    fillAllEvolutionGhosts();
    CouplingAndPotential coupling_and_potential(m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
	ModifiedGravityConstraints<FourDerivScalarTensorWithCouplingAndPotential>
        constraints(fdst, m_dx, m_p.center, m_p.G_Newton, m_bh_amr.m_rho_mean/*added for rho contrast*/, c_Ham,
                    Interval(c_Mom, c_Mom),c_Ham_abs_sum);
    RhoDiagnostics<FourDerivScalarTensorWithCouplingAndPotential>
        rho_diagnostics(fdst, m_dx, m_p.center);
    auto compute_pack = make_compute_pack(constraints, rho_diagnostics);

    BoxLoops::loop(compute_pack, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
}

void KerrBH4dSTLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                              const FArrayBox &current_state,
                                              const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(HamTaggingCriterion(m_dx, m_p.center_BS, m_p.rad), current_state_diagnostics,
                   tagging_criterion);
}

void KerrBH4dSTLevel::specificPostTimeStep()
{
    int min_level = 0;
    // No need to evaluate the diagnostics more frequently than every coarse
    // timestep, but must happen on every level (not just level zero or data
    // will not be populated on finer levels)
    bool calculate_diagnostics = at_level_timestep_multiple(min_level);

	bool first_step = (m_time == 0.);

    bool last_step = (m_time == m_p.stop_time);
    bool is_same_timestep = at_level_timestep_multiple(min_level);

    if (calculate_diagnostics)
    {
	fillAllGhosts();
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
    ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);
	ModifiedGravityConstraints<FourDerivScalarTensorWithCouplingAndPotential>
        constraints(fdst, m_dx, m_p.center, m_p.G_Newton, m_bh_amr.m_rho_mean/*added for rho contrast*/, c_Ham,
                    Interval(c_Mom, c_Mom),c_Ham_abs_sum);
    RhoDiagnostics<FourDerivScalarTensorWithCouplingAndPotential>
        rho_diagnostics(fdst, m_dx, m_p.center);

    RHSDiagnostics<FourDerivScalarTensorWithCouplingAndPotential,
                   ModifiedPunctureGauge, FourthOrderDerivatives>
        rhs_diagnostics(fdst, m_p.modified_ccz4_params, modified_puncture_gauge,
                        m_dx, m_p.sigma, m_bh_amr.m_K_mean, m_bh_amr.m_rho_mean, m_p.center, m_p.G_Newton);
    
    //
    auto compute_pack = make_compute_pack(constraints, rho_diagnostics, rhs_diagnostics);

    BoxLoops::loop(compute_pack, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
        if (m_level == min_level) //
        {
            AMRReductions<VariableType::diagnostic> amr_reductions_diag(m_bh_amr);
            double phys_vol = amr_reductions_diag.sum(c_sqrt_gam);
            double L2_Ham = amr_reductions_diag.norm(c_Ham);
            double L2_Mom = amr_reductions_diag.norm(c_Mom);
            //double Ham_abs = amr_reductions_diag.sum(c_Ham_abs_sum);
            //double Hamsum = amr_reductions_diag.sum(c_Ham);
	        //double rho_mean = amr_reductions_diag.sum(c_rho_scaled)/phys_vol;
            //pout() << "phyvol = " << phys_vol << endl;
            //pout() << "rho mean  = " << rho_mean << endl;
            //----
            //pout() << "PostTimeStep c_rho_contrast = " << amr_reductions_diag.sum(c_rho_contrast)<< endl;
            double wcc_max = amr_reductions_diag.max(c_weak_coupling_condition_GB);
            m_bh_amr.m_rho_mean = amr_reductions_diag.sum(c_rho_scaled) / phys_vol;
            double rho_mean_all = (amr_reductions_diag.sum(c_rho_phi) + amr_reductions_diag.sum(c_rho_g2) + amr_reductions_diag.sum(c_rho_g3) + amr_reductions_diag.sum(c_rho_GB))/ phys_vol;
            
            double rho_max = amr_reductions_diag.max(c_rho_scaled);
                pout() << "rho_max : " << rho_max << endl;
            BoxLoops::loop(
            Excision95Density<FourDerivScalarTensorWithCouplingAndPotential>(m_dx, m_p.center, rho_max),
             m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
             disable_simd());

            /*
            if (last_step){
                double rho_max = amr_reductions_diag.max(c_rho_scaled);
                pout() << "rho_max : " << rho_max << endl;
                BoxLoops::loop(
            Excision95Density<FourDerivScalarTensorWithCouplingAndPotential>(m_dx, m_p.center, rho_max),
             m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
             disable_simd());
            }
            */
        
           //m_bh_amr.m_S_mean = amr_reductions_diag.sum(c_S_scaled) / phys_vol;
            m_bh_amr.m_K_mean = - sqrt(3.0 * m_bh_amr.m_rho_mean);

            //---- calculate volume and mass of object
            double vol_obj = amr_reductions_diag.sum(c_sqrt_gam_exc);
            //double mass_obj = amr_reductions_diag.sum(c_rho_exc) * vol_obj;
            double mass_obj = amr_reductions_diag.sum(c_rho_exc);
            
            //----
	        AMRReductions<VariableType::evolution> amr_reductions_evo(m_bh_amr);

            double chi_mean = amr_reductions_evo.sum(c_chi) / phys_vol ;
            double lapse = amr_reductions_evo.sum(c_lapse) / phys_vol ;

            //----
          // BoxLoops::loop(SetValue(m_bh_amr.m_K_mean, Interval(c_K, c_K)),
          //             m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
            //----
            SmallDataIO constraints_file(m_p.data_path + "data_out",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"<chi>", "<rho>", "L2_Ham", "L2_Mom", "Vol_osc", "M_osc", "rho_max", "wcc max"});
            }
            constraints_file.write_time_data_line({chi_mean, m_bh_amr.m_rho_mean, L2_Ham, L2_Mom, vol_obj, mass_obj, rho_max, wcc_max});

        //Custom Extaction
            // set up an interpolator
            // pass the boundary params so that we can use symmetries if
            // applicable
        //    AMRInterpolator<Lagrange<4>> interpolator(
        //        m_bh_amr, m_p.origin, m_p.dx, m_p.boundary_params,
        //        m_p.verbosity);

            // this should fill all ghosts including the boundary ones according
            // to the conditions set in params.txt
        //    interpolator.refresh();

            // set up the query and execute it
        //    int num_points = 1;
        //    std::array<double, CH_SPACEDIM> extr_point = {33.25, 18.25, 33.25}; // specified point
        //    CustomExtraction extraction(c_rho_scaled, num_points, m_p.L, extr_point,
        //                               m_dt, m_time);
        //    extraction.execute_query(&interpolator, m_p.data_path + "rho_dat");
        //----
        }
    }
}
