// Copyright (C) 2004--2017 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cmath>
#include <cassert>
#include <gsl/gsl_math.h>

#include "MohrCoulombYieldStressEvan.hh"
#include "MohrCoulombYieldStress.hh"

#include "pism/hydrology/Hydrology.hh"
#include "pism/hydrology/hydrologyEvan.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {

//! \file PISMMohrCoulombYieldStressEvan.cc  Process model which computes pseudo-plastic yield stress for the subglacial layer.
/*! \file PISMMohrCoulombYieldStressEvan.cc
The output variable of this submodel is `tauc`, the pseudo-plastic yield stress
field that is used in the ShallowStressBalance objects.  This quantity is
computed by the Mohr-Coulomb criterion [\ref SchoofTill], but using an empirical
relation between the amount of water in the till and the effective pressure
of the overlying glacier resting on the till [\ref Tulaczyketal2000].

The "dry" strength of the till is a state variable which is private to
the submodel, namely `tillphi`.  Its initialization is nontrivial: either the
`-topg_to_phi`  heuristic is used or inverse modeling can be used.  (In the
latter case `tillphi` can be read-in at the beginning of the run.  Currently
`tillphi` does not evolve during the run.)

This submodel uses a pointer to a Hydrology instance to get the till (pore)
water amount, the effective water layer thickness.  The effective pressure is
derived from this.  Then the effective pressure is combined with tillphi to
compute an updated `tauc` by the Mohr-Coulomb criterion.

This submodel is inactive in floating areas.
*/


MohrCoulombYieldStressEvan::MohrCoulombYieldStressEvan(IceGrid::ConstPtr g,
                                               hydrology::Hydrology *hydro)
  : MohrCoulombYieldStress(g, hydro), m_hydrology(hydro) {


//  m_log->message(2,
//             "* Starting MohrCoulombYieldStressEvan ...\n");


  m_effective_pressure.create(m_grid, "effective_pressure", WITHOUT_GHOSTS);
  m_effective_pressure.set_attrs("internal",
                       "fraction of overburden pressure of subglacial hydrology system",
                       "1", "");

  m_sliding_mechanism.create(m_grid, "sliding_mechanism", WITHOUT_GHOSTS);
  m_sliding_mechanism.set_attrs("internal",
                       "0 for no sliding (i.e. ice free), 1 for sediment deformation, 2 for hydrological influence flow",
                       "1", "");

  m_till_cover_local.create(m_grid, "tillcover_local",
              WITHOUT_GHOSTS);
  m_till_cover_local.set_attrs("internal",
                 "copy of till cover from the hydrology model",
                 "1", "");



  m_velocity_temp.create(m_grid, "velocity_temp",
              WITHOUT_GHOSTS);
  m_velocity_temp.set_attrs("internal",
                 "copy of velocity from the hydrology model",
                 "m s-1", "");



  hydro_tauc.create(m_grid, "hydro_tauc",
              WITHOUT_GHOSTS);
  hydro_tauc.set_attrs("internal",
                 "tauc calculated from the hydrology model",
                 "Pa", "");


}

MohrCoulombYieldStressEvan::~MohrCoulombYieldStressEvan() {
  // empty
}


//! Initialize the pseudo-plastic till mechanical model.
/*!

*/
void MohrCoulombYieldStressEvan::init_impl() {

//  m_log->message(2,
//             "* Starting MohrCoulombYieldStressEvan::init_impl ...\n");

  MohrCoulombYieldStress::init_impl();

  m_effective_pressure.set(0.0);

}

void MohrCoulombYieldStressEvan::define_model_state_impl(const PIO &output) const {
  MohrCoulombYieldStress::define_model_state_impl(output);

  m_basal_yield_stress.define(output);
  m_sliding_mechanism.define(output);
  hydro_tauc.define(output);
}

void MohrCoulombYieldStressEvan::write_model_state_impl(const PIO &output) const {
  MohrCoulombYieldStress::write_model_state_impl(output);
  m_basal_yield_stress.write(output);
  m_sliding_mechanism.write(output);
  hydro_tauc.write(output);

}

//! Update the till yield stress for use in the pseudo-plastic till basal stress
//! model.  See also IceBasalResistancePlasticLaw.
/*!

 */
void MohrCoulombYieldStressEvan::update_impl(const YieldStressInputs &inputs) {

  m_log->message(2,
             "* Starting MohrCoulombYieldStressEvan::update_impl ...\n");
//
  bool slippery_grounding_lines = m_config->get_boolean("basal_yield_stress.slippery_grounding_lines");

  hydrology::hydrologyEvan* hydroEvan = dynamic_cast<hydrology::hydrologyEvan*>(m_hydrology);

  if ( not hydroEvan) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "MohrCoulombYieldStressEvan requires the Hydrology model hydrologyEvan");
  }

  const double high_tauc   = m_config->get_double("basal_yield_stress.ice_free_bedrock"),
               tillwat_max = m_config->get_double("hydrology.tillwat_max"),
               c0          = m_config->get_double("basal_yield_stress.mohr_coulomb.till_cohesion"),
               N0          = m_config->get_double("basal_yield_stress.mohr_coulomb.till_reference_effective_pressure"),
               e0overCc    = (m_config->get_double("basal_yield_stress.mohr_coulomb.till_reference_void_ratio") /
                              m_config->get_double("basal_yield_stress.mohr_coulomb.till_compressibility_coefficient")),
               delta       = m_config->get_double("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden"),
               tlftw       = m_config->get_double("basal_yield_stress.mohr_coulomb.till_log_factor_transportable_water"),
               max_effective_pressure_ratio = m_config->get_double("hydrology.maximum_effective_pressure_ratio"),
               K1          = m_config->get_double("basal_yield_stress.mohr_coulomb_evan.sliding_flow_factor"),
               K2          = m_config->get_double("basal_yield_stress.mohr_coulomb_evan.deformation_flow_factor"),
               rho_i = m_config->get_double("constants.ice.density"),
               g = m_config->get_double("constants.standard_gravity"),
               q = m_config->get_double("basal_resistance.pseudo_plastic.q");
        double m_pseudo_u_threshold = m_config->get_double("basal_resistance.pseudo_plastic.u_threshold", "m second-1");



  const IceModelVec2CellType &mask           = inputs.geometry->cell_type;
  const IceModelVec2S        &bed_topography = inputs.geometry->bed_elevation;
  const IceModelVec2S        &temp_thk = *m_grid->variables().get_2d_scalar("thk");

  IceModelVec::AccessList list{&m_tillwat, &m_till_phi, &m_basal_yield_stress, &mask,
      &bed_topography, &m_Po, &m_till_cover_local, &m_effective_pressure, &m_sliding_mechanism, &m_velocity_temp, &hydro_tauc, &temp_thk};

  if (hydroEvan) {
    hydroEvan->till_water_thickness(m_tillwat);
    hydroEvan->overburden_pressure(m_Po);
    hydroEvan->get_EffectivePressure(m_effective_pressure);
    hydroEvan->get_SedimentDistribution(m_till_cover_local);
    hydroEvan->update_velbase_mag(m_velocity_temp);

  }

  m_log->message(2,
             "* calculating basal strength for each cell ...\n");
  double seconds_in_year = 365.0*24.0*3600.0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

//  m_log->message(2,
//             "* %i %i %f %i %i %f\n", i, j, m_basal_yield_stress(i, j), mask.ocean(i, j), mask.ice_free(i, j), m_tillwat(i,j));

    if (mask.ocean(i, j)) {

      m_basal_yield_stress(i, j) = 0.0;

        m_sliding_mechanism(i,j) = 0;
    } else if (mask.ice_free(i, j)) {
      m_basal_yield_stress(i, j) = high_tauc;  // large yield stress if grounded and ice-free
      m_sliding_mechanism(i,j) = 0;
    } else { // grounded and there is some ice
      // user can ask that marine grounding lines get special treatment
      const double sea_level = 0.0; // FIXME: get sea-level from correct PISM source
      double water = m_tillwat(i,j); // usual case
      if (slippery_grounding_lines and
          bed_topography(i,j) <= sea_level and
          (mask.next_to_floating_ice(i,j) or mask.next_to_ice_free_ocean(i,j))) {
        water = tillwat_max;
      }
      double
        s    = water / tillwat_max,
        Ntil = N0 * pow(delta * m_Po(i,j) / N0, s) * pow(10.0, e0overCc * (1.0 - s));
        Ntil = std::min(m_Po(i,j), Ntil);

        m_basal_yield_stress(i, j) = c0 + Ntil * tan((M_PI/180.0) * m_till_phi(i, j));


//  m_log->message(2, 
//             "* %i %i %f ...\n", i, j, m_velocity_temp(i, j));

       // take into account sediment free areas

       m_basal_yield_stress(i, j) = m_basal_yield_stress(i, j) * m_till_cover_local(i, j) + temp_thk(i,j) * g * rho_i * (1.0 - m_till_cover_local(i, j));

       m_sliding_mechanism(i,j) = 1;

       if (hydroEvan) {

         // sliding multipliers used by Arnold and Sharp (2002)
 //        double K1 = 6.3376e-6;
 //        double K2 = 400.0;

          double K1_override = 5e-8;
          double K2_override = 400;
          double yield_stress_hydrology2;
          double yield_stress_hydrology3;
          double yield_stress_hydrology4;
         // based on equation A.4 in Arnold and Sharp (2002)
         // Note, the equation in Arnold and Sharp is messed up, I ended up converting it to the orignal form used in the paper by McInnes and Budd (1984), which
         // used value of "ice thickness above buoyancy". The values of K1 and K2 were designed for this, and the Arnold and Sharp equation used a value of K2 which
         // was not converted to the correct units.

         double z_star = m_effective_pressure(i,j) / (rho_i * g); //ice_thickness_above_buoyancy
         double yield_stress_hydrology;

//  m_log->message(2, 
//             "* %i %i %f %f ...\n", i, j, Ntil, m_effective_pressure(i,j));

         if (m_velocity_temp(i,j) > 0.0) {



 //          yield_stress_hydrology = (z_star+pow(z_star,2)/K2_override) / K1_override   * (pow(m_pseudo_u_threshold,q) * pow(m_velocity_temp(i,j),1.0-q));
           yield_stress_hydrology = m_effective_pressure(i,j);


         } else {
           yield_stress_hydrology = high_tauc; // to prevent errors
//		yield_stress_hydrology2 = high_tauc;
//            yield_stress_hydrology3 = high_tauc;
//            yield_stress_hydrology4 = high_tauc;
         }

         hydro_tauc(i,j) = yield_stress_hydrology;
/*
  m_log->message(2, 
             "*  %i %i %f %f  ...\n", i, j, z_star, (z_star) / K1);
 
//  m_log->message(2, 
//             "* > %i %i %f %f %f %f %f ...\n", i, j, m_effective_pressure(i,j), yield_stress_hydrology, yield_stress_hydrology2, yield_stress_hydrology3, yield_stress_hydrology4);
          
  m_log->message(2, 
             "* >> %i %i %f %f %f %f ...\n", i, j, (z_star+pow(z_star,2)/K2) / K1 * seconds_in_year,  m_pseudo_u_threshold,  m_velocity_temp(i, j)*seconds_in_year, pow(m_pseudo_u_threshold,q) * pow(m_velocity_temp(i,j),1.0-q));
  m_log->message(2,
             "* >>> %i %i %f %f %i\n", i, j, m_basal_yield_stress(i, j), yield_stress_hydrology, yield_stress_hydrology < m_basal_yield_stress(i, j));
*/
         // if the ice base is weaker than the sediments
         if (yield_stress_hydrology < m_basal_yield_stress(i, j)) {

           m_basal_yield_stress(i, j) = yield_stress_hydrology;
           m_sliding_mechanism(i,j) = 2;
         }



       }

/*
       // find ratio of effective pressure due to hydrology and overburden, limited by the user defined maximum ratio

         double pressure_ratio;
         if (m_Po(i, j) > 0.0) {
           pressure_ratio = m_effective_pressure(i, j) / m_Po(i, j);
         } else {
           pressure_ratio = 0.0;
         }

       pressure_ratio = std::min(pressure_ratio,max_effective_pressure_ratio);

       // reduce the effective basal yield stress by this ratio

       double k1 = 6.3e-8;
       double k2 = 400.0;

       double test_var = k1 * m_basal_yield_stress(i, j) / ( (1.0 - pressure_ratio) * m_Po(i, j) );

       m_basal_yield_stress(i, j) = m_basal_yield_stress(i, j) * (1.0 - pressure_ratio);

//  m_log->message(2,
//             "* %i %i %f %f %f ...\n", i, j, m_basal_yield_stress(i, j),test_var, pressure_ratio);

*/

    }


//  m_log->message(2,
//             "* %i %i %f ...\n", i, j, m_basal_yield_stress(i, j));
  }


  m_log->message(2,
             "* Finished MohrCoulombYieldStressEvan::update_impl ...\n");

  m_basal_yield_stress.update_ghosts();
}

} // end of namespace pism
