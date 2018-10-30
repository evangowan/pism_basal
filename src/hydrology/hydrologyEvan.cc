// Copyright (C) 2012-2015 PISM Authors
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

#include <iostream>
#include "Hydrology.hh"
#include "hydrologyEvan.hh"
#include "pism/util/Mask.hh"
#include "pism/util/error_handling.hh"
#include "hydrology_diagnostics.hh"
#include "pism/util/Vars.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_options.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/coupler/surface/TemperatureIndex.hh" // needed to grab srunoff variable
#include "pism/stressbalance/StressBalance.hh"



#include "pism/util/iceModelVec2T.hh"

#include "pism/util/io/PIO.hh"

#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"


namespace pism {
namespace hydrology {

hydrologyEvan :: hydrologyEvan(IceGrid::ConstPtr g, stressbalance::StressBalance *sb, surface::SurfaceModel *m_surface)
	: NullTransport(g)
   {

  m_log->message(2,
             "* Starting hydrologyEvan ...\n");

  m_stressbalance = sb;

  m_surfaceT = m_surface;

  unsigned int stencil_width = 1;

  m_melt_rate_local.create(m_grid, "melt_rate_local", WITHOUT_GHOSTS);
  m_melt_rate_local.set_attrs("internal",
                        "hydrology model workspace for melt_rate",
                        "m s-1", "");

  m_pressure_temp.create(m_grid, "pressure_temp", WITH_GHOSTS, stencil_width);
  m_pressure_temp.set_attrs("internal",
                        "pressure at the base",
                        "Pa", "");



  m_surface_gradient_temp.create(m_grid, "surface_gradient_temp", WITHOUT_GHOSTS);
  m_surface_gradient_temp.set_attrs("internal",
                        "temporary storage for the surface gradient",
                        "1", "");

  m_bed_gradient_temp.create(m_grid, "bed_gradient_temp", WITHOUT_GHOSTS);
  m_bed_gradient_temp.set_attrs("internal",
                        "temporary storage for the bed gradient",
                        "1", "");

  m_gradient_temp.create(m_grid, "gradient_temp", WITHOUT_GHOSTS);
  m_gradient_temp.set_attrs("internal",
                        "gradient at the base",
                        "Pa m-1", "");



  m_hydro_gradient.create(m_grid, "hydro_gradient", WITHOUT_GHOSTS);
  m_hydro_gradient.set_attrs("internal",
                        "pressure gradient at the base",
                        "Pa m-1", "");

  m_hydro_gradient_dir.create(m_grid, "hydro_gradient_dir", WITHOUT_GHOSTS);
  m_hydro_gradient_dir.set_attrs("internal",
                        "potential gradient at the base in directional components",
                        "Pa m-1", "");

  m_surface_gradient.create(m_grid, "hydro_surface_gradient", WITHOUT_GHOSTS);
  m_surface_gradient.set_attrs("internal",
                        "surface gradient in hydrology",
                        "1", "");

  m_surface_gradient_dir.create(m_grid, "hydro_surface_gradient_dir", WITHOUT_GHOSTS);
  m_surface_gradient_dir.set_attrs("internal",
                        "surface gradient in hydrology directional components",
                        "1", "");

  m_surface_elevation_temp.create(m_grid, "surface_elevation_temp", WITH_GHOSTS, stencil_width);
  m_surface_elevation_temp.set_attrs("model_state",
                       "temporary surface_elevation",
                       "m", "");

  m_bed_elevation_temp.create(m_grid, "bed_elevation_temp", WITH_GHOSTS, stencil_width);
  m_bed_elevation_temp.set_attrs("model_state",
                       "temporary bed_elevation",
                       "m", "");




  m_till_cover.create(m_grid, "tillcover", WITHOUT_GHOSTS);
  m_till_cover.set_attrs("model_state",
                       "fraction of surface that is covered in till",
                       "1", "");

  m_gradient_permutation.create(m_grid, "gradient_permutation", WITHOUT_GHOSTS);
  m_gradient_permutation.set_attrs("model_state",
                       "permutation array for sorting the gradient",
                       "1", "");


  m_total_input_ghosts.create(m_grid, "total_input", WITH_GHOSTS, stencil_width); // need ghosts here
  m_total_input_ghosts.set_attrs("internal",
                          "hydrology model workspace for total input rate into subglacial water layer with ghosts",
                          "m s-1", "");


  m_volume_water_flux.create(m_grid, "volume_water_flux", WITHOUT_GHOSTS);
  m_volume_water_flux.set_attrs("model_state",
                       "volume water flux through tunnel system",
                       "m3 s-1", "");

  m_tunnel_cross_section.create(m_grid, "tunnel_cross_section", WITHOUT_GHOSTS);
  m_tunnel_cross_section.set_attrs("internal",
                        "cross section area of a tunnel",
                        "m-2", "");

  m_velbase_mag.create(m_grid, "velbase_mag", WITHOUT_GHOSTS);
  m_velbase_mag.set_attrs("internal",
                        "ice sliding speed seen by subglacial hydrology",
                        "m s-1", "");
  m_velbase_mag.metadata().set_double("valid_min", 0.0);


  m_hydrology_effective_pressure.create(m_grid, "hydrology_effective_pressure", WITHOUT_GHOSTS);
  m_hydrology_effective_pressure.set_attrs("model_state",
                        "effective pressure due to hydrology",
                        "Pa", "");
  m_hydrology_effective_pressure.metadata().set_double("valid_min", 0.0);

  m_hydrology_fraction_overburden.create(m_grid, "fraction_overburden", WITHOUT_GHOSTS);
  m_hydrology_fraction_overburden.set_attrs("internal",
                        "fraction of effective pressure to overburden",
                        "1", "");
  m_hydrology_effective_pressure.metadata().set_double("valid_min", 0.0);



 m_hydrosystem.create(m_grid, "hydrology_type", WITHOUT_GHOSTS);
 m_hydrosystem.set_attrs("model_state",
                        "type of hydrology, 1 tunnels 2 cavity",
                        "1", "");
 m_hydrosystem.metadata().set_double("valid_min", 0.0);


 m_check_mask.create(m_grid, "check_mask", WITHOUT_GHOSTS);
 m_check_mask.set_attrs("internal",
                        "check for if a cell has been checked in the gradient sort",
                        "1", "");
 m_check_mask.metadata().set_double("valid_min", 0.0);


}

hydrologyEvan::~hydrologyEvan() {
}


void hydrologyEvan::init() {
  m_log->message(2,
             "* Initializing Evan's null-transport subglacial hydrology model ...\n");
  Hydrology::init();

  options::Real
    till_cover("-till_fraction_coverage", "fill value for fraction of surface covered by sediments",
                m_config->get_double("hydrology.till_fraction_coverage"));


  m_till_cover.set(till_cover);
  m_volume_water_flux.set(0.0);
  m_hydrosystem.set(0.0);

  // initialize the permutation array, after the first sorting it should not take long to sort



  int counter = 0;

  IceModelVec::AccessList list;
  list.add(m_gradient_permutation);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_log->message(2,
             "%i %i %i\n", i, j, counter);

    m_gradient_permutation(i, j) = counter;
    counter++;

  }

  m_log->message(2,
             "* Finished initializing the permutation array ...\n");

}


void hydrologyEvan::define_model_state_impl(const PIO &output) const {

  m_till_cover.define(output);
  m_hydrology_effective_pressure.define(output);
  m_total_input_ghosts.define(output);
  m_volume_water_flux.define(output);
  m_melt_rate_local.define(output);
  m_hydrosystem.define(output);
  m_velbase_mag.define(output);
  m_hydro_gradient.define(output);
  m_hydro_gradient_dir.define(output);
  m_tunnel_cross_section.define(output);
  m_hydrology_fraction_overburden.define(output);
  m_gradient_permutation.define(output);
}

void hydrologyEvan::write_model_state_impl(const PIO &output) const {
  m_Wtil.write(output);
  m_till_cover.write(output);
  m_hydrology_effective_pressure.write(output);
  m_total_input_ghosts.write(output);
  m_volume_water_flux.write(output);
  m_melt_rate_local.write(output);
  m_hydrosystem.write(output);
  m_velbase_mag.write(output);
  m_hydro_gradient.write(output);
  m_hydro_gradient_dir.write(output);
  m_tunnel_cross_section.write(output);
  m_hydrology_fraction_overburden.write(output);
  m_gradient_permutation.write(output);

}




void hydrologyEvan::get_SedimentDistribution(IceModelVec2S &result) {

  result.copy_from(m_till_cover);

}


void hydrologyEvan::get_EffectivePressure(IceModelVec2S &result) {

  result.copy_from(m_hydrology_effective_pressure);

}


// copied from Distributed hydrology
//! Update the the sliding speed |v_b| from ice quantities.
/*!
Calls a StressBalance method to get the vector basal velocity of the ice,
and then computes the magnitude of that.
 */
void hydrologyEvan::update_velbase_mag(IceModelVec2S &result) {
  // velbase_mag = |v_b|
  result.set_to_magnitude(m_stressbalance->advective_velocity());
}

void hydrologyEvan::update_surface_runoff(IceModelVec2S &result) {

  surface::TemperatureIndex *temp_index = dynamic_cast<surface::TemperatureIndex*>(m_surfaceT);
 //   hydrology::Routing *hydrology_routing = dynamic_cast<hydrology::Routing*>(m_hydrology);

  temp_index -> get_runoff_rate(m_melt_rate_local);
}


// potential gradient at the base of the ice sheet
void hydrologyEvan::potential_gradient(IceModelVec2V &result) {

  double rho_i       = m_config->get_double("constants.ice.density");
  double rho_w       = m_config->get_double("constants.fresh_water.density");
  double g           = m_config->get_double("constants.standard_gravity");
  double flotation_fraction = m_config->get_double("hydrology.floatation_fraction");

  double rho_i_g = rho_i * g;
  double density_ratio = rho_w / rho_i;

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  double u, v;


  IceModelVec::AccessList list;
  list.add(m_pressure_temp);
  list.add(m_surface_gradient_temp);
  list.add(m_bed_gradient_temp);
  list.add(m_gradient_temp);


  // grab the bed and surface gradients
  surface_gradient(m_surface_gradient_temp);
  bed_gradient(m_bed_gradient_temp);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Equation 6.13 in Cuffy and Paterson (2010)
    // The flotation fraction is the ratio of the pressure of water to the pressure of ice, and will influence the effect of the bed gradient on the total potential gradient
    // At a default, it is set to 0.8, which gives means the bed slope needs to be 2.7 times greater than the surface slope to have an equal influence on the direction of water flow.
    m_gradient_temp(i,j).u = rho_i_g * (flotation_fraction * m_surface_gradient_temp(i,j).u + (density_ratio - flotation_fraction) * m_bed_gradient_temp(i,j).u);
    m_gradient_temp(i,j).v = rho_i_g * (flotation_fraction * m_surface_gradient_temp(i,j).v + (density_ratio - flotation_fraction) * m_bed_gradient_temp(i,j).v);


  }

  result.copy_from(m_gradient_temp);

//  m_log->message(2,
//             "* Finished hydrologyEvan::potential_gradient ...\n");

}

// find the gradient of the surface
void hydrologyEvan::surface_gradient(IceModelVec2V &result) {


  const IceModelVec2S &surface_elevation = *m_grid->variables().get_2d_scalar("surface_altitude");

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  double u, v;


  IceModelVec::AccessList list;
  list.add(surface_elevation);
  list.add(m_surface_elevation_temp);

  m_surface_elevation_temp.copy_from(surface_elevation);
  m_surface_elevation_temp.update_ghosts();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

   
    // third order finite difference method for calculationg gradient, equations 3 and 4 in Skidmore (1989)
    result(i,j).u = ((m_surface_elevation_temp(i+1,j+1) + 2.0 * m_surface_elevation_temp(i+1,j) + m_surface_elevation_temp(i+1,j-1)) -
			    (m_surface_elevation_temp(i-1,j+1) + 2.0 * m_surface_elevation_temp(i-1,j) + m_surface_elevation_temp(i-1,j-1))) / (8.0 * dx);

    result(i,j).v = ((m_surface_elevation_temp(i+1,j+1) + 2.0 * m_surface_elevation_temp(i,j+1) + m_surface_elevation_temp(i-1,j+1)) -
			    (m_surface_elevation_temp(i+1,j-1) + 2.0 * m_surface_elevation_temp(i,j-1) + m_surface_elevation_temp(i-1,j-1))) / (8.0 * dy);

  }

}

// find the gradient of the bed
void hydrologyEvan::bed_gradient(IceModelVec2V &result) {

  const IceModelVec2S &bed_elevation = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  double u, v;

  IceModelVec::AccessList list;
  list.add(bed_elevation);
  list.add(m_bed_elevation_temp);

  m_bed_elevation_temp.copy_from(bed_elevation);
  m_bed_elevation_temp.update_ghosts();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // third order finite difference method for calculationg gradient, equations 3 and 4 in Skidmore (1989)
    result(i,j).u = ((m_bed_elevation_temp(i+1,j+1) + 2.0 * m_bed_elevation_temp(i+1,j) + m_bed_elevation_temp(i+1,j-1)) -
			    (m_bed_elevation_temp(i-1,j+1) + 2.0 * m_bed_elevation_temp(i-1,j) + m_bed_elevation_temp(i-1,j-1))) / (8.0 * dx);

    result(i,j).v = ((m_bed_elevation_temp(i+1,j+1) + 2.0 * m_bed_elevation_temp(i,j+1) + m_bed_elevation_temp(i-1,j+1)) -
			    (m_bed_elevation_temp(i+1,j-1) + 2.0 * m_bed_elevation_temp(i,j-1) + m_bed_elevation_temp(i-1,j-1))) / (8.0 * dy);

  }
}


//! Compute the total water input rate into the basal hydrology layer in the ice-covered region, allowing time-varying input from a file.
/*!
The user can specify the total of en- and supra-glacial drainage contributions
to subglacial hydrology in a time-dependent input file using option -hydrology_input_to_bed.
This method includes that possible input along with `bmelt` to get the total water
input into the subglacial hydrology.

This method crops the input rate to the ice-covered region.  It
also uses hydrology_const_bmelt if that is requested.

Call this method using the current \e hydrology time step.  This method
may be called many times per IceModel time step.  See update() method
in derived classes of Hydrology.


Modified to include a parameter for supraglacial melt
 */
void hydrologyEvan::get_input_rate(double hydro_t, double hydro_dt,
                               IceModelVec2S &result) {



//  m_log->message(2,
//             "* Starting hydrologyEvan::get_input_rate %f %f\n", hydro_t, hydro_dt);


  bool   use_const   = m_config->get_boolean("hydrology.use_const_bmelt");
  double const_bmelt = m_config->get_double("hydrology.const_bmelt");
  double const_water_from_surface = m_config->get_double("hydrology.fraction_from_surface"); // According to Caroline Clason, usually 70-95% of the surface melt goes to the base
  const double ice_density = m_config->get_double("constants.ice.density");



  const IceModelVec2S        &bmelt = *m_grid->variables().get_2d_scalar("bmelt");
  const IceModelVec2CellType &mask  = *m_grid->variables().get_2d_cell_type("mask");

  const IceModelVec2S &surface_elevation = *m_grid->variables().get_2d_scalar("surface_altitude");


  if (not m_hold_bmelt) {
    m_bmelt_local.copy_from(bmelt);
  }


  IceModelVec::AccessList list{&m_bmelt_local, &mask, &result, &surface_elevation};
// Need to grab the meltrate from PSTemperature Index
  list.add(m_melt_rate_local); 


  if(m_surfaceT) {
    update_surface_runoff(m_melt_rate_local);
  } else {

    m_melt_rate_local.set(0.0);

  }
  // cheat

 //   m_melt_rate_local.set(5.0);

  double seconds_in_year = 365.0*24.0*3600.0;
  
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

/*
  // cheat
   result(i,j) = 0.0;
   const double
    x = m_grid->x(i),
    y = m_grid->y(j); // hopefully the grid is always square

  // cheat
    double melt_factor = -log10(surface_elevation(i,j)) / 3.0;


    if(melt_factor > 1.0) {
      melt_factor = 1.0;
    }

    if(y < -50000.0 && x >= -50000 && x <= 50000.0 && hydro_t > 10000.0 * seconds_in_year) {
//     m_melt_rate_local(i,j) = 5.0 / seconds_in_year;
      m_melt_rate_local(i,j) =  100.*pow(10.0,melt_factor) / seconds_in_year;

    }

    if(y > 50000.0 && x >= -50000 && x <= 50000.0) {
 //    m_melt_rate_local(i,j) = 0.5 / seconds_in_year;

//      m_melt_rate_local(i,j) = pow(10.0,melt_factor) / seconds_in_year * 0.25;

    }
    // end cheat

*/

    if (mask.icy(i, j)) {

//      if( m_bmelt_local(i,j) < 0.1 / seconds_in_year) {
        result(i,j) = (use_const) ? const_bmelt : m_bmelt_local(i,j); // get the melt water from the base
//      } else {

 //       result(i,j) = (use_const) ? const_bmelt : 0.1 / seconds_in_year; // get the melt water from the base
 //     }

      result(i,j) += m_melt_rate_local(i,j) *  const_water_from_surface; // add on the meltwater from the surface




    } else {
      result(i,j) = 0.0;
    }



 // m_log->message(2,
 //            "* %i %i %e %e %f %f %f\n", i, j, m_melt_rate_local(i,j), result(i,j), surface_elevation(i,j), melt_factor, m_bmelt_local(i,j)*seconds_in_year);
  }

}


//! Determine the hydrology at the base
/*!

This is a non-conservative hydrology model, based on a combination of the sediment deformation
model that is default in PISM, and the tunnel/cavity hydrology model used by Arnold and Sharp (2002).
The model first calculates if the water input is sufficient to fill up the sediments in the grid
cell, and then if it does, determines whether there is enough water input to have a 
tunnel or cavity based drainage system. This takes into account the amount of the base
is covered by sediments (i.e. in places like the Canadian Shield, there is rock sticking out
that does not have water storage).  In this way, a region regarded as mostly rock is more likely
to develop a tunnel based system because there should be more water available to develop it.
 */
void hydrologyEvan::update_impl(double icet, double icedt) {



  m_log->message(2,
             "* Starting hydrologyEvan::update_impl ...\n");

  // if asked for the identical time interval as last time, then do nothing
//  if ((fabs(icet - m_t) < 1e-6) && (fabs(icedt - m_dt) < 1e-6)) {
//    return;
// }
  const double  m_t = icet / 365.0 / 24.0 / 3600.0;
//  m_dt = icedt;

  const double tillwat_max = m_config->get_double("hydrology.tillwat_max"),
               tillwat_decay_rate           = m_config->get_double("hydrology.tillwat_decay_rate"), // this C factor should probably be made to be spatially variable too, enhanced if sand, for now constant
               rho_w       = m_config->get_double("constants.fresh_water.density"),
               g           = m_config->get_double("constants.standard_gravity"),
		   tunnel_spacing = m_config->get_double("hydrology.tunnel_spacing"),
               channel_flow_constant = m_config->get_double("hydrology.channel_flow_constant"),
               bedrock_wavelength = m_config->get_double("hydrology.bedrock_wavelength"),
               bump_amplitude = m_config->get_double("hydrology.bedrock_bump_amplitude"),
               Glen_exponent = m_config->get_double("stress_balance.ssa.Glen_exponent"),
               cavity_area = m_config->get_double("hydrology.cavity_area"),
               cavity_spacing = m_config->get_double("hydrology.cavity_spacing"),
               mu = m_config->get_double("hydrology.power_bedrock"),
               arrhenius_parameter = m_config->get_double("flow_law.isothermal_Glen.ice_softness"),
               rho_i       = m_config->get_double("constants.ice.density"),
               latent_heat = m_config->get_double("constants.fresh_water.latent_heat_of_fusion"),
               shadowing_function = m_config->get_double("hydrology.shadowing_function"),
               max_effective_pressure_ratio = m_config->get_double("hydrology.maximum_effective_pressure_ratio")
  ;



  m_log->message(2,
             "* time: %f\n", m_t);

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy(); // hopefully the grid is always square

//  m_log->message(2,
//             "* dx: %f\n", dx);

  // grab the input rate, which will be the sum of the meltwater generated from basal heating, and meltwater from the surface
  get_input_rate(icet,icedt, m_total_input_ghosts); 


  if (tillwat_max < 0.0) {

 //       throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Hydrology::hydrologyEvan: hydrology_tillwat_max is negative\n"
 //                                     "This is not allowed.");
  }




  const IceModelVec2CellType &mask  = *m_grid->variables().get_2d_cell_type("mask");


  IceModelVec::AccessList list{&mask, &m_velbase_mag};
  list.add(m_Wtil);
  list.add(m_total_input_ghosts);
  list.add(m_till_cover);
  list.add(m_hydro_gradient);
  list.add(m_hydro_gradient_dir);
  list.add(m_surface_gradient);
  list.add(m_surface_gradient_dir);
  list.add(m_tunnel_cross_section);
  list.add(m_hydrology_effective_pressure);
  list.add(m_volume_water_flux);
  list.add(m_pressure_temp);
  list.add(m_hydrosystem);
  list.add(m_hydrology_fraction_overburden);
  list.add(m_gradient_permutation);
  list.add(m_check_mask);



  // first we need to find out if the water input is sufficient to fill up the till in the cell

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();


    if ( mask.ice_free(i,j) && ! mask.ocean(i,j)) {
      m_Wtil(i,j) = 0.0;

    } else if (mask.ocean(i,j)) {

	m_Wtil(i,j) = tillwat_max; // have to assume that sediments under water are saturated. important for when the ice advances, so it doesn't have to fill up.

    } else {

      // calculate the amount of water in the till
      double Wtill_before = m_Wtil(i,j);
      double Wtill_after;
      double before_wat=m_total_input_ghosts(i,j);
      if (m_till_cover(i,j) < 0.01) { // no till
        Wtill_after = Wtill_before;
      } else {
        Wtill_after = m_Wtil(i,j) + icedt * ( m_total_input_ghosts(i,j) - tillwat_decay_rate) / m_till_cover(i,j);
      }

      m_Wtil(i,j) = std::min(std::max(0.0, Wtill_after), tillwat_max);

      if(m_Wtil(i,j) < tillwat_max && m_till_cover(i,j) > 0.01) { // all of the water is taken up by the sediments
         m_total_input_ghosts(i,j) = 0.0;
      } else { // determine how much water was taken up by the sediments and subtract it from the total water

        double water_in_till = ((Wtill_before - Wtill_after)  / icedt + tillwat_decay_rate) * m_till_cover(i,j);
        m_total_input_ghosts(i,j) -= water_in_till;

      }


 




    }
  }


  // cheat

//  if (m_t < 5000.) {
 //   m_total_input_ghosts.set(0.0);
 // }



  m_total_input_ghosts.update_ghosts();

  // we're going to need the potential gradient
  potential_gradient(m_hydro_gradient_dir);
  surface_gradient(m_surface_gradient_dir);

  // find the magnitude of the potential gradient

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

     m_hydro_gradient(i,j) = sqrt(pow(m_hydro_gradient_dir(i,j).v,2.0) + pow(m_hydro_gradient_dir(i,j).u,2.0));
  }

  // sort the permutation array

  int i_perm, j_perm;

  int i_offset = m_grid->xs();
  int j_offset = m_grid->ys();
  int sub_width_i = m_grid->xm();
  int sub_width_j = m_grid->ym();

  int i_current, j_current, i_next, j_next, i_compare, j_compare, i_now, j_now;


  int counter = 0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();


      i_now = i;
      j_now = j;
      cell_coordinates(m_gradient_permutation(i,j), sub_width_i, sub_width_j, i_offset, j_offset, i_current, j_current);

      int backwards_count = counter - 1;

	bool found_back;

      if(backwards_count < 0) {
        found_back = true;
      } else {
        found_back = false;
      }

      while (! found_back) {


       cell_coordinates(double(backwards_count), sub_width_i, sub_width_j, i_offset, j_offset, i_next, j_next); // get the next permutation cell



        cell_coordinates(m_gradient_permutation(i_next,j_next), sub_width_i, sub_width_j, i_offset, j_offset, i_compare, j_compare); // convert permutation to cell coordinates

//  m_log->message(2,
 //            "* %i %i %i %i %f %f \n", i_now, j_now, i_next, j_next, m_hydro_gradient(i_current, j_current), m_hydro_gradient(i_compare, j_compare));

        if( m_hydro_gradient(i_current, j_current) < m_hydro_gradient(i_compare, j_compare)) { // swap if true

          double temp_permutation = m_gradient_permutation(i_next, j_next);
          m_gradient_permutation(i_next, j_next) = m_gradient_permutation(i_now, j_now);
          m_gradient_permutation(i_now, j_now) = temp_permutation;
          i_now = i_next;
          j_now = j_next;
          backwards_count--;

         } else { // found the natural position for the sort
           found_back = true;
         }

         if(backwards_count < 0) {
           found_back = true;
         }


      }

      counter++;
  }


  // The remaining water input rate is supplemented by adding the nearest upstream cell's water. Ideally you would
  // want to include the entire pathway up to the ice sheet dome, but that would require a lot of additional calculations
  // and crosstalk between processors. As a result, the amount of water is surely going to be underestimated (but
  // hopefully not by an amount that is important!)


  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double x_amount;
    double y_amount;
    double angle;

    if ( fabs(m_surface_gradient_dir(i,j).v) > 0.01 || fabs(m_surface_gradient_dir(i,j).u) > 0.01 ) { // if the gradient is not that much, then it should not transmit water
      angle = atan2(m_surface_gradient_dir(i,j).v, m_surface_gradient_dir(i,j).u);

      x_amount = m_surface_gradient_dir(i,j).u; //cos(angle);
      y_amount = m_surface_gradient_dir(i,j).v; //sin(angle);

//  m_log->message(2,
//             "* %i %i %f %f %f\n", i, j, x_amount, y_amount, fabs(x_amount) + fabs(y_amount));

      // normalize

      double sum_components = fabs(x_amount) + fabs(y_amount);
      x_amount = x_amount / sum_components;
      y_amount = y_amount / sum_components;

//  m_log->message(2,
//             "* %i %i %f %f %f\n", i, j, x_amount, y_amount, fabs(x_amount) + fabs(y_amount));

    } else {
      angle = 0.0;
      x_amount = 0.0;
      y_amount = 0.0;

    }

    int ghost_index_x = 1;
    int ghost_index_y = 1;

    if(x_amount < 0.0) {
      ghost_index_x = -1;
    }

    if(y_amount < 0.0) {
      ghost_index_y = -1;
    }


//  m_log->message(2,
//             "* %i %i %i %i %f %f %f %f %f\n", i, j, ghost_index_x, ghost_index_y, x_amount, y_amount, m_hydro_gradient_dir(i,j).v, m_hydro_gradient_dir(i,j).u, angle);

    // add the upstream water to the current cell


/////    m_total_input_ghosts(i,j) += m_total_input_ghosts(i+ghost_index_x,j) * fabs(x_amount) + m_total_input_ghosts(i,j+ghost_index_y) * fabs(y_amount);


//  m_log->message(2,
//             "* %i %i %f %f %f\n", i, j, m_total_input_ghosts(i,j), Wtill_before, Wtill_after);

  }

//  m_log->message(2,
//             "* finished calculating m_total_input_ghosts ...\n");


  // now we have the water input rate, it is possible to calculate the tunnel shape and the critical tunnel stability 

  const double bump_ratio = bump_amplitude / bedrock_wavelength;

  update_velbase_mag(m_velbase_mag);

//  m_log->message(2,
//             "* finished updating velocity ...\n");

  // need to grab the overburden pressure for the calculation of the hydrology scheme
  overburden_pressure(m_pressure_temp);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // pressure gradient magnitude
//    m_hydro_gradient(i,j) = sqrt(pow(m_hydro_gradient_dir(i,j).v,2.0) + pow(m_hydro_gradient_dir(i,j).u,2.0));

    // surface gradient magnitude
    m_surface_gradient(i,j) = sqrt(pow(m_surface_gradient_dir(i,j).v,2.0) + pow(m_surface_gradient_dir(i,j).u,2.0));

//  m_log->message(2,
//             "* m_hydro_gradient: %f\n", m_hydro_gradient(i,j));

//  m_log->message(2,
//             "* m_total_input_ghosts: %f\n", m_total_input_ghosts(i,j));

    // Water discharge
    // Arnold and Sharp (2002) assumed the volume water flux was the same through tunnels and cavities



    m_volume_water_flux(i,j) = m_total_input_ghosts(i,j) * pow(tunnel_spacing,2);

      if(m_volume_water_flux(i,j) > 1000.0) {

 // m_log->message(2,
  //           "* error: %i %i %f %e\n", i, j, m_volume_water_flux(i,j), m_total_input_ghosts(i,j));
  //         throw RuntimeError::formatted(PISM_ERROR_LOCATION, "total input is way too high for some reason.\n" );
      }

//  m_log->message(2,
//             "* m_volume_water_flux: %f\n", m_volume_water_flux(i,j));

    // Rothleisburger tunnel radius
    // equation A.10 from Arnold and Sharp (2002)

    if(m_surface_gradient(i,j) > 1.0e-8) {

      m_tunnel_cross_section(i,j) = pow(channel_flow_constant * pow(m_volume_water_flux(i,j),2.0) / (rho_w * g * m_surface_gradient(i,j)), 3.0/8.0);

//  m_log->message(2,
//             "* %i %i %f %f %f %f\n",i, j, m_volume_water_flux(i,j), pow(m_volume_water_flux(i,j),2.0), m_hydro_gradient(i,j), m_tunnel_cross_section(i,j));
    } else{
      m_tunnel_cross_section(i,j) = 0.0;
    }

//  m_log->message(2,
//             "* m_tunnel_cross_section: %i %i %f\n",i, j, m_tunnel_cross_section(i,j));

    // Tunnel effective pressure
    // Equation A.8 from Arnold and Sharp (2002)

    double effective_pressure_tunnel;
    if(m_total_input_ghosts(i,j) > 1.0e-12 &&  m_tunnel_cross_section(i,j) > 1.0e-8) {
      effective_pressure_tunnel = pow( (rho_w * g * m_surface_gradient(i,j) * m_volume_water_flux(i,j)) / 
                                         (rho_i * arrhenius_parameter * latent_heat * m_tunnel_cross_section(i,j)), (1.0 / Glen_exponent));
    } else {
      effective_pressure_tunnel = 0.0;
    }


    // Fowler suggested that the effective pressure become atmospheric, but I am setting it to be some fraction of the overburden
    if(effective_pressure_tunnel / m_pressure_temp(i,j) > max_effective_pressure_ratio ) {
      effective_pressure_tunnel = max_effective_pressure_ratio * m_pressure_temp(i,j);
    }


//  m_log->message(2,
//             "* effective_pressure_tunnel: %f\n", effective_pressure_tunnel);

    // Cavity effective pressure
    // Equation A.11 from Arnold and Sharp (2002) (note that the equation is wrong in the paper, see equation 4.16 in Fowler (1987))
    double effective_pressure_cavity;

    double number_of_cavities = cavity_spacing *  tunnel_spacing;

//    double number_of_cavities = tunnel_spacing / (

    if(m_total_input_ghosts(i,j) > 1e-12 &&  m_tunnel_cross_section(i,j) > 1.0e-8) {
      effective_pressure_cavity = shadowing_function * pow((  (rho_w * g * m_surface_gradient(i,j)) / (rho_i * arrhenius_parameter * latent_heat) * 
                                                              ( m_volume_water_flux(i,j) / (number_of_cavities * cavity_area) ) ), (1.0 / Glen_exponent));
    } else {
      effective_pressure_cavity = 0.0;
    }

//  m_log->message(2,
//             "* effective_pressure_cavity: %f\n", effective_pressure_cavity);

    // tunnel stability critical value
    // equation A.13 from Arnold and Sharp (2002)
    double critical_stability = pow((3.0 * Glen_exponent * m_tunnel_cross_section(i,j) / (cavity_area * cavity_spacing * tunnel_spacing)), ((4.0-mu)/mu));

//  m_log->message(2,
//            "* critical_stability: %f\n", critical_stability);

    // tunnel stability value
    // equation A.12 from Arnold and Sharp (2002)
    double tunnel_stability = bump_ratio * m_velbase_mag(i,j) / ( bedrock_wavelength * arrhenius_parameter * pow(effective_pressure_tunnel, Glen_exponent));

//  m_log->message(2,
//             "* tunnel_stability: %f\n", tunnel_stability);

    if(m_total_input_ghosts(i,j) < 1e-12) { // essentially no water available

      m_hydrology_effective_pressure(i,j) = 0.0;

      m_hydrosystem(i,j) = 0.;

    } else if(tunnel_stability < critical_stability) { // tunnels are stable, so take that as the effective pressure

      m_hydrology_effective_pressure(i,j) = effective_pressure_tunnel;
      m_hydrosystem(i,j) = 1.;

    } else { // cavity system

      m_hydrology_effective_pressure(i,j) = effective_pressure_cavity;
      m_hydrosystem(i,j) = 2.;

    }

    m_hydrology_fraction_overburden(i,j) = m_hydrology_effective_pressure(i,j) / m_pressure_temp(i,j);

//  m_log->message(2,
//             "* %i %i %f %f %f %f\n", i, j, effective_pressure_tunnel / m_pressure_temp(i,j), effective_pressure_cavity / m_pressure_temp(i,j), m_volume_water_flux(i,j), m_tunnel_cross_section(i,j));
//  m_log->message(2,
//             "* %i %i %f %f %f\n", i, j, tunnel_stability, critical_stability, m_hydro_gradient(i,j));

//  m_log->message(2,
//             "* %i %i %f %f %f\n", i, j, m_tunnel_cross_section(i,j), number_of_cavities * cavity_area, m_volume_water_flux(i,j));

  }


//  m_log->message(2,
//             "* finished hydrologyEvan::update_impl ...\n");
}


// returns the coordinate of a cell
void hydrologyEvan::cell_coordinates(double in_number, int number_i, int number_j, int i_offset, int j_offset, int& i, int& j){

  
  j = floor(in_number / double(number_i)); // number of rows


  i = rint(in_number - double(j) * double(number_i)); // column offset

  i = i + i_offset;
  j = j + j_offset;

//  m_log->message(2,
 //             "** %f %i %i %i\n", in_number, number_i, j_offset, j);

}


} // end of namespace hydrology
} // end of namespace pism
