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
#include <array>
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

  m_gradient_temp_u.create(m_grid, "gradient_temp_u", WITHOUT_GHOSTS);
  m_gradient_temp_u.set_attrs("internal",
                        "gradient at the base",
                        "Pa m-1", "");

  m_gradient_temp_v.create(m_grid, "gradient_temp_v", WITHOUT_GHOSTS);
  m_gradient_temp_v.set_attrs("internal",
                        "gradient at the base",
                        "Pa m-1", "");


  m_hydro_gradient.create(m_grid, "hydro_gradient", WITHOUT_GHOSTS);
  m_hydro_gradient.set_attrs("internal",
                        "pressure gradient at the base",
                        "Pa m-1", "");

  m_hydro_gradient_dir_u.create(m_grid, "hydro_gradient_dir_u", WITHOUT_GHOSTS);
  m_hydro_gradient_dir_u.set_attrs("internal",
                        "potential gradient at the base in directional components, u",
                        "Pa m-1", "");

  m_hydro_gradient_dir_v.create(m_grid, "hydro_gradient_dir_v", WITHOUT_GHOSTS);
  m_hydro_gradient_dir_v.set_attrs("internal",
                        "potential gradient at the base in directional components, v",
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

  m_total_input_ghosts_temp.create(m_grid, "total_input_temp", WITH_GHOSTS, stencil_width); // need ghosts here
  m_total_input_ghosts_temp.set_attrs("internal",
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


 m_processor_mask.create(m_grid, "processor_mask", WITHOUT_GHOSTS);
 m_processor_mask.set_attrs("internal",
                        "assigns a number for each processor",
                        "1", "");
 m_processor_mask.metadata().set_double("valid_min", 0.0);

 m_offset_mask_u.create(m_grid, "offset_mask_u", WITHOUT_GHOSTS);
 m_offset_mask_u.set_attrs("internal",
                        "assigns the grid offset for each processor",
                        "1", "");
 m_offset_mask_u.metadata().set_double("valid_min", 0.0);

 m_width_mask_u.create(m_grid, "width_mask_u", WITHOUT_GHOSTS);
 m_width_mask_u.set_attrs("internal",
                        "assigns the grid width for each processor",
                        "1", "");
 m_width_mask_u.metadata().set_double("valid_min", 0.0);

 m_offset_mask_v.create(m_grid, "offset_mask_v", WITHOUT_GHOSTS);
 m_offset_mask_v.set_attrs("internal",
                        "assigns the grid offset for each processor",
                        "1", "");
 m_offset_mask_v.metadata().set_double("valid_min", 0.0);

 m_width_mask_v.create(m_grid, "width_mask_v", WITHOUT_GHOSTS);
 m_width_mask_v.set_attrs("internal",
                        "assigns the grid width for each processor",
                        "1", "");
 m_width_mask_v.metadata().set_double("valid_min", 0.0);


// processor 0 memory

//m_processor_mask
//m_offset_mask
//m_width_mask
//m_total_input_ghosts
// m_gradient_permutation
// m_hydro_gradient

  m_processor_mask_p0 = m_processor_mask.allocate_proc0_copy();
  m_offset_mask_u_p0 = m_offset_mask_u.allocate_proc0_copy();
  m_offset_mask_v_p0 = m_offset_mask_v.allocate_proc0_copy();
  m_width_mask_u_p0 = m_width_mask_u.allocate_proc0_copy();
  m_width_mask_v_p0 = m_width_mask_v.allocate_proc0_copy();
  m_total_input_ghosts_p0 = m_total_input_ghosts.allocate_proc0_copy();
  m_total_input_ghosts_temp_p0 = m_total_input_ghosts_temp.allocate_proc0_copy();
  m_gradient_permutation_p0 =m_gradient_permutation.allocate_proc0_copy();
  m_hydro_gradient_p0 = m_hydro_gradient.allocate_proc0_copy();
  m_hydro_gradient_dir_u_p0 = m_hydro_gradient_dir_u.allocate_proc0_copy();
  m_hydro_gradient_dir_v_p0 = m_hydro_gradient_dir_v.allocate_proc0_copy();


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
  list.add(m_processor_mask);
  list.add(m_offset_mask_u);
  list.add(m_offset_mask_v);
  list.add(m_width_mask_u);
  list.add(m_width_mask_v);

  int i_offset = m_grid->xs();
  int j_offset = m_grid->ys();
  int sub_width_i = m_grid->xm();
  int sub_width_j = m_grid->ym();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_log->message(2,
             "%i %i %i\n", i, j, counter);

    m_gradient_permutation(i, j) = counter;
    m_processor_mask(i,j) = m_grid -> rank();
    m_offset_mask_u(i,j) = i_offset;
    m_offset_mask_v(i,j) = j_offset;
    m_width_mask_u(i,j) = sub_width_i;
    m_width_mask_v(i,j) = sub_width_j;
    counter++;

  }

  m_log->message(2,
             "* Finished initializing the permutation array ...\n");

}


void hydrologyEvan::define_model_state_impl(const PIO &output) const {

  m_till_cover.define(output);
  m_hydrology_effective_pressure.define(output);
  m_total_input_ghosts.define(output);
  m_total_input_ghosts_temp.define(output);
  m_volume_water_flux.define(output);
  m_melt_rate_local.define(output);
  m_hydrosystem.define(output);
  m_velbase_mag.define(output);
  m_hydro_gradient.define(output);
  m_hydro_gradient_dir_u.define(output);
  m_hydro_gradient_dir_v.define(output);
  m_tunnel_cross_section.define(output);
  m_hydrology_fraction_overburden.define(output);
  m_processor_mask.define(output);

}

void hydrologyEvan::write_model_state_impl(const PIO &output) const {
  m_Wtil.write(output);
  m_till_cover.write(output);
  m_hydrology_effective_pressure.write(output);
  m_total_input_ghosts_temp.write(output);
  m_total_input_ghosts.write(output);
  m_volume_water_flux.write(output);
  m_melt_rate_local.write(output);
  m_hydrosystem.write(output);
  m_velbase_mag.write(output);
  m_hydro_gradient.write(output);
  m_hydro_gradient_dir_u.write(output);
  m_hydro_gradient_dir_v.write(output);
  m_tunnel_cross_section.write(output);
  m_hydrology_fraction_overburden.write(output);
  m_processor_mask.write(output);
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
void hydrologyEvan::potential_gradient(IceModelVec2S &result_u, IceModelVec2S &result_v) {

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
  list.add(m_gradient_temp_u);
  list.add(m_gradient_temp_v);


  // grab the bed and surface gradients
  surface_gradient(m_surface_gradient_temp);
  bed_gradient(m_bed_gradient_temp);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Equation 6.13 in Cuffy and Paterson (2010)
    // The flotation fraction is the ratio of the pressure of water to the pressure of ice, and will influence the effect of the bed gradient on the total potential gradient
    // At a default, it is set to 0.8, which gives means the bed slope needs to be 2.7 times greater than the surface slope to have an equal influence on the direction of water flow.
    m_gradient_temp_u(i,j) = rho_i_g * (flotation_fraction * m_surface_gradient_temp(i,j).u + (density_ratio - flotation_fraction) * m_bed_gradient_temp(i,j).u);
    m_gradient_temp_v(i,j) = rho_i_g * (flotation_fraction * m_surface_gradient_temp(i,j).v + (density_ratio - flotation_fraction) * m_bed_gradient_temp(i,j).v);


  }

  result_u.copy_from(m_gradient_temp_u);
  result_v.copy_from(m_gradient_temp_v);

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
  list.add(m_total_input_ghosts_temp);
  list.add(m_till_cover);
  list.add(m_hydro_gradient);
  list.add(m_hydro_gradient_dir_u);
  list.add(m_hydro_gradient_dir_v);
  list.add(m_surface_gradient);
  list.add(m_surface_gradient_dir);
  list.add(m_tunnel_cross_section);
  list.add(m_hydrology_effective_pressure);
  list.add(m_volume_water_flux);
  list.add(m_pressure_temp);
  list.add(m_hydrosystem);
  list.add(m_hydrology_fraction_overburden);
  list.add(m_gradient_permutation);
  list.add(m_processor_mask);
  list.add(m_offset_mask_u);
  list.add(m_offset_mask_v);
  list.add(m_width_mask_u);
  list.add(m_width_mask_v);





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
  potential_gradient(m_hydro_gradient_dir_u,m_hydro_gradient_dir_v);
  surface_gradient(m_surface_gradient_dir);

  // find the magnitude of the potential gradient

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

     m_hydro_gradient(i,j) = sqrt(pow(m_hydro_gradient_dir_v(i,j),2.0) + pow(m_hydro_gradient_dir_u(i,j),2.0));
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


//m_processor_mask
//m_offset_mask
//m_width_mask
//m_total_input_ghosts
// m_gradient_permutation

  // find the routing of water, it is easiest done in a serial way, so everything is moved to one processor for this calculation

// Note: temporary cheat, uncomment this line when not cheating
//  m_total_input_ghosts_temp.copy_from(m_total_input_ghosts);

//  m_total_input_ghosts_temp.set(0.1);


  int num_i = m_grid->Mx();
  int num_j = m_grid->My();
  double seconds_in_year = 365.0*24.0*3600.0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

     if(i < num_i/2 && j < num_j / 2) {
       m_total_input_ghosts_temp(i,j) = 1.0 / seconds_in_year;
     } else {
       m_total_input_ghosts_temp(i,j) = 0.0;
     }

  }
 

  m_log->message(2,
             "* starting serial process ...\n");
  {
    m_processor_mask.put_on_proc0(*m_processor_mask_p0);
    m_offset_mask_u.put_on_proc0(*m_offset_mask_u_p0);
    m_offset_mask_v.put_on_proc0(*m_offset_mask_v_p0);
    m_width_mask_u.put_on_proc0(*m_width_mask_u_p0);
    m_width_mask_v.put_on_proc0(*m_width_mask_v_p0);
    m_total_input_ghosts.put_on_proc0(*m_total_input_ghosts_p0);
    m_total_input_ghosts_temp.put_on_proc0(*m_total_input_ghosts_temp_p0);
    m_gradient_permutation.put_on_proc0(*m_gradient_permutation_p0);
    m_hydro_gradient.put_on_proc0(*m_hydro_gradient_p0);
    m_hydro_gradient_dir_u.put_on_proc0(*m_hydro_gradient_dir_u_p0);
    m_hydro_gradient_dir_v.put_on_proc0(*m_hydro_gradient_dir_v_p0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
  m_log->message(2,
             "* in serial process ...\n");



        petsc::VecArray processor_mask_p0(*m_processor_mask_p0);
        petsc::VecArray offset_mask_u_p0(*m_offset_mask_u_p0);
        petsc::VecArray offset_mask_v_p0(*m_offset_mask_v_p0);
        petsc::VecArray width_mask_u_p0(*m_width_mask_u_p0);
        petsc::VecArray width_mask_v_p0(*m_width_mask_v_p0);
        petsc::VecArray total_input_ghosts_p0(*m_total_input_ghosts_p0);
        petsc::VecArray total_input_ghosts_temp_p0(*m_total_input_ghosts_temp_p0);
        petsc::VecArray gradient_permutation_p0(*m_gradient_permutation_p0);
        petsc::VecArray hydro_gradient_p0(*m_hydro_gradient_p0);
        petsc::VecArray hydro_gradient_p0_vec_u(*m_hydro_gradient_dir_u_p0);
        petsc::VecArray hydro_gradient_p0_vec_v(*m_hydro_gradient_dir_v_p0);


        double* processor_mask_vec =  processor_mask_p0.get();
        double* offset_mask_u_vec =  offset_mask_u_p0.get();
        double* offset_mask_v_vec =  offset_mask_v_p0.get();
        double* width_mask_u_vec =  width_mask_u_p0.get();
        double* width_mask_v_vec =  width_mask_v_p0.get();
        double* total_input_ghosts_vec =  total_input_ghosts_p0.get();
        double* total_input_ghosts_temp_vec =  total_input_ghosts_temp_p0.get();
        double* gradient_permutation_vec =  gradient_permutation_p0.get();
        double* hydro_gradient_vec =  hydro_gradient_p0.get();
        double* hydro_gradient_vec_u = hydro_gradient_p0_vec_u.get();
        double* hydro_gradient_vec_v = hydro_gradient_p0_vec_v.get();

        double pi = 3.14159265358979;

  m_log->message(2,
             "* switched up memory ...\n");

        int total_nodes = m_grid->xm() * m_grid->ym();

        int number_of_processors = m_grid->size();

        // temporary arrays to store information on each subdomain

        int processor_point_counter[number_of_processors] = {0}; // initialize to zero
        int processor_width_mask_u[number_of_processors];
        int processor_width_mask_v[number_of_processors];
        int processor_offset_mask_u[number_of_processors];
        int processor_offset_mask_v[number_of_processors];

        int max_points = 0;

        int processor, vector_index;

        // Fill up those arrays





  m_log->message(2,
             "* assigning first arrays ...\n");



        for (int j = 0; j < num_j; j++) {
          for (int i = 0; i < num_i; i++) {

            vector_index = j*num_i + i;

 // m_log->message(2,
 //            "* %i %i %i %i\n", i, j, num_i, num_j);

            processor = processor_mask_vec[vector_index];
            processor_point_counter[processor]++; // increment the number of points in that particular processor

            if (processor_point_counter[processor] > max_points) {
              max_points = processor_point_counter[processor];
            }

            processor_width_mask_u[processor] = width_mask_u_vec[vector_index];
            processor_width_mask_v[processor] = width_mask_v_vec[vector_index];
            processor_offset_mask_u[processor] = offset_mask_u_vec[vector_index];
            processor_offset_mask_v[processor] = offset_mask_v_vec[vector_index];

//  m_log->message(2,
//             "* %i %i %i %i %f %i %f %f\n", i, j, num_i, num_j, vector_index, processor_mask_vec[vector_index], offset_mask_u_vec[vector_index], offset_mask_v_vec[vector_index]);
          }
        }

  m_log->message(2,
             "* assigned first arrays ...\n");

        // read in the permutation and separate per processor
        double serial_permutation[number_of_processors][max_points];
        double gradient_storage[number_of_processors][max_points]; // used to reduce the amount of calculations
        int i_array[number_of_processors][max_points];
        int j_array[number_of_processors][max_points];

        for(int k = 0; k < number_of_processors; k++) {

          processor_point_counter[k] =-1; // first point will be zero, which will give the correct index
        }

        int i_temp, j_temp, permutation_index;

  m_log->message(2,
             "* assigning second arrays ...\n");

        for (int j = 0; j < num_j; j++) {
          for (int i = 0; i < num_i; i++) {

            vector_index = j*num_i + i;

//  m_log->message(2,
//            "* %i %i %i %f\n", i, j, vector_index), processor_mask_vec[vector_index];

            processor = processor_mask_vec[vector_index];
            processor_point_counter[processor]++; // increment the number of points in that particular processor



//  m_log->message(2,
//             "* %f %i %i %f %f %i %i\n", gradient_permutation_vec[vector_index], processor_width_mask_u[processor], processor_width_mask_v[processor], processor_offset_mask_u[processor], processor_offset_mask_v[processor], i_temp, j_temp);
            
            cell_coordinates(gradient_permutation_vec[vector_index], processor_width_mask_u[processor], processor_width_mask_v[processor], processor_offset_mask_u[processor], processor_offset_mask_v[processor], i_temp, j_temp);

            permutation_index = j_temp * num_j + i_temp;

            gradient_storage[processor][processor_point_counter[processor]] = hydro_gradient_vec[permutation_index]; // should be highest to lowest

            serial_permutation[processor][processor_point_counter[processor]] = permutation_index;

 // m_log->message(2,
 //            "* %i %i %f\n", processor, processor_point_counter[processor], hydro_gradient_vec[permutation_index]);

          }
        }

        // move the point count to a separate variable
        int max_point_count[number_of_processors];

        // rezero the point counter. C++ is not as nice as Fortran for this
        for (counter = 0; counter < number_of_processors; counter ++) {
          max_point_count[counter] = processor_point_counter[counter];
          processor_point_counter[counter] = 0;

        }


     // distribute the water
  m_log->message(2,
             "* distributing water ...\n");

        bool finished = false;

        int lowest_index, lowest_processor;
        double lowest_gradient;

        while (! finished) {

          // find the lowest index of the lowest gradient



          bool found_first = false;

//  m_log->message(2,
//            "* number proc: %i\n", number_of_processors);
          for (int processor_counter = 0; processor_counter < number_of_processors; processor_counter++) {



            if(processor_point_counter[processor_counter] <= max_point_count[processor_counter]) { // check if the processor still has points to check



              if(! found_first) {
                lowest_gradient = gradient_storage[processor_counter][processor_point_counter[processor_counter]];
                lowest_index = processor_point_counter[processor_counter];
                lowest_processor = processor_counter;
                found_first = true;
              } else {
 
                if( gradient_storage[processor_counter][processor_point_counter[processor_counter]] < lowest_gradient) { // use this as the next point
                  lowest_gradient = gradient_storage[processor_counter][processor_point_counter[processor_counter]];
                  lowest_index = processor_point_counter[processor_counter];
                  lowest_processor = processor_counter;
                }

              }
            }

          }


          if (found_first) { 


            if (gradient_storage[lowest_processor][lowest_index] > 1e-9) { // distribute water

              int index = serial_permutation[lowest_processor][lowest_index];

              cell_coordinates(index, num_i, num_j, 0, 0, i_temp, j_temp);

              // find angular direction of the maximum gradient

              double direction = atan2(-hydro_gradient_vec_v[index], -hydro_gradient_vec_u[index]);

              double min_direction = direction - pi / 2.0;
              if(min_direction < -pi) {
                min_direction = 2.0 * pi + min_direction;
              }




              double max_direction = direction + pi / 2.0;
              if(max_direction > pi) {
                max_direction =    max_direction - 2.0 * pi;
              }

              bool check_two = false;
              double max_direction_2, min_direction_2;
              if(min_direction > 0.0 && max_direction < 0.0) { // need two checks
                check_two = true;
                min_direction_2 = -pi;
                max_direction_2 = max_direction;
                max_direction = pi;
              }


              double water_store_multiplier[12];
              double total = 0;
              double current_dir, next_dir, angle_to_min_current, angle_to_min_next, angle_to_max_current, angle_to_max_next;

              for ( int dir_counter = 0; dir_counter < 12; dir_counter++){

                current_dir = -pi + double(dir_counter) * pi / 6.0;
                next_dir = -pi + double(dir_counter+1) * pi / 6.0;

                bool add_to_distribution = false;

                if (min_direction > current_dir && min_direction < next_dir) { // move current dir to min_direction
                  current_dir = min_direction;
                }

                if (max_direction > current_dir && max_direction < next_dir) { // move next_dir to max_direction
                  next_dir = max_direction;
                }

                if(check_two) {
                  if (max_direction_2 > current_dir && max_direction_2 < next_dir) { // move next_dir to max_direction
                   next_dir = max_direction_2;
                  }
                }

                bool distribute = false;

                if(current_dir >= min_direction && next_dir <= max_direction) { // distribute water
                  distribute = true;
                }

                if(check_two) {
                  if(current_dir >= min_direction_2 && next_dir <= max_direction_2) { // distribute water
                    distribute = true;
                  }
                }

                if(distribute) { // distribute water

                  // this is a cos^2 function, integrated, and scaled to have the sum under the curve be 1
                  water_store_multiplier[dir_counter] = 2.0/pi * ( ((next_dir - direction) / 2.0 + sin(2.0*(next_dir-direction)) / 4.0) - 
                                                                  ((current_dir - direction) / 2.0 + sin(2.0*(current_dir-direction)) / 4.0));


                  total = total + water_store_multiplier[dir_counter];

                } else {

                  water_store_multiplier[dir_counter] = 0.0;

                }


              }

              int neighbor_index =  index - 1; // goes left
              total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index] 
                                                            + water_store_multiplier[0]*total_input_ghosts_temp_vec[index]
                                                            + water_store_multiplier[11]*total_input_ghosts_temp_vec[index];

              neighbor_index =  index + (-1) + (-1) * num_i; // goes down left
              total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index] 
                                                            + water_store_multiplier[1]*total_input_ghosts_temp_vec[index];

              neighbor_index =  index +  (-1) * num_i; // goes down
              total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index] 
                                                            + water_store_multiplier[2]*total_input_ghosts_temp_vec[index]
                                                            + water_store_multiplier[3]*total_input_ghosts_temp_vec[index];

              neighbor_index =  index + (1) + (-1) * num_i; // goes down right
              total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index] 
                                                            + water_store_multiplier[4]*total_input_ghosts_temp_vec[index];

              neighbor_index =  index + 1; // goes right
              total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index] 
                                                            + water_store_multiplier[5]*total_input_ghosts_temp_vec[index]
                                                            + water_store_multiplier[6]*total_input_ghosts_temp_vec[index];

              neighbor_index =  index + (1) + (1) * num_i; // goes up right
              total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index] 
                                                            + water_store_multiplier[7]*total_input_ghosts_temp_vec[index];

              neighbor_index =  index +  (1) * num_i; // goes up
              total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index] 
                                                            + water_store_multiplier[8]*total_input_ghosts_temp_vec[index]
                                                            + water_store_multiplier[9]*total_input_ghosts_temp_vec[index];

              neighbor_index =  index + (-1) + (1) * num_i; // goes up left
              total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index] 
                                                            + water_store_multiplier[10]*total_input_ghosts_temp_vec[index];

            }

            processor_point_counter[lowest_processor]++;

          } else {

            finished = true;

          }

        }



      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();

    m_total_input_ghosts_temp.get_from_proc0(*m_total_input_ghosts_temp_p0);

  }


  m_total_input_ghosts.copy_from(m_total_input_ghosts_temp);


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

    // Water discharge
    // Arnold and Sharp (2002) assumed the volume water flux was the same through tunnels and cavities



    m_volume_water_flux(i,j) = m_total_input_ghosts(i,j) * pow(tunnel_spacing,2);


    // Rothleisburger tunnel radius
    // equation A.10 from Arnold and Sharp (2002)

    if(m_surface_gradient(i,j) > 1.0e-8) {

      m_tunnel_cross_section(i,j) = pow(channel_flow_constant * pow(m_volume_water_flux(i,j),2.0) / (rho_w * g * m_surface_gradient(i,j)), 3.0/8.0);


    } else{
      m_tunnel_cross_section(i,j) = 0.0;
    }


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


    // Cavity effective pressure
    // Equation A.11 from Arnold and Sharp (2002) (note that the equation is wrong in the paper, see equation 4.16 in Fowler (1987))
    double effective_pressure_cavity;

    double number_of_cavities = cavity_spacing *  tunnel_spacing;


    if(m_total_input_ghosts(i,j) > 1e-12 &&  m_tunnel_cross_section(i,j) > 1.0e-8) {
      effective_pressure_cavity = shadowing_function * pow((  (rho_w * g * m_surface_gradient(i,j)) / (rho_i * arrhenius_parameter * latent_heat) * 
                                                              ( m_volume_water_flux(i,j) / (number_of_cavities * cavity_area) ) ), (1.0 / Glen_exponent));
    } else {
      effective_pressure_cavity = 0.0;
    }


    // tunnel stability critical value
    // equation A.13 from Arnold and Sharp (2002)
    double critical_stability = pow((3.0 * Glen_exponent * m_tunnel_cross_section(i,j) / (cavity_area * cavity_spacing * tunnel_spacing)), ((4.0-mu)/mu));


    // tunnel stability value
    // equation A.12 from Arnold and Sharp (2002)
    double tunnel_stability = bump_ratio * m_velbase_mag(i,j) / ( bedrock_wavelength * arrhenius_parameter * pow(effective_pressure_tunnel, Glen_exponent));


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
