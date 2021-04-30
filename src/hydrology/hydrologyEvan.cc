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
#include "hydrology_diagnostics.hh"
#include "hydrologyEvan.hh"
#include "pism/util/Mask.hh"
#include "pism/util/error_handling.hh"
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

hydrologyEvan :: hydrologyEvan(IceGrid::ConstPtr g, stressbalance::StressBalance *sb)
	: NullTransport(g)
   {

  m_log->message(2,
             "* Starting hydrologyEvan ...\n");


  m_stressbalance = sb;

  if(m_stressbalance) {
  m_log->message(2,
             "* Stress balance model detected in the initialization\n");


  } else {
  m_log->message(2,
             "* Stress balance model not detected in the initialization\n");
  }


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


  m_surface_elevation_temp.create(m_grid, "surface_elevation_temp", WITH_GHOSTS, 2);
  m_surface_elevation_temp.set_attrs("internal",
                       "temporary surface_elevation",
                       "m", "");

  m_bed_elevation_temp.create(m_grid, "bed_elevation_temp", WITH_GHOSTS, 2);
  m_bed_elevation_temp.set_attrs("internal",
                       "temporary bed_elevation",
                       "m", "");




  m_till_cover.create(m_grid, "tillcover", WITHOUT_GHOSTS);
  m_till_cover.set_attrs("model_state",
                       "fraction of surface that is covered in till",
                       "1", "");
  m_till_cover.set_time_independent(true);

  m_gradient_permutation.create(m_grid, "gradient_permutation", WITHOUT_GHOSTS);
  m_gradient_permutation.set_attrs("internal",
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

  m_velbase_mag2.create(m_grid, "velbase_mag2", WITHOUT_GHOSTS);
  m_velbase_mag2.set_attrs("model_state",
                        "ice sliding speed seen by subglacial hydrology",
                        "m s-1", "");
  m_velbase_mag2.metadata().set_double("valid_min", 0.0);


  m_hydrology_effective_pressure.create(m_grid, "hydrology_effective_pressure", WITHOUT_GHOSTS);
  m_hydrology_effective_pressure.set_attrs("model_state",
                        "effective pressure due to hydrology",
                        "Pa", "");
  m_hydrology_effective_pressure.metadata().set_double("valid_min", 0.0);

  m_hydrology_fraction_overburden.create(m_grid, "fraction_overburden", WITHOUT_GHOSTS);
  m_hydrology_fraction_overburden.set_attrs("internal",
                        "fraction of effective pressure to overburden",
                        "1", "");
  m_hydrology_fraction_overburden.metadata().set_double("valid_min", 0.0);



 m_hydrosystem.create(m_grid, "hydrology_type", WITHOUT_GHOSTS);
 m_hydrosystem.set_attrs("model_state",
                        "type of hydrology, 0 dry 1 tunnels 2 cavity 3 overburden 4 almost floating",
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
  m_gradient_permutation_p0 = m_gradient_permutation.allocate_proc0_copy();
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

  const IceModelVec2S        &temp_thk = *m_grid->variables().get_2d_scalar("thk");
  double rho_i       = m_config->get_double("constants.ice.density");
  double rho_w       = m_config->get_double("constants.fresh_water.density");
  double g           = m_config->get_double("constants.standard_gravity");


  InputOptions opts = process_input_options(m_grid->com);

  if (opts.type == INIT_RESTART) {
    m_till_cover.read(opts.filename, opts.record);

    m_volume_water_flux.read(opts.filename, opts.record);
    m_hydrosystem.read(opts.filename, opts.record);

    m_velbase_mag2.read(opts.filename, opts.record);

    first_time = true;

  } else if (opts.type == INIT_BOOTSTRAP) {
    m_till_cover.regrid(opts.filename, OPTIONAL, till_cover);

    m_volume_water_flux.regrid(opts.filename, OPTIONAL, 0.0);
    m_hydrosystem.regrid(opts.filename, OPTIONAL, 0.0);

    m_velbase_mag2.regrid(opts.filename, OPTIONAL, 0.0);

    first_time = true;

  } else {
    m_till_cover.set(till_cover);

    m_volume_water_flux.set(0.0);
    m_hydrosystem.set(0.0);
    m_velbase_mag2.set(0.0);

    first_time = false;
  }


  // initialize the effective pressure



  IceModelVec::AccessList list;

  if (opts.type == INIT_RESTART) {
    m_hydrology_effective_pressure.read(opts.filename, opts.record);


  } else if (opts.type == INIT_BOOTSTRAP) {
    m_hydrology_effective_pressure.regrid(opts.filename, OPTIONAL, 36100800);


  } else {

    list.add(m_hydrology_effective_pressure);
    list.add(temp_thk);

    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        m_hydrology_effective_pressure(i,j) = temp_thk(i,j) * rho_i * g;

      }
    } catch (...) {
      loop.failed();
    }
    loop.check();

  }



  // initialize the permutation array, after the first sorting it should not take long to sort



  int counter = 0;

  list.add(m_gradient_permutation);
  list.add(m_processor_mask);
  list.add(m_offset_mask_u);
  list.add(m_offset_mask_v);
  list.add(m_width_mask_u);
  list.add(m_width_mask_v);

  list.add(temp_thk);




  int i_offset = m_grid->xs();
  int j_offset = m_grid->ys();
  int sub_width_i = m_grid->xm();
  int sub_width_j = m_grid->ym();


  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();


      m_gradient_permutation(i, j) = (double)counter;
      m_processor_mask(i,j) = m_grid -> rank();
      m_offset_mask_u(i,j) = i_offset;
      m_offset_mask_v(i,j) = j_offset;
      m_width_mask_u(i,j) = sub_width_i;
      m_width_mask_v(i,j) = sub_width_j;
      counter++;

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

//  m_log->message(2,
 //            "* Initializing Evan's null-transport subglacial hydrology model ...\n");



}


void hydrologyEvan::define_model_state_impl(const PIO &output) const {

  m_till_cover.define(output);
  m_hydrology_effective_pressure.define(output);
  m_total_input_ghosts.define(output);
  m_total_input_ghosts_temp.define(output);
  m_volume_water_flux.define(output);
  m_melt_rate_local.define(output);
  m_hydrosystem.define(output);
  m_velbase_mag2.define(output);
  m_hydro_gradient.define(output);
  m_hydro_gradient_dir_u.define(output);
  m_hydro_gradient_dir_v.define(output);
  m_tunnel_cross_section.define(output);
  m_hydrology_fraction_overburden.define(output);
//  m_processor_mask.define(output);  // does not need to be output
//  m_pressure_temp.define(output); // does not need to be output
//  m_gradient_permutation.define(output); // does not need to be output

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
  m_velbase_mag2.write(output);
  m_hydro_gradient.write(output);
  m_hydro_gradient_dir_u.write(output);
  m_hydro_gradient_dir_v.write(output);
  m_tunnel_cross_section.write(output);
  m_hydrology_fraction_overburden.write(output);
//  m_processor_mask.write(output);  // does not need to be output
//  m_pressure_temp.write(output); // does not need to be output
//  m_gradient_permutation.write(output); // does not need to be output
}


// Diagnostics


std::map<std::string, Diagnostic::Ptr> hydrologyEvan::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"hydrology_type",                     Diagnostic::Ptr(new hydrology_type(this))},
    {"volume_water_flux",                     Diagnostic::Ptr(new volume_water_flux(this))}
  };
  return combine(result, Hydrology::diagnostics_impl());
}


// Return the hydrology type.
void hydrologyEvan::get_hydrology_type(IceModelVec2S &result) const {
  result.copy_from(m_hydrosystem);
}

void hydrologyEvan::get_volume_water_flux(IceModelVec2S &result) const {
  result.copy_from(m_volume_water_flux);
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


  IceModelVec2S ice_vel_temp(m_grid, "ice_temp", WITHOUT_GHOSTS);;
  IceModelVec::AccessList list;
  list.add(ice_vel_temp);

  ice_vel_temp.set_to_magnitude(m_stressbalance->advective_velocity());

  // the velocity in PISM is defined in seconds, I want it to be m/yr
  double inverse_seconds_in_year = 1.0/ 365.0*24.0*3600.0;
  ice_vel_temp.scale(inverse_seconds_in_year);

  result.copy_from(ice_vel_temp);

}

// grabs the surface runoff
void hydrologyEvan::update_surface_runoff(IceModelVec2S &result) {

  m_surfaceT->runoff_rate(result);

}


// Calculate the potential gradient at the base of the ice sheet, based on Equation 6.13 in Cuffey and Paterson (2010).
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

  {
    IceModelVec::AccessList list;
    list.add(m_surface_gradient_temp);
    list.add(m_bed_gradient_temp);
    list.add(m_gradient_temp_u);
    list.add(m_gradient_temp_v);


    // grab the bed and surface gradients
    surface_gradient(m_surface_gradient_temp);
    bed_gradient(m_bed_gradient_temp);

    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        // Equation 6.13 in Cuffy and Paterson (2010)
        // The flotation fraction is the ratio of the pressure of water to the pressure of ice, and will influence the effect of the bed gradient on the total potential gradient
        // At a default, it is set to 0.8, which gives means the bed slope needs to be 2.7 times greater than the surface slope to have an equal influence on the direction of water flow.
        m_gradient_temp_u(i,j) = rho_i_g * (flotation_fraction * m_surface_gradient_temp(i,j).u + (density_ratio - flotation_fraction) * m_bed_gradient_temp(i,j).u);
        m_gradient_temp_v(i,j) = rho_i_g * (flotation_fraction * m_surface_gradient_temp(i,j).v + (density_ratio - flotation_fraction) * m_bed_gradient_temp(i,j).v);

      }
    } catch (...) {
      loop.failed();
    }
    loop.check();
  }

  result_u.copy_from(m_gradient_temp_u);
  result_v.copy_from(m_gradient_temp_v);


}

// find the gradient of the surface
void hydrologyEvan::surface_gradient(IceModelVec2V &result) {


  const IceModelVec2S &surface_elevation = *m_grid->variables().get_2d_scalar("surface_altitude");

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  double u, v;


  IceModelVec::AccessList list;
  list.add(m_surface_elevation_temp);
  list.add(result);

  m_surface_elevation_temp.copy_from(surface_elevation);
  m_surface_elevation_temp.update_ghosts();

  ParallelSection loop(m_grid->com);
  try {

    // Set gradient_grid_width to 3 for third order finite difference method from Skidmore (1989)
    // Set gradient_grid_width to 5 for 5x5 grid least squares
    int gradient_grid_width = 5; 
    double point_store3[3][3];
    double point_store5[5][5];
    double u, v;
    int i_check, j_check;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      for (int k=0; k < gradient_grid_width; k++) {
        for (int l=0; l < gradient_grid_width; l++) {

         i_check = i + k - ((gradient_grid_width-1)/2);
         i_check = low_check(i_check);
         i_check = high_i_check(i_check);     

         j_check = j + l - ((gradient_grid_width-1)/2);
         j_check = low_check(j_check);
         j_check = high_j_check(j_check); 


         if (k > 0 && k < 4 && l > 0 && l < 4) {

           point_store3[k-1][l-1] = m_surface_elevation_temp(i_check,j_check);
         }

         point_store5[k][l] = m_surface_elevation_temp(i_check,j_check);

        }
      }

      if(gradient_grid_width == 3) {
        finite_difference(point_store3, u, v);
      } else {

        gradient_five_point(point_store5, u, v);

      }
      result(i,j).u = u;
      result(i,j).v = v;

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

}

// find the gradient of the bed
void hydrologyEvan::bed_gradient(IceModelVec2V &result) {

  const IceModelVec2S &bed_elevation = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  double u, v;

  IceModelVec::AccessList list;
  list.add(result);
  list.add(m_bed_elevation_temp);

  m_bed_elevation_temp.copy_from(bed_elevation);
  m_bed_elevation_temp.update_ghosts();

  ParallelSection loop(m_grid->com);
  try {
    // Set gradient_grid_width to 3 for third order finite difference method from Skidmore (1989)
    // Set gradient_grid_width to 5 for 5x5 grid least squares
    int gradient_grid_width = 5;
    double point_store3[3][3];
    double point_store5[5][5];
    double u, v;
    int i_check, j_check;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();



      for (int k=0; k < gradient_grid_width; k++) {
        for (int l=0; l < gradient_grid_width; l++) {

         i_check = i + k - ((gradient_grid_width-1)/2);
         i_check = low_check(i_check);
         i_check = high_i_check(i_check);     

         j_check = j + l - ((gradient_grid_width-1)/2);
         j_check = low_check(j_check);
         j_check = high_j_check(j_check); 

         if (k > 0 && k < 4 && l > 0 && l < 4) {

           point_store3[k-1][l-1] = m_bed_elevation_temp(i_check,j_check);
         }

         point_store5[k][l] = m_bed_elevation_temp(i_check,j_check);

        }
      }

      if(gradient_grid_width == 3) {
        finite_difference(point_store3, u, v);
      } else {

        gradient_five_point(point_store5, u, v);

      }
      result(i,j).u = u;
      result(i,j).v = v;

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

int hydrologyEvan::high_i_check(int i) {

  int num_i = m_grid->Mx();
  int i_high;


  if(i >= num_i) {
    i_high = num_i;
  } else {
    i_high = i;
  }

  return i_high;
}

int hydrologyEvan::high_j_check(int j) {

  int num_j = m_grid->My();
  int j_high;


  if(j >= num_j) {
    j_high = num_j;
  } else {
    j_high = j;
  }

  return j_high;
}

int hydrologyEvan::low_check(int i) {

  int i_low;
      if(i < 0) {
        i_low = 0;
      } else {
        i_low = i;
      }
  return i_low;

}

// Third order finite difference method for calculationg gradient, equations 3 and 4 in Skidmore (1989)
// This was originally how I calculated the gradient, and it is still possible to use it if you change
// the grid size in the gradient calculations.
void hydrologyEvan::finite_difference(double point_array[3][3], double& u, double& v) {

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();


  u = ((point_array[2][2] + 2.0 * point_array[2][1] + point_array[2][0]) -
			      (point_array[0][2] + 2.0 * point_array[0][1] + point_array[0][0])) / (8.0 * dx);

  v = ((point_array[2][2] + 2.0 * point_array[1][2] + point_array[0][2]) -
			      (point_array[2][0] + 2.0 * point_array[1][0] + point_array[0][0])) / (8.0 * dy);

}


// Find the gradient via least squares.
// This is the default way to calculate gradient.
void hydrologyEvan::gradient_five_point(double point_array[5][5], double& u, double& v) {

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  int array_size = 5;
  int array_size_squared = array_size*array_size;
  int plane_variables = 3;
  int half_array = (5-1) / 2;

  double B_matrix[array_size_squared];
  double R_T_R[plane_variables][plane_variables];
  double R[array_size_squared][plane_variables];
  double R_T[plane_variables][array_size_squared];
  double R_T_R_inverse[plane_variables][plane_variables];

  double R_T_R_inverse_R_T[plane_variables][array_size_squared];

  double plane_solution[plane_variables];

  // create R and R transpose arrays, and B_matrix
  int counter = 0;
  for (int i = 0; i < array_size; i++) {
    for (int j = 0; j < array_size; j++) {

      R[counter][0] = 1;
      R[counter][1] = double(i-half_array)*dx;
      R[counter][2] = double(j-half_array)*dy;

      R_T[0][counter] = R[counter][0];
      R_T[1][counter] = R[counter][1];
      R_T[2][counter] = R[counter][2];

      B_matrix[counter] = point_array[i][j];

      counter++;
    }
  }

  // multiply R transpose and R
  for (int i = 0; i < plane_variables; i++) {
    for (int j = 0; j < plane_variables; j++) {

      double temp_storage = 0;

       for (int k = 0; k < array_size_squared; k++) {

         temp_storage = temp_storage + R_T[i][k] * R[k][j];

       }

       R_T_R[i][j] = temp_storage; 

    }
  }

  // find the inverse of the R_T_R matrix
  double det_R_T_R;
  {
    double a = R_T_R[0][0];
    double b = R_T_R[1][0];
    double c = R_T_R[2][0];
    double d = R_T_R[0][1];
    double e = R_T_R[1][1];
    double f = R_T_R[2][1];
    double g = R_T_R[0][2];
    double h = R_T_R[1][2];
    double i = R_T_R[2][2];



    R_T_R_inverse[0][0] = e * i - f * h;      // A
    R_T_R_inverse[0][1] = -(d * i - f * g);  // B
    R_T_R_inverse[0][2] = d * h - e * g;      // C
    R_T_R_inverse[1][0] = -(b * i - c * h);  // D
    R_T_R_inverse[1][1] = a * i - c * g;      // E
    R_T_R_inverse[1][2] = -( a * h - b * g);  // F
    R_T_R_inverse[2][0] = b * f - c * e;      // G
    R_T_R_inverse[2][1] = -( a * f - c * d);  // H
    R_T_R_inverse[2][2] = a * e - b * d;      // I

  // Rule of Sarrus for finding the determinant of a 3x3 matrix
    det_R_T_R = a * R_T_R_inverse[0][0] + b * R_T_R_inverse[0][1] + c * R_T_R_inverse[0][2];

  }

  // shouldn't need to check of the determinant is zero, because they are coordinates of a grid
  for (int i = 0; i < plane_variables; i++) {
    for (int j = 0; j < plane_variables; j++) {

      double before = R_T_R_inverse[i][j];
      R_T_R_inverse[i][j] = R_T_R_inverse[i][j] / det_R_T_R;

    }

  }

  // Multiply R_T_R_inverse_R by R_t
  for (int i = 0; i < plane_variables; i++) {
    for (int k = 0; k < array_size_squared; k++) {
    

      double temp_storage = 0;

      for (int j = 0; j < plane_variables; j++) {

         temp_storage = temp_storage + R_T_R_inverse[i][j] * R_T[j][k] ;

      }

      R_T_R_inverse_R_T[i][k] = temp_storage; 

    }
  }

  // find the variables for the equation of a plane
  for (int i = 0; i < plane_variables; i++) {

    double temp_storage = 0;

    for (int k = 0; k < array_size_squared; k++) {
    
         temp_storage = temp_storage + R_T_R_inverse_R_T[i][k] * B_matrix[k] ;

    }

    plane_solution[i] = temp_storage; 

  }

  // finally, output the gradient
  u = plane_solution[1];
  v = plane_solution[2];

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
  //  It should be updating the surface runoff
    update_surface_runoff(m_melt_rate_local);

  } else {
    //  Surface model not detected
    m_melt_rate_local.set(0.0);

  }

  // gets water from the surface and adds the water generated from melting at the base due to friction and geothermal heating
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.grounded_ice(i, j)) {

        result(i,j) = (use_const) ? const_bmelt : m_bmelt_local(i,j); // get the melt water from the base

        result(i,j) += m_melt_rate_local(i,j) *  const_water_from_surface; // add on the meltwater from the surface

      } else {
        result(i,j) = 0.0;
      }

    }

  } catch (...) {
    loop.failed();
  }
  loop.check();

}


//! Determine the hydrology at the base
/*!

This is a non-conservative hydrology model, based on a combination of the sediment deformation
model that is default in PISM, and the tunnel/cavity hydrology model by Schoof (2010).
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


  const double tillwat_max = m_config->get_double("hydrology.tillwat_max"),
               tillwat_decay_rate           = m_config->get_double("hydrology.tillwat_decay_rate"), // this C factor should probably be made to be spatially variable too, enhanced if sand, for now constant
               rho_w       = m_config->get_double("constants.fresh_water.density"),
               g           = m_config->get_double("constants.standard_gravity"),
		   tunnel_spacing = m_config->get_double("hydrology.tunnel_spacing"),
               Glen_exponent = m_config->get_double("stress_balance.ssa.Glen_exponent"),
               arrhenius_parameter = m_config->get_double("flow_law.isothermal_Glen.ice_softness"),
               rho_i       = m_config->get_double("constants.ice.density"),
               latent_heat = m_config->get_double("constants.fresh_water.latent_heat_of_fusion"),
               ice_thickness_threshold = m_config->get_double("hydrology.ice_thickness_threshold")

  ;


  const double
    dx = m_grid->dx(),
    dy = m_grid->dy(); // hopefully the grid is always square


  // grab the input rate, which will be the sum of the meltwater generated from basal heating, and meltwater from the surface
  get_input_rate(icet,icedt, m_total_input_ghosts);


  const IceModelVec2CellType &mask  = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec2S        &temp_thk = *m_grid->variables().get_2d_scalar("thk");


  // First, we need to find out if the water input is sufficient to fill up the till in the cell.
  // If there is excess water after filling the till, it is sent to the next step.

  {

    IceModelVec::AccessList list{&m_Wtil, &mask, &m_till_cover, &m_total_input_ghosts};

    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();


        if ( mask.ice_free(i,j) && ! mask.ocean(i,j)) {
          m_Wtil(i,j) = 0.0;

        } else if (mask.ocean(i,j)) {

             // We assume that sediments under water are saturated. This is important for when the ice advances, 
             // so it doesn't have to fill up.
	        m_Wtil(i,j) = tillwat_max; 

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
    } catch (...) {
      loop.failed();
    }
    loop.check();

  }




  // we're going to need the potential gradient
  potential_gradient(m_hydro_gradient_dir_u,m_hydro_gradient_dir_v);
  surface_gradient(m_surface_gradient_dir);

  // find the magnitude of the potential gradient
  {
    IceModelVec::AccessList list{&m_hydro_gradient, &m_hydro_gradient_dir_v, &m_hydro_gradient_dir_u, &m_total_input_ghosts, &temp_thk};
    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        // magnitude of the potential gradient
        m_hydro_gradient(i,j) = sqrt(pow(m_hydro_gradient_dir_v(i,j),2.0) + pow(m_hydro_gradient_dir_u(i,j),2.0));

        // if the ice thickness is really small, don't bother distributing the water
        if(temp_thk(i,j) <= ice_thickness_threshold) {
          m_total_input_ghosts(i,j) = 0.0;
        }
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();
  }

  m_total_input_ghosts.update_ghosts();

  // sort the permutation array, finding the order from lowest gradient to highest gradient
  {
    int i_perm, j_perm;

    int i_offset = m_grid->xs();
    int j_offset = m_grid->ys();
    int sub_width_i = m_grid->xm();
    int sub_width_j = m_grid->ym();

    int i_current, j_current, i_next, j_next, i_compare, j_compare, i_now, j_now;

    IceModelVec::AccessList list{&m_gradient_permutation, &m_hydro_gradient, &m_processor_mask};

    // I'm not entirely sure if this block is needed, but it has been so long since I coded it that I don't want to touch it
    {
      ParallelSection loop(m_grid->com);
      try {
        for (Points p(*m_grid); p; p.next()) {
          const int i = p.i(), j = p.j();

          cell_coordinates(m_gradient_permutation(i,j), sub_width_i, sub_width_j, i_offset, j_offset, i_current, j_current);

        }
      } catch (...) {
        loop.failed();
      }
      loop.check();
    }

    // This block is the meat of it, it sorts the permutation array for each processor
    int counter = 0;
    ParallelSection loop(m_grid->com);
    try {
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

          // get the next permutation cell
          cell_coordinates(double(backwards_count), sub_width_i, sub_width_j, i_offset, j_offset, i_next, j_next); 

          // convert permutation to cell coordinates
          cell_coordinates(m_gradient_permutation(i_next,j_next), sub_width_i, sub_width_j, i_offset, j_offset, i_compare, j_compare); 

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
    } catch (...) {
      loop.failed();
    }
    loop.check();

    // Again, I don't know if this block is necessary, but I don't want to change anything at this point.
    {
      ParallelSection loop(m_grid->com);
      try {
        for (Points p(*m_grid); p; p.next()) {
          const int i = p.i(), j = p.j();

          cell_coordinates(m_gradient_permutation(i,j), sub_width_i, sub_width_j, i_offset, j_offset, i_current, j_current);
        
        }
      } catch (...) {
        loop.failed();
      }
      loop.check();
    }
  }



  // In this section, we route the water, which it is easiest done in a serial process.
  // Everything is moved to one processor for this calculation

  m_total_input_ghosts_temp.copy_from(m_total_input_ghosts);

  int num_i = m_grid->Mx();
  int num_j = m_grid->My();

  {

    // this puts all the variables needed onto a single processor
    m_processor_mask.put_on_proc0(*m_processor_mask_p0);
    m_offset_mask_u.put_on_proc0(*m_offset_mask_u_p0);
    m_offset_mask_v.put_on_proc0(*m_offset_mask_v_p0);
    m_width_mask_u.put_on_proc0(*m_width_mask_u_p0);
    m_width_mask_v.put_on_proc0(*m_width_mask_v_p0);
    m_total_input_ghosts_temp.put_on_proc0(*m_total_input_ghosts_temp_p0);
    m_hydro_gradient.put_on_proc0(*m_hydro_gradient_p0);
    m_hydro_gradient_dir_u.put_on_proc0(*m_hydro_gradient_dir_u_p0);
    m_hydro_gradient_dir_v.put_on_proc0(*m_hydro_gradient_dir_v_p0);
    m_gradient_permutation.put_on_proc0(*m_gradient_permutation_p0);


    // This is the section of the code that is running on a single processor
    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {

        // creates vectors that can be used on a single processor 
        petsc::VecArray processor_mask_p0(*m_processor_mask_p0);
        petsc::VecArray offset_mask_u_p0(*m_offset_mask_u_p0);
        petsc::VecArray offset_mask_v_p0(*m_offset_mask_v_p0);
        petsc::VecArray width_mask_u_p0(*m_width_mask_u_p0);
        petsc::VecArray width_mask_v_p0(*m_width_mask_v_p0);
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
        double* total_input_ghosts_temp_vec =  total_input_ghosts_temp_p0.get();
        double* gradient_permutation_vec =  gradient_permutation_p0.get();
        double* hydro_gradient_vec =  hydro_gradient_p0.get();
        double* hydro_gradient_vec_u = hydro_gradient_p0_vec_u.get();
        double* hydro_gradient_vec_v = hydro_gradient_p0_vec_v.get();

        double pi = 3.14159265358979;

        int total_nodes = m_grid->xm() * m_grid->ym();

        int number_of_processors = m_grid->size();

        // temporary arrays to store information on each subdomain
        int processor_point_counter[number_of_processors];
        int processor_width_mask_u[number_of_processors];
        int processor_width_mask_v[number_of_processors];
        int processor_offset_mask_u[number_of_processors];
        int processor_offset_mask_v[number_of_processors];

        int max_points = 0;

        int processor, vector_index;

        // initialize to zero
        for (int proc_counter = 0; proc_counter < number_of_processors; proc_counter++) {
          processor_point_counter[proc_counter] = 0;

        }

        // Fill up those arrays
        for (int j = 0; j < num_j; j++) {
          for (int i = 0; i < num_i; i++) {

            vector_index = j*num_i + i;


            processor = processor_mask_vec[vector_index];
            processor_point_counter[processor]++; // increment the number of points in that particular processor

            if (processor_point_counter[processor] > max_points) {
              max_points = processor_point_counter[processor];
            }

            processor_width_mask_u[processor] = width_mask_u_vec[vector_index];
            processor_width_mask_v[processor] = width_mask_v_vec[vector_index];
            processor_offset_mask_u[processor] = offset_mask_u_vec[vector_index];
            processor_offset_mask_v[processor] = offset_mask_v_vec[vector_index];

          }
        }

        // read in the permutation arrays and separate per processor
        double serial_permutation[number_of_processors][max_points];
        double gradient_storage[number_of_processors][max_points]; // used to reduce the amount of calculations
        int i_array[number_of_processors][max_points];
        int j_array[number_of_processors][max_points];

        // Initialize the array counters, first point will be zero, which will give the correct index
        for(int k = 0; k < number_of_processors; k++) {

          processor_point_counter[k] =-1; 

        }

        int i_temp, j_temp, permutation_index;

        // this loop places the sorted gradient values from each processor into a storage array for sorting
        for (int j = 0; j < num_j; j++) {
          for (int i = 0; i < num_i; i++) {

            vector_index = j*num_i + i;

            processor = processor_mask_vec[vector_index];
            processor_point_counter[processor]++; // increment the number of points in that particular processor


            cell_coordinates(gradient_permutation_vec[vector_index], processor_width_mask_u[processor], processor_width_mask_v[processor], processor_offset_mask_u[processor], processor_offset_mask_v[processor], i_temp, j_temp);

            permutation_index = j_temp * num_i + i_temp;

            gradient_storage[processor][processor_point_counter[processor]] = hydro_gradient_vec[permutation_index]; // should be highest to lowest

            serial_permutation[processor][processor_point_counter[processor]] = permutation_index;

          }
        }

        // move the point count to a separate variable
        int max_point_count[number_of_processors];

        // rezero the point counter. C++ is not as nice as Fortran for this
        for (int counter = 0; counter < number_of_processors; counter ++) {
          max_point_count[counter] = processor_point_counter[counter];
          processor_point_counter[counter] = 0;

        }


        // distribute the water in the direction of largest gradient magnitude
        bool finished = false;

        int lowest_index, lowest_processor;
        double lowest_gradient;

        while (! finished) {

          // find the lowest index of the lowest gradient

          bool found_first = false;

          // This loop finds the smallest gradient value
          for (int processor_counter = 0; processor_counter < number_of_processors; processor_counter++) {

            // check if the processor still has points to check
            if(processor_point_counter[processor_counter] <= max_point_count[processor_counter]) { 

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

            int index = serial_permutation[lowest_processor][lowest_index];

            // distribute water if the gradient is significant enough
            if (gradient_storage[lowest_processor][lowest_index] > 1.0) {



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

              int max_ind = num_i * num_j;
              int neighbor_index =  index - 1; // goes left
              if ((neighbor_index >= 0) and (neighbor_index < max_ind) and (i_temp > 0)) {
                total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index]
                                                              + water_store_multiplier[0]*total_input_ghosts_temp_vec[index]
                                                              + water_store_multiplier[11]*total_input_ghosts_temp_vec[index];
              }

              neighbor_index =  index + (-1) + (-1) * num_i; // goes down left
              if ((neighbor_index >= 0) and (neighbor_index < max_ind) and (j_temp > 0) and (i_temp > 0)) {
                total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index]
                                                              + water_store_multiplier[1]*total_input_ghosts_temp_vec[index];
              }

              neighbor_index =  index +  (-1) * num_i; // goes down
              if ((neighbor_index >= 0) and (neighbor_index < max_ind) and (j_temp > 0)) {
                total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index]
                                                              + water_store_multiplier[2]*total_input_ghosts_temp_vec[index]
                                                              + water_store_multiplier[3]*total_input_ghosts_temp_vec[index];
              }

              neighbor_index =  index + (1) + (-1) * num_i; // goes down right
              if ((neighbor_index >= 0) and (neighbor_index < max_ind) and (j_temp > 0) and (i_temp < num_i)) {
                total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index]
                                                              + water_store_multiplier[4]*total_input_ghosts_temp_vec[index];
              }

              neighbor_index =  index + 1; // goes right
              if ((neighbor_index >= 0) and (neighbor_index < max_ind) and (i_temp < num_i)) {
                total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index]
                                                              + water_store_multiplier[5]*total_input_ghosts_temp_vec[index]
                                                              + water_store_multiplier[6]*total_input_ghosts_temp_vec[index];
              }

              neighbor_index =  index + (1) + (1) * num_i; // goes up right
              if ((neighbor_index >= 0) and (neighbor_index < max_ind) and (i_temp < num_i) and (j_temp < num_j)) {
                total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index]
                                                              + water_store_multiplier[7]*total_input_ghosts_temp_vec[index];
              }

              neighbor_index =  index +  (1) * num_i; // goes up
              if ((neighbor_index >= 0) and (neighbor_index < max_ind) and (j_temp < num_j)) {
                total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index]
                                                              + water_store_multiplier[8]*total_input_ghosts_temp_vec[index]
                                                              + water_store_multiplier[9]*total_input_ghosts_temp_vec[index];
              }

              neighbor_index =  index + (-1) + (1) * num_i; // goes up left
              if ((neighbor_index >= 0) and (neighbor_index < max_ind) and (i_temp > 0) and (j_temp < num_j)) {
                total_input_ghosts_temp_vec[neighbor_index] = total_input_ghosts_temp_vec[neighbor_index]
                                                              + water_store_multiplier[10]*total_input_ghosts_temp_vec[index];
              }

            } else { // set the storage to zero

              total_input_ghosts_temp_vec[index] = 0.;

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

  //  we need to grab the basal velocity, which is used to determine the hydrology regime
  if(first_time) {
    first_time = false;
  } else {
    update_velbase_mag(m_velbase_mag2);
  }

  // need to grab the overburden pressure for the calculation of the hydrology scheme
  {
    ParallelSection loop(m_grid->com);
    try {

      IceModelVec::AccessList list{&m_pressure_temp, &temp_thk, &m_total_input_ghosts, &mask};


      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        m_pressure_temp(i,j) = temp_thk(i,j) * rho_i * g;

        if( ! mask.grounded(i,j) ) {
          m_total_input_ghosts(i,j) = 0.0;
        }

      }
    } catch (...) {
      loop.failed();
    }
    loop.check();

  }

  m_hydrology_fraction_overburden.set(1.0);


  {
    ParallelSection loop(m_grid->com);
    try {

      IceModelVec::AccessList list{&m_surface_gradient, &m_surface_gradient_dir, &m_volume_water_flux, &m_total_input_ghosts, &m_tunnel_cross_section, &m_pressure_temp, &m_velbase_mag2, &m_hydrology_effective_pressure, &m_hydrology_fraction_overburden, &mask, &m_hydrosystem, &m_hydro_gradient};


      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        // surface gradient magnitude
        m_surface_gradient(i,j) = sqrt(pow(m_surface_gradient_dir(i,j).v,2.0) + pow(m_surface_gradient_dir(i,j).u,2.0));

        // Water discharge

        // let us use Schoof's parameterization

        double c1 = 1.0 / (rho_i * latent_heat);
        double c2 = 2.0 * arrhenius_parameter * pow(Glen_exponent,-Glen_exponent);
        double pi = 3.14159265359;
        double f = 0.1;
        double c3 = pow(2.0,(1.0/4.0)) * pow(pi+2.0,(1/2)) / (pow(pi,(1/4)) * pow(rho_w*f,(1/2)));
        double alpha = 5.0/4.0;
        double protrusion_height = 0.1;
        double psi_exponent = -1.0 / (2.0 * alpha);


        double effective_pressure_temp;

        // equation 3 in Schoof 2010, used to determine the hydrology system type
        double Qc = m_velbase_mag2(i,j) * protrusion_height / (c1 * (alpha-1.0) * m_hydro_gradient(i,j));
        m_volume_water_flux(i,j) = m_total_input_ghosts(i,j) * pow(dx,2) / ( dx / tunnel_spacing);


        m_hydrosystem(i,j) = 0.;
        if(m_volume_water_flux(i,j) <= 1e-12) { // essentially no water available

          effective_pressure_temp = m_pressure_temp(i,j);

          m_hydrosystem(i,j) = 0.;
        } else {

          // equation 2 in Schoof 2010


          effective_pressure_temp = pow(( c1 * m_volume_water_flux(i,j) * m_hydro_gradient(i,j) + m_velbase_mag2(i,j) * protrusion_height ) /
                                         ( c2 * pow(c3, -1.0 / alpha) * pow(m_volume_water_flux(i,j), 1.0/alpha) * pow(m_hydro_gradient(i,j), psi_exponent))
                                         , (1.0/Glen_exponent));


          if(m_volume_water_flux(i,j) > Qc) {
             m_hydrosystem(i,j) = 1.; // tunnels
          } else {
             m_hydrosystem(i,j) = 2.; // cavities
          }

        }

	  // prevents the effective pressure from exceeding the overburden presssure

        m_hydrology_effective_pressure(i,j) = effective_pressure_temp;
        if (m_hydrology_effective_pressure(i,j) > m_pressure_temp(i,j)) {
          m_hydrology_effective_pressure(i,j) = m_pressure_temp(i,j);
          m_hydrosystem(i,j) = 3.;
        }

        // prevents issues if the effective pressure is essentially zero

        if(mask.grounded_ice(i,j) && m_hydrology_effective_pressure(i,j) < 0.01 * m_pressure_temp(i,j)) {
         m_hydrology_effective_pressure(i,j) = 0.01 * m_pressure_temp(i,j);
         m_hydrosystem(i,j) = 4.;
        }

        if(m_pressure_temp(i,j) > 0.0) {
          m_hydrology_fraction_overburden(i,j) = m_hydrology_effective_pressure(i,j) / m_pressure_temp(i,j);
        }

      }

    } catch (...) {
      loop.failed();
    }
    loop.check();

  }


  m_log->message(2,
             "* finished hydrologyEvan::update_impl ...\n");
}


// This function returns the coordinate of a cell
void hydrologyEvan::cell_coordinates(double in_number, int number_i, int number_j, int i_offset, int j_offset, int& i, int& j){


  j = int(in_number) / number_i; // number of rows

  i = in_number - (j * number_i); // column offset

  i = i + i_offset;
  j = j + j_offset;

}




hydrology_type::hydrology_type(const hydrologyEvan *m)
  : Diag<hydrologyEvan>(m) {

  // set metadata:

  m_vars = {SpatialVariableMetadata(m_sys, "hydrology_type")};

  set_attrs("type of hydrology, 0 dry 1 tunnels 2 cavity 3 overburden 4 almost floating", "",
            "1", "1", 0);
}


IceModelVec::Ptr hydrology_type::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "hydrology_type", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  model->get_hydrology_type(*result);

  return result;
}




volume_water_flux::volume_water_flux(const hydrologyEvan *m)
  : Diag<hydrologyEvan>(m) {

  // set metadata:

  m_vars = {SpatialVariableMetadata(m_sys, "volume_water_flux")};

  set_attrs("volume water flux through tunnel system", "",
            "m3 s-1", "m3 s-1", 0);
}


IceModelVec::Ptr volume_water_flux::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "volume_water_flux", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  model->get_volume_water_flux(*result);

  return result;
}



} // end of namespace hydrology
} // end of namespace pism
