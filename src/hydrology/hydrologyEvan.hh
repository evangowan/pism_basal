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

#ifndef _PISMHYDROLOGYEVAN_H_
#define _PISMHYDROLOGYEVAN_H_

#include "Hydrology.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/Component.hh"
#include "pism/coupler/SurfaceModel.hh"

namespace pism {

class IceModelVec2T;

namespace stressbalance {
class StressBalance;
}

//! @brief Sub-glacial hydrology models and related diagnostics.
namespace hydrology {


// This is the same as NullTransport, except that the entire water thickness is saved instead of just the water in the till
class hydrologyEvan : public NullTransport {
public:
  hydrologyEvan(IceGrid::ConstPtr g, stressbalance::StressBalance *sb);
  virtual ~hydrologyEvan();
  virtual void init();


  virtual void write_model_state_impl(const PIO &output) const;
  virtual void define_model_state_impl(const PIO &output) const;


  virtual void get_SedimentDistribution(IceModelVec2S &result);
  virtual void potential_gradient(IceModelVec2S &result_u, IceModelVec2S &result_v);
  virtual void surface_gradient(IceModelVec2V &result);
  virtual void bed_gradient(IceModelVec2V &result);
  virtual void get_EffectivePressure(IceModelVec2S &result);
  virtual void get_hydrology_type(IceModelVec2S &result) const;
  virtual void get_volume_water_flux(IceModelVec2S &result) const;
  virtual void update_velbase_mag(IceModelVec2S &result);



protected:

  virtual void get_input_rate(double hydro_t, double hydro_dt, IceModelVec2S &result);
//  virtual MaxTimestep max_timestep_impl(double t);
  //! Solves an implicit step of a highly-simplified ODE.
  virtual void update_impl(double icet, double icedt);

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  virtual void update_surface_runoff(IceModelVec2S &result);


  IceModelVec2S m_till_cover;      // fraction of surface covered by till
  IceModelVec2S m_velbase_mag2; // velocity at the base of the ice sheet
  IceModelVec2S m_hydrology_effective_pressure; // effective pressure at the base due to hydrology

  // need to get basal sliding velocity (thus speed):
  stressbalance::StressBalance* m_stressbalance;


  IceModelVec2S m_volume_water_flux, m_hydrosystem;


  // diagnostic stuff

  friend class hydrology_type;

  friend class volume_water_flux;

//  virtual void hydrology_type_impl(IceModelVec2S &result) const = 0;


private:

  IceModelVec2S m_melt_rate_local, m_hydro_gradient, m_pressure_temp, m_total_input_ghosts, m_total_input_ghosts_temp, m_tunnel_cross_section,  m_surface_gradient, m_surface_elevation_temp,  m_bed_elevation_temp, m_hydrology_fraction_overburden, m_gradient_permutation, m_processor_mask, m_offset_mask_u, m_width_mask_u, m_offset_mask_v, m_width_mask_v, m_hydro_gradient_dir_u, m_hydro_gradient_dir_v, m_gradient_temp_u, m_gradient_temp_v;

  IceModelVec2V  m_surface_gradient_dir, m_surface_gradient_temp, m_bed_gradient_temp;// directional gradient

  void cell_coordinates(double in_number, int number_i, int number_j, int i_offset, int j_offset, int& i, int& j);

  int high_i_check(int i);
  int high_j_check(int j);
  int low_check(int i);
  void finite_difference(double point_array[3][3], double& u, double& v);
  void gradient_five_point(double point_array[5][5], double& u, double& v);

protected:

  petsc::Vec::Ptr m_offset_mask_u_p0, m_width_mask_u_p0, m_offset_mask_v_p0, m_width_mask_v_p0, m_processor_mask_p0, m_gradient_permutation_p0, m_total_input_ghosts_p0, m_total_input_ghosts_temp_p0, m_hydro_gradient_p0, m_hydro_gradient_dir_u_p0, m_hydro_gradient_dir_v_p0;

};


} // end of namespace hydrology
} // end of namespace pism


#endif /* _PISMHYDROLOGYEVAN_H_ */
