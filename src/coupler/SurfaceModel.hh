// Copyright (C) 2008-2017 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

#ifndef __PISMSurfaceModel_hh
#define __PISMSurfaceModel_hh

/*!
 * This file should contain the class definition and nothing else.
 * Implementations should go in separate files.
 */

#include "pism/util/Component.hh"

namespace pism {

namespace atmosphere {
class AtmosphereModel;
}

class IceModelVec2S;

//! @brief Surface models and modifiers: provide top-surface
//! temperature, mass flux, liquid water fraction, mass and thickness of the surface
//! layer.
namespace surface {

//! \brief The interface of PISM's surface models.
class SurfaceModel : public Component_TS {
public:
  SurfaceModel(IceGrid::ConstPtr g);
  virtual ~SurfaceModel();

  void init();

  void attach_atmosphere_model(atmosphere::AtmosphereModel *input);

  // the interface:
  void mass_flux(IceModelVec2S &result) const;

  void temperature(IceModelVec2S &result) const;
  void liquid_water_fraction(IceModelVec2S &result) const;

  void layer_mass(IceModelVec2S &result) const;
  void layer_thickness(IceModelVec2S &result) const;

  void runoff_rate(IceModelVec2S &result) const; // added by Evan

protected:
  virtual void init_impl();

  virtual void attach_atmosphere_model_impl(atmosphere::AtmosphereModel *input);

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual MaxTimestep max_timestep_impl(double my_t) const;

  virtual void layer_thickness_impl(IceModelVec2S &result) const;
  virtual void layer_mass_impl(IceModelVec2S &result) const;

  virtual void temperature_impl(IceModelVec2S &result) const = 0;
  virtual void liquid_water_fraction_impl(IceModelVec2S &result) const;

  virtual void mass_flux_impl(IceModelVec2S &result) const = 0;

  virtual void runoff_rate_impl(IceModelVec2S &result) const; // added by Evan

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;
  virtual std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const;
protected:
  atmosphere::AtmosphereModel *m_atmosphere;
};

} // end of namespace surface
} // end of namespace pism

#endif  // __PISMSurfaceModel_hh
