// Copyright (C) 2012-2017 PISM Authors
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


#include "hydrology_diagnostics.hh"

namespace pism {
namespace hydrology {

Hydrology_bwat::Hydrology_bwat(const Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "bwat")};
  set_attrs("thickness of transportable water in subglacial layer", "", "m", "m", 0);
}

IceModelVec::Ptr Hydrology_bwat::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bwat", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];
  model->subglacial_water_thickness(*result);
  return result;
}

Hydrology_bwp::Hydrology_bwp(const Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "bwp")};
  set_attrs("pressure of transportable water in subglacial layer", "", "Pa", "Pa", 0);
}


IceModelVec::Ptr Hydrology_bwp::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bwp", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];
  model->subglacial_water_pressure(*result);
  return result;
}


Hydrology_bwprel::Hydrology_bwprel(const Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "bwprel")};
  set_attrs("pressure of transportable water in subglacial layer as fraction of the overburden pressure", "",
            "", "", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}


IceModelVec::Ptr Hydrology_bwprel::compute_impl() const {
  double fill_value = m_fill_value;

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bwprel", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  IceModelVec2S Po;
  Po.create(m_grid, "Po_temporary", WITHOUT_GHOSTS);

  model->subglacial_water_pressure(*result);
  model->overburden_pressure(Po);

  IceModelVec::AccessList list{result.get(), &Po};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (Po(i,j) > 0.0) {
      (*result)(i,j) /= Po(i,j);
    } else {
      (*result)(i,j) = fill_value;
    }
  }

  return result;
}


Hydrology_effbwp::Hydrology_effbwp(const Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "effbwp")};
  set_attrs("effective pressure of transportable water in subglacial layer (overburden pressure minus water pressure)",
            "", "Pa", "Pa", 0);
}


IceModelVec::Ptr Hydrology_effbwp::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "effbwp", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  IceModelVec2S P;
  P.create(m_grid, "P_temporary", WITHOUT_GHOSTS);

  model->subglacial_water_pressure(P);
  model->overburden_pressure(*result);
  result->add(-1.0, P);  // result <-- result + (-1.0) P = Po - P

  return result;
}


Hydrology_hydrobmelt::Hydrology_hydrobmelt(const Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "hydrobmelt")};
  set_attrs("the version of bmelt seen by the hydrology model",
            "", "m s-1", "m year-1", 0);
}


IceModelVec::Ptr Hydrology_hydrobmelt::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "hydrobmelt", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];
  // the value reported diagnostically is merely the last value filled
  result->copy_from(model->m_bmelt_local);

  return result;
}


Hydrology_hydroinput::Hydrology_hydroinput(const Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "hydroinput")};
  set_attrs("total water input into subglacial hydrology layer",
            "", "m s-1", "m year-1", 0);
}


IceModelVec::Ptr Hydrology_hydroinput::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "hydroinput", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];
  // the value reported diagnostically is merely the last value filled
  result->copy_from(model->m_total_input);

  return result;
}


Hydrology_wallmelt::Hydrology_wallmelt(const Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "wallmelt")};
  set_attrs("wall melt into subglacial hydrology layer from (turbulent) dissipation of energy in transportable water",
            "", "m s-1", "m year-1", 0);
}


IceModelVec::Ptr Hydrology_wallmelt::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "wallmelt", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];
  model->wall_melt(*result);
  return result;
}


MCHydrology_ice_free_land_loss::MCHydrology_ice_free_land_loss(const Routing *m)
  : TSDiag<TSFluxDiagnostic, Routing>(m, "hydro_ice_free_land_loss") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("long_name",
                              "rate of liquid water loss from subglacial hydrology into "
                              "cells with mask as ice free land");
}

double MCHydrology_ice_free_land_loss::compute() {
  return model->boundary_mass_accounting().ice_free_land_loss;
}

MCHydrology_ocean_loss::MCHydrology_ocean_loss(const Routing *m)
  : TSDiag<TSFluxDiagnostic, Routing>(m, "hydro_ocean_loss") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("long_name",
                             "rate of liquid water loss from subglacial hydrology into "
                             "cells with mask as ocean");
}

double MCHydrology_ocean_loss::compute() {
  return model->boundary_mass_accounting().ocean_loss;
}

MCHydrology_negative_thickness_gain::MCHydrology_negative_thickness_gain(const Routing *m)
  : TSDiag<TSFluxDiagnostic, Routing>(m, "hydro_negative_thickness_gain") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("long_name",
                             "rate of non-conserving liquid water gain from subglacial "
                             "hydrology transportable water thickness coming out negative "
                             "during time step, and being projected up to zero");
}

double MCHydrology_negative_thickness_gain::compute() {
  return model->boundary_mass_accounting().negative_thickness_gain;
}

MCHydrology_null_strip_loss::MCHydrology_null_strip_loss(const Routing *m)
  : TSDiag<TSFluxDiagnostic, Routing>(m, "hydro_null_strip_loss") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("long_name",
                             "rate of liquid water loss from subglacial hydrology into "
                             "cells inside the null strip");
}

double MCHydrology_null_strip_loss::compute() {
  return model->boundary_mass_accounting().null_strip_loss;
}



// hydrologyEvan

hydrology_type::hydrology_type(const hydrologyEvan *m)
  : Diag<hydrologyEvan>(m) {

  // set metadata:

  m_vars = {SpatialVariableMetadata(m_sys, "hydrology_type")};

  set_attrs("type of hydrology, 0 dry 1 tunnels 2 cavity 3 overburden 4 almost floating", "",
            "1", "1", 0);
}


IceModelVec::Ptr hydrology_type::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "hydrology_type", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  model->hydrology_type(*result);

  return result;
}

} // end of namespace hydrology
} // end of namespace pism
