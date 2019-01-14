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

#include "Hydrology.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "hydrology_diagnostics.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace hydrology {

Hydrology::Hydrology(IceGrid::ConstPtr g)
  : Component_TS(g) {

  m_surfaceT = NULL;
  m_inputtobed = NULL;
  m_hold_bmelt = false;

  m_total_input.create(m_grid, "total_input", WITHOUT_GHOSTS);
  m_total_input.set_attrs("internal",
                          "hydrology model workspace for total input rate into subglacial water layer",
                          "m s-1", "");

  m_bmelt_local.create(m_grid, "basal_melt_rate_grounded", WITHOUT_GHOSTS);
  m_bmelt_local.set_attrs("internal",
                          "hydrology model workspace for basal_melt_rate_grounded",
                          "m s-1", "");

  // *all* Hydrology classes have layer of water stored in till as a state variable
  m_Wtil.create(m_grid, "tillwat", WITHOUT_GHOSTS);
  m_Wtil.set_attrs("model_state",
                   "effective thickness of subglacial water stored in till",
                   "m", "");
  m_Wtil.metadata().set_double("valid_min", 0.0);
}


Hydrology::~Hydrology() {
  // empty
}


void Hydrology::init() {

  options::String bmelt_file("-hydrology_bmelt_file",
                             "Read time-independent values for basal_melt_rate_grounded from a file;"
                             " replaces basal_melt_rate_grounded computed through conservation of energy");
  // itb = input_to_bed
  options::String itb_file("-hydrology_input_to_bed_file",
                           "A time- and space-dependent file with amount of water"
                           " (depth per time) which should be added to the amount of water"
                           " at the ice sheet bed at the given location at the given time;"
                           " adds to basal_melt_rate_grounded");

  options::Real itb_period_years("-hydrology_input_to_bed_period",
                                 "The period (i.e. duration before repeat), in years,"
                                 " of -hydrology_input_to_bed_file data", 0.0);

  options::Real itb_reference_year("-hydrology_input_to_bed_reference_year",
                                   "The reference year for periodizing the"
                                   " -hydrology_input_to_bed_file data", 0.0);

  // the following are IceModelVec pointers into IceModel generally and are read by code in the
  // update() method at the current Hydrology time

  if (bmelt_file.is_set()) {
    m_log->message(2,
               "  option -hydrology_bmelt_file seen; reading basal_melt_rate_grounded from '%s'.\n", bmelt_file->c_str());
    m_bmelt_local.regrid(bmelt_file, CRITICAL);
    m_hold_bmelt = true;
  }


  if (itb_file.is_set()) {
    m_inputtobed_period = itb_period_years;
    m_inputtobed_reference_time = units::convert(m_sys, itb_reference_year, "years", "seconds");

    unsigned int buffer_size = (unsigned int) m_config->get_double("climate_forcing.buffer_size");

    PIO nc(m_grid->com, "netcdf3", itb_file, PISM_READONLY);
    unsigned int n_records = nc.inq_nrecords("inputtobed", "", m_sys);

    // if -..._period is not set, make n_records the minimum of the
    // buffer size and the number of available records. Otherwise try
    // to keep all available records in memory.
    if (not itb_period_years.is_set()) {
      n_records = std::min(n_records, buffer_size);
    }

    if (n_records == 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "can't find 'inputtobed' in -hydrology_input_to_bed"
                                    " file with name '%s'",
                                    itb_file->c_str());
    }

    m_log->message(2,
               "    option -hydrology_input_to_bed_file seen ... creating 'inputtobed' variable ...\n");
    m_log->message(2,
               "    allocating buffer space for n = %d 'inputtobed' records ...\n", n_records);
    m_inputtobed = new IceModelVec2T;
    m_inputtobed->set_n_records(n_records);
    m_inputtobed->create(m_grid, "inputtobed");
    m_inputtobed->set_attrs("climate_forcing",
                            "amount of water (depth per time like basal_melt_rate_grounded)"
                            " which should be put at the ice sheet bed",
                            "m s-1", "");
    m_log->message(2,
               "    reading 'inputtobed' variable from file '%s' ...\n",
               itb_file->c_str());
    m_inputtobed->init(itb_file, m_inputtobed_period, m_inputtobed_reference_time);
  }

  InputOptions opts = process_input_options(m_grid->com);

  double tillwat_default = m_config->get_double("bootstrapping.defaults.tillwat");

  switch (opts.type) {
  case INIT_RESTART:
    m_Wtil.read(opts.filename, opts.record);
    break;
  case INIT_BOOTSTRAP:
    m_Wtil.regrid(opts.filename, OPTIONAL, tillwat_default);
    break;
  case INIT_OTHER:
  default:
    m_Wtil.set(tillwat_default);
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  regrid("Hydrology", m_Wtil);
}


std::map<std::string, Diagnostic::Ptr> Hydrology::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"bwat",       Diagnostic::Ptr(new Hydrology_bwat(this))},
    {"bwp",        Diagnostic::Ptr(new Hydrology_bwp(this))},
    {"bwprel",     Diagnostic::Ptr(new Hydrology_bwprel(this))},
    {"effbwp",     Diagnostic::Ptr(new Hydrology_effbwp(this))},
    {"hydrobmelt", Diagnostic::Ptr(new Hydrology_hydrobmelt(this))},
    {"hydroinput", Diagnostic::Ptr(new Hydrology_hydroinput(this))},
    {"wallmelt",   Diagnostic::Ptr(new Hydrology_wallmelt(this))},
    {"tillwat",    Diagnostic::wrap(m_Wtil)},
  };
  return result;
}

void Hydrology::define_model_state_impl(const PIO &output) const {
  m_Wtil.define(output);
}

void Hydrology::write_model_state_impl(const PIO &output) const {
  m_Wtil.write(output);
}

//! Update the overburden pressure from ice thickness.
/*!
Uses the standard hydrostatic (shallow) approximation of overburden pressure,
  \f[ P_0 = \rho_i g H \f]
Accesses H=thk from Vars, which points into IceModel.
 */
void Hydrology::overburden_pressure(IceModelVec2S &result) const {
  // FIXME issue #15
  const IceModelVec2S *thk = m_grid->variables().get_2d_scalar("thk");

  result.copy_from(*thk);  // copies into ghosts if result has them
  result.scale(m_config->get_double("constants.ice.density") * m_config->get_double("constants.standard_gravity"));
}


//! Return the effective thickness of the water stored in till.
void Hydrology::till_water_thickness(IceModelVec2S &result) const {
  result.copy_from(m_Wtil);
}


//! Set the wall melt rate to zero.  (The most basic subglacial hydrologies have no lateral flux or potential gradient.)
void Hydrology::wall_melt(IceModelVec2S &result) const {
  result.set(0.0);
}


/*!
Checks \f$0 \le W_{til} \le W_{til}^{max} =\f$hydrology_tillwat_max.
 */
void Hydrology::check_Wtil_bounds() {
  double tillwat_max = m_config->get_double("hydrology.tillwat_max");

  IceModelVec::AccessList list(m_Wtil);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_Wtil(i,j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Hydrology: negative till water effective layer thickness Wtil(i,j) = %.6f m\n"
                                      "at (i,j)=(%d,%d)", m_Wtil(i,j), i, j);
      }

      if (m_Wtil(i,j) > tillwat_max) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Hydrology: till water effective layer thickness Wtil(i,j) = %.6f m exceeds\n"
                                      "hydrology_tillwat_max = %.6f at (i,j)=(%d,%d)",
                                      m_Wtil(i,j), tillwat_max, i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

}


//! Compute the total water input rate into the basal hydrology layer in the ice-covered region, allowing time-varying input from a file.
/*!
The user can specify the total of en- and supra-glacial drainage contributions
to subglacial hydrology in a time-dependent input file using option -hydrology_input_to_bed.
This method includes that possible input along with `basal_melt_rate_grounded` to get the total water
input into the subglacial hydrology.

This method crops the input rate to the ice-covered region.  It
also uses hydrology_const_bmelt if that is requested.

Call this method using the current \e hydrology time step.  This method
may be called many times per IceModel time step.  See update() method
in derived classes of Hydrology.
 */
void Hydrology::get_input_rate(double hydro_t, double hydro_dt,
                               IceModelVec2S &result) {
  bool   use_const   = m_config->get_boolean("hydrology.use_const_bmelt");
  double const_bmelt = m_config->get_double("hydrology.const_bmelt");

  const IceModelVec2S        &bmelt = *m_grid->variables().get_2d_scalar("bmelt");
  const IceModelVec2CellType &mask  = *m_grid->variables().get_2d_cell_type("mask");

  if (not m_hold_bmelt) {
    m_bmelt_local.copy_from(bmelt);
  }

  IceModelVec::AccessList list{&m_bmelt_local, &mask, &result};
  if (m_inputtobed != NULL) {
    m_inputtobed->update(hydro_t, hydro_dt);
    m_inputtobed->interp(hydro_t + hydro_dt/2.0);
    list.add(*m_inputtobed);
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.icy(i, j)) {
      result(i,j) = (use_const) ? const_bmelt : m_bmelt_local(i,j);
      if (m_inputtobed != NULL) {
        result(i,j) += (*m_inputtobed)(i,j);
      }
    } else {
      result(i,j) = 0.0;
    }
  }
}

void Hydrology::passSufaceModelIn(surface::SurfaceModel *m_surface) {
  m_surfaceT = m_surface;


  if(m_surfaceT) {
  m_log->message(2,
             "* Surface model detected in the initialization\n");


  } else {
  m_log->message(2,
             "* Surface model not detected in the initialization\n");
  }
}

} // end of namespace hydrology
} // end of namespace pism
