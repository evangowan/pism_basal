/* Copyright (C) 2013, 2014, 2015, 2016, 2017 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */



#include "CalvingDepth.hh"

#include "pism/util/Mask.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_const.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"


#include "pism/util/Vars.hh"
#include "pism/util/iceModelVec2T.hh"


namespace pism {

//! @brief Calving and iceberg removal code.
namespace calving {

CalvingDepth::CalvingDepth(IceGrid::ConstPtr g)
  : Component(g) {

  m_old_mask.create(m_grid, "old_mask", WITH_GHOSTS, 1);

  m_calving_threshold_depth.create(m_grid, "calving_threshold_depth", WITHOUT_GHOSTS);

  m_calving_threshold_depth.set_attrs("diagnostic",
                                "threshold used by the 'calving at threshold for CalvingDepth' calving method",
                                "m",
                                ""); // no standard name

}

CalvingDepth::~CalvingDepth() {
  // empty
}


void CalvingDepth::init() {

  m_log->message(2, "* Initializing the 'calving at a threshold depth' mechanism...\n");


  double max_threshold = m_config->get_double("calving.depth_calving.threshold");

  m_calving_threshold_depth.set(max_threshold);


  m_log->message(2,
                  "  Maximum thickness threshold: %3.3f meters.\n", max_threshold);

}

/**
 * Updates ice cover mask and the ice thickness according to the
 * calving rule removing ice at the shelf front that is thinner than a
 * given threshold.
 *
 * @param[in,out] pism_mask ice cover mask
 * @param[in,out] ice_thickness ice thickness
 * @param[in,out] topg base topography
 *
 * @return 0 on success
 */
void CalvingDepth::update(IceModelVec2CellType &pism_mask,
                                IceModelVec2S &ice_thickness) {

  double max_threshold = m_config->get_double("calving.depth_calving.threshold");
  double fraction_depth = m_config->get_double("calving.depth_calving.fraction_depth");


  // this call fills ghosts of m_old_mask
  m_old_mask.copy_from(pism_mask);

  // need to grap the bed elevation to calculate the threshold
  const IceModelVec2S &bed_elevation = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  IceModelVec::AccessList list{&pism_mask, &ice_thickness, &m_old_mask, &m_calving_threshold_depth};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_calving_threshold_depth(i,j) = std::max(-bed_elevation(i,j)*fraction_depth,0.0);
    m_calving_threshold_depth(i,j) = std::min(m_calving_threshold_depth(i,j),max_threshold);

    if (m_old_mask.floating_ice(i, j)           &&
        m_old_mask.next_to_ice_free_ocean(i, j) &&
        ice_thickness(i, j) < m_calving_threshold_depth(i, j)) {
      ice_thickness(i, j) = 0.0;
      pism_mask(i, j)     = MASK_ICE_FREE_OCEAN;
    }
  }

  pism_mask.update_ghosts();
  ice_thickness.update_ghosts();
}

const IceModelVec2S& CalvingDepth::threshold() const {
  return m_calving_threshold_depth;
}

std::map<std::string, Diagnostic::Ptr> CalvingDepth::diagnostics_impl() const {
  return {{"calving_threshold_depth", Diagnostic::wrap(m_calving_threshold_depth)}};
}

} // end of namespace calving
} // end of namespace pism
