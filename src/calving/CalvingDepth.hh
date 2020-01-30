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


/*

This Calving methodology, which is modified from the Thickness Calving routine, is designed to create a threshold that
is a function water depth, rather than a static value (which was set to be 200 m). The rational is that in places
like Hudson Bay, where the water depth is often less than 300 m, it doesn't make sense to calve at 200 m.

- Evan Gowan

*/

#ifndef _PISMCALVINGDEPTH_H_
#define _PISMCALVINGDEPTH_H_

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace calving {

/*! \brief Calving mechanism removing the ice at the shelf front that
  has thickness below a given threshold. */
class CalvingDepth : public Component
{
public:
  CalvingDepth(IceGrid::ConstPtr g);
  virtual ~CalvingDepth();

  virtual void init();
  void update(IceModelVec2CellType &pism_mask, IceModelVec2S &ice_thickness);
  const IceModelVec2S& threshold() const;

protected:
  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;
  IceModelVec2S m_calving_threshold_depth;
  IceModelVec2CellType m_old_mask;
};

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMCALVINGDEPTH_H_ */
