// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#ifndef _PISMMOHRCOULOMBYIELDSTRESSEVAN_H_
#define _PISMMOHRCOULOMBYIELDSTRESSEVAN_H_

#include "YieldStress.hh"
#include "MohrCoulombYieldStress.hh"

#include "pism/util/iceModelVec.hh"

namespace pism {

class IceModelVec2CellType;

namespace hydrology {
class Hydrology;
}

//! @brief PISM's default basal yield stress model which applies the
//! Mohr-Coulomb model of deformable, pressurized till.
class MohrCoulombYieldStressEvan : public MohrCoulombYieldStress {
public:
  MohrCoulombYieldStressEvan(IceGrid::ConstPtr g, hydrology::Hydrology *hydro);
  virtual ~MohrCoulombYieldStressEvan();


protected:
  virtual void init_impl();
  virtual void update_impl(const YieldStressInputs &inputs);

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

protected:
//  IceModelVec2S m_till_phi, m_tillwat, m_Po;
  hydrology::Hydrology *m_hydrology;

  IceModelVec2S m_effective_pressure, m_till_cover_local, m_sliding_mechanism, m_velocity_temp, hydro_tauc;
};

} // end of namespace pism

#endif /* _PISMMOHRCOULOMBYIELDSTRESSEVAN_H_ */
