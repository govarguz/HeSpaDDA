/*
  Copyright (C) 2012,2013,2017(H)
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "log4espp.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "NodeGrid.hpp"

#include <boost/bind.hpp>
#include "mpi.hpp"

namespace espressopp {
  namespace storage {
    LOG4ESPP_LOGGER(NodeGrid::logger, "DomainDecomposition.NodeGrid");
  
    NodeGridIllegal::NodeGridIllegal()
      : std::invalid_argument("node grid dimensions have to be positive") {}

    NodeGrid::
    NodeGrid(const Int3D& grid,
         const longint nodeId,
         const Real3D& domainSize,
         const boost::python::list& _neiListx,
         const boost::python::list& _neiListy,
         const boost::python::list& _neiListz)
      : Grid(grid)
    {
      if (grid[0] <= 0 || grid[1] <= 0 || grid[2] <= 0) {
        throw NodeGridIllegal();
      }
      for (int i = 0; i < boost::python::len(_neiListx); ++i) {
        neiListx.push_back(boost::python::extract<int>(_neiListx[i]));
      }
      for (int i = 0; i < boost::python::len(_neiListy); ++i) {
        neiListy.push_back(boost::python::extract<int>(_neiListy[i]));
      }
      for (int i = 0; i < boost::python::len(_neiListz); ++i) {
        neiListz.push_back(boost::python::extract<int>(_neiListz[i]));
      }
    maxDomainSizeInCells[0]=neiListx[neiListx.size()-1];
    maxDomainSizeInCells[1]=neiListy[neiListy.size()-1];
    maxDomainSizeInCells[2]=neiListz[neiListz.size()-1];


    mapIndexToPosition(localNodePos, mpiWorld->rank());  // check this  |  mapIndexToPosition(nodePos, node);
    localBoxSize[0] = (static_cast<real>((neiListx[localNodePos[0]+1])-(neiListx[localNodePos[0]]))/static_cast<real>(maxDomainSizeInCells[0]))*domainSize[0];
    invLocalBoxSize[0] = 1.0/localBoxSize[0];

    localBoxSize[1] = (static_cast<real>((neiListy[localNodePos[1]+1])-(neiListy[localNodePos[1]]))/static_cast<real>(maxDomainSizeInCells[1]))*domainSize[1];
    invLocalBoxSize[1] = 1.0/localBoxSize[1];

    localBoxSize[2] = (static_cast<real>((neiListz[localNodePos[2]+1])-(neiListz[localNodePos[2]]))/static_cast<real>(maxDomainSizeInCells[2]))*domainSize[2];
    invLocalBoxSize[2] = 1.0/localBoxSize[2];

  smallestLocalBoxDiameter = std::min(std::min(localBoxSize[0], localBoxSize[1]), localBoxSize[2]);

  calcNodeNeighbors(nodeId);
    }

    longint NodeGrid::
    mapPositionToNodeClipped(const Real3D& pos) const
    {
      Int3D cpos;
      real kCellGridSize;
      kCellGridSize=localBoxSize[0]/(static_cast<real>((neiListx[localNodePos[0]+1])-(neiListx[localNodePos[0]]))/static_cast<real>(maxDomainSizeInCells[0]))/static_cast<real>(maxDomainSizeInCells[0]);

      for (int i = 0; i < 3; ++i) {
          cpos[i] = static_cast< int >(pos[i]/kCellGridSize);
        if (cpos[i] <= 0) {
          cpos[i] = 0;
        }
        else if (cpos[i] >= getGridSize(i)) {
          cpos[i] = getGridSize(i) - 1;
        }
        else {
            if (i==0){
            for (int j=0;j<neiListx.size(); ++j)
            {
            if (cpos[i]>(neiListx[j]) && cpos[i]<=(neiListx[j+1])) {
                cpos[i]=j;
                }
            else{
            }
            }
            }
            else if (i==1) {
            for (int k=0;k<neiListy.size(); ++k)
            {
            if (cpos[i]>(neiListy[k]) && cpos[i]<=(neiListy[k+1])) {
                cpos[i]=k;
                }
            else{
            }
            }
            }
            else{
                for (int l=0;l<neiListz.size(); ++l)
                {
                if (cpos[i]>(neiListz[l]) && cpos[i]<=(neiListz[l+1])) {
                    cpos[i]=l;
                    }
                else{
                }
                }
        }
        }
      }
      return mapPositionToIndex(cpos);
    }

    void NodeGrid::calcNodeNeighbors(longint node)
    {
      Int3D nPos;
  
      mapIndexToPosition(nodePos, node);

      LOG4ESPP_DEBUG(logger, "my position: "
		     << node << " -> "
		     << nodePos[0] << " "
		     << nodePos[1] << " "
		     << nodePos[2]);

      for(int dir = 0; dir < 3; ++dir) {
	for(int j = 0; j < 3; ++j) {
	  nPos[j] = nodePos[j];
	}

	// left neighbor in direction dir
	nPos[dir] = nodePos[dir] - 1;
	if(nPos[dir] < 0) {
	  nPos[dir] += getGridSize(dir);
	}
	nodeNeighbors[2*dir] = mapPositionToIndex(nPos);
	LOG4ESPP_DEBUG(logger, "left neighbor in dir " << dir << ": "
		       << getNodeNeighborIndex(2*dir) << " <-> "
		       << nPos[0] << " "
		       << nPos[1] << " "
		       << nPos[2]);

	// right neighbor in direction dir
	nPos[dir] = nodePos[dir] + 1;
	if(nPos[dir] >= getGridSize(dir)) {
	  nPos[dir] -= getGridSize(dir);
	}
	nodeNeighbors[2*dir + 1] = mapPositionToIndex(nPos);

	LOG4ESPP_DEBUG(logger, "right neighbor in dir " << dir << ": "
		       << getNodeNeighborIndex(2*dir+1) << " <-> "
		       << nPos[0] << " "
		       << nPos[1] << " "
		       << nPos[2]);

	// left or right boundary ?
	boundaries[2*dir] = (nodePos[dir] == 0) ? ToRight : 0;
	boundaries[2*dir + 1] = (nodePos[dir] == getGridSize(dir) - 1) ? ToLeft : 0;

	LOG4ESPP_DEBUG(logger, "boundaries in dir " << dir << ": "
		       << getBoundary(2*dir) << " "
		       << getBoundary(2*dir + 1));
      }
    }

  }
}
