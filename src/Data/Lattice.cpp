/*
 * Lattice.cpp
 *
 *  Created on: 02.02.2016
 *      Author: dschmidt
 */

#include "../Data/Lattice.h"

#include <iostream>

namespace FermiOwn {

Lattice::Lattice(size_t lenTime, size_t lenSpace, size_t numDim) :
Nt(lenTime),
Ns(lenSpace),
dim(numDim),
V( Nt*pow(Ns, dim-1) )
{
  makeNeighbourIndex();
  makeLineIndex();
}

Lattice::~Lattice(){};

/* =============================================================================
 * private functions
 * ============================================================================= */
void Lattice::makeNeighbourIndex() {
  nnIndex = std::vector< std::vector<size_t> >(V);
  for( size_t x = 0; x < V; x++ ) {
    for( size_t mu = 0; mu < dim ; mu++ ) {

      int dx;
      if( mu == 0 ) {
        dx = 1;
      } else {
        dx = Nt*pow(Ns,mu-1);
      }
      int blockLen = Nt * pow(Ns, mu);

      int nnM = x-dx;
      int nnP = x+dx;
      int line = x/blockLen;
      if( nnM < 0 || nnM/blockLen != line )
        nnM += blockLen;
      if( nnP/blockLen != line )
        nnP -= blockLen;
      nnIndex[x].push_back(nnM);
      nnIndex[x].push_back(nnP);
    }
  }
  #ifdef DEBUG
  std::cout << "NeighbourIndex: " << std::endl;
  for( auto point : nnIndex) {
    for( size_t nn : point ) {
      std::cout << nn << " ";
    }
    std::cout << std::endl;
  }
  #endif
}

void Lattice::makeLineIndex() {
  if( nnIndex.empty() )
  {
    std::cout << "Lattice cannot create LineIndex, because NeighbourIndex is empty. You have to create the NeighbourIndex first!" << std::endl;
    exit(1);
  }

  lineIndex = std::vector< std::vector<size_t> >(pow(Ns,dim-1));

  for( size_t x = 0; x < V; x++ ) {
    lineIndex[x/Nt].push_back(x);
  }

#ifdef DEBUG
  std::cout << "LineIndex: " << std::endl;
  for( auto line : lineIndex ) {
    for( auto point : line ) {
      std::cout << point << " ";
    }
    std::cout << std::endl;
  }
#endif
}

} // namespace FermiOwn
