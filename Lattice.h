/*
 * Lattice.h
 *
 *  Created on: 02.02.2016
 *      Author: dschmidt
 */

#ifndef LATTICE_H_
#define LATTICE_H_

/* Representation of the lattice, assuming the same number of points in every direction */

#include <cmath>
#include <cstddef>
#include <vector>

class Lattice {
public:
  Lattice(const size_t lenTime, const size_t lenSpace, const size_t numDim);
  ~Lattice();

  inline const size_t getDim() const;
  inline const size_t getTimeSize() const;
  inline const size_t getSpaceSize() const;
  inline const size_t getVol() const;

  inline const std::vector< size_t > getNeighbours( const size_t x ) const;
  inline const std::vector< std::vector<size_t> > & getLineIndex() const;


private:

  void makeNeighbourIndex();
  // creates a list of every point belonging to a line in time direction for every spatial point in the lattice
  void makeLineIndex();

  const size_t Nt;    // lattice size in time direction
  const size_t Ns;    // lattice size in space direction
  const size_t dim;   // number of dimensions
  const size_t V;     // total amount of lattice points
  std::vector< std::vector<size_t> > nnIndex;  // contains a list of nearest neighbours for every lattice point
  std::vector< std::vector<size_t> > lineIndex;
};

/* =============================================================================
 * implementation of inline functions
 * ============================================================================= */

const size_t Lattice::getDim() const {
  return dim;
}

const size_t Lattice::getSpaceSize() const {
	return Ns;
}

const size_t Lattice::getTimeSize() const {
	return Nt;
}

const size_t Lattice::getVol() const {
  return V;
}

const std::vector< size_t > Lattice::getNeighbours( const size_t x ) const {
	return nnIndex[x];
}

const std::vector< std::vector<size_t> > & Lattice::getLineIndex() const {
  return lineIndex;
}

#endif /* LATTICE_H_ */
