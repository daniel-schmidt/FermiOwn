//============================================================================
// Name        : FermiOwn.cpp
// Author      : Daniel Schmidt
// Version     :
// Copyright   : 
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>		// std::cout
#include <vector>		// std::equal
#include <algorithm>	// std::vector

#include "Lattice.h"

namespace FermiOwnTest {
	using namespace FermiOwn;

bool test_neighbours( const size_t x, const Lattice & lat, const std::vector<size_t>& solution ) {
	std::vector<size_t> testNeighbours = lat.getNeighbours( x );

	if( testNeighbours.size() == solution.size() ) {
		if( std::equal( solution.begin(), solution.end(), testNeighbours.begin() ) ) {
			std::cout << "Test passed." << std::endl;
			return true;
		} else {
			std::cout << "Content is not the same." << std::endl;
			return false;
		}
	} else {
		std::cout << "Size of neighbour list does not match solution." << std::endl;
		return false;
	}
}

bool test_neighbours_dir( const size_t x, const Lattice & lat, const std::vector<size_t>& solution, Lattice::Dir dir ) {
	std::vector<size_t> testNeighbours = lat.getNeighbours( x, dir );

	if( testNeighbours.size() == solution.size() ) {
		if( std::equal( solution.begin(), solution.end(), testNeighbours.begin() ) ) {
			std::cout << "Test passed." << std::endl;
			return true;
		} else {
			std::cout << "Content is not the same." << std::endl;
			return false;
		}
	} else {
		std::cout << "Size of neighbour list does not match solution. List is: " << std::endl;
		for( size_t x : testNeighbours ) {
			std::cout << x << std::endl;
		}
		return false;
	}
}

} // namespace FermiOwnTest

int main(void) {
	using namespace FermiOwn;
	using namespace FermiOwnTest;
	std::cout << "Testing a lattice object in 2 dimensions." << std::endl;

	size_t Nt = 2;
	size_t Ns = 3;
	size_t dim = 2;
	Lattice lat2d( Nt, Ns, dim);

	if( lat2d.getDim() != dim )
		std::cout << "Lattice returns wrong dimension." << std::endl;
	if( lat2d.getVol() != Nt*Ns)
		std::cout << "Lattice returns wrong volume." << std::endl;
	if( lat2d.getSpaceSize() != Ns)
		std::cout << "Lattice returns wrong spatial size." << std::endl;
	if( lat2d.getTimeSize() != Nt)
		std::cout << "Lattice returns wrong time size." << std::endl;

	size_t x = 0;
	std::vector<size_t> solution( {1,1,4,2} );
	test_neighbours(x, lat2d, solution);

	x = 1;
	solution = {0, 0, 5, 3};
	test_neighbours(x, lat2d, solution);

	x = 2;
	solution = {3, 3, 0, 4};
	test_neighbours(x, lat2d, solution);

	x = 5;
	solution = {4, 4, 3, 1};
	test_neighbours(x, lat2d, solution);

	std::cout << "Testing a lattice object in 3 dimensions." << std::endl;

	Nt = 3;
	Ns = 3;
	dim = 3;
	Lattice lat3d( Nt, Ns, dim);

	if( lat3d.getDim() != dim )
		std::cout << "Lattice returns wrong dimension." << std::endl;
	if( lat3d.getVol() != Nt*Ns*Ns)
		std::cout << "Lattice returns wrong volume." << std::endl;
	if( lat3d.getSpaceSize() != Ns)
		std::cout << "Lattice returns wrong spatial size." << std::endl;
	if( lat3d.getTimeSize() != Nt)
		std::cout << "Lattice returns wrong time size." << std::endl;

	x = 0;
	solution = {2,1,6,3,18,9};
	test_neighbours(x, lat3d, solution);

	x = 13;
	solution = {12,14,10,16,4,22};
	test_neighbours(x, lat3d, solution);

	x = 25;
	solution = {24,26,22,19,16,7};
	test_neighbours(x, lat3d, solution);

	solution = {24,22,16};
	test_neighbours_dir(x, lat3d, solution, Lattice::bwd);
	solution = {26,19,7};
	test_neighbours_dir(x, lat3d, solution, Lattice::fwd);

	return EXIT_SUCCESS;
}

