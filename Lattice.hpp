#pragma once
/* Lattice.hpp
 * by H. Schnait, 16.11.2020
 * 
 * Infinite (square) lattice for ising model simulations
 *
 * As the lattice is growing once the outer bound is reached (in the non-const getSpin() function),
 * the size of the lattice is stored as the "generation" of the lattice.
 * A lattice of "generation 0" is just one spin (of starting value "SPIN0").
 * A lattice of "generation 1" looks like this:
 * | + |
 * |+-+|
 * | + |
 * A lattice of "generation 2" looks like this:
 * |  -  |
 * | -+- |
 * |-+-+-|
 * | -+- |
 * |  -  |
 * and so on...
 * 
 * There are two types of site-indexing:
 * coord_t :
 * 	Carthesian coordinates (both in negative and positive direction) with the first spin being (0, 0).
 * 	Used for neighbor-detection. 1st number is x-direction (right being +), 2nd coordinate is y-direction (up being +)
 * indx_t :
 * 	Linear indexing used for array addressing, defined in the following way (for a lattice of generation 2):
 * 	Hex-numbers used in this comment for easier portrayal
 * 		|  5  |
 * 		| c16 |
 * 		|b4027|
 * 		| a38 |
 * 		|  9  |
 * coord_t and indx_t can be transformed into one another through coord2indx() and indx2coord()
*/






#include <vector>
#include <cstddef>
#include <iostream>
#include <cmath>
#include <array>

#define SPIN0 true

class InfiniteLattice {
public:
	typedef size_t indx_t;
	typedef std::array<indx_t, 4> neighbor_t;
	typedef std::pair<int,int> coord_t;
	

	InfiniteLattice() {	}
	InfiniteLattice(size_t initial_generaton) {
		for (size_t i = 0; i < initial_generaton; i++)
			IncreaseLatticeSize();
	}

	indx_t coord2indx(const coord_t c ) const {
		const size_t gen = abs(c.first) + abs(c.second);

		if ( c.first >= 0 && c.second > 0 ) {
			// Sector I
			const size_t pos_in_sector = (gen + c.first - c.second) / 2;
			return NumPointsUpToGen(gen-1) + pos_in_sector;
		}
		else if ( c.first > 0 && c.second <= 0 ) {
			// Sector II
			const size_t pos_in_sector = (gen - c.first - c.second) / 2;
			return NumPointsUpToGen(gen-1) + pos_in_sector + gen;
		}
		else if ( c.first <= 0 && c.second < 0 ) {
			// Sector III
			const size_t pos_in_sector = (gen - c.first + c.second) / 2;
			return NumPointsUpToGen(gen-1) + pos_in_sector + 2*gen;
		}
		else if ( c.first < 0 && c.second >= 0) {
			// Sector IV
			const size_t pos_in_sector = (gen + c.first + c.second) / 2;
			return NumPointsUpToGen(gen-1) + pos_in_sector + 3*gen;
		}
		else {
			// Point(0,0)
			return 0;
		}
	}
	
	coord_t indx2coord(const indx_t i ) const {
		if ( i == 0 )
			return {0, 0};

		const size_t gen = getGeneration(i);
		const size_t pos_in_gen = i - NumPointsUpToGen(gen-1);
		const size_t pos_in_sector = pos_in_gen % gen;

		switch (pos_in_gen / gen) {
			case 0: // Sector I
				return { pos_in_sector, gen-pos_in_sector };
			case 1: // Sector II
				return { gen-pos_in_sector, -pos_in_sector }; 
			case 2: // Sector III
				return { -pos_in_sector, -gen+pos_in_sector };
			case 3: // Sector IV
				return { -gen+pos_in_sector, pos_in_sector };
			default:
				// should not happen
				throw 0;
		}
		return {0,0};
	}

	size_t getGeneration(indx_t site) const { return (site==0) ? 0 : (0.5*((std::sqrt(2*site - 1)) - 1))+1;}
	size_t getGeneration() const { return _gen; }
	std::vector<neighbor_t> getNeighbors() const { return _neighbors; }
	neighbor_t getNeighbors(indx_t site) const { return _neighbors[site]; }

	/// getSpin is non-const, as it can increase the lattice size!
	bool getSpin(indx_t site) { 
		if ( site >= _spins.size() ) {
			IncreaseLatticeSize();
		}
		return _spins[site]; }
	bool getSpin(indx_t site) const {
		return _spins[site];	}


	void flipSpins(std::vector<indx_t> cluster) {
		for ( auto site : cluster ) {
			_spins[site] = !_spins[site];
		}
	}

	friend std::ostream& operator<<(std::ostream&, const InfiniteLattice&);
	friend std::ostream& operator<<(std::ostream&, const std::vector<neighbor_t>&);
	friend std::ostream& operator<<(std::ostream&, const neighbor_t&);

private:
	void IncreaseLatticeSize() {
		++_gen;

		_spins.reserve(NumPointsUpToGen(_gen));
		_neighbors.reserve(NumPointsUpToGen(_gen));

		
		const bool neelSpin { (_gen % 2 == 0 ) ? SPIN0 : !SPIN0 };
		for ( indx_t i = 0; i < _gen*4; ++i ) {
			_neighbors.push_back( get_neighbors(_spins.size()) );
			_spins.push_back(neelSpin);

		}
	}

	neighbor_t get_neighbors( indx_t site ) const {
		const auto site_coords = indx2coord(site);
		neighbor_t neighbors{
			coord2indx({site_coords.first, site_coords.second+1}),
			coord2indx({site_coords.first+1, site_coords.second}),
			coord2indx({site_coords.first, site_coords.second-1}),
			coord2indx({site_coords.first-1, site_coords.second})
		};

		return neighbors;
	}
	
	/// Number of total lattice sites needed for gen g
	size_t NumPointsUpToGen(const size_t g) const { return 2*g*g + 2*g + 1; }
	

	size_t _gen {0};

	std::vector<bool> _spins { SPIN0 };
	std::vector<neighbor_t> _neighbors { {1,2,3,4} };



	// For debugging only
	// void test_c2i(coord_t c) {
	// 	std::cout << "Point(" << c.first << "," << c.second << ") is equal to: " << coord2indx(c) << '\n';
	// }
	// void test_i2c(size_t i) {
	// 	std::cout << "Indx" <<  i << " is equal to: Point(" << indx2coord(i).first << "," << indx2coord(i).second <<  ")" << '\n';
	// }
};

#define SPIN1_CHAR "▓"
#define SPIN2_CHAR "░"
#define EMPTY_CHAR " "
std::ostream& operator<<(std::ostream& os, const InfiniteLattice& latt) {
	os << "Infinite lattice of generation " << latt._gen << ":\n";
	os << '\n';

	for ( size_t row = 0; row < latt._gen*2 + 1; row++ ) {
		for ( size_t col = 0; col < latt._gen*2 + 1; col++ ) {
			if ( col < abs(latt._gen - row) || col > latt._gen*2 - abs(latt._gen - row) ) {
				os << EMPTY_CHAR;
			}
			else {
				const bool spin = latt._spins[latt.coord2indx({col-latt._gen,row-latt._gen})];
				os << ( (spin == true) ? SPIN1_CHAR : SPIN2_CHAR);
				//os << ( (spin == true) ? '+' : '-');
			}
		}

		os << '\n';
	}
	os << '\n';

	
	
	
	return os;
}
std::ostream& operator<<(std::ostream& os, const InfiniteLattice::neighbor_t& neighbors) {
	for (auto neighbor : neighbors ) {
		os << neighbor << ", ";
	}

	return os;
}
std::ostream& operator<<(std::ostream& os, const std::vector<InfiniteLattice::neighbor_t>& neighbors_list) {
	os << "InfiniteLattice Neighbors:\n";

	for ( size_t i = 0; i < neighbors_list.size(); i++ ) {
		auto neighbors = neighbors_list[i];
		os << "SiteNo " << i << " neighbors: ";
		os << neighbors << '\n';
	}
	os << '\n';
	
	return os;
}

