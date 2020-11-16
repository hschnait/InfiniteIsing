#pragma once

#include "Lattice.hpp"
#include <set>
#include <algorithm>
#include <random>
#include <vector>

// template <typename return_t>
class WolffSolver {
public:
	typedef InfiniteLattice::indx_t indx_t ;
	typedef InfiniteLattice lattice_t;
	typedef std::vector<indx_t> cluster_t;
	typedef void (* onUpdate_t) (const WolffSolver*, const cluster_t&, const size_t);

	// typedef return_t (* Meassure_t) (const lattice_t&, const cluster_t&, const size_t);




	WolffSolver()  = delete;
	WolffSolver( InfiniteLattice& latt, double beta)
		: _latt { latt }, _beta{beta} {
		}

	cluster_t ClusterSearch( ) const {
		// TODO: Check & optimize random number generation
		std::random_device rd;
		std::mt19937 mt(rd());
		std::uniform_real_distribution<> RNG(0., 1.);


		cluster_t cluster { 0 };
		size_t currentIdx = 0;

		// Bonds are inserted with std::minmax(site1, site2)
		std::set< std::pair<indx_t, indx_t> > bonds_visited {};

		while ( currentIdx < cluster.size() ) {
			const auto current = cluster[currentIdx];
			for ( const auto neighbor : _latt.getNeighbors(current)) {
				if( _latt.getSpin(current) == _latt.getSpin(neighbor)) {
					
					// if neighbor not in cluster
					if(std::find(cluster.begin(), cluster.end(), neighbor) == cluster.end() ) {
						const auto bond = std::minmax(current, neighbor);
						// if bond not in bonds_visited
						if ( bonds_visited.find( bond ) == bonds_visited.end() ) {
							bonds_visited.insert( bond );
							if ( RNG(mt) < (1-exp(-2.*_beta)) ) {
								cluster.push_back(neighbor);
							}
						}
					}
				}
			}
			currentIdx++;
		}

		return cluster;
	}

	// return_t solve(size_t n_iterations) {
	std::vector<double> solve(size_t n_iterations) {
		n_iter = n_iterations;
		std::vector<double> G {};
		std::vector<size_t> G_count { 0 };

		for ( size_t i = 0; i < n_iterations; i++ ) {
			const auto cluster = ClusterSearch();
			_latt.flipSpins(cluster);

			for ( auto site : cluster ) {
				while(G_count.size() < site ) G_count.push_back(0);
				G_count[site]++;
			}

			for ( auto f : onUpdate )
				f(this, cluster, i);
		}
		n_iter = 0;

		G.reserve(G_count.size());
		for ( size_t i = 0; i < G_count.size(); i++)
			G.push_back(static_cast<double>(G_count[i]) / n_iterations);
		return G;
	}

	const lattice_t getLattice() const { return _latt; }
	const size_t getNIterations() const { return n_iter; }
	 
	std::vector<onUpdate_t> onUpdate {};
	// Meassure_t meassure = [](auto,auto,auto) { return_t tmp; return tmp;};
private: 

	size_t n_iter;
	double _beta;
	lattice_t& _latt;
};