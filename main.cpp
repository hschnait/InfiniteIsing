#include <iostream>
#include "Lattice.hpp"
#include "Wolff.hpp"
#include <iterator>
#include <fstream>
#include <unistd.h> // for usleep


/* For Color-CLI plot */
#define SPIN1_CHAR "▓"
#define SPIN2_CHAR "░"
#define CLUSTER_CHAR "\033[31m▒\033[0m" // red
#define EMPTY_CHAR " "
#define SLEEP_TIME 50000 // mu-seconds



void naive_cli(const WolffSolver* solver, const WolffSolver::cluster_t& cluster, const size_t i);
void color_cli(const WolffSolver* solver, const WolffSolver::cluster_t& cluster, const size_t i);
void print_percentage(const WolffSolver* l, const WolffSolver::cluster_t& c, const size_t i);

int main (int argc, char* argv[]) {
	if ( argc < 2 ) {
		std::cerr << "Usage: " << argv[0] << " <BETA>" << "\n";
		return 1;
	}

	const int n_iterations = 100000000;
	const double beta = atof(argv[1]);

	InfiniteLattice myLattice{0};

	WolffSolver mySolver {myLattice, beta}; //? beta_crit = 0.44068

	// UNCOMMENT this for a percentage-bar
	// mySolver.onUpdate.push_back(print_percentage);

	// UNCOMMENT this for simple b/w cli-plot (without sleep in between iterations)
	// mySolver.onUpdate.push_back(naive_cli);

	// UNCOMMENT this for color cli-plot with sleep between iterations
	mySolver.onUpdate.push_back(color_cli);


	auto G = mySolver.solve(n_iterations);

	std::cout << "Final lattice size: " << myLattice.getGeneration() << "\n";
	

	std::ofstream output_file("G_"+std::to_string(beta)+ ".dat");
	for ( size_t i = 0; i < G.size(); i++ ) {
		const auto c = myLattice.indx2coord(i);
		output_file << c.first << "\t" << c.second << "\t" << G[i] << "\n";
	}


	return 0;
}

void naive_cli(const WolffSolver* solver, const WolffSolver::cluster_t& cluster, const size_t i) {
	if ( i % 1 == 0 ) {
		for ( int j = 0; j < 2*solver->getLattice().getGeneration() + 4; j++) {
			std::cout << "\r\e[A";
		}
		std::cout << solver->getLattice();
		std::cout.flush();
	}

	// return 0.0;
}

void color_cli(const WolffSolver* solver, const WolffSolver::cluster_t& cluster, const size_t i) {
	std::cout << '\n';
	const auto& latt = solver->getLattice();
	const auto gen = latt.getGeneration();

	// Move cursor up
	for ( int j = 0; j < 2*latt.getGeneration() + 3; j++) {
		std::cout << "\r\e[A";
	}

	// Print lattice
	for ( size_t row = 0; row < gen*2 + 1; row++ ) {
		for ( size_t col = 0; col < gen*2 + 1; col++ ) {
			if ( col < abs(gen - row) || col > gen*2 - abs(gen - row) ) {
				std::cout << EMPTY_CHAR;
			}
			else {
				const auto indx = latt.coord2indx({col-gen,row-gen});
				// if not in cluster
				if( std::find(cluster.begin(), cluster.end(), indx) == cluster.end()) {
					const bool spin = latt.getSpin(indx);
					std::cout << ( (spin == true) ? SPIN1_CHAR : SPIN2_CHAR);
				}
				else {
					std::cout << CLUSTER_CHAR;
				}
			}
		}

		std::cout << '\n';
	}

	std::cout << '\n';
	std::cout.flush();
	usleep(SLEEP_TIME);
}

void print_percentage(const WolffSolver* l, const WolffSolver::cluster_t& c, const size_t i) {
	if ( i %( l->getNIterations() / 100) == 0 ) {
		std::cout << (i*100)/l->getNIterations() << "%\n";
	}
}