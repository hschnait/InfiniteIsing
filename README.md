# Ising model on infinite Lattice

For the Phase Transitions lecture by Prof. H.G. Evertz (WS2020).

## Compilation
Compile with 

    g++ main.cpp -o Ising -O3

and then run the code with

    ./Ising <BETA>


## Changes
The output of the program can be set in `main.cpp` by commenting in/out the `mySolver.onUpdate.push_back(...)` lines.
By default a color output with sleep in between iterations is set. For measuring G this should be changed of course...

The command line output works on my Linux machine, I've not tested it anywhere else...

