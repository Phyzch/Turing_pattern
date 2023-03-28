# Turing_pattern

Code for doing Turing pattern simulation.
This code do simulation of Turing pattern for pattern which has fixed position. 
Contact cz38@illinois.edu when you have any question.

## How to use the program

1. Link the project. Go to the project folder (use Cmake)

```bash
mkdir ./Release
cd ./Release
cmake .. -DCMAKE_BUILD_TYPE="Release"
make
```

2. example input file is given in ./input folder

## Structure of the program

### Class in the program

There are two classes : 

1. **Brusselator_reaction** : reaction in single box

2. **Brusselator_reaction_system** : reaction for 1d system. It has number of boxes  (**Brusselator_reaction** instance) where local reaction can happen (assumed well-mixed in each box). Boxes are connected by linked list: (pointer variable **previous** , **next** ) and particle $X_{1} , X_{2}, X_{3}$ can diffuse between boxes.

### Gillespie algorithm

This code using Gillespie algorithm to simulate extended 1d system.  

(Gillespie, Daniel T. (1977). "Exact Stochastic Simulation of Coupled Chemical Reactions". *The Journal of Physical Chemistry*. **81** (25): 2340–2361.)

The number of species is specified by $N$ and number of chemical reaction (including diffusion) is specified by $M$.

The initial concentration of particle X1_0, X2_0 and X3_0 is set in program as input, this is concentration of homogeneous state.

Parameter for simulation is based on choice made in [Self-organization and positioning of bacterial protein clusters | Nature Physics](https://www.nature.com/articles/nphys4155). See SI of above paper for reference.

**Speed up calculation** : Notice that in extended system, when diffusion or reaction happens, the chemical concentration to box not connected to the box that reaction happens will not change reaction rate at all. Therefore, when update reaction rate for the extended system, we only need to re-calculate reaction rate related to given box. This idea is implemented in the function: 

**Brusselator_reaction_system : System_evolve()**

The reaction constant is set as  Brusselator_reaction::c

The number of particles in given box is set as : Brusselator_reaction::S 

The reaction rate for each local reaction and diffusion is set as : Brusselator_reaction::a 
