# Documentation for NE-SCGLE libraries

This set of libraries serves the purpose of solving the equations presented by the NE-SCGLE formalism for the description of fluids dynamics undergoing a thermodynamic processes. The libraries has been built upon [GSL](https://www.gnu.org/software/gsl/) and are separated in a general purpose math library, a library for structural properties, and lastly, a library for dynamical properties.



In summary, given a time-dependent Helmholtz free energy of a fluid system, the NE-SCGLE formalism allow us to predict its dynamical and rheological properties. At such, the free energy is actually given by the arrangement of the system's particles, which under thermodynamical equilibrium conditions can actually be stated through the static structure factor of an homogeneous and isotropic system $S(k; \bar n , T...)$. Such structure factor can then be said to be the main input of the NE-SCGLE formalism. Thus, the computation of this quantity is the main aim of the structural library.

On the other hand, the dynamics library handles the computation of the NE-SCGLE set of equations. This is done by differentiating between two conditions: one for stationary systems where no process is given, which is handled by the equilibrium version of the theory (SCGLE); and another for non-equilibrium processes done through the general version of the theory (NE-SCGLE).



# Table of contents

- [Documentation for NE-SCGLE libraries](#documentation-for-ne-scgle-libraries)
- [Table of contents](#table-of-contents)
- [Pre-requisites](#pre-requisites)
- [Basic installation](#basic-installation)
- [Math](#math)
- [Structures](#structures)
  - [Approximations](#approximations)
    - [Hard Sphere (3D)](#hard-sphere-3d)
    - [Hard Sphere + Square Well (3D)](#hard-sphere--square-well-3d)
    - [Hard Sphere + Double Yukawa (3D)](#hard-sphere--double-yukawa-3d)
    - [Hard Sphere + Double Exponential (3D)](#hard-sphere--double-exponential-3d)
    - [Hard Disk (2D)](#hard-disk-2d)
- [Dynamics](#dynamics)
  - [SCGLE](#scgle)
  - [NE-SCGLE](#ne-scgle)
  - [Mean Density Heterogeneities](#mean-density-heterogeneities)

# Pre-requisites

# Basic installation

# Math

# Structures

## Approximations

### Hard Sphere (3D)

### Hard Sphere + Square Well (3D)

### Hard Sphere + Double Yukawa (3D)

### Hard Sphere + Double Exponential (3D)

### Hard Disk (2D)


# Dynamics

## SCGLE

## NE-SCGLE

## Mean Density Heterogeneities