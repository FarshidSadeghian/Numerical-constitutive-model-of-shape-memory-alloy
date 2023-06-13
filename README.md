# Numerical-constitutive-model-of-shape-memory-alloy
SMA_UM is a FORTRAN-coded thermomechanical constitutive model implementation for any finite element program.
Brinson's original one-dimensional phenomenological constitutive model for shape memory alloys is modified to predict asymmetric behavior in tension and compression.
We propose a method that divides the volume fraction of stress-induced martensite into two portions, one in tension and one in compression.
We implement the proposed model as a two-dimensional Euler–Bernoulli beam element in a user-defined material subroutine of the nonlinear finite element software ABAQUS/Standard.
