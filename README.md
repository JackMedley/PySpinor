# PySpinor
A python library for spinor calculations.

# Uses
This can be used for numerical computations of spinor calculations, for example, in calculating simple QED and QCD matrix element calculations.

To run one of the example codes you must first include this directory (where README.md is located) to the environment variable 'PYTHONPATH', for example by doing the following in a shell:

export PYTHONPATH=$pwd

The examples can then be run in the usual python way with, for example:

python Examples/gg2uux.py

# Aims
Aims include:
  1. Adding support for symbolic calculations,
  2. Supporting generation of matrix elements from diagrams via interfacing with MadGraph v5,
  3. Adding automated colour factor computations
  4. Adding a simple phase space generator for random phase space point evaluation of matrix elements