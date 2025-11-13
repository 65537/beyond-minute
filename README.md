# Affine Deligne-Lusztig varieties beyond the minute case
This is the software repository associated with the preprint <https://arxiv.org/abs/2511.08879> from Felix Schremmer and Eva Viehmann. The software provides mathematical calculations associated to admissible sets and affine Deligne-Lusztig varieties as described in our paper. This does not only serve as a tool to verify some of our claims made, but may also assist in identifying new cases worth studying.
## Installation
The requirement is an up to date installation of the [Sagemath](https://www.sagemath.org/) computer algebra system. Then, you may clone this github repository or manually download the program files into a folder of your choice.

## Usage
Call the Sagemath interpreter on one of the three program files, which are called `coweightsByDepth.py`, `geometricPositiveCoxeterType.py` and `L1BC_ING.py` respectively. The last file, named `sharedFunctions.py` should not be called directly. The first command line argument should always be the description of a connected Dynkin diagram, such as G2 or D6. You can call the files without any arguments to get an overview of the remaining arguments.

### Enumeration of coweights of depth <2

In the proof of Theorem 2.7, we rely on exhaustive computations of coweights of depth $<2$ for all root systems of rank $<20$. In order to verify our claims, one may make the following calls.

```
sage coweightsByDepth.py A1 2
sage coweightsByDepth.py A2 2
sage coweightsByDepth.py B2 2
sage coweightsByDepth.py G2 2

...

sage coweightsByDepth.py B4 2
sage coweightsByDepth.py C4 2
sage coweightsByDepth.py D4 2
sage coweightsByDepth.py F4 2

...
```

### Verification of geometric and positive Coxeter type

In the proof of Theorem 3.11 and Proposition 3.22, we claim that a number of admissible sets have to be checked individually to verify geometric resp. positive Coxeter type. For this, one may make the following calls.

```
sage geometricPositiveCoxeterType.py A1 omega[1],omega[1],omega[1]
sage geometricPositiveCoxeterType.py A4 omega[1],omega[2]
sage geometricPositiveCoxeterType.py A5 omega[3]
sage geometricPositiveCoxeterType.py A6 omega[3]
sage geometricPositiveCoxeterType.py A7 omega[3]
sage geometricPositiveCoxeterType.py C3 omega[3]
sage geometricPositiveCoxeterType.py D5 omega[5]
```

If one were to check, e.g., these properties for the entire admissible set associated to a minuscule cocharacter of a non-split group of type A2, the corresponding call would be as follows.
```
sage geometricPositiveCoxeterType.py A2 omega[1] --level iwahori --frobenius quasi-split
```

### Verification of (ING) and (L1BC)

In the proof of Theorem 4.8, we rely on checking the property (L1BC) in a finite number of exceptional cases. Explicitly, one may make the following calls.
```
sage L1BC_ING.py A4 omega[1]+omega[2]
sage L1BC_ING.py A5 omega[3]
sage L1BC_ING.py A6 omega[3]
sage L1BC_ING.py A7 omega[3]
```
The first such call also verifies the property (ING) and computes the maximal length as needed in the proof of Theorem 5.30.
