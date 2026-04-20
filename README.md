CLASS: Cosmic Linear Anisotropy Solving System  {#mainpage}
==============================================

 - Edit version of CLASS : [CLASS - Kerr-Cartan Cosmological Model](https://github.com/FlavioGomes-hub/class_kerr_cartan/blob/main/README.md#class---kerr-cartan-cosmological-model)

Authors: Julien Lesgourgues, Thomas Tram, Nils Schoeneberg

with several major inputs from other people, especially Benjamin
Audren, Simon Prunet, Jesus Torrado, Miguel Zumalacarregui, Francesco
Montanari, Deanna Hooper, Samuel Brieden, Daniel Meinert, Matteo Lucca, etc.

For download and information, see http://class-code.net


Compiling CLASS and getting started
-----------------------------------

(the information below can also be found on the webpage, just below
the download button)

Download the code from the webpage and unpack the archive (tar -zxvf
class_vx.y.z.tar.gz), or clone it from
https://github.com/lesgourg/class_public. Go to the class directory
(cd class/ or class_public/ or class_vx.y.z/) and compile (make clean;
make class). You can usually speed up compilation with the option -j:
make -j class. If the first compilation attempt fails, you may need to
open the Makefile and adapt the name of the compiler (default: gcc),
of the optimization flag (default: -O4 -ffast-math) and of the OpenMP
flag (default: -fopenmp; this flag is facultative, you are free to
compile without OpenMP if you don't want parallel execution; note that
you need the version 4.2 or higher of gcc to be able to compile with
-fopenmp). Many more details on the CLASS compilation are given on the
wiki page

https://github.com/lesgourg/class_public/wiki/Installation

(in particular, for compiling on Mac >= 10.9 despite of the clang
incompatibility with OpenMP).

To check that the code runs, type:

    ./class explanatory.ini

The explanatory.ini file is THE reference input file, containing and
explaining the use of all possible input parameters. We recommend to
read it, to keep it unchanged (for future reference), and to create
for your own purposes some shorter input files, containing only the
input lines which are useful for you. Input files must have a *.ini
extension. We provide an example of an input file containing a
selection of the most used parameters, default.ini, that you may use as a
starting point.

If you want to play with the precision/speed of the code, you can use
one of the provided precision files (e.g. cl_permille.pre) or modify
one of them, and run with two input files, for instance:

    ./class test.ini cl_permille.pre

The files *.pre are suppposed to specify the precision parameters for
which you don't want to keep default values. If you find it more
convenient, you can pass these precision parameter values in your *.ini
file instead of an additional *.pre file.

The automatically-generated documentation is located in

    doc/manual/html/index.html
    doc/manual/CLASS_manual.pdf

On top of that, if you wish to modify the code, you will find lots of
comments directly in the files.

Python
------

To use CLASS from python, or ipython notebooks, or from the Monte
Python parameter extraction code, you need to compile not only the
code, but also its python wrapper. This can be done by typing just
'make' instead of 'make class' (or for speeding up: 'make -j'). More
details on the wrapper and its compilation are found on the wiki page

https://github.com/lesgourg/class_public/wiki

Plotting utility
----------------

Since version 2.3, the package includes an improved plotting script
called CPU.py (Class Plotting Utility), written by Benjamin Audren and
Jesus Torrado. It can plot the Cl's, the P(k) or any other CLASS
output, for one or several models, as well as their ratio or percentage
difference. The syntax and list of available options is obtained by
typing 'pyhton CPU.py -h'. There is a similar script for MATLAB,
written by Thomas Tram. To use it, once in MATLAB, type 'help
plot_CLASS_output.m'

Developing the code
--------------------

If you want to develop the code, we suggest that you download it from
the github webpage

https://github.com/lesgourg/class_public

rather than from class-code.net. Then you will enjoy all the feature
of git repositories. You can even develop your own branch and get it
merged to the public distribution. For related instructions, check

https://github.com/lesgourg/class_public/wiki/Public-Contributing

Using the code
--------------

You can use CLASS freely, provided that in your publications, you cite
at least the paper `CLASS II: Approximation schemes <http://arxiv.org/abs/1104.2933>`. Feel free to cite more CLASS papers!

Support
-------

To get support, please open a new issue on the

https://github.com/lesgourg/class_public 

or me 

f.gomes.phys@proton.me
or
Reddit : u/CosmoFlavio

***

# CLASS - Kerr-Cartan Cosmological Model

This modified version of the **CLASS (Cosmological Linear Anisotropy Solvers System)** code implements the geometric corrections derived from the **Kerr-Cartan** metric. This model resolves the Hubble tension through a pure spin-torsion coupling factor, without the addition of scalar fields or exotic dark energy.

## Theoretical Foundations

The implementation is based on the assumption of a universe embedded in a quasi-extreme Kerr geometry ($a_{*} = 0.998$). Within the Einstein-Cartan-Sciama-Kibble (ECSK) formalism, the isotropization of the Kantowski-Sachs anisotropy isolates an invariant geometric factor:
$$\Gamma_{KC} = \frac{13}{12} \approx 1.0833$$

This factor alters the dynamics of Friedmann expansion by introducing a topological correction to the Hubble rate $H(z)$.

## Changes to the Source Code

Changes have been made to the cosmological background module (`source/background.c`) to ensure that the new metric is propagated consistently throughout all parts of the code (distances, proper time, and perturbations).

### Technical details (`source/background.c`):
The expansion rate $H$ and its derivative with respect to proper time have been scaled by the factor $\Gamma_{KC}$:
- **Standard:** $H = \sqrt{\rho_{tot}}$
- **Kerr-Cartan:** `pvecback[pba->index_bg_H] = sqrt(rho_tot - pba->K/a/a) * GAMMA_KC;`
- **Derivative:** `pvecback[pba->index_bg_H_prime] = (- (3./2.) * (rho_tot + p_tot) * a + pba->K/a) * GAMMA_KC;`

## Results and Numerical Validation (MCMC)

The model was validated using a Monte Carlo Markov Chain (MCMC) analysis with the **MontePython** sampler.

### "Gold" Run Configuration:
- **Data:** Pantheon+ (Type Ia supernovae).
- **Prior:** Absolute magnitude $M$ fixed to the SH0ES local calibration ($M = -19.24 \pm 0.04$).
- **Convergence:** Gelman-Rubin criterion $R-1 = 0.00024$ (full stability over 10,000 steps).

### Resolution of the Hubble Tension:
While the standard $\Lambda$CDM model fails to reconcile local and primordial measurements, the Kerr-Cartan model naturally converges to:
**$H_0 = 67.81 \pm 0.23$ km/s/Mpc**

This result is in perfect statistical agreement ($1\sigma$) with the Planck 2018 constraints ($H_0 \approx 67.4$), proving that the Hubble tension is an artifact resulting from the omission of spin-torsion topology in the standard metric.


## Use

1. **Compilation :**
   ```bash
   make clean
   make -j
   cd python
   python setup.py build
   python setup.py install
   ```

2. **Exécution :**
   Use the provided `.ini` files to test how the background changes. The torsion parameters are set to $\Gamma = 13/12$ by default in `source/background.c`.

## References

- **Full paper:** [Gomes, F. (2026). The Kerr-Cartan Cosmological Model. Zenodo.](https://doi.org/10.5281/zenodo.19570177)
- **Original code:** [CLASS Repository](https://github.com/lesgourg/class_public)

***
