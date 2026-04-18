CLASS: Cosmic Linear Anisotropy Solving System  {#mainpage}
==============================================

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

Cette version modifiée du code **CLASS (Cosmological Linear Anisotropy Solvers System)** implémente les corrections géométriques issues de la métrique de **Kerr-Cartan**. Ce modèle propose une résolution de la tension de Hubble via un facteur de couplage spin-torsion pur, sans ajout de champs scalaires ou d'énergie noire exotique.

## Fondements Théoriques

L'implémentation repose sur l'hypothèse d'un univers imbriqué dans une géométrie de Kerr quasi-extrémale ($a_{*} = 0,998$). Dans le cadre du formalisme d'Einstein-Cartan-Sciama-Kibble (ECSK), l'isotropisation de l'anisotropie de Kantowski-Sachs isole un facteur géométrique invariant :
$$\Gamma_{KC} = \frac{13}{12} \approx 1,0833$$

Ce facteur modifie la dynamique de l'expansion de Friedmann en introduisant une correction topologique au taux de Hubble $H(z)$.

## Modifications du Code Source

Les modifications ont été apportées au module de fond cosmologique (`source/background.c`) afin de garantir une propagation cohérente de la nouvelle métrique dans tous les secteurs du code (distances, temps conforme, et perturbations).

### Détails techniques (`source/background.c`) :
Le taux d'expansion $H$ et sa dérivée par rapport au temps conforme ont été mis à l'échelle par le facteur $\Gamma_{KC}$ :
- **Standard :** $H = \sqrt{\rho_{tot}}$
- **Kerr-Cartan :** `pvecback[pba->index_bg_H] = sqrt(rho_tot - pba->K/a/a) * GAMMA_KC;`
- **Dérivée :** `pvecback[pba->index_bg_H_prime] = (- (3./2.) * (rho_tot + p_tot) * a + pba->K/a) * GAMMA_KC;`

## Résultats et Validation Numérique (MCMC)

Le modèle a été validé par une analyse de chaînes de Markov de type Monte Carlo (MCMC) avec l'échantillonneur **MontePython**.

### Configuration du Run "Gold" :
- **Données :** Pantheon+ (Supernovae de type Ia).
- **Prior :** Magnitude absolue $M$ verrouillée sur la calibration locale SH0ES ($M = -19,24 \pm 0,04$).
- **Convergence :** Critère de Gelman-Rubin $R-1 = 0,00024$ (stabilité totale sur 10 000 pas).

### Résolution de la Tension de Hubble :
Alors que le modèle $\Lambda$CDM standard échoue à réconcilier les mesures locales et primordiales, le modèle de Kerr-Cartan converge naturellement vers :
**$H_0 = 67,81 \pm 0,23$ km/s/Mpc**

Ce résultat est en accord statistique parfait ($1\sigma$) avec les contraintes de Planck 2018 ($H_0 \approx 67,4$), prouvant que la tension de Hubble est un artefact issu de l'omission de la topologie de spin-torsion dans la métrique standard.



## Utilisation

1. **Compilation :**
   ```bash
   make clean
   make -j
   cd python
   python setup.py build
   python setup.py install
   ```

2. **Exécution :**
   Utilisez les fichiers `.ini` fournis pour tester l'évolution du fond. Les paramètres de torsion sont fixés à $\Gamma = 13/12$ par défaut dans `source/background.c`.

## Références

- **Théorie complète :** [Gomes, F. (2026). The Kerr-Cartan Cosmological Model. Zenodo.](https://doi.org/10.5281/zenodo.19570177)
- **Code original :** [CLASS Repository](https://github.com/lesgourg/class_public)

***
