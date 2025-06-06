[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# A stochastic radial basis function method for the global optimization of expensive functions

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper
[A stochastic radial basis function method for the global optimization of expensive functions](https://doi.org/10.1287/ijoc.1060.0182)
by R.G. Regis and C.A. Shoemaker.

This paper was awarded the _INFORMS Journal of Computing_ Test of Time Paper
Award for 2005-2009. We make the code available here for historical interest.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.1060.0182

https://doi.org/10.1287/ijoc.1060.0182.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Regis2007,
  author =        {R.G. Regis and C.A. Shoemaker},
  publisher =     {INFORMS Journal on Computing},
  title =         {{A stochastic radial basis function method for the global optimization of expensive functions}},
  year =          {2007},
  doi =           {10.1287/ijoc.1060.0182.cd},
  url =           {https://github.com/INFORMSJoC/1060.0182},
  note =          {Available for download at https://github.com/INFORMSJoC/1060.0182},
}
```


## Overview

The `MLMSRBF` folder contains MATLAB 9.4 codes for implementing the Multistart
Local Metric Stochastic RBF (MLMSRBF) method. This method can be used for
finding the global minimum of a computationally expensive objective function
defined over a hypercube, i.e., the only constraints are bound constraints on
all the decision variables and the difference between the upper bound and the
lower bound for each variable is constant. Ideally, the objective function
should be transformed so that each decision variable is between 0 and 1. The
MLMSRBF method is described in the paper:

R.G. Regis, C.A. Shoemaker. A stochastic radial basis function method for the
global optimization of expensive functions. INFORMS Journal on Computing, Vol.
19, No. 4, pp. 497-509, 2007.

Please cite this paper whenever you use these codes to generate results for your
research.

## Instructions

There are three main files: `LMSRBF.m`, `MLMSRBF.m` and `RunMLMSRBF.m`. The
first file implements the Local Metric Stochastic RBF (LMSRBF) algorithm
described in the above paper while the second file implements a multistart
version of LMSRBF called MLMSRBF. The third file is the interface file and it
contains information on the required input and the structure of the output when
running the MLMSRBF algorithm. I have also included two test functions (the 4-D
Shekel10 and 14-D Schoen14_100 functions) and two sample MATLAB m-files
(`SampleRun_Shekel10.m` and `SampleRun_Schoen14_100`) for running MLMSRBF on
these test functions.

To get started, type `help RunMLMSRBF` or type `SampleRun_Shekel10` on the
MATLAB prompt. This should work for MATLAB 9.4 or newer.

If you encounter any errors, please contact `rregis@sju.edu`.

Thank you very much for your interest.

Rommel G. Regis, Ph.D.
Professor
Department of Mathematics
Saint Joseph's University
Philadelphia, PA 19131, USA
