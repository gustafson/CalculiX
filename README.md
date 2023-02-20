# Capabilities

CalculiX Extras builds on the [CalculiX](https://www.calculix.de/) code base

-   It extends CalculiX to write results into
    [ExodusII](https://gsjaardema.github.io/seacas/) format.
    Post processing can be accomplished with several readers. A
    recommended post processor is [Paraview](https://www.paraview.org/).
	![ccx_exodusII](/assets/images/ccx_exodusII.png)
-   It adds CUDA based solvers (obsolete).
    -   Please note the solver extras were research level code which
		used early libraries interfacing CUDA GPU capabilities. Those
		libraries have been superseeded by more modern solver
		technology. Without support, there has been no justification
		for maintaining these solver extensions.  Hence, the solver
		extras may work, but must be used with due diligence.
    -   Currently, the solvers can be called for static analysis in
        mechanical models in CalculiX.
        -   [Cuda-Cusp](https://github.com/cusplibrary/cusplibrary)
            -   (Outstanding performance for appropriate models.)
        -   [SuiteSparse Cholmod](https://people.engr.tamu.edu/davis/suitesparse.html)
            -   (Modest performance at last testing in ccx version
                2.11. However, since then CUDA has been added to
                Cholmod and it has not been subsequently testing with
                CalculiX Extras)
-   See
	[PastiX4CalculiX](https://github.com/Dhondtguido/PaStiX4CalculiX)
	for the current most promising enhancements for solvers callable
	by CalculiX.

# Build Instructions

1.  Clone this repository
    -   The [main branch](https://github.com/gustafson/CalculiX) is
		the as-distributed CalculiX code.
	-   Most people will want to checkout the
		[exodus](https://github.com/gustafson/CalculiX/tree/exodus) branch
		or one of the distribution specific branches
		([opensuse](https://github.com/gustafson/CalculiX/tree/opensuse),
		[ubuntu](https://github.com/gustafson/CalculiX/tree/ubuntu))
3.  These instructions are limited. They assume you are able to build
	CalculiX from source. 
	-   The code contains a detailed
        [Makefile](https://github.com/gustafson/CalculiX/blob/exodus/ccx/src/Makefile)
        which can be modified for the local environment.
	-   The distribution specific branches contain automated scripts
        for building in that distribution's environment.
4.  Send bug reports to peter.gustafson@wmich.edu

# Usage

## Exodus

-   ExodusII output format is called on job execution using the command
    line.
```console
ccx -i jobname -o exo
```
-   Open the .exo file with [Paraview](https://www.paraview.org/).

## Solvers (Obsolete)

-   The solvers can be called using the *static keyword:
```console
*static, solver=cudacusp
*static, solver=cholmod
```

## Solver benchmarking (Obsolete)

A rigorous round of
[benchmarking](https://arc.aiaa.org/doi/pdf/10.2514/6.2014-0346) has
been completed. If you are unable to download the paper, please
contact me.


# Known Bugs and Functional Deficits

## Exodus

-   There are known limitations for output requests with ExodusII
    -   Several output requests to ExodusII format are untested.
    -   It is not currently possible to change output requests
        *between steps*.  This is due to a limitation in the exodus
        format.
    -   Node sets are partially implemented but element sets are not
        implemented

## Solvers (Obsolete)

-   The solvers have primarily been applied to static solid mechanics models.
-   The solvers have have not been actively tested since version 2.11.


# Developing for CalculiX


-   Doxygen was used to generate documentation which may help with
	development for CalculiX. The [documentation is
	here](https://doxygen.openaircraft.com/ccx-doxygen/index.html).
	Note also the user documentation for
	[ccx](https://www.openaircraft.com/ccx-doc/ccx/index.html)
	and
	[cgx](https://www.openaircraft.com/ccx-doc/cgx/index.html)

# Acknowledgments

-   Guido Dhondt and Klaus Wittig for their excellent open source
    finite element contribution.
-   The National Science Foundation (for hardware grants which
    supported this work)
-   Western Michigan University and my students
