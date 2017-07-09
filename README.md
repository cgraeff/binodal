# binodal

This code models a two phase system and determines the transition between a
hadron and a quark phase. The transition surface is known as the *binodal*.

Both phases are described by SU(2) NJL type models: the quark phase is described
by the NJLv model (M. Buballa, Physics Reports 407 (2005) 205-376) and the
hadron phase is described by the eNJL model (Pais, Menezes, and Providência
Physical Review C 93, 065805 (2016)).

## Obtaining this program
The latest version of this code is available at
[github.com/cgraeff/binodal](https://github.com/cgraeff/binodal).
A copy can be easily obtained if your system have `git` installed, just issue
the following command in a terminal:
```
git clone https://github.com/cgraeff/binodal.git
```
A submodule `libdatafun` contains some functions that are also used in other
projects. To be able to properly compile, after cloning, issue the following
commands
```
git submodule init
```
and
```
git submodule update
```

Also, you may download both the program and the submodule and place the contents
of `libdadafun` inside `binodal/src/libdatafun/`.

## Requisites

To build and run this code `make`, a C compiler (default is `gcc`) and the
GSL (GNU Scientific Library) must be installed on the system. For plotting
graphics, `gnuplot` is also required.

On Linux the installation varies from distribution to distribution, but generally
there are two packages, one for regular use and one for developing.
The development version is the one needed for compilation of this code.

On OSX, the location of GSL headers and libraries is assumed to be
`/usr/local/include` and `/usr/local/lib`. This is the default if GSL is
installed from Homebrew. To change that, edit the variables at `src/Makefile`.

## Build and run instructions

Basic build (generate `binodal` executable):
* Build with `make` in the root dir;
* Debug builds may be generated with `make debug`;

Running:
* When executed, `binodal` will calculate the binodal with the
  default parameterizations;
* The following options are available:
 * `-q par`: uses `par` quark parameters set;
 * `-h par`: uses `par` hadron parameters set;
 * `-t val`: uses `val` for temperature value. Must be a floating point
             (`double`) value;
 * `-l`: list available parameters set;
 * `-s`: supress progress information during execution;
 * `-d`: write results using a dir structure;
 * `-a`: run tests. Automatically sets `-d`;
 * `-u`: prints usage;

For easier running and testing, the following commands are provided (all to be
issued at the root dir):
* Run default parameterizations with `make run` (implies `-d`);
* Remove product files with `make clean`;
* Arguments may be passed with `make ARGS="-p Set" run`, where `-p Set`
  stand for example arguments;
* Run for many combinations of parameters sets at once using `make multirun` or
  `make ARGS="-t 10" multirun`, where `-t 10` stands for example arguments. The
  sets must be listed in the `QUARK_MULTIRUN_SETS` and `HADRON_MULTIRUN_SETS`
  variableS in the `Makefile` at the root dir;
* Run tests with `make tests`;

When running on the default tree (that is, on the cloned or downloaded dir), the
results can be plotted with included gnuplot scripts:
* Plot results with `make graph`;
* Plot results for multirun with `make mgraph`.
* Plot tests with `make tgraph`;

## Known limitations
* Due to numerical instabilities, the pressures are not equal at some points.
  This must be further investigated.

## Code structure

### Parameters

Model parameter sets must be declared in `ParametersSetup()` in `Parameters.c`.
After declaration, the set must be appendend to the list of sets using
`AppendQuarkParametersSetToList()` or `AppendHadronParametersSetToList()`.
Numerical are declared at `NumericalParameters()` in the same file.

To choose the parameter sets, the function `SetParametersSet()` must be used.
This is used before the calculation of the main results and if no parameter sets
are explicitly requested in the command line, the first declared set of each
type is used. Each test in `RunTests()` (in `Tests.c`) should declare a set
pair.

Once set using `SetParametersSet()`, use `parameters.*` to access parameter
values. See the sctructure in `Parameters.h` to see available parameters groups.
Each group is declared in the header file of the `.c` file where the parameter
is needed. The sub-parameters are accessed using dot notation, e.g.:
```
double a_value = parameters.quark.model.bare_mass;
```

(To see a list of all parameters, just look at the model sets declarations and
numerical parameters declarations in `Parameters.c`.)

### Commandline options

The commandline options set in `CommandlineOptions.*` are globally available using
`options.an_option`.

### Tests

Tests should be declared in the function `RunTests()` in `Tests.c`. This function is
executed when the executable is called with `-a` option, or using `make tests`.

### Files and paths

Paths can be easily set with `SetFilePath()`. After that, all files created with
`OpenFile()` will be created at that path. After a section of code that uses
`SetFilePath()`, the path should be reset to the working dir with
`SetFilePath(NULL)`. This is important as other sections of code may expect to
write at the working dir while using `OpenFile()` and this is accomplished with
a `NULL` path.
