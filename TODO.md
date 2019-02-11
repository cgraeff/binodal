# TODO

- Temporaly disable negative or positive semiplane calculation 
- Check `DetermineHadronPressureAndDensities()` results in the first point of
  the `GridLoop()`
  - The solutions for mass and densities are not unreasonable. Apparently,
    changing the temperature did't even change the results.
  - The calculation fails when determining the entropy, due to values of
    $bar{n}_p$ going to zero. This eventually results in a 'nan'in the
    determination of the integrand
    - Try to write proper cases to determine the value of the entropy, look 
      the results obtained to that effect in the notes, there is a version
      that works well, as far I can remember
    - We may also temporaly disable the calculation of the entropy to be able to
      test other code paths
  - The determination of roots appear to be worse than on the zero temperature
    case
    - We, however, are now calculating at T = 100 MeV, whick may be too high
- In `HadronMassAndDensitiesSolutionEquation()`:
  - Try to use only one code path for both T != 0 and T == 0
    - I doubt that this is worth it
- Solve TODOs in `CalculationOverview.md`; The descriptions must be only to the
  solutions of complex equations;
- Solve TODOs in the code;
- Document units of inputs and outputs of functions
- Document return status in functions that may fail to calculate results;
- Why there is are functions `MySign()` and `MyAdapterFunction()` in
  `QuarkPhaseEOS.c`?
  - `MySign()` may be substituted by stlib `copysgn()`
  - `MyAdapterFunction()` is used by
    `QuarkMassAndRenormChemPotSolutionBissection()`, but why the generic name,
    and why `QuarkMassAndRenormChemPotSolutionBissection()` is so weird?
    - This is a quick hack to be able so determine the solution that minimizes
      the quark thermodynamic potential when there are multiple solution. The
      `MyAdapterFunction()` exists just as glue code to make reusing the old
      code possible.
- Reevaluate what are the relevant information that must be in `Flow.md`. Less
  information will make keeping it up to date easier. A theorethical explanation
  is out of scope.
- In `SolveBinodalForBarionicAndIsovectorChemicalPotentialsGrid()`, when
  determining the vacuum value of the hadron energy density (this value is the
  same as the vacuum hadron thermodynamic potential), we just assume that the
  value of the hadron mass that minimizes the thermodynamic potential is the
  nucleon mass. This may not be true for a given parameterization and would be
  better to just calculate it, as is done for the quark case.
- Try to solve both finite and zero temperature in one case in the function
  `HadronMassAndDensitiesSolutionEquation()`
- Change variables/constants:
  - Function call will be simpler if there are less input parameters
  - A constant ought to be recalled from `parameters.*` directly to minimize
    the number of input parameters of functions
  - A parameter must be a variable whenever a function is expected to be called
    multiple times with different values of that parameter
  - In simple functions like the `FermiDirac` distributions, the temperature
    would be a constant parameter, but it is cleaner to just input it as a
    variable. This avoids having to include `Parameters.h` for those functions.
  - In external functions, there is no way to access parameters directly, so it
    is necessary to use them as variables or insert them in apropriate structs.
    
- In `HadronBarionicDensity()`:
  - Does it have a special case for zero mass?
  - Does it need a cutoff?
    - **NOTE:** This is important: are the integration of Fermi Dirac functions
                limited to the cutoff as an upper limit?
                
- In `FermiDiracDistributionIntegralFromScalarDensity()`:
  - Is the upper limit the cutoff?
    - **I'm assuming that it is**
