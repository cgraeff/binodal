# Calculation Overview

General flow of the program:
- Execution starts (file `main.c`):
  - check commandline parameters;
  - determine if the main loop or tests should be run;
    
- If tests should be run:
  - execute tests (file `Tests.c`);
    
- If normal run (binodal determination)
  - execute `SolveBinodalForBarionicAndIsovectorChemicalPotentialsGrid()`
      (file `Loop.c`);
    - determine hadron vacuum energy density: using the hadron vacuum mass 
      (the nucleon mass), zero temperature, and zero hadron chemical potential
      (in `HadronVacuumEnergyDensity()`). This result is a constant that is used
      to set the hadron thermodynamic potential minimum to zero. (Note that at
      zero hadron chemical potential, the hadron energy density and the hadron
      thermodynamic potetial are equal.) 
    - determine quark vacuum masses: determine the values of the quark masses
      that causes the thermodynamic potential to be a minimum (in 
      `QuarkVacuumMassDetermination()`);
    - determine the quark vacuum thermodynamic potential: determine the value of
      the quark thermodynamic potential using zero temperature and chemical
      potential, with quark masses equal to the vacuum values (in
      `QuarkThermodynamicPotential()`). This result, which is constant, is used
      to set the thermodynamic potential to zero at the minimum value.
      
    - Divide the hadron barionic chemical potential vs hadron isovector
      chemical potential phase space in a grid, looking for values of hadron
      and quark phases pressures (in `GridLoop()`) for which there is a phase
      coexistence (we look for points of the grid which are *binodal points*).
      The negative range of values of hadron isovector chemical potential is
      swept after the positive range;
      - [In `GridLoop()`], the grid is swept for all values of hadron barionic
        chemical potential for a given value of isovector chemical potential,
        then the hadron isovector chemical potential value is increased (in
        absolute value);
      - At each point of the grid, starting from initial guesses of hadron
        mass, proton density, neutron density, up quark mass, and down quark
        mass, determine pressure for hadron and quark phases (in function
        `BinodalPointCandidate()`, in `Binodal.c`). The guesses for the
        parameters are necessary as the process of solution involves
        root-finding algorithms;
        * [In `BinodalPointCandidate()`], *see below*.
      - The initial guesses of the parameters for a given grid point are the
        ones determined as proper solutions for the former grid point. For the
        initial point of the calculation, the guesses are input data and must
        be declared in `Parameters.c`;
      - The Gibbs' conditions for phase coexistence are
        $T^Q = T^H$, $\mu_B^Q = \mu_B^H$, and $P^Q = P^H$.
        The first condition is just a matter of choosing the same temperature
        for both phases; The second condition is used to determine the chemical
        potentials of the quarks; Determining if the point is a coexistence
        point (a binodal point) is just a matter of checking that the hadron and
        quark phases are equal (in pratice, however, we must check that the
        difference in the values of pressures is less than a given tolerance);
            
Determination of hadron and quark phases pressures
(in `BinodalPointCandidate()`):
- From values of barionic and isovector hadron chemical potentials, determine
  proton and neutron chemical potentials (in functions
  `ProtonChemicalPotential()` and `NeutronChemicalPotential()`);
- Determine hadron mass, proton density, neutron density, and hadron pressure
  (in `DetermineHadronPressureAndDensities()`);
  - [In in `DetermineHadronPressureAndDensities()`], from guesses for proton and
    neutron densities and chemical potentials, as well as a hadron mass guess,
    determine adequate hadron mass, proton and neutron densities (in function
    `HadronMassAndDensitiesSolution()`);
    - [In `HadronMassAndDensitiesSolution()`], we solve the gap 
      equation
      $M = m - 2G_s\rho_s + 2 G_{sv}\rho_s\rho^2 + 2G_{s\rho}\rho_s\rho_3^2$,
      determining the values of mass $M$, proton density $\rho^p$, and
      neutron density $\rho^n$ densities, for which the above relation is
      true, given values of proton and neutron chemical potentials. This is
      done with a multidimensional root-finder. The solution is particularly
      complex due to the fact that the three solutions may abruptly go to
      zero.                                                                     **TODO 2: try to change var. mappings**
  - Determine kinectic energy density (in `HadronKinecticEnergyDensity()`)
    - [In `HadronKinecticEnergyDensity()`], we solve the equation                                                            **TODO: check sol. in the code (T != 0)**
      $\epsilon_{\text{kin}} = 2 n_f n_c \int \frac{d^3p}{(2\pi)^3}
                                              \frac{p^2}{\sqrt{p^2 + M^2}
                                              (\eta_+ - \eta_-)$
  - Determine proton and neutron scalar densities (in `HadronScalarDensity()`)
    - [In `HadronScalarDensity()`], we solve
      $\rho_s = 2 n_f n_c \int \frac{d^3p}{(2\pi)^3}
                               \frac{M}{\sqrt{p^2 + m^2}}
                               (\eta_+ - eta_-)$                                
  - Determine the hadron thermodynamic potential
    (in `HadronThermodynamicPotential()`)
    - [In `HadronThermodynamicPotential()`], we solve 
      $\omega = \langle \mathcal{H}_{\text{MF}}\rangle
                - (\mu_p \rho_p + \mu_n \rho_n)
                -Ts$                         
      - Determine entropy (in `HadronEntropy()`) through
        $s = -2 n_c n_f \int \frac{d^3p}{(2\pi)^3}
                             [\eta_+ \ln\eta_+ + (1 - \eta_+)\ln(1 - \eta_+)
                              \eta_- \ln\eta_- + (1 - \eta_-)\ln(1 - \eta_-)]$
  - Determine hadron pressure (it's just $P = - \omega$)
- Using Gibbs' conditions, determine up and down quarks chemical potentials;
- Given the values of up and down chemical potentials and mass guesses,         **TODO 1: check T != 0 case**
  determine up and down quark masses, as well as quark pressure (in function
  `DetermineQuarkPressureAndMasses()`);
  - [In `DetermineQuarkPressureAndMasses()`], given the values of up and down   **TODO 1: check T != 0 case**
    chemical potentials, determine up and down masses and renormalized chemical
    potentials (in function `QuarkMassAndRenormChemPotSolutionBissection()`);
    - [In `QuarkMassAndRenormChemPotSolutionBissection()`], using up and down   **TODO 1: check T != 0 case**
      chemical potentials, we scan the values of mass searching for solutions
      to the equations
      $M = m - 4 n_f n_c G_s \int \frac{d^3p}{(2\pi)^3} \frac{M}{E_p}
                                  (1 - n_p - \bar{n}) = 0$
      and
      $\mu = \tilde\mu + 4 n_f N_c G_v int \frac{d^3p}{(2\pi)^3} (n - \bar{n})$.
      Note that both functions must be solved self-consistently. Moreover, if
      there are multiple solutions, we must choose the one that minimizes the
      quark thermodynamic potential. Currently, we do that by exploiting the
      fact that the quark masses are always the same and scanning the values of
      $F(M) = M - m + 4 n_f n_c G_s \int \frac{d^3p}{(2\pi)^3} \frac{M}{E_p}
                                         (1 - n_p - \bar{n})$
      in a given range of values of $M$. Whenever there is a sign change, we
      calculate the value of $M$ for which $F(M) = 0$. If there are multiple
      solutions, we choose the one that minimizes the quark thermodynamic
      potential. For each evaluation of $F(M)$, appropriate values of
      $\tilde\mu$ are calculated so that the second equation above is satisfied.
      (**NOTE: This probably is quite costly**)
  - Given up and down masses, chemical potentials, and renormalized chemical    **TODO 1: check T != 0 case, should have entropy**
    potentials, determine the quark thermodynamic potential (in 
    functions `QuarkThermodynamicPotential()`)
    - [In `QuarkThermodynamicPotential()`], we calculate
      $\omega = \omega_M + \frac{(M - m)^2}{4G_s}
                           - \frac{(\mu-\tilde\mu)^2}{4G_v}
                           - \omega_$.
      Here                                                                      **TODO 1: check T != 0 case**
      $\omega_M = - 2 n_f n_c
                    \int \frac{d^3p}{(2\pi)^3}
                         \{E_p + T \ln(1 + \exp{-\frac{E_p - \tilde\mu}{T}}
                               + T \ln(1 + \exp{-\frac{E_p + \tilde\mu}{T}} \}$
      and $\omega_0$ is the constant value needed so that the quark
      thermodynamic potential is zero at its minimum value for zero temperature
      and quark chemical potentials (as explained previously)
  - Given the quark thermodynamic potential, determine quark-phase pressure
    (again, the pressure is just $P = - \omega$.
- Determine quarks renormalized chemical potentials (in function                **TODO 1: Review this function, as I don't remember how it works and why the name**
  `QuarkSelfConsistentRenormChemPot()`);
- Determine up and down quark densities (in `QuarkDensity()`)                   **TODO 1: check T != 0 case**
  - [In `QuarkDensity()`], solve
    $n = 2 n_f n_c \int \frac{d^3p}{(2\pi)^3} (n_p - \bar{n}_p)$;
    
## Additional remarks

- We are define the Fermi Dirac distributions as
  $n_p = \frac{1}{\exp\left(\frac{(E_p - \mu)}{T}\right) + 1}$
  $\bar{n}_p = \frac{1}{\exp\left(\frac{(E_p + \mu)}{T}\right) + 1}$
  with
  $E_p = \sqrt{p^2 + M^2}$.
  Many texts adopt the distributions as
  $n_\pm = \left[ \exp\left(\frac{\pm E_p - \mu}{T}\right) + 1\right]^{-1}$,
  whose correspondence with the definition adopted here is
  $n_+ \equiv n_p$
  $n_- \equiv 1 - \bar{n}_p$.
  Also,
  $1 - n_p = \left[1+\exp\left(-\frac{(E_p-\mu)}{T}\right)\right]^{-1}$
  $1 - \bar{n}_p = \left[1+\exp\left(-\frac{(E_p+\mu)}{T}\right)\right]^{-1}$
  
   

        
