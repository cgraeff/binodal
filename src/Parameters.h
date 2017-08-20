//
//  Parameters.h
//  binodal
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#include "libdatafun/libdatafun.h"

#include "HadronPhaseEOS.h"
#include "QuarkPhaseEOS.h"

typedef struct _parameters
{
    // Parameters for the variables for which we are
    // performing the calculation of the EOS
    struct _VariableParameters{
        double temperature;                         // (MeV)
        double min_barionic_chemical_potential;     // (MeV)
        double max_barionic_chemical_potential;     // (MeV)
        double min_isovector_chemical_potential;    // (MeV)
        double max_isovector_chemical_potential;    // (MeV)
        int num_points;
    } variables;

    struct _hadron {
        HadronModelParameters model;
        HadronMassAndDensitiesSolutionParams mass_and_densities_solution;
    } hadron;

    struct _quark {
        QuarkModelParameters model;

        QuarkVacuumMassDeterminationParameters vacuum_mass_determination;
        QuarkRenormChemPotSolutionParameters renorm_chem_pot_solution;

        QuarkMassAndRenormChemPotSolParams mass_and_renorm_chem_pot_solution;

        QuarkThermodynamicPotMinDetermParams therm_pot_minimum;

        QuarkEntropyIntegrationParameters entropy_integration_params;
    } quark;

    IntegratorParameters fermi_dirac_integrals;
    IntegratorParameters therm_pot_free_gas_integral;

    UnidimensionalRootFindingParameters binodal_rootfinding_params;

} Parameters;

extern Parameters parameters;

void ParseModelParametersSets(void);
void SetParametersSet(char *quark_set_identifier,
                      char *hadron_set_identifier);

#endif /* Parameters_h */
