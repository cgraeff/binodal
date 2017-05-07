//
//  Parameters.c
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Constants.h"
#include "Parameters.h"
#include "CommandlineOptions.h"

// Chosen parameter set (globally accessible)
Parameters parameters;

static HadronModelParameters hadron_model_par_sets[256] = {0};
static int hadron_model_par_sets_count = 0;

static QuarkModelParameters quark_model_par_sets[256] = {0};
static int quark_model_par_sets_count = 0;

Parameters NewCopyOfParametersSetFromTemplate();
void AppendQuarkParametersSetToList(QuarkModelParameters a_set);
void AppendHadronParametersSetToList(HadronModelParameters a_set);
Parameters NumericalParameters();

void ParametersSetsSetup(void)
{
    // See header file for units and other relevant
    // information about the parameters below.

    // START DECLARATION OF QUARK PARAMETERS SETS:

    QuarkModelParameters pq;

    pq.parameters_set_identifier = "Buballa_1";
    pq.parameters_set_origin = "Set 1 from M. Buballa, "
                              "Nucl. Phys. A 611 (1996) 393-408";

    pq.bare_mass = 0.0;    // MeV
    pq.cutoff = 650.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.14 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.0;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "Buballa_2";

    pq.parameters_set_origin = "Set 2 from M. Buballa, "
                              "Nucl. Phys. A 611 (1996) 393-408";

    pq.bare_mass = 0.0;    // MeV
    pq.cutoff = 600.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.45 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.0;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "Buballa_3";

    pq.parameters_set_origin = "Set 3 from M. Buballa, "
                              "Nucl. Phys. A 611 (1996) 393-408";

    pq.bare_mass = 0.0;    // MeV
    pq.cutoff = 570.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.84 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.0;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "BuballaR_2";

    pq.parameters_set_origin = "Set 2 from M. Buballa, Physics "
                              "Reports 407 (2005) 205-376 with G_V = 0";

    pq.bare_mass = 5.6;    // MeV
    pq.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.44 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.0;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "BuballaR_2_GV";

    pq.parameters_set_origin = "Set 2 from M. Buballa, Physics "
                              "Reports 407 (2005) 205-376 with G_V = G_S";

    pq.bare_mass = 5.6;    // MeV
    pq.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.44 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "BuballaR_2_GV_0.25";

    pq.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports"
                              "407 (2005) 205-376 with G_V = 0.25 * G_S";

    pq.bare_mass = 5.6;    // MeV
    pq.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.44 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.25 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "BuballaR_2_GV_0.35";

    pq.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                              "407 (2005) 205-376 with G_V = 0.35 * G_S";

    pq.bare_mass = 5.6;    // MeV
    pq.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.44 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.35 * pq.G_S;

    AppendQuarkParametersSetToList(pq);


    ////

    pq.parameters_set_identifier = "BuballaR_2_GV_0.45";

    pq.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                              "407 (2005) 205-376 with G_V = 0.45 * G_S";

    pq.bare_mass = 5.6;    // MeV
    pq.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.44 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.45 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "BuballaR_2_GV_0.55";

    pq.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                              "407 (2005) 205-376 with G_V = 0.55 * G_S";

    pq.bare_mass = 5.6;    // MeV
    pq.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.44 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.55 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "BuballaR_2_GV_0.65";

    pq.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                                    "407 (2005) 205-376 with G_V = 0.65 * G_S";

    pq.bare_mass = 5.6;    // MeV
    pq.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.44 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.65 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "BuballaR_2_GV_0.75";

    pq.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                                    "407 (2005) 205-376 with G_V = 0.75 * G_S";

    pq.bare_mass = 5.6;    // MeV
    pq.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.44 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.75 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "PCP-0.0";

    pq.parameters_set_origin = "Set from PRD 94 094001, 2016";

    pq.bare_mass = 5.1;    // MeV
    pq.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.11 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.0 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "PCP-0.1";

    pq.parameters_set_origin = "Set from PRD 94 094001, 2016";

    pq.bare_mass = 5.1;    // MeV
    pq.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.11 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.1 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "PCP-0.2";

    pq.parameters_set_origin = "Set from PRD 94 094001, 2016";

    pq.bare_mass = 5.1;    // MeV
    pq.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.11 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.2 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "PCP-0.3";

    pq.parameters_set_origin = "Set from PRD 94 094001, 2016";

    pq.bare_mass = 5.1;    // MeV
    pq.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.11 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.3 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "PCP-0.4";

    pq.parameters_set_origin = "Set from PRD 94 094001, 2016";

    pq.bare_mass = 5.1;    // MeV
    pq.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.11 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.4 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    ////

    pq.parameters_set_identifier = "PCP-0.5";

    pq.parameters_set_origin = "Set from PRD 94 094001, 2016";

    pq.bare_mass = 5.1;    // MeV
    pq.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    pq.G_S = 2.11 * pow(CONST_HBAR_C / pq.cutoff, 2.0);
    pq.G_V = 0.5 * pq.G_S;

    AppendQuarkParametersSetToList(pq);

    //  DECLARATION OF HADRON PARAMETERS SETS

    HadronModelParameters ph;

    ph.parameters_set_identifier = "eNJL1";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
                               "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 4.855;         // (fm)^2
    ph.G_V = 4.65;          // (fm)^2
    ph.G_SV = -6.583;       // (fm)^8
    ph.G_RHO = 0.5876;      // (fm)^2
    ph.G_VRHO = 0.0;        // (fm)^8
    ph.G_SRHO = 0.0;        // (fm)^8
    ph.cutoff = 388.189;    // (MeV)
    ph.bare_mass = 0.0;     // (MeV)
    ph.nucleon_mass = 939.0; // (MeV)

  	AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL1OmegaRho1";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 4.855;
    ph.G_V = 4.65;
    ph.G_SV = -6.583;
    ph.G_RHO = 0.5976;
    ph.G_VRHO = -1.0;
    ph.G_SRHO = 0.0;
    ph.cutoff = 388.189;
    ph.bare_mass = 0.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL1OmegaRho2";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 4.855;
    ph.G_V = 4.65;
    ph.G_SV = -6.583;
    ph.G_RHO = 0.6476;
    ph.G_VRHO = -6.0;
    ph.G_SRHO = 0.0;
    ph.cutoff = 388.189;
    ph.bare_mass = 0.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL2";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 3.8;
    ph.G_V = 3.8;
    ph.G_SV = -4.228;
    ph.G_RHO = 0.6313;
    ph.G_VRHO = 0.0;
    ph.G_SRHO = 0.0;
    ph.cutoff = 422.384;
    ph.bare_mass = 0.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL2OmegaRho1";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 3.8;
    ph.G_V = 3.8;
    ph.G_SV = -4.288;
    ph.G_RHO = 0.6413;
    ph.G_VRHO = -1.0;
    ph.G_SRHO = 0.0;
    ph.cutoff = 422.384;
    ph.bare_mass = 0.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL3";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 1.93;
    ph.G_V = 3.0;
    ph.G_SV = -1.8;
    ph.G_RHO = 0.65;
    ph.G_VRHO = 0.0;
    ph.G_SRHO = 0.0;
    ph.cutoff = 534.815;
    ph.bare_mass = 0.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL3SigmaRho1";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 1.93;
    ph.G_V = 3.0;
    ph.G_SV = -1.8;
    ph.G_RHO = 0.0269;
    ph.G_VRHO = 0.0;
    ph.G_SRHO = 0.5;
    ph.cutoff = 534.815;
    ph.bare_mass = 0.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL1m";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 1.3833;
    ph.G_V = 1.781;
    ph.G_SV = -2.943;
    ph.G_RHO = 0.7;
    ph.G_VRHO = 0.0;
    ph.G_SRHO = 0.0;
    ph.cutoff = 478.248;
    ph.bare_mass = 450.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL1mSigmaRho1";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 1.3833;
    ph.G_V = 1.781;
    ph.G_SV = -2.943;
    ph.G_RHO = 0.0739;
    ph.G_VRHO = 0.0;
    ph.G_SRHO = 1.0;
    ph.cutoff = 478.248;
    ph.bare_mass = 450.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL2m";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 1.078;
    ph.G_V = 1.955;
    ph.G_SV = -2.74;
    ph.G_RHO = 0.75;
    ph.G_VRHO = 0.0;
    ph.G_SRHO = 0.0;
    ph.cutoff = 502.466;
    ph.bare_mass = 500.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);

  	////

    ph.parameters_set_identifier = "eNJL2mSigmaRho1";
    ph.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, "
							  "Phys. Rev. C 93, 065805 – Published 8 June 2016";

    ph.G_S = 1.078;
    ph.G_V = 1.955;
    ph.G_SV = -2.74;
    ph.G_RHO = -0.1114;
    ph.G_VRHO = 0.0;
    ph.G_SRHO = 1.0;
    ph.cutoff = 502.466;
    ph.bare_mass = 500.0;
    ph.nucleon_mass = 939.0;

    AppendHadronParametersSetToList(ph);


    // END DECLARATION OF PARAMETERS SETS

    // Verify that the set identifiers are unique
    for (int i = 0; i < quark_model_par_sets_count; i++){
        for (int j = 0; j < quark_model_par_sets_count; j++){

            if (i == j)
                continue;

            if(!strcasecmp(quark_model_par_sets[i].parameters_set_identifier,
                           quark_model_par_sets[j].parameters_set_identifier)){
                printf("Two parameters sets share the \"%s\" identifier,"
                       " but it should be unique.\n",
                       quark_model_par_sets[i].parameters_set_identifier);
                exit(EXIT_FAILURE);
            }
        }
    }

    for (int i = 0; i < hadron_model_par_sets_count; i++){
        for (int j = 0; j < hadron_model_par_sets_count; j++){

            if (i == j)
                continue;

            if(!strcasecmp(hadron_model_par_sets[i].parameters_set_identifier,
                           hadron_model_par_sets[j].parameters_set_identifier)){
                printf("Two parameters sets share the \"%s\" identifier,"
                       " but it should be unique.\n",
                       hadron_model_par_sets[i].parameters_set_identifier);
                exit(EXIT_FAILURE);
            }
        }
    }

    // If asked to, print parameters sets and exit
    if (options.list_available_parameterizations){
        for (int i = 0; i < quark_model_par_sets_count; i++){
            printf("Quark parameters set %s\n"
                   "\tOrigin: %s\n",
                   quark_model_par_sets[i].parameters_set_identifier,
                   quark_model_par_sets[i].parameters_set_origin);
        }
        for (int i = 0; i < hadron_model_par_sets_count; i++){
            printf("Hadron parameters set %s\n"
                   "\tOrigin: %s\n",
                   hadron_model_par_sets[i].parameters_set_identifier,
                   hadron_model_par_sets[i].parameters_set_origin);
        }
        exit(EXIT_SUCCESS);
    }

    // If there isn't at least one set, exit
    if (quark_model_par_sets_count == 0 || hadron_model_par_sets_count == 0){
        printf("There are not parameters set declared. Declare at least one set.\n");
        exit(EXIT_FAILURE);
    }
}

void SetParametersSet(char * quark_set_identifier, char * hadron_set_identifier)
{
    parameters = NumericalParameters();

    // If the identifier is null, return default case
    if (quark_set_identifier == NULL){
        parameters.quark_model = quark_model_par_sets[0];
    }
    else{
        for (int i = 0; i < quark_model_par_sets_count; i++){

            QuarkModelParameters p = quark_model_par_sets[i];

            if (!strcasecmp(p.parameters_set_identifier, quark_set_identifier)){
                parameters.quark_model = p;
            }
            else{
                printf("Parameters set %s unrecognized.\n"
                       "Use -l to list available parameterizations.\n",
                        quark_set_identifier);
                exit(EXIT_FAILURE);
            }
        }
    }

    if (hadron_set_identifier == NULL){
        parameters.hadron_model = hadron_model_par_sets[0];
    }
    else{
        for (int i = 0; i < hadron_model_par_sets_count; i++){

            HadronModelParameters p = hadron_model_par_sets[i];

            if (!strcasecmp(p.parameters_set_identifier, hadron_set_identifier)){
                parameters.hadron_model = p;
            }
            else{
                printf("Parameters set %s unrecognized.\n"
                       "Use -l to list available parameterizations.\n",
                        hadron_set_identifier);
                exit(EXIT_FAILURE);
            }
        }
    }

    return;
}

Parameters NumericalParameters()
{
    Parameters p;

    p.variables.temperature = 0.0; // (MeV)
    p.variables.proton_fraction = 0.5;

    // Low lower_bound but not zero, as it may be problematic if bare_mass == 0
    // upper_bound near the value of the nucleon mass.
    p.vacuum_mass_determination.lower_bound = 1.0E-3;
    p.vacuum_mass_determination.upper_bound = 1000.0;
    p.vacuum_mass_determination.abs_error = 1.0E-5;
    p.vacuum_mass_determination.rel_error = 1.0E-5;
    p.vacuum_mass_determination.max_iterations = 2000;

    p.q_renorm_chem_pot_finding.lower_bound = 1.0E-3;
    p.q_renorm_chem_pot_finding.upper_bound = 1000.0;
    p.q_renorm_chem_pot_finding.abs_error = 1.0E-5;
    p.q_renorm_chem_pot_finding.rel_error = 1.0E-5;
    p.q_renorm_chem_pot_finding.max_iterations = 2000;

    // In the following, the guesses are used for the first iteration.
    // Any subsequent iteration uses the values of mass
    // and renormalized chemical potential from the previous
    // iteration as new guesses. This shall work well for
    // the "stepping" evolution of a control variable.
    p.simultaneous_solution.max_iter = 8000;
    p.simultaneous_solution.barionic_density_guess = 0.10;
    p.simultaneous_solution.up_quark_mass_guess = 400.0; // (MeV)
    p.simultaneous_solution.down_quark_mass_guess = 400.0; // (MeV)
    p.simultaneous_solution.hadron_mass_guess = 939.0; // (MeV)
    p.simultaneous_solution.abs_error = 1.0E-8;
    p.simultaneous_solution.rel_error = 1.0E-8;

    p.fermi_dirac_integrals.lower_limit = 0.0; // (MeV)
    p.fermi_dirac_integrals.upper_limit = 2000.0; // (MeV)
    p.fermi_dirac_integrals.max_interval_num = 8000;
    p.fermi_dirac_integrals.integration_key = GSL_INTEG_GAUSS61;
    p.fermi_dirac_integrals.max_sub_interval = 8000;
    p.fermi_dirac_integrals.abs_error = 1.0E-10;
    p.fermi_dirac_integrals.rel_error = 1.0E-10;

    p.therm_pot_free_gas_integral.lower_limit = 0.0; // (MeV)
    p.therm_pot_free_gas_integral.upper_limit = 2000.0; // (MeV)
    p.therm_pot_free_gas_integral.max_interval_num = 8000;
    p.therm_pot_free_gas_integral.integration_key = GSL_INTEG_GAUSS61;
    p.therm_pot_free_gas_integral.max_sub_interval = 8000;
    p.therm_pot_free_gas_integral.abs_error = 1.0E-10;
    p.therm_pot_free_gas_integral.rel_error = 1.0E-10;

    /* pq.simultaneous_solution.entropy.lower_limit = 0.0; */
    /* pq.simultaneous_solution.entropy.upper_limit = parameters.model.cutoff; */
    /* pq.simultaneous_solution.entropy.abs_error = 1.0E-3; */
    /* pq.simultaneous_solution.entropy.rel_error = 1.0E-3; */
    /* pq.simultaneous_solution.entropy.max_sub_interval = 1000; */
    /* pq.simultaneous_solution.entropy.integration_key = GSL_INTEG_GAUSS61; */

    return p;
}

void AppendQuarkParametersSetToList(QuarkModelParameters a_set)
{
    // Append to list of parameterizations
    quark_model_par_sets[quark_model_par_sets_count] = a_set;
    quark_model_par_sets_count++;
}

void AppendHadronParametersSetToList(HadronModelParameters a_set)
{
    // Append to list of parameterizations
    hadron_model_par_sets[hadron_model_par_sets_count] = a_set;
    hadron_model_par_sets_count++;
}

