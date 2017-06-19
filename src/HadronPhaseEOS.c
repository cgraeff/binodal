//
//  ZeroTemperatureEOS.c
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-02-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>

#include "libdatafun/libdatafun.h"

#include "Parameters.h"
#include "Constants.h"
#include "Functions.h"

double HadronZeroedGapEquation(double mass,
                               void * params)
{
    hadron_gap_eq_input_params * p = (hadron_gap_eq_input_params *)params;

    double barionic_density = p->proton_density + p->neutron_density;
    double rho_3 = p->proton_density - p->neutron_density;

	double scalar_density = HadronScalarDensity(mass,
                                                p->proton_fermi_momentum,
                                                parameters.hadron_model.cutoff)
                            + HadronScalarDensity(mass,
                                                  p->neutron_fermi_momentum,
                                                  parameters.hadron_model.cutoff);

	double gap_1st_term = parameters.hadron_model.G_S * scalar_density;
	double gap_2nd_term = - parameters.hadron_model.G_SV
                            * scalar_density
                            * pow(barionic_density, 2.0);
    double gap_3rd_term = - parameters.hadron_model.G_SRHO
                            * scalar_density
                            * pow(rho_3, 2.0);

	return mass
           + 2.0 * CONST_HBAR_C
             * (gap_1st_term + gap_2nd_term + gap_3rd_term)
           - parameters.hadron_model.bare_mass;
}

double HadronSolveGapEquation(UnidimensionalRootFindingParameters rootfinding_params,
                              double proton_density,
                              double neutron_density,
                              double proton_fermi_momentum,
                              double neutron_fermi_momentum)
{
    hadron_gap_eq_input_params p;
    p.proton_density = proton_density;
    p.neutron_density = neutron_density;
    p.proton_fermi_momentum = proton_fermi_momentum;
    p.neutron_fermi_momentum= neutron_fermi_momentum;

    gsl_function func;
    func.function = &HadronZeroedGapEquation;
    func.params = &p;

    double return_result;
    int status = UnidimensionalRootFinder(&func,
                                          rootfinding_params,
                                          &return_result);

    if (status != 0){

        // If no root could be found, that means
        // we've reached chiral restoration
        return 0.0;
    }
    return return_result;
}

double HadronScalarDensity(double mass,
                           double fermi_momentum,
                           double cutoff)
{
    if (mass == 0.0){
        return 0.0;
    }

	return pow(CONST_HBAR_C, -3.0)
           * (mass / pow(M_PI, 2.0))
           * (F0(mass, fermi_momentum) - F0(mass, cutoff));
}

double HadronVacuumScalarDensity()
{
	return 2.0 * pow(CONST_HBAR_C, -3.0)
           * (parameters.hadron_model.nucleon_mass / pow(M_PI, 2.0))
           * (F0(parameters.hadron_model.nucleon_mass, 0.0)
              - F0(parameters.hadron_model.nucleon_mass,
                   parameters.hadron_model.cutoff));
}

double ProtonChemicalPotential(double proton_fermi_momentum,
							   double scalar_density,
							   double mass,
							   double barionic_density,
							   double proton_density,
							   double neutron_density)
{
	double E = sqrt(pow(mass, 2.0) + pow(proton_fermi_momentum, 2.0));

	double rho_3 = (proton_density - neutron_density);

	double rho_terms = (parameters.hadron_model.G_V * barionic_density
	                    + parameters.hadron_model.G_SV * barionic_density
                          * pow(scalar_density, 2.0)
						+ parameters.hadron_model.G_RHO * rho_3
						+ parameters.hadron_model.G_VRHO * pow(rho_3, 2.0)
                          * barionic_density
						+ parameters.hadron_model.G_VRHO * pow(barionic_density, 2.0)
                          * rho_3
                        + parameters.hadron_model.G_SRHO * pow(scalar_density, 2.0)
                          * rho_3)
                       * 2.0 * CONST_HBAR_C;

	return E + rho_terms;
}

double NeutronChemicalPotential(double neutron_fermi_momentum,
								double scalar_density,
								double mass,
								double barionic_density,
								double proton_density,
								double neutron_density)
{
	double E = sqrt(pow(mass, 2.0)
                    + pow(neutron_fermi_momentum, 2.0));

	double rho_3 = (proton_density - neutron_density);

	double rho_terms = (parameters.hadron_model.G_V * barionic_density
                        + parameters.hadron_model.G_SV * barionic_density
                          * pow(scalar_density, 2.0)
                        - parameters.hadron_model.G_RHO * rho_3
                        + parameters.hadron_model.G_VRHO * pow(rho_3, 2.0)
                          * barionic_density
                        - parameters.hadron_model.G_VRHO * pow(barionic_density, 2.0)
                          * rho_3
                        - parameters.hadron_model.G_SRHO * pow(scalar_density, 2.0)
                          * rho_3)
                       * 2.0 * CONST_HBAR_C;

	return E + rho_terms;
}


double HadronKinecticEnergyDensity(double mass,
                                   double proton_fermi_momentum,
                                   double neutron_fermi_momentum)
{
	double proton_kinectic_energy = (NUM_H_COLORS  / pow(M_PI, 2.0))
                                    * (F2(mass, proton_fermi_momentum)
                                       - F2(mass, parameters.hadron_model.cutoff));

	double neutron_kinectic_energy = (NUM_H_COLORS / pow(M_PI, 2.0))
                                     * (F2(mass, neutron_fermi_momentum)
                                        - F2(mass, parameters.hadron_model.cutoff));

	return (proton_kinectic_energy + neutron_kinectic_energy)
           / (pow(CONST_HBAR_C, 3.0));
}

double HadronVacuumKinecticEnergyDensity()
{
	return 2.0 * (NUM_H_COLORS / pow(M_PI, 2.0)) * (pow(CONST_HBAR_C, -3.0))
		   * (F2(parameters.hadron_model.nucleon_mass, 0)
              - F2(parameters.hadron_model.nucleon_mass,
                   parameters.hadron_model.cutoff));
}

double HadronVacuumEnergyDensity()
{
    double scalar_density_0 = HadronVacuumScalarDensity();

	return HadronVacuumKinecticEnergyDensity()
           + parameters.hadron_model.bare_mass * scalar_density_0
           - parameters.hadron_model.G_S * pow(scalar_density_0, 2.0) * CONST_HBAR_C;
}

double HadronEnergyDensity(double pressure,
	                	   double proton_chemical_potential,
					       double neutron_chemical_potential,
					       double proton_density,
					       double neutron_density)
{
	return - pressure
           + proton_chemical_potential * proton_density
           + neutron_chemical_potential * neutron_density;
}

double HadronPressure(double termodynamic_potential)
{
	return - termodynamic_potential;
}

double HadronThermodynamicPotential(double total_scalar_density,
                                    double barionic_density,
                                    double proton_density,
                                    double neutron_density,
                                    double proton_chemical_potential,
                                    double neutron_chemical_potential,
                                    double kinectic_energy_density)
{
	double rho_3 = proton_density - neutron_density;

	double omega = kinectic_energy_density
                   + parameters.hadron_model.bare_mass * total_scalar_density;

	omega += - proton_chemical_potential * proton_density
			 - neutron_chemical_potential * neutron_density;

	omega += (- parameters.hadron_model.G_S * pow(total_scalar_density, 2.0)
			  + parameters.hadron_model.G_V * pow(barionic_density, 2.0)
			  + parameters.hadron_model.G_SV
                * pow(total_scalar_density * barionic_density, 2.0)
			  + parameters.hadron_model.G_RHO * pow(rho_3, 2.0)
			  + parameters.hadron_model.G_VRHO
                * pow(barionic_density * rho_3, 2.0)
              + parameters.hadron_model.G_SRHO
                * pow(total_scalar_density * rho_3, 2.0))
             * CONST_HBAR_C;

	return omega;
}

double HPFermiMomentum(double density)
{
    return CONST_HBAR_C * pow(3.0 * pow(M_PI, 2.0) * density, 1.0 / 3.0);
}

