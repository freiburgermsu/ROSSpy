# import libraries
from matplotlib import pyplot 
from scipy.constants import nano, kilo, liter, day
from itertools import chain
from chempy.properties.water_density_tanaka_2001 import water_density
from pubchempy import get_compounds 
import subprocess
import datetime
import pandas
from math import pi, exp
import time
import os
import re


# calculation constants
grams_over_liters_h2o = 997.07 # @25 degrees celcius
grams_over_moles_h2o = 18.015 
kinematic_flow_velocity = 9.33E-7    #square meters / second
simulated_time_over_computational_time = 9.29    


# simulation constants
timesteps_over_cell = 1

# useful conditions and parameters
possible_answers = ['y', 'n']


class ROSSPkg():

    def __init__(self):
        # print the software introductory box
        print('\n\n')
        message = ('''* ROSS 1.1.1 *
        Reverse Osmosis Scaling Simulation
        by Andrew P. Freiburger, Ethan S. Chan, and Heather L. Buckley
        Summer 2021, Green Safe Water Lab, University of Victoria''')
        
        lines = message.split('\n')
        space = " " * indent
        width = max(map(len, lines))
        upper = f'╔{"═" * (width + indent * 2)}╗\n'
        middles = ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
        lower = f'╚{"═" * (width + indent * 2)}╝' 
        box = chain(upper, middles, lower)
        box_print = ''.join(box)
        print(box_print, '\n')

        
        # establish the general organization structures
        self.parameters = {}
        self.variables = {}
        self.results = {}

        
    def define_general(self, os, phreeqc_path, database_selection, simulation_type, simulation_title, indent = 1)
        '''
        Establish general conditions
        '''
        self.parameters['water_mw'] = get_compounds('water', 'name')[0].molecular_weight
        self.parameters['water_density'] = water_density()
        
        # parameterize the input file
        self.parameters['os'] =  os
        self.parameters['phreeqc_path'] = phreeqc_path
        self.parameters['database_selection'] = database_selection 
        self.parameters['simulation_type'] = simulation_type
        self.parameters['title'] = simulation_title
        
        title_line = 'TITLE\t %s' %(simulation_title)
        if os == 'Windows':
            self.parameters['database_path'] = phreeqc_path + '\\database\\%s.dat' %(database_selection)
            database_line = 'DATABASE %s' %(database_path)
            self.results['general_conditions'] = [database_line, title_line]
        else:
            self.results['general_conditions'] = [title_line]
        

    def transport(self, module_selection, module_characteristics = {}, quantity_of_modules, self.parameters['cells_per_module'], domain, output_perspective):
        '''
        Define the TRANSPORT block
        '''
        # parameterize the module 
        if module_selection == 'BW30-400':
            module_diameter =  201                   #mm
            permeate_tube_diameter =  29             #mm
            module_length =  1.016                   #m
            permeate_flow = 40                       #cubic meters / day
            max_feed_flow = 15.9                     #cubic meters / hour
            membrane_thickness = 250 * (constants.milli / constants.nano)   #mm
            feed_thickness = 0.7112                  #mm
            permeate_thickness = 0.3                 #mm
            polysulfonic_layer_thickness = 0.05      #mm 
            support_layer_thickness = 0.15           #mm
            repeated_membrane_thickness = 2 * membrane_thickness + feed_thickness + permeate_thickness + 2 * polysulfonic_layer_thickness + 2 * support_layer_thickness      #mm
            print('\nMembrane thickness:', '%s (mm)' %(repeated_membrane_thickness))

        elif module_selection == 'Custom':
            module_diameter = module_characteristics['diameter']                             #mm
            permeate_tube_diameter = module_characteristics['permeate_diameter']             #mm
            module_length = module_characteristics['length']                                 #m
            permeate_flow = module_characteristics['permeate_flow']                          #cubic meters / day
            max_feed_flow = module_characteristics['feed_flow']                              #cubic meters / day
            membrane_thickness = module_characteristics['membrane_thickness']                #mm
            feed_thickness = module_characteristics['feed_spacer_thickness']
            permeate_thickness = input('- What is the permeate spacer thickness? (mm) __ ')
            polysulfonic_layer_thickness = module_characteristics['polysulfonic_thickness']
            support_layer_thickness = module_characteristics['support_thickness']
            repeated_membrane_thickness = 2 * membrane_thickness + feed_thickness + permeate_thickness + 2 * polysulfonic_layer_thickness + 2 * support_layer_thickness
            print('Thickness of the repeated membrane unit (mm): %s' %(repeated_membrane_thickness))

        self.parameters['quantity_of_modules'] = quantity_of_modules
        self.parameters['cells_per_module'] = cells_per_module
        cell_length = module_length / self.parameters['cells_per_module']     #meters

        # calculate module properties
        module_cross_sectional_area = module_diameter**2 * math.pi / 4        #squared millimeters
        permeate_tube_cross_sectional_area = permeate_tube_diameter**2 * math.pi / 4     #squared millimeters
#             filtering_layer_thicknes = (module_diameter - permeate_thickness) / 2         #millimeters
        filtration_cross_sectional_area = module_cross_sectional_area - permeate_tube_cross_sectional_area        #squared millimeters
        feed_cross_sectional_area = (feed_thickness / repeated_membrane_thickness) * filtration_cross_sectional_area     #squared millimeters
        feed_volume = feed_cross_sectional_area * module_length * kilo**2   #cubic meters
        feed_mass = feed_volume * self.parameters['water_density'] * (1 / constants.liter) / constants.kilo    #kilograms, which assumes pure water for mass
        feed_moles = feed_mass * constants.kilo / self.parameters['water_mw'] 

        # calculate fluid flow characteristics
        feed_velocity = max_feed_flow / (feed_cross_sectional_area * kilo**2) / day     #meters / second
        reynolds_number = feed_velocity * module_length / kinematic_flow_velocity
        print('Reynold\'s number: ', reynolds_number)

        # calculate module cell characteristics
        self.parameters['feed_mass_cell'] = feed_mass / self.parameters['cells_per_module']      #kg
        feed_moles_cell = feed_moles / self.parameters['cells_per_module']    #moles

        # calculate simulation parameters that will populate the PHREEQC simulation 
        maximal_timestep = cell_length / feed_velocity * timesteps_over_cell         #seconds, from the Courant condition
        self.['permeate_removal_per_cell'] = maximal_timestep * permeate_flow / day * liter * self.parameters['water_density'] / self.parameters['water_mw'] / self.parameters['cells_per_module']      #moles / cell

        # define the transport black
        transport_line = '\nTRANSPORT'
        cells_line = '-cells\t\t\t%s' %(self.parameters['cells_per_module'])
        
        simulation_shifts = (2*self.parameters['cells_per_module']*self.parameters['quantity_of_modules'])
        shifts_line = '-shifts\t\t\t%s' %(simulation_shifts)
        lengths_line = '-lengths\t\t%s' %(cell_length)
        timestep_line = '-time_step\t\t%s\t# the Courant condition is satisfied with the cell_length of %s m and the feed velocity of %s m/s' %(maximal_timestep, cell_length, feed_velocity)
        initial_time_line = '-initial_time\t\t0'    
        boundary_conditions_line = '-boundary_conditions\tconstant\tconstant \t # Dirichlet boundary condition'
        
        if domain == 'single':
            domain_line = '-stagnant\t\t0\t0\t0\t0 \t # single domain\n#^stagnant cells\t^exchange factor\t^CP volume\t^bulk volume'
        elif domain == 'dual':
            domain_line = '-stagnant\t\t1\t1\t0.1\t0.9 \t # dual domain\n#^stagnant cells\t^exchange factor\t^CP volume\t^bulk volume'

        self.parameters['output_perspective'] = output_perspective
        if self.parameters['output_perspective'] == 'Scaling':
            punch_cells_line = '-punch_cells\t\t1-%s' %(self.parameters['cells_per_module'])
            punch_frequency_line = '-punch_frequency\t%s' %(self.parameters['cells_per_module'])
        elif self.parameters['output_perspective'] == 'Brine':
            punch_cells_line = '-punch_cells\t\t%s' %(self.parameters['cells_per_module'])
            punch_frequency_line = '-punch_frequency\t1'       

        # create the transport block
        self.results['transport_block'] = []
        self.results['transport_block'].extend((transport_line, cells_line, shifts_line, lengths_line, timestep_line, initial_time_line, boundary_conditions_line, domain_line, punch_cells_line, punch_frequency_line))


    def reaction(self, permeate_approach, permeate_efficiency, head_loss):
        '''
        Define the REACTION block
        '''
        cfs = []
        cell_moles = []
        reaction_parameters = []
        iteration = 0
        cumulative_cf = 1 
        self.results['reaction_block'] = []
        for module in range(quantity_of_modules):
            module_previous_moles_removed = 0
            if permeate_approach == 'Linear permeate':
                if iteration == 0:
                    initial_moles_cell = feed_moles_cell

                initial_moles_removed = permeate_removal_per_cell * 2 / (1 + math.exp(head_loss))
                final_moles_removed = initial_moles_removed * math.exp(head_loss)
                try:
                    removed_moles_slope = ((final_moles_removed - initial_moles_removed) / (self.parameters['cells_per_module'])) / permeate_efficiency
                except:
                    removed_moles_slope = 0
                average_moles_removed = (final_moles_removed + initial_moles_removed) / 2
                print('(Removed moles / cell) slope: ', removed_moles_slope)

                for cell in range(self.parameters['cells_per_module']):
                    removed_moles_in_cell = (cell * removed_moles_slope + initial_moles_removed)
                    reaction_parameters.append(removed_moles_in_cell)
                    module_previous_moles_removed += removed_moles_in_cell 

                cf = initial_moles_cell / (initial_moles_cell - module_previous_moles_removed)
                cumulative_cf *= cf

                initial_moles_cell -= module_previous_moles_removed
                previous_moles_removed = 0
                iteration += 1

            if permeate_approach == 'Linear CF':
                module_iteration = 0
                initial_cf = 1
                print('''\nConcentration factor (CF):
        The CF is the quotient of the effluent concentration and the influent concentration. The CF for single RO modules approximates 1.2, while that of in-series systems can reach 3-4.''')
                while True:
                    try:
                        final_cf = float(input('''- What is the final CF of your simulation? (1 < CF < 5) __ '''))
                        break
                    except:
                        print('ERROR: The passed value is not an accepted float.')
                        final_cf = float(input('''- What is the final CF of your simulation? (1 < CF < 5) __ '''))

                cf_slope = (final_cf - initial_cf) / self.parameters['cells_per_module']
                for cell in range(1,self.parameters['cells_per_module']+1):
                    cell_cf = cell * cf_slope + initial_cf
                    cfs.append(cell_cf)    

                for cf in cfs:
                    moles_to_be_removed =  feed_moles_cell - (feed_moles_cell / cf)
                    if module_iteration == 0:
                        initial_moles_cell = feed_moles_cell
                        reaction_parameters.append(moles_to_be_removed)
                    if module_iteration > 0:
                        module_previous_moles_removed += reaction_parameters[-1] 
                        reaction_parameter = moles_to_be_removed - module_previous_moles_removed
                        reaction_parameters.append(reaction_parameter)
                        initial_moles_cell -= moles_to_be_removed

                    module_iteration += 1

                cf = cfs[-1]
                cumulative_cf *= cf
                initial_moles_cell = feed_moles_cell - moles_to_be_removed   # moles_to_be_removed = the final quantity of moles is calculated with naming from the linear permeate flux method to be consistent  

            if self.parameters['simulation_type'] == 'transport':
                reaction_line.append('\n')
                for cell in range(1, self.parameters['cells_per_module']+1):
                    cell_number = cell + self.parameters['cells_per_module'] * module
                    reaction_line = 'REACTION %s' %(cell_number)
                    
                    if cell < self.parameters['cells_per_module']:
                        reaction_line += '\n\tH2O -1; %s' %(reaction_parameters[cell_number-1]) 
                    elif cell == self.parameters['cells_per_module']:
                        reaction_line += '''\n\tH2O -1; %s
        INCREMENTAL_REACTIONS \ttrue''' %(reaction_parameters[cell_number-1], module + 1)     

                    self.results['reaction_block'].append(reaction_line)

            elif self.parameters['simulation_type'] == 'evaporation':
                parameter_quantity = 15                          
                recursive_assymtote_multiplier = 1.335449219     # ??? arbitrary assignment of kg of water in the simulation?
                total_moles_removed = sum(reaction_parameters)
                initial_evaporation_parameter = total_moles_removed / recursive_assymptote_multiplier
                evaporation_reaction_parameters = ['0', initial_evaporation_parameter]
                for parameter in range(1, parameter_quantity):
                    evaporation_reaction_parameter = evaporation_reaction_parameters[parameter] * 1/4
                    evaporation_reaction_parameters.append(evaporation_reaction_parameter)

                # define the reaction block
                reaction_line = 'REACTION 1'
                reaction_line += '\n\tH2O -1; '
                reaction_line += ' '.join(evaporation_reaction_parameters) 
                self.parameters['reaction_block'] = [reaction_line]

            # the calculated reaction parameters will be added and printed to a generated PHREEQC input file
            final_solution_mass = initial_moles_cell * grams_over_moles_h2o / constants.kilo  #kg water mass
            final_cf_cell = self.parameters['feed_mass_cell'] / final_solution_mass
            print('Effluent module %s CF:' %(module + 1), final_cf_cell)

            if self.parameters['os'] == 'Windows':
                self.results['reaction_block'].append('#%s' %(permeate_approach))
                if permeate_approach == 'Linear permeate':
                    self.results['reaction_block'].append('''
        #Permeate efficiency parameter: %s
        #Head loss parameter: %s''' %(permeate_efficiency, head_loss))
                self.results['reaction_block'].append('''    #Effluent module %s:
        #Estimated CF: %s
        #Estimated solution mass: %s\n\n''' %(module + 1, cumulative_cf, final_solution_mass))


    def make_solutions():
        """
        Specify the SOLUTION block of the simulation.
        """
        global solution_block
        global elements
        global water_selection

        if database_selection == 'phreeqc':
            database_elements = ['Alkalinity', 'Alkalinity', 'Al', 'B', 'Ba', 'Br', 'C', 'C(4)', 'C(-4)', 'Ca', 'Cl', 'E', 'Fe', 'Fe(+3)', 'H', 'H(1)', 'H(0)', 'K', 'Li', 'Mg', 'Mn', 'Mn(+3)', 'Na', 'N', 'N(+3)', 'N(0)', 'N(-3)', 'P', 'F', 'Zn', 'Cd', 'Pb', 'Cu', 'Cu(+1)' 'O', 'O(-2)', 'S', 'S(6)', 'S(-2)', 'Si', 'Sr']

            database_species = ['HCO3', 'CaCO3', 'Al+3', 'B(OH)3', 'Ba+2', 'Br-', 'CO3-2', 'CO3-2', 'CH4', 'Ca+2', 'Cl-', 'e-', 'Fe+2', 'Fe+3', 'H+', 'H+', 'H2', 'K+', 'Li+', 'Mg+2', 'Mn+2', 'Mn+3', 'Na+', 'NO3-', 'NO2-', 'N2', 'NH4+', 'PO4-3', 'F-', 'Zn+2', 'Cd+2', 'Pb+2', 'Cu(+2)', 'Cu(+1)', 'H2O', 'H2O', 'SO4-2', 'SO4-2', 'HS-', 'H4SiO4', 'Sr+2']

        elif database_selection == 'pitzer':
            database_elements = ['Alkalinity', 'Alkalinity', 'B', 'Ba', 'Br', 'C', 'C(4)', 'Ca', 'Cl', 'E', 'Fe', 'H', 'H(1)', 'K', 'Li', 'Mg', 'Mn', 'Na', 'O', 'O(-2)', 'S', 'S(6)', 'Si', 'Sr']
            database_species = ['HCO3', 'CaCO3', 'B(OH)3', 'Ba+2', 'Br-', 'CO3-2', 'CO3-2', 'Ca+2', 'Cl-', 'e-', 'Fe+2', 'H+', 'H+', 'K+', 'Li+', 'Mg+2', 'Mn+2', 'Na+', 'H2O', 'H2O', 'SO4-2', 'SO4-2', 'H4SiO4', 'Sr+2']

        long_element_species = ['B(OH)3', 'CaCO3', 'CO3-2', 'SO4-2', 'H4SiO4']
        elements = []

        # create the solution line of the input file
        solution_block = []
        announcement = '\nParameterize the feed solution:'
        print(announcement, '\n', '='*len(announcement))
        solution_description = input('''- What is a description of your feed solution? __ 
        Default is none\t''')

        if self.parameters['simulation_type'] == 'transport':
            initial_solution_line = '\nSOLUTION 0\t%s' % (solution_description)
        elif self.parameters['simulation_type'] == 'evaporation':
            initial_solution_line = '\nSOLUTION 1\t%s' % (solution_description)
        solution_block.append(initial_solution_line)

        #=============================================================================
        # determine which predefined solution should be simulated
        water_selections = ['Red Sea', 'Mediterranean Sea', 'German Basin', 'Marcellus Appalachian Basin', 'Bakken Formation', 'Michigan Basin', 'Palo Duro Basin', 'Western Pennsylvannia Basin', 'Custom']   
        for selection in water_selections:
            print('< %s >' %(selection))

        print('The elemental descriptions of the predefined water bodies are simplified to comform with the most restrictive PHREEQC database, which is the pitzer database. Elements beyond the pitzer database must be manually added to the SOLUTION blocks of the input file. __ ')
        water_selection = input('- Which water geochemistry would you like to simulate?  __ ')
        while water_selection not in water_selections:
            print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
            water_selection = input('- Which water geochemistry would you like to simulate?')

        if water_selection != 'Custom':       
            if water_selection == 'Red Sea':
                elements = ['Mn', 'Fe', 'B', 'Cl', 'Na', 'S(6)', 'Ca', 'K', 'Mg', 'Sr', 'Ba', 'Li']
                concentrations = ['0.000306', '0.006281', '1.344', '24756', '16417.2', '9500', '774', '301', '1646', '8.3', '0.011', '0.228']
                units = []
                i = 0
                while i <= len(elements):
                    units.append('ppm')
                    i += 1
                element_comments = ['#average of Al-Taani et al., 2014 and Longinelli and Craig, 1967.', '#Sabine et al., 2002 as reported for the Indian Ocean', '#Al-Taani et al., 2014 // 4.00 is the default (?)', '#Al-Taani et al., 2014 for the Northern Gulf of Aqaba', '#Al-Taani et al., 2014 for the Northern Gulf of Aqaba', '#Al-Taani et al., 2014 ', '#https://www.lenntech.com/composition-seawater.htm in the Red Sea, and Longinelli and Craig, 1967', '#https://www.lenntech.com/composition-seawater.htm in the Red Sea, and Longinelli and Craig, 1967 describes [Na]=15834', '#Longinelli and Craig, 1967 and Llyod, 1967.', '#Abdel-Aal et al., 2015', '#Abdel-Aal et al., 2015', '#Abdel-Aal et al., 2015', '#Bernat, Church, and Allegre, 1972 from the Mediterranean', '#Bernat, Church, and Allegre, 1972 from the Mediterranean', '#Stoffyn-Egli and Mackenzie, 1984 for the Mediterranean Sea']

                # parameterize the environmental conditions
                temperature_line = 'temp \t 24.5 \t #average of Al-Taani et al., 2014 and Longinelli and Craig, 1967.'
                ph_line = 'pH \t 8.22 charge'
                pe_line = 'pe \t 0.2679 \t   #Al-Taani et al., 2014 // 4.00 is the default (?)'
                solution_unit = 'ppm'
                unit_line = 'units \t %s' %(solution_unit)
                solution_block.extend([temperature_line, ph_line, pe_line, unit_line])

            #=========================================================================================               

            elif water_selection == 'Mediterranean Sea':
                elements = ['Mn', 'Fe', 'B', 'Cl', 'Na', 'S(6)', 'Ca', 'K', 'Mg', 'Sr', 'Ba', 'Li']
                concentrations = ['0.000734', '0.0005566', '5.98', '22590', '15120', '9500', '610', '410', '1655', '9.0', '0.011', '0.228']
                units = []
                i = 0
                while i <= len(elements):
                    units.append('ppm')
                    i += 1
                element_comments = ['#El Sayed, Aminot, and Kerouel, 1994', '#Sherrell and Boyle, 1988; El Sayed, Aminot, and Kerouel, 1994; Sarthou and Jeandel, 2001','#Lee et al., 2009 for the Atlantic Ocean','#https://www.lenntech.com/composition-seawater.htm; Sherrell and Boyle, 1988; El Sayed, Aminot, and Kerouel, 1994 ; reduced from 24090 ppm', '#https://www.lenntech.com/composition-seawater.htm; Sherrell and Boyle, 1988; El Sayed, Aminot, and Kerouel, 1994 ; increased from 13410 ppm','#Longinelli and Craig, 1967 and Llyod, 1967 for the Red Sea', '#Abdel-Aal et al., 2015; Culkin and Cox, 1966; Kumgalz, 1982', '#Abdel-Aal et al., 2015; Culkin and Cox, 1966', '#Abdel-Aal et al., 2015; Culkin and Cox, 1966', '#Bernat, Church, and Allegre, 1972; Culkin and Cox, 1966','#Bernat, Church, and Allegre, 1972 ', '#Stoffyn-Egli and Mackenzie, 1984']

                # parameterize the environmental conditions
                temperature_line = 'temp \t 18.5 \t #Sherrell and Boyle, 1988; El Sayed, Aminot, and Kerouel, 1994.'
                ph_line = 'pH \t 8.22 charge'
                solution_unit = 'ppm'
                unit_line = 'units \t %s' %(solution_unit)
                solution_block.extend([temperature_line, ph_line, unit_line])

            #=========================================================================================   

            elif water_selection == 'Palo Duro Basin':
                elements = ['Mn', 'Fe', 'B', 'Cl', 'Na', 'S(6)', 'Ca', 'K', 'Mg', 'Sr', 'Ba', 'Li']
                concentrations = ['3.7', '11.12', '2.77', '138109.3', '80130', '1220', '5960', '369','1320', '99.5', '962', '30.33']
                units = []
                i = 0
                while i <= len(elements):
                    units.append('ppm')
                    i += 1
                element_comments = ['#Fisher and Kreitler, 1987','#Dresel and Rose, 2010, decreased proportionally to match the TDS values','#Kloppman et al., 2001, decreased proportionally to match the TDS values','#Fisher and Kreitler, 1987, Cl + Br + I','#Fisher and Kreitler, 1987; increased from 78400', '#Fisher and Kreitler, 1987','#Fisher and Kreitler, 1987','#Fisher and Kreitler, 1987','#Fisher and Kreitler, 1987','#Fisher and Kreitler, 1987','#Chapman et al., 2012','#Dresel and Rose, 2010, decreased proportionally to match the TDS values']

                # parameterize the environmental conditions
                temperature_line = 'temp \t 18.5 \t #Fisher and Kreitler, 1987.'
                alkalinity_line = 'Alkalinity 131 \t #Fisher and Kreitler, 1987'
                ph_line = 'pH \t 5.6 \t#Fisher and Kreitler, 1987'
                solution_unit = 'ppm'
                unit_line = 'units \t %s' %(solution_unit)
                solution_block.extend([temperature_line, alkalinity_line, ph_line, unit_line])

            #=========================================================================================  

            elif water_selection == 'Marcellus Appalachian Basin':
                elements = ['Mn', 'Fe', 'B', 'Cl', 'Na', 'S(6)', 'Ca', 'K', 'Mg', 'Sr', 'Ba', 'Li']
                concentrations = ['3000', '26.6', '20', '81900', '39630', '45', '8786', '350','841', '2415', '962', '95']
                units = []
                i = 0
                while i <= len(elements):
                    units.append('ppm')
                    i += 1
                element_comments = ['#Haluszczak, Rose, and Kump, 2013 [estimated from another Marcellus publication]', '#Chapman et al., 2012', '#Haluszczak, Rose, and Kump, 2013 [reported average form another Marcellus publication]','#Chapman et al., 2012','#Chapman et al., 2012; increased from 32800', '#Haluszczak, Rose, and Kump, 2013 [estimated from another Marcellus publication]', '#Chapman et al., 2012', '#Haluszczak, Rose, and Kump, 2013 [estimated from another Marcellus publication]', '#Chapman et al., 2012','#Chapman et al., 2012','#Chapman et al., 2012', '#Haluszczak, Rose, and Kump, 2013 [reported average from another Marcellus publication]']

                # parameterize the environmental conditions
                temperature_line = 'temp \t 24 \t #Dresel and Rose, 2010'
                alkalinity_line = 'Alkalinity 71 \t #Haluszczak, Rose, and Kump, 2013 [reported average from another Marcellus publication]'
                ph_line = 'pH \t 8.22'
                solution_unit = 'ppm'
                unit_line = 'units \t %s' %(solution_unit)
                solution_block.extend([temperature_line, alkalinity_line, ph_line, unit_line])


            #=========================================================================================       

            elif water_selection == 'Michigan Basin':
                elements = ['Mn', 'B', 'Cl', 'Na', 'S(6)', 'Ca', 'K', 'Mg', 'Sr', 'Ba', 'Li']
                concentrations = ['9.780', '283', '206062', '73500', '150', '35600', '1620', '6300', '980', '17500', '24']
                units = []
                i = 0
                while i <= len(elements):
                    units.append('ppm')
                    i += 1
                element_comments = ['#Peterman et al., 2017', '#Wilson and Long, 1993, esimated from adjacent samples', '#Cl + Br + I, Wilson and Long, 1993; increased from 183062','#Wilson and Long, 1993','#Wilson and Long, 1993', '#Wilson and Long, 1993','#Wilson and Long, 1993', '#Wilson and Long, 1993', '#Wilson and Long, 1993','#Peterman et al., 2017','#Wilson and Long, 1993']

                # parameterize the environmental solutions
                temperature_line = 'temp \t 33.8 \t #Wilson and Long, 1993'
                ph_line = 'pH \t 5.6 charge\t#Wilson and Long, 1993'
                solution_unit = 'ppm'
                unit_line = 'units \t %s' %(solution_unit)
                solution_block.extend([temperature_line, ph_line, unit_line])

            #=========================================================================================    

            elif water_selection == 'Western Pennsylvannia Basin':
                elements = ['Fe', 'Mn', 'B', 'Cl', 'Na', 'S(6)', 'Ca', 'K', 'Mg', 'Sr', 'Ba', 'Li']
                concentrations = ['90', '15.8', '3.934', '124144', '58500', '4', '19000', '193', '2520', '1470', '815', '45.5']
                units = []
                i = 0
                while i <= len(elements):
                    units.append('ppm')
                    i += 1
                element_comments = ['#Dresel and Rose, 2010','#Dresel and Rose, 2010','#Kloppman et al., 2001, decreased proportionally to match the TDS values','#Cl + Br + I, Dresel and Rose, 2010; increased from 124144','#Dresel and Rose, 2010','#Dresel and Rose, 2010, estimated from the adjacent sample values', '#Dresel and Rose, 2010','#Dresel and Rose, 2010','#Dresel and Rose, 2010','#Dresel and Rose, 2010','#Dresel and Rose, 2010','#Dresel and Rose, 2010']

                # parameterize the environmental solutions
                temperature_line = 'temp \t 24 \t #Dresel and Rose, 2010'
                alkalinity_line = 'Alkalinity 90 as HCO3-'
                ph_line = 'pH \t 5.6 \t#Wilson and Long, 1993'
                solution_unit = 'ppm'
                unit_line = 'units \t %s' %(solution_unit)
                solution_block.extend([temperature_line, alkalinity_line, ph_line, unit_line])

            #=========================================================================================  

            elif water_selection == 'German Basin':
                elements = ['Mn', 'B', 'Cl', 'Na', 'S(6)', 'Ca', 'K', 'Mg', 'Sr', 'Ba', 'Li']
                concentrations = ['9.780', '5.54', '198765', '121000', '5870', '1270', '1290', '1340', '24.19', '17500', '0.23']
                units = []
                i = 0
                while i <= len(elements):
                    units.append('ppm')
                    i += 1
                element_comments = ['#Peterman et al., 2017','#Kloppman et al., 2001','#Cl + Br, Kloppman et al., 2001; increased from 183065','#Kloppman et al., 2001','#Kloppman et al., 2001', '#Kloppman et al., 2001','#Kloppman et al., 2001','#Kloppman et al., 2001','#Kloppman et al., 2001','#Peterman et al., 2017','#Kloppman et al., 2001']

                # parameterize the environmental solutions
                temperature_line = 'temp \t 33.8 \t #Wilson and Long, 1993'
                ph_line = 'pH \t 6.8 charge\t#Kloppman et al., 2001' 
                solution_unit = 'ppm'
                unit_line = 'units \t %s' %(solution_unit)
                solution_block.extend([temperature_line, ph_line, unit_line])          

            #=========================================================================================  

            elif water_selection == 'Bakken Formation':
                elements = ['Mn', 'B', 'Cl', 'Na', 'S(6)', 'Ca', 'K', 'Mg', 'Sr', 'Ba', 'Li']
                concentrations = ['9.780', '283', '185979', '84010', '274', '18000', '5570','1620', '1210', '17500', '52.7']
                units = []
                i = 0
                while i <= len(elements):
                    units.append('ppm')
                    i += 1
                element_comments = ['#Peterman et al., 2017','#Peterman et al., 2017','#Cl + Br + F, Peterman et al., 2017','#Peterman et al., 2017, elevated from 84000 ppm to match the reported charge balance -1.3','#Peterman et al., 2017', '#Peterman et al., 2017','#Peterman et al., 2017','#Peterman et al., 2017','#Peterman et al., 2017','#Peterman et al., 2017','#Peterman et al., 2017']

                # parameterize the environmental solutions
                temperature_line = 'temp \t 33.8 \t #Wilson and Long, 1993'
                ph_line = 'pH \t 6.22 charge\t#Dresel and Rose, 2010' 
                solution_unit = 'ppm'
                unit_line = 'units \t %s' %(solution_unit)
                solution_block.extend([temperature_line, ph_line, unit_line])

            #=========================================================================================  

            # parameterize the environmental conditions
            for element in elements:
                element_index = elements.index(element)
                database_index = database_elements.index(element)

                element_concentration = concentrations[element_index]
                element_unit = units[element_index]
                elemental_formula = database_species[database_index]
                element_comment = element_comments[element_index]

                if len(element_concentration) > 3:
                    if elemental_formula in long_element_species:
                        element_line = '%s\t%s %s\tas %s\t%s'  %(element, element_concentration, element_unit, elemental_formula, element_comment)  
                    else:
                        element_line = '%s\t%s %s\tas %s\t\t%s'  %(element, element_concentration, element_unit, elemental_formula,element_comment)  
                elif len(element_concentration) <= 3:
                    if elemental_formula in long_element_species:
                        element_line = '%s\t%s %s\t\tas %s\t%s'  %(element, element_concentration, element_unit, elemental_formula, element_comment)    
                    else:
                        element_line = '%s\t%s %s\t\tas %s\t\t%s'  %(element, element_concentration, element_unit, elemental_formula, element_comment)    

                solution_block.append(element_line)

            # parameterize the inital water mass
            water_mass = self.parameters['feed_mass_cell']
            if water_selection == 'Bakken formation':
                water_line = '-water \t%s\t#TDS=300 ppthousand [before fudging]' %(water_mass)
            elif water_selection == 'German Basin':
                water_line = '-water \t%s\t#TDS=314 ppthousand [before fudging]' %(water_mass)
            else:
                water_line = '-water \t%s' %(water_mass)

            solution_block.append(water_line)

            #=========================================================================================  

        elif water_selection == 'Custom':

            # create the temperature line of the input file
            solution_temperature = input('''- What is the solution temperature? 
            Default = 25 (K) __  ''') or str(25)
            temperature_line = 'temp \t %s' %(solution_temperature)

            # create the pH line of the input file
            solution_ph = input('''- What is the solution pH?
            Default = 7.0  __ ''') or str(7)
            error = 'yes'
            while error == 'yes':
                try:
                    float(solution_ph) == True
                    error = 'no'
                except:
                    print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                    solution_ph = input('- What is the solution pH?')
            ph_charge_balance = input('''Should the pH be charge balanced? 
            Default = < n >  ;  < y > or < n > __ ''') or 'n'
            if ph_charge_balance == 'y':
                solution_ph += '\tcharge'
            while ph_charge_balance not in possible_answers:
                print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                ph_charge_balance = input('''- Should the pH be charge balanced? < y > or < n > __ ''')

            ph_line = 'pH \t %s' %(solution_ph)

            # create the pe line of the input file
            solution_pe = input('''- What is the solution pe (-log(electron activity))? 
            Default = 4.0  __  ''') or str(4)
            while not solution_pe.isnumeric():
                print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                solution_pe = input('''- What is the solution pe (-log(electron activity))? 
                Default = 4.0 __  ''') or str(4)

            pe_charge_balance = input('''- Should the pe phase be balanced? 
            Default = < n >  ;  < y > or < n > __ ''') or 'n'
            pe_charge = ''
            if pe_charge_balance == 'y':
                pe_phase = input('- What is the balanced phase name?')
                pe_phase_index = input('''- What is the phase saturation index?
                Default = 0  __ ''') or str(0)
                pe_charge = ' %s %s' %(pe_phase, pe_phase_index)

            pe_line = 'pe \t ' + solution_pe + pe_charge

            # create the redox line of the input file      
            solution_redox = input('''- Is a special pe calculated via a redox reaction? 
            Default = < n >  ; < y > or < n > __ ''') or 'n'
            if solution_redox == 'y':
                redox_couple = input('''- What is the redox couple?
                e.g. O(-2)/O(0) for oxygen with (-2) and (0) oxidation states\t''')
                redox_line = 'redox \t\t ' + redox_couple

            else:
                redox_line = ''

            # create the units line of the input file        
            unit_numerators = ['mg', 'ug', 'mmol', 'umol'] 
            unit_denominators = ['L', 'kgs', 'kgw'] 
            unit_concentration = ['ppt', 'ppm', 'ppb']
            print('Numerators:')
            for unit in unit_numerators:
                print('< %s >' %(unit)) 
            print('Denominators:')
            for unit in unit_denominators:
                print('< %s >' %(unit)) 
            print('Concentrations:')
            for unit in unit_concentration:
                print('< %s >' %(unit))
            solution_units = input('''- What are the concentration units of your solution?
            Default = mmol/kgw\t''') or 'mmol/kgw'
            if re.search('(\/)', solution_units):
                numerator, denominator = solution_units.split('/')
                numerator = numerator.strip()
                denominator = denominator.strip()
                while numerator not in unit_numerators or denominator not in unit_denominators:
                    print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                    solution_units = input('- What are the concentration units of your solution?   __ ')
            unit_line = 'unit\t ' + solution_units

            # create the units line of the input file        
            solution_density = input('''- What is the solution density? 
            Default = 1 (kg/L)  __  ''') or str(1)
            while not solution_pe.isnumeric():
                print('''ERROR: The entry is beyond the accepted values. 
                Please provide an accepted value''')
                solution_density = input('''- What is the solution density? 
                Default = 1 (kg/L)  __  ''') or str(1)
            density_line = 'density\t ' + solution_density

            # add the created lines to the SOLUTION block
            solution_block.extend([temperature_line, ph_line, pe_line, redox_line, unit_line, density_line])

            #creating the elemental lines of the input file     
            elemental_addition = 'y'
            element_lines = []
            elements = []
            while elemental_addition == 'y':
                print('''Element options for the < %s > database: \n< element >\t\tcompound''' %(database_selection))
                for possible_element in database_elements:
                    if possible_element not in ('Alkalinity', 'C(4)', 'H(1)', 'O(-2)', 'S(6)'):
                        print('< %s >' %(possible_element), '\t\t\t', database_species[database_elements.index(possible_element)]) 
                    if possible_element in ('Alkalinity', 'C(4)', 'H(1)', 'O(-2)', 'S(6)'):
                        print('< %s >' %(possible_element), '\t\t', database_species[database_elements.index(possible_element)])                 
                element = input('''- What element will you simulate?
                Type < done > when you have finished parameterizing elements.   __ ''')
                while element not in database_elements and element != 'done':
                    print('''The element is not supported by the %s database. A Solution_Master must be defined for the %s element''' %(database_selection, element))  
                    element = input('- What element will you simulate?  __ ')

                if element == 'done':
                    break

                elements.append(element)

                element_concentration = input('- What is the elemental concentration?')
                while not element_concentration.isnumeric():
                    print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                    element_concentration = input('- What is the elemental concentration?   __ ')

                special_element_unit = input('''- Are the units different than is previously defined?
                Default = < n >  ; < y > or < n > __ ''') or 'n'
                if special_element_unit == 'y':
                    element_unit = input('''- What is the element unit?
                    numerator = < mg >, < ug >, < mmol >, < umol > 
                    Denominator = < /L >, < /kgs> (kilograms), < /kgw > (kilograms of water) OR < ppt > (part per thousand), < ppm >, < ppb >
                    Default = mmol/kgw  __  ''') or 'mmol/kgw'
                elif special_element_unit == 'n':
                    element_unit = ''

                elemental_formula_selection = input('''- Is the formula different than the provided compound?
                < y > or < n > ; Default = < n >  __  ''') or 'n'
                if elemental_formula_selection == 'y':
                    elemental_formula = input('- What is the molecular formula of the compound?  __ ')
                elif elemental_formula_selection == 'n':
                    elemental_formula = database_species[database_elements.index(element)]

                possible_elemental_comment = input('''- Does the elemental information have associated commentary or a reference? 
        < y > or < n > __ ''')
                element_comment = ''
                if possible_elemental_comment == 'y':
                    elemental_comment = input('- What is the comment or reference?   __ ')
                    element_comment = '#%s' %(elemental_comment)
                elif possible_elemental_comment == 'n':
                    pass
                else:
                    while possible_elemental_comment not in ('y' or 'n'):
                        print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                        possible_elemental_comment = input('''- Does the elemental information have associated commentary or a reference? 
        < y > or < n > __ ''')

                if len(element_concentration) > 3 or element == 'Alkalinity':
                    if elemental_formula in long_element_species:
                        element_line = '%s\t%s %s\tas %s\t%s'  %(element, element_concentration, element_unit, elemental_formula, element_comment)  
                    else:
                        element_line = '%s\t%s %s\tas %s\t\t%s'  %(element, element_concentration, element_unit, elemental_formula, element_comment)  
                elif len(element_concentration) <= 3:
                    if elemental_formula in long_element_species:
                        element_line = '%s\t%s %s\t\tas %s\t%s'  %(element, element_concentration, element_unit, elemental_formula, element_comment)    
                    else:
                        element_line = '%s\t%s %s\t\tas %s\t\t%s'  %(element, element_concentration, element_unit, elemental_formula, element_comment)   

                solution_block.append(element_line)


            #creating the elemental lines of the input file   
            solution_isotope = input('''- Are exceptional elemental isotopes used in the simulation? 
            Default = < n >  ; < y > or < n >  __ ''') or 'n'                      
            if solution_isotope == 'y':
                isotope_mass = input('''- What is the isotopic mass? (amu)\t''')
                isotope_element = input('''- What is the elemental symbol of the isotope?
                e.g. < C > or < Si >  __  ''')
                isotope_line = '-isotope\t' + isotope_mass + isotope_element
                solution_block.append(isotope_line)


            #defining the mass of water
            water_mass = input('''- What is the mass of water in the simulation?
            Default = 1 (kg)  __  ''') or str(1)  
            while not water_mass.isnumeric():
                print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                water_mass = input('''- What is the mass of water in the simulation?
                Default = 1 (kg)  __  ''') or str(1)  

            water_line = '-water \t ' + water_mass
            solution_block.append(water_line)


            #looping to the next element in the simulation
            elemental_addition = input('''- Will you parameterize another element? < y > or < n > __ ''')
            while elemental_addition not in possible_answers:
                print('ERROR: Provide one of the provided answers')
                elemental_addition = input('''- Will you parameterize another element? < y > or < n > __ ''')

        else:
            print('''ERROR: The water selection exceeds the possible options. Examine the passed entry.''')

        #=========================================================================================  

        #parameterize the initial module solution
        if self.parameters['simulation_type'] == 'transport':
            feed_solution_line = '\nSOLUTION 1-%s\tInitial solution in the RO module' %(self.parameters['cells_per_module'] * quantity_of_modules) 
            solution_block.extend([feed_solution_line,
                                  'temp \t 25',
                                  'units \t ppm'])

            for element in elements:
                element_concentration = 0
                element_line = '%s\t%s'  %(element, element_concentration)    
                solution_block.append(element_line)

            water_line = '-water \t %s' %(self.parameters['feed_mass_cell'])
            solution_block.append(water_line)


    def make_equilibrium_phases():
        """
        Specify the EQUILIBRIUM_PHASES block of the simulation.
        """
        global equilibrium_phases_block
        global minerals
        global mineral_formulas
        global selected_minerals

        if database_selection == 'pitzer':

            minerals = ['Akermanite', 'Anhydrite', 'Anthophyllite', 'Antigorite', 'Aragonite', 'Arcanite', 'Artinite', 'Barite', 'Bischofite', 'Bloedite', 'Borax', 'Boric_acid,s', 'Brucite', 'Burkeite', 'Calcite', 'Carnallite', 'Celestite', 'Chalcedony', 'Chrysotile', 'Diopside', 'Dolomite', 'Enstatite', 'Epsomite', 'Forsterite', 'Gaylussite', 'Glaserite', 'Glauberite', 'Goergeyite', 'Gypsum', 'Halite', 'Hexahydrite', 'Huntite', 'K2B4O7:4H2O', 'KB5O8:4H2O', 'Kainite', 'Kalicinite', 'Kieserite', 'Labile_S', 'Leonhardite', 'Leonite', 'Magnesite', 'MgCl2_2H2O', 'MgCl2_4H2O', 'Mirabilite', 'Misenite', 'NaB5O8:5H2O', 'NaBO2:4H2O', 'Nahcolite', 'Natron', 'Nesquehonite', 'Pentahydrite', 'Pirssonite', 'Polyhalite', 'Portlandite', 'Quartz', 'Schoenite', 'Sepiolite', 'Sepiolite(d)', 'SiO2(a)', 'Sylvite', 'Syngenite', 'Talc', 'Teepleite', 'Thenardite', 'Trona']

            mineral_formulas = ['Ca2Mg[Si2O7]', 'CaSO4', '☐Mg2Mg5Si8O22(OH)2', 'Mg48Si34O85(OH)62', 'CaCO3', 'K2SO4', 'Mg2(CO3)(OH)2·3H2O', 'BaSO4', 'MgCl2·6H2O', 'Na2Mg(SO4)2·4H2O', 'Na2B4O5(OH)4•8(H2O)', 'B(OH)3', 'Mg(OH)2', 'Na6(CO3)(SO4)2', 'CaCO3', 'KMgCl3•6(H2O)', 'SrSO4', 'SiO2', 'Mg3Si2O5(OH)4', 'CaMgSi2O6', 'CaMg(CO3)2', 'MgSiO3', 'MgSO4•7(H2O)', 'Mg2SiO4', 'Na2Ca(CO3)2•5(H2O)', 'NaK3(SO4)2', 'Na2Ca(SO4)2', 'K2Ca5(SO4)6•(H2O)', 'CaSO4•2(H2O)', 'NaCl', 'MgSO4•6(H2O)', 'CaMg3(CO3)4', 'K2B4O7•4H2O', 'KB5O8•4H2O', 'MgSO4•KCl•3(H2O)', 'KHCO3', 'MgSO4•(H2O)', 'Na4Ca(SO4)3•2H2O', 'MgSO4•4(H2O)', 'K2Mg(SO4)2•4(H2O)', 'MgCO3', 'MgCl2:2H2O', 'MgCl2•4H2O', 'Na2SO4•10(H2O)', 'K8H6(SO4)7', 'NaB5O8•5H2O', 'NaBO2•4H2O', 'NaHCO3', 'Na2CO3•10(H2O)', 'Mg(HCO3)(OH)•2(H2O)', 'MgSO4•5(H2O)', 'Na2Ca(CO3)2•2(H2O)', 'K2Ca2Mg(SO4)4•2(H2O)', 'Ca(OH)2', 'SiO2', 'K2Mg(SO4)2•6(H2O)', 'Mg4Si6O15(OH)2•6(H2O)', 'Mg4Si6O15(OH)2•6(H2O)', 'SiO2', 'KCl', 'K2Ca(SO4)2•(H2O)', 'Mg3Si4O10(OH)2', 'Na2B(OH)4Cl', 'Na2SO4', 'Na3(CO3)(HCO3)•2(H2O)']     


        elif database_selection == 'phreeqc':

            minerals = ['Al(OH)3(a)', 'Albite', 'Anhydrite', 'Anorthite', 'Aragonite', 'Barite', 'Ca-Montmorillonite', 'Calcite', 'Celestite', 'Chalcedony', 'Chlorite(14A)', 'Chrysotile', 'Dolomite', 'Fe(OH)3(a)', 'FeS(ppt)', 'Fluorite', 'Gibbsite', 'Goethite', 'Gypsum', 'Halite', 'Hausmannite', 'Hematite', 'Hydroxyapatite', 'Illite', 'K-feldspar', 'K-mica', 'Kaolinite', 'Mackinawite', 'Manganite', 'Pyrite', 'Pyrochroite', 'Pyrolusite', 'Quartz', 'Rhodochrosite', 'Sepiolite', 'Sepiolite(d)', 'SiO2(a)', 'Siderite', 'Strontianite', 'Sulfur', 'Sylvite', 'Talc', 'Vivianite', 'Witherite']

            mineral_formulas = ['Al(OH)3', 'NaAlSi3O8', 'CaSO4', 'CaAl2Si2O8', 'CaCO3', 'BaSO4', 'Ca0.165Al2.33Si3.67O10(OH)2', 'CaCO3', 'SrSO4', 'SiO2', 'Mg5Al2Si3O10(OH)8', 'Mg3Si2O5(OH)4', 'CaMg(CO3)2', 'Fe(OH)3', 'FeS', 'CaF2', 'Al(OH)3', 'FeO(OH)', 'CaSO4•2(H2O)', 'NaCl', 'Mn(II)Mn(III)2O4', 'Fe2O3', 'Ca5(PO4)3(OH)', 'K0.6Mg0.25Al2.3Si3.5O10(OH)2', 'KAlSi3O8', 'KAl3Si3O10(OH)2', 'Al2Si2O5(OH)4', 'FeS', 'MnO(OH)', 'FeS2', 'Mn(OH)2', 'MnO2', 'SiO2', 'MnCO3', 'Mg4Si6O15(OH)2•6(H2O)', 'Mg4Si6O15(OH)2•6(H2O)', 'SiO2', 'FeCO3', 'SrCO3', 'S8', 'KCl', 'Mg3Si4O10(OH)2', 'Fe3(PO4)2•8(H2O)', 'BaCO3']    


        short_mineral_names = ['Barite', 'Gypsum','Halite', 'Natron', 'Quartz', 'Talc','Trona', 'Borax', 'Albite', 'K-mica','Illite', 'Pyrite', 'Sulfur',]

        long_mineral_names = ['Anthophyllite', 'Hexahydrite', 'Leonhardite', 'Nesquehonite', 'Pentahydrite', 'Portlandite','Sepiolite(d)', 'Boric_acid,s', 'K2B4O7:4H2O', 'NaB5O8:5H2O', 'Rhodochrosite', 'Strontianite','Hydroxyapatite', 'Chlorite(14A)', 'Mackinawite', 'Hausmannite', 'Pyrochroite']

        cleaned_elements = []
        for element in elements:
            element = re.sub('(\(.\))', '', element)
            element.strip()
            cleaned_elements.append(element)

        remaining_characters = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '(', ')', 'I', '•']
        permitted_auxillary_reagents = ['OH', 'H2O','O']

        undescribed_minerals = []
        described_minerals = []
        formulas = []
        described_mineral_formulas = []
        for mineral in minerals:
            mineral_index = minerals.index(mineral)
            mineral_formula = mineral_formulas[mineral_index]
            formulas.append(mineral_formula)

            for element in cleaned_elements:
                if re.search(element, mineral_formula):
                    mineral_formula = re.sub(element, '', mineral_formula)

            for reagent in permitted_auxillary_reagents:
                if re.search(reagent, mineral_formula):
                    mineral_formula = re.sub(reagent, '', mineral_formula)

            if all(character in remaining_characters for character in mineral_formula):
                described_minerals.append(mineral)
                described_mineral_formulas.append(formulas[-1])


        #============================================================================
        # create the initial line of the equilibrium phases block
        equilibrium_phases_block = []
        announcement = '\nParameterize the scaling minerals:'
        print(announcement, '\n', '='*len(announcement))

        if self.parameters['simulation_type'] == 'transport':
            equilibrium_phases_number = input('''- Over what cell(s) will scaling occur? < # > or < #-# > 
            Default of all = 1 - %s\t''' %(self.parameters['cells_per_module']*quantity_of_modules)) or '1-{}'.format(self.parameters['cells_per_module']*quantity_of_modules)
        elif self.parameters['simulation_type'] == 'evaporation':
            equilibrium_phases_number = '1'

        equilibrium_phases_description = input('''- Do the scaling conditions have a special description?
        Default is none\t''')

        equilibrium_phases_line = '\nEQUILIBRIUM_PHASES %s\t%s' % (equilibrium_phases_number, 
                                                                   equilibrium_phases_description)
        equilibrium_phases_block.append(equilibrium_phases_line)

        # determine which predefined set of minerals should be simulated
        mineral_selections = ['All', 'All but a few', 'Build from scratch']
        print('''\nMineral options for the < %s > database and the parameters from < %s > elements:
    mineral\t\t\tcompound''' %(database_selection, water_selection))  
        print('--------------------------------------')
        for possible_mineral in described_minerals:
            if possible_mineral in short_mineral_names:
                print('%s' %(possible_mineral), '\t\t\t', described_mineral_formulas[described_minerals.index(possible_mineral)])  
            elif possible_mineral == 'Ca-Montmorillonite':
                print('%s' %(possible_mineral), '\t', described_mineral_formulas[described_minerals.index(possible_mineral)])  
            else:
                print('%s' %(possible_mineral), '\t\t', described_mineral_formulas[described_minerals.index(possible_mineral)])
        for selection in mineral_selections:
            print('\n< %s >' %(selection))

        mineral_selection = input('''- Which selection of minerals would you like to simulate?
        Default = 'All'  __ ''')
        while mineral_selection not in mineral_selections:
            print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
            mineral_selection = input('- Which selection of minerals would you like to simulate?   __ ')

        #=============================================================================
        # establish and edit all minerals for a simulation
        if mineral_selection in ['All','All but a few']:
            if mineral_selection == 'All':
                mineral_saturation_states = []
                minerals_initial_moles = []
                selected_minerals = described_minerals
                i = 0
                while i <= len(minerals):
                    mineral_saturation_states.append(0)
                    minerals_initial_moles.append(0)
                    i += 1

                possible_saturation_states = input('''- Do any minerals possess a non-zero saturation index? < y > or < n > 
        Default = < n >  __ ''') or 'n'
                while possible_saturation_states not in possible_answers:
                    print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                    possible_saturation_states = input('''- Do any minerals possess a non-zero saturation index? < y > or < n > 
        Default = < n >  __ ''') or 'n'
                if possible_saturation_states == 'y':
                    print('\nYour selected minerals:')
                    for mineral in minerals:
                        print('< %s >' %(mineral))

                    nonzero_mineral_saturation = input('''- Which mineral possesses an existing saturation index?
                    Type < done > when you have completed the parameterization.   __ ''') 
                    while nonzero_mineral_saturation != 'done' and nonzero_mineral_saturation not in selected_minerals:
                        print('''ERROR: The entered mineral is not in the list of accepted minerals. Please provide an accepted value''')
                        nonzero_mineral_saturation = input('- Which mineral possesses an existing saturation index?   __ ')   
                    while nonzero_mineral_saturation != 'done':   
                        mineral_saturation_state = input('- What is the saturation index?   __ ')
                        while not mineral_saturation_state.isnumeric():
                            print('''ERROR: The entered value is not an accepted number. Please provide an accepted value''')
                            mineral_saturation_state = input('- What is the saturation index?   __ ')

                        mineral_saturation_states[minerals.index(nonzero_mineral_saturation)] = mineral_saturation_state
                        nonzero_mineral_saturation = input('- Which mineral possesses an existing saturation index?   __ ')
                        while nonzero_mineral_saturation != 'done' and nonzero_mineral_saturation not in selected_minerals:
                            print('''ERROR: The entered mineral is not in the list of accepted minerals. Please provide an accepted value''')
                            nonzero_mineral_saturation = input('- Which mineral possesses an existing saturation index? __ ')   

                elif possible_saturation_states == 'n':
                    pass


                preexisting_minerals = input('''- Do any minerals currently exist in the module? 
        Default = 'n'  ; < y > or < n > __ ''')
                if preexisting_minerals == 'y':
                    print('\nYour selected minerals:')
                    for mineral in selected_minerals:
                        print(mineral)

                    existing_mineral = input('''- Which mineral possesses an existing saturation index?
                    Type < done > when you have finished parameterizing the simulation. __ ''')    
                    while existing_mineral != 'done':  
                        while existing_mineral != 'done' and existing_mineral not in selected_minerals:
                            print('''ERROR: The entered mineral is not in the list of accepted minerals. Please provide an accepted value''')
                            existing_mineral = input('- Which mineral possesses an existing saturation index?  __ ')

                        mineral_inital_moles = input('- What is the initial concentration?  __ ')
                        while not mineral_inital_moles.isnumeric():
                            print('''ERROR: The entered value is not an accepted number. Please provide an accepted value''')
                            mineral_inital_moles = input('What is the initial concentration?  __ ')

                        minerals_initial_moles[selected_minerals.index(existing_mineral)] = mineral_inital_moles
                        existing_mineral = input('- Which mineral possesses an existing saturation index?')

                elif preexisting_minerals == 'n':
                    pass

                else:
                    while preexisting_minerals not in possible_answers:
                        print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                        preexisting_minerals = input('''- Do any minerals currently previously exist in the module? < y > or < n > __ ''')


        #=============================================================================
            # establish and systematically remove all but a few minerals for simulation
            elif mineral_selection == 'All but a few':
                mineral_saturation_states = []
                minerals_initial_moles = []
                selected_minerals = described_minerals
                i = 0
                while i <= len(selected_minerals):
                    mineral_saturation_states.append(0)
                    minerals_initial_moles.append(0)
                    i += 1
                print('\nAll possible minerals with the %s database:' %(database_selection))
                for mineral in selected_minerals:
                    print('< %s >' %(mineral))

                removed_mineral = input('''- Which mineral should be removed?
                    One must be selected at a time. Enter < done > when you are finished.  __ ''')
                while removed_mineral != 'done':
                    while removed_mineral != 'done' and removed_mineral not in selected_minerals:
                        print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                        removed_mineral = input('''- Which mineral should be removed?
                        Only one can be removed at a time. Enter < done > when you are finished.  __ ''')

                    removed_mineral_index = selected_minerals.index(removed_mineral)
                    del selected_minerals[removed_mineral_index], mineral_saturation_states[removed_mineral_index], minerals_initial_moles[removed_mineral_index]
                    removed_mineral = input('- Which mineral should be removed?   __ ')

                # include non-zero saturation indices    
                possible_saturation_states = input('''- Do any minerals possess a non-zero saturation index? < y > or < n > __ ''')
                if possible_saturation_states == 'y':
                    print('\nYour selected minerals:')
                    for mineral in selected_minerals:
                        print('< %s >' %(mineral))

                    nonzero_mineral_saturation = input('''- Which mineral possesses an existing saturation index?
                    One must be selected at a time. Enter < done > when you are finished.   __ ''')    
                    while nonzero_mineral_saturation != 'done':
                        while removed_mineral != 'done' and nonzero_mineral_saturation not in selected_minerals:
                            print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                            removed_mineral = input('''Which mineral possesses an existing saturation index?
                            Only one can be removed at a time.Enter < done > when you are finished.  __ ''')

                        mineral_saturation_state = input('- What is the saturation index?  __ ')
                        mineral_saturation_states[selected_minerals.index(nonzero_mineral_saturation)] = mineral_saturation_state
                        nonzero_mineral_saturation = input('- Which mineral possesses an existing saturation index?   __ ')

                while possible_saturation_states not in possible_answers:
                    print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                    possible_saturation_states = input('''- Do any minerals possess a non-zero saturation index?
                    < y > or < n > __ ''')

                # include previously existing mineral moles
                preexisting_minerals = input('''- Do any minerals previously exist in the module?
                < y > or < n >   __ ''')
                if preexisting_minerals == 'y':
                    print('\nYour selected minerals:')
                    for mineral in selected_minerals:
                        print(mineral)

                    existing_mineral = input('''- Which mineral possesses an existing saturation index?
                    One must be selected at a time.Enter < done > when you are finished.   __ ''')  
                    while existing_mineral != 'done':
                        while existing_mineral != 'done' and existing_mineral not in selected_minerals:
                            print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                            removed_mineral = input('''- Which mineral possesses an existing saturation index?
                            Only one can be removed at a time.Enter < done > when you are finished.   __ ''')

                        mineral_inital_moles = input('- What is the saturation index?')
                        minerals_initial_moles[selected_minerals.index(existing_mineral)] = mineral_inital_moles
                        existing_mineral = input('- Which mineral possesses an existing saturation index?')     

                while preexisting_minerals not in possible_answers:
                    print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                    preexisting_minerals = input('''- Do any minerals previously exist in the module? < y > or < n > __ ''')

            # add the elemental lines to the EQUILIBRIUM_PHASES block of the code
            for mineral_index in range(len(selected_minerals)):
                if selected_minerals[mineral_index] in short_mineral_names:
                    equilibrium_phases_line = '%s\t\t%s\t%s' %(selected_minerals[mineral_index],
                                                               mineral_saturation_states[mineral_index], 
                                                               minerals_initial_moles[mineral_index])
                elif selected_minerals[mineral_index] not in short_mineral_names:
                    equilibrium_phases_line = '%s\t%s\t%s' %(selected_minerals[mineral_index], 
                                                             mineral_saturation_states[mineral_index], 
                                                             minerals_initial_moles[mineral_index]) 

                equilibrium_phases_block.append(equilibrium_phases_line)        

        #=============================================================================
        # establish a simulation with user-defined minerals
        elif mineral_selection == 'Build from scratch':
            mineral_saturation_states = []
            minerals_initial_moles = []
            selected_minerals = []

            print('''Mineral options for the < %s > database: \n< mineral >\t\tcompound''' %(database_selection))
            for possible_mineral in described_minerals:
                if possible_mineral in long_mineral_names or possible_mineral == 'Ca-Montmorillonite':
                    print('< %s >' %(possible_mineral), '\t', mineral_formulas[described_minerals.index(possible_mineral)])
                elif possible_mineral in short_mineral_names:
                    print('< %s >' %(possible_mineral), '\t\t', sorted_mineral_formulas[described_minerals.index(possible_mineral)])  
                else:
                    print('< %s >' %(possible_mineral), '\t\t', sorted_mineral_formulas[described_minerals.index(possible_mineral)])


            mineral = input('''- What mineral will you simulate?
            Type < done > when you have parameterized the last mineral   __ ''')
            while mineral not in minerals:
                print('''The mineral is not supported by the %s database.A Solution_Master must be defined for the %s mineral''' %(database_selection, mineral))  
                mineral = input('- What mineral will you simulate?\t')
            while mineral != 'done':
                selected_minerals.append(mineral)

                #include non-zero saturation indices    
                mineral_saturation_state = input('''- What is the saturation index of the mineral?
                Default = 0   __ ''') or str(0)
                while not mineral_saturation_state.isnumeric():
                    print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                    mineral_saturation_state = input('''- What is the saturation index of the mineral?
                    Default = 0   __ ''')

                #include previously existing mineral moles
                preexisting_mineral_moles = input('''- What is the existing quantity of moles for the mineral in the module?
                Default = 0   __ ''') or str(0)
                while not preexisting_mineral_moles.isnumeric():
                    print('''ERROR: The entry is beyond the accepted values. Please provide an accepted value''')
                    preexisting_mineral_moles = input('''- What is the existing quantity of moles for the mineral in the module?
                    Default = 0   __ ''') or str(0)

                #add the mineral information to the EQUILIBRIUM_PHASES block
                if mineral in short_mineral_names:
                    equilibrium_phases_line = '%s\t\t%s\t%s' %(mineral, 
                                                               mineral_saturation_state,
                                                               preexisting_mineral_moles)
                elif mineral not in short_mineral_names:
                    equilibrium_phases_line = '%s\t%s\t%s' %(mineral,
                                                             mineral_saturation_state,
                                                             preexisting_mineral_moles) 
                equilibrium_phases_block.append(equilibrium_phases_line)
                mineral = input('- What mineral will you simulate?   __ ')

        else:
            print('ERROR: The passed EQUILIBRIUM_PHASES scope is beyond the alotted options.')


    def make_selected_output():
        '''
        Specify the output file after a PHREEQC simulation
        '''
        global selected_output_block
        global selected_output_file_name

        announcement = '\nParameterize the output file:'
        print(announcement, '\n', '='*len(announcement))
        first_line = '\nSELECTED_OUTPUT'
        print('''The inclusion of "Scaling" or "Brine" and "phreeqc" or "pitzer" in the selected output file name facilitates data processing. 
    Exclude the extension of the file name.''')
        count = 0
        proposed_selected_output_file_name = '{}_{}, {}, {}, {}_{}'.format(datetime.date.today(), water_selection, self.parameters['simulation_type'], database_selection, self.parameters['output_perspective'], count) 
        while os.path.exists('{}.txt'.format(proposed_selected_output_file_name)):
            count += 1
            proposed_selected_output_file_name = '{}_{}, {}, {}_{}'.format(datetime.date.today(), water_selection, self.parameters['simulation_type'], database_selection, self.parameters['output_perspective'], count) 

        selected_output_file_name = input('''- What is the name for your output file? 
        Default = {}  __'''.format(proposed_selected_output_file_name)) or proposed_selected_output_file_name

        if not re.search('scaling|brine', selected_output_file_name, flags=re.IGNORECASE):
            selected_output_file_name += ', {}'.format(self.parameters['output_perspective'])

        selected_output_file_name += '.txt'
        file_name_line = '-file\t\t\t%s' %(selected_output_file_name)

        reaction_line = '-reaction\t\ttrue'
        temperature_line = '-temperature\t\ttrue'

        minerals_line = ''
        for mineral in selected_minerals:
             minerals_line += ' %s' %(mineral)
        saturation_indices_line = '-saturation_indices\t%s' %(minerals_line)
        equilibrium_phases_line = '-equilibrium_phases\t%s' %(minerals_line)

        elements_line = ''
        for element in elements:
             elements_line += ' %s' %(element)
        total_elements_line = '-totals\t\t\t' + elements_line    

        ph_line = '-pH\t\t\ttrue'
        solution_line = '-solution'
        time_line = '-time\t\t\ttrue'
        distance_line = '-distance\t\ttrue'
        simulation_line = '-simulation\t\ttrue'
        high_precision_line = '-high_precision\ttrue'
        step_line = '-step'
        water_line = '-water'

        selected_output_block = []
        selected_output_block.extend((first_line, file_name_line, reaction_line, temperature_line, total_elements_line, saturation_indices_line, equilibrium_phases_line, ph_line, time_line, distance_line, simulation_line, high_precision_line, solution_line, step_line,water_line))


    def export():
        """
        View and export the PHREEQC input file
        """
        global complete_lines
        global input_file_name

        complete_lines = chain(general_conditions, solution_block, equilibrium_phases_block, self.results['reaction_block'], selected_output_block, self.results['transport_block']) 

        # create the export input file
        file_number = 0
        proposed_file_name = '%s_PHREEQC input file, %s, %s_%s.pqi' %(datetime.date.today(), water_selection, 
                                                                        self.parameters['output_perspective'], file_number)
        if not os.path.exists(proposed_file_name):
            input_file_name = proposed_file_name

        elif os.path.exists(proposed_file_name):
            while os.path.exists('%s_PHREEQC input file, %s, %s, %s_%s.pqi' %(datetime.date.today(), water_selection, self.parameters['simulation_type'], self.parameters['output_perspective'], file_number)):
                file_number += 1
                input_file_name = '%s_PHREEQC input file, %s, %s, %s_%s.pqi' %(datetime.date.today(), water_selection, self.parameters['simulation_type'], self.parameters['output_perspective'], file_number)

        # printing and exporting the input file
        printing_block = input('''- Would you like to view your input file? < y > or < n > __ ''')
        if printing_block == 'y':
            print('\n\nInput file:\n')

        input_file = open(input_file_name,'a')
        for line in complete_lines:
            if printing_block == 'y':
                print(line)
            input_file.write(line + '\n')

        input_file.close()


    def execute():
        '''
        Execute the input file through the batch PHREEQC software 
        '''
        global database_selection
        global input_file_name

        announcement = '\nExecute the input file:'
        print(announcement, '\n', '='*len(announcement))
        #subprocess.call("start powershell", shell=True)
        working_directory = os.getcwd()
        phreeqc_path = 'C:\\Program Files\\USGS\\phreeqc-3.6.2-15100-x64'
        if input_selection == 'n':
            database_selection = input('''- What database do you select?
            < pitzer > or < phreeqc >  __ ''')
            input_file_name = input('- What is the input file name?')
            input_path = os.path.join(working_directory, input_file_name)
            while not os.path.exists(input_path):
                print('''ERROR: The input file path {} does not exist. Provide a valid input file path.'''.format(input_path))
                input_file_name = input('What is the input file name?')
                input_path = os.path.join(working_directory, input_file_name)

            output_file_name = re.sub('(?<=\.)(.+)' , 'pqo', input_file_name)
            database_path = ''

        else:
            input_path = os.path.join(working_directory, input_file_name)
            output_file_name = re.sub('pqi', 'pqo', input_file_name)
            database_path = os.path.join(phreeqc_path, 'database\\%s.dat' %(database_selection))        

        output_path = os.path.join(working_directory, output_file_name)
        bat_path = os.path.join(phreeqc_path, 'bin\\phreeqc.bat')
        print('\ninput file path: {}'.format(input_path))

        proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
        command = str.encode("\"" + bat_path + "\" \"" + input_path + "\" \"" + output_path + "\"  \"" + database_path + "\"\n") 
        proc.stdin.write(command)
        proc.stdin.close()  
        proc.wait()

        if os.path.exists(output_path):
            if input_selection == 'y':
                if os.path.exists(selected_output_file_name):
                    print('The execution is complete.')
                else:
                    print('ERROR: The SELECTED_OUTPUT file is missing. The simulation failed to simulate.')
            else:
                print('The execution is complete.')
        else:
            print('\nERROR: The simulation failed to execute.')


    def process_selected_output():
        """
        Interpreting the PHREEQC SELECTED_OUTPUT file and conducting the plotting functions
        """
        global selected_output_file_name
        global self.parameters['simulation_type']
        global simulation_perspective
        global initial_solution_mass
        global selected_output_file
        global graphical_selection
        global self.parameters['output_perspective']
        global simulation_cf
        global csv_data   
        global database


        announcement = '\nProcess the selected_output file:'
        print(announcement, '\n', '='*len(announcement))

        # determining the appropriate variables
        if input_selection == 'y':
            selected_output_file = selected_output_file_name
            simulation_perspective = self.parameters['output_perspective']
            database = database_selection 

        else:       
            if execute_selection == 'y':
                input_file = open(input_file_name, 'r')
                for line in input_file:
                    if re.search('(-file\s+)', line):
                        selected_output_file = re.sub('(-file\s+)', '', line)
                        selected_output_file = re.sub('(\n)', '', selected_output_file)

                    if re.search('(TRANSPORT)', line):
                        self.parameters['simulation_type'] = 'transport'
                    else:
                        self.parameters['simulation_type'] = 'evaporation'

            else:  
                selected_output_file = input('''- What is the name and\or path of the simulation file?
                Include the extension, like < .txt > __ ''')
                selected_output_file_name = re.sub('(?<=\.)(.+)', '', selected_output_file)
                while not os.path.exists(selected_output_file):
                    print('ERROR: The simulation file is missing from the current directory.')  
                    selected_output_file = input('What is the name and\or path of the simulation file?')  

                for line in selected_output_file:
                    if re.search('(TRANSPORT)', line):
                        self.parameters['simulation_type'] = 'transport'
                    else:
                        self.parameters['simulation_type'] = 'evaporation'

            #determining the scope of the simulation data
            if re.search('(Scaling)', selected_output_file, flags=re.IGNORECASE):
                simulation_perspective = 'Scaling'

            elif re.search('(Brine)', selected_output_file, flags=re.IGNORECASE):
                simulation_perspective = 'Brine'

            else:     
                simulation_perspective = input('''- Is the output file representative of a simulation for scaling or brine?
                < Scaling > or < Brine >''')
                while simulation_perspective != 'Brine' and simulation_perspective != 'Scaling':
                    print('''ERROR: The value is not one of the options. Select one of the choices to proceed.''')  
                    simulation_perspective = input('- Is the output file representative of a simulation for scaling or brine?')

            self.parameters['output_perspective'] = ''

        graphical_selection = input('''- Would you like to view the effluent brine or the module scaling from your {} simulation?
        < Brine > or < Scaling >'''.format(self.parameters['output_perspective']))
        while graphical_selection != 'Brine' and graphical_selection != 'Scaling':
            print('''ERROR: The value is not one of the options. Select one of the choices to proceed.''')  
            graphical_selection = input('''- Would you like to view brine over time in the module < Brine > , 
            or would your like to view scaling over distance in the module < Scaling > ? __ ''')

        # preparing the SELECTED_OUTPUT file into a dataframe
        selected_output = open(selected_output_file, 'r')
        original_data = pandas.read_table(selected_output, sep = '\t')
        csv_data = pandas.DataFrame(original_data)
        for column in csv_data.columns:
            new_column = column.strip()
            csv_data.rename(columns={column:new_column}, inplace = True)

        initial_solution_mass = csv_data.at[0, 'mass_H2O']
        final_solution_mass = csv_data['mass_H2O'].iloc[-1]
        simulation_cf = initial_solution_mass / final_solution_mass

        # conducting the appropriate visualization function
        if graphical_selection == 'Brine':
            make_brine_plot()
        elif graphical_selection == 'Scaling':
            make_scaling_plot()
        else:
            print('''ERROR: Reconfigure the SELECTED_OUTPUT simulation file name and\or the above visualizaiton parameters.''')


    def make_brine_plot():
        """
        Generate plots of the elemental concentrations from effluent brine in the PHREEQC SELECTED_OUTPUT file  
        """
        global export_name
        global plot_title 
        global figure

        elements = []
        for column in csv_data.columns:
            if re.search('([A-Z][a-z]?(\(\d\))?){1}$', column) and not re.search('(_|H2O|pH)', column):
                elements.append(column)

        csv_data.drop(csv_data.index[:3], inplace=True)

        # plot the brine concentrations figure
        unit = 'mol/kgw'
        concentration_array = []
        pyplot.figure(figsize = (17,10))
        plot_title = (input('''- What is the title of the plot?
        Default = Effluent brine elemental concentrations ____ ''')) or 'Effluent brine elemental concentrations'
        pyplot.title(plot_title, fontsize = 'xx-large')
        pyplot.xlabel('time (s)', fontsize = 'x-large')
        pyplot.ylabel('concentration (%s)' %(unit), fontsize = 'x-large')
        pyplot.grid(True)

        for element in elements:  
            concentration_serie = []
            time_serie = []
            initial_solution_time = 0
            for index, row in csv_data.iterrows():
                if csv_data.at[index, 'Cl'] == 0:
                    initial_solution_time += 1
                    #print('yes')
                else:
                    concentration_serie.append(csv_data.at[index, element])
                    time_serie.append(csv_data.at[index, 'time'] - initial_solution_time * (csv_data.at[index, 'time'] - csv_data.at[index-1, 'time']))

            pyplot.plot(time_serie,concentration_serie)

        plot_caption = '''\n\nBrine Figure:\n%s 
        The effluent concentrations of each existing element in the brine. Brine plots from brine data rapidly reach a steady state elemental concentration. Brine plots from scaling data consist of vertical concentrations that represent the concentrations at each distance throughout the RO module at the specified time, where the low end represents the influent concentration while the high end represents the effluent concentration.''' %('='*len('Brine Figure'))

        pyplot.legend(elements, loc='best', title = 'non-zero elements', fontsize = 'x-large')
        pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(simulation_cf), wrap=True, horizontalalignment='left', fontsize=12)
        pyplot.yscale('log')
        figure = pyplot.gcf()
        print('\nClose the figure to proceed.')
        print(plot_caption)
        pyplot.show()
        export_figure = input('''Will you export the figure? < n > or < y >  __  ''')
        while export_figure not in possible_answers:
            print('ERROR: The entered mineral is not among the options.')  
            export_figure = input('''\nWill you export the figure? < n > or < y >  __  ''')

        # create a complementary concentration data table for the scaling figure 
        loop_iteration = 1
        table_view = {}
        average_concentrations_table = pandas.DataFrame()
        index_elements = []
        if simulation_perspective == 'Scaling':
            for element in elements:
                quantity_of_steps_index = 0
                average_iteration = 0
                time_serie = []            
                time_averages = {}
                for index, row in csv_data.iterrows():
                    if csv_data.at[index, 'time'] == 0:
                        time_serie.append(csv_data.at[index,element])                 
                        quantity_of_steps_index += 1                    

                    elif csv_data.at[index-1,'soln'] == quantity_of_steps_index:       
                        #process the complete time serie
                        try:
                            average_concentration = sum(time_serie) / len(time_serie)
                        except:
                            average_concentration = 0

                        #print(average_concentration)
                        table_view['Time (s): %s' %(average_iteration * quantity_of_steps_index)] = average_concentration
                        average_iteration += 1
                        #print('mid-End')   

                        #begin the new time serie
                        time_serie = []
                        time_averages = {}
                        time_serie.append(csv_data.at[index,element])
                        #print('middle')    

                    elif index == len(csv_data[element]) + 2:       
                        time_serie.append(csv_data.at[index,element])            
                        try:
                            average_concentration = sum(time_serie) / len(time_serie)
                        except:
                            average_concentration = 0
                        table_view['Time (s): %s' %(round(average_iteration * quantity_of_steps_index), 1)] = average_concentration
                        average_iteration += 1
                        #print('end-End')  
                        index_elements.append(element)
                        average_concentrations_table = average_concentrations_table.append(table_view, ignore_index=True)
                        loop_iteration += 1

                    else:
                        #print('middle')      
                        time_serie.append(csv_data.at[index,element])

            average_concentrations_table.index = index_elements
            average_concentrations_table.index.name = 'Elements'
            dataframe_title = 'Average elemental molal concentrations of the feed water in the RO module for each %s seconds of simulation:' %(quantity_of_steps_index)


        # create a complementary concentration data table for the brine figure 
        elif simulation_perspective == 'Brine':
            total_time = csv_data['time'].iloc[-1]
            for element in elements:  
                concentration_serie = []
                time_serie = []
                for index, row in csv_data.iterrows():
                    if csv_data.at[index, 'Cl'] != 0:
                        concentration_serie.append(csv_data.at[index,element])      
                average_concentration = sum(concentration_serie) / len(concentration_serie)
                table_view['%s' %(element)] = average_concentration

            average_concentrations_table = average_concentrations_table.append(table_view, ignore_index=True)
            average_concentrations_table.rename(index = {0:'Concentrations (molal)'}, inplace = True)
            dataframe_title = 'Average elemental molal concentrations of the feed water in the RO module over %s seconds of simulation:' %(total_time)

        print('\n\n\n',dataframe_title,'\n%s'%('='*len(dataframe_title)))
        print(average_concentrations_table)

        # export the output graphic
        if export_figure == 'y':
            export_filename_progenitor = re.sub('(\.\w+)', '', selected_output_file_name)
            if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
                export_name = '{}, {}'.format(export_filename_progenitor, simulation_perspective)
            else:
                export_name = export_filename_progenitor

            export_option = input('''- Would you like to export the figure? < y > or < n > __ ''')
            if export_option == 'y':
                export_plot()
            while export_option not in possible_answers:
                print('''ERROR: The value is not one of the options.''')  
                export_option = input('''- Would you like to export the figure?
                < y > or < n >''')

    def make_scaling_plot():
        """
        Generate plots of scaling along the module distance in the PHREEQC SELECTED_OUTPUT file  
        """
        global initial_solution_mass
        global individual_plots
        global mineral_formulas
        global export_name
        global plot_title
        global mineral
        global minerals
        global figure

        if input_selection == 'n':
            minerals = ['Akermanite', 'Al(OH)3(a)', 'Albite', 'Anhydrite', 'Anorthite', 'Anthophyllite', 'Antigorite', 'Aragonite', 'Arcanite', 'Artinite', 'Barite', 'Bischofite', 'Bloedite', 'Borax', 'Boric_acid,s', 'Brucite', 'Burkeite', 'Ca-Montmorillonite', 'Calcite', 'Carnallite', 'Celestite', 'Chalcedony', 'Chlorite(14A)', 'Chrysotile', 'Diopside', 'Dolomite', 'Enstatite', 'Epsomite', 'Fe(OH)3(a)', 'FeS(ppt)', 'Fluorite', 'Forsterite', 'Gaylussite', 'Gibbsite', 'Glaserite', 'Glauberite', 'Goergeyite', 'Goethite', 'Gypsum', 'Halite', 'Hausmannite', 'Hematite', 'Hexahydrite', 'Huntite', 'Hydroxyapatite', 'Illite', 'K-feldspar', 'K-mica', 'K2B4O7:4H2O', 'KB5O8:4H2O', 'Kainite', 'Kalicinite', 'Kaolinite', 'Kieserite', 'Labile_S', 'Leonhardite', 'Leonite', 'Mackinawite', 'Magnesite', 'Manganite', 'MgCl2_2H2O', 'MgCl2_4H2O', 'Mirabilite', 'Misenite', 'NaB5O8:5H2O', 'NaBO2:4H2O', 'Nahcolite', 'Natron', 'Nesquehonite', 'Pentahydrite', 'Pirssonite', 'Polyhalite', 'Portlandite', 'Pyrite', 'Pyrochroite', 'Pyrolusite', 'Quartz', 'Rhodochrosite', 'Schoenite', 'Sepiolite', 'Sepiolite(d)', 'SiO2(a)', 'Siderite', 'Strontianite', 'Sulfur', 'Sylvite', 'Syngenite', 'Talc', 'Teepleite', 'Thenardite', 'Trona', 'Vivianite', 'Witherite']

            mineral_formulas = ['Ca2Mg[Si2O7]', 'Al(OH)3', 'NaAlSi3O8', 'CaSO4', 'CaAl2Si2O8', '☐Mg2Mg5Si8O22(OH)2', 'Mg48Si34O85(OH)62', 'CaCO3', 'K2SO4', 'Mg2(CO3)(OH)2·3H2O', 'BaSO4', 'MgCl2·6H2O', 'Na2Mg(SO4)2·4H2O', 'Na2B4O5(OH)4•8(H2O)', 'B(OH)3', 'Mg(OH)2', 'Na6(CO3)(SO4)2', 'Ca0.165Al2.33Si3.67O10(OH)2', 'CaCO3', 'KMgCl3•6(H2O)', 'SrSO4', 'SiO2', 'Mg5Al2Si3O10(OH)8', 'Mg3Si2O5(OH)4', 'CaMgSi2O6', 'CaMg(CO3)2', 'MgSiO3', 'MgSO4•7(H2O)', 'Fe(OH)3', 'FeS', 'CaF2', 'Mg2SiO4', 'Na2Ca(CO3)2•5(H2O)', 'Al(OH)3', 'NaK3(SO4)2', 'Na2Ca(SO4)2', 'K2Ca5(SO4)6•(H2O)', 'FeO(OH)', 'CaSO4•2(H2O)', 'NaCl', 'Mn(II)Mn(III)2O4', 'Fe2O3', 'MgSO4•6(H2O)', 'CaMg3(CO3)4', 'Ca5(PO4)3(OH)', 'K0.6Mg0.25Al2.3Si3.5O10(OH)2', 'KAlSi3O8', 'KAl3Si3O10(OH)2', 'K2B4O7•4H2O', 'KB5O8•4H2O', 'MgSO4•KCl•3(H2O)', 'KHCO3', 'Al2Si2O5(OH)4', 'MgSO4•(H2O)', 'Na4Ca(SO4)3•2H2O', 'MgSO4•4(H2O)', 'K2Mg(SO4)2•4(H2O)', 'FeS', 'MgCO3', 'MnO(OH)', 'MgCl2:2H2O', 'MgCl2•4H2O', 'Na2SO4•10(H2O)', 'K8H6(SO4)7', 'NaB5O8•5H2O', 'NaBO2•4H2O', 'NaHCO3', 'Na2CO3•10(H2O)', 'Mg(HCO3)(OH)•2(H2O)', 'MgSO4•5(H2O)', 'Na2Ca(CO3)2•2(H2O)', 'K2Ca2Mg(SO4)4•2(H2O)', 'Ca(OH)2', 'FeS2', 'Mn(OH)2', 'MnO2', 'SiO2', 'MnCO3', 'K2Mg(SO4)2•6(H2O)', 'Mg4Si6O15(OH)2•6(H2O)', 'Mg4Si6O15(OH)2•6(H2O)', 'SiO2', 'FeCO3', 'SrCO3', 'S8', 'KCl', 'K2Ca(SO4)2•(H2O)', 'Mg3Si4O10(OH)2', 'Na2B(OH)4Cl', 'Na2SO4', 'Na3(CO3)(HCO3)•2(H2O)', 'Fe3(PO4)2•8(H2O)', 'BaCO3']      

        # the complete list of all minerals is created
        csv_minerals = []
        for column in csv_data.columns:
            if re.search('([A-Z].{3,})', column) and not re.search('(\(|\_|\:)', column):
                csv_minerals.append(column)

        if self.parameters['simulation_type'] == 'transport':
            csv_data.drop(csv_data.index[:3], inplace=True)

        plot_title = (input('''- What is the title of the plot? 
            Default = Scaling throughout the RO module  ___ ''')) or 'Scaling throughout the RO module'



        # all of the non-zero minerals are identified and the chemical formulas are sorted into a list
        non_zero_minerals = []
        for mineral in csv_minerals:
            for row in csv_data[mineral]:
                if row != 0 and mineral not in non_zero_minerals:
                    non_zero_minerals.append(mineral)

        quantity_nonzero_minerals = len(non_zero_minerals)      

        non_zero_mineral_formulas = []
        for mineral in non_zero_minerals:
            mineral_index = minerals.index(mineral)
            mineral_formula = mineral_formulas[mineral_index]
            non_zero_mineral_formulas.append(mineral_formula)      

        # plot the simulation depending upon the simulation perspective
        unit = 'moles'
        if simulation_perspective == "Brine":
            individual_plots = 'n'
            pyplot.figure(figsize = (17,10))
            pyplot.title(plot_title, fontsize = 'xx-large')
            pyplot.xlabel('Simulation time (s)', fontsize = 'x-large')
            pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  
            experimental_loop = []
            formula_index = 0
            for mineral in non_zero_minerals:
                mineral_serie = []
                time_serie = []
                for index, row in csv_data.iterrows():
                    mineral_serie.append(csv_data.at[index, mineral]) 
                    time = csv_data.at[index, 'time']
                    time_serie.append(time)

                pyplot.plot(time_serie,mineral_serie)
                pyplot.scatter(time_serie,mineral_serie)

                experimental_loop.append('%s [%s]' %(mineral,non_zero_mineral_formulas[formula_index]))   

            pyplot.legend(experimental_loop, loc='best', fontsize = 'x-large')
            pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(simulation_cf), wrap=True, horizontalalignment='left', fontsize=12)
            figure = pyplot.gcf()
            print('\nClose the figure to proceed.\n')
            pyplot.show()

            export_filename_progenitor = re.sub('(\.\w+)', '', selected_output_file_name)
            if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
                export_name = '{}, {}'.format(export_filename_progenitor, simulation_perspective)
            else:
                export_name = export_filename_progenitor

            export_option = input('''- Would you like to export the figure? < y > or < n > __ ''')
            while export_option not in possible_options:
                print('''ERROR: The value is not one of the options.''')  
                export_option = input('''Would you like to export the figure? < y > or < n > __ ''')  

            if export_option == 'y':
                export_plot()

        elif simulation_perspective == 'Scaling':
            if quantity_nonzero_minerals < 2:
                individual_plots = 'n'
            elif quantity_nonzero_minerals >= 2:
                individual_plots = 'y'
            print('\nQuantity of precipitated minerals: %s' %(quantity_nonzero_minerals))
            #print('Number of timesteps per mineral: %s')
            individual_plots = input('''- Would you like to plot each mineral on a separate figure?
        Suggestion = < %s >  ;  < y > or < n >  __ '''  %(individual_plots)) or individual_plots
            while individual_plots not in possible_answers:
                print('ERROR: The entered value is not accepted.')  
                individual_plots = input('- Would you like to plot each mineral on a separate figure?')

            if individual_plots == 'y':
                viewing_plots = input('''How many mineral figures would you like to view?
        < All >, < A few >, or < None >''')
                while viewing_plots != 'All' and exporting_plots != 'A few' and exporting_plots != 'None':
                    print('ERROR: The entered value is not accepted.')  
                    viewing_plots = input('- How many mineral figures would you like to view?')

                if viewing_plots == ('All' or 'A few'):
                    viewing_minerals = []
                    if viewing_plots == 'All':
                        viewing_minerals = non_zero_minerals

                    elif exporting_plots == 'A few':
                        for mineral in non_zero_minerals:
                            print('< %s >' %(mineral))

                        viewing_mineral = input('''- Which minerals will you view?
            Type < done > when you are finished.''')
                        while viewing_mineral not in non_zero_minerals:
                            print('ERROR: The entered mineral is not among the options.')  
                            viewing_mineral = input('- Which minerals will you view?')
                        while viewing_mineral != 'done':
                            viewing_minerals.append(exporting_mineral)
                            viewing_mineral = input('- Which minerals will you view?')
                            while viewing_mineral not in non_zero_minerals:
                                print('ERROR: The entered mineral is not among the options.')  
                                viewing_mineral = input('- Which minerals will you view?')

                    for mineral in viewing_minerals:
                        pyplot.figure(figsize = (17,10))
                        pyplot.title(plot_title, fontsize = 'xx-large')

                        if self.parameters['simulation_type'] == 'transport':
                            pyplot.xlabel('Midpoint module distance (m)', fontsize = 'x-large')
                            pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  
                            experimental_loop = []

                            iteration = 0
                            distance_serie = []
                            time_serie = []
                            quantity_of_steps_index = 0   
                            for index, row in csv_data.iterrows():
                                if csv_data.at[index, 'time'] == 0:
                                    time_serie.append(csv_data.at[index, mineral]) 
                                    distance_serie.append(csv_data.at[index, 'dist_x'])
                                    quantity_of_steps_index += 1   
                                    time = 0

                                elif csv_data.at[index-1, 'soln'] == quantity_of_steps_index:
                                    experimental_loop.append('%s [%s] ; time (s): %s' 
                                                             %(mineral, 
                                                                mineral_formulas[minerals.index(mineral)], 
                                                                round(time, 2)))
                                    pyplot.plot(distance_serie,time_serie)
                                    distance_serie = []
                                    time_serie = []
                                    time_serie.append(csv_data.at[index, mineral])
                                    distance_serie.append(csv_data.at[index, 'dist_x'])
                                    time = csv_data.at[index, 'time']

                                elif index == len(csv_data[mineral]) + 2:   
                                    experimental_loop.append('%s [%s] ; time (s): %s' 
                                                             %(mineral, 
                                                                mineral_formulas[minerals.index(mineral)], 
                                                                round(time, 2)))
                                    pyplot.plot(distance_serie,time_serie)

                                else:
                                    time_serie.append(csv_data.at[index, mineral])
                                    distance_serie.append(csv_data.at[index, 'dist_x'])
                                    iteration += 1

                        elif self.parameters['simulation_type'] == 'evaporation':
                            pyplot.xlabel('Concentration Factor (CF)', fontsize = 'x-large')
                            pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  

                            experimental_loop = []
                            cf_series = []
                            concentration_series = []
                            data_length = len(csv_data['mass_H2O'])
                            for index, row in csv_data.iterrows():
                                if index < data_length:
                                    if csv_data.at[index, 'step'] >= 1:
                                        concentration_series.append(csv_data.at[index, mineral]) 
                                        solution_mass = csv_data.at[index, 'mass_H2O']
                                        cf_series.append(initial_solution_mass / solution_mass)  

                                    elif index > 1:
                                        print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')

                            concentration_series.append(csv_data.at[index, mineral]) 
                            solution_mass = csv_data.at[index, 'mass_H2O']
                            cf_series.append(initial_solution_mass / solution_mass)  

                            experimental_loop.append('%s [%s]' %(mineral, mineral_formulas[minerals.index(mineral)]))
                            pyplot.plot(cf_series,concentration_series)                    


                        pyplot.legend(experimental_loop, loc='best', fontsize = 'x-large')
                        pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(simulation_cf), wrap=True, horizontalalignment='left', fontsize=12)
                        figure = pyplot.gcf()
                        print('\nClose the figure to proceed.')
                        pyplot.show()
                        export_figure = input('''Will you export the figure? < n > or < y >  __  ''')
                        while export_figure not in possible_answers:
                            print('ERROR: The entered mineral is not among the options.')  
                            export_figure = input('''\nWill you export the figure? < n > or < y >  __  ''')

                        # export the direct figures
                        if export_figure == 'y':
                            export_filename_progenitor = re.sub('(\.\w+)', '', selected_output_file)
                            if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
                                export_name = '{}, {}'.format(export_filename_progenitor, simulation_perspective)
                            else:
                                export_name = export_filename_progenitor

                            export_plot()


            elif individual_plots == 'n':
                pyplot.figure(figsize = (17,10))
                pyplot.title(plot_title, fontsize = 'xx-large')
                pyplot.xlabel('Midpoint module distance (m)', fontsize = 'x-large')
                pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  
                experimental_loop = []
                for mineral in non_zero_minerals:
                    if self.parameters['simulation_type'] == 'transport':
                        iteration = 0
                        distance_serie = []
                        time_serie = []
                        quantity_of_steps_index = 0   
                        for index, row in csv_data.iterrows():
                            if csv_data.at[index, 'time'] == 0:
                                time_serie.append(csv_data.at[index, mineral]) 
                                distance_serie.append(csv_data.at[index, 'dist_x'])
                                quantity_of_steps_index += 1   
                                time = 0

                            elif csv_data.at[index-1, 'soln'] == quantity_of_steps_index:
                                experimental_loop.append('%s [%s] ; time (s): %s' %(mineral,                                                       mineral_formulas[minerals.index(mineral)], round(time, 2)))
                                pyplot.plot(distance_serie,time_serie)
                                distance_serie = []
                                time_serie = []
                                time_serie.append(csv_data.at[index, mineral])
                                distance_serie.append(csv_data.at[index, 'dist_x'])
                                time = csv_data.at[index, 'time']

                            elif index == len(csv_data[mineral]) + 2:   
                                experimental_loop.append('%s [%s] ; time (s): %s' %(mineral, mineral_formulas[minerals.index(mineral)], round(time, 2)))
                                pyplot.plot(distance_serie,time_serie)

                            else:
                                time_serie.append(csv_data.at[index, mineral])
                                distance_serie.append(csv_data.at[index, 'dist_x'])
                                iteration += 1


                    elif self.parameters['simulation_type'] == 'evaporation':
                        pyplot.xlabel('Concentration Factor (CF)', fontsize = 'x-large')
                        pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  

                        experimental_loop = []
                        iteration = 0
                        mass_series = []
                        concentration_series = []
                        quantity_of_steps_index = 0  
                        initial_solution_mass = csv_data.at[0, 'mass_H2O']
                        for index, row in csv_data.iterrows():
                            try:
                                if csv_data.at[index+1, 'mass_H2O']:
                                    if csv_data.at[index, 'step'] >= 1:
                                        concentration_series.append(csv_data.at[index, mineral]) 
                                        solution_mass = csv_data.at[index, 'mass_H2O']
                                        mass_series.append(initial_solution_mass / solution_mass)   

                                    else:
                                        print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')


                            except:
                                concentration_series.append(csv_data.at[index, mineral]) 
                                solution_mass = csv_data.at[index, 'mass_H2O']
                                mass_series.append(initial_solution_mass / solution_mass)  

                                experimental_loop.append('%s [%s]' %(mineral, mineral_formulas[minerals.index(mineral)]))
                                pyplot.plot(mass_series,concentration_series)

                                mass_series = []
                                concentration_series = []
                                concentration_series.append(csv_data.at[index, mineral]) 
                                solution_mass = csv_data.at[index, 'mass_H2O']
                                mass_series.append(initial_solution_mass / solution_mass)         


                pyplot.legend(experimental_loop, loc='best', fontsize = 'x-large')
                pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(simulation_cf), wrap=True, horizontalalignment='left', fontsize=12)
                figure = pyplot.gcf()
                print('\nClose the figure to proceed.')
                pyplot.show()
                export_figure = input('''Will you export the figure? < n > or < y >  __  ''')
                while export_figure not in possible_answers:
                    print('ERROR: The entered mineral is not among the options.')  
                    export_figure = input('''\nWill you export the figure? < n > or < y >  __  ''')

                # export the output graphic
                if export_figure == 'y':
                    export_filename_progenitor = re.sub('(\.\w+)', '', selected_output_file)
                    if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
                        export_name = '{}, {}'.format(export_filename_progenitor, simulation_perspective)
                    else:
                        export_name = export_filename_progenitor

                    export_option = input('''Would you like to export the figure? < y > or < n > __ ''')
                    while export_option not in possible_answers:
                        print('''ERROR: The value is not one of the options.''')  
                        export_option = input('''Would you like to export the figure? < y > or < n > __ ''')       

                    if export_option == 'y':
                        export_plot()


    def export_plot():
        """
        Export the plots to the current working directory  
        """

        if graphical_selection == 'Brine' or (graphical_selection == 'Scaling' and individual_plots == 'n'):
            export_file_name = input('''- What is the name of your export figure?
            Omit < . > and < \ > in the name.
            Default = %s''' %(export_name)) or export_name
            while re.search('(\.|\\))', export_name):
                print('''ERROR: Remove < . > and < \ > from the figure name.''')  
                export_file_name = input('- What will be the name of your export figure?') 
        else:
            export_file_name = input('''- What is the name of your export figure?
            Omit < . > and < \ > in the name.
            Default = %s, %s''' %(export_name, mineral)) or '%s, %s' %(export_name, mineral)
            while re.search('(\.|\\))',export_name):
                print('''ERROR: Remove < . > and < \ > from the figure name.''')  
                export_file_name = input('- What will be the name of your export figure?')

        available_formats = ['jpg', 'png', 'svg']
        export_format = input('''- What will be the format of your export figure?
        Select from < %s >, < %s >, and < %s >.
        Default = < jpg >''' %('jpg', 'png', 'svg')) or 'jpg'
        while export_format not in available_formats:
            print('''ERROR: Select from < %s >, < %s >, and < %s >.''' %('jpg', 'png', 'svg'))  
            export_format = input('- What will be the format of your export figure?')       

        file_number = 0
        if not os.path.exists('%s.%s' %(export_file_name, export_format)):
            figure.savefig('%s.%s' %(export_file_name, export_format))
        elif os.path.exists('%s.%s' %(export_file_name, export_format)):
            while os.path.exists('%s_%s.%s' %(export_file_name, file_number, export_format)):
                file_number += 1
            figure.savefig('%s_%s.%s' %(export_file_name, file_number, export_format))


    def choose_your_path():
        """
        The software user directs the code functions based upon the use case  
        """
        global visualize_selection
        global input_selection
        global execute_selection

        welcome()

        input_selection = input('''- Will you create an input file? < y > or < n > ___ ''')
        while input_selection not in possible_answers:
            print('ERROR: provide < y > or < n > as your answer.')
            input_selection = input('''- Will you create an input file? < y > or < n > ___ ''')

        execute_selection = input('''- Will you execute a simulation?
        < y > or < n >''')
        while execute_selection not in possible_answers:
            print('ERROR: provide < y > or < n > as your answer.')
            execute_selection = input('''- Will you execute a simulation? < y > or < n > ___ ''')

        visualize_selection = input('''- Will you visualize simulation results? < y > or < n > ___ ''')
        while visualize_selection not in possible_answers:
            print('ERROR: provide < y > or < n > as your answer.')
            visualize_selection = input('''- Will you visualize simulation results? < y > or < n > ___ ''')

        while input_selection == execute_selection == visualize_selection == 'n':
            print('ERROR: one of the software options equate < y >.')
            simulation_purpose = input('''- Will you create an input file, execute a simulation, or visualize simulation results?
            < input >, < execute >, or < visualize > __ ''')
            while simulation_purpose not in ['input', 'execute', 'visualize']:
                print('ERROR: provide the purpose of the simulation, or exit the software.')
                simulation_purpose = input('''- Will you create an input file, execute a simulation, or visualize simulation results?
                < input >, < execute >, or < visualize > __ ''')

            if simulation_purpose == 'input':
                input_selection = 'y'

            elif simulation_purpose == 'execute':
                execute_selection = 'y'

            elif simulation_purpose == ' visualize':
                visualize_selection = 'y'

        if visualize_selection == input_selection == execute_selection == 'y':  
            print('\nSimulate and visualize a created input\n', '+'*len('Simulate and visualize a created input'))
            make_general_conditions()
            make_reactive_transport()
            make_solutions()
            make_equilibrium_phases()
            make_selected_output()
            export()
            execute()
            process_selected_output()

        elif execute_selection == visualize_selection == 'y' and input_selection == 'n':
            print('\nSimulation and visualization\n', '+'*len('Simulation and visualization'))
            execute()
            process_selected_output()

        elif input_selection == execute_selection == 'y' and visualize_selection == 'n':
            print('\nSimulate a created input\n', '+'*len('Simulate a created input'))
            make_general_conditions()
            make_reactive_transport()
            make_solutions()
            make_equilibrium_phases()
            make_selected_output()
            export()
            execute()
            process_selected_output()

        elif visualize_selection == 'y' and input_selection == execute_selection == 'n':
            print('\nVisualize\n', '+'*len('Visualize'))
            process_selected_output()

        elif input_selection == 'y' and visualize_selection == execute_selection == 'n':
            print('\nInput creation\n', '+'*len('Input creation'))
            make_general_conditions()
            make_reactive_transport()
            make_solutions()
            make_equilibrium_phases()
            make_selected_output()
            export()

        elif execute_selection == 'y' and input_selection == visualize_selection == 'n':
            print('\nSimulation\n', '+'*len('Simulation'))
            execute()

        final_message = '\n\nThe simulation is complete.'
        print('%s\n%s' %(final_message, '='*len(final_message)))
        
        
        
    def define(self):
        pass
    
    
    def calculate(self):
        if self.parameters['simulation_type'] == 'transport':
            transport()
    
    def execute(self):
        pass
    
    def simulate(self):
        pass