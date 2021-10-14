# import libraries
from matplotlib import pyplot 
from to_precision import auto_notation
from scipy.constants import nano, kilo, milli, centi, liter, minute, day, hour
from itertools import chain
from chempy.properties.water_density_tanaka_2001 import water_density
import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
from pubchempy import get_compounds 
from chemicals import periodic_table
import subprocess
import datetime
import pandas
from math import pi, exp, ceil
from glob import glob
import zipfile
import time
import json
import os
import re

# calculation constants
simulated_time_over_computational_time = 9.29    


class ROSSPkg():
    def __init__(self, verbose = False):       
        # establish the general organization structures
        self.parameters = {}
        self.variables = {}
        self.results = {}
        self.results['figures'] = {}
        self.verbose = verbose
        
    def define_general(self, phreeqc_path, database_selection, simulation_title, simulation = 'scaling', operating_system = 'windows', simulation_type = 'transport'):
        '''Establish general conditions'''
        self.parameters['water_mw'] = float(get_compounds('water', 'name')[0].molecular_weight)
        self.parameters['water_grams_per_liter'] = water_density()
        
        # parameterize the input file
        self.parameters['os'] =  operating_system
        self.parameters['phreeqc_path'] = phreeqc_path
        self.parameters['simulation_type'] = simulation_type
        self.parameters['simulation'] = simulation
        self.parameters['root_path'] = os.path.join(os.path.dirname(__file__), '..')
        database_path = os.path.join(self.parameters['root_path'], 'databases', f'{database_selection}.json') 
        print(database_path)
        
        title_line = 'TITLE\t %s' %(simulation_title)
        if operating_system == 'Windows':
            database_line = 'DATABASE %s' %(database_path)
            self.results['general_conditions'] = [database_line, title_line]
        else:
            self.results['general_conditions'] = [title_line]
            
        # establish the database content
        self.parameters['database_selection'] = database_selection 
        database = json.load(open(database_path, 'r'))
        self.elements = database['elements']
        self.minerals = database['minerals']

    def transport(self, simulation_time, module_characteristics = {}, simulation_perspective = None, cells_per_module = 12, domain = 'dual', parameterized_timestep = None, kinematic_flow_velocity = None):
        '''Define the TRANSPORT block'''
        self.parameters['simulation_time'] = simulation_time
        self.parameters['simulation_perspective'] = simulation_perspective
        if self.parameters['simulation_perspective'] is None:
            if self.parameters['simulation'] == 'scaling':
                self.parameters['simulation_perspective'] = 'all_distance'
            elif self.parameters['simulation'] == 'brine':
                self.parameters['simulation_perspective'] = 'all_time'
        
        # assign default RO module dimensions 
        self.parameters['module_diameter_mm'] =  201                 
        self.parameters['permeate_tube_diameter_mm'] =  29            
        self.parameters['module_length_m'] =  1.016                 
        self.parameters['permeate_flow_m3_per_day'] = 40                  
        self.parameters['max_feed_flow_m3_per_hour'] = 15.9                  
        self.parameters['membrane_thickness_mm'] = 250 * (nano / milli)   
        self.parameters['feed_thickness_mm'] = 0.8636
        self.parameters['active_m2'] = 37
        self.parameters['permeate_thickness_mm'] = 0.3               
        self.parameters['polysulfonic_layer_thickness_mm'] = 0.05     
        self.parameters['support_layer_thickness_mm'] = 0.15           

        # assign parameterized module dimensions
        for parameter in module_characteristics:
            self.parameters[parameter] = module_characteristics[parameter]

        # calculate module properties
        self.parameters['repeated_membrane_winding_mm'] = 2 * self.parameters['membrane_thickness_mm'] + self.parameters['feed_thickness_mm'] + self.parameters['permeate_thickness_mm'] + 2 * self.parameters['polysulfonic_layer_thickness_mm'] + 2 * self.parameters['support_layer_thickness_mm']     
        self.parameters['cells_per_module'] = cells_per_module
        self.variables['cell_meters'] = self.parameters['module_length_m'] / self.parameters['cells_per_module']        
        self.parameters['active_m2_cell'] = self.parameters['active_m2'] / self.parameters['cells_per_module']
        
        module_cross_sectional_area = self.parameters['module_diameter_mm']**2 * pi / 4        #squared millimeters
        permeate_tube_cross_sectional_area = self.parameters['permeate_tube_diameter_mm']**2 * pi / 4     #squared millimeters
        filtration_cross_sectional_area = (module_cross_sectional_area - permeate_tube_cross_sectional_area) * milli**2         #squared meters
        feed_cross_sectional_area = (self.parameters['feed_thickness_mm'] / self.parameters['repeated_membrane_winding_mm']) * filtration_cross_sectional_area       #squared meters
        self.variables['feed_cubic_meters'] = feed_cross_sectional_area * self.parameters['module_length_m'] 
        self.variables['feed_kg'] = self.variables['feed_cubic_meters'] / liter * self.parameters['water_grams_per_liter'] * milli    
        self.variables['feed_moles'] = self.variables['feed_kg'] * kilo / self.parameters['water_mw'] 

        # calculate fluid flow characteristics
        if not kinematic_flow_velocity:
            kinematic_flow_velocity = 9.33E-7    #square meters / second
        feed_velocity = self.parameters['max_feed_flow_m3_per_hour'] / (feed_cross_sectional_area) / hour     #meters / second
        reynolds_number = feed_velocity * (self.parameters['module_diameter_mm'] - self.parameters['permeate_tube_diameter_mm']) * milli**2 / kinematic_flow_velocity
        self.variables['Reynold\'s number'] = reynolds_number

        # calculate module cell characteristics
        self.variables['feed_kg_cell'] = self.variables['feed_kg'] / self.parameters['cells_per_module']   
        self.variables['feed_moles_cell'] = self.variables['feed_moles'] / self.parameters['cells_per_module']   

        # calculate simulation timestep that adheres to the Courant condition   
        self.parameters['timestep'] = self.parameters['module_length_m'] / feed_velocity  # seconds
        if parameterized_timestep:
            self.parameters['timestep'] = parameterized_timestep
            
        courant_timestep = self.variables['cell_meters'] / feed_velocity
        if self.parameters['timestep'] < courant_timestep:
            self.parameters['timestep'] = courant_timestep
            
        self.parameters['permeate_moles_per_cell'] = (self.parameters['permeate_flow_m3_per_day'] / day / liter * self.parameters['water_grams_per_liter'] / self.parameters['water_mw']) * (self.parameters['timestep'] / self.parameters['cells_per_module'])      #moles / (cell * self.parameters['timestep'])

        # define the transport black
        transport_line = '\nTRANSPORT'
        cells_line = '-cells\t\t\t{}'.format(self.parameters['cells_per_module'])
        
        self.simulation_shifts = ceil(simulation_time / self.parameters['timestep']) #(self.parameters['cells_per_module']*self.parameters['quantity_of_modules'])
        shifts_line = f'-shifts\t\t\t{self.simulation_shifts}'
        lengths_line = '-lengths\t\t{}'.format(self.variables['cell_meters'])
        timestep_line = '-time_step\t\t{}\t# this satisfies the Courant condition with a feed velocity of {} m/s'.format(self.parameters['timestep'], auto_notation(feed_velocity, 4))
        initial_time_line = '-initial_time\t\t0'    
        boundary_conditions_line = '-boundary_conditions\tconstant\tconstant \t # Dirichlet boundary condition'
        
        if domain == 'single':
            domain_line = '-stagnant\t\t0\t\t0\t\t\t0\t\t0 \t # single domain\n#\t\t\t^stagnant cells\t^exchange factor\t^CP volume\t^bulk volume'
        elif domain == 'dual':
            domain_line = '-stagnant\t\t1\t\t1\t\t\t0.1\t\t0.9 \t # dual domain\n#\t\t\t^stagnant cells\t^exchange factor\t^CP volume\t^bulk volume'
        
        if self.parameters['simulation_perspective'] == 'all_distance':
            punch_cells_line = '-punch_cells\t\t1-{}'.format(self.parameters['cells_per_module'])
            punch_frequency_line = f'-punch_frequency\t{self.simulation_shifts}'
        elif self.parameters['simulation_perspective'] == 'all_time':
            punch_cells_line = '-punch_cells\t\t{}'.format(self.parameters['cells_per_module'])
            punch_frequency_line = '-punch_frequency\t1'     
            

        # create the transport block
        self.results['transport_block'] = []
        self.results['transport_block'].extend((transport_line, cells_line, shifts_line, lengths_line, timestep_line, initial_time_line, boundary_conditions_line, domain_line, punch_cells_line, punch_frequency_line))

        # print simulation parameters
        if self.verbose:
            print('\nMembrane thickness (mm):', (self.parameters['repeated_membrane_winding_mm']))
            print('cell length (m): ', self.variables['cell_meters'])
            print('feed velocity (m/s): ', feed_velocity) 
            print('feed_cross_sectional_area (m^2): ', feed_cross_sectional_area)
            print('permeate_removal_per_cell', self.parameters['permeate_moles_per_cell'])
            print('active_cm_squared_cell', (self.parameters['active_m2_cell'] / centi**2))

    def reaction(self, quantity_of_modules = 1, permeate_approach = 'linear_permeate', permeate_efficiency = 1, head_loss = -0.15, final_cf = 2):
        '''Define the REACTION block'''
        # establish parameters
        self.parameters['permeate_approach'] = permeate_approach
        self.parameters['quantity_of_modules'] = quantity_of_modules
        
        cfs = []
        cell_moles = []
        reaction_parameters = []
        iteration = 0
        cumulative_cf = 1 
        self.results['reaction_block'] = []
        for module in range(self.parameters['quantity_of_modules']): 
            if permeate_approach == 'linear_permeate':
                initial_moles_removed = self.parameters['permeate_moles_per_cell'] * 2 / (1 + exp(head_loss))
                final_moles_removed = initial_moles_removed * exp(head_loss)
                removed_moles_slope = ((final_moles_removed - initial_moles_removed) / (self.parameters['cells_per_module'])) / permeate_efficiency
                average_moles_removed = (final_moles_removed + initial_moles_removed) / 2

                for cell in range(self.parameters['cells_per_module']):
                    removed_moles_in_cell = (cell * removed_moles_slope + initial_moles_removed)
                    reaction_parameters.append(removed_moles_in_cell)
                    
                moles_removed = sum(reaction_parameters)
                cf = self.variables['feed_moles'] / (self.variables['feed_moles'] - moles_removed)
                self.variables[f'module {module}'] = {'cumulative_cf': cf, 'permeate (moles/cell)': removed_moles_slope}
                
                moles_remaining = self.variables['feed_moles'] - moles_removed
                iteration += 1
                
                if auto_notation(average_moles_removed,14) != auto_notation(self.parameters['permeate_moles_per_cell'], 14):
                    print('--> ERROR: Inconsistent REACTION calculations.', )
                    print('average_moles_removed', average_moles_removed)
                    print('permeate_moles_per_cell', self.parameters['permeate_moles_per_cell'])
               
                if self.verbose:
                    print('\n')
                    print('permeate (moles/cell) ', removed_moles_slope)
                    print('moles_remaining', moles_remaining)
                    print('moles_removed', moles_removed)

            if permeate_approach == 'linear_cf':
                module_iteration = 0
                initial_cf = 1

                cf_slope = (final_cf - initial_cf) / self.parameters['cells_per_module']
                for cell in range(self.parameters['cells_per_module']):
                    cell_cf = (cell+1) * cf_slope + initial_cf
                    cfs.append(cell_cf)    

                for cf in cfs:
                    moles_to_be_removed =  self.parameters['feed_moles_cell'] - (self.parameters['feed_moles_cell'] / cf)
                    if module_iteration == 0:
                        moles = self.variables['feed_moles']
                        reaction_parameters.append(moles_to_be_removed)
                    if module_iteration > 0:
                        module_previous_moles_removed += reaction_parameters[-1] 
                        reaction_parameter = moles_to_be_removed - module_previous_moles_removed
                        reaction_parameters.append(reaction_parameter)
                        moles -= moles_to_be_removed

                    module_iteration += 1

                cf = cfs[-1]
                cumulative_cf *= cf
                moles = self.parameters['feed_moles'] - moles_to_be_removed   # moles_to_be_removed = the final quantity of moles is calculated with naming from the linear permeate flux method to be consistent  

            if self.parameters['simulation_type'] == 'transport':
                self.results['reaction_block'].append('\n')
                for cell in range(1, self.parameters['cells_per_module']+1):
                    cell_number = (cell) + self.parameters['cells_per_module'] * module
                    reaction_line = f'REACTION {cell_number}'
                    
                    reaction_index = cell_number-1
                    if (cell) < self.parameters['cells_per_module']:
                        reaction_line += f'\n\tH2O -1; {reaction_parameters[reaction_index]}' 
                    elif (cell) == self.parameters['cells_per_module']:
                        reaction_line += f'''\n\tH2O -1; {reaction_parameters[reaction_index]}
        INCREMENTAL_REACTIONS \ttrue'''   
                           
                    self.results['reaction_block'].append(reaction_line)

            elif self.parameters['simulation_type'] == 'evaporation':
                parameter_quantity = 15                          
                recursive_assymtote_multiplier = 1.335449219     # ??? arbitrary assignment of kg of water in the simulation?
                moles_removed = sum(reaction_parameters)
                initial_evaporation_parameter = moles_removed / recursive_assymptote_multiplier
                evaporation_reaction_parameters = ['0', initial_evaporation_parameter]  # ???
                for parameter in range(1, parameter_quantity):
                    evaporation_reaction_parameter = evaporation_reaction_parameters[parameter] * 1/4
                    evaporation_reaction_parameters.append(evaporation_reaction_parameter)

                # define the reaction block
                reaction_line = 'REACTION 1'
                reaction_line += '\n\tH2O -1; '
                reaction_line += ' '.join(evaporation_reaction_parameters) 
                self.parameters['reaction_block'] = [reaction_line]

            # the calculated reaction parameters will be added and printed to a generated PHREEQC input file
            final_solution_mass = moles_remaining * self.parameters['water_mw'] * milli  #kg water mass
            final_cf_cell = self.variables['feed_kg'] / final_solution_mass

            if self.parameters['os'] == 'windows':
                self.results['reaction_block'].append('#%s' %(permeate_approach))
                if permeate_approach == 'linear permeate':
                    self.results['reaction_block'].append(f'''
        #Permeate efficiency parameter: {permeate_efficiency}
        #Head loss parameter: {head_loss}''')

                self.results['reaction_block'].append(f'''    #Effluent module {module + 1}:
        #Estimated CF: {auto_notation(cf, 4)}
        #Estimated solution mass: {final_solution_mass}\n\n''')

            if self.verbose:
                print('Effluent module %s CF:' %(module + 1), final_cf_cell)

    def solutions(self, water_selection = '', custom_water_parameters = {}, solution_description = '', parameterized_alkalinity = False, parameterized_ph_charge = True):
        """Specify the SOLUTION block of the simulation."""
        # create the solution line of the input file
        self.results['solution_block'] = []
        
        if water_selection != '':
            self.parameters['water_selection'] = water_selection
        if self.parameters['simulation_type'] == 'transport':
            initial_solution_line = '\nSOLUTION 0\t%s' % (solution_description)
        elif self.parameters['simulation_type'] == 'evaporation':
            initial_solution_line = '\nSOLUTION 1\t%s' % (solution_description)
        self.results['solution_block'].append(initial_solution_line)

        #=============================================================================
        # determine which predefined solution should be simulated

        elements_lines = []
        self.parameters['solution_elements'] = []
        temperature = ph = alkalinity = pe = None
        if water_selection != '':       
            # import the predefined water body
            water_file_path = os.path.join(self.parameters['root_path'], 'water_bodies', f'{water_selection}.json')
            water_body = json.load(open(water_file_path))
            
            for content, information in water_body.items():
                if content == 'element':
                    for element, information2 in information.items():
                        self.parameters['solution_elements'].append(element)
                        if element in self.elements:
                            self.parameters['solution_elements'].append(element)
                            conc = information2['concentration (ppm)']
                            ref = information2['reference']
                            elements_lines.append(f'{element}\t{conc}\t#{ref}')
                        else:
                            print('\n--> ERROR: The {} element is not accepted by the {} database')
                                
                elif content == 'temperature':
                    temperature = information['celcius']
                    temperature_reference = information['reference']
                elif content == 'pe':
                    pe = information['value']
                    pe_reference = information['reference']
                elif content == 'Alkalinity':
                    alkalinity = information['value']
                    alkalinity_reference = information['reference'] 
                elif content == 'pH':
                    ph = information['value']
                    ph_reference = information['reference']

        if custom_water_parameters != {}:
            for content, information in custom_water_parameters.items():
                if content == 'element':
                    for element, information2 in information.items():
                        self.parameters['solution_elements'].append(element)
                        conc = information2['concentration (ppm)']
                        ref = information2['reference']
                        if element in self.elements:                          
                            if len(str(conc)) <= 3:
                                elements_lines.append(f'{element}\t\t{conc}\t#{ref}')
                            else:
                                elements_lines.append(f'{element}\t{conc}\t#{ref}')
                        else:
                            print('\n--> ERROR: The {} element is not accepted by the {} database'.format(element, self.parameters['database_selection']))
                                    
                # create the temperature line of the input file
                elif content == 'temperature':                    
                    temperature = custom_water_parameters['temperature']['value']
                    temperature_reference = custom_water_parameters['temperature']['reference']
                elif content == 'pe':       
                    pe = custom_water_parameters['pe']['value']
                    pe_reference = custom_water_parameters['pe']['reference']
                elif content == 'Alkalinity':
                    alkalinity = custom_water_parameters['Alkalinity']['value']
                    alkalinity_reference = custom_water_parameters['Alkalinity']['reference'] 
                elif content == 'pH':
                    ph = custom_water_parameters['ph']['value']
                    ph_reference = custom_water_parameters['ph']['reference']
                    
        # parameterize the lines of the SOLUTIONS block
        temperature_line = ''
        if temperature is not None:
            temperature_line = f'temp \t {temperature} \t #{temperature_reference}.'
        pe_line = ''
        if pe is not None:
            pe_line = f'pe \t\t {pe} \t   #{pe_reference} // 4.00 is the default (?)'     

        alkalinity_line = ''
        ph_line = ''
        if ph is not None:
            if parameterized_ph_charge and not parameterized_alkalinity:
                ph_line = f'pH \t\t {ph} charge #{ph_reference}'
                alkalinity_line = ''      
            elif parameterized_alkalinity and not parameterized_ph_charge:
                if alkalinity is not None:
                    ph_line = f'pH \t\t {ph} #{ph_reference}'
                    alkalinity_line = f'Alkalinity \t {alkalinity} #{alkalinity_reference}'
                   
        unit_line = 'units \t ppm' 
        elements_line = '\n'.join(elements_lines)
        if water_selection == 'Bakken formation':
            water_line = '-water \t{}\t#TDS=300 ppthousand [before fudging]'.format(self.variables['feed_kg'])
        elif water_selection == 'German Basin':
            water_line = '-water \t{}\t#TDS=314 ppthousand [before fudging]'.format(self.variables['feed_kg'])
        else:
            water_line = '-water \t{}'.format(self.variables['feed_kg'])

        self.results['solution_block'].extend([temperature_line, ph_line, pe_line, alkalinity_line, unit_line, elements_line, water_line])

        #parameterize the initial module solution
        if self.parameters['simulation_type'] == 'transport':
            total_cells = self.parameters['cells_per_module'] * self.parameters['quantity_of_modules']
            feed_solution_line = f'\nSOLUTION 1-{total_cells}\tInitial solution in the RO module'
            self.results['solution_block'].extend([feed_solution_line,'temp \t 25','units \t ppm'])

            for element in self.parameters['solution_elements']:
                element_line = f'{element}\t0'    
                self.results['solution_block'].append(element_line)

            water_line = '-water \t{}'.format(self.variables['feed_kg'])
            self.results['solution_block'].append(water_line)


    def equilibrium_phases(self, block_comment = '', ignored_minerals = [], existing_parameters = {}):
        """Specify the EQUILIBRIUM_PHASES block of the simulation."""
        # define mineral sizes for later spacing
        short_mineral_names = ['Barite', 'Gypsum','Halite', 'Natron', 'Quartz', 'Talc','Trona', 'Borax', 'Albite', 'K-mica','Illite', 'Pyrite', 'Sulfur',]
        long_mineral_names = ['Anthophyllite', 'Hexahydrite', 'Leonhardite', 'Nesquehonite', 'Pentahydrite', 'Portlandite','Sepiolite(d)', 'Boric_acid,s', 'K2B4O7:4H2O', 'NaB5O8:5H2O', 'Rhodochrosite', 'Strontianite','Hydroxyapatite', 'Chlorite(14A)', 'Mackinawite', 'Hausmannite', 'Pyrochroite']
        
        # determine the set of possible minerals 
        self.variables['described_minerals'] = {}
        for mineral in self.minerals:
            mineral_formula = self.minerals[mineral]['formula']
            original_formula = mineral_formula
            
            # remove entities is an ordered fashion
            mineral_formula = re.sub('(H2O|OH|CO3)', '', mineral_formula)
            if 'S(6)' in self.parameters['solution_elements']:
                mineral_formula = re.sub('(SO4)', '', mineral_formula)
            mineral_formula = re.sub('([0-9()•:.])', '', mineral_formula)
            mineral_elements = re.findall('[A-Z][a-z]?', mineral_formula)
            
            for element in mineral_elements:
                if all(element in self.parameters['solution_elements'] for element in mineral_elements):
                    self.variables['described_minerals'][mineral] = self.minerals[mineral]
                elif element not in self.parameters['solution_elements']:
                    if self.verbose:
                        print('--> The {} element is absent in the solution to describe the mineral < {} / {} >.'.format(element, original_formula, mineral))
                        
        # define the equilibrium_phases block
        self.results['equilibrium_phases_block'] = []
        if self.parameters['simulation_type'] == 'transport':
            equilibrium_phases_number = '1-{}'.format(self.simulation_shifts)
        elif self.parameters['simulation_type'] == 'evaporation':
            equilibrium_phases_number = '1'
            
        equilibrium_phases_line = f'\nEQUILIBRIUM_PHASES {equilibrium_phases_number}\t{block_comment}'
        self.results['equilibrium_phases_block'].append(equilibrium_phases_line)

        # define the equilibrium_phases lines for the code block
        for possible_mineral in self.variables['described_minerals']:
            if possible_mineral not in ignored_minerals:
                if possible_mineral in short_mineral_names:
                    mineral_line = f'{possible_mineral}\t\t' 
                elif possible_mineral == 'Ca-Montmorillonite':
                    mineral_line = f'{possible_mineral}'
                else:
                    mineral_line = f'{possible_mineral}\t'

                if possible_mineral in existing_parameters:
                    for key, value in existing_parameters[possible_mineral].items():
                        if key == 'saturation':
                            mineral_saturation = value['saturation']
                            mineral_line += f'\t{mineral_saturation}'
                        if key == 'initial_moles':
                            initial_moles = value['initial_moles']
                            mineral_line += f'\t{initial_moles}'
                else:
                    mineral_line += f'\t0\t0'
                    
                self.results['equilibrium_phases_block'].append(mineral_line)     
                
        if self.verbose:
            print(self.variables['described_minerals'])
            
    def selected_output(self, output_filename = None):
        '''Specify the output file after a PHREEQC simulation'''
        # create parameter lines 
        water_selection = ''
        if 'water_selection' in self.parameters:
            water_selection = self.parameters['water_selection']
            
        if output_filename is None:
            count = 0
            selected_output_file_name = '_'.join([str(x) for x in [datetime.date.today(), water_selection, self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], count]]) 
            while os.path.exists(f'{selected_output_file_name}.txt'):
                count += 1
                selected_output_file_name = '_'.join([str(x) for x in [datetime.date.today(), water_selection, self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], count]]) 
        else:
            selected_output_file_name = output_filename

        selected_output_file_name += '.txt'

        minerals_line = ' '.join([mineral for mineral in self.variables['described_minerals']])
        elements_line = ' '.join([element for element in self.parameters['solution_elements']])

        # define parameter lines
        first_line = '\nSELECTED_OUTPUT'
        file_name_line = f'-file\t\t\t{selected_output_file_name}'
        reaction_line = '-reaction\t\ttrue'
        temperature_line = '-temperature\t\ttrue'
        total_elements_line = '-totals\t\t\t' + elements_line    
        saturation_indices_line = f'-saturation_indices\t{minerals_line}'
        equilibrium_phases_line = f'-equilibrium_phases\t{minerals_line}'
        ph_line = '-pH\t\t\ttrue'
        solution_line = '-solution'
        time_line = '-time\t\t\ttrue'
        distance_line = '-distance\t\ttrue'
        simulation_line = '-simulation\t\ttrue'
        high_precision_line = '-high_precision\ttrue'
        step_line = '-step'
        water_line = '-water'

        # establish the selected_output_block
        self.results['selected_output_block'] = []
        self.results['selected_output_block'].extend((first_line, file_name_line, reaction_line, temperature_line, total_elements_line, saturation_indices_line, equilibrium_phases_line, ph_line, time_line, distance_line, simulation_line, high_precision_line, solution_line, step_line,water_line))

    def export(self, input_path = None, output_path = None, print_block = True):
        """View and export the PHREEQC input file"""
        # define the simulation input path 
        water_selection = ''
        if 'water_selection' in self.parameters:
            water_selection = self.parameters['water_selection']
        if input_path is None:
            file_number = 0
            if self.parameters['permeate_approach'] == 'linear_permeate':
                permeate_approach_name = 'LinPerm'
            elif self.parameters['permeate_approach'] == 'linear_cf':
                permeate_approach_name = 'LinCF'

            simulation_name = '_'.join([str(x) for x in [datetime.date.today(), 'ROSS', water_selection, self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], self.parameters['simulation_perspective'], permeate_approach_name, file_number]])
            while os.path.exists(simulation_name):
                file_number += 1
                simulation_name = '_'.join([str(x) for x in [datetime.date.today(), 'ROSS', water_selection, self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], self.parameters['simulation_perspective'], permeate_approach_name, file_number]])

            self.parameters['input_file_name'] = 'input.pqi'
            working_directory = os.getcwd()
            self.simulation_path = os.path.join(working_directory, simulation_name)
            os.mkdir(self.simulation_path)
            self.parameters['input_path'] = os.path.join(self.simulation_path, self.parameters['input_file_name'])
        else:
            self.parameters['input_path'] = input_path
            
        # comment the corresponding simulation in the input file
        simulation_line = '# {}'.format(self.simulation_path)
        self.results['solution_block'].insert(0, simulation_line)
        self.results['complete_lines'] = chain(self.results['general_conditions'], self.results['solution_block'], self.results['equilibrium_phases_block'], self.results['reaction_block'], self.results['selected_output_block'], self.results['transport_block']) 
            
        # printing and exporting the input file
        with open(self.parameters['input_path'],'w') as input_file:
            for line in self.results['complete_lines']:
                if print_block:
                    print(line)
                input_file.write(line + '\n')
            
        self.input_file = open(self.parameters['input_path'],'r').read()
        
        # define the simulation output path 
        if output_path is None:
            self.parameters['output_file_name'] = 'selected_output.pqo'
            self.parameters['output_path'] = os.path.join(self.simulation_path, self.parameters['output_file_name'])
        else:
            self.parameters['output_path'] = output_path    
            
        
        # define a table of parameters
        parameters = {'parameter':[], 'value':[]}
        parameters['parameter'].append('simulation')
        parameters['value'].append(self.simulation_path)
        for parameter in self.parameters:
            parameters['parameter'].append(parameter)
            parameters['value'].append(self.parameters[parameter])
        parameters_table = pandas.DataFrame(parameters)
        if print_block:
            display(parameters_table)
        
        parameters_file_name = 'parameters.csv'
        parameters_path = os.path.join(self.simulation_path, parameters_file_name)
        parameters_table.to_csv(parameters_path)
        
        # define a table of variables
        variables = {'variable':[], 'value':[]}
        variables['variable'].append('simulation_path')
        variables['value'].append(self.simulation_path)
        for variable in self.variables:
            variables['variable'].append(variable)
            variables['value'].append(self.variables[variable])
        variables_table = pandas.DataFrame(variables)
        if print_block:
            display(variables_table)
        
        variables_file_name = 'variables.csv'
        variables_path = os.path.join(self.simulation_path, variables_file_name)
        variables_table.to_csv(variables_path)

    def execute(self, print_output = True, simulated_to_real_time = 7):
        '''Execute a PHREEQC input file '''
        database_path = os.path.join(self.parameters['phreeqc_path'], 'database\\{}.dat'.format(self.parameters['database_selection']))

        def run(input_file, first=False):
            phreeqc = phreeqc_mod.IPhreeqc()                 
            phreeqc.load_database(r"%s" %(database_path))
            phreeqc.run_string(input_file)
            
            # define the conc dictionary
            output = phreeqc.get_selected_output_array()

            header = output[0]
            conc = {}
            for head in header:
                conc[head] = []
            for row in output[1:]:
                for col, head in enumerate(header):
                    conc[head].append(row[col])
                    
            return phreeqc, conc

        def main(input_file, print_output):
            import timeit

            def measure_time(func, *args, **kwargs):
                start = timeit.default_timer()
                phreeqc, conc = func(*args, **kwargs)
                return phreeqc, conc, timeit.default_timer() - start

            phreeqc, conc, run_time = measure_time(run, input_file, print_output)
            
            
            # export the simulation results
            fobj = open(self.parameters['output_path'], 'w')
            headers = conc.keys()
            self.results['csv_data'] = pandas.DataFrame(conc, columns = headers)
            if print_output:
                pandas.set_option('display.max_columns', None)
                display(self.results['csv_data'])
            fobj.write(self.results['csv_data'].to_string())
            
            self.variables['run_time (s)'] = run_time

        # communicate progress to the user
        estimated_time = ceil(self.parameters['simulation_time'] / simulated_to_real_time)
        estimated_completion = datetime.datetime.now() + datetime.timedelta(seconds = estimated_time)
        unit = 'seconds'
        if estimated_time > 600:
            if estimated_time > 7200:
                if estimated_time > 2e5:
                    estimated_time /= day
                    unit = 'days'
                else:
                    estimated_time /= hour
                    unit = 'hours'
            else:
                estimated_time /= minute
                unit = 'minutes'
                        
        print(f'\nEstimated completion in {estimated_time} {unit} by {estimated_completion} local time.')
        
        # execute the simulation
        main(self.input_file, print_output)
        if self.verbose:
            print('run_time (s):', self.variables['run_time (s)'])
        
        # verify that the PHREEQC executed and generated the appropriate files
        if not os.path.exists(self.parameters['output_path']):
            print('\nERROR: The simulation failed to execute.')
            
        return self.results['csv_data']

    def process_selected_output(self, selected_output_path = None, plot_title = '', title_font = 'xx-large', label_font = 'x-large', plot_caption = '', table_title = None, export_format = 'svg', individual_plots = True):
        """Interpreting the PHREEQC SELECTED_OUTPUT file and conducting the plotting functions"""
        
        databases = [database for database in glob('./databases/*.json')]
        databases = [re.search('(\w+)(?=.json)', database).group() for database in databases]

        if 'csv_data' not in self.results:
            # determining the appropriate variables
            if selected_output_path is not None:
                # define the simulation perspective
                self.parameters['simulation'] = 'scaling'
                if re.search('(brine)', selected_output_path, flags=re.IGNORECASE):
                    self.parameters['simulation'] = 'brine'

                # define the simulation type
                self.parameters['simulation_type'] = 'transport'
                for line in selected_output_path:
                    if re.search('(evaporation)', line, re.IGNORECASE):
                        self.parameters['simulation_type'] = 'evaporation'

                # define the database contents
                self.parameters['database_selection'] = 'llnl'
                for database in databases:
                    names = database.split('_')
                    if all(re.search(name, selected_output_path, re.IGNORECASE) for name in names):
                        self.parameters['database_selection'] = database                   
                        break                    

                database_json = json.load(open('./databases/{}.json'.format(self.parameters['database_selection'])))
                self.minerals = database_json['minerals']
                self.elements = database_json['elements']

                # define the simulation
                self.parameters['selected_output_file_name'] = re.search('(\w+)(?=\.)', selected_output_path).group()

            else:
                working_directory = os.getcwd()
                selected_output_path = os.path.join(working_directory, self.parameters['selected_output_file_name'])                     

            # preparing the SELECTED_OUTPUT file into a dataframe
            selected_output = open(selected_output_path, 'r')
            original_data = pandas.read_table(selected_output, sep = '\t')
            self.results['csv_data'] = pandas.DataFrame(original_data)
            for column in self.results['csv_data'].columns:
                new_column = column.strip()
                self.results['csv_data'].rename(columns={column:new_column}, inplace = True)

        self.variables['initial_solution_mass'] = self.results['csv_data'].at[0, 'mass_H2O']
        final_solution_mass = self.results['csv_data']['mass_H2O'].iloc[-1]
        self.variables['simulation_cf'] = self.variables['initial_solution_mass'] / final_solution_mass
                                 
        # conducting the appropriate visualization function
        if self.parameters['simulation'] == 'brine':
            data = self.brine_plot(pyplot, plot_title, title_font, label_font, plot_caption, table_title, export_format)
        elif self.parameters['simulation'] == 'scaling':
            data = self.scaling_plot(plot_title, title_font, label_font, plot_caption, table_title, individual_plots, export_format)
        else:
            print('--> ERROR: The < simulation_perspective > parameter is unpredicted.')
        return data
                                 
    def brine_plot(self, pyplot, plot_title, title_font, label_font, plot_caption, table_title, export_format, x_label_number = 6):
        """Generate plots of the elemental concentrations from effluent brine in the PHREEQC SELECTED_OUTPUT file  """
        # determine the minerals in the simulation      
        columns = []
        for column in self.results['csv_data'].columns:
            if re.search('([A-Z][a-z]?(?:\(\d\))?(?=\(mol\/kgw\)))', column) and not re.search('(_|H2O|pH)', column):
                columns.append(column)
        self.results['csv_data'].drop(self.results['csv_data'].index[:3], inplace=True)

        # parse the brine concentrations from the raw data
        pyplot.figure(figsize = (17,10))
        pyplot.figtext(0.2, 0, 'Final CF: {}'.format(auto_notation(self.variables['simulation_cf'], 4)), wrap=True, horizontalalignment='left', fontsize=12)  
        non_zero_elements = []
        loop_iteration = 1
        time = initial_solution_time = 0
        data = {} 
        concentration_serie = []
        total_time = self.results['csv_data']['time'].iloc[-1]
        
        for element in columns:  
            stripped_element = re.search('([A-Z][a-z]?(?:\(\d\))?(?=\(mol\/kgw\)))', element).group()
            non_zero_elements.append(stripped_element)
            concentration_serie = []
            if self.parameters['simulation_perspective'] == 'all_time':
                time_serie = []
                data[element] = {}
                for index, row in self.results['csv_data'].iterrows():
                    if all(row[element] > 1e-12 for element in columns):
                        concentration_serie.append(row[element])
                        time_serie.append(auto_notation(row['time'], 4)) # - initial_solution_time * self.parameters['timestep'])
                        data[element][row['time']] = row[element]
                    else:
                        initial_solution_time += 1
                pyplot.plot(time_serie,concentration_serie)

                # define the brine figure
                plotted_time = total_time - initial_solution_time*self.parameters['timestep']
                if table_title is None:
                    table_title = f'Molal concentrations of feed elements over {auto_notation(plotted_time, 3)} seconds' 
                    
                x_location = []
                x_label_distance = float(time_serie[-1])/x_label_number
                for x in range(x_label_number):
                    index = ceil(x * x_label_distance)
                    time = time_serie[index]
                    x_location.append(time)
                x_location.append(time_serie[-1])
                pyplot.xticks(x_location)
                x_label = 'Time (s)'
                if plot_title == '':
                    plot_title = 'Effluent brine concentrations'
                
            elif self.parameters['simulation_perspective'] == 'all_distance':
                distance_serie = []
                quantity_of_steps_index = average_iteration = 0
                data[element] = {}
                for index, row in self.results['csv_data'].iterrows():
                    if row['time'] == 0:            
                        quantity_of_steps_index += 1                    
                    elif self.results['csv_data'].at[index-1,'soln'] == quantity_of_steps_index:   
                        pyplot.plot(distance_serie,concentration_serie)
                        concentration_serie.append(row[element])
                        distance_serie.append(row['dist_x'])
                    elif index == len(self.results['csv_data'][element]) + 2:
                        pyplot.plot(distance_serie,concentration_serie)   
                        average_iteration += 1
                        loop_iteration += 1
                        
                        # define the dataframe
                        for x in distance_serie:
                            sigfig_x = auto_notation(x, 3)
                            index = distance_serie.index(x)
                            data[element][f'{sigfig_x} m'] = auto_notation(concentration_serie[index], 3)
                    else:
                        concentration_serie.append(row[element])
                        distance_serie.append(row['dist_x'])                       
                        
                # plot parameters
                x_label = 'Distance (m)'
                if plot_title == '':
                    plot_title = 'Brine concentrations along the module'
                
        # define the brine plot
        pyplot.grid(True)
        pyplot.title(plot_title, fontsize = title_font)   
        pyplot.xlabel(x_label, fontsize = label_font)
        pyplot.ylabel('Concentration (molal)', fontsize = label_font)
        pyplot.yscale('log')
        pyplot.legend(non_zero_elements, loc='best', title = 'non-zero elements', fontsize = 'x-large')
        figure = pyplot.gcf()
        pyplot.show()
        self.export_plot(figure, plot_title, export_format = export_format)

        # defining the datatable of brine concentrations
        concentrations_table = pandas.DataFrame(data)
        concentrations_table.index.name = x_label
        if table_title is None:
            table_title = 'Average elemental molal concentrations of the feed water in the RO module for each %s seconds of simulation:' %(quantity_of_steps_index)
        concentrations_table.to_csv(os.path.join(self.simulation_path, 'brine_concentrations.csv'))
        
        if self.verbose:
            print('\n\n\n',table_title,'\n%s'%('='*len(table_title)))
            print(concentrations_table)
            
        return concentrations_table
                                 

    def scaling_plot(self, plot_title, title_font, label_font, plot_caption, table_title, individual_plots, export_format = 'svg'):
        """Generate plots of scaling along the module distance in the PHREEQC SELECTED_OUTPUT file  """
        # define reused functions
        def series_creation():
            time = iteration = quantity_of_steps_index = 0   
            legend_entry = []
            distance_serie = []
            scaling_serie = []
            for index, row in self.results['csv_data'].iterrows():
                if row['time'] == 0:
                    quantity_of_steps_index += 1
                elif self.results['csv_data'].at[index-1, 'soln'] == quantity_of_steps_index:
                    if time != 0:
                        legend_entry.append(f'{mineral} [{mineral_formula}] ; {auto_notation(time, 3)} sec')
                        pyplot.plot(distance_serie,scaling_serie)
                        
                        distance_serie = [] 
                        scaling_serie = []
                        grams_area = (float(row[mineral]) * self.minerals[mineral]['mass']) / (self.parameters['active_m2_cell'])
                        scaling_serie.append(auto_notation(grams_area, 3))
                        distance_serie.append(row['dist_x'])
                    time = row['time']

                elif index == len(self.results['csv_data'][mineral]) + 2:   
                    legend_entry.append(f'{mineral} [{mineral_formula}] ; {auto_notation(time, 3)} sec')
                    pyplot.plot(distance_serie,scaling_serie)
                    
                    # define the dataframe
                    for x in distance_serie:
                        sigfig_x = auto_notation(x, 3)
                        index = distance_serie.index(x)
                        data[mineral] = {}
                        data[mineral][f'{sigfig_x} m'] = auto_notation(scaling_serie[index], 3)
                        
                else:
                    grams_area = (float(row[mineral]) * self.minerals[mineral]['mass']) / (self.parameters['active_m2_cell'])
                    scaling_serie.append(grams_area)
                    distance_serie.append(row['dist_x'])
                    iteration += 1
                    
            if self.verbose:
                print(mineral, self.minerals[mineral])
                                 
            return legend_entry, data
                                 
        def illustrate(pyplot, legend_entries):                     
            pyplot.legend(legend_entry, loc='best', fontsize = 'x-large')
            pyplot.figtext(0.2, 0, 'Final CF: {}'.format(auto_notation(self.variables['simulation_cf'], 4)), wrap=True, horizontalalignment='left', fontsize=12)
            figure = pyplot.gcf()
            pyplot.show()
                        
            return figure
                                 
        # the complete list of all minerals is created
        if self.parameters['simulation_type'] == 'transport':
            self.results['csv_data'].drop(self.results['csv_data'].index[:3], inplace=True)
            
        csv_minerals = []
        for column in self.results['csv_data'].columns:
            if re.search('([A-Z].{3,})', column) and not re.search('[(_:]', column):
                csv_minerals.append(column)

        # all of the non-zero minerals are identified and the chemical formulas are sorted into a list
        non_zero_minerals = set()
        for mineral in csv_minerals:
            for value in self.results['csv_data'][mineral]:
                if value != 0:
                    non_zero_minerals.add(mineral)   
        if non_zero_minerals == set():
            print('No scaling occurred.')
            return None
        
        # define a dictionary of precipitated minerals
        self.variables['precipitated_minerals'] = {}
        for mineral in self.minerals:
            if mineral in non_zero_minerals:
                self.variables['precipitated_minerals'][mineral] = self.minerals[mineral]

        # plot the simulation depending upon the simulation perspective
        unit = 'moles'
        data = {}
        scaling_data = pandas.DataFrame(data)
        if self.parameters['simulation_perspective'] == "all_time":
            legend_entry = []
            formula_index = 0
            for mineral in non_zero_minerals:
                mineral_formula = self.minerals[mineral]['formula']  
                mineral_serie = scaling_serie = []
                for index, row in self.results['csv_data'].iterrows():
                    mineral_serie.append(row[mineral]) 
                    time = row['time']
                    scaling_serie.append(time)

                pyplot.plot(scaling_serie,mineral_serie)
                pyplot.scatter(scaling_serie,mineral_serie)
                legend_entry.append(f'{mineral} [{mineral_formula}]')   

                # export the figure
                figure = illustrate(pyplot, legend_entry)
                self.results['figures'][self.parameters['selected_output_file_name']] = {'figure':figure, 'title':plot_title}
                self.export_plot(figure, plot_title, export_format)

        elif self.parameters['simulation_perspective'] == 'all_distance':
            legend_entries = []
            if individual_plots:
                for mineral in non_zero_minerals:
                    print(mineral)
                    mineral_formula = self.minerals[mineral]['formula']
                    pyplot.figure(figsize = (17,10))
                    pyplot.title(plot_title, fontsize = title_font)

                    if self.parameters['simulation_type'] == 'transport':
                        pyplot.xlabel('Midpoint module distance (m)', fontsize = label_font)
                        pyplot.ylabel('Quantity (g/m^2)', fontsize = label_font) 
                        legend_entry, data = series_creation()
                        data_df = pandas.DataFrame(data)
                        scaling_data = scaling_data.append(data_df)

                    elif self.parameters['simulation_type'] == 'evaporation':
                        pyplot.xlabel('Concentration Factor (CF)', fontsize = label_font)
                        pyplot.ylabel('Quantity (%s)' %(unit), fontsize = label_font)  

                        legend_entry = []
                        cf_series = []
                        concentration_series = []
                        data_length = len(self.results['csv_data']['mass_H2O'])
                        for index, row in self.results['csv_data'].iterrows():
                            if index < data_length:
                                if self.results['csv_data'].at[index, 'step'] >= 1:
                                    concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                                    solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                                    cf_series.append(self.variables['initial_solution_mass'] / solution_mass)  
                                elif index > 1:
                                    print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')

                        concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                        solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                        cf_series.append(self.variables['initial_solution_mass'] / solution_mass)  

                        legend_entry.append(f'{mineral} [{mineral_formula}]')
                        pyplot.plot(cf_series,concentration_series)                    
                    
                    legend_entries.append(legend_entry)
                    
                    # export the figure
                    figure = illustrate(pyplot, legend_entries)
                    scaling_data.index.name = 'scale (g/m^2)'
                    scaling_data.to_csv('scaling_data.csv')
                    self.export_plot(figure, plot_title, mineral, export_format)

            elif not individual_plots:
                legend_entry = []
                for mineral in non_zero_minerals:
                    if self.parameters['simulation_type'] == 'transport':
                        legend_entry = series_creation()

                    elif self.parameters['simulation_type'] == 'evaporation':
                        pyplot.xlabel('Concentration Factor (CF)', fontsize = label_font)
                        pyplot.ylabel('Quantity (%s)' %(unit), fontsize = label_font)  

                        legend_entry = []
                        iteration = 0
                        cf_series = []
                        concentration_series = []
                        self.variables['initial_solution_mass'] = self.results['csv_data'].at[0, 'mass_H2O']
                        for index, row in self.results['csv_data'].iterrows():
                            try:
                                if self.results['csv_data'].at[index+1, 'mass_H2O']:
                                    if self.results['csv_data'].at[index, 'step'] >= 1:
                                        concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                                        solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                                        cf_series.append(self.variables['initial_solution_mass'] / solution_mass)   
                                    else:
                                        print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')

                            except:
                                concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                                solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                                cf_series.append(self.variables['initial_solution_mass'] / solution_mass)  

                                legend_entry.append(f'{mineral} [{mineral_formula}]')
                                pyplot.plot(mass_series,concentration_series)

                                cf_series = []
                                concentration_series = []
                                concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                                solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                                cf_series.append(self.variables['initial_solution_mass'] / solution_mass)         
                                                         
                # export the figure
                figure = illustrate(pyplot, legend_entries)
                if export_figure:
                    self.export_plot(figure, plot_title, export_format = export_format)
                       
        return scaling_data

    def export_plot(self, figure, plot_title, mineral = None, individual_plots = False, export_name = None, export_format = 'svg'):
        """Export the plots to the current working directory  """
        # define the output name
        if export_name is None:
            if self.parameters['simulation'] == 'scaling':
                export_name = mineral
            if self.parameters['simulation'] == 'brine':
                export_name = 'brine'
            
        self.results['figures'][export_name] = {'figure':figure, 'title':plot_title}
                                
        # export the plot
        file_number = 0
        figure_path = os.path.join(self.simulation_path, export_name)
        if not os.path.exists('{}.{}'.format(figure_path, export_format)):
            self.results['figures'][export_name]['figure'].savefig('{}.{}'.format(figure_path, export_format))
        elif os.path.exists('{}.{}'.format(figure_path, export_format)):
            while os.path.exists('{}_{}.{}'.format(figure_path, file_number, export_format)):
                file_number += 1
            figure.savefig('{}_{}.{}'.format(figure_path, file_number, export_format))


    def input_file(self, operating_system, phreeqc_path, database_selection, simulation_type, simulation_title, water_selection, quantity_of_modules = 1, module_characteristics = {}, simulation = 'scaling', domain = 'dual', permeate_approach = 'linear permeate', permeate_efficiency = 1, head_loss = -0.15, final_cf = 2, custom_water_parameters = {}, ignored_minerals = [], existing_parameters = {}, export_figure = True):
        """Concisely create an input file of the software """
        self.define_general(operating_system, phreeqc_path, database_selection, simulation_type, simulation_title)
        self.transport(module_characteristics = module_characteristics, quantity_of_modules = quantity_of_modules, domain = domain, simulation = simulation)
        self.reaction(permeate_approach = permeate_approach, permeate_efficiency = permeate_efficiency, head_loss = head_loss, final_cf = final_cf)
        self.solutions(water_selection = water_selection, custom_water_parameters = custom_water_parameters)
        self.equilibrium_phases(ignored_minerals = ignored_minerals, existing_parameters = existing_parameters)
        self.selected_output()
        self.export()
        
            
    def complete_simulation():
        self.execute()
        self.process_selected_output(export_figure = export_figure)