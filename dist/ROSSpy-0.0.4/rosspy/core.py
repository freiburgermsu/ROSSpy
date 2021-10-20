# import libraries
from matplotlib import pyplot 
from to_precision import auto_notation
from scipy.constants import nano, kilo, milli, centi, liter, minute, day, hour
from itertools import chain
from chempy.properties.water_density_tanaka_2001 import water_density
from pubchempy import get_compounds 
from chemicals import periodic_table
import subprocess
import datetime
import pandas
from math import pi, exp, ceil
from glob import glob
import zipfile
import json
import os
import re


def time_determination(time):
    unit = 'seconds'
    if time > 600:
        if time > 7200:
            if time > 2e5:
                time /= day
                unit = 'days'
            else:
                time /= hour
                unit = 'hours'
        else:
            time /= minute
            unit = 'minutes'
    return auto_notation(time, 3), unit


class ROSSPkg():
    def __init__(self, operating_system = 'windows', verbose = False, jupyter = False):      
        # establish the general organization structures
        self.parameters = {}
        self.variables = {}
        self.results = {}
        self.results['figures'] = {}
        self.verbose = verbose
        if operating_system == 'mac':
            import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
        elif operating_system == 'windows':
            import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
        else:
            print(f'--> ERROR: The operating system {operating_system} is not supported.')
        self.phreeqc_mod = phreeqc_mod
        self.parameters['os'] =  operating_system
        
    def define_general(self, database_selection, simulation = 'scaling', domain = 'single', domain_phase = None, quantity_of_modules = 1, simulation_type = 'transport', simulation_title = None):
        '''Establish general conditions'''
        self.parameters['water_mw'] = float(get_compounds('water', 'name')[0].molecular_weight)
        self.parameters['water_grams_per_liter'] = water_density()
        
        # parameterize the input file
        self.parameters['simulation_type'] = simulation_type
        self.parameters['simulation'] = simulation
        self.parameters['quantity_of_modules'] = quantity_of_modules
        self.parameters['domain'] = domain
        self.parameters['domain_phase'] = domain_phase
        self.parameters['root_path'] = os.path.join(os.path.dirname(__file__))
        database_path = os.path.join(self.parameters['root_path'], 'databases', f'{database_selection}.json') 
        
        title_line = f'TITLE\t {simulation_title}'
        if self.parameters['os'] == 'Windows':
            database_line = f'DATABASE {database_path}'
            self.results['general_conditions'] = [database_line, title_line]
        else:
            self.results['general_conditions'] = [title_line]
            
        # establish the database content
        self.parameters['database_selection'] = database_selection 
        database = json.load(open(database_path, 'r'))
        self.elements = database['elements']
        self.minerals = database['minerals']

    def transport(self, simulation_time, simulation_perspective = None, module_characteristics = {}, cells_per_module = 12, parameterized_timestep = None, kinematic_flow_velocity = None, exchange_factor = 1e10):
        '''Define the TRANSPORT block'''
        self.parameters['simulation_time'] = simulation_time
        self.parameters['exchange_factor'] = exchange_factor
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
        
        # exit the function for evaporation simulations
        self.results['transport_block'] = []
        if self.parameters['simulation_type'] == 'evaporation':
            return None

        # define the transport black
        transport_line = '\nTRANSPORT'
        cells_line = '-cells\t\t\t{}'.format(self.parameters['cells_per_module'] * self.parameters['quantity_of_modules'])
        
        self.simulation_shifts = ceil(simulation_time / self.parameters['timestep']) #(self.parameters['cells_per_module']*self.parameters['quantity_of_modules'])
        shifts_line = f'-shifts\t\t\t{self.simulation_shifts}'
        lengths_line = '-lengths\t\t{}'.format(self.variables['cell_meters'])
        timestep_line = '-time_step\t\t{}\t# this satisfies the Courant condition with a feed velocity of {} m/s'.format(self.parameters['timestep'], auto_notation(feed_velocity, 4))
        initial_time_line = '-initial_time\t\t0'    
        boundary_conditions_line = '-boundary_conditions\tconstant\tconstant \t # Dirichlet boundary condition'
        
        # define the domain-dependent parameters
        domain_line = ''
        if self.parameters['domain'] == 'dual':
            domain_line = f'-stagnant\t\t1\t\t{exchange_factor}\t\t\t0.1\t\t0.9 \t # dual domain\n#\t\t\t^stagnant cells\t^exchange factor\t^CP volume\t^bulk volume'
        
        if self.parameters['domain'] == 'single' or (self.parameters['domain'] == 'dual' and self.parameters['domain_phase'] == 'mobile'):
            first_cell = 1
            final_cell = self.parameters['cells_per_module'] * self.parameters['quantity_of_modules']
            if self.parameters['simulation_perspective'] == 'all_distance':
                punch_cells_line = '-punch_cells\t\t{}-{}'.format(first_cell, final_cell)
                punch_frequency_line = f'-punch_frequency\t{self.simulation_shifts}'
            elif self.parameters['simulation_perspective'] == 'all_time':
                punch_cells_line = '-punch_cells\t\t{}'.format(final_cell)
                punch_frequency_line = '-punch_frequency\t1'    
            
        elif self.parameters['domain'] == 'dual' and self.parameters['domain_phase'] == 'immobile':
            first_cell = self.parameters['cells_per_module']* self.parameters['quantity_of_modules']+2
            final_cell = self.parameters['cells_per_module']* self.parameters['quantity_of_modules']*2+1
            if self.parameters['simulation_perspective'] == 'all_distance':
                punch_cells_line = '-punch_cells\t\t{}-{}'.format(first_cell, final_cell)
                punch_frequency_line = f'-punch_frequency\t{self.simulation_shifts}'
            elif self.parameters['simulation_perspective'] == 'all_time':
                punch_cells_line = '-punch_cells\t\t{}'.format(final_cell)
                punch_frequency_line = '-punch_frequency\t1'
        
        # create the transport block
        self.results['transport_block'].extend((transport_line, cells_line, shifts_line, lengths_line, timestep_line, initial_time_line, boundary_conditions_line, domain_line, punch_cells_line, punch_frequency_line))

        # print simulation parameters
        if self.verbose:
            print('\nMembrane thickness (mm):', (self.parameters['repeated_membrane_winding_mm']))
            print('cell length (m): ', self.variables['cell_meters'])
            print('feed velocity (m/s): ', feed_velocity) 
            print('feed_cross_sectional_area (m^2): ', feed_cross_sectional_area)
            print('permeate_removal_per_cell', self.parameters['permeate_moles_per_cell'])
            print('active_cm_squared_cell', (self.parameters['active_m2_cell'] / centi**2))

    def reaction(self, permeate_approach = 'linear_permeate', permeate_efficiency = 1, head_loss = -0.15, final_cf = 2):
        '''Define the REACTION block'''        
        # establish parameters
        self.parameters['permeate_approach'] = permeate_approach
        
        cfs = []
        cell_moles = []
        reaction_parameters = []
        self.results['reaction_block'] = []
        cumulative_cf = 1 
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
                module_iteration = module_previous_moles_removed = moles_removed = 0
                initial_cf = 1

                cf_slope = (final_cf - initial_cf) / self.parameters['cells_per_module']
                for cell in range(self.parameters['cells_per_module']):
                    cell_cf = (cell+1) * cf_slope + initial_cf
                    cfs.append(cell_cf)    

                for cf in cfs:
                    moles_to_be_removed = self.variables['feed_moles'] - (self.variables['feed_moles'] / cf)
                    if module_iteration == 0:
                        moles = self.variables['feed_moles']
                        reaction_parameters.append(moles_to_be_removed)
                    if module_iteration > 0:
                        module_previous_moles_removed += reaction_parameters[-1] 
                        reaction_parameter = moles_to_be_removed - module_previous_moles_removed
                        reaction_parameters.append(reaction_parameter)
                        moles_removed -= moles_to_be_removed

                    module_iteration += 1

                cf = cfs[-1]
                cumulative_cf *= cf
                moles_remaining = self.variables['feed_moles'] - moles_removed  

            if self.parameters['simulation_type'] == 'transport':
                self.results['reaction_block'].append('\n')
                for cell in range(1, self.parameters['cells_per_module']+1):
                    if self.parameters['domain'] == 'single':
                        cell_number = cell + self.parameters['cells_per_module'] * module
                        reaction_index = cell_number-1
                    elif self.parameters['domain'] == 'dual':
                        cell_number = cell + self.parameters['cells_per_module'] * (module + self.parameters['quantity_of_modules']) + 1
                        reaction_index = (cell + self.parameters['cells_per_module'] * module)-1
                        
                    reaction_line = f'REACTION {cell_number}'
                    if cell < self.parameters['cells_per_module']:
                        reaction_line += f'\n\tH2O -1; {reaction_parameters[reaction_index]}' 
                    elif cell == self.parameters['cells_per_module']:
                        reaction_line += f'''\n\tH2O -1; {reaction_parameters[reaction_index]}
        INCREMENTAL_REACTIONS \ttrue'''   
                           
                    self.results['reaction_block'].append(reaction_line)

            elif self.parameters['simulation_type'] == 'evaporation':
                parameter_quantity = 15                          
                recursive_asymptote_multiplier = 1.335449219     # ??? arbitrary assignment of kg of water in the simulation?
                moles_removed = sum(reaction_parameters)
                initial_evaporation_parameter = moles_removed / recursive_asymptote_multiplier
                evaporation_reaction_parameters = ['0', initial_evaporation_parameter]  # ???
                for parameter in range(1, parameter_quantity):
                    evaporation_reaction_parameter = evaporation_reaction_parameters[parameter] * 1/4
                    evaporation_reaction_parameters.append(evaporation_reaction_parameter)

                # define the reaction block
                reaction_line = 'REACTION 1'
                reaction_line += '\n\tH2O -1; '
                reaction_line += ' '.join([str(x) for x in evaporation_reaction_parameters]) 
                reaction_line += ';\nINCREMENTAL_REACTIONS \ttrue'
                self.results['reaction_block'].append(reaction_line)

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
        
        if water_selection == '' and custom_water_parameters == {}:
            print('--> ERROR: Feed geochemistry is not provided.')
        
        self.parameters['water_selection'] = water_selection
        if water_selection == '':
            self.parameters['water_selection'] = solution_description
            
        initial_number = 0
        if self.parameters['simulation_type'] == 'evaporation':
            initial_number = 1
            
        initial_solution_line = '\nSOLUTION {}\t{}'.format(initial_number, self.parameters['water_selection'])
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
                            print('\n--> ERROR: The {} element is not accepted by the {} database'.format(element, self.parameters['database_selection']))
                                
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
            if self.parameters['domain'] == 'dual':
                total_cells = total_cells*2+1
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
            total_cells = self.parameters['cells_per_module'] * self.parameters['quantity_of_modules']
            if self.parameters['domain'] == 'dual':
                total_cells = total_cells*2+1
            equilibrium_phases_number = '1-{}'.format(total_cells)     
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
                    mineral_line = possible_mineral
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
        if output_filename is None:
            count = 0
            selected_output_file_name = '_'.join([str(x) for x in [datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], count]]) 
            while os.path.exists(f'{selected_output_file_name}.txt'):
                count += 1
                selected_output_file_name = '_'.join([str(x) for x in [datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], count]]) 
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

    def export(self, input_path = None, output_path = None):
        """View and export the PHREEQC input file"""        
        # define the simulation input path 
        file_number = 0
        if self.parameters['simulation_type'] == 'evaporation':
            simulation_name = '_'.join([str(x) for x in [datetime.date.today(), 'ROSS', self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], self.parameters['simulation_perspective'], file_number]])
            while os.path.exists(simulation_name):
                file_number += 1
                simulation_name = '_'.join([str(x) for x in [datetime.date.today(), 'ROSS', self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], self.parameters['simulation_perspective'], file_number]])
                
        elif self.parameters['simulation_type'] == 'transport':
            if self.parameters['permeate_approach'] == 'linear_permeate':
                permeate_approach_name = 'LinPerm'
            elif self.parameters['permeate_approach'] == 'linear_cf':
                permeate_approach_name = 'LinCF'

            simulation_name = '_'.join([str(x) for x in [datetime.date.today(), 'ROSS', self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], self.parameters['simulation_perspective'], permeate_approach_name, file_number]])
            while os.path.exists(simulation_name):
                file_number += 1
                simulation_name = '_'.join([str(x) for x in [datetime.date.today(), 'ROSS', self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], self.parameters['simulation_perspective'], permeate_approach_name, file_number]])

        if input_path is None:
            self.parameters['input_file_name'] = 'input.pqi'
            working_directory = os.getcwd()
            self.simulation_path = os.path.join(working_directory, simulation_name)
            os.mkdir(self.simulation_path)
            self.parameters['input_path'] = os.path.join(self.simulation_path, self.parameters['input_file_name'])
        else:
            self.parameters['input_path'] = input_path
            
        # comment the corresponding simulation in the input file
        simulation_line = f'# {self.simulation_path}'
        self.results['solution_block'].insert(0, simulation_line)
        self.results['complete_lines'] = chain(self.results['general_conditions'], self.results['solution_block'], self.results['equilibrium_phases_block'], self.results['reaction_block'], self.results['selected_output_block'], self.results['transport_block']) 
            
        # printing and exporting the input file
        with open(self.parameters['input_path'],'w') as input_file:
            for line in self.results['complete_lines']:
                if self.verbose:
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
        parameters['parameter'].append('simulation_path')
        parameters['value'].append(self.simulation_path)
        for parameter in self.parameters:
            parameters['parameter'].append(parameter)
            parameters['value'].append(self.parameters[parameter])
        parameters_table = pandas.DataFrame(parameters)
        if self.verbose:
            if jupyter:
                display(parameters_table)
            else:
                print(parameters_table)
        
        parameters_path = os.path.join(self.simulation_path, 'parameters.csv')
        parameters_table.to_csv(parameters_path)
        
        # define a table of variables
        variables = {'variable':[], 'value':[]}
        variables['variable'].append('simulation_path')
        variables['value'].append(self.simulation_path)
        for variable in self.variables:
            variables['variable'].append(variable)
            variables['value'].append(self.variables[variable])
        variables_table = pandas.DataFrame(variables)
        if self.verbose:
            if jupyter:
                display(variables_table)
            else:
                print(variables_table)
        
        variables_path = os.path.join(self.simulation_path, 'variables.csv')
        variables_table.to_csv(variables_path)

    def execute(self, simulated_to_real_time = 9.29):
        '''Execute a PHREEQC input file '''
        database_path = os.path.join(self.parameters['root_path'], 'databases','{}.dat'.format(self.parameters['database_selection']))

        def run(input_file, first=False):
            phreeqc = self.phreeqc_mod.IPhreeqc()                 
            phreeqc.load_database(database_path)
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

        def main(input_file):
            import timeit

            def measure_time(func, *args, **kwargs):
                start = timeit.default_timer()
                phreeqc, conc = func(*args, **kwargs)
                return phreeqc, conc, timeit.default_timer() - start

            phreeqc, conc, run_time = measure_time(run, input_file)            
            
            # export the simulation results
            fobj = open(self.parameters['output_path'], 'w')
            headers = conc.keys()
            self.results['csv_data'] = pandas.DataFrame(conc, columns = headers)
            if self.verbose:
                if jupyter:
                    pandas.set_option('display.max_columns', None)
                    display(self.results['csv_data'])
                else:
                    print(self.results['csv_data'])
            fobj.write(self.results['csv_data'].to_string())
            
            self.variables['run_time (s)'] = float(run_time)

        # communicate progress to the user
        estimated_time = ceil(self.parameters['simulation_time'] / simulated_to_real_time)
        if self.parameters['simulation_type'] == 'transport':
            estimated_time = ceil(self.parameters['simulation_time'] / simulated_to_real_time + self.simulation_shifts**0.7)
        
        estimated_completion = datetime.datetime.now() + datetime.timedelta(seconds = estimated_time)
        estimated_time, unit = time_determination(estimated_time)
        print(f'\nEstimated completion in {estimated_time} {unit} by {estimated_completion} local time.')
        
        # execute the simulation
        main(self.input_file)
        if self.verbose:
            print('run_time (s):', self.variables['run_time (s)'])
        
        # verify that the PHREEQC executed and generated the appropriate files
        if not os.path.exists(self.parameters['output_path']):
            print('\nERROR: The simulation failed to execute.')
            
        return self.results['csv_data']

    def process_selected_output(self, selected_output_path = None, plot_title = None, title_font = 'xx-large', label_font = 'x-large', plot_caption = '', legend_title = None, x_label_number = 6, export_name = None, export_format = 'svg', individual_plots = None):
        """Interpreting the PHREEQC SELECTED_OUTPUT file and conducting the plotting functions"""
        databases = [re.search('(\w+)(?=.json)', database).group() for database in glob('./databases/*.json')]
        
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
                self.parameters['database_selection'] = 'pitzer'
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

        self.variables['initial_solution_mass'] = self.results['csv_data']['mass_H2O'][0]
        if self.parameters['simulation_type'] == 'transport':
            self.results['csv_data'].drop(self.results['csv_data'].index[:3], inplace=True)
            self.variables['initial_solution_mass'] = self.results['csv_data']['mass_H2O'][3]
        self.variables['final_solution_mass'] = self.results['csv_data']['mass_H2O'].iloc[-1]
        self.variables['simulation_cf'] = self.variables['initial_solution_mass'] / self.variables['final_solution_mass']
        self.variables['final_time'] = self.results['csv_data']['time'].iloc[-1]
                                 
        # conducting the appropriate visualization function
        if self.parameters['simulation'] == 'brine':
            if plot_title is None:
                plot_title = 'Effluent brine concentrations'
            data = self.brine_plot(pyplot, plot_title, title_font, label_font, plot_caption, export_format, x_label_number)
        elif self.parameters['simulation'] == 'scaling':
            data = self.scaling_plot(plot_title, title_font, label_font, plot_caption, individual_plots, legend_title, x_label_number, export_name, export_format)
        else:
            print('--> ERROR: The < simulation_perspective > parameter is unpredicted.')
        return data
                                 
    def brine_plot(self, pyplot, plot_title, title_font, label_font, plot_caption, export_format, x_label_number = 6):
        """Generate plots of the elemental concentrations from effluent brine in the PHREEQC SELECTED_OUTPUT file"""
        # determine the minerals in the simulation      
        columns = []
        for column in self.results['csv_data'].columns:
            if re.search('([A-Z][a-z]?(?:\(\d\))?(?=\(mol\/kgw\)))', column) and not re.search('(_|H2O|pH)', column):
                columns.append(column)

        # parse the brine concentrations from the raw data
        pyplot.figure(figsize = (17,10))
        non_zero_elements = []
        time = initial_solution_time = 0
        data = {} 
        concentration_serie = []
        
        # plot parameters
        if self.parameters['simulation_perspective'] == 'all_distance':
            x_label = 'Distance (m)'
            if plot_title == '':
                plot_title = 'Brine concentrations along the module'
        elif self.parameters['simulation_perspective'] == 'all_time':
            x_label = 'Time (s)'
            if plot_title == '':
                plotted_time = self.variables['final_time'] - initial_solution_time*self.parameters['timestep']
                plot_title = 'Effluent brine concentration' 
        
        for element in columns:  
            stripped_element = re.search('([A-Z][a-z]?(?:\(\d\))?(?=\(mol\/kgw\)))', element).group()
            non_zero_elements.append(stripped_element)
            concentration_serie = []
            if self.parameters['simulation_perspective'] == 'all_time':
                time_serie = []
                data[element] = {}
                for index, row in self.results['csv_data'].iterrows():
                    if all(row[element] > 1e-16 for element in columns):
                        concentration_serie.append(row[element])
                        time = auto_notation(row['time'], 3)
                        time_serie.append(time) # - initial_solution_time * self.parameters['timestep'])
                        data[element][time] = auto_notation(row[element], 4)
                    else:
                        initial_solution_time += 1
                pyplot.plot(time_serie,concentration_serie)
                                    
            elif self.parameters['simulation_perspective'] == 'all_distance':
                distance_serie = []
                quantity_of_steps_index = 0
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
                        
                        # define the dataframe
                        for x in distance_serie:
                            sigfig_x = auto_notation(x, 3)
                            index = distance_serie.index(x)
                            data[element][sigfig_x] = auto_notation(concentration_serie[index], 4)
                    else:
                        concentration_serie.append(row[element])
                        distance_serie.append(row['dist_x'])                       
                
        # define the brine plot
        if self.parameters['simulation_perspective'] == 'all_time':
            x_location = []
            index_space = len(time_serie)/x_label_number
            for x in range(x_label_number):
                index = int(index_space * x)
                time = time_serie[index]
                x_location.append(time)
            x_location.append(time_serie[-1])
            pyplot.xticks(x_location)
        legend = non_zero_elements
        legend_title = 'non-zero elements'
        pyplot.yscale('log')
        mineral = export_name = None
        self.illustrate(pyplot, legend, legend_title, plot_title, mineral, export_name, label_font, title_font, export_format)

        # defining the datatable of brine concentrations
        concentrations_table = pandas.DataFrame(data)
        concentrations_table.index.name = x_label
        concentrations_table.to_csv(os.path.join(self.simulation_path, 'brine_concentrations.csv'))
        if self.verbose:
            if jupyter:
                display(concentrations_table)
            else:
                print(concentrations_table)
        return concentrations_table            

    def scaling_plot(self, plot_title, title_font, label_font, plot_caption, individual_plots, legend_title = None, x_label_number = 6, export_name = None, export_format = 'svg'):
        """Generate plots of scaling along the module distance in the PHREEQC SELECTED_OUTPUT file  """
        # define reused functions
        def legend_determination(time, mineral, mineral_formula):
            time, unit = time_determination(time)
            if self.parameters['simulation_perspective'] == "all_time":
                if self.parameters['individual_plots']:
                    return time
                elif not self.parameters['individual_plots']:
                    return f'{mineral} [{mineral_formula}] {time} {unit}'
            if self.parameters['simulation_perspective'] == "all_distance":
                if not self.parameters['individual_plots']:
                    if len(set([t for t in self.results['csv_data']['time']])) == 2:
                        return f'{mineral} [{mineral_formula}]'
                    else:
                        return f'{mineral} [{mineral_formula}] {time} {unit}'
                elif self.parameters['individual_plots']:
                    return auto_notation(time, 3)
        
        def all_time_plotting(plot_title):
            def time_serie(mineral, scaling_data):
                legend_entries = []
                scaling_serie = []
                data = {}
                data[f'{mineral} (mol)'] = {}
                time_serie = []
                for index, row in self.results['csv_data'].iterrows():
                    if row[mineral] > 1e-12:
                        scaling_serie.append(row[mineral])
                        time_serie.append(row['time']) # - initial_solution_time * self.parameters['timestep'])
                        data[f'{mineral} (mol)'][auto_notation(row['time'], 4)] = auto_notation(row[mineral], 4)
                        legend_entries.append(legend_determination(row['time'], mineral, mineral_formula))
                        
                # defining the plot
                x_location = []
                index_space = len(time_serie)/x_label_number
                for x in range(x_label_number):
                    index = int(index_space * x)
                    time = time_serie[index]
                    x_location.append(time)
                x_location.extend([time_serie[0], time_serie[-1]])
                x_location = [ceil(float(x)) for x in x_location] 
                pyplot.xticks(x_location, x_location)        
                pyplot.plot(time_serie,scaling_serie)
                
                # assigning the data dictionary
                data_df = pandas.DataFrame(data)
                scaling_data = scaling_data.merge(data_df, how = 'outer', left_index = True, right_index = True)
                return legend_entries, scaling_data
            
            scaling_data = pandas.DataFrame({})
            custom_plot = False
            if plot_title is not None:
                custom_plot = True
            
            for mineral in self.variables['precipitated_minerals']:
                pyplot.figure(figsize = (17,10))
                mineral_formula = self.minerals[mineral]['formula']
                if self.parameters['simulation_type'] == 'transport':
                    legend_entry, scaling_data = time_serie(mineral, scaling_data)
                elif self.parameters['simulation_type'] == 'evaporation':
                    cf_series = []
                    concentration_series = []
                    data_length = len(self.results['csv_data']['mass_H2O'])
                    for index, row in self.results['csv_data'].iterrows():
                        if index < data_length:
                            if row['step'] >= 1:
                                concentration_series.append(row[mineral]) 
                                cf_series.append(self.variables['initial_solution_mass'] / row['mass_H2O'])  
                            elif index > 1:
                                print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')

                    concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                    cf_series.append(initial_solution_mass / row['mass_H2O'])  
                    pyplot.plot(cf_series,concentration_series)                    

                # export the scaling figure
                if not custom_plot:
                    plot_title = f'Module end scaling of {mineral} [{mineral_formula}]'
                legend_title = legend = None
                self.illustrate(pyplot, legend, legend_title, plot_title, mineral, export_name, label_font, title_font, export_format)
            return scaling_data
            
        def all_distance_plotting(plot_title):
            def distance_serie(mineral,):          
                mineral_formula = self.minerals[mineral]['formula'] 
                time = quantity_of_steps_index = 0   
                legend_entries = []
                distance_serie = []
                scaling_serie = []
                data = {}
                data[f'{mineral} (g/m^2)'] = {}
                for index, row in self.results['csv_data'].iterrows():
                    if row['time'] == 0:
                        quantity_of_steps_index += 1
                    elif self.results['csv_data'].at[index-1, 'soln'] == quantity_of_steps_index:
                        if time != 0:
                            legend_entries.append(legend_determination(time, mineral, mineral_formula))

                            pyplot.plot(distance_serie,scaling_serie)
                            grams_area = (float(row[mineral]) * self.minerals[mineral]['mass']) / (self.parameters['active_m2_cell'])
                            scaling_serie.append(auto_notation(grams_area, 3))
                            distance_serie.append(row['dist_x'])
                        time = row['time']

                    elif index == len(self.results['csv_data'][mineral]) + 2:  
                        legend_entries.append(legend_determination(time, mineral, mineral_formula)) 
                        pyplot.plot(distance_serie,scaling_serie)
                        # define the dataframe
                        for x in distance_serie:
                            sigfig_x = auto_notation(x, 3)
                            index = distance_serie.index(x)
                            data[f'{mineral} (g/m^2)'][sigfig_x] = auto_notation(scaling_serie[index], 3)

                    else:
                        grams_area = (float(row[mineral]) * self.minerals[mineral]['mass']) / (self.parameters['active_m2_cell'])
                        scaling_serie.append(grams_area)
                        distance_serie.append(row['dist_x'])

                if self.verbose:
                    print('mineral', mineral,' ', self.minerals[mineral]) 
                return legend_entries, data
        
            scaling_data = pandas.DataFrame({})
            pyplot.figure(figsize = (17,10))
            legend = []
            for mineral in self.variables['precipitated_minerals']: 
                if self.parameters['simulation_type'] == 'transport':
                    legend_entries, data = distance_serie(mineral)
                    data_df = pandas.DataFrame(data)
                    scaling_data = scaling_data.merge(data_df, how = 'outer', left_index = True, right_index = True)
                    legend.extend(legend_entries)
                elif self.parameters['simulation_type'] == 'evaporation':
                    mineral_formula = self.minerals[mineral]['formula'] 
                    legend.append(f'{mineral} [{mineral_formula}]')
                    data = {}
                    data[mineral] = {}
                    cf_series = []
                    scaling_series = []
                    for index, row in self.results['csv_data'].iterrows():
                        if index != len(self.results['csv_data']['mass_H2O']):
                            if row['step'] >= 1:
                                cf = self.variables['initial_solution_mass'] / row['mass_H2O']
                                scale = row[mineral]
                                scaling_series.append(scale) 
                                cf_series.append(cf)   
                                data[mineral][auto_notation(cf, 6)] = scale

                    pyplot.plot(cf_series,scaling_series)
                    data_df = pandas.DataFrame(data)
                    scaling_data = scaling_data.merge(data_df, how = 'outer', left_index = True, right_index = True)
                    if plot_title is None:
                        plot_title = 'Evaporation precipiation from the {}'.format(self.parameters['water_selection'])
                        
            #pyplot.yscale('log')
            if plot_title is None:
                time, units = time_determination(self.variables['final_time'])
                plot_title = 'Scaling from the {} after {} {}'.format(self.parameters['water_selection'],time,units)
            
            return legend, plot_title, scaling_data

        # all of the non-zero minerals are identified and the chemical formulas are sorted into a list
        csv_minerals = []
        non_zero_minerals = set()
        self.variables['precipitated_minerals'] = {}
        for column in self.results['csv_data'].columns:
            if re.search('([A-Z][a-z]{2,})', column) and not re.search('[_(]', column):
                mineral = re.search('([A-Z][a-z]{2,})', column).group()
                csv_minerals.append(mineral)
                for value in self.results['csv_data'][mineral]:
                    if value != 0:
                        non_zero_minerals.add(mineral)   
                        if column in self.minerals:
                            self.variables['precipitated_minerals'][mineral] = self.minerals[mineral]
        
        if non_zero_minerals == set():
            print('No scaling occurred.')
            return None            

        # plot the simulation depending upon the simulation perspective        
        if self.parameters['simulation_perspective'] == "all_time":
            if individual_plots is None:
                self.parameters['individual_plots'] = True
            if self.parameters['individual_plots']:
                scaling_data = all_time_plotting(plot_title)
            elif not self.parameters['individual_plots']:
                legend, plot_title, scaling_data = all_distance_plotting(plot_title)
                mineral = None
                self.illustrate(pyplot, legend, legend_title, plot_title, mineral, export_name, label_font, title_font, export_format)

        elif self.parameters['simulation_perspective'] == 'all_distance':
            if individual_plots is None:
                self.parameters['individual_plots'] = False    
            if not self.parameters['individual_plots']:
                legend, plot_title, scaling_data = all_distance_plotting(plot_title)
                legend_title = 'scale'
                mineral = None
                self.illustrate(pyplot, legend, legend_title, plot_title, mineral, export_name, label_font, title_font, export_format)
            elif self.parameters['individual_plots']:
                scaling_data = all_time_plotting(plot_title)
                                
        # finalize the output data
        x_label, y_label = self.determine_labels()
        scaling_data.index.name = x_label
        scaling_data.sort_index(inplace = True)
        scaling_data.to_csv(os.path.join(self.simulation_path, 'scaling_data.csv'))
        
        return scaling_data
    
    def determine_labels(self): 
        # determine the y-axis label
        if self.parameters['simulation'] == 'scaling':
            y_label = 'Quantity (g/m^2)'
        elif self.parameters['simulation'] == 'brine':
            y_label = 'Concentration (molal)'
        
        # determine the x-axis label
        if self.parameters['simulation_type'] == 'transport':
            if self.parameters['simulation_perspective'] == 'all_distance':
                x_label = 'Distance (m)'
            elif self.parameters['simulation_perspective'] == 'all_time':
                x_label = 'Time (s)'
        elif self.parameters['simulation_type'] == 'evaporation':
            x_label = 'Concentration Factor (CF)'
        return x_label, y_label
    
    def illustrate(self, pyplot, legend, legend_title, plot_title, mineral, export_name, label_font, title_font, export_format):
        # apply the attributes of the figure
        x_label, y_label = self.determine_labels()
        pyplot.grid(True)
        pyplot.title(plot_title, fontsize = title_font)
        pyplot.xlabel(x_label, fontsize = label_font)
        pyplot.ylabel(y_label, fontsize = label_font)  
        if legend is not None:
            pyplot.legend(legend, title = legend_title, loc='best', title_fontsize = 'x-large', fontsize = 'large')
        pyplot.title(plot_title, fontsize = title_font)   
        pyplot.figtext(0.2, 0, 'Final CF: {}'.format(auto_notation(self.variables['simulation_cf'], 4)), wrap=True, horizontalalignment='left', fontsize=12)
        
        figure = pyplot.gcf()
        pyplot.show()
        self.export_plot(figure, plot_title, mineral, export_name, export_format)

    def export_plot(self, figure, plot_title, mineral = None, export_name = None, export_format = 'svg'):
        """Export the plots to the current working directory  """
        # define the output name
        if export_name is None:
            if self.parameters['simulation'] == 'scaling':
                if self.parameters['individual_plots']:
                    export_name = mineral
                else:
                    export_name = 'all_minerals'
            if self.parameters['simulation'] == 'brine':
                export_name = 'brine'
            
        self.results['figures'][export_name] = {'figure':figure, 'title':plot_title}
                                
        # export the plot
        file_number = 0
        figure_path = os.path.join(self.simulation_path, export_name)
        if not os.path.exists('{}.{}'.format(figure_path, export_format)):
            figure.savefig('{}.{}'.format(figure_path, export_format))
        elif os.path.exists('{}.{}'.format(figure_path, export_format)):
            while os.path.exists('{}_{}.{}'.format(figure_path, file_number, export_format)):
                file_number += 1
            figure.savefig('{}_{}.{}'.format(figure_path, file_number, export_format))

    def test(self):
        self.define_general(database_selection = 'pitzer')
        self.transport(simulation_time = 200)
        self.reaction()
        self.solutions(water_selection = 'red_sea')
        self.equilibrium_phases()
        self.selected_output()
        self.export()
        self.execute()
        self.process_selected_output()