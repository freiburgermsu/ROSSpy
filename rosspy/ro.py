from scipy.constants import kilo, milli, centi, liter, minute, day, hour
from chempy.properties.water_density_tanaka_2001 import water_density
from numpy import array, log10, seterr
from math import pi, ceil, inf
from collections import OrderedDict
from matplotlib import pyplot 
from itertools import chain
from warnings import warn
from pprint import pprint
from glob import glob
import datetime
import sigfig 
import pandas
import chemw
import timeit, json, os, re

seterr(divide = 'ignore') 

def sigfigs_conversion(num, sigfigs_in = 2):
    return sigfig.round(num, sigfigs=sigfigs_in, notation = 'sci', warn = False)

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
    return int(float(sigfigs_conversion(time, 3))), unit

def isnumber(num):
    try:
        float(num)
        return True
    except:
        try:
            int(num)
            return True
        except:
            return False

class ROSSPkg():
    def __init__(self, database_selection, simulation = 'scaling', simulation_type = 'transport', operating_system = 'windows',  export_content = True, domain_phase = None, quantity_of_modules = 1, simulation_title = None, verbose = False, printing = True, jupyter = False):      
        # establish the general organization structures and values
        self.parameters, self.variables, self.results, self.results['figures'] = {}, {}, {}, {}
        self.printing, self.verbose, self.jupyter, self.export_content = printing, verbose, jupyter, export_content
        self.figure_title = None
        
        # define calculation values
        self.chem_mw = chemw.ChemMW(printing = False,)
        self.water_mw = float(self.chem_mw.mass('H2O'))
        self.water_gL = water_density()
        
        # import PHREEQpy
        self.parameters['os'] =  operating_system
        if operating_system == 'windows':
            import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
        elif operating_system == 'unix':
            import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
        else:
            self._error(f'The operating system {operating_system} is not supported.', 'type')
        self.phreeqc_mod = phreeqc_mod
        
        # define paths
        self.parameters['root_path'] = os.path.join(os.path.dirname(__file__))
        self.parameters['database_path'] = os.path.join(self.parameters['root_path'], 'databases',f'{database_selection}.dat')
        
        # define options
        self.databases = [re.search(r'(?<=processed_databases.)(.+)(?=.json)', database).group() for database in glob(
                os.path.join(self.parameters['root_path'], 'processed_databases', '*.json')
                )]
        self.feed_sources = [os.path.basename(feed).split('.')[0] for feed in glob(
                os.path.join(self.parameters['root_path'], 'water_bodies', '*.json')
                )]
        with open(os.path.join(self.parameters['root_path'], 'ro_module.json')) as module:
            self.ro_modules = json.load(module)
        
        # define simulation parameters
        self.parameters['database_selection'] = database_selection 
        self.parameters['quantity_of_modules'] = quantity_of_modules
        self.parameters['simulation_type'] = simulation_type
        self.parameters['simulation'] = simulation
        self._define_database()
        
        self.parameters['domain'] = 'single'
        if domain_phase in ['mobile', 'immobile']:
            self.parameters['domain_phase'] = domain_phase
            self.parameters['domain'] = 'dual'
        elif domain_phase is not None:
            self._error(f'The {domain_phase} domain phase must be either < mobile > or < immobile >.', 'value')

        # create the general_conditions block of the input file
        title_line = f'TITLE\t {simulation_title}'
        database_line = 'DATABASE {}'.format(self.parameters['database_path'])
        self.results['general_conditions'] = [database_line, title_line]
            
        
    def _define_database(self,):
        database_json = os.path.join(self.parameters['root_path'], 'processed_databases', self.parameters['database_selection']+'.json') 
        with open(database_json, 'r') as db:
            database = json.load(db)
        self.elements = database['elements']
        self.minerals = database['minerals']

    def _transport(self, simulation_time, simulation_perspective = None, module_characteristics = {}, ro_module = 'BW30-400', cells_per_module = 12, coarse_timestep = True, kinematic_flow_velocity = None, exchange_factor = 1e5):
        '''Define the TRANSPORT block'''
        self.parameters['simulation_time'], self.parameters['exchange_factor'] = simulation_time, exchange_factor
        self.parameters['coarse_timestep'], self.parameters['simulation_perspective'] = coarse_timestep, simulation_perspective
        if self.parameters['simulation_perspective'] is None:
            if self.parameters['simulation'] == 'scaling':
                self.parameters['simulation_perspective'] = 'all_distance'
            elif self.parameters['simulation'] == 'brine':
                self.parameters['simulation_perspective'] = 'all_time'
        
        # assign default RO module dimensions 
        if ro_module not in self.ro_modules:
            self._error(f'The parameterized ro_module {ro_module} is not defined in the ro_module.json parameter file ({self.ro_modules.keys()}))', 'value')
        self.ro_module = self.ro_modules[ro_module]
        for parameter in module_characteristics:
            self.ro_module[parameter] = module_characteristics[parameter]

        # calculate module properties
        self.parameters['repeated_membrane_winding_mm'] = float(
                2*self.ro_module['membrane_thickness_mm']['value'] + self.ro_module['feed_thickness_mm']['value'] + 
                self.ro_module['permeate_thickness_mm']['value'] + 2*self.ro_module['polysulfonic_layer_thickness_mm']['value'] +
                2*self.ro_module['support_layer_thickness_mm']['value']   
                )  
        if self.parameters['repeated_membrane_winding_mm'] == 0:
            warn('--> ERROR: The module dimensions are not sensible.')
        self.parameters['cells_per_module'] = cells_per_module
        self.variables['cell_meters'] = self.ro_module['module_length_m']['value'] / self.parameters['cells_per_module']        
        self.parameters['active_m2_cell'] = self.ro_module['active_m2']['value'] / self.parameters['cells_per_module']
        
        module_cross_sectional_area = self.ro_module['module_diameter_mm']['value']**2 * pi/4        #squared millimeters
        permeate_tube_cross_sectional_area = self.ro_module['permeate_tube_diameter_mm']['value']**2 * pi/4     #squared millimeters
        filtration_cross_sectional_area = (module_cross_sectional_area - permeate_tube_cross_sectional_area) * milli**2         #squared meters
        feed_cross_sectional_area = (
                self.ro_module['feed_thickness_mm']['value']/self.parameters['repeated_membrane_winding_mm']
                ) * filtration_cross_sectional_area       #squared meters
        self.variables['feed_cubic_meters'] = feed_cross_sectional_area * self.ro_module['module_length_m'] ['value']
        self.variables['feed_kg'] = self.variables['feed_cubic_meters']/liter * self.water_gL*milli    
        self.variables['feed_moles'] = self.variables['feed_kg']*kilo / self.water_mw 

        # calculate fluid flow characteristics
        if not kinematic_flow_velocity:
            kinematic_flow_velocity = 9.33E-7    #square meters / second
        feed_velocity = self.ro_module['max_feed_flow_m3_per_hour']['value']/hour / feed_cross_sectional_area     #meters / second
        self.variables['Reynold\'s number'] = feed_velocity * (
                self.ro_module['module_diameter_mm']['value'] - self.ro_module['permeate_tube_diameter_mm']['value']
                )*milli**2 / kinematic_flow_velocity

        # calculate module cell characteristics
        self.variables['feed_kg_cell'] = self.variables['feed_kg'] / self.parameters['cells_per_module']   
        self.variables['feed_moles_cell'] = self.variables['feed_moles'] / self.parameters['cells_per_module']  
        self.parameters['timestep'] = self.ro_module['module_length_m']['value']/feed_velocity
        if not self.parameters['coarse_timestep']:
            self.parameters['timestep'] = self.variables['cell_meters']/feed_velocity
        self.parameters['permeate_moles_per_cell'] = self.ro_module['permeate_flow_m3_per_hour']['value']/hour/liter * self.water_gL/self.water_mw * self.variables['cell_meters']/feed_velocity     #moles/cell
        
        # exit the function for evaporation simulations
        self.results['transport_block'] = []
        if self.parameters['simulation_type'] == 'evaporation':
            return None

        # define the transport black
        self.parameters['simulated_cells'] = self.parameters['cells_per_module']*self.parameters['quantity_of_modules']
        transport_line = '\nTRANSPORT'
        cells_line = f'-cells\t\t\t{self.parameters["simulated_cells"]}'
        
        self.simulation_shifts = ceil(simulation_time / self.parameters['timestep']) 
        shifts_line = f'-shifts\t\t\t{self.simulation_shifts}'
        lengths_line = '-lengths\t\t{}'.format(self.variables['cell_meters'])
        timestep_line = '-time_step\t\t{}\t# this satisfies the Courant condition with a feed velocity of {} m/s'.format(self.parameters['timestep'], sigfigs_conversion(feed_velocity, 4))
        initial_time_line = '-initial_time\t\t0'    
        boundary_conditions_line = '-boundary_conditions\tconstant\tflux \t # Dirichlet and Cachy boundary conditions, respectively'
        
        # define the domain-dependent parameters
        domain_line = ''
        if self.parameters['domain'] == 'dual':
            domain_line = f'''-stagnant\t\t 1 \t\t{exchange_factor}\t\t\t{0.05}\t\t{0.95} \t # dual domain
                            #\t\t\t^stagnant cells\t^exchange factor\t^CP porosity\t^bulk porosity''' # porosity here means volume fraction
        
        first_cell, final_cell = 1, self.parameters['simulated_cells']
        if self.parameters['domain'] == 'single' or (self.parameters['domain'] == 'dual' and self.parameters['domain_phase'] == 'mobile'):
            if self.parameters['simulation_perspective'] == 'all_distance':
                punch_cells_line = '-punch_cells\t\t{}-{}'.format(first_cell, final_cell)
                punch_frequency_line = f'-punch_frequency\t{self.simulation_shifts}'
            elif self.parameters['simulation_perspective'] == 'all_time':
                punch_cells_line = '-punch_cells\t\t{}'.format(final_cell)
                punch_frequency_line = '-punch_frequency\t1'    
            
        elif self.parameters['domain'] == 'dual' and self.parameters['domain_phase'] == 'immobile':
            first_cell, final_cell = self.parameters['simulated_cells']+2, self.parameters['simulated_cells']*2+1
            if self.parameters['simulation_perspective'] == 'all_distance':
                punch_cells_line = '-punch_cells\t\t{}-{}'.format(first_cell, final_cell)
                punch_frequency_line = f'-punch_frequency\t{self.simulation_shifts}'
            elif self.parameters['simulation_perspective'] == 'all_time':
                punch_cells_line = '-punch_cells\t\t{}'.format(final_cell)
                punch_frequency_line = '-punch_frequency\t1'
        
        # create the transport block
        self.results['transport_block'].extend([transport_line, cells_line, shifts_line, lengths_line, timestep_line, initial_time_line, boundary_conditions_line, domain_line, punch_cells_line, punch_frequency_line])

        # print simulation parameters
        if self.verbose:
            print('\nMembrane thickness (mm):', (self.parameters['repeated_membrane_winding_mm']))
            print('cell length (m): ', self.variables['cell_meters'])
            print('feed velocity (m/s): ', feed_velocity) 
            print('feed_cross_sectional_area (m^2): ', feed_cross_sectional_area)
            print('permeate_removal_per_cell', self.parameters['permeate_moles_per_cell'])
            print('active_cm_squared_cell', (self.parameters['active_m2_cell']/centi**2))

    def _reaction(self, final_cf = None, permeate_efficiency = 1, head_loss = 0.1, evaporation_steps = 15):     
        '''Define the REACTION blocks'''
        self.parameters['permeate_approach'] = 'linear_permeate'
        if isnumber(final_cf):
            self.parameters['permeate_approach'] = 'linear_cf'
        elif final_cf is not None:
            self._error(f'The final_cf parameter value < {final_cf} > is not a number.', 'value')
    
        # determine the permeate fluxes
        reaction_parameters, self.results['reaction_block'] = [], []
        self.cumulative_cf = 1 
        module_previous_moles_removed = 0 
        for module in range(self.parameters['quantity_of_modules']): 
            if self.parameters['permeate_approach'] == 'linear_permeate':
                initial_moles_removed = self.parameters['permeate_moles_per_cell']*2 / (1+(1-head_loss))
                final_moles_removed = initial_moles_removed * (1-head_loss)
                removed_moles_slope = ((final_moles_removed - initial_moles_removed)/self.parameters['cells_per_module']) / permeate_efficiency
                average_moles_removed = (final_moles_removed + initial_moles_removed) / 2

                for cell in range(self.parameters['cells_per_module']):
                    removed_moles_in_cell = (cell * removed_moles_slope + initial_moles_removed)
                    reaction_parameters.append(removed_moles_in_cell)
                moles_removed = sum(reaction_parameters)
                cf = self.cumulative_cf = self.variables['feed_moles'] / (self.variables['feed_moles'] - moles_removed)
                
                if sigfigs_conversion(average_moles_removed,14) != sigfigs_conversion(self.parameters['permeate_moles_per_cell'], 14):
                    warn(f'--> ERROR: Inconsistent REACTION calculations. The two definitions of the average permeate flux per cell < {average_moles_removed} > & < {self.parameters["permeate_moles_per_cell"]} > do not equate.')

            elif self.parameters['permeate_approach'] == 'linear_cf':
                module_iteration = moles_removed = 0
                initial_cf = 1
                cfs = []
                cf_slope = (final_cf - initial_cf) / self.parameters['simulated_cells']
                for cell in range(self.parameters['cells_per_module']*(module), self.parameters['cells_per_module']*(module+1)):
                    cell_cf = (cell+1) * cf_slope + initial_cf
                    cfs.append(cell_cf)  
                for cf in cfs:
                    if self.verbose:
                        print('\n',)
                    moles_to_be_removed = self.variables['feed_moles'] - (self.variables['feed_moles'] / cf)
                    if module_iteration == 0:
                        if len(reaction_parameters) > 0:
                            if self.verbose:
                                print('to be', moles_to_be_removed, 
                                      '\tlast parameter', reaction_parameters[-1], 
                                      '\tprevious removal', sum(reaction_parameters)
                                      )
                            module_previous_moles_removed += reaction_parameters[-1]
                        else:
                            if self.verbose:
                                print('to be', moles_to_be_removed, )
                        if module == 0:
                            reaction_parameters.append(moles_to_be_removed)
                        else:
                            reaction_parameters.append(moles_to_be_removed - module_previous_moles_removed)
                    else:
                        module_previous_moles_removed += reaction_parameters[-1] 
                        if moles_to_be_removed < module_previous_moles_removed:
                            warn(f'--> ERROR: The reaction parameter is negative: {moles_to_be_removed - module_previous_moles_removed}.')
                        reaction_parameter = moles_to_be_removed - module_previous_moles_removed
                        if self.verbose:
                            print('to be', moles_to_be_removed, '\tlast parameter', reaction_parameters[-1], '\tthis parameter', reaction_parameter, '\tprevious removal', module_previous_moles_removed)
                        reaction_parameters.append(reaction_parameter)

                    moles_removed = sum(reaction_parameters)
                    module_iteration += 1

                    # verify the CF calculations and efffects in the system
                    measured_cf = self.variables['feed_moles'] / (self.variables['feed_moles']-moles_removed)
                    consistency = sigfig.round(measured_cf, 5) == sigfig.round(cf, 5)
                    if not consistency:
                        warn(f'--> ERROR: The measured cf {measured_cf} is dissimilar from the target cf {cf}.') 
                    else:
                        if self.verbose:
                            print('cf consistency between the measured and predicted: ', consistency)

                removed_moles_slope = ((reaction_parameters[-1] - reaction_parameters[0])/self.parameters['cells_per_module']) / permeate_efficiency
                self.cumulative_cf = cf = cfs[-1]
                
            moles_remaining = self.variables['feed_moles'] - moles_removed  
            self.variables[f'module {module}'] = {'cf': cf, 'permeate (moles/cell)': removed_moles_slope}

            if self.parameters['simulation_type'] == 'transport':
                for cell in range(1, self.parameters['cells_per_module']+1):
                    if self.parameters['domain'] == 'single':
                        cell_number = cell + self.parameters['cells_per_module']*module
                        reaction_index = cell_number-1
                    elif self.parameters['domain'] == 'dual':
                        cell_number = cell + 1 + self.parameters['cells_per_module'] * (module + self.parameters['quantity_of_modules'])
                        reaction_index = (cell + self.parameters['cells_per_module']*module)-1
                        
                    reaction_line = f'\nREACTION {cell_number}'
                    if cell < self.parameters['cells_per_module']:
                        reaction_line += f'\n\tH2O -1; {reaction_parameters[reaction_index]}' 
                    elif cell == self.parameters['cells_per_module']:
                        reaction_line += f'''\n\tH2O -1; {reaction_parameters[reaction_index]}'''   
                           
                    self.results['reaction_block'].append(reaction_line)
                    
                # the calculated reaction parameters will be added and printed to a generated PHREEQC input file
                final_solution_mass = float(moles_remaining*self.water_mw*milli)  #kg water mass
                if self.parameters['os'] == 'windows':
                    self.results['reaction_block'].append('# {}'.format(self.parameters['permeate_approach']))
                    if self.parameters['permeate_approach'] == 'linear_permeate':
                        self.results['reaction_block'].append(f'''
    #Permeate efficiency parameter: {permeate_efficiency*100}%
    #Effluent relative pressure: {(1-head_loss)*100}%''')

                    self.results['reaction_block'].append(f'''    
    #Effluent module {module + 1}:
#Estimated CF: {sigfigs_conversion(cf, 4)}
#Estimated final water mass: {final_solution_mass}\n\n''')

            elif self.parameters['simulation_type'] == 'evaporation':
                # define the reaction block
                self.results['reaction_block'].append(f'''
                REACTION 1
                    H2O -1; {moles_removed} in {evaporation_steps} step;
                    INCREMENTAL_REACTIONS   true
                ''')

            if self.printing:
                print(f'Effluent module {module + 1} CF: {cf}')
            if self.verbose:
                print('\n')
                print('permeate (moles/cell) ', removed_moles_slope)
                print('moles_remaining', moles_remaining)
                print('moles_removed', moles_removed)

                
    def reactive_transport(self, simulation_time, simulation_perspective = None, final_cf = None, module_characteristics = {}, ro_module = 'BW30-400', permeate_efficiency = 1, head_loss = 0.1, evaporation_steps = 15, cells_per_module = 12, coarse_timestep = True, kinematic_flow_velocity = None, exchange_factor = 1e5):
        self._transport(simulation_time, simulation_perspective, module_characteristics, ro_module, cells_per_module, coarse_timestep, kinematic_flow_velocity, exchange_factor)
        self._reaction(final_cf, permeate_efficiency, head_loss, evaporation_steps)

                       
    def _described_minerals(self,):
        # determine the set of possible minerals 
        self.variables['described_minerals'] = {}
        for mineral in self.minerals:
            mineral_formula = self.minerals[mineral]['formula']
            original_formula = mineral_formula
            
            # remove entities is an ordered fashion
            mineral_formula = re.sub('(H2O|OH|CO3)', '', mineral_formula)
            if 'S(6)' in self.parameters['solution_elements']:
                mineral_formula = mineral_formula.replace('SO4', '')
            mineral_formula = re.sub('([0-9()•:.])', '', mineral_formula)
            mineral_elements = re.findall('[A-Z][a-z]?', mineral_formula)
            
            for element in mineral_elements:
                if all(element in self.parameters['solution_elements'] for element in mineral_elements):
                    self.variables['described_minerals'][mineral] = self.minerals[mineral]
                elif element not in self.parameters['solution_elements']:
                    if self.verbose:
                        print(f'--> The {element} element of < {mineral} | {original_formula} > is undefined in the feed.')
                        

    def _define_elements(self):
        elements_lines, undefined_elements = [], []
        self.predicted_effluent = {}
        for element, information2 in self.water_body['element'].items():
            self.parameters['solution_elements'].append(element)
            self.predicted_effluent[element] = information2['concentration (ppm)']*self.cumulative_cf
            
            ref = form = ''
            if 'reference' in information2:
                ref = information2['reference']
            if 'form' in information2:
                form = 'as '+information2['form']

            if element in self.elements:                          
                elements_lines.append(f'{element}\t\t{information2["concentration (ppm)"]}\t{form}\t#{ref}')
            else:
                undefined_elements.append(element)
        if undefined_elements != []:
            print('\n--> The {} elements are not accepted by the {} database'.format(undefined_elements, self.parameters['database_selection']))
        return elements_lines


    def feed_geochemistry(self, water_selection = '', water_characteristics = {}, solution_description = '', ignored_minerals = [], existing_scale = {}, parameterized_ph_charge = True):
        # create the solution line of the input file
        if water_selection not in self.feed_sources and water_characteristics == {}:
            self._error(f'The {water_selection} does not have a corresponding parameter file. Either select an existing feed water {self.feed_sources}, or define a custom feed water with the water_characteristics dictionary.', 'value')
        
        self.parameters['water_selection'] = water_selection
        if water_selection == '':
            self.parameters['water_selection'] = solution_description
            
        initial_number = 0
        if self.parameters['simulation_type'] == 'evaporation':
            initial_number = 1
        
        self.results['solution_block'] = []
        self.results['solution_block'].append(f'''
SOLUTION {initial_number}\t{self.parameters["water_selection"]}'''
                    )

        #========================= FEED IONS =================================           
        self.parameters['solution_elements'] = []
        temperature = ph = alkalinity = pe = None
        self.water_body = {}
        if water_selection in self.feed_sources:  
            with open(os.path.join(self.parameters['root_path'], 'water_bodies', f'{water_selection}.json')) as water_body:
                self.water_body = json.load(water_body)
        if water_characteristics != {}:
            if self.water_body == {}:
                self.water_body = water_characteristics
            else:
                for content, information in water_characteristics.items():
                    if content == 'element':
                        for element, information2 in information.items():
                            if 'concentration (ppm)' in information2:
                                self.water_body['element'][element]['concentration (ppm)'] = information2['concentration (ppm)']
                            if 'reference' in information2:
                                self.water_body['element'][element]['reference'] = information2['reference']
                            if 'form' in information2:
                                self.water_body['element'][element]['form'] = information2['form']
                    elif content == 'temperature (C)':     
                        if 'value' in information:
                            self.water_body['element']['temperature (C)']['value'] = information['value']
                        if 'reference' in information:
                            self.water_body['element']['temperature (C)']['reference'] = information['reference']
                    elif content == 'pe':
                        if 'value' in information:
                            self.water_body['element']['pe']['value'] = information['value']
                        if 'reference' in information:
                            self.water_body['element']['pe']['reference'] = information['reference']
                    elif content == 'Alkalinity':
                        if 'value' in information:
                            self.water_body['element']['Alkalinity']['value'] = information['value']
                        if 'reference' in information:
                            self.water_body['element']['Alkalinity']['reference'] = information['reference']
                        if 'form' in information:
                            self.water_body['element']['Alkalinity']['form'] = information['form'] 
                    elif content == 'pH':
                        if 'value' in information:
                            self.water_body['element']['pH']['value'] = information['value']
                        if 'reference' in information:
                            self.water_body['element']['pH']['reference'] = information['reference']
            
        for content, information in self.water_body.items():
            alkalinity_form = temperature_reference = pe_reference = alkalinity_reference = ph_reference = ''
            if content == 'element':
                elements_lines = self._define_elements()     
                elements_line = '\n'.join(elements_lines)
            elif content == 'temperature (C)':
                temperature = information['value']
                if 'reference' in information:
                    temperature_reference = information['reference']
            elif content == 'pe':
                pe = information['value']
                if 'reference' in information:
                    pe_reference = information['reference']
            elif content == 'Alkalinity':
                alkalinity = information['value']
                if 'reference' in information:
                    alkalinity_reference = information['reference'] 
                if 'form' in information:
                    alkalinity_form = 'as {}'.format(information['form'])    
            elif content == 'pH':
                ph = information['value']
                if 'reference' in information:
                    ph_reference = information['reference']
            
        # parameterize the lines of the SOLUTIONS block
        temperature_line = pe_line = alkalinity_line = ph_line = ''
        if temperature is not None:
            temperature_line = f'temp \t {temperature} \t #{temperature_reference}.'
        if pe is not None:
            pe_line = f'pe \t\t {pe} \t   #{pe_reference} // 4.00 is the default (?)'     
        if alkalinity:
            alkalinity_line = f'Alkalinity \t {alkalinity} {alkalinity_form} #{alkalinity_reference}'
        if ph is not None:
            ph_line = f'pH \t\t {ph} #{ph_reference}'
            if parameterized_ph_charge and not alkalinity:
                ph_line = f'pH \t\t {ph} charge #{ph_reference}'
                alkalinity_line = ''
                   
        unit_line = 'units \t ppm' 
        if water_selection == 'Bakken formation':
            water_line = '-water \t{}\t#TDS=300 per mille [before fudging]'.format(self.variables['feed_kg'])
        elif water_selection == 'German Basin':
            water_line = '-water \t{}\t#TDS=314 per mille [before fudging]'.format(self.variables['feed_kg'])
        else:
            water_line = f'-water \t{self.variables["feed_kg"]}'

        self.results['solution_block'].extend([temperature_line, ph_line, pe_line, alkalinity_line, unit_line, elements_line, water_line])

        #parameterize the initial module solution
        if self.parameters['simulation_type'] == 'transport':
            total_cells = self.parameters['simulated_cells']
            if self.parameters['domain'] == 'dual':
                total_cells = self.parameters['simulated_cells']*2+1
            feed_solution_line = f'\nSOLUTION 1-{total_cells}\tInitial solution in the RO module'
            water_line = '-water \t{}'.format(self.variables['feed_kg'])
            self.results['solution_block'].extend([feed_solution_line,'temp \t 25','units \t ppm',water_line])

        #========================= POTENTIAL SCALANTS =================================
        short_mineral_names = ['Barite', 'Gypsum','Halite', 'Natron', 'Quartz', 'Talc','Trona', 'Borax', 'Albite', 'K-mica','Illite', 'Pyrite', 'Sulfur',]
#        long_mineral_names = ['Anthophyllite', 'Hexahydrite', 'Leonhardite', 'Nesquehonite', 'Pentahydrite', 'Portlandite','Sepiolite(d)', 'Boric_acid,s', 'K2B4O7:4H2O', 'NaB5O8:5H2O', 'Rhodochrosite', 'Strontianite','Hydroxyapatite', 'Chlorite(14A)', 'Mackinawite', 'Hausmannite', 'Pyrochroite']

        # define the equilibrium_phases block
        self.results['equilibrium_phases_block'] = []
        equilibrium_phases_number = '1'
        if self.parameters['simulation_type'] == 'transport':
            equilibrium_phases_number = '1-'+str(total_cells)
        self.results['equilibrium_phases_block'].append(f'\nEQUILIBRIUM_PHASES {equilibrium_phases_number}')

        # define the equilibrium_phases lines for the code block
        self._described_minerals()
        for possible_mineral in self.variables['described_minerals']:
            if possible_mineral not in ignored_minerals:
                if possible_mineral in short_mineral_names:
                    mineral_line = f'{possible_mineral}\t\t' 
                elif possible_mineral == 'Ca-Montmorillonite':
                    mineral_line = possible_mineral
                else:
                    mineral_line = f'{possible_mineral}\t'

                if possible_mineral in existing_scale:
                    for key, value in existing_scale[possible_mineral].items():
                        if key == 'saturation':
                            mineral_line += f'\t{value["saturation"]}'
                        if key == 'initial_moles':
                            mineral_line += f'\t{value["initial_moles"]}'
                else:
                    mineral_line += f'\t0\t0'
                    
                self.results['equilibrium_phases_block'].append(mineral_line)     
                
        if self.verbose:
            pprint(self.variables['described_minerals'])
                        
    def _selected_output(self, output_filename):
        # create parameter lines             
        if output_filename is None:
            count = 0
            self.parameters['selected_output_file_name'] = '-'.join([str(x) for x in [
                    datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], count
                    ]]) 
            while os.path.exists(self.parameters['selected_output_file_name']+'.txt'):
                count += 1
                self.parameters['selected_output_file_name'] = '-'.join([str(x) for x in [
                        datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], count
                        ]]) 
        else:
            self.parameters['selected_output_file_name'] = output_filename
        self.parameters['selected_output_file_name'] += '.txt'

        # establish the selected_output_block
        minerals_line = ' '.join([mineral for mineral in self.variables['described_minerals']])
        elements_line = ' '.join([element for element in self.parameters['solution_elements']])
        self.results['selected_output_block'] = [f'''
SELECTED_OUTPUT
-file\t\t{self.parameters["selected_output_file_name"]}
-reaction\t\ttrue
-temperature\ttrue
-totals\t\t{elements_line}
-saturation_indices\t{minerals_line}
-equilibrium_phases\t{minerals_line}
-pH\t\t\ttrue
-alkalinity\ttrue
-solution
-time\t\ttrue
-distance\t\ttrue
-simulation\t\ttrue
-high_precision\ttrue
-charge_balance\ttrue
-ionic_strength\ttrue
-step
-water
        ''']
        
    def _define_paths(self, simulation_name = None, simulation_directory = None):
        # define the simulation name 
        simulation_number = 1
        if simulation_name is None:
            if self.parameters['simulation_type'] == 'evaporation':
                simulation_name = '-'.join([str(x).replace(' ', '_') for x in [
                        datetime.date.today(), 'ROSSpy', self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], self.parameters['simulation_perspective']
                        ]])
            elif self.parameters['simulation_type'] == 'transport':
                permeate_approach_name = 'LinPerm'
                if self.parameters['permeate_approach'] == 'linear_cf':
                    permeate_approach_name = 'LinCF'
                simulation_name = '-'.join([str(x).replace(' ', '_') for x in [
                        datetime.date.today(), 'ROSSpy', self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['simulation'], self.parameters['simulation_perspective'], permeate_approach_name
                        ]])
            
        directory = os.getcwd()
        if simulation_directory is not None:
            directory = simulation_directory
            
        while os.path.exists(os.path.join(directory, simulation_name)):
            simulation_number += 1
            simulation_name = re.sub(r'(-\d+$)', '', simulation_name)
            simulation_name = '-'.join([simulation_name, str(simulation_number)])
            
        # define the simulation path 
        self.simulation_path = os.path.join(directory, simulation_name)
        if self.export_content:
            if not os.path.exists(self.simulation_path):
                os.mkdir(self.simulation_path)
        
        # define the input file export path
        self.parameters['input_path'] = os.path.join(self.simulation_path, 'input.pqi')                       
        self.parameters['output_path'] = os.path.join(self.simulation_path, 'selected_output.csv')
        self.parameters['simulation_path'] = self.variables['simulation_path'] = self.simulation_path
        
        return directory, simulation_name
        
    def parse_input(self, input_file_path, water_selection = None, active_m2 = None):        
        # identify the database in the name
        with open(input_file_path, 'r') as file:
            self.input_file = file.read()
            
        for db in self.databases:
            if db in input_file_path:
                self.parameters['database_selection'] = db
            
        # parse the input_file       
        input_df = pandas.read_table(input_file_path, names = ['content'], sep='\n').squeeze()
        self.parameters['simulation_type'], self.parameters['permeate_approach'] = 'evaporation', 'linear_permeate'
        self.parameters['domain'], self.parameters['water_selection'] = 'single', water_selection
        
        self.predicted_effluent, self.parameters['solution_elements'] = {}, []
        permeate_moles = 0
        for index, row in input_df.iteritems():
            if 'DATABASE' in row:
                row = row.replace('DATABASE ','')
                self.parameters['database_selection'] = os.path.basename(row).replace('.dat', '')
            elif re.search(r'(^[A-Z][a-z \(]?[\s\d])', row):
                element = re.search(r'(^[A-Z][a-z]?)', row).group()
                conc = re.search(r'(?<!\()([\d\.]+)', row).group(1)
                if float(conc) > 0:
                    self.parameters['solution_elements'].append(element)
                    self.predicted_effluent[element] = float(conc)
            elif re.search(r'-water\s+\d', row):
                water_mass = float(re.sub(r'(-water\s+)', '', row))
            elif 'linear_cf' in row:
                self.parameters['permeate_approach'] = 'linear_cf'
            elif 'EQUILIBRIUM_PHASES' in row:
                row = row.replace('EQUILIBRIUM_PHASES ','') 
                row = float(re.sub(r'(\-\d+)', '', row))
                if row != 1:
                    self.parameters['domain'] = 'dual'
                    self.parameters['domain_phase'] = 'Mobile'
            elif 'H2O -1; ' in row:
                moles = re.search(r'([0-9]+\.[0-9]+)', row).group()
                permeate_moles += float(moles)
            elif 'TRANSPORT' in row:
                self.parameters['simulation_type'] = 'transport'
            elif '-cells' in row:
                self.parameters['cells_per_module'] = float(re.sub(r'-cells\s+', '', row))
            elif '-shifts' in row:
                self.simulation_shifts = float(re.sub(r'-shifts\s+', '', row))
            elif '-file' in row:
                self.parameters['selected_output_file_name'] = re.sub(r'-file\s+', '', row)
            elif '-time_step' in row:
                self.parameters['timestep'] = float(re.sub(r'-time_step\s+', '', row).split('\t')[0])
                self.parameters['simulation_time'] = self.parameters['timestep'] * self.simulation_shifts
            elif '-punch_cells' in row:
                row = re.sub(r'-punch_cells\s+','',row)
                row = float(re.sub(r'(-\d+)', '', row))
                if row != 1:
                    self.parameters['domain'] = 'dual'
                    self.parameters['domain_phase'] = 'Immobile'
            elif '-punch_frequency' in row:
                row = float(re.sub(r'-punch_frequency\s+','',row))
                self.parameters['simulation_perspective'] = 'all_distance'               
                if row == 1:
                    self.parameters['simulation_perspective'] = 'all_time'
                    
        # open the respective database
        self.parameters['database_path'] = os.path.join(self.parameters['root_path'], 'databases', self.parameters['database_selection']+'.dat')
        self._define_database()
        self._described_minerals()
                    
        # predict the effluent concentrations
        self.cumulative_cf = water_mass/(water_mass - permeate_moles*self.water_mw*milli)
        for element in self.predicted_effluent:
            self.predicted_effluent[element] = self.predicted_effluent[element]*self.cumulative_cf
                
        # define the simulation folder and parameters
        if self.parameters['domain'] == 'dual':
            self.figure_title = f'{self.parameters["domain_phase"]} phase {self.parameters["simulation"]} from the {self.parameters["water_selection"]} after {sigfigs_conversion(self.parameters["simulation_time"])}'
        self.parameters['active_m2'] = active_m2
        if self.parameters['active_m2'] is None:
            self.parameters['active_m2'] = 37
        self.parameters['active_m2_cell'] = self.parameters['active_m2']/self.parameters['cells_per_module']
        
        
    def execute(self, simulation_name = None, selected_output_path = None, simulation_directory = None, figure_title = None, title_font = 'xx-large', label_font = 'xx-large', x_label_number = 6, export_name = None, export_format = 'svg', scale_ions = True, define_paths = True, selected_output_filename = None):
        '''Execute a PHREEQC input file '''
        def run(input_file, first=False):
            try:
                phreeqc = self.phreeqc_mod.IPhreeqc()  
            except:
                raise ModuleNotFoundError('The IPHREEQC module has not yet been installed. The IPHREEQC module (https://water.usgs.gov/water-resources/software/PHREEQC/index.html) must be installed before ROSSpy can be used, which is detailed in the ROSSpy docs: https://rosspy.readthedocs.io/en/latest/?badge=latest#installation .')
            phreeqc.load_database(self.parameters['database_path'])
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
            def measure_time(func, *args, **kwargs):
                start = timeit.default_timer()
                phreeqc, conc = func(*args, **kwargs)
                return phreeqc, conc, timeit.default_timer() - start

            phreeqc, conc, run_time = measure_time(run, input_file)            
            
            # export the simulation results
            headers = conc.keys()
            self.selected_output = pandas.DataFrame(conc, columns = headers)
            if self.verbose:
                if self.jupyter:
                    pandas.set_option('display.max_columns', None)
                    display(self.selected_output)
                else:
                    print(self.selected_output)
            self.variables['run_time (s)'] = float(run_time)
       
        # construct the SELECTED_OUTPUT file
        if define_paths:
            self._define_paths(simulation_name, simulation_directory)
        self._selected_output(selected_output_filename)
        if 'solution_block' in self.results:
            self.results['complete_lines'] = chain(
                    [f'# {self.simulation_path}'], self.results['general_conditions'], self.results['solution_block'], self.results['equilibrium_phases_block'], self.results['reaction_block'], self.results['selected_output_block'], self.results['transport_block']
                     ) 
            self.input_file = '\n'.join([line for line in self.results['complete_lines']])       
            
        if self.export_content:
            with open(self.parameters['input_path'], 'w') as input:
                input.write(self.input_file)
          
        # communicate estimated completion to the user
        estimated_time = (self.parameters['simulation_time']/(2*self.parameters['timestep'])) + 10*(self.parameters['quantity_of_modules'])
        if self.parameters['database_selection'] == 'sit':
            estimated_time *= 14
        if 'coarse_timestep' in self.parameters:
            if self.parameters['coarse_timestep']:
                estimated_time /= 6
        estimated_completion = datetime.datetime.now() + datetime.timedelta(seconds = estimated_time)
        estimated_time, unit = time_determination(estimated_time)
        print(f'\nEstimated completion in {estimated_time} {unit}: {estimated_completion} local time.')
        
        # execute the simulation or parse the selected_output file
        if selected_output_path is None:
            main(self.input_file)
            print('run_time (s):', self.variables['run_time (s)'])
            selected_output_path = os.path.join(self.simulation_path, self.parameters['selected_output_file_name'])  
            
            # verify that the PHREEQC executed and generated the appropriate files
            if self.export_content:
                self.selected_output.to_csv(self.parameters['output_path'])
                if not os.path.exists(self.parameters['output_path']):
                    self._error(f'The {self.parameters["output_path"]} output file does not exist. The simulation likely failed to execute.', 'index')
        else:
            # define the simulation perspective
            self.parameters['simulation'] = 'scaling'
            if re.search('brine', selected_output_path, flags=re.IGNORECASE):
                self.parameters['simulation'] = 'brine'

            # define the simulation type
            self.parameters['simulation_type'] = 'transport'
            for line in selected_output_path:
                if re.search('evaporation', line, re.IGNORECASE):
                    self.parameters['simulation_type'] = 'evaporation'

            # define the database contents
            self.parameters['database_selection'] = 'pitzer'
            for database in self.databases:
                names = database.split('_')
                if all(re.search(name, selected_output_path, re.IGNORECASE) for name in names):
                    self.parameters['database_selection'] = database                   
                    break                    

            # define the simulation
            self.parameters['selected_output_file_name'] = re.search(r'(\w+)(?=\.)', selected_output_path).group()
                                  
            # preparing the SELECTED_OUTPUT file into a dataframe
            sep = '\t'
            if not 'pqo' in selected_output_path:
                sep = ','
            with open(selected_output_path, 'r') as selected_output:              
                self.selected_output = pandas.read_table(selected_output, sep = sep)
                
            for column in self.selected_output.columns:
                new_column = column.strip()
                self.selected_output.rename(columns={column:new_column}, inplace = True)
            
        # process the data
        if figure_title is not None:
            self.figure_title = figure_title

        self.variables['initial_solution_mass'] = float(self.selected_output['mass_H2O'][0])
        if self.parameters['simulation_type'] == 'transport':
            self.selected_output.drop(self.selected_output.index[:3], inplace=True)
            self.variables['initial_solution_mass'] = float(self.selected_output['mass_H2O'][3])
        self.variables['final_solution_mass'] = float(self.selected_output['mass_H2O'].iloc[-1])
        self.variables['simulation_cf'] = self.variables['initial_solution_mass'] / self.variables['final_solution_mass']
        self.variables['final_time'] = float(self.selected_output['time'].iloc[-1])
                                 
        # conducting the appropriate visualization function
        if self.parameters['simulation'] == 'brine':
            self.processed_data = self._brine_plot(title_font, label_font, export_format, x_label_number, export_name)
            csv_export_name = 'brine_concentrations.csv'
        elif self.parameters['simulation'] == 'scaling':
            self.processed_data = self._scaling_plot(title_font, label_font, x_label_number, export_name, export_format)
            csv_export_name = 'scaling_data.csv'
            if scale_ions and self.processed_data is not None:
                self._ion_proportions()
        else:
            self._error('The < simulation_perspective > perspective is not supported.', 'type')
            
        # export the simulation data
        if self.export_content:            
            if self.processed_data is not None:
                self.processed_data.to_csv(os.path.join(self.simulation_path, csv_export_name))
            
            # predicted effluent concentrations
            effluent_concentrations = pandas.DataFrame(list(self.predicted_effluent.items()), columns = ['elements', 'ppm'])
            effluent_concentrations.index = effluent_concentrations['elements']
            del effluent_concentrations['elements']
            effluent_concentrations.to_csv(os.path.join(self.simulation_path, 'effluent_predictions.csv'))
            
            # parameters
            parameters = {'parameter':[], 'value':[]}
            for parameter in self.parameters:
                parameters['parameter'].append(parameter)
                parameters['value'].append(self.parameters[parameter])
                
            parameters_table = pandas.DataFrame(parameters)
            parameters_table.to_csv(os.path.join(self.simulation_path, 'parameters.csv'))

            # variables
            variables = {'variable':[], 'value':[]}
            for variable in self.variables:
                variables['variable'].append(variable)
                variables['value'].append(self.variables[variable])
                
            variables_table = pandas.DataFrame(variables)
            variables_table.to_csv(os.path.join(self.simulation_path, 'variables.csv'))
                    
        if self.verbose:
            print(self.input_file)
            if self.jupyter:
                display(variables_table)
                display(parameters_table)
            else:
                print(variables_table)
                print(parameters_table)
                
        return self.processed_data
    
    def _ion_proportions(self):
        # calculate the mass of each scalant
        mineral_elements = {}
        for column in self.processed_data:                   
            mineral = re.search('([A-Za-z]+)', column).group()
            mineral_elements[mineral] = {}
            self.chem_mw.mass(self.minerals[mineral]['formula'])
            mineral_elements[mineral]['proportions'] = self.chem_mw.proportions

        # calculate the precipitated mass of each element
        self.elemental_masses = {}
        for mineral in mineral_elements:
            index, unit, x_axis = 0, '(g/m^2)', '(m)'
            if self.parameters['simulation_type'] == 'evaporation':
                unit = '(g)'
            if self.parameters['simulation_perspective'] == 'all_time':
                x_axis = '(s)'
            for scale in self.processed_data[f'{mineral} {unit}']:  
                index_value = float(self.processed_data.index[index])
                if self.parameters['simulation_perspective'] == 'all_time':
                    index_value = int(index_value)
                df_index = f'{index_value} {x_axis}'
                if df_index not in self.elemental_masses:
                    self.elemental_masses[df_index] = {}
                    self.elemental_masses[df_index]['ion (g/m^2)'] = {}
                for element in mineral_elements[mineral]['proportions']:
                    precipitated_mass = mineral_elements[mineral]['proportions'][element] * float(scale)
                    if element in self.elemental_masses[df_index]['ion (g/m^2)']:
                        self.elemental_masses[df_index]['ion (g/m^2)'][element] += precipitated_mass
                    elif element in self.elemental_masses[df_index]:
                        self.elemental_masses[df_index][element] += precipitated_mass
                    else:
                        self.elemental_masses[df_index]['ion (g/m^2)'][element] = precipitated_mass
                index += 1

        # print and export the results
        if self.printing:
            print('ionic masses\n', '='*len('ionic masses  '))
            pprint(OrderedDict(self.elemental_masses))
            
        if self.export_content:
            with open(os.path.join(self.simulation_path, 'scale_ions.json'), 'w') as output:
                json.dump(self.elemental_masses, output, indent = 4)
                
    def _error(self, error, error_type):
        try:
            export_path = self.simulation_path 
        except:
            export_path = os.getcwd()
            
        with open(os.path.join(export_path, 'error.txt'), 'w') as output:
            output.write(error)
            output.close()
            
        if error_type == 'index':
            raise IndexError(error)
        elif error_type == 'value':
            raise ValueError(error)
        elif error_type == 'type':
            raise TypeError(error)
            
                                 
    def _brine_plot(self, title_font, label_font, export_format, x_label_number, export_name, log_scale = True):
        """Generate plots of the elemental concentrations from effluent brine in the PHREEQC SELECTED_OUTPUT file"""
        # determine the minerals in the simulation      
        columns = []
        for column in self.selected_output.columns:
            if re.search(r'([A-Z][a-z]?(?:\(\d\))?(?=\(mol\/kgw\)))', column) and not re.search('(_|H2O|pH)', column):
                columns.append(column)

        # plot parameters
        if self.parameters['domain'] == 'dual':
            title_end = 'after '+sigfigs_conversion(self.parameters['simulation_time'])
            if self.parameters['simulation_perspective'] == 'all_time':
                title_end = 'over time'     
            self.figure_title = f'{self.parameters["domain_phase"].capitalize()} phase {self.parameters["simulation"]} from the {self.parameters["water_selection"]} {title_end}'
        
        pyplot.figure(figsize = (17,10))
        non_zero_elements, non_zero_columns, data = [], [], {} 
        insufficient_elements = set()
        for element in columns:
            stripped_element = re.search(r'([A-Z][a-z]?(?:\(\d\))?(?=\(mol\/kgw\)))', element).group()
            if self.selected_output[element].iloc[-1] > 0: 
                non_zero_elements.append(stripped_element)
                non_zero_columns.append(element)
                if max(self.selected_output[element]) < 1E-16:
                    insufficient_elements.add(element)
                
            concentration_serie, x_serie, data[element] = [], [], {}
            for index, row in self.selected_output.iterrows():
                if self.parameters['simulation_perspective'] == 'all_time':
                    if not all(row[element] > 1e-16 for element in non_zero_columns):
                        continue
                    time = float(sigfigs_conversion(row['time'], 3))
                    x_serie.append(time)                        
                    data[element][time] = row[element]
                elif self.parameters['simulation_perspective'] == 'all_distance':       
                    if row['time'] == 0:
                        continue
                    distance = float(sigfigs_conversion(row['dist_x'], 3))
                    x_serie.append(distance)                        
                    data[element][distance] = row[element]
                concentration_serie.append(row[element])
                
            pyplot.plot(x_serie,concentration_serie, label = stripped_element)
                        
        if self.parameters['simulation_perspective'] == 'all_time' and len(x_serie) == 0:
            self._error(f'The {insufficient_elements} elements remain below the 1E-16 molal threshold, and thus the figure could not be constructed.', 'index')
                        
        # define the brine plot
        x_label = 'Distance (m)'
        if self.figure_title is None:
            self.figure_title = f'Brine concentrations along the module after {sigfigs_conversion(self.variables["final_time"])} seconds'
        if self.parameters['simulation_perspective'] == 'all_time':
            x_label = 'Time (s)'
            if self.figure_title is None:
                self.figure_title = 'Effluent brine concentration over time' 
                
        if self.parameters['simulation_perspective'] == 'all_time':
            x_location = [x_serie[0]]
            index_space = len(x_serie)/x_label_number
            for x in range(x_label_number):
                index = int(index_space * x)
                x_location.append(x_serie[index])
            x_location.append(x_serie[-1])
            pyplot.xticks(x_location)
        self._illustrate(pyplot, 'non-zero ions', export_name, label_font, title_font, export_format, log_scale)

        # defining the datatable of brine concentrations
        concentrations_table = pandas.DataFrame(data)
        concentrations_table.index.name = x_label
        if self.printing:
            if self.jupyter:
                display(concentrations_table)
            else:
                print(concentrations_table)
        return concentrations_table            

    def _scaling_plot(self, title_font, label_font, x_label_number, export_name, export_format = 'svg', log_scale = None):
        """Generate plots of scaling along the module distance in the PHREEQC SELECTED_OUTPUT file"""
        def evaporation():
            scaling_data = pandas.DataFrame({})
            pyplot.figure(figsize = (17,10))
            data = {}
            for mineral in self.variables['precipitated_minerals']: 
                mineral_formula = self.minerals[mineral]['formula']
                data[f'{mineral} (g)'], cf_series, scaling_series = {}, [], []
                for index, row in self.selected_output.iterrows():
                    if index != len(self.selected_output['mass_H2O']):
                        if row['step'] >= 1:
                            cf = self.variables['initial_solution_mass'] / row['mass_H2O']
                            scale_mass = row[mineral] * float(self.minerals[mineral]['mass'])
                            scaling_series.append(scale_mass) 
                            cf_series.append(cf)   
                            data[f'{mineral} (g)'][float(sigfigs_conversion(cf, 6))] = scale_mass            

                pyplot.plot(cf_series,scaling_series, label = f'{mineral} [{mineral_formula}]')
                data_df = pandas.DataFrame(data)
                scaling_data = scaling_data.merge(data_df, how = 'outer', left_index = True, right_index = True)
                    
            if self.figure_title is None:
                self.figure_title = 'Evaporation scaling from the {}'.format(self.parameters['water_selection'])                    
            return scaling_data
        
        def time_serie(mineral, mineral_formula, molar_mass_area):
            for index, row in self.selected_output.iterrows():
                if row[mineral] > 1e-12:
                    x_serie.append(float(row['time'])) # - initial_solution_time * self.parameters['timestep'])
                    grams_area = float(row[mineral])*molar_mass_area
                    scaling_serie.append(float(sigfigs_conversion(grams_area, 3)))
                    
            # defining the plot
            x_location = []
            index_space = len(x_serie)/x_label_number
            for x in range(x_label_number):
                index = int(index_space * x)
                time = x_serie[index]
                x_location.append(time)
            x_location.extend([x_serie[0], x_serie[-1]])
            x_location = [int(x) for x in x_location] 
            pyplot.xticks(x_location, x_location)        
            pyplot.plot(x_serie,scaling_serie, label = f'{mineral} [{mineral_formula}]')
            
            data_df = pandas.DataFrame(
                    columns = [f'{mineral} (g/m^2)'],
                    data = [sigfigs_conversion(y, 4) for y in array(scaling_serie)],
                    index = [sigfigs_conversion(str(x), 3) for x in array(x_serie)]
                    )
            
            return data_df
            
        def distance_serie(mineral, mineral_formula, molar_mass_area):          
            for index, row in self.selected_output.iterrows():
                grams_area = float(row[mineral])*molar_mass_area
                if row['time'] != 0:
                    scaling_serie.append(float(sigfigs_conversion(grams_area, 3)))
                    x_serie.append(row['dist_x'])
                     
            # define the figure and dataframe
            pyplot.plot(x_serie,scaling_serie, label = f'{mineral} [{mineral_formula}]')
            data_df = pandas.DataFrame(
                    columns = [f'{mineral} (g/m^2)'],
                    data = [sigfigs_conversion(y, 4) for y in array(scaling_serie)],
                    index = [sigfigs_conversion(x, 3) for x in array(x_serie)]
                    )
                    
            if self.printing:
                if len(mineral) < 7:
                    print(f'{mineral}\t\t{self.minerals[mineral]}') 
                elif len(mineral) >= 7:
                    print(f'{mineral}\t{self.minerals[mineral]}') 
            return data_df
        
        # all of the non-zero minerals are identified and the chemical formulas are sorted into a list
        csv_minerals, non_zero_minerals, self.variables['precipitated_minerals'] = [], set(), {}
        maximum_scaling, minimum_scaling = -inf, 0
        for column in self.selected_output.columns:
            if re.search('([A-Z][a-z]{2,})', column) and not re.search('[_(]|(?:Metal)', column):
                mineral = re.search('([A-Z][a-z]{2,})', column).group()
                csv_minerals.append(mineral)
                
                # the spread of scaling is determined to select either a linear or logarithmic y-axis
                if max(self.selected_output[column]) != 0:
                    maximum_scaling = max(max([y for y in self.selected_output[column] if y>0]), maximum_scaling)
                    minimum_scaling = min(min([y for y in self.selected_output[column] if y>0]), minimum_scaling)
                    non_zero_minerals.add(mineral)   
                    if column in self.minerals:
                        self.variables['precipitated_minerals'][mineral] = self.minerals[mineral]
        
        if non_zero_minerals == set():
            warn('No scaling occurred.')
            return None          
        else:
            if log_scale is None:
                log_scale = False
                if log10(maximum_scaling)-log10(minimum_scaling) > 1:
                    log_scale = True

        # finalize the output data
        scaling_data = pandas.DataFrame({})
        pyplot.figure(figsize = (17,10))
        for mineral in self.variables['precipitated_minerals']: 
            scaling_serie, x_serie = [], []
            mineral_formula = self.minerals[mineral]['formula'] 
            molar_mass_area = float(self.minerals[mineral]['mass']) / self.parameters['active_m2_cell']
            
            if self.parameters['simulation_type'] == 'evaporation':
                data_df = evaporation()
            elif self.parameters['simulation_perspective'] == "all_time":
                data_df = time_serie(mineral, mineral_formula, molar_mass_area)
            elif self.parameters['simulation_perspective'] == 'all_distance':
                data_df = distance_serie(mineral, mineral_formula, molar_mass_area)
            
            scaling_data = scaling_data.merge(data_df, how = 'outer', left_index = True, right_index = True)
        
        if self.figure_title is None:
            time, units = time_determination(self.variables['final_time'])
            self.figure_title = f'Scaling from the {self.parameters["water_selection"]} after {sigfigs_conversion(time)} {units}'
            
        self._illustrate(pyplot, 'scale', export_name, label_font, title_font, export_format, log_scale)
        x_label, y_label = self._determine_labels()
        scaling_data.index.name = x_label
        return scaling_data
    
    def _determine_labels(self): 
        # determine the y-axis label
        if self.parameters['simulation'] == 'scaling':
            y_label = 'Mass concentration (g/m^2)'
            if self.parameters['simulation_type'] == 'evaporation':
                y_label = 'Mass (g)'
        elif self.parameters['simulation'] == 'brine':
            y_label = 'Concentration (molal)'
        
        # determine the x-axis label
        if self.parameters['simulation_type'] == 'transport':
            x_label = 'Distance (m)'
            if self.parameters['simulation_perspective'] == 'all_time':
                x_label = 'Time (s)'
        elif self.parameters['simulation_type'] == 'evaporation':
            x_label = 'Concentration Factor (CF)'
        return x_label, y_label
    
    def _illustrate(self, pyplot, legend_title, export_name, label_font, title_font, export_format, log_scale):
        def export_plot(figure, export_name = None, export_format = 'svg'):
            if export_name is None:
                export_name = 'all_minerals'
                if self.parameters['simulation'] == 'brine':
                    export_name = 'brine'
            self.results['figures'][export_name] = {'figure':figure, 'title':self.figure_title}

            # export the plot
            file_number, figure_path = 0, os.path.join(self.simulation_path, export_name)
            if not os.path.exists(f'{figure_path}.{export_format}'):
                figure.savefig(f'{figure_path}.{export_format}')
            elif os.path.exists(f'{figure_path}.{export_format}'):
                while os.path.exists(f'{figure_path}_{file_number}.{export_format}'):
                    file_number += 1
                figure.savefig(f'{figure_path}_{file_number}.{export_format}')
                
        # apply the attributes of the figure
        x_label, y_label = self._determine_labels()
        pyplot.grid(True)
        pyplot.title(self.figure_title, fontsize = title_font)
        pyplot.xlabel(x_label, fontsize = label_font)
        pyplot.ylabel(y_label, fontsize = label_font)  
        pyplot.legend(title = legend_title, loc='best', title_fontsize = 'x-large', fontsize = 'large')  
        if self.parameters['simulation_type'] == 'transport':
            pyplot.figtext(0.2, 0.07, f'Final CF: {float(sigfigs_conversion(self.variables["simulation_cf"], 4))}', 
                                                   wrap=True, horizontalalignment='left', fontsize=12)
        if log_scale:
            pyplot.yscale('log')
        
        figure = pyplot.gcf()
        if self.printing:
            pyplot.show()
        if self.export_content:
            export_plot(figure, export_name, export_format)

    def test(self):
        self.reactive_transport(simulation_time = 200)
        self.feed_geochemistry(water_selection = 'red_sea')
        self.execute(simulation_name = f'rosspy_test_{self.parameters["os"]}')