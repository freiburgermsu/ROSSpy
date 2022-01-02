# print the welcome message
from itertools import chain
print('\n\n')
message = ('''* iROSSpy: Interactive Reverse Osmosis Scaling Simulation in Python *
Andrew Philip Freiburger, Green Safe Water Lab, University of Victoria''')

indent = 1
lines = message.split('\n')
space = " " * indent
width = max(map(len, lines))

upper = f'╔{"═" * (width + indent * 2)}╗\n'
middles = ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
lower = f'╚{"═" * (width + indent * 2)}╝' 

box = chain(upper, middles, lower)
box_print = ''.join(box)
print(box_print, '\n')

print('The program is loading. It may take a few minutes...\n\n')

# import libraries
from scipy.constants import nano, milli, minute
import subprocess
import rosspy
import pandas
import time, os, re

# create the input file
class iROSSpy():
    def __init__(self):
        self.working_directory = os.getcwd()
        
        # announcements
        announcement = 'Parameterize initial details:'
        print(announcement, '\n', '='*len(announcement))
        
        # define operating system
        operating_systems = ['windows', 'unix']
        operating_system =  input('''- Is your computer operating system windows or unix?
        < {} > or < {} > ; Default = < windows >  __ '''.format(operating_systems[0], operating_systems[1])) or 'windows'
        while operating_system not in operating_systems:
            print('''--> ERROR: The operating system is not one of the accepted options.''')
            operating_system =  input('''- Is your computer operating system windows or unix?
            < {} > or < {} >  __ '''.format(operating_systems[0], operating_systems[1])) 
    
        # define verbosity
        verbose =  input('''- Your simulation be verbose?
        < True > or < False > ; Default = < False >  __ ''') or False
        while verbose not in ['True', 'False', False]:
            print('''--> ERROR: The verbosity must be an True or False.''')    
            verbose =  input('''- Your simulation be verbose?
            < True > or < False > ; Default = < False >  __ ''') or False
        verbose = bool(verbose)
                             
        # define the maximal argument path length 
        self.max_argument_length = len(r"C:\Users\Andrew Freiburger\Dropbox\My PC (DESKTOP-M302P50)\Documents\UVic Civil Engineering\PHREEQC\ROSS\examples\scaling\scale_validation\2021-12-18-ROSSpy-red_sea-transport-pitzer-scaling-all_distance-LinPerm\2021-12-28-ROSSpy--transport-pitzer-scaling-all_")
            
        # initiate the ROSSpy object
        self.ross = rosspy.ROSSPkg(operating_system, verbose, jupyter = False)

    def define_general(self,):
        # announcement
        announcement = '\nParameterize general simulation conditions:'
        print(announcement, '\n', '='*len(announcement))       
                             
        # define database
        for database in self.ross.databases:
            print(f'< {database} >') 
        database_selection = input('''- What database do you select?
        Default = < pitzer >  __ ''') or 'pitzer'
        while database_selection not in ross.databases:
            print('''The database is not current defined in this interface.
            Only < pitzer > and < phreeqc > are defined for Windows and < pitzer > is defined for Macintosh. __ ''')    
            database_selection = input('- What database do you select?')
        
        # define simulation
        simulations = ['scaling', 'brine']
        simulation = input('''- Will you simulate scaling or brine?
        Default = < scaling >  __ ''') or 'scaling'
        while simulation not in simulations:
            print('''--> ERROR: Only scaling or brine formation simulations can be simulated through ROSSpy.''')    
            simulation = input('''- Will you simulate scaling or brine formation?
            Default = < scaling >  __ ''') or 'scaling'
            
        # define domain
        domains = ['single', 'dual']
        domain = input('''- Will you simulate the single or dual domain?
        Default = < single >  __ ''') or 'single'
        while domain not in domains:
            print('''--> ERROR: Only the single or dual domains can be simulated through ROSSpy.''')    
            domain = input('''- Will you simulate the single or dual domain?
            Default = < single >  __ ''') or 'single'
        
        # define domain_phase
        if domain == 'dual': 
            domain_phases = ['mobile', 'immobile']
            domain_phase = input('''- Will you evaluate the mobile or immobile phase?
            Default = < immobile >  __ ''') or 'immobile'
            while domain_phase not in domain_phases:
                print('''--> ERROR: Only the mobile or immobile phases can be simulated through ROSSpy.''')    
                domain_phase = input('''- Will you simulate the single or dual domain?
                Default = < single >  __ ''') or 'immobile'
        else:
            domain_phase = None
                             
        # define quantity_of_modules
        quantity_of_modules = input('''- How many in-series RO modules will you simulate?
        Default = < 1 >  __ ''') or 1  
        while True:
            try:
                quantity_of_modules = int(quantity_of_modules)
                break
            except:            
                print('''--> ERROR: The quantity of modules must be an integer.''')    
                quantity_of_modules = input('''- How many in-series RO modules will you simulate?
                Default = < 1 >  __ ''') or 1
        
        # define simulation_perspective
        simulation_types = ['transport','evaporation']
        simulation_type = input('''- Will you simulate evaporation or transport?
        Default = < transport >  __ ''') or 'transport'
        while simulation_type not in simulation_types:
            print('''--> ERROR: Only transport or evaporation can be simulated through ROSSpy.''')  
            simulation_type = input('''- Will you simulate evaporation or transport?
            Default = < transport >  __ ''') or 'transport'   
        
        # define simulation title
        simulation_title = input('- What is the title of your simulation?   __ ')   
                             
        # apply the parameters to the ROSSpy instance
        self.ross.define_general(database_selection = database_selection, simulation = simulation, domain = domain, domain_phase = domain_phase, quantity_of_modules = quantity_of_modules, simulation_type = simulation_type, simulation_title = simulation_title)


    def transport(self,):
        # announcement
        announcement = '\nParameterize transport conditions:'
        print(announcement, '\n', '='*len(announcement))                                  
                             
        # define simulation_time
        simulation_time = input('''- How many minutes will you simulate?
        Default = < 3 > minutes __ ''') or 3
        while True:
            try:
                simulation_time = int(simulation_time)
                break
            except:            
                print('''--> ERROR: The simulation_time must be an integer.''')    
                simulation_time = input('''- How many minutes will you simulate?
                Default = < 3 > minutes __ ''') or 3        
        simulation_time *= minute
        
        # define simulation_perspective
        simulation_perspectives = ['all_time', 'all_distance']
        if self.ross.parameters['simulation'] == 'scaling':
            default = 'all_time'
        if self.ross.parameters['simulation'] == 'brine':
            default = 'all_distance'
        simulation_perspective = input(f'''- Will you simulate all_time or all_distance?
        Default = < {default} >  __ ''') or default
        while simulation_perspective not in simulation_perspectives:
            print('''--> ERROR: Only {} can be simulated.'''.format([x for x in simulation_perspectives]))    
            simulation_perspective = input(f'''- Will you simulate all_time or all_distance?
            Default = < {default} >  __ ''') or default
            
        # define module_characteristics
        custom_module = input('''- You will simulate a custom module?
        < True > or < False > ; Default = < False >  __ ''') or False
        while custom_module not in ['True', 'False', False]:
            print('''--> ERROR: Only True or False are accepted responses.''')    
            custom_module = input('''- You will simulate a custom module?
            < True > or < False > ; Default = < False >  __ ''') or False
        custom_module = bool(custom_module)
        
        module_characteristics = {}
        if custom_module:
            default_characteristics = {
                'module_diameter_mm':201,
                'permeate_tube_diameter_mm':29,
                'module_length_m':1.016,
                'permeate_flow_m3_per_day':40,
                'max_feed_flow_m3_per_hour':15.9,
                'membrane_thickness_mm':250 * (nano / milli),
                'feed_thickness_mm':0.8636,
                'active_m2':37,
                'permeate_thickness_mm':0.3,
                'polysulfonic_layer_thickness_mm':0.05,
                'support_layer_thickness_mm':0.15
            }
            for characteristic in default_characteristics:
                module_characteristics[characteristic] = input(''' What is the value of the {} module characteric?
                default = {}  ___ '''.format(characteristic, default_characteristics[characteristic])) or default_characteristics[characteristic]
                while True:
                    try:
                        module_characteristics[characteristic] = int(module_characteristics[characteristic])
                        break
                    except:            
                        print('''--> ERROR: Only integers are accepted for the module_characteristic.''')    
                        module_characteristics[characteristic] = input(''' What is the value of the {} characteric?
                        default = {}'''.format(characteristic, default_characteristics[characteristic])) or default_characteristics[characteristic]
        
        # define cells_per_module
        cells_per_module = input('''- Into how many cells is your simulated module discretized?
        Default = < 12 > __ ''') or 12     
        while True:
            try:
                cells_per_module = int(cells_per_module)
                break
            except:            
                print('''--> ERROR: The cells_per_module must be an integer.''')    
                cells_per_module = input('''- Into how many cells is your simulated module discretized?
                Default = < 12 >  __ ''') or 12  
                             
        # define exchange_factor
        exchange_factor = ''
        if self.ross.parameters['domain'] == 'dual':
            exchange_factor = input('''- What is the (1/sec) exchange rate between the mobile and immobile phases of the simulated dual domain?
            Default = < 1e10 >  __ ''') or 1e10     
            while True:
                try:
                    exchange_factor = int(exchange_factor)
                    break
                except:            
                    print('''--> ERROR: The exchange_factor must be an integer.''')    
                    exchange_factor = input('''- What is the (1/sec) exchange rate between the mobile and immobile phases of the simulated dual domain?
                    Default = < 1e10 >  __ ''') or 1e10                                                        
                             
        # apply the parameters to the ROSSpy instance                 
        self.ross.transport(simulation_time = simulation_time, simulation_perspective = simulation_perspective, module_characteristics = module_characteristics, cells_per_module = cells_per_module, parameterized_timestep = None, kinematic_flow_velocity = None, exchange_factor = exchange_factor)


    def reaction(self,):
        # announcement
        announcement = '\nParameterize reaction conditions:'
        print(announcement, '\n', '='*len(announcement))      
        
        # define permeate_approach
        permeate_approaches = ['linear_permeate', 'linear_cf']
        for approach in permeate_approaches:
            print(f'< {approach} >')
        permeate_approach = input('''- Which permeate_approach will be simulated?
        Default = < linear_permeate >  __ ''') or 'linear_permeate'
        while permeate_approach not in permeate_approaches:
            print('''--> ERROR: Only {} are accepted parameters.'''.format([x for x in permeate_approaches]))    
            permeate_approach = input('''- Which permeate_approach will be simulated?
            Default = < linear_permeate >  __ ''') or 'linear_permeate'
                                                  
        # define permeate_efficiency
        permeate_efficiency = input('''- What is the efficiency of the simulated module?
        Default = < 1 >  __ ''') or 1 
        while True:
            try:
                permeate_efficiency = int(permeate_efficiency)
                break
            except:            
                print('''--> ERROR: The cells_per_module must be an integer.''')    
                permeate_efficiency = input('''- What is the efficiency of the simulated module?
                Default = < 1 >  __ ''') or 1 
                         
        # define head_loss
        head_loss = input('''- What is the relative pressure of the effluent to the feed?
        Default = < 0.89 >  __ ''') or 0.89
        while True:
            try:
                head_loss = float(head_loss)
                break
            except:            
                print('''--> ERROR: The head_loss must be a number.''')    
                head_loss = input('''- What is the relative pressure of the effluent to the feed?
                Default = < 0.89 >  __ ''') or 0.89
                         
        # define final_cf
        final_cf = None
        if permeate_approach == 'linear_cf':
            final_cf = input('''- What is the effluent concentration relative to the feed concentration?
            Default = < 2 >  __ ''') or 2
            while True:
                try:
                    final_cf = float(final_cf)
                    break
                except:            
                    print('''--> ERROR: The final_cf must be a number.''')    
                    final_cf = input('''- What is the effluent concentration relative to the feed concentration?
                    Default = < 2 >  __ ''') or 2
                         
        # apply the parameters to the ROSSpy instance      
        self.ross.reaction(permeate_approach = permeate_approach, permeate_efficiency = permeate_efficiency, head_loss = head_loss, final_cf = final_cf)


    def solutions(self,):
        # announcement
        announcement = '\nParameterize solutions:'
        print(announcement, '\n', '='*len(announcement))      
        
        # define water_selection
        water_files = glob(os.path.join(self.parameters['root_path'], 'water_bodies','*.json'))
        self.water_bodies = [re.search('([a-z\_]+)(?=\.)', file).group() for file in water_files].append('custom')
        for water in self.water_bodies:
            print(f'< {water} >')   
        water_selection = input('''- Which water body will you simulate?  __ ''')
        while water_selection not in self.water_bodies:
            print('''--> ERROR: One of the printed water_bodies may be parameterized.''')    
            water_selection = input('''- Which water body will you simulate?  __ ''')
        
        solution_description = None
        water_characteristics = {}
        if water_selection == 'custom':
            solution_description = input('- What is the description of the custom water source?')
            for element in self.ross.elements:
                if element == 'Alkalinity':
                    continue
                water_characteristics[element] = {}
                water_characteristics[element]['concentration (ppm)'] = input(f'''- What is ppm concentration of {element} in the feed?''')
                while True:
                    try:
                        water_characteristics[element]['concentration (ppm)'] = float(water_characteristics[element]['concentration (ppm)'])
                        break
                    except:            
                        print('''--> ERROR: Only numbers are accepted for elemental concentrations.''')    
                        water_characteristics[element]['concentration (ppm)'] = input(f'''- What is ppm concentration of {element} in the feed?''') 
                water_characteristics[element]['reference'] = input(f'- What is the reference for the {element} concentration?')
                         
            # define pH
            water_characteristics['pH'] = {}
            water_characteristics['pH']['value'] = input('''- What is the pH of the customized solution? ___ ''')
            while True:
                try:
                    water_characteristics['pH']['value'] = float(water_characteristics['pH']['value'])
                    break
                except:            
                    print('''--> ERROR: The pH must be a number.''')    
                    water_characteristics['pH']['value'] = input('''- What is the pH of the customized solution? ___ ''')
            water_characteristics['pH']['reference'] = input('- What is the reference for this pH?')

            # define temperature
            water_characteristics['temperature'] = {}
            water_characteristics['temperature']['value'] = input('''- What is the temperature of the customized solution? ___ ''')
            while True:
                try:
                    water_characteristics['temperature']['value'] = float(water_characteristics['temperature']['value'])
                    break
                except:            
                    print('''--> ERROR: The temperature must be a number.''')
                    water_characteristics['temperature']['value'] = input('''- What is the temperature of the customized solution? ___ ''')
            water_characteristics['temperature']['reference'] = input('- What is the reference for this temperature?')

            # define pe
            water_characteristics['pe'] = {}
            water_characteristics['pe']['value'] = input('''- What is the pe of the customized solution? ___ ''')
            while True:
                try:
                    water_characteristics['pe']['value'] = float(water_characteristics['pe']['value'])
                    break
                except:            
                    print('''--> ERROR: The pe must be a number.''')
                    water_characteristics['pe']['value'] = input('''- What is the pe of the customized solution? ___ ''')
            water_characteristics['pe']['reference'] = input('- What is the reference for this pe?')
                         
        # define pH charge balance
        ph_charge_balance = input('''- You will charge balance the pH?
        < True > or < False > ; Default = True ___ ''') or True
        water_characteristics['Alkalinity'] = {}
        if ph_charge_balance is True:
            water_characteristics['Alkalinity']['value'] = water_characteristics['Alkalinity']['reference'] = None
            parameterized_alkalinity = False
            parameterized_ph_charge = True
        else:
            water_characteristics['Alkalinity']['value'] = input('''- What is the alkalinity in ((eq of CaCO3)/(Kg of water))?''')
            water_characteristics['Alkalinity']['reference'] = input('''- What is the reference for the alkalinity?''')
            parameterized_alkalinity = True
            parameterized_ph_charge = False
            while True:
                try:
                    water_characteristics['Alkalinity']['value'] = float(water_characteristics['Alkalinity']['value'])
                    break
                except:            
                    print('''--> ERROR: The Alkalinity must be a number.''')
                    water_characteristics['Alkalinity']['value'] = input('''- What is the Alkalinity of the customized solution? ___ ''')
                         
        # apply the parameters to the ROSSpy instance      
        self.ross.solutions(water_selection = water_selection, water_characteristics = water_characteristics, solution_description = solution_description, parameterized_alkalinity = parameterized_alkalinity, parameterized_ph_charge = parameterized_ph_charge)


    def equilibrium_phases(self,):
        # announcement
        announcement = '\nParameterize evaluated scale:'
        print(announcement, '\n', '='*len(announcement))      
        
        # define ignored minerals
        self.ross.described_minerals()
        ignored_minerals = []
        ignore_minerals = input('''- Will you ignore any minerals?
        < y/n > ; Default = n  ___ ''') or 'n'
        if ignore_minerals == 'y':
            for mineral in self.ross.parameters['described_minerals']:
                ignore_mineral = input('''- Will you ignore {}?
                < y/n > ; Default = n  ___ ''') or 'n'
                if ignore_mineral == 'y':
                    ignored_minerals.append(mineral)
                         
        # define existing scale        
        existing_scale = input('''- Scale already exists in the simulated module?
        < True > or < False > ; Default = False  ___ ''') or False
        while existing_scale not in ['True', 'False', False]:
            print('--> ERROR: ')
            existing_scale = input('''- Scale already exists in the simulated module?
            < True > or < False > ; Default = False  ___ ''') or False
        existing_scale = bool(existing_scale)
        existing_parameters = {}
        if existing_scale:
            for mineral in self.ross.variables['described_minerals']:
                if mineral not in ignored_minerals:
                    existing_parameters['saturation'] = input(f'''- What is the saturation index of {mineral}?''')
                    while True:
                        try:
                            existing_parameters['saturation'] = float(existing_parameters['saturation'])
                            break
                        except:            
                            print('''--> ERROR: The saturation must be a number.''')
                            existing_parameters['saturation'] = input(f'''- What is the saturation index of {mineral}?''')
                         
                    existing_parameters['initial_moles'] = input(f'''- What is the existing moles of {mineral}?''')
                    while True:
                        try:
                            existing_parameters['initial_moles'] = float(existing_parameters['initial_moles'])
                            break
                        except:            
                            print('''--> ERROR: The initial_moles must be a number.''')
                            existing_parameters['initial_moles'] = input(f'''- What is the existing moles of {mineral}?''')

        self.ross.equilibrium_phases(block_comment = '', ignored_minerals = ignored_minerals, existing_parameters = existing_parameters)


    def selected_output(self,):
        # announcement
        announcement = '\nParameterize the output content:'
        print(announcement, '\n', '='*len(announcement))      
        
        # define output_filename
        selected_output_file_name = self.ross.selected_output_file_name(output_filename = None)
        output_filename = input(f'''- What is the selected_output_file_name?
        Default = {selected_output_file_name}''') or selected_output_file_name
                         
        self.ross.selected_output(output_filename = output_filename)


    def export(self, external_file = False):
        # announcement
        announcement = '\nParameterize the export:'
        print(announcement, '\n', '='*len(announcement))      
        
        default_simulation_name, default_input_path, default_output_path = self.ross.define_paths()
                         
        # define simulation_name
        simulation_name = input(f'''- What is the simulation_name?
        Default = {default_simulation_name}''') or default_simulation_name          
                         
        # define input_path
        input_path = input(f'''- What is the input_path?
        Default = {default_input_path}''') or default_input_path    
                         
        # define output_path
        output_path = input(f'''- What is the output_path?
        Default = {default_output_path}''') or default_output_path    
                         
        self.ross.export(simulation_name = simulation_name, input_path = input_path, output_path = output_path, external_file = external_file)

    def parse_input(self,):
        # announcement
        announcement = '\nParse the existing input file:'
        print(announcement, '\n', '='*len(announcement))      
        
        # load the input file
        input_file_path = input('What is the input file path?')
        while len(input_file_path) > self.max_argument_length:
            excessive_characters = self.ross.parameters['input_path'][self.max_argument_length:]
            print(f'''--> ERROR: The {excessive_characters} characters of the input file path {input_file_path} are beyond the limit of characters. Provide a valid input file path.''')
            input_file_path = input('What is the input file path?')
        while not os.path.exists(input_file_path):
            print(f'''--> ERROR: The input file path {input_file_path} does not exist. Provide a valid input file path.''')
            input_file_path = input('What is the input file path?')
                  
        # define simulation
        simulations = ['scaling', 'brine']
        simulation = input('''- Will you simulate scaling or brine formation?
        Default = < scaling >  __ ''') or 'scaling'
        while simulation not in simulations:
            print('''--> ERROR: Only scaling or brine formation simulations can be simulated through ROSSpy.''')    
            simulation = input('''- Will you simulate scaling or brine formation?
            Default = < scaling >  __ ''') or 'scaling'
                  
        # define water_selection
        solution_description = input('- What is a 1-2 word description of the simulated water source?')
                  
        # define simulation_name
        simulation_name = input(f'''- What is the simulation_name?
        Default = Date-ROSSpy-water_body-simulation_type-database-scaling/brine-perspective-permeate_approach''') or None      
                  
        # define active_feed_area
        active_m2 = input(f'''- What is the active_m2?
        Default = 37''') or 37
        while True:
            try:
                active_m2 = float(active_m2)
                break
            except:            
                print('''--> ERROR: The active_m2 must be a number.''')
                active_m2 = input(f'''- What is the active_m2?
                Default = 37''') or 37
                      
        self.ross.parse_input(input_file_path = input_file_path, simulation = simulation, water_selection = solution_description, simulation_name = simulation_name, active_m2 = active_m2)
        
    def execute(self,):
        # print the announcement
        announcement = '\nExecute the input file:'
        print(announcement, '\n', '='*len(announcement))

        # execute the PHREEQC batch software
#         bat_path = os.path.join(os.getcwd(), 'phreeqc.bat')
        bat_path = r'phreeqc.bat'
        database_path = r'{}'.format(self.ross.parameters['database_path'])
        input_path = r'{}'.format(self.ross.parameters['input_path'])        
        
        output_path = r'{}'.format(self.ross.parameters['output_path'])
        if len(output_path) > self.max_argument_length:
            output_path = re.sub('(input.pqi)', 'output.pqo', input_path)
            print(f'''\n\n--> ERROR: The output file path was abridged to {output_path} to maintain validity as an argument for the batch PHREEQC software.\n\n''')            
        
        for path in [bat_path, input_path, os.path.dirname(output_path), database_path]:
            if not os.path.exists(path):
                print(f'-> ERROR: The < {path} > path does not exist\n')
        
        proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
        command = str.encode(bat_path + " \"" + input_path + "\" \"" + output_path + "\" \"" + database_path + "\"\n") 
        proc.stdin.write(command)
        proc.stdin.close()  
        proc.wait()
        

#         self.raw_data = self.ross.execute(simulated_to_real_time = 9.29)
        selected_output_path = os.path.join(os.path.dirname(__file__), self.ross.selected_output_file_name)
        self.ross.results['csv_data'] = pandas.read_table(open(selected_output_path), sep='\t')
        for column in self.ross.results['csv_data'].columns:
            new_column = column.strip()
            self.ross.results['csv_data'].rename(columns={column:new_column}, inplace = True)
            if re.search('Unnamed', column):
                del self.ross.results['csv_data'][column]
        self.selected_output_new_path = os.path.join(self.ross.simulation_path, 'selected_output.csv')
        self.ross.results['csv_data'].to_csv(self.selected_output_new_path)

# execute and process the input file
    def process_selected_output(self,):
        # announcement
        announcement = '\nProcess the output:'
        print(announcement, '\n', '='*len(announcement))     
        
        selected_output_path = input(f'''- What is the selected_output_path?
        Default = {self.selected_output_new_path} ____ ''') or self.selected_output_new_path
        while not os.path.exists(selected_output_path):
            selected_output_path = input(f'''- What is the selected_output_path?
            Default = {self.selected_output_new_path} ____ ''') or self.selected_output_new_path
        
        # define plot_title
        plot_title = input('''- What is the title of the plot?
        Default = None ____ ''')
                  
        # define the fonts
        possible_fonts = ['xx-small','x-small','small', 'medium', 'large', 'x-large', 'xx-large']
        for font in possible_fonts:
            print(f'< {font} >')
        title_font = input('''- What is the title font?
        Default = xx-large''') or 'xx-large'
        while title_font not in possible_fonts:
            print(f'''--> ERROR: Only one of the {possible_fonts} fonts can be parameterized.''')    
            title_font = input('''- What is the title font?
            Default = xx-large''') or 'xx-large'
                  
        label_font = input('''- What is the title font?
        Default = x-large''') or 'x-large'
        while label_font not in possible_fonts:
            print(f'''--> ERROR: Only one of the {possible_fonts} fonts can be parameterized.''')    
            label_font = input('''- What is the title font?
            Default = x-large''') or 'x-large'
                  
        # define x_label_number
        x_label_number = input('''- How many x-axis ticks are desired for the plot?
        Default = 6 ____ ''') or 6
        while True:
            try:
                x_label_number = float(x_label_number)
                break
            except:            
                print('''--> ERROR: The x_label_number must be a number.''')
                x_label_number = input('''- How many x-axis ticks are desired for the plot?
                Default = 6 ____ ''') or 6
                
        # define individual_plots
        default_individual_plots = False
        if self.ross.parameters['simulation_perspective'] == 'all_time':
            default_individual_plots = True
        individual_plots = input(f'''- Will each mineral be individually plotted?
        < True > or < False > ; Default = {default_individual_plots} ____ ''') or default_individual_plots
        while individual_plots not in ['True', 'False', True, False]:
            print(f'''--> ERROR: Only < True > or < False > are supported.''')    
            individual_plots = input(f'''- Will each mineral be individually plotted?
            < True > or < False > ; Default = {default_individual_plots} ____ ''') or default_individual_plots
        individual_plots = bool(individual_plots)   
        
        # define export_name
        if individual_plots:
            default_figure_name = 'mineral_names'
        else:
            default_figure_name = 'all_minerals'
        export_name = input(f'''- What is the export name for the figure(s)?
        Default = < {default_figure_name} > ____ ''') or default_figure_name
                  
        # define export_name
        figure_formats = ['svg', 'pdf', 'png', 'jpeg', 'jpg', 'eps']
        for format in figure_formats:
            print(f'< {format} >')
        export_format = input(f'''- What is the export name for the figure(s)?
        Default = svg ____ ''') or 'svg'
        while export_format not in figure_formats:
            print(f'''--> ERROR: Only one of the {figure_formats} formats is supported.''')    
            export_format = input('''- What is the title font?
            Default = svg ___ ''') or 'svg'
            
        self.processed_data = self.ross.process_selected_output(selected_output_path = selected_output_path, plot_title = plot_title, title_font = title_font, label_font = label_font, x_label_number = x_label_number, export_name = export_name, export_format = export_format, individual_plots = individual_plots)
        

def conduct_iROSSpy():
    iross = iROSSpy()
    
    # execute iROSSpy
    create_input = input('''- You will create an input file?
    Default = < True >''') or True
    while create_input not in ['True', 'False', True]:
        print(f'''--> ERROR: Only < True > or < False > are supported.''')    
        create_input = input('''- You will create an input file?
        < True > or < False > ; Default = < True >''') or True
    if create_input == 'False':
        create_input = False
    else:
        create_input = True
    
    # create or import the simulation input file
    if create_input:
        iross.define_general()
        iross.transport()
        iross.reaction()
        iross.solutions()
        iross.equilibrium_phases()
        iross.selected_output()
        iross.export()
    else:
        iross.parse_input()
    
    # execute and process the simulation
    iross.execute()
    iross.process_selected_output()
    

conduct_iROSSpy()