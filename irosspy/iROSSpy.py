# ============= WELCOME MESSAGE =================
from itertools import chain
print('\n\n')
message = ('''* iROSSpy v1: Interactive Reverse Osmosis Scaling Simulation in Python * 
Andrew Philip Freiburger, Ethan Sean Chan; University of Victoria, 2022''')

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

# ============= CODE =================
from scipy.constants import nano, milli, minute
import subprocess
import rosspy
import pandas
import time, json, os, re

class iROSSpy():
    def __init__(self):
        self.working_directory = os.getcwd()
        
        # announcements
        announcement = 'Parameterize initial details:'
        print(announcement, '\n', '='*len(announcement))
        
        # execute either a predefined input file or define an input file
        self.create_input = input('''- You will create an input file?
        Default = < True >''') or True
        while self.create_input not in ['True', 'False', True]:
            print(f'''--> ERROR: Only < True > or < False > are supported.''')    
            self.create_input = input('''- You will create an input file?
            < True > or < False > ; Default = < True >''') or True

        if self.create_input == 'False':
            self.create_input = False
        else:
            self.create_input = True
            
        # define operating system
        operating_systems = ['windows', 'unix']
        self.operating_system =  input('''- Is your computer operating system windows or unix?
        < {} > or < {} > ; Default = < windows >  __ '''.format(operating_systems[0], operating_systems[1])) or 'windows'
        while self.operating_system not in operating_systems:
            print('''--> ERROR: The operating system is not one of the accepted options.''')
            self.operating_system =  input('''- Is your computer operating system windows or unix?
            < {} > or < {} >  __ '''.format(operating_systems[0], operating_systems[1]))  or 'windows'

        # define verbosity
        self.verbose =  input('''- Will your simulation be verbose?
        < True > or < False > ; Default = < False >  __ ''') or False
        while self.verbose not in ['True', 'False', False]:
            print('''--> ERROR: The verbosity must be an True or False.''')    
            self.verbose =  input('''- Your simulation be verbose?
            < True > or < False > ; Default = < False >  __ ''') or False
        if self.verbose == 'True':
            self.verbose = True
                             
        # define the maximal argument path length 
        self.max_argument_length = len(r"C:\Users\Andrew Freiburger\Dropbox\My PC (DESKTOP-M302P50)\Documents\UVic Civil Engineering\PHREEQC\ROSS\examples\scaling\scale_validation\2021-12-18-ROSSpy-red_sea-transport-pitzer-scaling-all_distance-LinPerm\2021-12-28-ROSSpy--transport-pitzer-scaling-all_")
                
        self.ross = rosspy.ROSSPkg('pitzer', operating_system = self.operating_system)

    def define_general(self,):
        # announcement
        announcement = '\nParameterize general simulation conditions:'
        print(announcement, '\n', '='*len(announcement))   
                             
        # define database
        for database in self.ross.databases:
            print(f'< {database} >') 
        database_selection = input('''- What database do you select?
        Default = < pitzer >  __ ''') or 'pitzer'
        while database_selection not in self.ross.databases:
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
#        domains = ['single', 'dual']
#        domain = input('''- Will you simulate the single or dual domain?
#        Default = < single >  __ ''') or 'single'
#        while domain not in domains:
#            print('''--> ERROR: Only the single or dual domains can be simulated through ROSSpy.''')    
#            domain = input('''- Will you simulate the single or dual domain?
#            Default = < single >  __ ''') or 'single'
#        
        # define domain_phase
#        if domain == 'dual': 
#            domain_phases = ['mobile', 'immobile']
#            domain_phase = input('''- Will you evaluate the mobile or immobile phase?
#            Default = < immobile >  __ ''') or 'immobile'
#            while domain_phase not in domain_phases:
#                print('''--> ERROR: Only the mobile or immobile phases can be simulated through ROSSpy.''')    
#                domain_phase = input('''- Will you simulate the single or dual domain?
#                Default = < single >  __ ''') or 'immobile'
#        else:
#            domain_phase = None
                             
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

        # initiate the ROSSpy object
        self.ross = rosspy.ROSSPkg(database_selection, simulation, simulation_type, self.operating_system, quantity_of_modules = quantity_of_modules, simulation_title = simulation_title, verbose = self.verbose)


    def reactive_transport(self,):
        # announcement
        announcement = '\nParameterize reactive transport conditions:'
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
        custom_module = input('''- You will customize the module characteristics?
        < True > or < False > ; Default = < False >  __ ''') or False
        while custom_module not in ['True', 'False', False]:
            print('''--> ERROR: Only True or False are accepted responses.''')    
            custom_module = input('''- You will customize the module characteristics?
            < True > or < False > ; Default = < False >  __ ''') or False
        if custom_module== 'True':
            custom_module = True
        
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
                        module_characteristics[characteristic] = float(module_characteristics[characteristic])
                        break
                    except:            
                        print('''--> ERROR: The provided value is not a number.''')    
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
#        exchange_factor = ''
#        if self.ross.parameters['domain'] == 'dual':
#            exchange_factor = input('''- What is the (1/sec) exchange rate between the mobile and immobile phases of the simulated dual domain?
#            Default = < 1e10 >  __ ''') or 1e10     
#            while True:
#                try:
#                    exchange_factor = int(exchange_factor)
#                    break
#                except:            
#                    print('''--> ERROR: The exchange_factor must be an integer.''')    
#                    exchange_factor = input('''- What is the (1/sec) exchange rate between the mobile and immobile phases of the simulated dual domain?
#                    Default = < 1e10 >  __ ''') or 1e10                                                        
                                    
        # define permeate_approach
        permeate_approaches = ['linear_permeate', 'linear_cf']
        print('\n')
        for approach in permeate_approaches:
            print(f'< {approach} >')
        permeate_approach = input('''- Which permeate_approach will be simulated?
        Default = < linear_permeate >  __ ''') or 'linear_permeate'
        while permeate_approach not in permeate_approaches:
            print(f'''--> ERROR: Only {permeate_approaches} are accepted parameters.''')    
            permeate_approach = input('''- Which permeate_approach will be simulated?
            Default = < linear_permeate >  __ ''') or 'linear_permeate'

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
                                                  
        # define permeate_efficiency
        permeate_efficiency = input('''- What is the 0 < x < 1 efficiency of the simulated module?
        Default = < 1 >  __ ''') or 1 
        while True:
            try:
                permeate_efficiency = float(permeate_efficiency)
                break
            except:            
                print('''--> ERROR: The cells_per_module must be an integer.''')    
                permeate_efficiency = input('''- What is the 0 < x < 1 efficiency of the simulated module?
                Default = < 1 >  __ ''') or 1 
                         
        # define head_loss
        head_loss = input('''- What is the  0 < x < 1 relative pressure of the effluent to the feed?
        Default = < 0.89 >  __ ''') or 0.89
        while True:
            try:
                head_loss = float(head_loss)
                break
            except:            
                print('''--> ERROR: The head_loss must be a number.''')    
                head_loss = input('''- What is the  0 < x < 1 relative pressure of the effluent to the feed?
                Default = < 0.89 >  __ ''') or 0.89
                         
        # apply the parameters to the ROSSpy instance                 
        self.ross.reactive_transport(simulation_time, simulation_perspective, final_cf, module_characteristics, permeate_efficiency, head_loss, cells_per_module = cells_per_module)


    def feed_geochemistry(self,):
        # announcement
        announcement = '\nParameterize feed_geochemistry:'
        print(announcement, '\n', '='*len(announcement))      
        
        # define water_selection
        options = self.ross.feed_sources
        options.append('custom')
        for water in options:
            print(f'< {water} >')   
        water_selection = input('''- Which water body will you simulate?
                                Type < custom > to define a feed water. __ ''')
        while water_selection not in options:
            print('''--> ERROR: One of the printed water_bodies may be parameterized.''')    
            water_selection = input('''- Which water body will you simulate?
                                    Type < custom > to define a feed water. __ ''')
        
        solution_description = None
        water_characteristics = {}
        self.ross.parameters['solution_elements'] = []
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
                if water_characteristics[element]['concentration (ppm)'] != 0:
                    self.ross.parameters['solution_elements'].append(element)
                                         
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
        < True > or < False > ; Default = False ___ ''') or False
        if ph_charge_balance == 'True':
            ph_charge_balance = True
        water_characteristics['Alkalinity'] = {}
        if ph_charge_balance is True:
            water_characteristics['Alkalinity']['value'] = water_characteristics['Alkalinity']['reference'] = None
            parameterized_ph_charge = True
        else:
            water_characteristics['Alkalinity']['value'] = input('''- What is the alkalinity in ((eq of CaCO3)/(Kg of water))?''')
            water_characteristics['Alkalinity']['reference'] = input('''- What is the reference for the alkalinity?''')
            parameterized_ph_charge = False
            while True:
                try:
                    water_characteristics['Alkalinity']['value'] = float(water_characteristics['Alkalinity']['value'])
                    break
                except:            
                    print('''--> ERROR: The Alkalinity must be a number.''')
                    water_characteristics['Alkalinity']['value'] = input('''- What is the alkalinity in ((eq of CaCO3)/(Kg of water))?  ___''')

        # define the solution elements and possible minerals
        if water_selection != 'custom':
            path = os.path.join(self.ross.parameters['root_path'], 'water_bodies', f'{water_selection}.json')
            feed_dict= json.load(open(path))
            self.ross._define_elements(feed_dict['element'])
        self.ross._described_minerals()                         
        
        # define ignored minerals
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
        if existing_scale == 'True':
            existing_scale = True
        
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

        # dexecute the ROSSpy function
        self.ross.feed_geochemistry(water_selection, water_characteristics, solution_description, ignored_minerals, existing_parameters, parameterized_ph_charge)
        
        
    def parse_input(self,):
        # announcement
        announcement = '\nParse the existing input file:'
        print(announcement, '\n', '='*len(announcement))      
        
        # load the input file
        input_file_path = input('What is the input file path?')
#        while len(input_file_path) > self.max_argument_length:
#            excessive_characters = self.ross.parameters['input_path'][self.max_argument_length:]
#            print(f'''--> ERROR: The {excessive_characters} characters of the input file path {input_file_path} are beyond the limit of characters. Provide a valid input file path.''')
#            input_file_path = input('What is the input file path?')
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
        water_selection = input('- What is a 1-2 word description of the simulated water source?')
                                   
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
                      
        # initiate the ROSSpy object
        self.ross = rosspy.ROSSPkg('pitzer', simulation, )
        self.ross.parse_input(input_file_path, water_selection, active_m2)
        
        
    def execute(self,):
        # announcement
        announcement = '\nExecute the simulation:'
        print(announcement, '\n', '='*len(announcement))      
        
#        # define output_filename
#        selected_output_file_name = self.ross.selected_output_file_name(output_filename = None)
#        output_filename = input(f'''- What is the selected_output_file_name?
#        Default = {selected_output_file_name}''') or selected_output_file_name  
        
        default_simulation_path, default_simulation_name = self.ross._define_paths()
                         
        # define simulation_name
        simulation_name = input(f'''- What is the simulation_name?
        Default = {default_simulation_name}''') or default_simulation_name
                         
        # define output_path
        output_path = input(f'''- What is the simulation directory?
        Default = {default_simulation_path}''') or default_simulation_path       
        
        # define plot_title
        plot_title = input('''- What is the title of the plot?
        The default is a description of the simulation parameters ____ ''')
                  
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
                  
        label_font = input('''- What is the label font?
        Default = x-large''') or 'x-large'
        while label_font not in possible_fonts:
            print(f'''--> ERROR: Only one of the {possible_fonts} fonts can be parameterized.''')    
            label_font = input('''- What is the label font?
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
        if individual_plots == 'False':
            individual_plots = False
        elif individual_plots == 'True':
            individual_plots = True 
        
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
            
        # define the selected_output section
        self.ross._selected_output(simulation_name)
        
        # execute the PHREEQC batch software
        bat_path = os.path.join(os.path.dirname(__file__), 'phreeqc.bat')
        #bat_path = r'./phreeqc.bat'
        database_path = r'{}'.format(self.ross.parameters['database_path'])
        input_path = r'{}'.format(self.ross.parameters['input_path'])        
        output_path = re.sub('(input.pqi)', 'output.pqo', input_path)
#        if len(output_path) > self.max_argument_length:
#            output_path = re.sub('(input.pqi)', 'output.pqo', input_path)
#            print(f'''\n\n--> ERROR: The output file path was abridged to {output_path} to maintain validity as an argument for the batch PHREEQC software.\n\n''')            
        
        for path in [bat_path, input_path, os.path.dirname(output_path), database_path]:
            unicode = ''.join([ch if ord(ch) < 128 else '\t\t' for ch in path])
            if len(path) != len(unicode):
                print(f'-> ERROR:', unicode)
            if not os.path.exists(path):
                print(f'-> ERROR: The < {path} > path does not exist\n')
        
        proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
        command = str.encode(bat_path + " \"" + input_path + "\" \"" + output_path + "\" \"" + database_path + "\"\n") 
        proc.stdin.write(command)
        proc.stdin.close()  
        proc.wait()   

        # process the simulation results
        directory = os.path.dirname(__file__)
        so_filename = os.path.basename(self.ross.parameters['selected_output_file_name'])
        selected_output_path = os.path.join(directory, so_filename)
        self.ross.selected_output = pandas.read_table(open(selected_output_path), sep='\t')
        for column in self.selected_output.columns:
            new_column = column.strip()
            self.ross.selected_output.rename(columns={column:new_column}, inplace = True)
            if re.search('Unnamed', column):
                del self.ross.selected_output[column]
        self.ross.selected_output.to_csv(os.path.join(self.ross.simulation_path, 'selected_output.csv'))
        
        # visualize the simulation results, while bypassing the API execution operations
        self.ross.execute(simulation_name, selected_output_path, output_path, plot_title, title_font, label_font, x_label_number, export_name, export_format, individual_plots)
        

def conduct_iROSSpy():
    iross = iROSSpy()
    
    # create or import the simulation input file
    if iross.create_input:
        iross.define_general()
        iross.reactive_transport()
        iross.feed_geochemistry()
    else:
        iross.parse_input()
    
    # execute and process the simulation
    iross.execute()
    
    # closing input
    options = ['another', 'close']
    continue_simulation = input('''Will you perform another simulation, or will you close iROSSpy?
    < another > or < close > ''')
    if continue_simulation == 'close':
        pass
    elif continue_simulation == 'another':
        conduct_iROSSpy()
    while continue_simulation not in options:
        print(f'--> ERROR: The input is not one of the accepted {options} options. Provide an accepted option.')
        continue_simulation = input('''Will you perform another simulation, or will you close iROSSpy?
    < another > or < close > ''')
            

conduct_iROSSpy()              