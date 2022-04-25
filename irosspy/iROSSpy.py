# ============= WELCOME MESSAGE =================
from itertools import chain
print('\n\n')
message = ('''* iROSSpy I: interactive Reverse Osmosis Scaling Simulation in Python, I *
Andrew Philip Freiburger, 2022''')
lines = message.split('\n')
width = max(map(len, lines))

upper = f'╔{"═" * (width+2)}╗\n'
lower = f'╚{"═" * (width+2)}╝' 
middles = []
for line in lines:
    space = " " * int((width-len(line)+2)/2)
    middles.append(f'║{space}{line}{space}║\n')

print(''.join(chain(upper, middles, lower)), '\n')

print('The program is loading. It may take a minute...\n\n')

# ============= CODE =================
from scipy.constants import minute
from warnings import warn
#import subprocess
import rosspy
#import pandas
#import json, re
import os 

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
        
max_argument_length = len(r"C:\Users\Andrew Freiburger\Dropbox\My PC (DESKTOP-M302P50)\Documents\UVic Civil Engineering\PHREEQC\ROSS\examples\scaling\scale_validation\2021-12-18-ROSSpy-red_sea-transport-pitzer-scaling-all_distance-LinPerm\2021-12-28-ROSSpy--transport-pitzer-scaling-all_")

class iROSSpy():
    def __init__(self):
        # announcements
        announcement = '\nParameterize initial details:'
        print(announcement, '\n', '='*len(announcement))
            
        # define operating system
        operating_systems = ['windows', 'unix']
        my_binder =  input(f'''- You are operating in MyBinder? (boolean)
        Default = < True >  __ ''') or True
        if my_binder in [True, 'True']:
            operating_system = 'unix'
        else:            
            operating_system =  input(f'''- Is your computer operating system windows or unix?
            < {operating_systems[0]} > or < {operating_systems[1]} > ; Default = < windows >  __ ''') or 'windows'
            while operating_system not in operating_systems:
                warn('''The operating system is not one of the accepted options.''')
                operating_system =  input(f'''- Is your computer operating system windows or unix?
                < {operating_systems[0]} > or < {operating_systems[1]} > ; Default = < windows >  __ ''') or 'windows'
            
        # define database
        self.ross = rosspy.ROSSPkg('pitzer', operating_system = operating_system)
        for database in self.ross.databases:
            print(f'< {database} >') 
        database_selection = input('''- What database do you select?
        Default = < pitzer >  __ ''') or 'pitzer'
        while database_selection not in self.ross.databases:
            warn('''The database is not current defined in this interface.
            Only < pitzer > and < phreeqc > are defined for Windows and < pitzer > is defined for Macintosh. __ ''')    
            database_selection = input('- What database do you select?')
        
        # define simulation
        simulations = ['scaling', 'brine']
        simulation = input('''- Will you simulate < scaling > or < brine >?
        Default = < scaling >  __ ''') or 'scaling'
        while simulation not in simulations:
            warn('''Only < scaling > or < brine > generation can be simulated through ROSSpy.''')    
            simulation = input('''- Will you simulate < scaling > or < brine > formation?
            Default = < scaling >  __ ''') or 'scaling'
            
        # define simulation_perspective
        simulation_types = ['transport','evaporation']
        simulation_type = input('''- Will you simulate evaporation or transport?
        Default = < transport >  __ ''') or 'transport'
        while simulation_type not in simulation_types:
            print('''--> ERROR: Only transport or evaporation can be simulated through ROSSpy.''')  
            simulation_type = input('''- Will you simulate evaporation or transport?
            Default = < transport >  __ ''') or 'transport'   

        # define verbosity
        export_content =  input('''- Will you export the simulation content? (boolean)
        Default = < True >  __ ''') or True
        while export_content not in ['False', 'True', True]:
            warn('''The verbosity must be < True > or < False >.''')    
            export_content =  input('''- Will your simulation be verbose? (boolean)
            Default = < True >  __ ''') or True
        if export_content == 'False':
            export_content = False
                          
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
        while not isnumber(quantity_of_modules):
            print('''--> ERROR: The quantity of modules must be an integer.''')    
            quantity_of_modules = input('''- How many in-series RO modules will you simulate?
            Default = < 1 >  __ ''') or 1
        
        # define simulation title
        simulation_title = input('- What is the title of your simulation?   __ ')   

        # define verbosity
        verbose =  input('''- Will your simulation be verbose? (boolean)
        Default = < False >  __ ''') or False
        while verbose not in ['True', 'False', False]:
            warn('''The verbosity must be < True > or < False >.''')    
            verbose =  input('''- Will your simulation be verbose? (boolean)
            Default = < False >  __ ''') or False
        if verbose == 'True':
            verbose = True

        # initiate the ROSSpy object
        self.ross = rosspy.ROSSPkg(database_selection, simulation, simulation_type, operating_system, export_content, None, quantity_of_modules, simulation_title, verbose, printing = True, jupyter = True)


    def reactive_transport(self,):
        # announcement
        announcement = '\nParameterize reactive transport conditions:'
        print(announcement, '\n', '='*len(announcement))                                  
                             
        # define simulation_time
        simulation_time = input('''- How many minutes will you simulate?
        Default = < 3 > minutes __ ''') or 3
        while not isnumber(simulation_time):
            warn('''The simulation_time must be numerical.''')    
            simulation_time = input('''- How many minutes will you simulate?
            Default = < 3 > minutes __ ''') or 3        
        simulation_time = float(simulation_time) * minute
        
        # define simulation_perspective
        simulation_perspectives = ['all_time', 'all_distance']
        if self.ross.parameters['simulation'] == 'scaling':
            default = 'all_distance'
        if self.ross.parameters['simulation'] == 'brine':
            default = 'all_time'
        simulation_perspective = input(f'''- Will you simulate all_time or all_distance?
        Default = < {default} >  __ ''') or default
        while simulation_perspective not in simulation_perspectives:
            warn(f'''Only {[x for x in simulation_perspectives]} can be simulated.''')
            simulation_perspective = input(f'''- Will you simulate all_time or all_distance?
            Default = < {default} >  __ ''') or default
            
        # define module_characteristics
        ro_module = input('''What RO module will you simulate?
        Default = BW30-400 ''') or 'BW30-400'
        while ro_module not in self.ross.ro_modules:
            warn(f'The parameterized ro_module {ro_module} is not defined in the ro_module.json parameter file ({self.ross.ro_modules.keys()}))')
            ro_module = input('''What RO module will you simulate?
            Default = BW30-400 ''') or 'BW30-400'
        
        custom_module = input('''- You will customize the module characteristics? (boolean)
        Default = < False >  __ ''') or False
        while custom_module not in ['True', 'False', False]:
            warn('''Only < True > or < False > are accepted.''')    
            custom_module = input('''- You will customize the module characteristics? (boolean)
            Default = < False >  __ ''') or False
        if custom_module== 'True':
            custom_module = True
            
        module_characteristics = {}
        if custom_module:
            for characteristic in self.ross.ro_module:
                module_characteristics[characteristic] = {}
                module_characteristics[characteristic]['value'] = input(f''' What is the value of the {characteristic} module characteric?
                default = {self.ross.ro_module[characteristic]['value']}  ___ ''') or self.ross.ro_module[characteristic]['value']
                while not isnumber(module_characteristics[characteristic]['value']):
                    warn('''The provided value is not a number.''')    
                    module_characteristics[characteristic]['value'] = input(f''' What is the value of the {characteristic} module characteric?
                    default = {self.ross.ro_module[characteristic]['value']}  ___ ''') or self.ross.ro_module[characteristic]['value']
                         
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
            warn(f'''Only {permeate_approaches} are accepted parameters.''')    
            permeate_approach = input('''- Which permeate_approach will be simulated?
            Default = < linear_permeate >  __ ''') or 'linear_permeate'

        # define final_cf
        final_cf = None
        if permeate_approach == 'linear_cf':
            final_cf = input('''- What is the effluent CF (the effluent concentration relative to the feed concentration)?
            Default = < 2 >  __ ''') or 2
            while not isnumber(final_cf):
                warn('''The final_cf must be numerical.''')    
                final_cf = input('''- What is the effluent concentration relative to the feed concentration?
                Default = < 1.3 >  __ ''') or 1.3
                                                                      
        # define head_loss
        head_loss = input('''- What is the relative decline in pressure over the module? (0 < x < 1)
        Default = < 0.1 >  __ ''') or 0.1
        while not isnumber(head_loss):
            warn('''--> ERROR: The head_loss must be a number.''')    
            head_loss = input('''- What is the relative decline in pressure over the module? (0 < x < 1)
            Default = < 0.1 >  __ ''') or 0.1
                        
        # define cells_per_module
        cells_per_module = input('''- Into how many cells is your module discretized?
        Default = < 12 > __ ''') or 12     
        while not isnumber(cells_per_module):
            warn('''The cells_per_module must be numerical.''')    
            cells_per_module = input('''- Into how many cells is your simulated module discretized?
            Default = < 12 >  __ ''') or 12  
                        
        # define the timestep resolution
        coarse_timestep = input('''Will your simulation use a < coarse > or < fine > timestep?
        Default = < coarse >''') or 'coarse'
        while coarse_timestep not in ['fine', 'coarse']:
            warn('''The verbosity must be < coarse > or < fine >.''')    
            coarse_timestep = input('''Will your simulation use a < coarse > or < fine > timestep?
            Default = < coarse >''') or 'coarse'
                                     
        # apply the parameters to the ROSSpy instance                 
        self.ross.reactive_transport(simulation_time, simulation_perspective, final_cf, module_characteristics, ro_module, 1, head_loss, cells_per_module = cells_per_module, coarse_timestep = coarse_timestep)


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
                                Type < custom > to define feed geochemistry. __ ''')
        while water_selection not in options:
            warn('''Select one of the printed water_bodies.''')    
            water_selection = input('''- Which water body will you simulate?
                                    Type < custom > to define feed geochemistry. __ ''')
        
        solution_description, parameterized_ph_charge = None, True
        self.ross.parameters['solution_elements'], water_characteristics = [], {}
        if water_selection == 'custom':
            solution_description = input('- What is the description of the custom water source?')
            for element in self.ross.elements:
                if element == 'Alkalinity':
                    continue
                water_characteristics[element] = {}
                water_characteristics[element]['concentration (ppm)'] = input(f'''- What is ppm concentration of {element} in the feed?''')
                while not isnumber(water_characteristics[element]['concentration (ppm)']):
                    warn('''Elemental concentrations must be numeric.''')    
                    water_characteristics[element]['concentration (ppm)'] = input(f'''- What is ppm concentration of {element} in the feed?''') 
                water_characteristics[element]['reference'] = input(f'- What is the reference for {element} = {water_characteristics[element]["concentration (ppm)"]} ppm?')
                self.ross.parameters['solution_elements'].append(element)
                                         
            # define pH
            water_characteristics['pH'] = {}
            water_characteristics['pH']['value'] = input('''- What is the pH of the customized solution? ___ ''')
            while not isnumber(water_characteristics['pH']['value']):
                warn('''The pH must be numeric.''')    
                water_characteristics['pH']['value'] = input('''- What is the pH of the customized solution? ___ ''')
            water_characteristics['pH']['reference'] = input('- What is the reference for this pH?')

            # define temperature
            water_characteristics['temperature'] = {}
            water_characteristics['temperature']['value'] = input('''- What is the temperature of the customized solution? ___ ''')
            while not isnumber(water_characteristics['temperature']['value']):
                warn('''The temperature must be numberic.''')
                water_characteristics['temperature']['value'] = input('''- What is the temperature of the customized solution? ___ ''')
            water_characteristics['temperature']['reference'] = input('- What is the reference for this temperature?')

            # define pe
            water_characteristics['pe'] = {}
            water_characteristics['pe']['value'] = input('''- What is the pe of the customized solution? ___ ''')
            while not isnumber(water_characteristics['pe']['value']):
                print('''--> ERROR: The pe must be a number.''')
                water_characteristics['pe']['value'] = input('''- What is the pe of the customized solution? ___ ''')
            water_characteristics['pe']['reference'] = input('- What is the reference for this pe?')
                         
            # define pH charge balance
            parameterized_ph_charge = input('''- charge balanced pH and defined alkalinity are exclusive. You will charge balance the pH? (boolean)
            Default = True ___ ''') or True
            if parameterized_ph_charge == 'False':
                parameterized_ph_charge = False
            water_characteristics['Alkalinity'] = {}
            if parameterized_ph_charge is True:
                water_characteristics['Alkalinity']['value'] = water_characteristics['Alkalinity']['reference'] = None
            else:
                water_characteristics['Alkalinity']['value'] = input('''- What is the alkalinity in (eq of CaCO3/kg of water)?''')
                water_characteristics['Alkalinity']['reference'] = input(f'''- What is the reference for the alkalinity = {water_characteristics['Alkalinity']['value']} (eq of CaCO3/Kg of water)?''')
                while not isnumber(water_characteristics['Alkalinity']['value']):
                    warn('The Alkalinity must be numberic.')
                    water_characteristics['Alkalinity']['value'] = input('''- What is the alkalinity in ((eq of CaCO3)/(Kg of water))?  ___''')

        # define ignored minerals
        self.ross._described_minerals()   
        ignored_minerals = []
        ignore_minerals = input('''- You will ignore a mineral(s)? (boolean)
        Default = False  ___ ''') or False
        while ignore_minerals not in ['True', 'False', False]:
            warn('''The ignore_minerals must be < True > or < False >.''')    
            ignore_minerals = input('''- You will ignore a mineral(s)? (boolean)
            Default = False  ___ ''') or False
        if ignore_minerals == 'True':
            ignore_minerals = True
        if ignore_minerals is True:
            for mineral in self.ross.parameters['described_minerals']:
                ignore_mineral = input(f'''- You will ignore {mineral}? (boolean)
                Default = False  ___ ''') or False
                if ignore_mineral == 'True':
                    ignore_mineral = True
                if ignore_mineral is True:
                    ignored_minerals.append(mineral)
                         
        # define existing scale        
        existing = input('''- Scale already exists in the simulated module? (boolean)
        Default = False  ___ ''') or False
        while existing not in ['True', 'False', False]:
            warn('''Only < True > and < False > are accepted.''')    
            existing = input('''- Scale already exists in the simulated module? (boolean)
            Default = False  ___ ''') or False
        if existing == 'True':
            existing = True
        
        existing_scale = {}
        if existing is True:
            for mineral in self.ross.variables['described_minerals']:
                if mineral not in ignored_minerals:
                    existing_scale['saturation'] = input(f'''- What is the saturation index of {mineral}?''')
                    while not isnumber(existing_scale['saturation']):
                        warn('''The saturation must be numberic.''')
                        existing_scale['saturation'] = input(f'''- What is the saturation index of {mineral}?''')
                         
                    existing_scale['initial_moles'] = input(f'''- What is the existing moles of {mineral}?''')
                    while not isnumber(existing_scale['initial_moles']):
                        warn('''The saturation must be numberic.''')
                        existing_scale['initial_moles'] = input(f'''- What is the existing moles of {mineral}?''')

        # dexecute the ROSSpy function
        self.ross.feed_geochemistry(water_selection, water_characteristics, solution_description, ignored_minerals, existing_scale, parameterized_ph_charge)
               
        
def execute(ross):
    # announcement
    announcement = '\nExecute the simulation:'
    print(announcement, '\n', '='*len(announcement))      
    
    # define simulation name and path
    default_simulation_path, default_simulation_name = ross._define_paths()
    simulation_name = input(f'''- What is the simulation_name?
    Default = {default_simulation_name}''') or default_simulation_name
    simulation_directory = input(f'''- What is the simulation directory?
    Default = {default_simulation_path}''') or default_simulation_path   
    
    define_paths = True
    if simulation_name == default_simulation_name and simulation_directory == default_simulation_path:
        define_paths = False
    
    # define plot_title
    figure_title = input('''- What is the title of the figure?
    The default is a description of the simulation parameters ____ ''')
              
    # define the fonts
    possible_fonts = ['xx-small','x-small','small', 'medium', 'large', 'x-large', 'xx-large']
    for font in possible_fonts:
        print(f'< {font} >')
    title_font = input('''- What is the title font?
    Default = xx-large''') or 'xx-large'
    while title_font not in possible_fonts:
        warn(f'''Only one of the {possible_fonts} fonts can be parameterized.''')    
        title_font = input('''- What is the title font?
        Default = xx-large''') or 'xx-large'
              
    label_font = input('''- What is the label font?
    Default = x-large''') or 'x-large'
    while label_font not in possible_fonts:
        warn(f'''Only one of the {possible_fonts} fonts can be parameterized.''')    
        label_font = input('''- What is the label font?
        Default = x-large''') or 'x-large'
              
    # define x_label_number
    x_label_number = input('''- How many x-axis ticks are desired for the plot?
    Default = 6 ____ ''') or 6
    while not isnumber(x_label_number):
        warn('''The x_label_number must be numberic.''')
        x_label_number = input('''- How many x-axis ticks are desired for the plot?
        Default = 6 ____ ''') or 6
                          
    # define export_name
    figure_formats = ['svg', 'pdf', 'png', 'jpeg', 'jpg', 'eps']
    for format in figure_formats:
        print(f'< {format} >')
    export_format = input(f'''- What is the export name for the figure(s)?
    Default = svg ____ ''') or 'svg'
    while export_format not in figure_formats:
        warn(f'''Only one of the {figure_formats} formats is supported.''')    
        export_format = input('''- What is the title font?
        Default = svg ___ ''') or 'svg'
        
    scale_ions = input('''- You will determine the ionic proportion of the scale? (boolean)
    Default = True  ___ ''') or True
    while scale_ions not in ['False', 'True', True]:
        warn('''Only < True > and < False > are accepted.''')    
        scale_ions = input('''- You will determine the ionic proportion of the scale? (boolean)
        Default = True  ___ ''') or True
    if scale_ions == 'False':
        scale_ions = False
    
    # define the selected_output section
#    ross._selected_output(simulation_name)
    
    # execute the PHREEQC batch software
#        bat_path = os.path.join(os.path.dirname(__file__), 'phreeqc.bat')
##        bat_path = 'phreeqc.bat'
#        database_path = r'{}'.format(self.ross.parameters['database_path'])
#        input_path = r'{}'.format(self.ross.parameters['input_path'])        
#        output_path = re.sub('(input.pqi)', 'output.pqo', input_path)
##        if len(output_path) > max_argument_length:
##            output_path = re.sub('(input.pqi)', 'output.pqo', input_path)
##            print(f'''\n\n--> ERROR: The output file path was abridged to {output_path} to maintain validity as an argument for the batch PHREEQC software.\n\n''')            
#        
#        for path in [bat_path, input_path, os.path.dirname(output_path), database_path]:
#            unicode = ''.join([ch if ord(ch) < 128 else '\t\t' for ch in path])
#            if len(path) != len(unicode):
#                print(f'-> ERROR:', unicode)
#            if not os.path.exists(path):
#                raise FileNotFoundError(f'The < {path} > path does not exist\n')
#        
#        proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
#        command = str.encode(bat_path + " \"" + input_path + "\" \"" + output_path + "\" \"" + database_path + "\"\n") 
#        proc.stdin.write(command)
#        proc.stdin.close()  
#        proc.wait()   

    # process the simulation results
#    selected_output_path = os.path.join(os.path.dirname(__file__), os.path.basename(ross.parameters['selected_output_file_name']))
#    ross.selected_output = pandas.read_table(open(selected_output_path), sep='\t')
#    for column in ross.selected_output.columns:
#        new_column = column.strip()
#        ross.selected_output.rename(columns={column:new_column}, inplace = True)
#        if re.search('Unnamed', column):
#            del ross.selected_output[column]
#    ross.selected_output.to_csv(os.path.join(ross.simulation_path, 'selected_output.csv'))
    
    # visualize the simulation results, while bypassing the API execution operations
    ross.execute(simulation_name, None, simulation_directory, figure_title, title_font, label_font, x_label_number, None, export_format, scale_ions, define_paths)
        
def conduct_iROSSpy():
    # execute either a predefined input file or define an input file
    create_input = input('''- You will create an input file? (boolean)
    Default = < True >''') or True
    while create_input not in ['True', 'False', True]:
        warn(f'''--> ERROR: Only < True > or < False > are supported.''')    
        create_input = input('''- You will create an input file?
        < True > or < False > ; Default = < True >''') or True
    if create_input == 'False':
        create_input = False
    else:
        create_input = True
            
    # create or import the simulation input file
    if create_input is True:
        iross = iROSSpy()
        iross.reactive_transport()
        iross.feed_geochemistry()
        ross = iross.ross
    else:
        # announcement
        announcement = '\nParse the existing input file:'
        print(announcement, '\n', '='*len(announcement))      
        
        # load the input file
        input_file_path = input('What is the input file path?')
#        while len(input_file_path) > max_argument_length:
#            excessive_characters = self.ross.parameters['input_path'][max_argument_length:]
#            print(f'''--> ERROR: The {excessive_characters} characters of the input file path {input_file_path} are beyond the limit of characters. Provide a valid input file path.''')
#            input_file_path = input('What is the input file path?')
        while not os.path.exists(input_file_path):
            warn(f'''The input file path {input_file_path} does not exist.''')
            input_file_path = input('What is the input file path?')
                  
        # define simulation
        simulation = input('''- Will you simulate < scaling > or < brine > formation?
        Default = < scaling >  __ ''') or 'scaling'
        while simulation not in ['scaling', 'brine']:
            warn('''Only < scaling > or < brine > can be simulated through ROSSpy.''')    
            simulation = input('''- Will you simulate < scaling > or < brine > formation?
            Default = < scaling >  __ ''') or 'scaling'
                  
        # define water_selection
        water_selection = input('- What is a 1-2 word description of the simulated water source?')
                                   
        # define active_feed_area
        active_m2 = input(f'''- What is the active_m2?
        Default = 37''') or 37
        while not isnumber(active_m2):
            warn('''The active_m2 must be numberic.''')
            active_m2 = input(f'''- What is the active_m2?
            Default = 37''') or 37
                      
        # initiate the ROSSpy object
        ross = rosspy.ROSSPkg('pitzer', simulation)
        ross.parse_input(input_file_path, water_selection, active_m2)
    
    # execute and process the simulation
    execute(ross)
    
    # closing input
#    options = ['another', 'close']
#    continue_simulation = input('''Will you perform another simulation, or will you close iROSSpy? (boolean)''')
#    if continue_simulation in [False, 'False']:
#        pass
#    elif continue_simulation in [True, 'True']:
#        conduct_iROSSpy()
#    while continue_simulation not in [True, 'True', False, 'False']:
#        warn(f'The input {continue_simulation} is not one of the accepted {options} options.')
#        continue_simulation = input('''You will perform another simulation?
#    < another > or < close > ''')
            

conduct_iROSSpy()              