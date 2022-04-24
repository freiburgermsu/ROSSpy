from scipy.constants import nano, milli
import sigfig 
import rosspy
import pandas
import shutil, os, re

def isnumber_type(obj):
    if type(obj) is float or type(obj) is int:
        return True
    
    
# standard test conditions
simulation_time = 100
module_characteristics = {
    'module_diameter_mm': {
        'value':160
        },
    'permeate_tube_diameter_mm':{
            'value':9
            },  
    'module_length_m':{
            'value':1
            },  
    'permeate_flow_m3_per_hour':{
            'value':45
            },   
    'max_feed_flow_m3_per_hour':{
            'value':10
            },
    'membrane_thickness_mm':{
            'value':300 * (nano / milli)
            },  
    'feed_thickness_mm':{
            'value':0.98
            },
    'active_m2':{
            'value':34
            },
    'permeate_thickness_mm':{
            'value':0.31
            },
    'polysulfonic_layer_thickness_mm':{
            'value':0.05
            },     
    'support_layer_thickness_mm':{
            'value':0.15
            }
}   
    
def test_init():
    ross = rosspy.ROSSPkg('pitzer', verbose = False)

    # ensure that the dictionaries are created
    for dic in [ross.parameters, ross.variables, ross.results, ross.results['figures']]:
        assert type(dic) is dict
    
    # ensure that the booleans are created
    for boo in [ross.verbose, ross.jupyter,ross.printing]:
        assert type(boo) is bool
        
    # ensure that the booleans are created
    for lis in [ross.databases, ross.feed_sources,ross.results['general_conditions']]:
        assert type(lis) is list
        
    assert int(sigfig.round(ross.water_gL, 2)) == 1.0e3
    assert ross.water_mw == 18.0153
    assert isnumber_type(ross.parameters['quantity_of_modules'])
    
    for param in ['os', 'simulation_type', 'simulation', 'simulation_type', 'domain', 'database_selection']:
        assert type(ross.parameters[param]) is str

    for path in ['database_path','root_path']:
        assert os.path.exists(ross.parameters[path])

    for dic in [ross.elements, ross.minerals,]:
        assert type(dic) is dict

def test_reactive_transport():
    ross = rosspy.ROSSPkg('pitzer', verbose = False)
    ross.reactive_transport(simulation_time, module_characteristics = module_characteristics)

    # affirm qualities of the simulation
    for param in ['module_diameter_mm', 'permeate_tube_diameter_mm', 'module_length_m', 'permeate_flow_m3_per_hour', 'max_feed_flow_m3_per_hour', 'membrane_thickness_mm', 'feed_thickness_mm', 'active_m2', 'permeate_thickness_mm', 'polysulfonic_layer_thickness_mm', 'support_layer_thickness_mm']:
        assert ross.ro_module[param] == module_characteristics[param]   

    for param in ['repeated_membrane_winding_mm', 'simulation_time', 'exchange_factor', 'cells_per_module', 'active_m2_cell', 'timestep', 'permeate_moles_per_cell']:
        assert isnumber_type(ross.parameters[param])

    for var in ['cell_meters', 'feed_cubic_meters', 'feed_kg', 'feed_moles', 'Reynold\'s number', 'feed_kg_cell', 'feed_moles_cell']:
        assert isnumber_type(ross.variables[var])

    assert type(ross.parameters['simulation_perspective']) is str
    assert type(ross.results['transport_block']) is list
    for line in ross.results['reaction_block']:
        assert type(line) is str   
        if line not in ['', '\n', None, ' ', '\n\tH2O -1; \d+', '#[A-Za-z\s]+']:
            assert re.search('REACTION|\#',line)

def test_feed_geochemistry():
    water_selection = 'michigan_basin'
    ignored_minerals = ['gypsum']
    
    ross = rosspy.ROSSPkg('pitzer', verbose = False)    
    ross.reactive_transport(simulation_time)
    ross.feed_geochemistry(water_selection, ignored_minerals = ignored_minerals)

    # affirm qualities of the simulation
    for line in ross.results['solution_block']:
        assert type(line) is str   

    assert type(ross.parameters['solution_elements']) is list
    for param in ['water_selection']:
        assert type(ross.parameters[param]) is str
    
    # affirm qualities of the simulation
    for var in ['cell_meters', 'feed_cubic_meters', 'feed_kg', 'feed_moles', 'Reynold\'s number', 'feed_kg_cell', 'feed_moles_cell']:
        assert isnumber_type(ross.variables[var])
        
    assert type(ross.variables['described_minerals']) is dict

    for line in ross.results['equilibrium_phases_block']:
        assert type(line) is str    

def test_execute():
    water_selection = 'michigan_basin'
    
    ross = rosspy.ROSSPkg('pitzer', verbose = False)    
    ross.reactive_transport(simulation_time)
    ross.feed_geochemistry(water_selection)
    ross.execute()

    # affirm qualities of the simulation
    for block in [ross.results['complete_lines'], ross.results['selected_output_block']]:
        for line in block:
            assert type(line) is str   

    # affirm qualities of the simulation
    for path in [ross.parameters['simulation_path'], ross.simulation_path, os.path.join(ross.simulation_path, 'parameters.csv'), os.path.join(ross.simulation_path, 'variables.csv'), os.path.join(ross.simulation_path, 'effluent_predictions.csv')]: 
        assert os.path.exists(path)
        
    for file in ['parameters.csv', 'variables.csv', 'effluent_predictions.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
        
    for df in [ross.selected_output, ross.processed_data]:
        assert type(df) is pandas.core.frame.DataFrame
    assert type(ross.variables['run_time (s)']) is float
    shutil.rmtree(ross.simulation_path)   
    
def test_parse_input():
    simulation = 'scaling'
    water_selection = 'palo_duro_basin'
    input_file_path = os.path.join(os.path.dirname(__file__), '2020-09-03_APF_Palo Duro Basin-BW30-400_PE=100%_1.1.pqi')
    
    ross = rosspy.ROSSPkg('pitzer', simulation, verbose = False)
    ross.parse_input(input_file_path, water_selection)
    ross.execute()
    
    # affirm the execution of the simulation
    for file in ['all_minerals.svg', 'parameters.csv', 'scaling_data.csv', 'selected_output.csv', 'variables.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
        
    shutil.rmtree(ross.simulation_path)      

def test_all_distance_brine(): #!!! This displays an erroneous plot and x-axis domain
    water_selection = 'michigan_basin'
    simulation = 'brine'
    simulation_perspective = 'all_distance'
    
    ross = rosspy.ROSSPkg('pitzer', simulation, verbose = False)    
    ross.reactive_transport(simulation_time, simulation_perspective)
    ross.feed_geochemistry(water_selection)
    ross.execute()

    # affirm qualities of the simulation
    for var in ['initial_solution_mass', 'final_solution_mass','simulation_cf','final_time']:
        print(var)
        assert isnumber_type(ross.variables[var])

    for file in ['brine.svg', 'input.pqi', 'parameters.csv', 'brine_concentrations.csv', 'selected_output.csv', 'variables.csv', 'effluent_predictions.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
    
    shutil.rmtree(ross.simulation_path)   
    
def test_all_time_brine():
    water_selection = 'michigan_basin'
    simulation = 'brine'

    ross = rosspy.ROSSPkg('pitzer', simulation, verbose = False)    
    ross.reactive_transport(simulation_time)
    ross.feed_geochemistry(water_selection)
    ross.execute()

    # affirm qualities of the simulation
    for var in ['initial_solution_mass', 'final_solution_mass','simulation_cf','final_time']:
        assert isnumber_type(ross.variables[var])
        
    for file in ['brine.svg', 'input.pqi', 'parameters.csv', 'brine_concentrations.csv', 'selected_output.csv', 'variables.csv', 'effluent_predictions.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
        
    shutil.rmtree(ross.simulation_path) 

def test_all_distance_scaling():
    water_selection = 'michigan_basin'
    
    ross = rosspy.ROSSPkg('pitzer', verbose = False)    
    ross.reactive_transport(simulation_time)
    ross.feed_geochemistry(water_selection)
    ross.execute()

    # affirm qualities of the simulation
    for var in ['precipitated_minerals']:
        assert ross.variables[var]

    for file in ['all_minerals.svg', 'input.pqi', 'parameters.csv', 'scaling_data.csv', 'selected_output.csv', 'variables.csv', 'effluent_predictions.csv', 'scale_ions.json']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
        
    shutil.rmtree(ross.simulation_path) 

def test_all_time_scaling():
    water_selection = 'michigan_basin'
    simulation_perspective = 'all_time'
    
    ross = rosspy.ROSSPkg('pitzer', verbose = False)    
    ross.reactive_transport(simulation_time, simulation_perspective)
    ross.feed_geochemistry(water_selection)
    ross.execute()

    # affirm qualities of the simulation
    for var in ['initial_solution_mass', 'final_solution_mass', 'simulation_cf', 'final_time']:
        assert ross.variables[var]

    for file in ['all_minerals.svg', 'input.pqi', 'parameters.csv', 'scaling_data.csv', 'selected_output.csv', 'variables.csv', 'effluent_predictions.csv', 'scale_ions.json']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
        
    shutil.rmtree(ross.simulation_path) 
    
def test_test():    
    ross = rosspy.ROSSPkg('pitzer', verbose = False)
    ross.test()
    
    # affirm the execution of the simulation
    for file in ['all_minerals.svg', 'input.pqi', 'parameters.csv', 'scaling_data.csv', 'selected_output.csv', 'variables.csv', 'effluent_predictions.csv', 'scale_ions.json']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
        
    shutil.rmtree(ross.simulation_path)   