from scipy.constants import nano, milli
from to_precision import auto_notation
from glob import glob
import rosspy
import pandas
import os
import re

def isnumber(obj):
    if type(obj) is float or type(obj) is int:
        return True
    
def test_init():
    ross = rosspy.ROSSPkg(verbose = False)

    assert type(ross.parameters) is dict
    assert type(ross.variables) is dict
    assert type(ross.results) is dict
    assert type(ross.verbose) is bool

def test_define_general():
    ross = rosspy.ROSSPkg(verbose = False)

    # execute the general_conditions function of ROSSpy
    database_selection = 'pitzer'
    ross.define_general(database_selection)

    # affirm qualities of the simulation
    assert int(auto_notation(ross.parameters['water_mw'], 2)) == 18
    assert int(auto_notation(ross.parameters['water_grams_per_liter'], 2)) == 1.0e3
    for param in ['os', 'simulation_type', 'simulation', 'domain', 'root_path', 'database_selection']:
        assert type(ross.parameters[param]) is str
        
    assert type(ross.parameters['quantity_of_modules']) is int

    for path in ['root_path']:
        assert os.path.exists(ross.parameters[path])

    assert type(ross.elements) is dict
    assert type(ross.minerals) is dict

def test_transport():
    ross = rosspy.ROSSPkg(verbose = False)

    # execute the transport function of ROSSpy
    simulation_time = 100
    module_characteristics = {
        'module_diameter_mm':160,
        'permeate_tube_diameter_mm':9,  
        'module_length_m':1,  
        'permeate_flow_m3_per_day':45,   
        'max_feed_flow_m3_per_hour':10,
        'membrane_thickness_mm':300 * (nano / milli),  
        'feed_thickness_mm':0.98,
        'active_m2':34,
        'permeate_thickness_mm':0.31,
        'polysulfonic_layer_thickness_mm':0.05,     
        'support_layer_thickness_mm':0.15
    }   
    database_selection = 'pitzer'
    ross.define_general(database_selection)
    ross.transport(simulation_time, module_characteristics = module_characteristics)

    # affirm qualities of the simulation
    for param in ['module_diameter_mm', 'permeate_tube_diameter_mm', 'module_length_m', 'permeate_flow_m3_per_day', 'max_feed_flow_m3_per_hour', 'membrane_thickness_mm', 'feed_thickness_mm', 'active_m2', 'permeate_thickness_mm', 'polysulfonic_layer_thickness_mm', 'support_layer_thickness_mm']:
        assert ross.parameters[param] == module_characteristics[param]   

    for param in ['simulation_time', 'exchange_factor', 'module_diameter_mm', 'permeate_tube_diameter_mm', 'module_length_m', 'permeate_flow_m3_per_day', 'max_feed_flow_m3_per_hour', 'membrane_thickness_mm', 'feed_thickness_mm', 'active_m2', 'permeate_thickness_mm', 'polysulfonic_layer_thickness_mm', 'support_layer_thickness_mm', 'repeated_membrane_winding_mm', 'cells_per_module', 'active_m2_cell', 'timestep', 'permeate_moles_per_cell']:
        assert isnumber(ross.parameters[param])

    for var in ['cell_meters', 'feed_cubic_meters', 'feed_kg', 'feed_moles', 'Reynold\'s number', 'feed_kg_cell', 'feed_moles_cell']:
        assert isnumber(ross.variables[var])

    for param in ['domain', 'simulation_perspective']:
        assert type(ross.parameters[param]) is str

def test_reaction():
    ross = rosspy.ROSSPkg(verbose = False)
    database_selection = 'pitzer'
    simulation_time = 100
    
    ross.define_general(database_selection)
    ross.transport(simulation_time)
    ross.reaction()

    # affirm qualities of the simulation
    for line in ross.results['reaction_block']:
        assert type(line) is str   
        if line not in ['', '\n', None, ' ']:
            print(line)
            assert re.search('REACTION|\#',line) 

def test_solutions():
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    simulation_time = 100
    
    ross = rosspy.ROSSPkg(verbose = False)    
    ross.define_general(database_selection)
    ross.transport(simulation_time)
    ross.reaction()
    ross.solutions(water_selection)

    # affirm qualities of the simulation
    for line in ross.results['solution_block']:
        assert type(line) is str   

    assert type(ross.parameters['solution_elements']) is list
    for param in ['water_selection']:
        assert type(ross.parameters[param]) is str

def test_equilibrium_phases():
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    ignored_minerals = ['gypsum']
    simulation_time = 100
    
    ross = rosspy.ROSSPkg(verbose = False)
    ross.define_general(database_selection)
    ross.transport(simulation_time)
    ross.reaction()
    ross.solutions(water_selection) 
    ross.equilibrium_phases(None, ignored_minerals)

    # affirm qualities of the simulation
    for var in ['cell_meters', 'feed_cubic_meters', 'feed_kg', 'feed_moles', 'Reynold\'s number', 'feed_kg_cell', 'feed_moles_cell']:
        assert isnumber(ross.variables[var])

    for line in ross.results['equilibrium_phases_block']:
        assert type(line) is str   

def test_selected_output():
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    simulation_time = 100
    
    ross = rosspy.ROSSPkg(verbose = False)
    ross.define_general(database_selection)
    ross.transport(simulation_time)
    ross.reaction()
    ross.solutions(water_selection) 
    ross.equilibrium_phases(None)
    ross.selected_output()

    # affirm qualities of the simulation
    for line in ross.results['selected_output_block']:
        assert type(line) is str   

def test_export():
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    simulation_time = 100
    
    ross = rosspy.ROSSPkg(verbose = False)
    ross.define_general(database_selection)
    ross.transport(simulation_time)
    ross.reaction()
    ross.solutions(water_selection) 
    ross.equilibrium_phases(None)
    ross.selected_output()
    ross.export()

    # affirm qualities of the simulation
    for path in [ross.parameters['input_path'], ross.simulation_path]: 
        assert os.path.exists(path)

    assert type(ross.parameters['input_file_name']) is str

    for line in ross.results['complete_lines']:
        assert type(line) is str
        
    for file in ['parameters.csv', 'variables.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
    
def test_parse_input():
    simulation = 'scaling'
    water_selection = 'michigan_basin'
    input_file_path = '2020-09-03_APF_Palo Duro Basin-BW30-400_PE=100%_1.1.phr'
    
    ross = rosspy.ROSSPkg(verbose = False)
    ross.parse_input(input_file_path, simulation, water_selection)
    ross.execute()
    ross.process_selected_output()
    
    # affirm the execution of the simulation
    for file in ['all_minerals.svg', 'parameters.csv', 'scaling_data.csv', 'selected_output.pqo', 'variables.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))

def test_execute():
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    simulation_time = 100
    
    ross = rosspy.ROSSPkg(verbose = False)
    ross.define_general(database_selection)
    ross.transport(simulation_time)
    ross.reaction()
    ross.solutions(water_selection) 
    ross.equilibrium_phases(None)
    ross.selected_output()
    ross.export()
    ross.execute()

    # affirm qualities of the simulation
    assert os.path.exists(ross.parameters['output_path'])
    assert os.path.exists(os.path.join(ross.simulation_path, 'parameters.csv'))
    assert os.path.exists(os.path.join(ross.simulation_path, 'variables.csv'))
    assert type(ross.results['csv_data']) is pandas.core.frame.DataFrame
    assert type(ross.variables['run_time (s)']) is float

def test_process_selected_output_all_distance_brine():
    ross = rosspy.ROSSPkg(verbose = False)
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    simulation_time = 100
    simulation = 'brine'
    simulation_perspective = 'all_distance'
    
    ross.define_general(database_selection, simulation)
    ross.transport(simulation_time, simulation_perspective)
    ross.reaction()
    ross.solutions(water_selection) 
    ross.equilibrium_phases(None)
    ross.selected_output()
    ross.export()
    ross.execute()
    ross.process_selected_output()

    # affirm qualities of the simulation
    for var in ['initial_solution_mass', 'final_solution_mass','simulation_cf','final_time']:
        ross.variables[var]

    for file in ['brine.svg', 'input.pqi', 'parameters.csv', 'brine_concentrations.csv', 'selected_output.pqo', 'variables.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))

def test_process_selected_output_all_time_brine():
    ross = rosspy.ROSSPkg(verbose = False)
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    simulation_time = 100
    simulation = 'brine'

    ross.define_general(database_selection, simulation)
    ross.transport(simulation_time)
    
    ross.reaction()
    ross.solutions(water_selection) 
    ross.equilibrium_phases(None)
    ross.selected_output()
    ross.export()
    ross.execute()
    ross.process_selected_output()

    # affirm qualities of the simulation
    for var in ['initial_solution_mass', 'final_solution_mass','simulation_cf','final_time']:
        ross.variables[var]
        
    for file in ['brine.svg', 'input.pqi', 'parameters.csv', 'brine_concentrations.csv', 'selected_output.pqo', 'variables.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))

def test_process_selected_output_all_distance_scaling():
    ross = rosspy.ROSSPkg(verbose = False)
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    simulation_time = 100
    
    ross.define_general(database_selection)
    ross.transport(simulation_time)
    ross.reaction()
    ross.solutions(water_selection) 
    ross.equilibrium_phases(None)
    ross.selected_output()
    ross.export()
    ross.execute()
    ross.process_selected_output()

    # affirm qualities of the simulation
    for var in ['precipitated_minerals']:
        ross.variables[var]

    for file in ['all_minerals.svg', 'input.pqi', 'parameters.csv', 'scaling_data.csv', 'selected_output.pqo', 'variables.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))

def test_process_selected_output_all_time_scaling():
    ross = rosspy.ROSSPkg(verbose = False)
    database_selection = 'pitzer'
    water_selection = 'michigan_basin'
    simulation_time = 100
    simulation_perspective = 'all_time'
    
    ross.define_general(database_selection)
    ross.transport(simulation_time, simulation_perspective)
    ross.reaction()
    ross.solutions(water_selection) 
    ross.equilibrium_phases(None)
    ross.selected_output()
    ross.export()
    ross.execute()
    ross.process_selected_output()

    # affirm qualities of the simulation
    for var in ['initial_solution_mass', 'final_solution_mass', 'simulation_cf', 'final_time']:
        ross.variables[var]

    for file in ['Barite.svg', 'Halite.svg', 'input.pqi', 'parameters.csv', 'scaling_data.csv', 'selected_output.pqo', 'variables.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))
    
def test_test():    
    ross = rosspy.ROSSPkg(verbose = False)
    ross.test()
    
    # affirm the execution of the simulation
    for file in ['all_minerals.svg', 'input.pqi', 'parameters.csv', 'scaling_data.csv', 'selected_output.pqo', 'variables.csv']:
        assert os.path.exists(os.path.join(ross.simulation_path, file))