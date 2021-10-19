from chemicals import periodic_table
from glob import glob
from numpy import nan
import pandas
import json
import re


# define elemental masses
elemental_masses = {}
for element in periodic_table:
    elemental_masses[element.symbol] = element.MW

def database_parsing(db):
    print(f'\n\n\n{db}')
    
    database = pandas.read_table(db, sep='\n')
    start_master = False
    if database.columns == ['SOLUTION_MASTER_SPECIES']:
        start_master = True
    database.columns = ['content']

    elements_rows = []
    minerals_rows = []
    elemental_parsing = False
    mineral_parsing = False
    for index, row in database.iterrows():
        if (re.search('SOLUTION_MASTER_SPECIES', row['content']) or start_master) and not elemental_parsing :
            while not re.search('SOLUTION_SPECIES', database.at[index, 'content']):
                split_row = database.at[index, 'content'].split()
                if all(not re.search('^#', entity) for entity in split_row):
                    elements_rows.append(split_row)
                index+= 1
            elemental_parsing = True

        if re.search('PHASES', row['content']) and not mineral_parsing:
            loop = False
            while not re.search('PITZER|EXCHANGE_MASTER_SPECIES|SURFACE_MASTER_SPECIES', database.at[index, 'content']):
                if not loop:
                    minerals_rows.append(['phases', 'formula'])
                    loop = True
                    index += 1
                    continue

                if re.search('(^\w+\s*\d*$)',database.at[index, 'content']):
                    reactants = database.at[index+1, 'content'].split(' = ')[0]
                    if all('#' not in entity for entity in reactants):
                        formula = reactants.split('+')[0].strip()
                        name = database.at[index, 'content']
                        name = re.sub('(\s+\d*)', '', name)
                        minerals_rows.append([name, formula])
                index+= 1
                
                if index == len(database):
                    break
            mineral_parsing = True

    # define the elements content for the database
    elements = pandas.DataFrame(elements_rows)
    elements.fillna(' ')
#     elements.columns = elements.iloc[0]
    elements.drop([0], inplace = True)
    for column in elements:
        nan_entries = 0
        alphanumeric_entries = 0
        for entry in elements[column]:
#             if entry in [None, ' ', nan]:
            if entry is not None:
                if re.search('[a-z]|[0-9]', entry, re.IGNORECASE) and entry is not None:
    #             if entry is str or entry is float or entry is int:
                    alphanumeric_entries += 1
                else:
                    nan_entries += 1
            else:
                nan_entries += 1
        if nan_entries > alphanumeric_entries and len(elements.columns) > 5:
            print('deleted column: ', column)
            del elements[column]
    
    print(elements)
    elements.columns = ['elements', 'species', 'alk', 'gfw_formula', 'element_gfw']
#     elements.rename(columns = {'SOLUTION_MASTER_SPECIES':'elements'}, inplace = True)
    print(elements)
#         elements = elements.iloc[pandas.RangeIndex(len(elements)).drop([x for x in range(4)])]
    elements_list = list(elements['elements'])

    # define the minerals content for the database
    minerals = pandas.DataFrame(minerals_rows)
    minerals.columns = minerals.iloc[0]
    minerals = minerals.drop(0)
    print(minerals)
    mineral_list = list(minerals['phases'])
    formula_list = list(minerals['formula'])
    
    return elements, minerals

def parse_stoich(formula, ch_number):
    ch_no = ch_number
    stoich = ''
    skip_minus = False
    if formula[ch_no] == '.':
        stoich = '0.'
        ch_no += 1
        skip_minus = True
    while re.search('[0-9\.]', formula[ch_no]):                    
        stoich += formula[ch_no]
        if ch_no == len(formula)-1:
            break
#             if not re.search('[0-9\.]', formula[ch_no+1]):
#                 break
        ch_no += 1

    skips = ch_no - ch_number
    if skip_minus:
        skips -= 1
    stoich = float(stoich)
    if re.search('(\.0$)', str(stoich)):
        stoich = int(stoich)
    return skips, stoich

def group_parsing(formula, ch_number, final):
    group_mass = 0
    group_mass2 = 0
    ch_no2 = ch_number + 1
    skip_characters = 0
    double = False
    final_mass2 = 0
    group = True
    while formula[ch_no2] != ')' and ch_no2 != len(formula)-1:
        print('ch_no2', ch_no2)
        print('character', formula[ch_no2])
        if skip_characters > 0:
            skip_characters -= 1
            ch_no2 += 1
            continue
        if re.search('[\(\s]', formula[ch_no2]):
            ch_no2 += 1
            double = True
            continue
            
        if double:
            skips, mass = parse_mineral_formula(formula, ch_no2, final, group = group)
            group_mass2 += mass
            skip_characters += skips
            ch_no2 += 1
    #                     print(ch_no2)
            if formula[ch_no2] == ')':
                print('group_mass2', group_mass2)
                skips, stoich = parse_stoich(formula, ch_no2 + 1)
                print('group_stoich2', stoich)
                final_mass2 = group_mass2 * stoich
                skip_characters += skips 
                ch_no2 += 1
                double = False
        
        else:
            print(formula[ch_no2])
            skips, mass = parse_mineral_formula(formula, ch_no2, final, group = group)
            group_mass += mass
            skip_characters += skips
            ch_no2 += 1
#                     print(ch_no2)

    print('group_mass', group_mass)
    if ch_no2 == len(formula)-1:
        stoich = 1
        print('group_stoich1', stoich)
    elif re.search('[A-Z:\(]', formula[ch_no2+1]):
        stoich = 1
        print('group_stoich2', stoich)
    else:
        skips, stoich = parse_stoich(formula, ch_no2 + 1)
        print('group_stoich3', stoich)
        
    final_mass = group_mass * stoich + final_mass2
    skip_characters = ch_no2 - ch_number + skips 
    print('skips', skips)
    print('skip_characters_group', skip_characters)
#     if double:
#         skip_characters -= 1
    return skip_characters, final_mass

def parse_mineral_formula(formula, ch_number, final = False, group = False):
    stoich = 0
    skips = 0
    if re.search('[ +)]',formula[ch_number]):
        skip_characters = 0
        mass = 0
        return skip_characters, mass
    elif re.search('[A-Z]',formula[ch_number]):
        if final or len(formula) == 1:
            element = formula[ch_number]
            print('\n', element)
            stoich = 1

            print('elemental_mass', elemental_masses[element])
            print('stoich', stoich)
            mass = stoich * elemental_masses[element] 
            skip_characters = 0
            print('mass', mass)
            return skip_characters, mass

        elif re.search('[a-z]', formula[ch_number+1]):
            element = formula[ch_number] + formula[ch_number+1]
            print('\n', element)            
            if ch_number+1 != len(formula)-1:
                if re.search('[0-9]', formula[ch_number+2]):
                    skips, stoich = parse_stoich(formula, ch_number+2) # float(formula[ch_number+2])
#                     if group:
#                         skips = 0
                elif re.search('[A-Z\(:]', formula[ch_number+2]):
                    stoich = 1
                elif formula[ch_number+2] == '.':
                    skips, stoich = parse_stoich(formula, ch_number+2) # float(formula[ch_number+2])
                    skips += len('.')
                elif formula[ch_number+2] == ' ':
                    stoich = 1
                    skips = 1
                else:
                    print('--> ERROR: The mineral formula {} is unpredictable.'.format(formula))
            else:
                stoich = 1
                
            print('elemental_mass1', elemental_masses[element])
            print('stoich1', stoich)
            mass = stoich * elemental_masses[element] 
            skip_characters = skips + 1
            print('mass1', mass)
            
            print('skips1', skip_characters)

            return skip_characters, mass

        elif re.search('[0-9]', formula[ch_number+1]):
            skips, stoich = parse_stoich(formula, ch_number+1)

            element = formula[ch_number]
            print('\n', element)
            print('elemental_mass2', elemental_masses[element])
            print('stoich2', stoich)
            mass = stoich * elemental_masses[element] 
            skip_characters = skips
            print('mass2', mass)
            return skip_characters, mass

        elif formula[ch_number+1] == '.':
            skips, stoich = parse_stoich(formula, ch_number+1)

            element = formula[ch_number]
            print('\n', element)
            print('elemental_mass2', elemental_masses[element])
            print('stoich2', stoich)
            mass = stoich * elemental_masses[element] 
            skip_characters = skips+1
            print('mass2', mass)
            return skip_characters, mass

        elif re.search('[A-Z():+ ]', formula[ch_number+1]):
            element = formula[ch_number]
            print('\n', element)
            stoich = 1

            print('elemental_mass3', elemental_masses[element])
            print('stoich3', stoich)
            mass = stoich * elemental_masses[element] 
            skip_characters = 0
            print('mass3', mass)
#                 print('skips', skip_characters)
            return skip_characters, mass
            

    elif re.search(':',formula[ch_number]):
        print('element', formula[ch_number])
        skips = stoich = space = back_space = 0
        if re.search('[0-9]', formula[ch_number+1]):
            skips, stoich = parse_stoich(formula, ch_number+1)
        if re.search('[( ]',formula[ch_number+1+skips]):
            space = 1
            back_space = 1
        if formula[ch_number+1+skips+space:ch_number+4+skips+space] == 'H2O':
            print('subbed portion', formula[ch_number+1+skips+space:ch_number+4+skips+space])
            skip_characters = len('H2O')+skips+back_space
            water_mass = elemental_masses['H'] * 2 + elemental_masses['O']
            print('water_mass', water_mass)
            print('water_stoich', stoich)
            mass = float(stoich) * water_mass
            return skip_characters, mass
        elif formula[ch_number+1+skips+space:ch_number+4+skips+space] == 'H\+':
            print('subbed portion', formula[ch_number+1+skips+space:ch_number+4+skips+space])
            skip_characters = len('H+')+skips
            proton_mass = elemental_masses['H']
            mass = float(stoich) * proton_mass
            return skip_characters, mass
        elif re.search('[A-Z]',formula[ch_number+skips+1]):
            skip_characers, mass = group_parsing(formula, ch_number+skips, final)
            group_mass = mass * stoich
            return skip_characers, group_mass
        else:
            print(f'--> ERROR: The {formula} formula is not predictable.')
            return 0, 0
      
    elif formula[ch_number-1] == '.':
        return 0, 0
      
    elif re.search('[0-9]', formula[ch_number]):
        skips, stoich = parse_stoich(formula, ch_number)
        return skips, stoich
    
        
def mineral_masses(db, minerals):
    # evaluate all of the described minerals
    for mineral in minerals:
        formula = minerals[mineral]['formula']
        print(f'\n\n{db} {mineral} {formula}\n', '='*len(f'{mineral} {formula}'), '\n')
        skip_characters = 0
        minerals[mineral]['mass'] = 0
        mass = 0
        first = True
        double = triple = False
        
        for ch_number in range(len(formula)):
            print('total_mass', minerals[mineral]['mass'])
            print('ch_number', ch_number)
            print('character', formula[ch_number])
            if skip_characters > 0:
                print('skip_characters', skip_characters)
                skip_characters -= 1
                continue

            final = False
            if ch_number == len(formula)-1:
                print('final')
                final = True
            if formula[ch_number] == '(':
                skip_characters, mass = group_parsing(formula, ch_number, final)
                minerals[mineral]['mass'] += mass 
                if triple:
                    skip_characters += 1 
                    if mineral == 'Berthierine_ISGS':
                        skip_characters += 3
#                     if mineral == 'Glauconite':
#                         skip_characters += 1
                    
                if double:
                    if mineral in ['Boltwoodite', 'Corkite']:
                        skip_characters -= 1
                    if mineral == 'Glauconite':
                        skip_characters -= 5
                    if mineral in ['Saponite_SapCa', 'Vermiculite_SO']:
                        skip_characters -= 9
                    print('second')
                    triple = True
                    double = False
                    
                if first:
                    if mineral in ['Brochantite', 'Borax', 'Antlerite', 'Corkite', 'Kasolite', 'Phosgenite', 'Tsumebite', 'Artinite', 'Jaffeite'] or (mineral in ['Burkeite', 'Dawsonite'] and re.search('sit', db)):
                        skip_characters -= 1
                    if mineral == 'Berthierine_ISGS':
                        skip_characters -= 10
                    if mineral == 'Glauconite':
                        skip_characters -= 5
                    if mineral in ['Saponite_SapCa', 'SmectiteMX80', 'Vermiculite_SO']:
                        skip_characters -= 9
                    first = False
                    double = True
                    print('first')
                elif mineral == 'Berthierine_ISGS':
                    skip_characters -= 4
                        
            else:
                print('element', formula[ch_number])
                skip_characters, mass = parse_mineral_formula(formula, ch_number, final)
                minerals[mineral]['mass'] += mass

#         if mineral in ['Hydroxyapatite']:
#             minerals[mineral]['mass'] += elemental_masses['O']
                
        print('{} mass: {}'.format(mineral, minerals[mineral]['mass']))
    
    return minerals
        
    
def database_json_creation(database, elements, minerals):
    database_json = {'elements': {}, 'minerals': {}}
    
    # create the elements JSON
    for index, element in elements.iterrows():
        database_json['elements'][element['elements']] = {'charge_specie': element['alk'], 'gfw_formula':element['gfw_formula'], 'element_gfw':element['element_gfw']}
    
    # create the mienrals JSON
    for index, mineral in minerals.iterrows():
        mineral['formula'] = re.sub('Cyanide|Cyanate', 'CN', mineral['formula'])
        if re.search('PHASES', mineral['phases']):
            continue
        database_json['minerals'][mineral['phases']] = {'formula': mineral['formula'], 'mass': ''}
       
    database_json['minerals'] = mineral_masses(database, database_json['minerals'])
        
    # export the JSON files
    database_json_name = re.sub('.dat$', '.json', database)
    with open(database_json_name, 'w') as output:
        json.dump(database_json, output, indent = 4)