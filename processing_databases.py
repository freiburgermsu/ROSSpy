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
                        formula = reactants.split(' + ')[0].strip()
                        name = database.at[index, 'content']
                        name = re.sub('\s*\d*', '', name)
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
        
def mineral_masses(minerals):
    def parse_stoich(formula, ch_number):
        ch_no = ch_number
        stoich = ''
        while re.search('[0-9\.]', formula[ch_no]):                    
            stoich += formula[ch_no]
            if ch_no == len(formula)-1:
                break
            ch_no += 1
        skips = ch_no - ch_number
        return skips, float(stoich)

    def parse_mineral_formula(formula, ch_number, final):
        stoich = 0
        skips = 0
        if re.search('[A-Z]',formula[ch_number]):
            if final:
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
                if re.search('[0-9]', formula[ch_number+2]):
                    skips, stoich = parse_stoich(formula, ch_number+2) # float(formula[ch_number+2])
                elif re.search('[A-Z\(]', formula[ch_number+2]):
                    stoich = 1
                else:
                    print('--> ERROR: The mineral formula {} is unpredictable.'.format(formula))

                print('elemental_mass1', elemental_masses[element])
                print('stoich1', stoich)
                mass = stoich * elemental_masses[element] 
                skip_characters = skips + 1
                print('mass1', mass)
                return skip_characters, mass

            elif re.search('[0-9]', formula[ch_number+1]):
                skips, stoich = parse_stoich(formula, ch_number+1)

                element = formula[ch_number]
                print('\n', element)
                stoich = float(stoich)

                print('elemental_mass2', elemental_masses[element])
                print('stoich2', stoich)
                mass = stoich * elemental_masses[element] 
                skip_characters = skips  + 1
                print('mass2', mass)
                return skip_characters, mass

            elif re.search('[A-Z\(]', formula[ch_number+1]):
                element = formula[ch_number]
                print('\n', element)
                stoich = 1

                print('elemental_mass3', elemental_masses[element])
                print('stoich3', stoich)
                mass = stoich * elemental_masses[element] 
                skip_characters = 0
                print('mass3', mass)
                return skip_characters, mass

        if re.search('[0-9]',formula[ch_number]):
            stoich = formula[ch_number]
            subbed = re.sub('H|2|O', '', formula[ch_number+1:ch_number+4])
            if subbed == '':
                skip_characters = 3
                water_mass = elemental_masses['H'] * 2 + elemental_masses['O']
                mass = float(stoich) * water_mass
                return skip_characters, mass
            else:
                print(f'--> ERROR: The {formula} formula is not predictable.')

    # evaluate all of the described minerals
    for mineral in minerals:
        formula = minerals[mineral]['formula']
        print(f'\n\n{mineral} {formula}\n', '='*len(f'{mineral} {formula}'), '\n')
        skip_characters = 0
        minerals[mineral]['mass'] = 0
        mass = 0
        for ch_number in range(len(formula)):
            print('mass', minerals[mineral]['mass'])
            print('ch_number', ch_number)
            print('character', formula[ch_number])
            if skip_characters > 0:
                skip_characters -= 1
                continue

            final = False
            if ch_number == len(formula)-1:
                final = True

            if formula[ch_number] == ':':
                continue
            elif formula[ch_number] == '(':
                group_mass = 0
                ch_no2 = ch_number + 1
                while formula[ch_no2] != ')':
                    if skip_characters > 0:
                        skip_characters -= 1
                        ch_no2 += 1
                        continue

                    skips, mass = parse_mineral_formula(formula, ch_no2, final = final)
                    group_mass += mass
                    skip_characters += skips
                    ch_no2 += 1

                ch_no2 += 1
                print('group_mass', group_mass)
                skips, stoich = parse_stoich(formula, ch_no2)
                print('group_stoich', stoich)
                minerals[mineral]['mass'] += group_mass * stoich
                skip_characters = ch_no2 - ch_number + skips - 1

                print(skip_characters)
            else:
                skip_characters, mass = parse_mineral_formula(formula, ch_number, final = final)
                minerals[mineral]['mass'] += mass

        print('{} mass: {}'.format(mineral, minerals[mineral]['mass']))
    
    return minerals
        
    
def database_json_creation(database, elements, minerals):
    database_json = {'elements': {}, 'minerals': {}}
    
    # create the elements JSON
    for index, element in elements.iterrows():
        database_json['elements'][element['elements']] = {'charge_specie': element['alk'], 'gfw_formula':element['gfw_formula'], 'element_gfw':element['element_gfw']}
    
    # create the mienrals JSON
    for index, mineral in minerals.iterrows():
        mineral['formula'] = re.sub('Cyanide', 'CN', mineral['formula'])
        database_json['minerals'][mineral['phases']] = {'formula': mineral['formula'], 'mass': ''}
       
    database_json['minerals'] = mineral_masses(database_json['minerals'])
        
    # export the JSON files
    database_json_name = re.sub('.dat$', '.json', database)
    with open(database_json_name, 'w') as output:
        json.dump(database_json, output, indent = 4)