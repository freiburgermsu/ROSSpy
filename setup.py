# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst', encoding='utf-8') as file:
    readme = file.read()

with open('LICENSE') as file:
    license = file.read()
    
setup(
  name = 'ROSSpy',      
  package_dir = {'ro':'rosspy'},
  packages = find_packages(),
  package_data = {
	'rosspy':[
            'water_bodies/*',
            'databases/*',
            'processed_databases/*', 
            'test/*', 
            'test/*/*', 
            'ro_module.json'
            ],
  },
  version = '0.1.7',
  license = license,
  description = "Software for predicting the brine concentrations and scaling quantities after RO desalination.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/ROSSpy',   
  keywords = ['desalination', 'reactive transport', 'geochemistry', 'sustainability', 'civil engineering', 'fluid mechanics'],
  install_requires = [
      'matplotlib', 
      'chempy', 
      'scipy', 
      'chemw', 
      'pubchempy', 
      'pandas', 
      'phreeqpy', 
      'sigfig'
      ],
  project_urls={
      'Documentation': 'https://rosspy.readthedocs.io/en/latest/index.html',
      'Issues': 'https://github.com/freiburgermsu/ROSSpy/issues',
  },
  classifiers={
      'Programming Language :: Python :: 3.7',
      'Programming Language :: Python :: 3.8',
      'Programming Language :: Python :: 3.9',
      'Programming Language :: Python :: 3.10',
      'Operating System :: POSIX :: Linux',
      'Operating System :: Microsoft :: Windows :: Windows 10',
      'Operating System :: MacOS',
      
      }
)

