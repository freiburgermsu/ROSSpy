# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('rosspy/README_pypi.md') as file:
    readme = file.read()
    
setup(
  name = 'ROSSpy',      
  package_dir = {'core':'rosspy'},
  packages = find_packages(),
  package_data = {
	'rosspy':['water_bodies/*','databases/*', 'test/*', 'test/*/*', 'README_pypi.md'],
  },
  version = '0.0.4',
  license = 'GNU',
  description = "Software for predicting the brine and scaling consequences of a Desalination system.", 
  long_description = readme,
  long_description_content_type = "text/markdown",
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/ROSS',   
  keywords = ['desalination', 'reactive transport', 'geochemistry'],
  install_requires = ['matplotlib', 'chemicals', 'chempy', 'to_precision', 'scipy', 'pubchempy', 'pandas', 'phreeqpy']
)