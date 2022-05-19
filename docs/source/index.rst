Simulate Scale Formation and Brine Concentration during Reverse Osmosis Desalination
---------------------------------------------------------------------------------------------------------------------

|PyPI version| |docs| |Downloads| |docs| |License| 

.. |PyPI version| image:: https://img.shields.io/pypi/v/rosspy.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/ROSSpy/
   :alt: PyPI version

.. |Downloads| image:: https://pepy.tech/badge/rosspy
   :target: https://pepy.tech/project/rosspy
   :alt: Downloads

.. |Actions Status| image:: https://github.com/freiburgermsu/rosspy/workflows/Test%20ROSSpy/badge.svg
   :target: https://github.com/freiburgermsu/rosspy/actions
   :alt: Actions Status

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. |MyBinder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/freiburgermsu/rosspy/main?labpath=irosspy%2Firosspy.ipynb
   :alt: MyBinder
   
.. |docs| image:: https://readthedocs.org/projects/rosspy/badge/?version=latest
   :target: https://rosspy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
   

Desalinating ocean water is crucial for meeting the 6th UN Sustainable Development Goal of univeralizing potable water. Reverse Osmosis (RO) is the leading desalination technology, although, it remains hindered by membrane scaling, which lessens its energy efficiency and economic practicality. The geochemistry of mineral scaling evades experimental methods, and existing software programs to simulate scaling geochemistry -- e.g. French Creek -- are insufficient to explore all relevant variables and are furthermore inaccessible to many researchers. 

`Reverse Osmosis Scaling Software in Python (ROSSpy) <https://pypi.org/project/ROSSpy/>`_ satisfies this niche in RO research by implementating a one-dimesnional RO model through `PHREEQpy <https://pypi.org/project/phreeqpy/>`_, which is the Python version of `PHREEQC <https://www.usgs.gov/software/phreeqc-version-3>`_. The ``examples/scaling/scaling_validation`` directory of the `ROSSpy GitHub <https://github.com/freiburgermsu/ROSSpy>`_ details validation studies and sensitivity analyses of ROSSPy via Notebook examples. We encourage users and developers to critique and improve the open-source (`MIT License <https://opensource.org/licenses/MIT>`_) ROSSpy package by creating new `GitHub issues <https://github.com/freiburgermsu/ROSSpy/issues>`_ or emailing afreiburger@uvic.ca.

++++++++++++++++++++++
Theory
++++++++++++++++++++++

The ROSSpy framework represents RO desalination as a 1D reactive transport model of the membrane-solution interface in the RO feed channel. The feed solution is represented as a single, homogeneous, solution. The inlet boundary is defined by the Dirichlet condition, where the inlet concentrations are assumed to be independent of desalination, with the feed as an infinite reservoir. The outlet boundary is defined by the Cachy condition, where the effluent concentration is assumed to be dependent upon the reactive transport processes within the RO module. 


.. note::

   This project is under active development.


++++++++++++++++++++++
Installation
++++++++++++++++++++++

The IPHREEQC module must be installed before ROSSpy can be executed, since this is the source of geochemical calculations and data. The appropriate version of IPHREEQC can be installed from the `USGS <https://water.usgs.gov/water-resources/software/PHREEQC/index.html>`_ . Linux installation may require addition steps::

    wget https://water.usgs.gov/water-resources/software/PHREEQC/iphreeqc-3.7.3-15968.tar.gz
    tar -xzvf iphreeqc-3.7.3-15968.tar.gz
    cd iphreeqc-3.7.3-15968
    ./configure
    make
    make check
    sudo make install

Then either use the full path to the shared libray when making an instance of phreeqc_dll.IPhreeqc::

 phreeqc = phreeqpy.iphreeqc.phreeqc_dll.IPhreeqc\('/full/path/to/libiphreeqc.so'\)


or copy the shared object into phreeqpy/iphreeqc replacing the existing one, such as::

 sudo cp /usr/local/lib/libiphreeqc.so  /path/to/site-packages/phreeqpy/iphreeqc/libiphreeqc.so.0.0.0


ROSSpy per se is then installed in a command prompt, Powershell, Terminal, or Anaconda Command Prompt via ``pip``::

 pip install rosspy
   
   

Contents
--------

.. toctree::

   rosspy_usage
   rosspy_api
   parameter_files
   irosspy
