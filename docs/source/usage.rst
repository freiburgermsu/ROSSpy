Usage
=====

.. autosummary::
   :toctree: generated


A multitude of values are stored within the ``ROSSpy`` object, and can be subsequently used in a workflow. The complete list of content within the ``ROSSpy`` object can be identified and printed through the built-in ``dir()`` function in the following example sequence:

.. code-block:: python

 # conduct a ROSSpy simulation
 from rosspy import ROSSPkg
 ross = ROSSPkg(database_selection, simulation)
 ross.reactive_transport(simulation_time, simulation_perspective, final_cf)
 ross.feed_geochemistry(water_selection, water_characteristics)
 ross.execute()
 
 # evaluate the ROSSpy simulation contents
 print(dir(ross))

Accessible content
----------------------

The following list highlights stored content in the ``ROSSpy`` object after a simulation:

- *selected_output* & *processed_data* ``DataFrame``: `Pandas DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that possesses the raw and processed simulation data, respectively, from the PHREEQ simulation.
- *elemental_masses* ``dict``: A dictionary of the mass for each ion that constitutes precipitated scale. This is only determined for ``scaling`` simulations.
- *databases* ``list``: The available databases in the ``rosspy/databases/`` directory.
- *feed_sources* ``list``: The available feed waters in the ``rosspy/water_bodies/`` directory.
- *elements* & *minerals* ``dict``: Dictionaries of the elements and minerals, respectively, that are defined by the selected database.
- *cumulative_cf* ``float``: The final effluent CF after the entire simulation
- *input_file* ``str``: The complete PHREEQ input file of the simulation.
- *predicted_effluent* ``dict``: The predicted effluent concentrations of each ion that is defined in the feed.
- *parameters* & *variables* ``dict``: Dictionaries with the simulation parameters stored as key-value pairs.
- *results* ``dict``: A dictionary with the simulation results and each block of the simulation.
- *simulation_shifts* ``float``: The number of simulation shifts.
- *water_mw* & *water_gL* ``float``: The molecular weight and density of water, respectively.
- *chem_mw* ``ChemMW``: The ``ChemMW`` object from the `ChemW module <https://pypi.org/project/ChemW/>`_, which allows users to calculate the molecular weight from a string of any chemical formula. The formatting specifications are detailed in the README of the ChemW module. 