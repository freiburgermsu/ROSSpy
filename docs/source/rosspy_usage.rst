Usage
=====

A complete ROSSpy simulation can be designed and executed through the following example sequence:

.. code-block:: python

 # conduct a ROSSpy simulation
 from rosspy import ROSSPkg
 ross = ROSSPkg(database_selection, simulation)
 ross.reactive_transport(simulation_time, simulation_perspective, final_cf)
 ross.feed_geochemistry(water_selection, water_characteristics)
 ross.execute()

Accessible content
----------------------

A multitude of values are stored within the ``ROSSpy`` object, and can be subsequently used in post-processing. The complete list of these values can be reviewed through the built-in ``dir()`` function, which is highlighted in the following list:

- *selected_output* & *processed_data* ``DataFrame``: `Pandas DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that possesses the raw and processed simulation data, respectively, from the PHREEQ simulation.
- *elemental_masses* ``dict``: A dictionary of the mass for each ion that constitutes scale. This is only determined for ``scaling`` simulations.
- *databases* & *feed_sources* ``list``: The available databases and feed waters in the ``rosspy/databases/`` and ``rosspy/water_bodies/`` directories, respectively.
- *elements* & *minerals* ``dict``: Dictionaries of the elements and minerals that are defined by the selected database, respectively.
- *cumulative_cf* ``float``: The final effluent CF after the entire simulation
- *ro_module* & *water_body* ``dict``: The complete sets of parameterized values for the simulated RO module and feed water, respectively.
- *input_file* ``str``: The complete PHREEQ input file that is executed for the simulation.
- *predicted_effluent* ``dict``: The predicted effluent concentrations of each feed ion.
- *parameters* & *variables* ``dict``: Dictionaries with the simulation parameters and calculated variables, respectively.
- *results* ``dict``: A dictionary with the simulation results and each block of the input file.
- *simulation_shifts* ``float``: The number of simulation shifts.
- *water_mw* & *water_gL* ``float``: The molecular weight and density of water, respectively.
- *chem_mw* ``ChemMW``: The ``ChemMW`` object from the `ChemW module <https://pypi.org/project/ChemW/>`_, which allows users to calculate the molecular weight from any chemical formula or chemical common name. The formatting specifications are detailed in the README of the ChemW module. 