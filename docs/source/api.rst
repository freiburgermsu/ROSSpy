ROSSpy API
--------------

.. autosummary::
   :toctree: generated


++++++++++++++++++++++
Installation
++++++++++++++++++++++

ROSSpy is installed in a command prompt, Powershell, Terminal, or Anaconda Command Prompt via ``pip``::

 pip install rosspy

++++++++++++++++++++++
__init__
++++++++++++++++++++++

The simulation environment is defined:

.. code-block:: python

 import rosspy
 ross = rosspy.ROSSPkg(database_selection, simulation = 'scaling', simulation_type = 'transport', 
 domain_phase = None, quantity_of_modules = 1, simulation_title = None, verbose = False, printing = True, jupyter = False)

- *database_selection* ``str``: specifies which PHREEQ database file -- ``Amm``, ``ColdChem``, ``core10``, ``frezchem``, ``iso``, ``llnl``, ``minteq``, ``minteq.v4``, ``phreeqc``, ``pitzer``, ``sit``, ``Tipping_Hurley``, or ``wateq4f`` -- will be imported and used to execute the simulation. These databases were all processed via the ``PHREEQdb`` object of the `ChemW module <https://pypi.org/project/ChemW/>`_, in this specific Notebook: `here <https://github.com/freiburgermsu/ChemW/blob/main/examples/PHREEQ/PHREEQ%20databases.ipynb>`_.
- *simulation* ``str``: specifies whether the ``scaling`` or ``brine`` of the simulation will be evaluated.
- *simulation_type* ``str``: specifies whether RO reactive transport ``transport``, or the geochemistry of ``evaporation``, will be simulated with the parameterized feed solution.
- *operating_system* ``str``: specifies whether the user is using a ``windows`` or ``unix`` system, which directs importing PHREEQpy and commenting in the ``PQI`` PHREEQ input files.
- *domain_phase* ``str``: specifies whether the ``mobile`` (i.e. bulk solution) or the ``immobile`` (i.e. the CP solution layer) will be evaluated for dual domain simulations. Parameterizing an argument other than ``None`` implicitly signifies that simulation of the dual-domain model, as opposed to the default single-domain model.  
- *quantity_of_modules* ``int``: specifies the number of in-series RO modules that will be simulated.
- *simulation_title* ``str``: specifies the title of the simulation, which is only observed in the PHREEQC ``PQI`` input file.
- *verbose* ``bool``: specifies whether simulation details and calculated values will be printed. This is valuable for trobuleshooting.
- *printing* ``bool``: specifies whether simulation results will be printed. 
- *jupyter* ``bool``: specifies whether the simulation is being conducted in a Jupyter Notebook, which allows ``display()`` to illustrate data tables and figures.

++++++++++++++++++++++
reactive_transport()
++++++++++++++++++++++

The spatiotemporal conditions and the permeate flux, or evaporation rate, are defined for reactive transport:

.. code-block:: python

 ross.reactive_transport(simulation_time, simulation_perspective = None, final_cf = None, 
 module_characteristics = {}, permeate_efficiency = 1, head_loss = 0.89, evaporation_steps = 15, 
 timestep = None, cells_per_module = 12, kinematic_flow_velocity = None, exchange_factor = 1e5)

- *simulation_time* ``float``: specifies the total simulated time in seconds.
- *simulation_perspective* ``str``: specifies whether the simulation data is parsed to view the end of the module over the entire simulated time "all_time" or to view the entire module distance at the final timestep ``all_distance``. The ``None`` parameter defaults to ``all_time`` for brine simulations and ``all_distance`` for scaling simulations.
- *final_cf* ``float``: specifies the effluent CF in the simulated RO system. The ``None`` parameter indicates that the ``linear_permeate`` permeate flux method will be simulated, while any numerical value indicates that a ``linear_cf`` permeate flux method will be simulated. 
- *module_characteristics* ``dict``: specifies custom RO specifications that diverge from those of the default DOW FILMTEC BW30-400 RO module. The expected dictionary keys are: ``module_diameter_mm``, ``permeate_tube_diameter_mm``, ``module_length_m``, ``permeate_flow_m3_per_day``, ``max_feed_flow_m3_per_hour``, ``membrane_thickness_mm``, ``feed_thickness_mm``, ``active_m2``, ``permeate_thickness_mm``, ``polysulfonic_layer_thickness_mm``, ``support_layer_thickness_mm``. Any specifications that are not defined in the parameterized dictionary defaults to that specification of the FILMTEC BW30-400 module. The dictionary values are all floats in the units from the end of the corresponding key:

.. code-block:: json

 {
 "active_m2": 37,
 "permeate_thickness_mm": 0.3,
 "polysulfonic_layer_thickness_mm": 0.05
 }


- *permeate_efficiency* ``float``: specifies the 0<=PE<=1 proportion of calculated permeate flux that is simulated: e.g. ``PE=1`` denotes a perfectly operational module and ``PE=0.5`` denotes a 50% operational module, etc. 
- *head_loss* ``float``: specifies the 0<=HL<=1 proportion of effluent pressure relative to the influent. The `default value of 0.89 <https://doi.org/10.1063/1.3109795>`_ corresponds to an 11% pressure drop.
- *cells_per_module* ``int``: specifies the quantity of cells into which the RO module is discretized. This controls the resolution of data over the distance of the module, and thus is only consequential for ``simulation_perspective = "all_distance"`` simulations.
- *kinematic_flow_velocity* ``float``: specifies the kinetic flow velocity of the feed solution. The ``None`` parameter defaults to 9.33E-7 (m^2/sec).

++++++++++++++++++++++
feed_geochemistry()
++++++++++++++++++++++

The feed geochemistry is defined from either a parameter file in the ``rosspy/water_bodies`` directory or in a dictionary argument, from which the potential scalants are determined:

.. code-block:: python

 ross.feed_geochemistry(water_selection = '', water_characteristics = {}, solution_description = '', ignored_minerals = [], existing_parameters = {}, parameterized_ph_charge = True)

- *water_selection* ``str``: specifies a feed water from the *rosspy/water_bodies* directory, where default parameter files for natural waters -- the ``red_sea`` and the ``mediterranean_sea`` -- and produced waters of fracking oil wells -- the ``bakken_formation``, ``marcellus_appalachian_basin``, ``michigan_basin``, ``north_german_basin``, ``palo_duro_basin``, and ``western_pennsylvania_basin`` -- are provided. Parameter files for other feed waters can be created by emulating the syntax of these default files and storing the created file in the aforementioned directory.
- *water_characteristics* ``dict``: defines the geochemistry and conditions of a custom feed solution. The expected keys are: ``element``, ``temperature (C)``, ``pe``, ``Alkalinity``, and ``pH``. Each value of these keys is itself a dictionary, with the keys of ``value`` for the numerical value and ``reference`` to denote an experimental citation for the numerical value. The ``element`` key deviates slightly from this organization, by using another dictionary layer for each ion in the feed. The keys here, for each ion, are: ``concentration (ppm)`` for its ppm concentration, ``form`` for the mineral form or charge-state of the ion, and ``reference`` with the same aforementioned purpose. The following dictionary illustrates this organization:

.. code-block:: json

 {
    "element": {
        "Mn": {
            "concentration (ppm)": 0.000734,
            "reference": "El Sayed, Aminot, and Kerouel, 1994"
        }, 
        "Si": {
            "concentration (ppm)": 95,
            "reference": "Haluszczak, Rose, and Kump, 2013",
            "form": "SiO2"
        }
    },
    "temperature (C)": {
        "value": 24,
        "reference": "Dresel and Rose, 2010"
    }
 }

- *solution_description* ``str``: briefly describes the solution, which will be used in the simulation folder name in the absence of a parameterized *water_selection*.
- *ignored_minerals* ``list``: defines the minerals that will be excluded from the determined set of minerals that can potentially precipitate from the parameterized feed ions.
- *existing_parameters* ``dict``: specifies pre-existing equilibria conditions that influence the geochemical calculations of PHREEQ. The expected keys are the referenced mineral names, with values of ``saturation`` and ``initial_moles`` that correspond to the pre-existing saturation index and the initial moles, respectively, of the simulated mineral.
- *parameterized_ph_charge* ``bool``: specifies whether the pH will be charged balance, where ``True`` prevents the parameterization of alkalinity in the feed solution. 



++++++++++++++++++++++
parse_input()
++++++++++++++++++++++

This function is used to parse and execute pre-existing input file:

.. code-block:: python

 ross.parse_input(input_file_path, simulation, water_selection = None, simulation_name = None, active_m2= None)

- *input_file_path* ``str``: specifies the path of the existing input file that will be parsed and executed. 
- *simulation* ``str``: specifies whether ``scaling`` or ``brine`` will be processed from the simulation. 
- *water_selection* ``str``: describes the simulated feed water. 
- *simulation_name* ``str``: specifies the name of the simulation folder to which all of the simulation files will be exported, where ``None`` defaults to a naming scheme that is designed by the software with unique simulation details. 
- *active_m2* ``float``: defines the area of active filtration in the simulated RO module, where ``None`` defaults to 37 from the standard FILMTEC BW30-400 module. 


++++++++++++++++++++++
execute()
++++++++++++++++++++++

The input file is executed through PHREEQ:

.. code-block:: python

 processed_data = ross.execute(simulation_name = None, selected_output_path = None, simulation_path = None, 
 plot_title = None, title_font = 'xx-large', label_font = 'x-large', x_label_number = 6, export_name = None, 
 export_format = 'svg', individual_plots = None, scale_ions = True, selected_output_filename = None,)

- *simulation_name* ``str``: specifies the name of a PHREEQ output file.
- *selected_output_path* ``str``: specifies the path of a simulation output file that will be processed into data tables and figures. This imported file can be independent of executing ROSSpy, and thus can be used to process old data. This parameter must be ``None`` to execute PHREEQ input files.
- *simulation_path* ``str``: The path to where the simulation content will be saved, where ``None`` signifies the current working directory.
- *plot_title* ``str``: specifies the title of the simulation figure, where ``None`` defaults to customized titles that incorporate unique simulation details: e.g. ``scaling`` or ``brine``, the water body, and the total simulation time.
- *title_font* & *label_font* ``str``: these specify the fonts of the figure title and axis labels, respectively, in terms of MatPlotLib font specifications: ``xx-small``, ``x-small``, ``small``, ``medium``, ``large``, ``x-large``, or ``xx-large``. 
- *x_label_number* ``int``: quantifies the ticks along the x-axis of the simulation figure.
- *export_name* ``str``: specifies the export name of the simulation figure. The default name for ``brine`` simulations is ``brine`` . The default names for ``scaling`` simulations, depending upon a ``True`` or ``False`` value of the *individual_plots* argument, is an individual mineral name (e.g. ``Gypsum``) or ``all_minerals``, respectively.
- *export_format* ``str``: specifies the format of the exported simulation figure, from the MatPlotLib options: ``svg``, ``pdf``, ``png``, ``jpeg``, ``jpg``, or ``eps``. The default is ``svg``, which is a lossless format that is highly customizable in software like `Inkscape <https://inkscape.org/>`_.
- *individual_plots* ``bool``: specifies whether each mineral of ``scaling`` simulations are plotted individually, or whether each scalant is plotted in a combined single figure. The ``None`` parameter defaults to ``True`` for the "all_time" *simulation_perspective* and ``False`` otherwise.
- *scale_ions* ``bool``: specifies whether the scale from ``scaling`` simulations will be refined into quantities of individual ions that constitute the mineral scale. This information of ionic quantities is exported as a JSON file to the simulation folder. The default value is ``True``.

- *simulation_name* ``str``: specifies the name of the simulation folder to which simulation content will be exported. The ``None`` parameter assigns a default name for the simulation folder, which follows the format of **today's_date-ROSSpy-water_selection-simulation_type-database_selection-simulation-simulation_perspective-#**. 
- *input_path* & *output_path* ``str``: specifies the directory path to where the input file will be exported, where ``None`` defaults to "input.pqi" and "selected_output.csv", respectively, in the current working directory. 

**Returned** *processed_data* ``DataFrame``: `Pandas DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that possesses the processed simulation data.

++++++++++++++++++++++
test()
++++++++++++++++++++++

ROSSpy can be tested with a simple built-in ``test()`` function, which can be executed through these three lines:

.. code-block:: python

 import rosspy
 ross = rosspy.ROSSPkg(database_selection, simulation)
 ross.test()

The ``Test()`` function executes a predefined sample simulation to exemplify ROSSpy with a simple use case.

