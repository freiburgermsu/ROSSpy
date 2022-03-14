ROSSpy parameter files
-----------------------

Module parameters may be more succinctly provided through JSON files, which are automatically imported by the code, than through function arguments; although, the function arguments can complement parameters from these JSON files. Each JSON file pertains to a distinct category of simulation parameters, which are individually detailed in the following sections.


ro_module
+++++++++++

The dimensions and specifications of the simulated RO module can be provided through the ``ro_module.json`` file, in addition to module characteristics that are provided through the dictionary argument of the function. The default entries of this file reflect the DOW FILMTEC BW30-400 RO module; however, other modules can be defined in this file by emulating the structure and content of the default entry in the file. Each characteristic of the RO module is expressed with the respective units in its key name of the JSON file.

.. code-block:: json

 {
    "BW30-400":{
        "module_diameter_mm": {
            "value": 201
        },
        "permeate_tube_diameter_mm": {
            "value": 29
        },
        "module_length_m": {
            "value": 1.016
        },
        "permeate_flow_m3_per_hour": {
            "value": 1.667
        },
        "max_feed_flow_m3_per_hour": {
            "value": 15.9
        },
        "feed_thickness_mm": {
            "value": 0.8636
        },
        "active_m2": {
            "value": 37
        },
        "permeate_thickness_mm": {
            "value": 0.3
        },
        "membrane_thickness_mm": {
            "value": 0.25
        },
        "polysulfonic_layer_thickness_mm": {
            "value": 0.05
        },
        "support_layer_thickness_mm": {
            "value": 0.15
        }      
    }
 }

- *module_diameter_mm* ``dict``: specifies the total diameter of the RO module.
- *permeate_tube_diameter_mm* ``dict``: specifies the diameter of the permeate tube within the RO module.
- *module_length_m* ``dict``: specifies the total length of the RO module.
- *permeate_flow_m3_per_hour* & *max_feed_flow_m3_per_hour* ``dict``: specifies the flow of permeate and the maximum flow of feed through the RO module.
- *feed_thickness_mm* ``dict``: defines the thickness of the feed spacer through which the feed solution passes in the RO module.
- *active_m2* ``dict``: defines the total, cumulative, filtration area of the RO module.
- *permeate_thickness_mm* ``dict``: defines the thickness of the permeate spacer.
- *membrane_thickness_mm* & *polysulfonic_layer_thickness_mm* & *support_layer_thickness_mm* ``dict``: define the thickensses of each layer in the composite filtration membrane: the polyamide layer that filters the feed, and the polysulfonic and support layers that provide resiliency to the membrane structure.
		
Other key:value pairs may be introduced alongside the ``value`` key to store metadata about the parameterized module:


water_bodies
+++++++++++++
 
The feed geochemistry of the simulated water source may be defined through a JSON file, in the ``water_bodies`` folder within the ``rosspy`` package directory, in addition to defining module characteristics via the dictionary argument of the ``feed_geochemistry`` function. The default files in this folder embody curated experimental data from both natural and produced water sources; however, other water sources can be defined by emulating the structure and content of the default file options. Each characteristic of the RO module is expressed with the respective units in its key name of the JSON file.
      
.. code-block:: json
		
 {
    "element": {
        "Mn": {
            "concentration (ppm)": 3000,
            "reference": "Haluszczak, Rose, and Kump, 2013 [estimated from another Marcellus publication]"
        },
        "Fe": {
            "concentration (ppm)": 26.6,
            "reference": "Chapman et al., 2012"
        },
        "B": {
            "concentration (ppm)": 20,
            "reference": "Haluszczak, Rose, and Kump, 2013 [reported average form another Marcellus publication]"
        },
        "Cl": {
            "concentration (ppm)": 81900,
            "reference": "Chapman et al., 2012"
        },
        "Na": {
            "concentration (ppm)": 32800,
            "reference": "Chapman et al., 2012"
        },
        "S(6)": {
            "concentration (ppm)": 45,
            "reference": "Haluszczak, Rose, and Kump, 2013 [estimated from another Marcellus publication]"
        },
        "Ca": {
            "concentration (ppm)": 8786,
            "reference": "Chapman et al., 2012"
        },
        "K": {
            "concentration (ppm)": 350,
            "reference": "Haluszczak, Rose, and Kump, 2013 [estimated from another Marcellus publication]"
        },
        "Mg": {
            "concentration (ppm)": 841,
            "reference": "Chapman et al., 2012"
        },
        "Sr": {
            "concentration (ppm)": 2415,
            "reference": "Chapman et al., 2012"
        },
        "Ba": {
            "concentration (ppm)": 962,
            "reference": "Chapman et al., 2012"
        },
        "Li": {
            "concentration (ppm)": 95,
            "reference": "Haluszczak, Rose, and Kump, 2013 [reported average from another Marcellus publication]"
        }
    },
    "temperature (C)": {
        "value": 24,
        "reference": "Dresel and Rose, 2010"
    },
    "pe": {
        "value": null,
        "reference": null
    },
    "Alkalinity": {
        "value": 71,
        "reference": "Haluszczak, Rose, and Kump, 2013 [reported average from another Marcellus publication]",
	"form": "CaCO3"
    },
    "pH": {
        "value": 7,
        "reference": "Haluszczak, Rose, and Kump, 2013 [estimated from another Marcellus publication]"
    }
 } 
	  
- *element* ``dict``: specifies all of the elements that are present in the simulated water source, with sub-dictionaries of their concentrations and metadata. Some of these elements will not be amenable with some databases, nevertheless, the elements can be parameterized and the ROSSpy code will not simulate these elements.
- *temperature (C)*, *pe*, *Alkalinity*, & *pH* ``dict``: specify conditions and characteristics of the feed solution, with sub-directories of their respective value, chemical form where it is applicable, and metadata.