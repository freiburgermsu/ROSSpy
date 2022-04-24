ROSSpy parameter files
-----------------------

Simulation parameters may be more succinctly provided through JSON files, which are imported by the code through identified function arguments, than through dictionaries that are passed as function arguments. Argument parameters can be used synergistically by supplanting specific characteristics that are different in the parameter files. 

ro_module
+++++++++++

The RO module characteristics can be provided through ``ro_module.json``. The default entry embodies the DOW FILMTEC BW30-400 RO module; however, other modules can be defined by emulating the structure of the default entry, where each characteristic possesses the respective units in the key name:

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
- *permeate_flow_m3_per_hour* & *max_feed_flow_m3_per_hour* ``dict``: specify the permeate and maximum feed flows through the RO module.
- *feed_thickness_mm* ``dict``: defines the thickness of the feed spacer through which the feed passes.
- *active_m2* ``dict``: defines the total filtration area of the RO module.
- *permeate_thickness_mm* ``dict``: defines the thickness of the permeate spacer.
- *membrane_thickness_mm* & *polysulfonic_layer_thickness_mm* & *support_layer_thickness_mm* ``dict``: define the thickensses of each layer in the composite filtration membrane: the polyamide layer that filters the feed, and the polysulfonic and support layers that provide resiliency to the membrane structure, respectively.
		
Other key:value sub-dictionaries may be introduced as metadata for the parameters.


water_bodies
+++++++++++++
 
The feed geochemistry may be defined as a parameter file, in the ``water_bodies`` folder within the ``rosspy`` package directory, that supplements feed characteristics as a dictionary through the ``feed_geochemistry`` function. The default files in this folder embody curated experimental data from both natural and produced water sources, which can be emulated for constructing files for other water sources, which each characteristic is expressed with the respective units in its key name:
      
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
	  
- *element* ``dict``: specifies all of the elements that are present in the feed, with sub-dictionaries of their concentrations and metadata. Some of these elements will not be amenable with some databases, which ROSSpy will simply ignore when defining the input file for an incompatible database.
- *temperature (C)*, *pe*, *Alkalinity*, & *pH* ``dict``: specify conditions and characteristics of the feed solution, with sub-directories of their respective value, chemical formula, and optionally metadata.