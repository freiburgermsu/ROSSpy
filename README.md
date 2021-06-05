# Reverse Osmosis Scaling Simulation (ROSS)

Desalination is an increasingly critical technology for human sustainability. Reverse Osmosis (RO) is a prominent desalination technology, yet, some of the geochemical phenomena like brine formation and mineral scaling hinder RO efficacy. We developed a RO scaling simulation (ROSS) software that effectively evaluates brine formation and mineral scaling from desalination. The source code and executable files for this software, which incorporates batch PHREEQC into the executable files, is provided to support community research and user customization.


# Input file generation 

The software first facilitates user parameterization of simulation conditions through intuitive command-line prompts.
The functions and parameterized code blocks for PHREEQC file are detailed in the fillowing sub-sections. 


## make_general_conditions

The "make_general_conditions" function notebook defines initial simulation conditions


## make_solutions

The "make_solutions" function defines the feed solution for the SOLUTION block of PHREEQC. The geochemical characteristics like elemental concentrations, alkalinity, and pH are parameterized with either user-defined quantities or from predefined  options like the Mediterranean and Red Seas or American produced groundwaters that are sourced from literature. 


## make_equilibrium_phases

The "make_equilibrium_phases" function defines mineral scaling for the EQUILIBRIUM_PHASES block. The user selects minerals from an alphabetized list of mineral names and chemical formulae that are defined by the user-selected database. The default set of displayed minerals are those that may be assembled from the defined set of elements in the "make_solutions" function. 


## make_reactive_transport

The "make_reactive_transport" function defines the reactive transport conditions for the REACTION and TRANSPORT blocks. Users may selected predefined parameters of the Dow-DuPont FILMTEC BW30-400 RO module, or users may assign custom spatiotemporal variables for the module characteristics like permeate flux, feed flow rate, membrane thickness, and module length, et cetera. The function includes numerous intermediary values and calculates final PHREEQC parameters from the user inputs.


## make_selected_output

The "make_selected_output" function defines the output file and data for the SELECTED_OUTPUT block. The selected set of simulation data is output to the working directory as a tab-delimited file that subsequently processed for simulation results. 


## export

The "export" function prints and exports the generated PHREEQC input file to the working directory. The export file is established with a predefined naming structure that incorporates simulation details, which prevents overwriting files, facilitates user estimation of the contained simulation, and is used by subsequent data processing code to generate the data figures. 


# PHREEQC file execution
The "execute" function executes PHREEQC for either previously created input files or the input file from the aforementioned code. The simulation progress and errors are printed in the same terminal window as the other components of the software, which facilitates user interpretation of the software.   


# Output data processing

The PHREEQC simulation output data may either be designed to view scaling over the module distance or brine formation over time. The two simulation perspectives are differentiated by the function, either through the filename of the output data file or through user specification. The two perspectives are then differentially processed into either figures of mineral scaling over modular distance or elemental concentration over time.


## process_selected_output

The "process_selected_output" function imports and interprets the SELECTED_OUTPUT data file from the working directory. The function then executes the appropriate function that processes the data into data figures. 


## make_brine_plot

The "make_brine_plot" function generates a figure of elemental concentration over time. The function first scans the elemental data for all non-zero elemental concentrations and then plots the concentrations on a logarithmic y-axis via the matplotlib Python library. The illustrated time commences after the concentration of the most concentrated ion, chloride, becomes non-zero, which represents the initial solution being flushed from the RO module. The generated figure of elemental concentrations contains a title, axes labels, a caption, and a ticker of the concentration factor (CF) from the corresponding simulation, which is a validation of the calculation accuracy. The function finally yields a data table of the average elemental concentrations, which augments identifying and quantifying the elemental concentrations of the figure. The elemental concentrations may be illustrated from scaling simulations as well, which may be contrasted with the figure of scaling from the simulation, however, these plots lack the data table of elemental concentrations. The function concludes by displaying the figure(s) and executing the "export_plot" function to export selected figures. 


## make_scaling_plot

The "make_scaling_plot" function generates figures of scaling over module distance. The function first quantifies the number of minerals that are precipitated during the simulation. The user then selects whether the minerals should be plotted simultaneously on a single figure, or whether the each mineral should be separately plotted; the latter is suggested to convey clear and organized simulation results. The figure legends depict the mineral name and chemical formula, and the time series every plot in the figure, which thereby illustrates scale accumulation over time. The function concludes by displaying the figure(s) and executing the "export_plot" function to export selected figures. 


## export_plot

The "export_plot" function exports the desired brine and scaling figures to the working directory. A default image format and name are provided for the figures which contain details that describe the essence of the corresponding simulation and prevent overwriting files, although, the user can provide any arbitrary export filename and image format from either PNG, JPG, and SVG. 
