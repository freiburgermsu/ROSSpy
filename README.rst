Simulate Scale Formation and Brine Concentration during Reverse Osmosis Desalination
---------------------------------------------------------------------------------------------------------------------

|PyPI version| |Downloads| |docs| |License| 

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


`Reverse Osmosis Scaling Software in Python (ROSSpy) <https://pypi.org/project/ROSSpy/>`_ offers an open-source API to simulate the reactive transport geochemistry of Reverse Osmosis desalination. ROSSpy essentially translates user specifications of an RO system into `PHREEQpy <https://pypi.org/project/phreeqpy/>`_, which is the Python version of `PHREEQC <https://www.usgs.gov/software/phreeqc-version-3>`_. The ``examples/scaling/scaling_validation`` directory of the `ROSSpy GitHub <https://github.com/freiburgermsu/ROSSpy>`_ details the numerous functions and accuracy of ROSSPy via Notebook examples. We encourage users and developers to critique and improve ROSSpy, as an open-source (`MIT License <https://opensource.org/licenses/MIT>`_) library, through `GitHub issues <https://github.com/freiburgermsu/ROSSpy/issues>`_.

The complete documentation is provided by `ReadTheDocs <https://rosspy.readthedocs.io/en/latest/index.html>`_.


++++++++++++++++++++++
Installation
++++++++++++++++++++++

ROSSpy is installed in a command prompt, Powershell, Terminal, or Anaconda Command Prompt via ``pip``::

 pip install rosspy

The IPHREEQC module must then be installed, since this is the source of geochemical calculations and data for ROSSpy. The appropriate version of IPHREEQC can be installed from the `USGS <https://water.usgs.gov/water-resources/software/PHREEQC/index.html>`_ . **Linux** installation may require addition steps::

    wget https://water.usgs.gov/water-resources/software/PHREEQC/iphreeqc-3.7.3-15968.tar.gz
    tar -xzvf iphreeqc-3.7.3-15968.tar.gz
    cd iphreeqc-3.7.3-15968
    ./configure
    make
    make check
    sudo make install
    pip show phreeqpy
    sudo cp /usr/local/lib/libiphreeqc.so  /path/to/site-packages/phreeqpy/iphreeqc/libiphreeqc.so.0.0.0
