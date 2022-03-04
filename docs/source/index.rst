Simulating Scale Formation during Reverse Osmosis Desalination
---------------------------------------------------------------------------------------------------------------------

|PyPI version| |Downloads| |License|

.. |PyPI version| image:: https://img.shields.io/pypi/v/rosspy.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/ROSSpy/
   :alt: PyPI version

.. |Actions Status| image:: https://github.com/freiburgermsu/rosspy/workflows/Test%20ROSSpy/badge.svg
   :target: https://github.com/freiburgermsu/rosspy/actions
   :alt: Actions Status

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/rosspy
   :target: https://pepy.tech/project/rosspy
   :alt: Downloads

Desalination is an unavoidable technology for meeting the 6th UN Sustainable Development Goal of providing potable for all people. Reverse Osmosis (RO) is the leading desalination technology, although, it can be further improved in energy efficiency and economic practicality by mitigating membrane fouling like mineral scaling. The geochemistry of mineral scaling is generally inaccessible to physical experimentation, and existing software programs to simulate scaling geochemistry -- e.g. French Creek -- are esoteric and/or financially expensive. 

`Reverse Osmosis Scaling Software in Python (ROSSpy) <https://pypi.org/project/ROSSpy/>`_ seeks to satisfy this void. ROSSpy essentially translates user specifications of an RO system into `PHREEQpy <https://pypi.org/project/phreeqpy/>`_, which is the Python version of `PHREEQC <https://www.usgs.gov/software/phreeqc-version-3>`_. The ``examples/scaling/scaling_validation`` directory of the `ROSSpy GitHub <https://github.com/freiburgermsu/ROSSpy>`_ details the numerous functions and accuracy of ROSSPy via Notebook examples. We encourage users and developers to critique and improve ROSSpy, as an open-source (`MIT License <https://opensource.org/licenses/MIT>`_) library, through `GitHub issues <https://github.com/freiburgermsu/ROSSpy/issues>`_.

++++++++++++++++++++++
Theory
++++++++++++++++++++++

The ROSSpy framework represents RO desalination as a 1D reactive transport model of the membrane-solution interface in the RO feed channel. The feed solution can be represented either by a single, homogeneous, solution domain. The inlet boundary is defined by the Dirichlet condition, where the feed is assumed to be an infinite reservoir. The outlet boundary is defined by the Cachy condition, where the effluent is assumed to be dependent upon the reactive transport processes within the RO module. 


.. note::

   This project is under active development.
   
   

Contents
--------

.. toctree::

   rosspy_usage
   rosspy_api
