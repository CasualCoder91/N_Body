==============
DATABASE
==============

SQLite Database is used to store simulations including all relevant parameters, stars with their respective positions and velocities as well as data resulting from analysis.
Simulations, which have been performed, can be loaded from the database to perform analysis like energy vs time or average star velocity vs time.
The separation of simulation and analysis is the main reason this database exists, since performing analysis during simulation would lead to ugly code or low performance.
Another reason being that loading data from files is inefficient depending on the queried information.

Entity Relationship Diagram
---------------------------

.. image:: ..Images/ERD.png

Database Interface
------------------

Interactions with the database are handled by the Database class

.. doxygenclass:: Database
   :members:
