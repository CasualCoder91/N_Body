========
Database
========

A SQLite Database is used to store simulations including all relevant parameters, stars with their respective positions and velocities as well as data resulting from analysis.
SQLite supports both C++ and python. And is therefore used in most data related parts of the project. Generated data can be loaded and plotted directly from the database.
Previously performed simulations can be loaded from the database in order to carry out analysis like energy vs time or average star velocity vs time.

Due to bad performance, no foreign key constraints where set for position.id_star and velocity.id_star.
A multi-column index is used to improve the execution time of queries containing both velocity.timestep and velocity.id_star.
NULL values are avoided where possible. However, at many stages during execution, some information is unknown, but creating entries is still required.
Observed stars contain no information about their mass or GCA phase space coordinates, some cannot be mapped to simulated stars, etc.


Entity Relationship Diagram
---------------------------

.. image:: Images/ERD.*