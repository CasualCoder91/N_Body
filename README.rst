Download and Installation
=========================

Documentation
-------------

.. image:: https://readthedocs.org/projects/n-body/badge/?version=latest
   :target: https://n-body.readthedocs.io/en/latest/?badge=latest
   :alt: Here be dragons!

`https://n-body.readthedocs.io/en/latest/?badge=latest
<https://n-body.readthedocs.io/en/latest/?badge=latest>`_

Clone repository
----------------
git clone https://github.com/CasualCoder91/N_Body.git

Windows
-------

Navigate into the build directory and use Microsoft Build Engine to build the application
 ``cd build
 
 cmake --config Release .. 
 
 "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\MSBuild\Current\Bin\MSBuild.exe" ALL_BUILD.vcxproj /p:Configuration=Release``
