
hgauss
=======

|Build-Status| |Win-Build-Status| |Codacy-Grade|

Cephes library.

Install
-------

The recommended way of installing it is via conda_

.. code:: bash

    conda install -c conda-forge hgauss

A second installation option would be to download the latest source and to
build it by yourself.
On Linux or macOS systems it can be as simple as

.. code:: bash

    bash <(curl -fsSL https://raw.githubusercontent.com/limix/hgauss/master/install)

On Windows systems you might want to have a look at the nmake_ command as a
replacement for ``make``.
Assuming you are at folder ``C:\projects\hgauss`` you might want to try

.. code:: dos

    mkdir build
    cd build
    cmake .. -G "NMake Makefiles" ^
             -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE
    nmake
    nmake test
    nmake install

Authors
-------

* `Danilo Horta`_

License
-------

This project is licensed under the MIT License - see the `license file`_ for
details.


.. |Build-Status| image:: https://travis-ci.org/limix/hgauss.svg?branch=master
    :target: https://travis-ci.org/limix/hgauss

.. |Win-Build-Status| image:: https://ci.appveyor.com/api/projects/status/enlpm6cfltjy7h3m/branch/master?svg=true
    :target: https://ci.appveyor.com/project/Horta/hgauss/branch/master

.. |Codacy-Grade| image:: https://api.codacy.com/project/badge/Grade/11201b1fc0d34d38bedeb218ad54dc6d
    :target: https://www.codacy.com/app/danilo.horta/hgauss?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=limix/hgauss&amp;utm_campaign=Badge_Grade

.. _conda: http://conda.pydata.org/docs/index.html

.. _license file: https://raw.githubusercontent.com/limix/hgauss/master/LICENSE.txt

.. _Danilo Horta: https://github.com/horta

.. _nmake: https://msdn.microsoft.com/en-us/library/dd9y37ha.aspx
