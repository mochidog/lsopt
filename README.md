# lsopt
 Non-linear least square optimizer based on K. Madsen, H.B. Nielsen, O.Tingleff, "Methods for Non-Linear Least Squares Problems", 2014

Installation
=======================
This library requires clapack (tested version 3.2.1). A quick way to install clapack on a Linux machine is out-line as following:

1. Download clapack-3.2.1-CMAKE.tgz from Netlib.org and save to a/temp/dir
1. cd a/temp/dir; tar xvf clapack-3.2.1-CMAKE.tgz
1. cd clapack-3.2.1-CMAKE
1. Make the library using: cmake .; make
1. Create a folder to hold the library: sudo mkdir /usr/local/clapack-3.2.1
1. cd /usr/local/clapack-3.2.1
1. sudo cp -r a/temp/dir/INCLUDE/ .
1. sudo cp a/temp/dir/BLAS/SRC/libblas.a .
1. sudo cp a/temp/dir/SRC/liblapack.a .
1. sudo cp a/temp/dir/F2CLIBS/libf2c.a .
