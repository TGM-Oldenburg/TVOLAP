[W]eighted [O]ver[L]ap [A]dd in [P]artitions 
==============================

Reference: Jaeger, Bitzer, Simmer, Blau: "Echtzeitfähiges binaurales Rendering mit Bewegungssensoren von 3D-Brillen", Tagungsband DAGA 2017 - 43. Jahrestagung für Akustik, 6.-9. März 2017, Kiel, ISBN 978-3-939296-12-6


This project implements a weighted overlap add routine which works in partitions. Therefore, TVOLAP has low overall latency and is both efficient and suitable for time-variant filtering processes.


``TVOLAP.cpp`` and ``TVOLAP.h`` implement the processing routine as C++ class, with no external dependencies except the stdc++11 library. The TVOLAP processing class itself is licensed under the LGPLv3 (see ``COPYING.LESSER.txt`` for a copy of the license). The fast fourier transform routine (``fft.cpp`` and ``fft.h``) is available under the MIT License (see end of file for the copyrights and a copy of the license).

``TVOLAP.m`` and ``testTVOLAP.h`` in the MATLAB_Octave subfolder implement the processing routine as Octave class and show its usage. The implementation is also compatible with MATLAB (Tested with MATLAB r2016a and Octave 4.0.0). This implementation is available under the MIT license.

``pyTVOLAP.py`` and ``testTVOLAP.py`` in the Python subfolder implement the processing routine as Python class and show its usage. Tested with Python 2.7 and 3.6. This implementation is available under the MIT license.

``processExample.cpp`` stores a processing example of the convolution process by writing the resulting test signal to a .wav file, so that it can be easily loaded into MATLAB / Audacity a.s.o. for visualization purposes. This example is licensed under the LGPLv3 (see ``COPYING.LESSER.txt`` for a copy of the license).


``QtExample.cpp`` with its GUI class ``QtMainwindow.cpp`` / ``QtMainwindow.h`` directly visualizes a processing example by plotting the results via Qt and qcustomplot. This example is licensed under the GPLv3 (see ``COPYING.txt`` for a copy of the license).


Please have a look at the example codes to get an idea on how to use the TVOLAP class. Project configuring is done via CMake, there is a confiogure variable ``BUILD_QT_EXAMPLE`` which toggles the examples that will be build. Set it via the cmake GUI or command line option (-DBUILD_QT_EXQTAMPLE = ON / OFF) while configuring.


Build
-----

Run Cmake without errors and build with your favourite compiler. For building the Qt example, run Cmake with -DBUILD_QT_EXAMPLE=ON, or toggle the option in the Cmake GUI if you use it. 


Tested setups:

- Windows 7, 8 and 10 with MinGW-w64-i686 (32Bit), MinGW-w64-x86_64 (64Bit) compiler
- Linux Ubuntu 16.04, Debian 8.7 with GCC

on x86-32 (Win32) and x86-64 processor architecture.


Output of the build is a shared library libTVOLAP and the test executable(s).


Functionality
------------

The examples show a processing of the TVOLAP-class. A square wave is convolved with an inverse, delayed, decayed sine as impulse response. The Qt example additionally visualizes the results in a plot. The MATLAB-example ``testTVOLAP.m`` shows the same signal processing results.


If you like to use the TVOLAP class in a published project, you have two options: 

- Copy the relevant source code (``TVOLAP.cpp``, ``TVOLAP.h``, ``fft.cpp`` and ``fft.h``), include it in your project (or generate libTVOLAP static library and link against this) and include us as author of this (and only this) program part. Include the reference. This information should be clearly visible in your release. You have to copyleft your sources / license your software under a GPL.
- Link dynamically to the generated shared library libTVOLAP, include ``TVOLAP.h`` and make clearly visible, that this processing functionality is provided by us. Include the reference. Then, you do not have to licence your program under a GPL.

If there are any questions please feel free to contact me: Email: hagenvontronje1@gmx.de
