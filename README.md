[T]ime [V]ariant [O]ver[L]ap [A]dd in [P]artitions 
==============================

Reference: Jaeger, Bitzer, Simmer, Blau: "Echtzeitfähiges binaurales Rendering mit Bewegungssensoren von 3D-Brillen", Tagungsband DAGA 2017 - 43. Jahrestagung für Akustik, 6.-9. März 2017, Kiel, ISBN 978-3-939296-12-6


This project implements a weighted overlap add routine which works in partitions. Therefore, TVOLAP has low overall latency and is both efficient and suitable for time-variant filtering processes.


``TVOLAP.cpp`` and ``TVOLAP.h`` implement the processing routine as C++ class, with no external dependencies except the stdc++11 library. The TVOLAP processing class itself is licensed under the LGPLv3 (see ``COPYING.LESSER.txt`` for a copy of the license). The fast fourier transform routine (``fft.cpp`` and ``fft.h``) is available under the MIT License (see end of file for the copyrights and a copy of the license).

``TVOLAP.m`` and ``testTVOLAP.h`` in the MATLAB_Octave subfolder implement the processing routine as Octave class and show its usage. The implementation is also compatible with MATLAB (Tested with MATLAB r2016a and Octave 4.0.0). This implementation is available under the MIT license.

``pyTVOLAP.py`` and ``testTVOLAP.py`` in the Python subfolder implement the processing routine as Python class and show its usage. The example depends on bastibes soundfile library (https://github.com/bastibe/SoundFile). Tested with Python 2.7 and 3.6. This implementation is available under the MIT license.

``processExample.cpp`` stores a processing example of the convolution process by writing the resulting test signal to a .wav file, so that it can be easily loaded into MATLAB / Audacity a.s.o. for visualization purposes. This example is licensed under the LGPLv3 (see ``COPYING.LESSER.txt`` for a copy of the license).

Please have a look at the example codes to get an idea on how to use the TVOLAP class. C++ Project configuration is done via CMake. Project configuration is tested for Visual Studio 14 and Eclipse Neon with GNU / MinGW / MSVC compilers under Windows and Ubuntu. Other setups (like Debian / OSX) might work as well, but were not tested yet.


Build
-----

Run CMake without any errors and build the library and test binary with your favourite compiler.

Tested development toolchains:

- Windows 7, 8 and 10 with MinGW-w64-i686 (32Bit), MinGW-w64-x86_64 (64Bit) 6.3.0 compiler via Eclipse Neon and QtCreator, MSVC (32+64Bit) compiler via Visual Studio 14.
- Linux Ubuntu 16.04, Debian 8.7 with GCC via Eclipse Neon and QtCreator.

on x86-32 (Win32) and x86-64 processor architecture.

Output of the build is a shared library libTVOLAP and the test executable testTVOLAP.


Functionality
------------

The examples show a processing of the TVOLAP-class. 8Channel White noise is convolved with power complementary, switching bandpass filters, designed in frequency domain. The output is written to a .wav file, so you can easily visualize and/or play it. The Octave/MATLAB, as well as the Python example ``testTVOLAP.m`` / ``testTVOLAP.py`` in their subdirectories generate the same signal processing results.


If you like to use the TVOLAP class in a published project, you have two options: 

- Copy the relevant source code (``TVOLAP.cpp``, ``TVOLAP.h``, ``fft.cpp`` and ``fft.h``), include it in your project (or generate libTVOLAP static library and link against this) and include us as author of this (and only this) program part. Include the reference. This information should be clearly visible in your release. You have to copyleft your sources / license your software under a GPL.
- Link dynamically to the generated shared library libTVOLAP, include ``TVOLAP.h`` and make clearly visible, that this processing functionality is provided by us. Include the reference. Then, you do not have to licence your program under a GPL.

If there are any questions, please feel free to contact me: Email: hagenvontronje1@gmx.de
