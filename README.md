[W]eighted [O]ver[L]ap [A]dd in [P]artitions 
==============================

Reference: Jaeger, Bitzer, Simmer, Blau: "Echtzeitfähiges binaurales Rendering mit Bewegungssensoren von 3D-Brillen", DAGA 2017, Kiel

This project implements a weighted overlap add routine which works in partitions and therefore is both efficient and suitable for time-variant filtering processes.

``WOLAP.cpp`` and ``WOLAP.h`` implement the processing routine as C++ class, with no external dependencies except the stdc++11 library. The WOLAP processing class itself plus the fast fourier transform routine (``fft.cpp`` and ``fft.h``) are licenced under a LGPL (see ``COPYING.LESSER.txt`` for a copy of the licence).

``processExample.cpp`` stores a processing example of the convolution process by writing the resulting test signal to a .wav file, so that it can be easily visualized in MATLAB / Audacity a.s.o.). This example is licenced under a LGPL (see ``COPYING.LESSER.txt`` for a copy of the licence).

``QtExample.cpp`` with its GUI class ``QtMainwindow.cpp`` / ``QtMainwindow.h`` visualize a processing example by plotting the results via Qt and qcustomplot. This example is licenced under a GPL (see ``COPYING.txt`` for a copy of the licence).

Please see the example codes to get an idea on how to use the WOLAP class. There is also a MATLAB subfolder which shows how to implement WOLAP there. Tha MATLAB-example is not encapsulated in functions and scripts / objects yet.


Build
-----

Run Cmake without errors and build with your favourite compiler. For building the Qt example, run Cmake with -DBUILD_QT_EXAMPLE=ON, or toggle the option in the Cmake GUI if you use it. 

Tested setups:

- Windows 7, 8 and 10 with MinGW-w64-i686 (32Bit), MinGW-w64-x86_64 (64Bit) compiler
- Linux Ubuntu 16.04, Debian 8.7 with GCC

on x86-32 (Win32) and x86-64 processor architecture.


Functionality
------------

The examples show an processing of the WOLAP-class. A square wave is convolved with an inverse, delayed, decayed sine as impulse response. The Qt example additionally visualized the results in a plot. The MATLAB-example ``testWOLAP.m`` shows the same signal processing results.

If you like to use the WOLAP class in a published project, you have two options: 

- Copy the relevant source code (``WOLAP.cpp``, ``WOLAP.h``, ``fft.cpp`` and ``fft.h``), include it in your project (or generate libWOLAP static library and link against this) and include us as author of this (and only this) program part. Include the reference. This information should be clearly visible in your release. You have to licence your software under an GPL.
- Link dynamically to the generated shared library libWOLAP, include ``WOLAP.h`` and make clearly visible, that this processing functionality is provided by us. Include the reference. Then, you do not have to licence your program under an GPL.

If there are any questions please feel free to contact me: Email: hagenvontronje1@gmx.de, Tel. 015129140392