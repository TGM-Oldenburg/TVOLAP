[W]eighted [O]ver[L]ap [A]dd in [P]artitions 
==============================

This project implements a weighted overlap add routine which works in partitions and therefore is both efficient and suitable for time-variant filtering processes.

``WOLAP.cpp`` and ``WOLAP.h`` implement the processing routine as C++ class, with no external dependencies except the stdc++11 library. The WOLAP processing class itself plus the fast fourier transform routine (``fft.cpp`` and ``fft.h``) are licenced under an LGPL (see ``COPYING.LESSER.txt`` for a copy of the licence).

``processExample.cpp`` show a processing example with writing the result to a .wav file, so that it can easily be visualized in MATLAB / Audacity a.s.o.). This example is licenced under an LGPL (see ``COPYING.txt`` for a copy of the licence).

``QtExample.cpp`` plus its GUI class ``QtMainwindow.cpp`` / ``QtMainwindow.h`` show a processing example with plotting of the results. This example is licenced under an GPL (see ``COPYING.txt`` for a copy of the licence).

Please see the example codes to get an idea on how to use the WOLAP class. There is also a MATLAB subfolder which shows how to implement WOLAP there. Tha MATLAB-example is not encapsulated in functions and scripts or objects yet.


Build
-----

Run cmake without errors and build with your favourite compiler. Tested setups: 

- Windows7, 8 and 10 with MinGW-w64-i686 (32Bit), MinGW-w64-x86_64 (64Bit)
- Linux Ubuntu 16.04, Debian 8.7 with GCC

on x86-32 (Win32) and x86-64 processor architecture.


Functionality
------------

The examples show an processing of the WOLAP-class. A square wave is convolved with an inverse, delayed, decayed sine as impulse response. The Qt example additionally visualized the results in a plot. The MATLAB-example ``testWOLAP.m`` shows the same signal processing results.

If there are any questions please feel free to contact me: Email: hagenvontronje1@gmx.de, Tel. 015129140392