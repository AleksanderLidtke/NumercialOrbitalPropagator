"""
12/19/2013
Author: Joshua Milas
Python Version: 3.3.2

The NRLMSISE-00 model 2001 ported to python
Based off of Dominik Brodowski 20100516 version available here
http://www.brodo.de/english/pub/nrlmsise/

This is the test program, and the output should be compaired to
/* -------------------------------------------------------------------- */
/* ---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ---------- */
/* -------------------------------------------------------------------- */

/* This file is part of the NRLMSISE-00  C source code package - release
 * 20041227
 *
 * The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
 * Doug Drob. They also wrote a NRLMSISE-00 distribution package in
 * FORTRAN which is available at
 * http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
 *
 * Dominik Brodowski implemented and maintains this C version. You can
 * reach him at mail@brodo.de. See the file "DOCUMENTATION" for details,
 * and check http://www.brodo.de/english/pub/nrlmsise/index.html for
 * updated releases of this package.
 */
"""

import time
from nrlmsise_00_header import *
from nrlmsise_00 import *

output = [nrlmsise_output() for _ in range(17)]
Input = [nrlmsise_input() for _ in range(17)]
flags = nrlmsise_flags()
aph = ap_array()

for i in range(7):
    aph.a[i]=100
flags.switches[0] = 0
for i in range(1, 24):
    flags.switches[i]=1
    
for i in range(17):
    Input[i].doy=172;
    Input[i].year=0; #/* without effect */
    Input[i].sec=29000;
    Input[i].alt=400;
    Input[i].g_lat=60;
    Input[i].g_long=-70;
    Input[i].lst=16;
    Input[i].f107A=150;
    Input[i].f107=150;
    Input[i].ap=4;
	
Input[1].doy=81;
Input[2].sec=75000;
Input[2].alt=1000;
Input[3].alt=100;
Input[10].alt=0;
Input[11].alt=10;
Input[12].alt=30;
Input[13].alt=50;
Input[14].alt=70;
Input[16].alt=100;
Input[4].g_lat=0;
Input[5].g_long=0;
Input[6].lst=4;
Input[7].f107A=70;
Input[8].f107=180;
Input[9].ap=40;

Input[15].ap_a = aph
Input[16].ap_a = aph

#evaluate 0 to 14
for i in range(15):
    gtd7(Input[i], flags, output[i])

#/* evaluate 15 and 16 */
flags.switches[9] = -1
for i in range(15, 17):
    gtd7(Input[i], flags, output[i])

#/* output type 1 */
for i in range(17):
    print '\n'
    for j in range(9):
        print'%E ' % output[i].d[j]
    print '%E ' % output[i].t[0]
    print '%E ' % output[i].t[1]
    #/* DL omitted */

#/* output type 2 */
for i in range(3):
    print '\n'
    print "\nDAY   "
    for j in range(5):
        print "         %3i" % Input[i*5+j].doy
    print "\nUT    "
    for j in range(5):
        print "       %5.0f" % Input[i*5+j].sec
    print "\nALT   "
    for j in range(5):
        print "        %4.0f" % Input[i*5+j].alt
    print "\nLAT   "
    for j in range(5):
        print "         %3.0f" % Input[i*5+j].g_lat
    print "\nLONG  "
    for j in range(5):
        print "         %3.0f" % Input[i*5+j].g_long
    print "\nLST   "
    for j in range(5):
        print "       %5.0f" % Input[i*5+j].lst
    print "\nF107A "
    for j in range(5):
        print "         %3.0f" % Input[i*5+j].f107A
    print "\nF107  "
    for j in range(5):
        print "         %3.0f" % Input[i*5+j].f107

    print '\n\n'
    
    print "\nTINF  "
    for j in range(5):
        print "     %7.2f" % output[i*5+j].t[0]
    print "\nTG    "
    for j in range(5):
        print "     %7.2f" % output[i*5+j].t[1]
    print "\nHE    "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[0]
    print "\nO     "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[1]
    print "\nN2    "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[2]
    print "\nO2    "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[3]
    print "\nAR    "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[4]
    print "\nH     "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[6]
    print "\nN     "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[7]
    print "\nANM   "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[8]
    print "\nRHO   "
    for j in range(5):
        print "   %1.3e" % output[i*5+j].d[5]
    print '\n'