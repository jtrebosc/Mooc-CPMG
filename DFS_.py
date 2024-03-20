#print "execute ssDFS"
import os
from math import *
import cmath
import argparse
import sys


description = """Generates a shape pulse with Double Frequency Sweep (DFS)
The sweep consists of symmetrical with respect to carier frequency. Direction of sweep can be out high to low or low to high frequency.
DFS parameter can be read from dataset parameter or from commnd line argument.
"""

class dummy():
    """A class to mimick argparse object in case argparse is not available."""
    def __init__(self):
        self.spname = False # -s
        self.spname_file = False # -f
        self.end = False # -r
        self.begin = False # -o
        self.resolution = False # -n
        self.pulse_length = False # -p
        self.sweep_direction = False # -d
        self.wurst_cutoff = False # -c
# input parameters
# shape_name
# shape_filename
# rotation_speed
# offset
# resolution
# pulse_length
# sweep_direction
# wurst_cutoff

try : 
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-s', '--spname',
                help='Shape name number to update. Default: 8 for SPNAM 8.')
    parser.add_argument('-f', '--spname_file', 
                help='File where to store shape. Default: ssDFS.')
    parser.add_argument('-n', '--resolution',
                help='Resolution of the shape in ns. Default: value stored in CNST 3.')
    parser.add_argument('-b', '--begin',
                help='Start of frequency sweep in kHz. Default: value stored in CNST 1.')
    parser.add_argument('-e', '--end',
                help='End of frequency sweep in kHz. Default: value stored in CNST 2.')
    parser.add_argument('--pulse_name', '-pn', 
                help='Name number of the pulse. Default: 8 for P 8.')
    parser.add_argument('--pulse_length', '-pl', 
                help='Length of the shape. Default: value stored in pulse_name (e.g. P 8).')
    parser.add_argument('--wurst_cutoff', '-c', 
                help='WURST cutoff value for shape enveloppe. Default: value stored in CNST 2.')

    args  =  parser.parse_args(sys.argv[1:])
except ImportError:
    if len(sys.argv) > 1:
        MSG("Argparse module not found!\n Arguments won't be processed")
    args = dummy()
except SystemExit:
    MSG(""" Script is exiting : either you asked for help or there is an argument error.
    Check console for additional information
    """  + parser.format_help() )
    EXIT()
    
if args.spname is None:
    SPNAME = "SPNAM 8"
else:
    SPNAME = "SPNAM " + args.spname.strip()

if args.spname_file is None:
    spname_file = "DFS"
else:
    spname_file = args.spname_file
    
if args.begin is None:
    offset1 = float(GETPAR("CNST 1"))*1000
else:
    offset1 = float(args.begin)*1000.
    
if args.resolution is None:
    res = float(GETPAR("CNST 3"))/1000.
else:
    res = float(args.resolution)/1000.0

if args.end is None:
    offset2 = float(GETPAR("CNST 2"))*1000.
else:
    offset2 = float(args.end)*1000.
    
if args.pulse_name is None:
    pulse_name = "P 8"
else:
    pulse_name = "P " + args.pulse_name.strip()

if args.pulse_length is None:
    pl = float(GETPAR(pulse_name))/1e6
else:
    pl = float(args.pulse_length)/1e6
    
  
if args.wurst_cutoff is None:
    WURSTcutoff = 80.
else:
    WURSTcutoff = float(args.wurst_cutoff)


wavDir=os.getenv('XWINNMRHOME')+"/exp/stan/nmr/lists/wave/user/"
ShapeFull = wavDir+spname_file


def writeShape(amp,ph,filename,description,parameters):
    # print filename
    import datetime
    now=datetime.datetime.now()
    today=now.strftime("%d/%m/%Y")
    time=now.strftime("%H:%M:%S")
    npoints=len(amp)
    header="""##TITLE= %s
##USAGE= %s
##JCAMP-DX= 5.00 $$ Bruker JCAMP library
##DATA TYPE= Shape Data
##ORIGIN= Bruker Analytik GmbH
##DATE= %s
##TIME= %s
##$SHAPE_PARAMETERS= %s
##MINX= -1.000000e+02
##MAXX= 1.000000e+02
##MINY= 0.000000e+00
##MAXY= 0.000000e+00
##$SHAPE_EXMODE= None
##$SHAPE_TOTROT= 0.000000e+00
##$SHAPE_BWFAC= 0.000000e+00
##$SHAPE_INTEGFAC= 7.460936e-01
##$SHAPE_MODE= 1
##NPOINTS= %d
##XYPOINTS= (XY..XY)
""" % (filename,description,today,time,parameters,npoints)
    FILE=open(filename,'w')
    FILE.write(header)
    for i in range(npoints):
	FILE.write(', '.join((str(amp[i]),str(ph[i]))))
	FILE.write("\n")
    FILE.write("##END")
    FILE.close()

# number of points defining shape is depending on max offset and pulse length to avoid aliasing
minres=int(1e6/(2.0*offset1)/0.025)*0.025
NPT=int(round(pl*1e6/res))
if res>minres : 
    MSG(str(minres))
    MSG("%f ns (CNST3) resolution is not enough : min is %f ns" %(res*1000,minres*1000))
    sys.exit()

off1=offset1
off2=-offset1
span = offset1 - offset2

ph=[]
amp=[]

for i in range(NPT):
    t=(i+0.5)*pl/NPT
    # amplitude en pourcents
    amp1=50*(1-abs(sin(-pi/2+i*pi/NPT))**WURSTcutoff)
    # la phase en radians
    ph1=2*pi*(off1*t-span/2*t*t/pl)
    ph2=2*pi*(off2*t+span/2*t*t/pl)
    comp=amp1*cmath.exp(1j*ph1)+amp1*cmath.exp(1j*ph2)
    P=float(atan2(comp.imag, comp.real))/pi*180
    A=abs(comp)
    amp.append(A)
    ph.append(P)
    # print amp[i],ph[i]

descrip="double freq sweep-%d, sweep from %6.0f kHz to %6.0f kHz in %6.3f us" % (WURSTcutoff, offset1/1000, offset2/1000, pl*1e6)
parameters="Type : DFS ; cutoff:%d ; start sweep:%6.0f kHz ; end sw:%6.0f kHz ; pulse length:%6.3f us" % (WURSTcutoff, offset1/1000., offset2/1000., pl*1e6)
writeShape(amp,ph,ShapeFull,descrip,parameters) 
print("DFS shape written with %d points" % (NPT))
print(parameters)
PUTPAR(SPNAME,spname_file)
