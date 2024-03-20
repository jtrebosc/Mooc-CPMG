#print "execute ssDFS"
import os
from math import *
import cmath
import argparse
import sys


description = """Generates a shape pulse with spinning-sideband Double Frequency Sweep (ssDFS)
The sweep consists of symmetrical with respect to carier frequency with frequency sweep of width 
equal to rotation, centered at offset. Direction of sweep can be out high to low or low to high frequency.
"""

class dummy():
    """A class to mimick argparse object in case argparse is not available."""
    def __init__(self):
        self.spname = False # -s
        self.spname_file = False # -f
        self.rotation = False # -r
        self.offset = False # -o
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
    parser.add_argument('-r', '--rotation',
                help='Rotation speed in Hz. Default: value stored in CNST 31.')
    parser.add_argument('-o', '--offset',
                help='Offset of the center of the sweep. Default: value stored in CNST 1.')
    parser.add_argument('--pulse_name', '-pn', 
                help='Name number of the pulse. Default: 8 for P 8.')
    parser.add_argument('--pulse_length', '-pl', 
                help='Length of the shape. Default: value stored in pulse_name (e.g. P 8).')
    parser.add_argument('--sweep_direction', '-d', 
                help='Direction of the sweep : +1 low to high, -1: high to low. Default: -1.')
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
    spname_file = "ssDFS"
else:
    spname_file = args.spname_file
    
if args.offset is None:
    offset = float(GETPAR("CNST 1"))*1000
else:
    offset = float(args.offset)*1000.
    
if args.resolution is None:
    res = float(GETPAR("CNST 3"))/1000.
else:
    res = float(args.resolution)/1000.0

if args.rotation is None:
    nuR = float(GETPAR("CNST 31"))
else:
    nuR = float(args.rotation)
    
if args.pulse_name is None:
    pulse_name = "P 8"
else:
    pulse_name = "P " + args.pulse_name.strip()

if args.pulse_length is None:
    pl = float(GETPAR(pulse_name))/1e6
else:
    pl = float(args.pulse_length)/1e6
    
if args.sweep_direction is None:
    sweepDir = -1 # high to low frequency sweep
else:
    sweepDir = int(args.sweep_direction)
    
if args.wurst_cutoff is None:
    WURSTcutoff = float(GETPAR("CNST 2"))
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
minres=int(1e6/(2.0*offset)/0.025)*0.025
NPT=int(round(pl*1e6/res))
if res>minres : 
    MSG(str(minres))
    MSG("%f ns (CNST3) resolution is not enough : min is %f ns" %(res*1000,minres*1000))
    sys.exit()

# calculate the sweep from 0 to nuR
#NnuR=int(round(2*offset/nuR))
sweepRate=nuR/pl
span=nuR
off1=offset+span/2.0
off2=-offset-span/2.0

ph=[]
amp=[]

for i in range(NPT):
    t=(i+0.5)*pl/NPT
    # amplitude en pourcents
    amp1=50*(1-abs(sin(-pi/2+i*pi/NPT))**WURSTcutoff)
    # amp2=40*(1-abs(sin(-pi/2+i*pi/NPT))**WURSTcutoff)
    # la phase en radians
    ph1=2*pi*(off1*t-span/2*t*t/pl)
    ph2=2*pi*(off2*t+span/2*t*t/pl)
    """	ph1=2*pi*(sweepDir*0.5*sweepRate*((i-NPT/2.0)*pl/NPT)**2 + offset*(i-NPT/2.0)*pl/NPT)
    ph3=2*pi*(sweepDir*0.5*sweepRate*((i-NPT/2.0)*pl/NPT)**2 + offset*(i-NPT/2.0)*pl/NPT)
    # le meme balayage mais NnuR*nuR=int(2*offset/nuR)*nuR plus loin
    ph2=2*pi*(sweepDir*0.5*sweepRate*((i-NPT/2.0)*pl/NPT)**2 + (offset - NnuR*nuR)*(i-NPT/2.0)*pl/NPT)
    # faire la somme des deux shapes :
    #amp2=amp1
    #amp1=0
    """
    comp=amp1*cmath.exp(1j*ph1)+amp1*cmath.exp(1j*ph2)
    P=float(atan2(comp.imag, comp.real))/pi*180
    A=abs(comp)
    amp.append(A)
    ph.append(P)
    # print amp[i],ph[i]

descrip="ss double freq sweep-%d, sweep of -%6.0f Hz at offsets centered at %6.0f and %6.0f kHz in %6.3f us" % (WURSTcutoff,nuR,(off1-span/2)/1000.0,(off2+span/2)/1000.0,pl*1e6)
parameters="Type : ssDFS ; cutoff:%d ; sweep:%6.0f Hz ; offset:%6.0f kHz ; pulse length:%6.3f us" % (WURSTcutoff,sweepDir*nuR,offset/1000.,pl*1e6)
writeShape(amp,ph,ShapeFull,descrip,parameters)
print("ssDFS shape writen with %d points" % (NPT))
print(parameters)
PUTPAR(SPNAME,spname_file)
