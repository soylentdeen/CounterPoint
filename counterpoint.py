import scipy
import numpy
import os
import MoogTools
import AstroUtils
import sys


def sort_file(name, weak=False):
    data = open(name, 'r').readlines()
    wl = []
    for line in data:
        wl.append(float(line[0:10]))

    order = numpy.argsort(wl)
    out = open(name, 'w')
    if weak:
        out.write("Weak Lines\n")
    for i in order:
        out.write(data[i])
    out.close()


def generateLineList(b_dir, prefix, wl_start, wl_stop, Bfield, **kwargs):
    ''' This subroutine generates Zeeman-split MOOG line lists (both 
    strong and weak) '''

    staging_dir = 'stage_'+prefix+'/'
    print "Feature : ", prefix, ' ', b_dir

    try:
        os.makedirs(os.path.join(staging_dir,'Parfiles',b_dir))
        os.makedirs(os.path.join(staging_dir,'Output', b_dir))
        os.makedirs(os.path.join(staging_dir,'Linelists', b_dir))
    except:
        pass
    MoogTools.write_par_file(wl_start, wl_stop, staging_dir, b_dir, prefix,
            mode='gridstokes', strongLines=True, **kwargs)

    # Load in configuration file
    config = AstroUtils.parse_config('counterpoint.cfg')
    weak_file = config['weak_file']
    strong_file = config['strong_file']
    molecules = config['molecules']
    VALD_list = config['VALD_file']
    gf_corrections = config['gf_file']

    strongLines, weakLines = MoogTools.parse_VALD(VALD_list, strong_file,
            wl_start, wl_stop, Bfield, gf_corrections)

    #     CO
    weakLines = numpy.append(weakLines, MoogTools.parse_HITRAN(
             molecules+'05_HITEMP2010new.par', wl_start, wl_stop, Bfield,
             gf_corrections, weedout=2.5))
    #     OH
    weakLines = numpy.append(weakLines, MoogTools.parse_HITRAN(
             molecules+'13_HITEMP2010.par', wl_start, wl_stop, Bfield,
             gf_corrections))

    #     CN
    weakLines = numpy.append(weakLines, MoogTools.parse_Plez_CN(
             molecules+'CN_Plez_linelist.dat', wl_start, wl_stop, Bfield,
             gf_corrections))

    sfn = os.path.join(staging_dir, 'Linelists', b_dir,
            prefix+'_strong_linelist.stokes')
    outfile = open(sfn, 'w')

    for strongLine in strongLines:
        strongLine.dump(out=outfile, mode='MOOGSTOKES')
    outfile.close()
    sort_file(sfn)

    wfn = os.path.join(staging_dir, 'Linelists', b_dir,
            prefix+'_weak_linelist.stokes')
    weakfile = open(wfn, 'w')
    for weakLine in weakLines:
        weakLine.dump(out=weakfile, mode='MOOGSTOKES')

    weakfile.close()
    sort_file(wfn, weak=True)

    if Bfield==0.0:
        MoogTools.write_par_file(wl_start, wl_stop, staging_dir, b_dir,
                prefix, mode='gridsyn', strongLines=True, **kwargs)

        sfn = os.path.join(staging_dir, 'Linelists', b_dir,
                prefix+'_strong_linelist.scalar')
        outfile = open(sfn, 'w')

        for strongLine in strongLines:
            strongLine.dump(out=outfile, mode='MOOGSCALAR')
        outfile.close()
        sort_file(sfn)

        wfn = os.path.join(staging_dir, 'Linelists', b_dir,
                prefix+'_weak_linelist.scalar')
        outfile = open(wfn, 'w')

        for weakLine in weakLines:
            weakLine.dump(out=outfile, mode='MOOGSCALAR')
        outfile.close()
        sort_file(wfn, weak=True)


parameters = AstroUtils.parse_config(sys.argv[-1])
if type(parameters["BFields"])==str:
    Bfields = numpy.array(parameters["BFields"].split(','), dtype = float).tolist()
else:
    Bfields = numpy.array(parameters["BFields"]).tolist().tolist()
prefixes = numpy.array(parameters["prefixes"].split(','), dtype=str)
if type(parameters["wl_starts"])== str:
    wl_starts = numpy.array(parameters["wl_starts"].split(','), dtype=float)
    wl_stops = numpy.array(parameters["wl_stops"].split(','), dtype=float)
else:
    wl_starts = numpy.array([parameters["wl_starts"]], dtype=float)
    wl_stops = numpy.array([parameters["wl_stops"]], dtype=float)
if type(parameters["temps"])==str:
    temps = numpy.array(parameters["temps"].split(','), dtype=int).tolist()
else:
    temps = numpy.array([parameters["temps"]]).tolist()
if type(parameters["gravs"])==str:
    gravs = numpy.array(parameters["gravs"].split(','), dtype=int).tolist()
else:
    gravs = numpy.array([parameters["gravs"]]).tolist()


for B in Bfields:
    for feature in zip(prefixes, wl_starts, wl_stops):
        prefix = feature[0].strip()
        wl_start = feature[1]
        wl_stop = feature[2]
        B_dir = 'B_'+str(B)+'kG'
        generateLineList(B_dir, prefix, wl_start*10000.0,
                wl_stop*10000.0, B/10.0, temps=temps, gravs=gravs)
