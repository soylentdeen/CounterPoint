import scipy
import numpy
import os
import MoogTools
import AstroUtils
import sys

def write_par_file(wl_start, wl_stop, stage_dir, b_dir, prefix,
        temps=None, gravs=None, **kwargs):
    outfile_name = os.path.join(stage_dir,'Parfiles',b_dir,'batch.par')
    pf = open(outfile_name, 'w')

    stronglines = True

    pf.write('gridstokes\n')
    pf.write('terminal      \'x11\'\n')
    pf.write('summary_out   \'../../Output/'+b_dir+'/summary.out\'\n')
    pf.write('smoothed_out  \'../../Output/'+b_dir+'/junk.out\'\n')
    pf.write('standard_out  \'../../Output/'+b_dir+'/out1\'\n')
    pf.write('atmos_dir     \'/home/deen/Data/Atmospheres/MARCS/\'\n')
    pf.write('out_dir       \'../../Output/'+b_dir+'/\'\n')
    pf.write('lines_in      \'../../Linelists/'+b_dir+'/'+prefix
            +'_weak_linelist.out\'\n')
    if stronglines:
        pf.write('stronglines_in      \'../../Linelists/'+b_dir+'/'+prefix+
                 '_strong_linelist.out\'\n')
        pf.write('strong          1\n')
    pf.write('atmosphere      1\n')
    pf.write('molecules       2\n')
    pf.write('lines           1\n')
    pf.write('damping         0\n')
    pf.write('freeform        0\n')
    pf.write('ncells          695\n')
    pf.write('nrings          23\n')
    #if (theta != -99.0):
    #    pf.write('flux/int        1\n')
    #else:
    #    pf.write('flux/int        0\n')
    pf.write('plot            1\n')
    pf.write('synlimits\n')
    pf.write('                '+str(wl_start)+' '+str(wl_stop)+' 0.05 1.00\n')
    pf.write('plotpars        1\n')
    pf.write('                '+str(wl_start)+' '+str(wl_stop)+' 0.05 1.01\n')
    pf.write('                0.0 0.000 0.000 1.00\n')
    pf.write('                g 0.000 0.00 0.00 0.00 0.00\n')
    
    run_number = 1

    if (not temps):
        temps = range(2500, 4100, 100)+range(4250, 6250, 250)
    if (not gravs):
        gravs = range(300, 550, 50)
    
    for T in temps:
        for G in gravs:
            pf.write('RUN            '+str(run_number)+'\n')
            pf.write('stokes_out   \''+prefix+
                    '_MARCS_T'+str(T)+'G'+str(G)+'\'\n')
            pf.write('hardpost_out   \'../../Output/'+b_dir+'/dummy.ps\'\n')
            pf.write('model_in       \'MARCS_T'+
                    str(T)+'_G'+str(G/100.0)+'_M0.0_t2.0.md\'\n')
            pf.write('abundances     1  1\n')
            pf.write('    12      0.0\n')
            run_number += 1
    pf.close()

def sort_file(name):
    data = open(name, 'r').readlines()
    wl = []
    for line in data:
        wl.append(float(line[0:10]))

    order = numpy.argsort(wl)
    out = open(name, 'w')
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
    config = AstroUtils.parse_config('cp.cfg')
    weak_file = config['weak_file']
    strong_file = config['strong_file']
    molecules = config['molecules']
    VALD_list = config['VALD_file']

    strongLines, weakLines = MoogTools.parse_VALD(VALD_list, strong_file,
            molecules, wl_start, wl_stop, Bfield)

    '''
    if ( len(weakLines) == 0 ):
        wave = (wl_start+wl_stop)/2.0
        species = 26.0
        exp_pot = 9.99
        loggf = -9.99
        Jlow = 1.0
        Jup = 1.0
        glow = 0.1
        gup = 0.1
        VdW = 0.0
        weakLines.append(MOOGTools.spectral_Line(wave, species, exp_pot,
            loggf, Jlow=Jlow, Jup=Jup, glow=glow, gup=gup, VdW=VdW))
        weakLines[-1].zeeman_splitting(Bfield)
    '''

    sfn = os.path.join(staging_dir, 'Linelists', b_dir,
            prefix+'_strong_linelist.stokes')
    outfile = open(sfn, 'w')

    for strongLine in strongLines:
        strongLine.dump(out=outfile, mode='MOOG')
    outfile.close()
    sort_file(sfn)

    wfn = os.path.join(staging_dir, 'Linelists', b_dir,
            prefix+'_weak_linelist.stokes')
    weakfile = open(wfn, 'w')
    for weakLine in weakLines:
        weakLine.dump(out=weakfile, mode='MOOG')

    weakfile.close()
    sort_file(wfn)

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
        sort_file(sfn)


parameters = AstroUtils.parse_config(sys.argv[-1])
if type(parameters["BFields"])==str:
    Bfields = numpy.array(parameters["BFields"].split(','), dtype = float)
else:
    Bfields = numpy.array(parameters["BFields"])
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
