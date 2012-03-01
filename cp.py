import scipy
import numpy
import os
import MOOGTools

def write_par_file(wl_start, wl_stop, stage_dir, b_dir, prefix, theta=-99.0):
    outfile_name = os.path.join(stage_dir,'Parfiles',b_dir,prefix+'.par')
    pf = open(outfile_name, 'w')

    stronglines = True

    pf.write('gridsyn\n')
    pf.write('terminal      \'x11\'\n')
    pf.write('summary_out   \'../../Output/'+b_dir+'/summary.out\'\n')
    pf.write('smoothed_out   \'../../Output/'+b_dir+'/junk.out\'\n')
    pf.write('standard_out  \'../../Output/'+b_dir+'/out1\'\n')
    pf.write('lines_in      \'../../Linelists/'+b_dir+'/'+prefix+'_weak_linelist.out\'\n')
    if stronglines:
        pf.write('stronglines_in      \'../../Linelists/'+b_dir+'/'+prefix+'_strong_linelist.out\'\n')
        pf.write('strong          1\n')
    pf.write('atmosphere      1\n')
    pf.write('molecules       2\n')
    pf.write('lines           1\n')
    pf.write('damping         0\n')
    pf.write('freeform        0\n')
    if (theta != -99.0):
        pf.write('flux/int        1\n')
    else:
        pf.write('flux/int        0\n')
    pf.write('plot            1\n')
    pf.write('synlimits\n')
    pf.write('                '+str(wl_start)+' '+str(wl_stop)+' 0.05 1.00\n')
    pf.write('plotpars        1\n')
    pf.write('                '+str(wl_start)+' '+str(wl_stop)+' 0.05 1.01\n')
    pf.write('                0.0 0.000 0.000 1.00\n')
    pf.write('                g 0.00 0.00 0.00 0.00 0.00\n')
    
    run_number = 1

    temps = range(2500, 4100, 100)+range(4250, 6250, 250)
    #temps = range(5250, 6250, 250)
    gravs = range(300, 550, 50)
    #temps = range(3800, 4100, 100)
    #gravs = range(300, 600, 100)
    #temps = [4000]
    #gravs = [400]

    for T in temps:
        #if T < 4000:
        #    gravs = range(300, 600, 50)
        #else:
        #    gravs = range(300, 550, 50)
        for G in gravs:
            for angle in theta:
                pf.write('RUN            '+str(run_number)+'\n')
                pf.write('stokes_out   \'../../Output/'+b_dir+'/'+prefix+'_THETA_'+str(angle)+'_MARCS_T'+str(T)+'G'+str(G)+'.moog\'\n')
                pf.write('mu              %10.4f\n' % numpy.cos(numpy.radians(angle)))
                pf.write('hardpost_out   \'../../Output/'+b_dir+'/dummy.ps\'\n')
                pf.write('model_in       \'../../Atmospheres/MARCS/MARCS_T'+str(T)+'_G'+str(G/100.0)+'_M0.0_t1.0.md\'\n')
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


def generateLineList(b_dir, prefix, wl_start, wl_stop, Bfield):
    ''' This subroutine generates a Zeeman-split MOOG line list '''
    mus = {74.5:1.3, 62.4:1.0895, 53.3:0.9303, 45.0:0.7854, 36.7:0.6405, 27.6:0.4813, 15.5:0.2706}  #degress:radians
    #Sets up the directories
    staging_dir = 'stage_'+prefix+'/'
    print "Feature : ", prefix, ' ', b_dir
    try:
        os.makedirs(os.path.join(staging_dir,'Parfiles',b_dir))
        os.makedirs(os.path.join(staging_dir,'Output',b_dir))
        os.makedirs(os.path.join(staging_dir,'Linelists', b_dir))
    except:
        pass
    write_par_file(wl_start, wl_stop, staging_dir, b_dir, prefix, mus)

    
    if (Bfield != -99):
        # These are the data files used for the atomic/molecular line data and QM values
        weak_file = '/home/deen/Code/python/MOOG/VALD_lines/atomic_corrected.zeeman'
        strong_file = '/home/deen/Code/python/MOOG/VALD_lines/strongLines.dat'
        molecules = '/home/deen/Code/python/MOOG/VALD_lines/molecular_corrected.dat'

        weakLines = []

        # Reads in the "weak" lines and determined the amount of Zeeman splitting
        for line in open(weak_file, 'r').readlines():
            l = line.split()
            if ( (float(l[0]) > wl_start) & (float(l[0]) < wl_stop) ):
                wave = float(l[0])
                species = float(l[1])
                exp_pot = float(l[2])
                loggf = float(l[3])
                Jlow = float(l[4])
                Jup = float(l[5])
                glow = float(l[6])
                gup = float(l[7])
                VdW = float(l[9])
                weakLines.append(MOOGTools.spectral_Line(wave, species, exp_pot, loggf, Jlow=Jlow, Jup=Jup, glow=glow, gup=gup,
                VdW=VdW))
                weakLines[-1].zeeman_splitting(Bfield)

        # Reads in the molecular lines.  Since Lande g factors for molecules are much weaker than atomic lines,
        # we do not determine Zeeman splitting to save computation time (and shorten the line lists)
        for line in open(molecules, 'r').readlines():
            l = line.split()
            if ( (float(l[0]) > wl_start) & (float(l[0]) < wl_stop) ):
                wave = float(l[0])
                species = float(l[1])
                exp_pot = float(l[2])
                loggf = float(l[3])
                DE = float(l[4])
                weakLines.append(MOOGTools.spectral_Line(wave, species, exp_pot, loggf, DissE=DE))

        # Sorts the weak lines by wavelength
        weakLines.sort()

        # Prepares to read in the strong lines
        strongLines = []

        # For each line in the strong line file, read it in (if within the wavelength bounds) and apply Zeeman splitting
        for line in open(strong_file, 'r').readlines():
            l = line.split()
            if ( (float(l[0]) > wl_start) & (float(l[0]) < wl_stop) ):
                wave = float(l[0])
                species = float(l[1])
                exp_pot = float(l[2])
                loggf = float(l[3])
                Jlow = float(l[4])
                Jup = float(l[5])
                glow = float(l[6])
                gup = float(l[7])
                VdW = float(l[9])
                strongLines.append(MOOGTools.spectral_Line(wave, species, exp_pot, loggf, Jlow=Jlow, Jup=Jup, glow=glow, gup=gup, VdW=VdW))
                strongLines[-1].zeeman_splitting(Bfield)

        # Sorts the strong lines by wavelength
        strongLines.sort()

        # Maximum number of strong lines (must match the constant in MOOG)
        strong_line_limit = 400

        zp = 'straight'

        sfn = os.path.join(staging_dir,'Linelists', b_dir, prefix+'_strong_linelist.out')
        outfile = open(sfn, 'w')

        for strongLine in strongLines:
            strongLine.dump(out=outfile, zeeman=zp, mode='MOOG')
        outfile.close()
        sort_file(sfn)

        wfn = os.path.join(staging_dir, 'Linelists', b_dir, prefix+'_weak_linelist.out')
        weakfile = open(wfn, 'w')
        for weakLine in weakLines:
            if not(weakLine in strongLines):
                weakLine.dump(out=weakfile, zeeman=zp, mode='MOOG')

        weakfile.close()
        sort_file(wfn)


Bfields = numpy.arange(0, 4.5, 0.5)
prefixes = ['alpha', 'bravo', 'charlie', 'delta', 'echo']
wl_starts = [1.15, 1.48, 2.17, 2.23, 2.22]
wl_stops = [1.22, 1.52, 2.22, 2.31, 2.23]

for B in Bfields:
    for feature in [zip(prefixes, wl_starts, wl_stops)[-1]]:
        prefix = feature[0]
        wl_start = feature[1]
        wl_stop = feature[2]
        B_dir = 'B_'+str(B)+'kG'
        generateLineList(B_dir, prefix, wl_start*10000.0, wl_stop*10000.0, B/10.0)
