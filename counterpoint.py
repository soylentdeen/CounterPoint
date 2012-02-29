import scipy
import numpy
import os
import MOOGTools

def write_par_file(wl_start, wl_stop, datadir, prefix, index, zeeman_prefix):
    outfile_name = os.path.join(datadir, zeeman_prefix)+'/'+
                   prefix+str(index)+'.par'
    pf = open(outfile_name, 'w')

    stronglines = True

    pf.write('gridsyn\n')
    pf.write('terminal      \'x11\'\n')
    pf.write('summary_out   \'./output/summary.out\'\n')
    pf.write('standard_out  \'./output/out1\'\n')
    pf.write('lines_in      \'./'+prefix+'weak_linelist_'+str(index)+'.out\'\n')
    if stronglines:
        pf.write('stronglines_in      \'./'+prefix+
                 'strong_linelist_'+str(index)+'.out\'\n')
        pf.write('strong          1\n')
    pf.write('atmosphere      1\n')
    pf.write('molecules       2\n')
    pf.write('lines           1\n')
    pf.write('damping         0\n')
    pf.write('freeform        0\n')
    pf.write('flux/int        0\n')
    pf.write('plot            1\n')
    pf.write('synlimits\n')
    pf.write('                '+str(wl_start)+' '+str(wl_stop)+' 0.05 1.00\n')
    pf.write('plotpars        1\n')
    pf.write('                '+str(wl_start)+' '+str(wl_stop)+' 0.02 1.01\n')
    pf.write('                0.0 0.000 0.000 1.00\n')
    pf.write('                g 0.250 0.00 0.00 0.00 0.00\n')
    
    run_number = 1

    temps = range(2500, 4100, 100)+range(4250, 6250, 250)
    gravs = range(300, 550, 50)
    
    for T in temps:
        #if T < 4000:
        #    gravs = range(300, 600, 50)
        #else:
        #    gravs = range(300, 550, 50)
        for G in gravs:
            pf.write('RUN            '+str(run_number)+'\n')
            pf.write('smoothed_out   \'./output/'+prefix+
                    'MARCS_T'+str(T)+'G'+str(G)+'_'+str(index)+'.moog\'\n')
            pf.write('hardpost_out   \'./output/dummy.ps\'\n')
            pf.write('model_in       \'../../../../atmospheres/MARCS/MARCS_T'+
                    str(T)+'_G'+str(G/100.0)+'_M0.0_t1.0.md\'\n')
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


def generateLineList(datadir, prefix, wl_start, wl_stop, Bfield):
    weak_file = '/home/deen/Code/python/StarFormation/MOOG/linelists/VALD_lines/atomic_corrected.zeeman'
    strong_file = '/home/deen/Code/python/StarFormation/MOOG/linelists/VALD_lines/strongLines.dat'
    molecules = '/home/deen/Code/python/StarFormation/MOOG/linelists/VALD_lines/molecular_corrected.dat'

    weakLines = []

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
            weakLines.append(MOOGTools.spectral_Line(wave, species, exp_pot,
                loggf, Jlow=Jlow, Jup=Jup, glow=glow, gup=gup, VdW=VdW))
            weakLines[-1].zeeman_splitting(Bfield)

    for line in open(molecules, 'r').readlines():
        l = line.split()
        if ( (float(l[0]) > wl_start) & (float(l[0]) < wl_stop) ):
            wave = float(l[0])
            species = float(l[1])
            exp_pot = float(l[2])
            loggf = float(l[3])
            DE = float(l[4])
            weakLines.append(MOOGTools.spectral_Line(wave, species, exp_pot,
                loggf, DissE=DE))

    weakLines.sort()

    strongLines = []

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
            strongLines.append(MOOGTools.spectral_Line(wave, species, exp_pot,
                loggf, Jlow=Jlow, Jup=Jup, glow=glow, gup=gup, VdW=VdW))
            strongLines[-1].zeeman_splitting(Bfield)

    strongLines.sort()

    strong_line_limit = 400

    # Longitudinal line list:
    zeeman_polarization = ['FULL']
    #zeeman_polarization = ['LONG', 'TRANS', 'MU_'+str(mu)]

    for zp in zeeman_polarization:
        try:
            os.makedirs(os.path.join(datadir,zp))
        except:
            pass
        try:
            os.makedirs(os.path.join(datadir,zp,'output'))
        except:
            pass

        i = 0
        j = 0
        strong_file_counter = 0
        sfn = os.path.join(datadir,zp)+'/'+prefix+'strong_linelist_'+
              str(strong_file_counter)+'.out'
        outfile = open(sfn, 'w')
        strong_line_counter = 0
        breakpoints = []
        breakpoints.append(weakLines[0].wl)

        while i < len(strongLines):
            if (strong_line_counter + 
                    len(strongLines[i].zeeman[zp][0])) < strong_line_limit:
                strongLines[i].dump(out=outfile, zeeman=zp, mode='MOOG')
                strong_line_counter += len(strongLines[i].zeeman[zp][0])
                i = i + 1
            else:
                breakpoints.append(numpy.mean([strongLines[i].wl,
                    strongLines[i-1].wl]))
                j = i
                while ( (strong_line_counter < strong_line_limit) &
                        (j < len(strongLines) ) ):
                    strongLines[j].dump(out=outfile, mode='MOOG')
                    strong_line_counter += 1
                    j += 1
                outfile.close()
                sort_file(sfn)

                wfn = os.path.join(datadir,zp)+'/'+prefix+'weak_linelist_'+
                          str(strong_file_counter)+'.out'
                weakfile = open(wfn, 'w')
                bm = scipy.where( (weakLines > breakpoints[-2]) &
                        (weakLines < breakpoints[-1]) )[0]
                if len(bm) == 0:
                    weakfile.write('%10.3f%10s%10.3f%10.3f\n' %
                            (numpy.mean([breakpoints[-2], breakpoints[-1]]),
                            26.0, 6.0, -8.0))
                for k in bm:
                    if not(weakLines[k] in strongLines):
                        weakLines[k].dump(out=weakfile, zeeman=zp, mode='MOOG')

                weakfile.close()
                sort_file(wfn)
            

                write_par_file(breakpoints[-2], breakpoints[-1], datadir,
                        prefix, strong_file_counter, zp)

                strong_file_counter += 1
                sfn=os.path.join(datadir,zp)+'/'+prefix+'strong_linelist_'+
                        str(strong_file_counter)+'.out'
                outfile = open(sfn, 'w')
                strongLines[i-1].dump(out=outfile, mode='MOOG')
                strong_line_counter = 1
 
        outfile.close()
        breakpoints.append(weakLines[-1].wl)
        write_par_file(breakpoints[-2], breakpoints[-1], datadir, prefix,
                strong_file_counter,zp)

        wfn = os.path.join(datadir,zp)+'/'+prefix+'weak_linelist_'+
                str(strong_file_counter)+'.out'
        weakfile = open(wfn, 'w')

        for k in range(len(weakLines)):
            if not(weakLines[k] in strongLines):
                weakLines[k].dump(out=weakfile, zeeman=zp, mode='MOOG')
        
        weakfile.close()
        sort_file(wfn)

#datadir = raw_input("Enter output directory :")
#prefix = raw_input("Enter prefix :")
#wl_start = float(raw_input("Enter starting wavelength: "))*10000.0
#wl_stop = float(raw_input("Enter stopping wavelength :"))*10000.0
#Bfield = float(raw_input("Enter Magnetic Field (in Gauss) :"))*1e-4

Bfields = numpy.arange(0, 4.5, 0.5)
prefixes = ['f1_', 'f2_', 'f3_', 'f4_']
wl_starts = [1.15, 1.48, 2.17, 2.24]
wl_stops = [1.22, 1.52, 2.23, 2.31]

for B in [0.0, 2.0]:
    for feature in zip(prefixes, wl_starts, wl_stops):
        prefix = feature[0]
        wl_start = feature[1]
        wl_stop = feature[2]
        datadir = 'test_B_'+str(B)+'kG'
        generateLineList(datadir, prefix, wl_start*10000.0,
                wl_stop*10000.0, B/10.0)
