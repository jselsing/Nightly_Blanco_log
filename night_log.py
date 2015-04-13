#!/usr/bin/env python
#######################
# Script to generate nightly log a Blanco 4m telescope, using DECam
########################


import sys, argparse, random, numpy as np
import json, math, os, subprocess, glob

########################

def split(tmp):

    length = len(tmp[0])
    tmp_split = []           
    for i in range(length):
        tmp_split.append([str(k[i]) for k in tmp])

    return tmp_split

########################

def float_convert(num):

    try:
        #num = float("{0:.3f}".format(float(num)))
	num = float("%.3f" % float(num))
    except:
	num = np.NaN

    return num




########################
def run_inventory(args):
    filename = args.date+"/inventory.txt"
    filt = '{if ($1 == "MJD") exit; if ($1 == "#expnum") found=1; if (!found || $8 != "obj") next; if (!found || $9 == "pointing") next; if (found) print}'



    first = ['~/'+args.propid+'/run_kentool.csh inventory '+args.day+'']
#    first = ['cat', 'inventory.txt']    

    cat = subprocess.Popen(first,
                            stdout=subprocess.PIPE,
                            shell=True)

    awk = subprocess.Popen(['awk', filt],
                            stdin=cat.stdout,
                            stdout=open(filename, "w"),
                            )

    from time import sleep
    sleep(10.05)

    data = np.genfromtxt(filename, dtype=None)

    return data

########################
def run_seeingall(args, expnum):
    filt = '{print $3}'

    seeing = []
    for i, k in enumerate(expnum):
        first = ['~/'+args.propid+'/run_kentool.csh seeingall '+str(k)+'']
        cat = subprocess.Popen(first,
                               stdout=subprocess.PIPE,
                               shell=True
                               )
    
        grep = subprocess.Popen(['grep', 'Avg Seeing'],
                               stdin=cat.stdout,
                               stdout=subprocess.PIPE
                               )
        awk = subprocess.Popen(['awk', filt],
                                stdin=grep.stdout,
                                stdout=subprocess.PIPE
                                )                          
    
        from time import sleep
        sleep(0.05)
        see, tmp = awk.communicate()
        seeing.append(see)


    return seeing  

########################
def read_from_fits(keyword, args, expnum, column):
    filt = '{print $'+str(column)+'}'
    

    keyword_val = []
    for i, k in enumerate(expnum):
        first = 'dfits /home4/images/fits/'+args.propid+'/DECam_00'+k+'.fits.fz'
        cat = subprocess.Popen(first, 
                               stdout=subprocess.PIPE,
                               shell=True
                               )
    
        grep = subprocess.Popen(['grep', keyword],
                               stdin=cat.stdout,
                               stdout=subprocess.PIPE
                               )
        awk = subprocess.Popen(['awk', filt],
                                stdin=grep.stdout,
                                stdout=subprocess.PIPE
                                )                          
    
        from time import sleep
        sleep(0.05)
        val, tmp = awk.communicate()
        keyword_val.append(val)


    return keyword_val  

########################

def run_psc(args, expnum):
    
    filt = '{if ($1 == "Clouds") found=1; if (found) print}'

    cloud = []
    sky = []
    teff = []
    for i, k in enumerate(expnum):
        first = ['~/'+args.propid+'/run_kentool.csh psc '+str(k)+'']
        cat = subprocess.Popen(first,
                               stdout=subprocess.PIPE,
                               shell=True
                               )
    
        awk = subprocess.Popen(['awk', filt],
                               stdin=cat.stdout,
                               stdout=subprocess.PIPE
                               )                
                    
        psc, tmp = awk.communicate()
        psc_coll = ' '.join(psc.split())
        cloud.append(psc_coll.split(' ')[2])
        sky.append(psc_coll.split(' ')[14])
        teff.append(psc_coll.split(' ')[16])


    return cloud, sky, teff



########################

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('--propid', type=str, default='2015A-0121', help = 'Your propsal id for the observing night')
    parser.add_argument('--date', type=str, default='20150319', help = 'The date of your observation, ie. 20150319')
    parser.add_argument('--day', type=str, default='', help = 'Day id to supply to kent tools, i.e -1')
    

    args = parser.parse_args(argv)

    if not args.propid:
        print('You need to supply a propid')
    
    if not args.date:
        print('You need to supply a date')



    print('Generating night log for:')
    print('PROPID: %s') % args.propid
    print('DATE: %s') % args.date
    print('Which is %s days ago') %args.day    


    
    print('Getting exposure ids belonging to propid')
    expid_tot = []
    for i in glob.glob('/home4/images/fits/'+args.propid+'/DECam_00*.fits.fz'):
        expid_tot.append(i[-14:-8])

        
    print('Running kentools: inventory')
    inv = run_inventory(args)    
    expnum, ra, dec, ut, filt, exp, secz, type_obs, object = split(inv)
    
    
    
    #Matching exposure numbers from running inventory against exposure numbers from proposal id 
    expnumber = [k for i, k in zip(expnum, expnum) if i in expid_tot]
    ra = [k for i, k in zip(expnum, ra) if i in expid_tot]
    dec = [k for i, k in zip(expnum, dec) if i in expid_tot]
    ut = [k for i, k in zip(expnum, ut) if i in expid_tot]
    filt = [k for i, k in zip(expnum, filt) if i in expid_tot]
    exp = [k for i, k in zip(expnum, exp) if i in expid_tot]
    secz = [k for i, k in zip(expnum, secz) if i in expid_tot]
    type_obs = [k for i, k in zip(expnum, type_obs) if i in expid_tot]
    object = [k for i, k in zip(expnum, object) if i in expid_tot]


    expnum = expnumber
    
    
    print('Reading header values:')
    lskypow = read_from_fits('LSKYPOW', args, expnum, 3)
    lskypow = np.array([float_convert(k[:-1]) for k in lskypow])
    lskyhot = read_from_fits('LSKYHOT', args, expnum, 3)
    lskyhot = np.array([float_convert(k[:-1]) for k in lskyhot])
    gskyhot = read_from_fits('GSKYHOT', args, expnum, 3)
    gskyhot = np.array([float_convert(k[:-1]) for k in gskyhot])
    msurtemp = read_from_fits('MSURTEMP', args, expnum, 2)
    msurtemp = np.array([float_convert(k[:-1]) for k in msurtemp])
    domelow = read_from_fits('DOMELOW', args, expnum, 3)
    domelow = np.array([float_convert(k[:-1]) for k in domelow])
    tempdiff = []
    for i,k in zip(domelow, msurtemp):
        try:
            tempdiff.append(float(i) - float(k) )
        except:
            print('Tempdiff is nan')
            tempdiff.append(np.NaN)          
    date_obs = read_from_fits('DATE-OBS', args, expnum, 2)
    date_obs = np.array([k[1:-18] for k in date_obs])    



    print('Reading dumped QR file')
    try:
        qr = np.genfromtxt(glob.glob(args.date+'/QR*')[0], dtype=None)
        expid, mjd, exptime, filter, nccds, fwhm, ellipticity, bkg, zd, ha = split(qr)
        expnumber = [i for i in expid if i in expnum]
        bkg = [float(k) for i, k in zip(expid, bkg) if i in expnum]
        fwhm = [float(k) for i, k in zip(expid, fwhm) if i in expnum]
        ellip = [float(k) for i, k in zip(expid, ellipticity) if i in expnum]
        zd = [float(k) for i, k in zip(expid, zd) if i in expnum]
    except:
        print('QR file should be dumped from terminal and placed in folder')
 


    print('Reading dumped qcInv file')
    try:
        qc = np.genfromtxt(glob.glob(args.date+'/*qcinv')[0], dtype=None, skip_footer = 1, delimiter='/t', missing_values='    ', filling_values=np.NaN)
        qc2 = []    
        for i in qc:
            list_qc = ' '.join(i.split()).split(' ')
            length = len(list_qc)
            if length == 12:
                qc2.append(list_qc)
            elif length == 8:
                list_missing = list_qc
                list_out = list_missing[:-1] + [np.NaN for i in range(12 - 8)] + list_missing[-1:]
                qc2.append(list_out)
        expid,  ra_tmp,  dec_tmp,  ut_tmp,  fil_tmp,  time_tmp, secz_tmp,  psf,  sky,  cloud,  teff, Object = split(qc2)
        expnumber = [i for i in expid if i in expnum]
        psf = [float(k) for i, k in zip(expid, psf) if i in expnum]
        sky = [float(k) for i, k in zip(expid, sky) if i in expnum]
        cloud = [float(k) for i, k in zip(expid, cloud) if i in expnum]
        teff = [float(k) for i, k in zip(expid, teff) if i in expnum]
    except:
        print('Run qcInv from godb and place file in folder')



    expnum = expnumber
    print('Generating log for exposure numbers:')    
    print expnum
    
    
    #    # In case you want to run seeingall and psc
    #    print('Running kentools: seeingall')
    #    seeall = run_seeingall(args, expnum)
    #
    #    seeing = np.array([float(k[:-1]) for k in seeall])
    #    
    #    #Saving to .txt file
    #    dt = [("seeing", type(seeing))]
    #    data = np.array(seeing, dtype=dt)
    #    np.savetxt(args.date+"/seeingall.txt", data, header="seeing", fmt = ['%s'] )      
    #
    #
    #    print('Running kentools: psc')
    #    cloud, sky, teff = run_psc(args, expnum)
    #
    #
    #    cloud = [float(k) for k in cloud]
    #    sky = [float(k) for k in sky]
    #    teff = [float(k) for k in teff]
    #    seeing = [float(k) for k in seeing]
          
    dt = [("expnum", type(expnum)), ("ra", type(ra)), ("dec", type(dec)), 
          ("filt", type(filt)), ("exp", type(exp)) , ("secz", type(secz)),
          ("object", type(object)), ("psf", type(psf)), ("zd", type(zd)),
          ("bkg", type(bkg)),
          ("fwhm", type(fwhm)), ("ellip", type(ellip)), ("cloud", type(cloud)),
          ("sky", type(sky)), ("teff", type(teff)), ("lskypow", type(lskypow)),
          ("lskyhot", type(lskyhot)), ("gskyhot", type(gskyhot)), ("tempdiff", type(tempdiff)),
          ("date_obs", type(date_obs)), ("ut", type(ut))]          
          
    data = np.array(zip(expnum, ra, dec, filt, exp, secz, object,
                        psf, zd, bkg, fwhm, ellip, cloud, sky, teff, lskypow,
                        lskyhot, gskyhot, tempdiff, date_obs, ut), dtype=dt)

    np.savetxt(args.date+"/summary.txt", data,
               header="expnum ra             dec          filt      exp      secz   object  seeing      zd    bkg     fwhm   ellip  cloud      sky   teff   lskypow lskyhot gskyhot tempdiff  date-obs       time",
               fmt = ['%s', '%s', '%s', '%s', '%4s', '%4s', '%s', '%.2f', '%.2f', '%.2f', '%.3f', '%.3f', '%.2f', '%.2f', '%.2f', '%4s', '%4s', '%4s', '%4s', '%s', '%s'],
               delimiter = '\t' )    
    

#########################


if __name__ == '__main__':

    main(argv = sys.argv[1:])

