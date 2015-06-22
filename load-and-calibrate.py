import my_vlbatasks
from my_vlbatasks import *
import ParallelTask
import sys, os
import ConfigParser, ast
from datetime import datetime
from astropy.time import Time

if len(sys.argv) < 2:
    print "You need to supply a config file. Check out ~/manuals/configs/load-and-calibrate.cfg \n"
    sys.exit()

start = Time(datetime.now(), scale='utc').isot.replace(':','-')
config = ConfigParser.ConfigParser() 
config.readfp(open(sys.argv[1]))

# all confguration info is stored as dictionaries
basics = ast.literal_eval(config.get('basics', 'basics'))
tasks = ast.literal_eval(config.get('tasks', 'tasks'))
fitldd = ast.literal_eval(config.get('fitld', 'fitld'))
cal = ast.literal_eval(config.get('cal', 'cal'))

outdir = basics['outdir']
indir = basics['indir']


AIPS.userno =basics['userno']
queue = ParallelTask.ParallelQueue()
go = AIPSTask('go')
go.qcreate = 1

if not os.path.exists(outdir):
    os.system('mkdir -p ' + outdir)
log_file = outdir + '/' + basics['log_file_base'] + '.' + str(start) + '.aipslog'
step_file = log_file + '.steps'
step_log = open(step_file, 'a')
step_log.write("Started at " + str(start) + ' .\n')
AIPS.log = open(log_file, 'a')


uvdata = AIPSUVData(fitldd['uvname'],'UVDATA',1,1)
if uvdata.exists():
    max_sn = maxtab(uvdata,'SN')
    max_cl = maxtab(uvdata,'cl')
    max_bp = maxtab(uvdata,'bp')
    max_fg = maxtab(uvdata,'fg')
    if any([max_sn,max_bp,max_fg]) > 0 or max_cl > 1:
        if interaction.yesno("Delete existing calibration tables? " \
                             "(Can keep a backup in next line. No will abort pipeline)"):
            if interaction.yesno("Do you want to keep a backup of those tables?"):
                out_tasav=tasav(uvdata)
                step_log.write("tasav(%s)" %(uvdata))
                out_tasav_file=outdir + '/' + out_tasav.name + '.TASAV.' + \
                                str(start) + '.fits'
                fittp(out_tasav,out_tasav_file)
                step_log.write("fittp(%s, %s)" %(out_tasav,out_tasav_file))            
            uvdata.zap_table('SN',-1)
            uvdata.zap_table('BP',-1)
            uvdata.zap_table('FG',-1)
            for i in range(2,max_cl+1):
                uvdata.zap_table('CL',i)
            step_log.write('Zapped all SN[1,%s], BP[1,%s], FG[1,%s], CL[2:%s]' %(max_sn,max_bp,max_fg,max_cl))
        else:
            sys.exit()
        
if tasks['load_data']:
    if uvdata.exists():
        if interaction.yesno("Delete existing UV dataset " + uvdata.name + \
                             "? (No will abort pipeline)"):
            uvdata.clrstat()
            uvdata.zap()
        else:
            sys.exit()
    fitld(indir+'/'+fitldd['raw_data_prefix'], uvdata, ncount=fitldd['n_files'],
          doconcat=fitldd['doconcat'], clint=fitldd['clint'], sources=fitldd['sources'],
          digicor=fitldd['digicor'], wtthresh=fitldd['wtthresh'])
    step_log.write("fitld(" + indir+'/'+fitldd['raw_data_prefix'] +', ' + str(uvdata) +',' 
                   + 'ncount=' + str(fitldd['n_files']) +', doconcat=' + str(fitldd['doconcat']) +
                   ', clint=' + str(fitldd['clint']) + ', sources=' + str(fitldd['sources']) + 
                   ', digicor=' + str(fitldd['digicor']) + ', wtthresh=' + str(fitldd['wtthresh'])
                   + '\n')

if tasks['inspect_data']:
    if not uvdata.exists():
        raise ValueError('You have not even loaded the data yet...')
    listr(uvdata,outdir + '/' + uvdata.name + '.listr.txt')
    prtan(uvdata,outdir + '/' + uvdata.name + '.prtan.txt')
    if interaction.yesno('Just put out listr.txt and prtan.txt.' \
                         'Do you want to continue? (Will reload the config file. ' \
                         'No will abort here.)'):
        print 'Going on...'
        basics = ast.literal_eval(config.get('basics', 'basics'))
        tasks = ast.literal_eval(config.get('tasks', 'tasks'))
        fitldd = ast.literal_eval(config.get('fitld', 'fitld'))
        cal = ast.literal_eval(config.get('cal', 'cal'))

        outdir = basics['outdir']
        indir = basics['indir']
        if interaction.yesno('Do you want to reload the data?'\
                             '(No will work with what is here.)'):
            uvdata.clrstat()
            uvdata.zap()
            fitld(indir+'/'+fitldd['raw_data_prefix'], uvdata, ncount=fitldd['n_files'],
                  doconcat=fitldd['doconcat'], clint=fitldd['clint'], sources=fitldd['sources'],
                  digicor=fitldd['digicor'], wtthresh=fitldd['wtthresh'])
            step_log.write("fitld(" + indir+'/'+fitldd['raw_data_prefix'] +', ' + 
                           str(uvdata) +',' + 'ncount=' + str(fitldd['n_files']) +
                           ', doconcat=' + str(fitldd['doconcat']) +
                           ', clint=' + str(fitldd['clint']) + 
                           ', sources=' + str(fitldd['sources']) + 
                           ', digicor=' + str(fitldd['digicor']) + 
                           ', wtthresh=' + str(fitldd['wtthresh'])
                           + '\n')
    else:
        sys.exit()

if tasks['load_evn_tbls']:
    evn_tbls = AIPSUVData(fitldd['evn_tbls_name'],'TASAV',1,1)
    if evn_tbls.exists():
        if interaction.yesno("Delete existing EVN tables " + evn_tbls.name + \
                             "? (No will abort pipeline)"):
            evn_tbls.clrstat()
            evn_tbls.zap()
        else:
            sys.exit()
    fitld(indir+'/'+fitldd['evn_tbls_raw_prefix'], evn_tbls)
    step_log.write("fitld(" + indir+'/'+fitldd['evn_tbls_raw_prefix']
                   + ', ' + str(evn_tbls) + ')\n')

if tasks['interactive']:
    man_fring_src = raw_input("man_fring_src: ")
    man_fring_timer = raw_input("man_fring_timer (as d h m s d h m s): ")
    bp_cal = raw_input("bp_cal: ")
    phase_ref_srcs = raw_input("phase_ref_srcs (as src1,src2,src3,...): ")
    main_phase_ref_src = raw_input("main_phase_ref_src: ")
    target = raw_input("target: ")
    refant = raw_input("refant (needs to be the two-letter station code, e.g. wb): ")
    search = raw_input("Antennas to search for extensive search in fring (refant,ant1,...): ")
    n_tec = int(raw_input("number of tec-files for tecor: "))
else:
    man_fring_src = cal['man_fring']['man_fring_src']
    man_fring_timer = cal['man_fring']['man_fring_timer']
    bp_cal = cal['bpass']['bp_src']
    phase_ref_srcs = cal['fring']['phase_ref_srcs']
    main_phase_ref_src = cal['fring']['main_phase_ref_src']
    target = basics['target']
    refant = cal['refant']
    search = cal['search']
    n_tec = cal['a_priori']['n_tec_files']

if isinstance(refant, basestring):
    refant = uvdata.antennas.index(refant.upper()) + 1
if isinstance(search, basestring):
    search = [ant.strip() for ant in search.split(',')]
    search = [uvdata.antennas.index(ant.upper()) + 1 for ant in search]
    if not refant in search:
        search.insert(0,refant)
if isinstance(man_fring_timer, basestring):
    man_fring_timer = map(int,man_fring_timer.split())

man_fring_timer = fix_timerange(man_fring_timer,cal['man_fring']['solint'])
phase_ref_srcs = [s.strip() for s in phase_ref_srcs.split(',')]

# do apriori calibration: apply Tsys cal from JIVE, run TECOR, parallactic angle correction
evn_tbls = AIPSUVData(fitldd['evn_tbls_name'],'TASAV',1,1)
if not evn_tbls.exists():
    raise ValueError('EVN tables ' + uvdata.name + ' do not exists.')

tacop(evn_tbls, 'SN',1,uvdata,0)
step_log.write("Copy EVN Pipeline tables\n" +
               "TACOP(SN 1" + 'from ' + str(evn_tbls) + 'to ' + str(uvdata) + ')\n')
clcal(uvdata, opcode='', interpol='', snver=1, invers=1, inCL=1, outCL=2, refant=0)
step_log.write("Apply them\n" +
               'CLCAL(' + str(uvdata) + ", opcode='', interpol='', snver=1, invers=1," +
               " inCL=1, outCL=2, refant=0)\n")
pang(uvdata, inCL=2, outCL=2)
step_log.write("Parralacti angle correction:\n"+
               "pang(" + str(uvdata) +", inCL=2, outCL=2)\n")

yr, day = aipstime2doy(uvdata.header['date_obs'])
tecfile1 = 'CODG' + day + '0.' + yr[-2:] + 'I' 
for n in range(0, n_tec+1):
    tecfile = 'CODG' + str(int(day) + n) + '0.' + yr[-2:] + 'I'
    if not os.path.exists(indir + '/' + tecfile):
        print 'Downloding tecfile ' + tecfile
        os.system('cd ' + indir + '; wget ftp://ftp.unibe.ch/aiub/CODE/' + yr
                  + '/' + tecfile + '.Z; gunzip ' + tecfile + '.Z')

tecor(uvdata, tecfile=indir + '/' + tecfile1, inCL=2, outCL=3,
      nfiles=n_tec)
step_log.write("Ionospheric Corrections:\n"+
               "tecor(" + str(uvdata) + ', ' + tecfile1 +', inCL=2, outCL=3,' +
               "n_files=" + str(n_tec) + "\n")

# run manfring 
manfring(uvdata, calsrc=man_fring_src, 
         timerang=man_fring_timer,
         solintmins=cal['man_fring']['solint'], 
         inttimesecs=basics['inttimesecs'], 
         refant=refant, snversion=2, clversion = 3)

step_log.write('Ran man_fring\n'+
               "manfring(" +  str(uvdata) + ', calsrc=' + man_fring_src +
               "timerang=" + str(man_fring_timer) + ", solintmins=" + str(cal['man_fring']['solint']) +
               "inttimesecs=" + str(basics['inttimesecs']) + ', refant=' + str(refant) +
               ", snversion=2, clversion = 3)"
               + "\n")

### apply manfring solts to all sources
clcal(uvdata, opcode='', interpol='ambg', snver=2, invers=2, inCL=3, outCL=4, refant=refant,
      sources=[''], calsour=[man_fring_src])
step_log.write("Apply manfring cal\n" +
               'CLCAL(' + str(uvdata) + ", opcode='', interpol='ambg', snver=2, invers=2," +
               " inCL=3, outCL=4, refant=" +str(refant) + 
               "sources='', " + "calsour=" + man_fring_src + ")\n")

# run fring on bp_calibrator
fring(uvdata, calsrcs=bp_cal, solintmins=cal['fring']['solint'],
      inttimesecs=basics['inttimesecs'], refant=refant,snversion=3,
      clversion=4, doband=-1,snr=cal['fring']['snr_cutoff'],
      sumrrll=cal['fring']['sumrrll'],
      sumifs=cal['fring']['sumifs'], sumfreq=cal['fring']['sumfreq'])
step_log.write('Ran FRING on bp_calibrator:\n'
               'fring(' + str(uvdata) + ', calsrcs=' + bp_cal +
               ', solintmins=' + str(cal['fring']['solint']) +
               ', inttimesecs=' + str(basics['inttimesecs']) +
               ', refant=' + str(refant)+', snversion=3'+
               ', clversion=4' +
               ', doband=-1' + 
               ', bpver=0' + 
               ', sumrrll=' + str(cal['fring']['sumrrll']) +
               ', sumifs=' + str(cal['fring']['sumifs']) +
               ', sumfreq=' + str(cal['fring']['sumfreq']))

## apply fringslts to bp_calibrator
clcal(uvdata, opcode='', interpol='ambg', snver=3, invers=3, inCL=4, outCL=5, refant=refant,
      sources=[bp_cal], calsour=[bp_cal])
step_log.write("Apply fring-cal to bp_calibrator\n" +
               'CLCAL(' + str(uvdata) + ", opcode='', interpol='ambg'," + 
               'snver=3, invers=3,' +
               " inCL=4, outCL=5, refant=" +str(refant) + 
               "sources="  + bp_cal + ", calsour=" + bp_cal + ")\n")

# plot the fring solutions with possm
outplot = outdir + '/' + 'fring.' + bp_cal + '.possm.ps'
plot_possm(uvdata, outplot, sources=bp_cal, clversion=5)

# run bpass on bp_calibrator
bpass(uvdata, bp_cal, solint=-1,
      solty=cal['bpass']['solty'], clversion=5,
      refant=refant,amps_only=cal['bpass']['amps_only'],
      phase_only=cal['bpass']['phase_only'])
step_log.write('Running Bpass:\n' +
               'bpass(' + str(uvdata) + ', srcname=' + bp_cal +
               ', solty=' + cal['bpass']['solty'] + ', clversion=5, ' +
               'refant=' + str(refant) +
               ' ,amps_only=' + str(cal['bpass']['amps_only']) +
               ', phase_only=' + str(cal['bpass']['phase_only']))

# plot the BP-table with possm
outplot = outdir + '/' + 'BP-table.' + bp_cal + '.possm.ps'
plot_possm(uvdata, outplot, sources=bp_cal, bpversion=1, plotbptable=True)


# run fring on phase-ref sources with Bpass applied, phase-ref src can be target itself
# get dict of scans for each source first
scans = {}
for src in phase_ref_srcs:
    d = get_source_obs_times(uvdata,src)
    scans[src] = [d[key] for key in d.keys()]
fring_parallel(n_jobs=8, uvdata=uvdata, dictionary_of_sources_and_scans=scans,
               solintmins=cal['fring']['solint'],
               inttimesecs=basics['inttimesecs'], refant=refant, search=search, flag=-1,
               clversion=4,
               bpver=1, snr=cal['fring']['snr_cutoff'],
               sumrrll=cal['fring']['sumrrll'],
               sumifs=cal['fring']['sumifs'], sumfreq=cal['fring']['sumfreq'])
step_log.write('Ran fring_parallel:\n'
               'fring_parallel(n_jobs='+ str(n_jobs)+ ', ' + str(uvdata) + 
               ', solintmins=' + str(cal['fring']['solint']) +
               ', inttimesecs=' + str(basics['inttimesecs']) +
               ', refant=' + str(refant)+ ', search=' + str(search) + ', snversion=4'+
               ', clversion=4' +
               ', bpver=1' + ', flag=-1' + 
               ', sumrrll=' + str(cal['fring']['sumrrll']) +
               ', sumifs=' + str(cal['fring']['sumifs']) +
               ', sumfreq=' + str(cal['fring']['sumfreq']))

# apply fring fitting slts to fringe fitters themselves
for src in phase_ref_srcs:
    clcal(uvdata, opcode='', interpol='ambg', snver=4, invers=4, inCL=4, outCL=5, refant=refant,
          sources=[src], calsour=[src])
    step_log.write("Apply fring-cal to phase_ref_srcs\n" +
                   'CLCAL(' + str(uvdata) + ", opcode='', interpol='ambg', snver=4, invers=4," +
                   " inCL=4, outCL=5, refant=" +str(refant) + 
                   "sources="  + src + ", calsour=" + src + ")\n")
    # plot the fring solutions with possm
    outplot = outdir + '/' + 'fring+bpass.' + src + '.possm.ps'
    plot_possm(uvdata, outplot, sources=src, clversion=5, bpversion=1)

# and then the solts from the phase_ref_src to the target (if target not in phase_ref_srcs)
if target not in phase_ref_srcs:
    clcal(uvdata, opcode='', interpol='ambg', snver=4, invers=4, inCL=4, outCL=5, 
          refant=refant,
          sources=[target], calsour=main_phase_ref_src)
    step_log.write("Apply fring-cal to target\n" +
                   'CLCAL(' + str(uvdata) + ", opcode='', interpol='ambg', snver=3," + 
                   'invers=3,' +
                   " inCL=4, outCL=5, refant=" +str(refant) + 
                   "sources="  + target + ", calsour=" + 
                   main_phase_ref_src + ")\n")
    

stop = Time(datetime.now(), scale='utc').isot.replace(':','-')
step_log.write('Finished at ' + stop)
step_log.close()
