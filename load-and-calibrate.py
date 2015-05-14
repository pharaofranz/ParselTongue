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

start = Time(datetime.now(), scale='utc').isot
config = ConfigParser.ConfigParser() 
config.readfp(open(sys.argv[1]))

# all confguration info is stored as dictionaries
basics = ast.literal_eval(config.get('basics', 'basics'))
tasks = ast.literal_eval(config.get('tasks', 'tasks'))
fitldd = ast.literal_eval(config.get('fitld', 'fitld'))
cal = ast.literal_eval(config.get('cal', 'cal'))

AIPS.userno =basics['userno']
queue = ParallelTask.ParallelQueue()
go = AIPSTask('go')
go.qcreate = 1

if not os.path.exists(basics['outdir']):
    os.system('mkdir -p ' + basics['outdir'])
log_file = basics['outdir'] + '/' + basics['log_file_base'] + start + '.aipslog'
step_file = log_file + '.steps'
step_log = open(step_file, 'a')
step_log.write("Started at " + start + ' .\n')
AIPS.log = open(log_file, 'a')

uvdata = AIPSUVData(fitldd['uvname'],'UVDATA',1,1)
if tasks['load_data']:
    if uvdata.exists():
        if interaction.yesno("Delete existing UV dataset " + uvdata.name + \
                             "? (No will abort pipeline)"):
            uvdata.clrstat()
            uvdata.zap()
        else:
            sys.exit()
    fitld(basics['indir']+'/'+fitldd['raw_data_prefix'], uvdata, ncount=fitldd['n_files'],
          doconcat=fitldd['doconcat'], clint=fitldd['clint'], sources=fitldd['sources'],
          digicor=fitldd['digicor'], wtthresh=fitldd['wtthresh'])
    step_log.write("fitld(" + basics['indir']+'/'+fitldd['raw_data_prefix'] +', ' + str(uvdata) +',' 
                   + 'ncount=' + str(fitldd['n_files']) +', doconcat=' + str(fitldd['doconcat']) +
                   ', clint=' + str(fitldd['clint']) + ', sources=' + str(fitldd['sources']) + 
                   ', digicor=' + str(fitldd['digicor']) + ', wtthresh=' + str(fitldd['wtthresh'])
                   + '\n')

if tasks['load_evn_tbls']:
    evn_tbls = AIPSUVData(fitldd['evn_tbls_name'],'TASAV',1,1)
    if evn_tbls.exists():
        if interaction.yesno("Delete existing EVN tables " + evn_tbls.name + \
                             "? (No will abort pipeline)"):
            evn_tbls.clrstat()
            evn_tbls.zap()
        else:
            sys.exit()
    fitld(basics['indir']+'/'+fitldd['evn_tbls_raw_prefix'], evn_tbls)
    step_log.write("fitld(" + basics['indir']+'/'+fitldd['evn_tbls_raw_prefix']
                   + ', ' + str(evn_tbls) + ')\n')
    
if tasks['do_apriori']:
    evn_tbls = AIPSUVData(fitldd['evn_tbls_name'],'TASAV',1,1)
    tecfile = cal['a_priori']['tecfile']
    if not evn_tbls.exists():
        raise ValueError('EVN tables ' + uvdata.name + ' do not exists.')
    if not os.path.exists(basics['indir'] + '/' + tecfile):
        yr, day = aipstime2doy(uvdata.header['date_obs'])
        tecfile = 'CODG' + day + '0.' + yr[-2:] + 'I'
        print 'Downloding tecfile ' + tecfile
        os.system('cd ' + basics['indir'] + '; wget ftp://ftp.unibe.ch/aiub/CODE/' + yr
                  + '/' + tecfile + '.Z; gunzip ' + tecfile + '.Z')
    uvdata.zap_table('SN', -1)
    for i in range(2, maxtab('CL',uvdata)+1):
        uvdata.zap_table('CL',i)
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
    tecor(uvdata, tecfile=basics['indir'] + '/' + tecfile, inCL=2, outCL=3,
          n_files=cal['a_priori']['n_files'])
    step_log.write("Ionospheric Corrections:\n"+
                   "tecor(" + str(uvdata) + ', ' + tecfile +', inCL=2, outCL=3,' +
                   "n_files=" + str(cal['a_priori']['n_tec_files']) + "\n")

if tasks['do_man_fring']:
    manfring(uvdata, calsrc=cal['man_fring']['man_fring_src'], 
             timerang=cal['man_fring']['man_fring_time'], 
             solintmins=cal['man_fring']['solint'], 
             inttimesecs=basics['inttimesecs'], 
             refant=cal['refant'], snversion=2, clversion = 3)
    step_log.write('Ran man_fring\n'+
                   "manfring(" +  str(uvdata) + ', calsrc=' + cal['man_fring']['man_fring_src'] +
                   "timerang=" + str(timerang) + ", solintmins=" + str(cal['man_fring']['solint']) +
                   "inttimesecs=" + str(basics['inttimesecs']) + ', refant=' + str(cal['refant']) +
                   ", snversion=2, clversion = 3)"
                   + "\n")
    clcal(uvdata, opcode='', interpol='amb', snver=2, invers=2, inCL=3, outCL=4, refant=cal['refant'],
          sources=[basics['target']], calsour=[cal['man_fring_src']])
    step_log.write("Apply manfring cal\n" +
                   'CLCAL(' + str(uvdata) + ", opcode='', interpol='amb', snver=2, invers=2," +
                   " inCL=3, outCL=4, refant=" +str(cal['refant']) + 
                   "sources=" str(basics['target']) + "calsour=" + str(cal['man_fring_src']) + ")\n")

if tasks['do_bpass']:
    bpass(uvdata, cal['bpass']['bp_src'], solint=cal['bpass']['solint'], 
          solty=cal['bpass']['solty'], clversion=-1,
          refant=basics['refant'],amps_only=cal['bpass']['amps_only'],
          phase_only=cal['bpass']['phase_only'])
    step_log.write('Running Bpass:\n' +
                   'bpass(' + str(uvdata) + ', srcname=' + cal['bpass']['bp_src'] +
                   ', solint=' + str(cal['bpass']['solint']) + 
                   ', solty=' + cal['bpass']['solty'] + ', clversion=-1, ' +
                   'refant=' + str(basics['refant']) +
                   ' ,amps_only=' + str(cal['bpass']['amps_only']) +
                   ', phase_only=' + str(cal['bpass']['phase_only']))

if tasks['do_fring']:
    fring(uvdata, calsrc=cal['fring']['fring_src'], solintmins=cal['fring']['solint'],
          inttimesecs=basics['inttimesecs'], refant=basics['refant'],
          clversion=cal['fring']['cl_table'], doband=cal['fring']['bo_band'],
          bpver=cal['fring']['bp_tble'], sumrrll=cal['fring']['sumrrll'],
          sumifs=cal['fring']['sumifs'], sumfreq=cal['fring']['sumfreq'])
    step_log.write('Running FRING:\n'
                   'fring(' + str(uvdata) + ', calsrc=' + cal['fring']['fring_src']+
                   ', solintmins=' + str(cal['fring']['solint']) +
                   ', inttimesecs=' + str(basics['inttimesecs']) +
                   ', refant=' + str(basics['refant'])+
                   ', clversion=' + str(cal['fring']['cl_table'])+
                   ', doband=' + str(cal['fring']['bo_band'])+
                   ', bpver=' + str(cal['fring']['bp_tble'])+
                   ', sumrrll=' + str(cal['fring']['sumrrll']) +
                   ', sumifs=' + str(cal['fring']['sumifs']) +
                   ', sumfreq=' + str(cal['fring']['sumfreq']))
)
if tasks['apply_fring']:
    caldata = AIPSUVData('bin17','UVDATA',1,1)
    clcal(caldata, opcode='cali', interpol='self', snver=2,invers=2, inCL=3,outCL=4, refant=3,\
          sources=['B1957+20'], calsour=['B1957+20'])

stop = Time(datetime.now(), scale='utc').isot
step_log.write('Finished at ' + stop)
step_log.close()
