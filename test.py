#import my_vlbatasks
import Utilities, AIPS, AIPSData, AIPSTask
import ParallelTask
from my_vlbatasks import fring_parallel

AIPS.userno = 52

AIPS.proxies.append('http://franz@202.8.37.186:8000')
AIPS.proxies.append('http://franz@202.8.37.186:8000')
AIPS.proxies.append('http://franz@202.8.37.186:8000')
AIPS.proxies.append('http://franz@202.8.37.186:8000')
AIPS.proxies.append('http://franz@202.8.37.186:8000')
AIPS.proxies.append('http://franz@202.8.37.186:8000')
AIPS.proxies.append('http://franz@202.8.37.186:8000')
AIPS.proxies.append('http://franz@202.8.37.186:8000')

rdisk1=AIPS.AIPSDisk(AIPS.proxies[1],1)
rdisk2=AIPS.AIPSDisk(AIPS.proxies[2],1)
rdisk3=AIPS.AIPSDisk(AIPS.proxies[3],1)
rdisk4=AIPS.AIPSDisk(AIPS.proxies[4],1)
rdisk5=AIPS.AIPSDisk(AIPS.proxies[5],1)
rdisk6=AIPS.AIPSDisk(AIPS.proxies[6],1)
rdisk7=AIPS.AIPSDisk(AIPS.proxies[7],1)
rdisk8=AIPS.AIPSDisk(AIPS.proxies[8],1)

AIPS.disks.append(rdisk1)
AIPS.disks.append(rdisk2)
AIPS.disks.append(rdisk3)
AIPS.disks.append(rdisk4)
AIPS.disks.append(rdisk5)
AIPS.disks.append(rdisk6)
AIPS.disks.append(rdisk7)
AIPS.disks.append(rdisk8)

data1=AIPSData.AIPSUVData('bin14_a','UVDATA',len(AIPS.disks)-8,1,AIPS.userno)
data2=AIPSData.AIPSUVData('bin14_a','UVDATA',len(AIPS.disks)-7,1,AIPS.userno)
data3=AIPSData.AIPSUVData('bin14_a','UVDATA',len(AIPS.disks)-6,1,AIPS.userno)
data4=AIPSData.AIPSUVData('bin14_a','UVDATA',len(AIPS.disks)-5,1,AIPS.userno)
data5=AIPSData.AIPSUVData('bin14_a','UVDATA',len(AIPS.disks)-4,1,AIPS.userno)
data6=AIPSData.AIPSUVData('bin14_a','UVDATA',len(AIPS.disks)-3,1,AIPS.userno)
data7=AIPSData.AIPSUVData('bin14_a','UVDATA',len(AIPS.disks)-2,1,AIPS.userno)
data8=AIPSData.AIPSUVData('bin14_a','UVDATA',len(AIPS.disks)-1,1,AIPS.userno)

uvdata_list = [data1]
# it makes no difference whether all FRING instances work on data1 or on data[1-8] individually
# uvdata_list = [data1, data2, data3, data4, data5, data6, data7, data8]
scans={'B1957+20':[[0, 5, 35, 48, 0, 5, 49, 22],[0, 5, 49, 24, 0, 5, 53, 54],[0, 5, 56, 9, 0, 6, 7, 55], [0, 6, 7, 57, 0, 6, 14, 16], [0, 6, 16, 31, 0, 6, 26, 26], [0, 6, 26, 28, 0, 6, 34, 38], [0, 6, 36, 53, 0, 6, 43, 30], [0, 6, 43, 32, 0, 6, 54, 58]]}

fring_parallel(n_jobs=8,uvdata_list=uvdata_list,srcs_and_scans=scans,1, 4, refant=3,snr=4.0)


# just for your reference, below is the defintion of fring_parallel()
def fring_parallel(n_jobs, uvdata_list, srcs_and_scans, solintmins,
                   inttimesecs, refant, search=[0], flag=-1, clversion=-1,
                   bpver=-1, snr=5.0, sumrrll=False, sumifs=False, sumfreq=False):
    if not isinstance(dictionary_of_sources_and_scans,dict):
        raise TypeError('dictionary_of_sources_and_scans is not a dictionary.')
    for src in dictionary_of_sources_and_scans.keys():
        if not all([len(t)==8 for t in dictionary_of_sources_and_scans[src]]):
            raise ValueError('Some of your times are messed up.')
    
    queue = ParallelTask.ParallelQueue()
    go = AIPSTask('go')
    start_sn = maxtab(uvdata_list[0],'SN')
    snversion = start_sn
    jobcounter = 0
    for source in srcs_and_scans.keys():
        for timerang in srcs_and_scans[source]:
            jobcounter += 1
            snversion += 1
            uvdata = uvdata_list[0]
            if len(uvdata_list) == n_jobs:
                uvdata = uvdata_list[jobcounter-1]
            uvdata.clrstat()
            fring = AIPSTask('fring')
            fring.indata = uvdata
            fring.outdisk = fring.indisk
            fring.timerang[1:] = timerang
            fring.docal = -1
            if clversion > -1:
                fring.docal = 1
                fring.gainuse = clversion
            fring.doband = -1
            if bpver > -1:                
                fring.doband = 3
                fring.bpver = bpver
            fring.calsour[1] = source
            fring.flagver = flag
            fring.refant = refant
            fring.search[1:] = search
            fring.aparm[1:] = [2,0,0,0,0,0,snr,0,1,0]
            fring.dparm[1:] = [0,0,0,inttimesecs,0,0,0,0]
            if sumrrll:
                fring.aparm[3] = 1
            if sumfreq:
                fring.aparm[4] = 1
            if sumifs:
                fring.aparm[5] = 1
            fring.solint = solintmins
            fring.snver = snversion
            fring.doapply = -1
            queue.queue(fring)
            if jobcounter == n_jobs:
                queue.go()
                queue = ParallelTask.ParallelQueue()
                jobcounter = 0
    if not jobcounter == 0:
        queue.go()

    start_sn += 1
    last_sn = maxtab(uvdata_list[0],'SN')
    clcal(uvdata, opcode='MERG', interpol='', snver=start_sn, invers=last_sn, 
          inCL=start_sn, outCL=last_sn+1,
          refant=refant, sources=[''], calsour=[''], dobtween=-1)
    for tbl in range(start_sn, last_sn+1):
        uvdata.zap_table('SN',tbl)
    
    tacop(indataset=uvdata, tabletype='SN', inver=maxtab(uvdata,'SN'), outdataset=uvdata,
          outver=start_sn, ncount=1)
    uvdata.zap_table('SN', last_sn+1)
