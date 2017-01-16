#
# CHARIS.py -- CHARIS personality for g2cam instrument interface
#
# Eric Jeschke (eric@naoj.org)
#
"""This file implements a simulator for a simulated instrument (CHARIS).
"""
import logging
import math
import sys, os, time
import re
import threading
from datetime import datetime, timedelta

import numpy
import select
import SocketServer
import subprocess

import astropy.io.fits as pyfits
from astropy.time import Time

# gen2 base imports
from g2base import Bunch, Task

# g2cam imports
from g2cam.Instrument import BASECAM, CamError, CamCommandError
from g2cam.util import common_task


# Value to return for executing unimplemented command.
# 0: OK, non-zero: error
unimplemented_res = 0

class HeaderServer(SocketServer.TCPServer):
    timeout = 1.0

    def shutdown(self):
        self.__shutdown_request = True
        
    def XXhandle_timeout(self):
        if self.__shutdown_request:
            return False
        return True

    def XXhandle_request(self):
        """Handle one request, possibly blocking.

        Respects self.timeout.
        """

        fd_sets = SocketServer._eintr_retry(select.select, [self], [], [], self.timeout)
        if not fd_sets[0]:
            stillOK = self.handle_timeout()
            if not stillOK:
                return
        self._handle_request_noblock()

class HeaderQueryHandler(SocketServer.BaseRequestHandler):
    def handle(self):
        req = self.request.recv(1024).strip()
        try:
            hdr = self.server.localFunc(req)
        except Exception as e:
            logging.warn('failed to running command %s: %s' % (req, e))
            hdr = pyfits.Header()
        self.request.sendall(hdr.tostring())
        
class HeaderTask(Task.Task):
    def __init__(self, port, localFunc, timeout=0.5):
        self.port = port
        self.localFunc = localFunc
        self.timeout = timeout
        self.handlerClass = HeaderQueryHandler
        
        super(HeaderTask, self).__init__()
        
    def stop(self):
        self.ev_quit.set()

    def execute(self):
        self.logger.info('Starting waiter')

        server = HeaderServer(('', self.port), self.handlerClass)
        server.timeout = self.timeout
        server.localFunc = self.localFunc
        
        while not self.ev_quit.isSet():
            try:
                server.handle_request()
            except Exception, e:
                self.logger.error("Error invoking fn: %s" % str(e))

        self.logger.info('Stopping periodic interval task')


class CHARISError(CamCommandError):
    pass

class CHARIS(BASECAM):

    def __init__(self, logger, env, ev_quit=None):

        super(CHARIS, self).__init__()

        self.logger = logger
        self.env = env
        # Convoluted but sure way of getting this module's directory
        self.mydir = os.path.split(sys.modules[__name__].__file__)[0]

        if not ev_quit:
            self.ev_quit = threading.Event()
        else:
            self.ev_quit = ev_quit

        # Holds our link to OCS delegate object
        self.ocs = None

        # We define our own modes that we report through status
        # to the OCS
        self.mode = 'default'

        # Thread-safe bunch for storing parameters read/written
        # by threads executing in this object
        self.param = Bunch.threadSafeBunch()

        # Interval between status packets (secs)
        self.param.status_interval = 10.0


    #######################################
    # INITIALIZATION
    #######################################

    def read_header_list(self, inlist):
        # open input list and check if it exists
        try:
            fin = open(inlist)
        except IOError:
            raise CHARISError("cannot open %s" % str(inlist))

        StatAlias_list = []
        FitsKey_list = []
        FitsType_list = []
        FitsDefault_list = []
        FitsComment_list = []
        header_num = 0
        for line in fin:
            if not line.startswith('#'):
                param = re.split('[\s!\n\t]+',line[:-1])
                StatAlias = param[0]
                FitsKey = param[1]
                FitsType = param[2]
                if FitsType == 'string':
                    FitsDefault = param[3]
                elif FitsType == 'float':
                    FitsDefault = float(param[3])
                elif FitsType == 'int':
                    FitsDefault = long(param[3])

                FitsComment = ""
                for i in range(len(param)):
                    if i >= 4:
                        if i == 4:
                            FitsComment = param[i]
                        else:
                            FitsComment = FitsComment + ' ' + param[i]
                
                StatAlias_list.append(StatAlias)
                FitsKey_list.append(FitsKey)
                FitsType_list.append(FitsType)
                FitsDefault_list.append(FitsDefault)
                FitsComment_list.append(FitsComment)
                header_num += 1 

        header = zip(StatAlias_list, FitsKey_list, FitsType_list, FitsDefault_list, FitsComment_list)

        fin.close()

        return header

    def init_stat_dict(self, header):
        statusDict = {}
        for i in range(len(header)):
            if header[i][0] != 'NA':
                statusDict[header[i][0]] = header[i][3]
        return statusDict

    def initialize(self, ocsint):
        '''Initialize instrument.
        '''
        super(CHARIS, self).initialize(ocsint)
        self.logger.info('***** INITIALIZE CALLED *****')
        # Grab my handle to the OCS interface.
        self.ocs = ocsint

        # Get instrument configuration info
        self.obcpnum = self.ocs.get_obcpnum()
        self.insconfig = self.ocs.get_INSconfig()

        # Thread pool for autonomous tasks
        self.threadPool = self.ocs.threadPool

        # For task inheritance:
        self.tag = 'charis'
        self.shares = ['logger', 'ev_quit', 'threadPool']

        # Get our 3 letter instrument code and full instrument name
        self.inscode = self.insconfig.getCodeByNumber(self.obcpnum)
        self.insname = self.insconfig.getNameByNumber(self.obcpnum)

        # Figure out our status table name.
        if self.obcpnum == 9:
            # Special case for SUKA.  Grrrrr!
            tblName1 = 'OBCPD'
        else:
            tblName1 = ('%3.3sS%04.4d' % (self.inscode, 1))

        self.stattbl1 = self.ocs.addStatusTable(tblName1,
                                                ['status', 'mode', 'count',
                                                 'time'])

        # read telescope status
        self.tel_header = self.read_header_list("header_telescope_20160917.txt")
        self.statusDictTel = self.init_stat_dict(self.tel_header)
        
        # AO188 status dictionary
        self.ao_header = self.read_header_list("header_ao188+lgs_obs_20110425.txt")
        self.statusDictAO = self.init_stat_dict(self.ao_header)
        self.ao_stat = 1
        
        # Add other tables here if you have more than one table...

        # Establish initial status values
        self.stattbl1.setvals(status='ALIVE', mode='LOCAL', count=0)

        # Handles to periodic tasks
        self.status_task = None
        self.power_task = None

        # Lock for handling mutual exclusion
        self.lock = threading.RLock()


    def start(self, wait=True):
        super(CHARIS, self).start(wait=wait)

        self.logger.info('CHARIS STARTED.')

        # Start auto-generation of status task
        t = common_task.IntervalTask(self.putstatus,
                                     self.param.status_interval)
        self.status_task = t
        t.init_and_start(self)

        # Start task to monitor summit power.  Call self.power_off
        # when we've been running on UPS power for 60 seconds
        t = common_task.PowerMonTask(self, self.power_off, upstime=60.0)
        #self.power_task = t
        #t.init_and_start(self)

        # Start header generating task
        self.logger.info('HeaderTask: %s' % (HeaderTask))        
        t = HeaderTask(6666, self.localCmd)
        t.init_and_start(self)

    def localCmd(self, cmdStr):
        args = cmdStr.split()
        cmd = args[0]
        cmdargs = args[1:]
        self.logger.debug('cmd=%s arg=%s' % (cmd, args))
        
        if cmd == 'hdr':
            return self.return_new_header(*cmdargs)
        elif cmd == 'seqno':
            return self.reqframes(*cmdargs)
    
    def stop(self, wait=True):
        super(CHARIS, self).stop(wait=wait)

        # Terminate status generation task
        if self.status_task is not None:
            self.status_task.stop()

        self.status_task = None

        # Terminate power check task
        if self.power_task is not None:
            self.power_task.stop()

        self.power_task = None

        self.logger.info("CHARIS STOPPED.")


    #######################################
    # INTERNAL METHODS
    #######################################

    def execCmd(self, cmdStr, subtag=None, callback=None):
        self.logger.debug('execIng: %s', cmdStr)
        proc = subprocess.Popen([cmdStr], shell=True, bufsize=1,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ret = []
        while True:
            l = proc.stdout.readline()
            if not l or len(ret) > 1000:
                break
            ret.append(l)
            if callback is not None:
                callback(subtag, l)
                
            if False and subtag is not None:
                self.ocs.setvals(subtag, cmd_str='dribble: %s' % (l))

            self.logger.debug('exec ret: %s', l)
            
        err = proc.stderr.read()
        self.logger.warn('exec stderr: %s', err)

        return ret

    def execOneCmd(self, actor, cmdStr, timelim=60.0, subtag=None, callback=None):
        """ Execute a command using the actorcore oneCmd wrapper.
        """

        return self.execCmd('oneCmd.py %s --level=i --timelim=%0.1f %s' % (actor, timelim, cmdStr),
                            subtag=subtag, callback=callback)
        
    def dispatchCommand(self, tag, cmdName, args, kwdargs):
        self.logger.debug("tag=%s cmdName=%s args=%s kwdargs=%s" % (
            tag, cmdName, str(args), str(kwdargs)))

        params = {}
        params.update(kwdargs)
        params['tag'] = tag

        try:
            # Try to look up the named method
            method = getattr(self, cmdName)

        except AttributeError, e:
            result = "ERROR: No such method in subsystem: %s" % (cmdName)
            self.logger.error(result)
            raise CamCommandError(result)

        return method(*args, **params)

    def update_header_stat(self):
        """ Update the external data feeding our headers. """

        self.logger.info('updating telescope info')
        self.ocs.requestOCSstatus(self.statusDictTel)
        if self.ao_stat == 1:
            self.logger.info('updating AO188 info')
            self.ocs.requestOCSstatus(self.statusDictAO)

        self.logger.info('fetching header...')
        hdr = self.fetch_header('here', 9999, 1, 4.5, 'planetX', datetime.utcnow())
        phdu = pyfits.PrimaryHDU(header=hdr)
        
        self.logger.info('writing file...')
        phdu.writeto('/tmp/foo_2.fits', clobber=True, checksum=True)
        
    def return_new_header(self, frameid, mode, itime, fullHeader=True):
        """ Update the external data feeding our headers and generate one. """

        self.logger.info('updating telescope info')
        self.ocs.requestOCSstatus(self.statusDictTel)
        if self.ao_stat == 1:
            self.logger.info('updating AO188 info')
            self.ocs.requestOCSstatus(self.statusDictAO)

        self.logger.info('fetching header...')
        hdr = self.fetch_header('here', frameid, mode, float(itime),
                                'planetX', datetime.utcnow(),
                                fullHeader=fullHeader)
        
        return hdr
        
    def fetch_header(self, path, frameid, mode, itime, target, utc_start,
                     fullHeader=True):

        hdr = pyfits.Header()

        # Date and Time
        utc_end = utc_start + timedelta(seconds=float(itime))
        hst_start = utc_start - timedelta(hours=10)
        hst_end = hst_start + timedelta(seconds=float(itime))

        date_obs_str = utc_start.strftime('%Y-%m-%d')
        utc_start_str = utc_start.strftime('%H:%M:%S.%f')[:-3]
        utc_end_str = utc_end.strftime('%H:%M:%S.%f')[:-3]
        hst_start_str = hst_start.strftime('%H:%M:%S.%f')[:-3]
        hst_end_str = hst_end.strftime('%H:%M:%S.%f')[:-3]

        hdr.set('DATE-OBS',date_obs_str, "Observation start date (yyyy-mm-dd)")
        hdr.set('UT', utc_start_str, "HH:MM:SS.SS typical UTC at exposure")
        hdr.set('UT-STR', utc_start_str, "HH:MM:SS.SS UTC at exposure start")
        hdr.set('UT-END', utc_end_str, "HH:MM:SS.SS UTC at exposure end")
        hdr.set('HST', hst_start_str, "HH:MM:SS.SS typical HST at exposure")
        hdr.set('HST-STR', hst_start_str, "HH:MM:SS.SS HST at exposure start")
        hdr.set('HST-END', hst_end_str, "HH:MM:SS.SS HST at exposure end")

        # calculate MJD
        t_start = Time(utc_start, scale='utc')
        t_end = Time(utc_end, scale='utc')
        hdr.set('MJD',t_start.mjd, "Modified Julian Day at typical time")
        hdr.set('MJD-STR',t_start.mjd, "Modified Julian Day at exposure start")
        hdr.set('MJD-END',t_end.mjd, "Modified Julian Day at exposure end")

        # Local sidereal time
        # longitude = 155.4761
        # gmst = utcDatetime2gmst(utc_start)
        # lst = gmst2lst(longitude, hour=gmst.hour, minute=gmst.minute, second=(gmst.second + gmst.microsecond/10**6),
        #                lonDirection="W", lonUnits="Degrees")
        # lst_fmt = '%02d:%02d:%06.3f' % (lst[0],lst[1],lst[2])
        # hdr.set('LST',lst_fmt, "HH:MM:SS.SS typical LST at exposure")
        # hdr.set('LST-STR',lst_fmt, "HH:MM:SS.SS LST at exposure start")

        # gmst = utcDatetime2gmst(utc_end)
        # lst = gmst2lst(longitude, hour=gmst.hour, minute=gmst.minute, second=(gmst.second + gmst.microsecond/10**6),
        #                lonDirection="W", lonUnits="Degrees")
        # lst_fmt = '%02d:%02d:%06.3f' % (lst[0],lst[1],lst[2])
        # hdr.set('LST-END',lst_fmt, "HH:MM:SS.SS LST at exposure end")

        # frame ID
        hdr.set('FRAMEID',frameid, "Image sequential number")

        # readout mode
        # if mode == 0:
        #     if itime == 0:
        #         hdr.set('DATA-TYP','BIAS', "Type / Characteristics of this data")
        #     else:
        #         hdr.set('DATA-TYP','OBJECT', "Type / Characteristics of this data")
        # elif mode == 1:
        #         hdr.set('DATA-TYP','DARK', "Type / Characteristics of this data")
        # else:
        #     hdr.set('DATA-TYP','UNKNOWN')

        # exposure time 
        hdr.set('EXPTIME',float(itime), "Total integration time of the frame (sec)")

        if fullHeader is False:
            return hdr
        
        # object name 
        hdr.set('OBJECT', target, "Target description")

        # detector temperature (TBD)
        # hdr.set('DET-TMP', 0.0, "Detector temperature (K)")

        # Telescope header 
        for i in range(len(self.tel_header)):
            if self.tel_header[i][0] == 'NA':
                hdr.set(self.tel_header[i][1], self.tel_header[i][3], self.tel_header[i][4])
            else:
                hdr.set(self.tel_header[i][1], self.statusDictTel[self.tel_header[i][0]], self.tel_header[i][4])

        # WCS parameters 
        imrpad = float(self.statusDictAO['AON.IMR.PAD']) # degrees

        # convert ra/dec to degrees 
        ra = self.statusDictTel['FITS.SBR.RA']
        dec = self.statusDictTel['FITS.SBR.DEC']
        ra_param = ra.split(":")
        ra_deg = 15.0*(float(ra_param[0]) + float(ra_param[1])/60.0 + float(ra_param[2])/3600.0)  
        dec_param = dec.split(":")
        if dec_param[0].find("-") == -1:
            dec_deg = float(dec_param[0]) + float(dec_param[1])/60.0 + float(dec_param[2])/3600.0
        else:
            dec_deg = float(dec_param[0]) - float(dec_param[1])/60.0 - float(dec_param[2])/3600.0

        sin_pa = math.sin(imrpad * math.pi / 180.0)
        cos_pa = math.cos(imrpad * math.pi / 180.0)
        if False:
            pixscale = 1.54321e-5 # degrees
            cd1_1 = pixscale * bin * sin_pa
            cd1_2 = -pixscale * bin * cos_pa
            cd2_1 = pixscale * bin * cos_pa
            cd2_2 = pixscale * bin * sin_pa
            cdelta1 = pixscale
            cdelta2 = pixscale
            pc001001 = -sin_pa
            pc001002 = cos_pa
            pc002001 = -cos_pa
            pc002002 = -sin_pa
            hdr.set('OBS-MOD', 'IMAG')

            hdr.set('CRVAL1', ra_deg, 'Physical value of the reference pixel X')
            hdr.set('CRVAL2', dec_deg, 'Physical value of the reference pixel Y')
            hdr.set('CRPIX1', 2140.0, 'Reference pixel in X (pixel)')
            hdr.set('CRPIX2', 1064.0, 'Reference pixel in Y (pixel)')
            hdr.set('CTYPE1', 'RA---TAN', 'Units used in both CRVAL1 and CDELT1')
            hdr.set('CTYPE2', 'DEC--TAN', 'Units used in both CRVAL2 and CDELT2')
            hdr.set('CUNIT1', 'degree', 'Units used in both CRVAL1 and CDELT1')
            hdr.set('CUNIT2', 'degree', 'Units used in both CRVAL2 and CDELT2')
            hdr.set('CDELT1', cdelta1, 'Size projected into a detector pixel X')
            hdr.set('CDELT2', cdelta2, 'Size projected into a detector pixel Y')
            hdr.set('PC001001', pc001001, 'Pixel Coordinate translation matrix')
            hdr.set('PC001002', pc001002, 'Pixel Coordinate translation matrix')
            hdr.set('PC002001', pc002001, 'Pixel Coordinate translation matrix')
            hdr.set('PC002002', pc002002, 'Pixel Coordinate translation matrix')
            hdr.set('CD1_1', cd1_1, 'Pixel Coordinate translation matrix')
            hdr.set('CD1_2', cd1_2, 'Pixel Coordinate translation matrix')
            hdr.set('CD2_1', cd2_1, 'Pixel Coordinate translation matrix')
            hdr.set('CD2_2', cd2_2, 'Pixel Coordinate translation matrix')


        hdr['COMMENT'] = "------------------------------------------------------------------------"
        hdr['COMMENT'] = "---------------- Parameters for AO188/LGS -- ---------------------------"
        hdr['COMMENT'] = "------------------------------------------------------------------------"

        # AO188 header 
        for i in range(len(self.ao_header)):
            if self.ao_header[i][0] == 'NA':
                hdr.set(self.ao_header[i][1], self.ao_header[i][3], self.ao_header[i][4])
            else:
                hdr.set(self.ao_header[i][1], self.statusDictAO[self.ao_header[i][0]], self.ao_header[i][4])
                          
        return hdr
            
    #######################################
    # INSTRUMENT COMMANDS
    #######################################

    def obcp_mode(self, motor='OFF', mode=None, tag=None):
	"""One of the commands that are in the SOSSALL.cd
        """
        self.mode = mode

    def sleep(self, tag=None, sleep_time=0):

        itime = float(sleep_time)

        # extend the tag to make a subtag
        subtag = '%s.1' % tag

        # Set up the association of the subtag in relation to the tag
        # This is used by integgui to set up the subcommand tracking
        # Use the subtag after this--DO NOT REPORT ON THE ORIGINAL TAG!
        self.ocs.setvals(tag, subpath=subtag)

        # Report on a subcommand.  Interesting tags are:
        # * Having the value of float (e.g. time.time()):
        #     task_start, task_end
        #     cmd_time, ack_time, end_time (for communicating systems)
        # * Having the value of str:
        #     cmd_str, task_error

        self.ocs.setvals(subtag, task_start=time.time(),
                         cmd_str='Sleep %f ...' % itime)

        self.logger.info("\nSleeping for %f sec..." % itime)
        while int(itime) > 0:
            self.ocs.setvals(subtag, cmd_str='Sleep %f ...' % itime)
            sleep_time = min(1.0, itime)
            time.sleep(sleep_time)
            itime -= 1.0

        self.ocs.setvals(subtag, cmd_str='Awake!')
        self.logger.info("Woke up refreshed!")
        self.ocs.setvals(subtag, task_end=time.time())


    def obcp_mode(self, motor='OFF', mode=None, tag=None):
	"""One of the commands that are in the SOSSALL.cd
        """
        self.mode = mode

    def fits_file(self, motor='OFF', frame_no=None, target=None, template=None, delay=0,
                  tag=None):
	"""One of the commands that are in the SOSSALL.cd.
        """

        self.logger.info("fits_file called...")

	if not frame_no:
	    return 1

        # TODO: make this return multiple fits files
	if ':' in frame_no:
	    (frame_no, num_frames) = frame_no.split(':')
	    num_frames = int(num_frames)
        else:
            num_frames = 1

        # Check frame_no
        match = re.match('^(\w{3})(\w)(\d{8})$', frame_no)
        if not match:
            raise CHARISError("Error in frame_no: '%s'" % frame_no)

        inst_code = match.group(1)
        frame_type = match.group(2)
        # Convert number to an integer
        try:
            frame_cnt = int(match.group(3))
        except ValueError, e:
            raise CHARISError("Error in frame_no: '%s'" % frame_no)

        statusDict = {
            'FITS.CRS.PROP-ID': 'None',
            'FITS.CRS.OBSERVER': 'None',
            'FITS.CRS.OBJECT': 'None',
            }
        try:
            res = self.ocs.requestOCSstatus(statusDict)
            self.logger.debug("Status returned: %s" % (str(statusDict)))

        except CHARISError, e:
            return (1, "Failed to fetch status: %s" % (str(e)))

        # Iterate over number of frames, creating fits files
        frame_end = frame_cnt + num_frames
        framelist = []

        while frame_cnt < frame_end:
            # Construct frame_no and fits file
            frame_no = '%3.3s%1.1s%08.8d' % (inst_code, frame_type, frame_cnt)
            if template is None:
                fits_f = pyfits.HDUList(pyfits.PrimaryHDU())
            else:
                templfile = os.path.abspath(template)
                if not os.path.exists(templfile):
                    raise CHARISError("File does not exist: %s" % (templfile))

                fits_f = pyfits.open(templfile)

            hdu = fits_f[0]
            updDict = {'FRAMEID': frame_no,
                       'EXP-ID': frame_no,
                       }

            self.logger.info("updating header")
            for key, val in updDict.items():
                hdu.header.update(key, val)

            subaruCards = self.return_new_header(frame_cnt, 'blank', 0.0)
            hdu.header.extend(subaruCards)

            fitsfile = '/tmp/%s.fits' % frame_no
            try:
                os.remove(fitsfile)
            except OSError:
                pass
            fits_f.writeto(fitsfile, output_verify='ignore')
            fits_f.close()

            # Add it to framelist
            framelist.append((frame_no, fitsfile))

            frame_cnt += 1

        # self.logger.debug("done exposing...")

        # If there was a non-negligible delay specified, then queue up
        # a task for later archiving of the file and terminate this command.
        if delay:
            if type(delay) == type(""):
                delay = float(delay)
            if delay > 0.1:
                # Add a task to delay and then archive_framelist
                self.logger.info("Adding delay task with '%s'" % \
                                 str(framelist))
                t = common_task.DelayedSendTask(self.ocs, delay, framelist)
                t.initialize(self)
                self.threadPool.addTask(t)
                return 0

        # If no delay specified, then just try to archive the file
        # before terminating the command.
        self.logger.info("Submitting framelist '%s'" % str(framelist))
        self.ocs.archive_framelist(framelist)


    def grism(self, tag=None, pos=None):
        # extend the tag to make a subtag
        subtag = '%s.1' % tag
        self.ocs.setvals(tag, subpath=subtag)

        # Report on a subcommand.  Interesting tags are:
        # * Having the value of float (e.g. time.time()):
        #     task_start, task_end
        #     cmd_time, ack_time, end_time (for communicating systems)
        # * Having the value of str:
        #     cmd_str, task_error

        self.ocs.setvals(subtag, task_start=time.time(),
                         cmd_str='NOT really Moving grism stage....')

        self.logger.info("Would set the grism position to %s" % (pos))
        time.sleep(2)

        self.ocs.setvals(subtag, cmd_str='Moved! (liar!)')
        self.logger.info("Woke up refreshed!")
        self.ocs.setvals(subtag, task_end=time.time(), cmd_str="Done.")

    def shutter(self, tag=None, pos=None):
        # extend the tag to make a subtag
        subtag = '%s.1' % tag
        self.ocs.setvals(tag, subpath=subtag)

        # Report on a subcommand.  Interesting tags are:
        # * Having the value of float (e.g. time.time()):
        #     task_start, task_end
        #     cmd_time, ack_time, end_time (for communicating systems)
        # * Having the value of str:
        #     cmd_str, task_error

        self.ocs.setvals(subtag, task_start=time.time(),
                         cmd_str="moving shutter to %s" % (pos))


        self.ocs.setvals(subtag, cmd_str='Moved! ')
        self.ocs.setvals(subtag, task_end=time.time())

        self.logger.info("Would set the shutter position to %s" % (pos))
        
    def filter(self, tag=None, name=None):
        # extend the tag to make a subtag
        subtag = '%s.1' % tag
        self.ocs.setvals(tag, subpath=subtag)

        # Report on a subcommand.  Interesting tags are:
        # * Having the value of float (e.g. time.time()):
        #     task_start, task_end
        #     cmd_time, ack_time, end_time (for communicating systems)
        # * Having the value of str:
        #     cmd_str, task_error

        self.ocs.setvals(subtag, task_start=time.time(),
                         cmd_str="moving filter to %s" % (name))

        ret = self.execOneCmd('charis', 'filter '+name, timelim=90)
        self.ocs.setvals(subtag, cmd_str='Moved! ')
        self.ocs.setvals(subtag, task_end=time.time())

        self.logger.info("Would set the filter to %s" % (name))
        
    def ramp(self, tag=None, exptype='TEST', exptime=0.0, nreset=1, nread=2, target=None):

        # extend the tag to make a subtag
        subtag = '%s.1' % tag
        self.ocs.setvals(tag, subpath=subtag)

        # Report on a subcommand.  Interesting tags are:
        # * Having the value of float (e.g. time.time()):
        #     task_start, task_end
        #     cmd_time, ack_time, end_time (for communicating systems)
        # * Having the value of str:
        #     cmd_str, task_error

        if nread > 0 and exptime > 0:
            self.ocs.setvals(subtag, task_error='Either exptime OR nread can be set. Not both.')
            return

        # frames = self.reqframes(type='A9')
        
        self.ocs.setvals(subtag, task_start=time.time(),
                         cmd_str="Taking a %s ramp(%s, %s)" % (exptype, nreset, nread))

        if nread > 0:
            readArg = "nread=%d" % (nread)
            timelim = (nread + nreset)*1.5 + 20
        else:
            readArg = "itime=%0.1f" % (exptime)
            timelim = nreset*1.5 + exptime + 20
            
        self.execOneCmd('hx', 'ramp nreset=%d %s' % (nreset, readArg),
                        timelim=timelim, subtag=subtag)
        self.logger.info("Would take a %s ramp(%s, %s)" % (exptype, nreset, nread))

        self.ocs.setvals(subtag, cmd_str='Moved! ')
        self.ocs.setvals(subtag, task_end=time.time(),
                         cmd_str="Done with ramp")

        
    def putstatus(self, target="ALL"):
        """Forced export of our status.
        """
	# Bump our status send count and time
	self.stattbl1.count += 1
	self.stattbl1.time = time.strftime("%4Y%2m%2d %2H%2M%2S",
                                           time.localtime())

        self.ocs.exportStatus()


    def getstatus(self, target="ALL"):
        """Forced import of our status using the normal status interface.
        """
	ra, dec, focusinfo, focusinfo2 = self.ocs.requestOCSstatusList2List(['STATS.RA',
                                                      'STATS.DEC',
                                                      'TSCV.FOCUSINFO',
                                                      'TSCV.FOCUSINFO2'])

        self.logger.info("Status returned: ra=%s dec=%s focusinfo=%s focusinfo2=%s" % (ra, dec, focusinfo, focusinfo2))


    def getstatus2(self, target="ALL"):
        """Forced import of our status using the 'fast' status interface.
        """
	ra, dec = self.ocs.getOCSstatusList2List(['STATS.RA',
                                                  'STATS.DEC'])

        self.logger.info("Status returned: ra=%s dec=%s" % (ra, dec))

    def view_file(self, path=None, num_hdu=0, tag=None):
        """View a FITS file in the OCS viewer.
        """
        self.ocs.view_file(path, num_hdu=num_hdu)


    def view_fits(self, path=None, num_hdu=0, tag=None):
        """View a FITS file in the OCS viewer
             (sending entire file as buffer, no need for astropy).
        """
        self.ocs.view_file_as_buffer(path, num_hdu=num_hdu)


    def reqframes(self, num=1, type="A"):
        """Forced frame request.
        """
        framelist = self.ocs.getFrames(num, type)

        # This request is not logged over DAQ logs
        self.logger.info("framelist: %s" % str(framelist))

        return framelist

    def kablooie(self, motor='OFF'):
	"""Generate an exception no matter what.
        """
        raise CHARISError("KA-BLOOIE!!!")


    def defaultCommand(self, *args, **kwdargs):
        """This method is called if there is no matching method for the
        command defined.
        """

        # If defaultCommand is called, the cmdName is pushed back on the
        # argument tuple as the first arg
        cmdName = args[0]
        self.logger.info("Called with command '%s', params=%s" % (cmdName,
                                                                  str(kwdargs)))

        res = unimplemented_res
        self.logger.info("Result is %d\n" % res)

        return res

    def power_off(self, upstime=None):
        """
        This method is called when the summit has been running on UPS
        power for a while and power has not been restored.  Effect an
        orderly shutdown.  upstime will be given the floating point time
        of when the power went out.
        """
        res = 1
        try:
            self.logger.info("!!! POWERING DOWN !!!")
            # res = os.system('/usr/sbin/shutdown -h 60')

        except OSError, e:
            self.logger.error("Error issuing shutdown: %s" % str(e))

        self.stop()

        self.ocs.shutdown(res)
