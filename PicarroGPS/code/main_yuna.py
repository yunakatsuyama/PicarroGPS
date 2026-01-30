# Created by Jakob Borchardt, University of Bremen
# bjakob@iup.physik.uni-bremen.de
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read buffer files and create according kml files

.. currentmodule:: PicarroGPS.code.main

.. autosummary::
    :toctree: api/

    read_config
    strfrommsg
    track_gps
    track_to_kml
    write_msg
    export_colorbar
"""

import pdb

import py
import os
import shutil
import pathlib
import time
from contextlib import ExitStack
import collections
import datetime
import configparser
import pynmea2
import numpy as np
import pandas
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
#from ..code import set_systime as sst
#from ..code import kml as KML
from . import set_systime as sst
from . import kml as KML
import warnings
import glob
import subprocess

filename = '/Users/mameyuna/stu/Romania_data/LGR 2026-01-23 Software Backup/PythonPrograms/PicarroGPSconfig.cfg'
#def read_config(filename='PicarroGPSconfig.cfg'):
def read_config(filename = filename):
    """Reads config file and returns config dictionary.

    Parameters
    ----------
    filename : :class:`str <python:str>`
        relative or absolute path to file and filename of the config file.

    Returns
    -------
    :class:`dict <python:dict>`
        Dictionary containing the config values.

    Notes
    -----
    Conversion of the different fields should take place in this routine


    """

    config = configparser.ConfigParser(allow_no_value=True)
    file = os.path.abspath(filename)
    config.read(file)
    #files_read = config.read(file)
    #print("Config file read:", files_read)
    #print("Sections found:", config.sections())

    return config


def strfrommsg(msg):
    """Read longitude, latitude and altitude from pynmea message and return
    them as a string

    .. note:: Deprecated in 0.1

        Deprecated as functionality is special and not reusable. Will be removed
        in future releases

    Parameters
    ----------
    msg : `pynmea2 message`
        Pynmea2 gps message

    Returns
    -------
    str : :class:`str <python:str>`
        String of the format 'lat,lon,alt' or empty string if wrong message

    """

    with warnings.catch_warnings():
            warnings.filterwarnings('always', category=DeprecationWarning)
            warnings.warn(
                'The strfrommsg function is deprecated as its functionality is'
                'incorporated in the standard program.', DeprecationWarning)

    if msg.sentence_type == 'GGA':
        string = ','.join([str(msg.longitude),
                           str(msg.latitude),
                           str(msg.altitude)])
    else:
        string = ''

    return string


def track_gps(filename='FlightTrack'):
    """Uses a gps mouse via serial port to log the position to a file

    The tracking is aborted by pressing <ctrl> + <c>

    Parameters
    ----------
    filename : :class:`str <python:str>`
        Name prefix of the file to which the gps messages are logged. The
        output file will be named filename_gps.dat


    """

    try:
        portdict = sst.get_COMconfig()
    except KeyError:
        portdict = {'port': None}

    stream = pynmea2.NMEAStreamReader()

    with sst.get_gpscom(**portdict) as gps_ser:
        with open(filename + '_gps.dat', 'w') as gpsfile:
            while 1:
                try:
                    data = gps_ser.readline()
                    for msg in stream.next(data.decode()):
                        gpsfile = write_msg(msg, gpsfile)
                except (pynmea2.nmea.SentenceTypeError, ValueError,
                        AttributeError):
                    pass
                except KeyboardInterrupt:
                    print('Keyboard Interrupt detected, closing files and '
                          'serial connection stopping tracking of gps')
                    break


def track_to_kml(filename='Picarro', timelag=16, vrange=[380.0,420.0,1.80,2.0,0.5,1.5],
                 reprocess=False, reprocessfolder='', skipfiles=[0]):
    """Read buffer files and write them to a kml file with gps coordinates.

    This function reads in special buffer files and extracts concentration and
    gps information. Additionally, 5 kml files are produced: ch4, co2 and h2o
    concentration, the current position and the flight path.

    The concentrations are displayed using a colorbar. The colorbar range can
    be changed in the config file. The colorbars are written to additional kml
    files in the beginning.

    The program is terminated when an :class:`KeyboardInterruptError <python:KeyboardInterruptError>`
    is detected.

    .. note:: Parameters may change!

        Parameters may change or be reduced in the future, as the values should be
        read in from the config file.

    Parameters
    ----------
    filename : :class:`str <python:str>`
        Name prefix for the kml files. Should only be changed when you know
        what you are doing, as the live view depends on fixed kml file names!
    timelag : :class:`int <python:int>`
        Timelag between sampling the signal with the aircraft and appearance
        of the signal at the instrument. For remote sensing this should be
        0, otherwise it should be optimized for the current instrument.

        With the FU Berlin Cessna and a LosGatos GGA the value was around 11,
        for a Picarro it was around 4-5.
    vrange : :class:`list <python:list>`
        List of :class:`floating point numbers <python:float>`. Must have at
        least 6 entries. These are lower and upper boundaries for ch4, co2 and
        h2o colorbars, respectively.
    reprocess : :class:`boolean <python:bool>` or similar
        Flag if reprocess folder should be searched for buffer files and if
        there are any, reprocess them. This can be any representation of `True`
        or `False`.
    reprocessfolder : :class:`str <python:str>` or path
        Folder to which the buffer files are copied and where
        the function searches for files at the beginning if `reprocess` is set.
    skipfiles : :class:`int <python:int>`
        Integer defining the number of files to be skipped when reprocessing.

    See also
    --------
    :func:`PicarroGPS.code.set_systime.get_COMconfig` : get port configuration
    :func:`PicarroGPS.code.set_systime.get_gpscom` : set up serial connection
    """

    config = read_config(filename)
    print(f"here config{config}")
    print("CONFIG SECTIONS:", config.sections())


    reprocess = config['Paths'].getboolean('reprocess')

    reprocessfolder = config['Paths']['reprocessfolder']
    skipfiles = [int(sfr) for sfr in config['Paths']['skipfiles'].split()]
    max_points = int(config['Paths']['max_len']) + 17  # 16 styles plus header are additional children
    print(f"config check {config['Paths']['iconfolder']}")
    try:
        timelag = config['Picarro'].getint('timelag')
        vrange = [float(num) for num in config['Picarro']['vrange'].split()]
        filenames = config['Picarro']['filenames'].split()
    except:
        pass

    if config['Picarro'].getboolean('dryflag'):
        if config['Picarro'].getboolean('coflag'):
            columns = [0,5,7,8,10,11,12,13,14]
            datanames = ['0','1','2','3','8','4','5','6','7']
        else:
            columns = [5,7,8,9,10,11,12,0]
            datanames = ['0','1','2','3','4','5','6','7']
    else:
        if config['Picarro'].getboolean('coflag'):
            columns = [0,4,6,8,9,11,12,13,14]
            datanames = ['0','1','2','3','8','4','5','6','7']
        else:
            columns = [0,4,6,8,9,10,11,12]
            datanames = ['0','1','2','3','4','5','6','7']


    curr_pos = KML.kmlinit(filename + '_pos')
    flighttrack = KML.kmlinit(filename + '_track')
    concCH4 = KML.kmlinit(filename + '_CH4')
    concCO2 = KML.kmlinit(filename + '_CO2')
    concH2O = KML.kmlinit(filename + '_H2O')

    iconfolder = config['Paths']['iconfolder']
    print(f"here config iconfolder {config['Paths']['iconfolder']}")
    KML.add_style(curr_pos, stid='pyst1',
                  iconStyle={'scale': '1.5',
                             'color': 'ff80c0ff',
                             'Icon': {'href': iconfolder + '/icon46.png'}}
                  )
    KML.add_style(flighttrack, stid='pyst1',
                  lineStyle={'color': 'ffffffff', 'width': '5'})
    ch4_cmap, nbins = KML.add_colorscheme(concCH4, nbins=config['Picarro'].getint('nbins'))
    export_colorbar(ch4_cmap, (vrange[2],vrange[3]), 'CH4 in ppm', config)
    co2_cmap, nbins = KML.add_colorscheme(concCO2, nbins=config['Picarro'].getint('nbins'))
    export_colorbar(co2_cmap, (vrange[0],vrange[1]), 'CO2 in ppm', config)
    h2o_cmap, nbins = KML.add_colorscheme(concH2O, nbins=config['Picarro'].getint('nbins'))
    export_colorbar(h2o_cmap, (vrange[4],vrange[5]), 'H2O in %', config)

    KML.add_Placemark(curr_pos, placemark={'name': 'Current Position',
                                           'styleUrl': 'pyst1',
                                           'Point': {'coordinates': '',
                                                     'altitudeMode': 'absolute'}})
    KML.add_Placemark(flighttrack,
                      placemark={'name': 'Flight Track',
                                 'styleUrl': 'pyst1',
                                 'LineString': {'coordinates': '\n        ',
                                                'altitudeMode': 'absolute'},
                                 'tesselate': '1'})

    conckml = [concCO2, concCH4, concH2O]
    vrangearr = [vrange[0:2], vrange[2:4], vrange[4:6]]

    if config['Picarro'].getboolean('coflag'):
        concCO = KML.kmlinit(filename + '_CO')
        co_cmap, nbins = KML.add_colorscheme(concCO, nbins=config['Picarro'].getint('nbins'))
        export_colorbar(co_cmap, (vrange[6],vrange[7]), 'CO in ppm', config)
        conckml += [concCO]
        vrangearr.append(vrange[6:8])


    coordsque = collections.deque(maxlen=timelag)

    remotefolder = config['Paths']['remotefolder']
    backupfolder = config['Paths']['backuppath']

    try:
        os.mkdir(backupfolder)
    except FileExistsError:
        pass

    if reprocess and len(reprocessfolder):
        reprofiles = os.listdir(reprocessfolder)
        for repfile in reprofiles[skipfiles[0]:]:
            repfilepath = pathlib.Path(reprocessfolder + repfile)
            if os.stat(repfilepath).st_size > 0:
                datin = pandas.read_csv(repfilepath, delim_whitespace=True,
                                               header=None,
                                               usecols=columns,
                                               names=datanames)
                data = [np.mean(datin[str(ii+1)]) for ii in range(len(datanames)-1)]
                if data[0] > 0.0:
                    coor = ','.join([str(data[4]), str(data[3]), str(data[5])])

                    coordsque.appendleft(coor)
                    curr_pos.Document.Placemark.Point.replace(
                            curr_pos.Document.Placemark.Point.coordinates,
                            KML.objectify.fromstring(KML.dicttokmlstring(
                                {'coordinates': coor}, curr_pos.nsmap)
                            ))
                    KML.add_coordinate(flighttrack, data[3], data[4], data[5])
                    concentrations = data[0:3]
                    if config['Picarro'].getboolean('coflag'):
                        concentrations.append(data[-1])
                    if len(coordsque) == timelag:
                        coorout = coordsque.pop()
                        for kml, conc, delval in zip(conckml,
                                                     concentrations,
                                                     vrangearr):
                            stid = KML.selectkmlstyle(conc, delval[0], delval[1])
                            KML.add_Placemark(kml, placemark={'description': ' '.join(['Concentration ', str(conc),
                                                                                       '\nAltitude ', (coorout.split(','))[2],
                                                                                       '\nMeasNr ', str(data[6]),
                                                                                       '\nMeasTime ', str(datin['0'][0]),
                                                                                       '\nGPS ', '  '.join((coorout.split(','))[:2])]),
                                                              'styleUrl': stid,
                                                              'Point': {'coordinates': coorout,
                                                                        'altitudeMode': 'absolute'}
                                                             })
                            children = kml.Document.getchildren()
                            if len(children) > max_points:
                                kml.Document.remove(children[18])
                        print('Updated Measurement Nr {0:05d}'.format(int(data[6])))
                        print(len(conckml[0].Document.getchildren()))
#                        print(' '.join(['Concentration ', str(conc),
#                                        '\nAltitude ', (coorout.split(','))[2],
#                                        '\nMeasNr ', str(data[6]),
#                                        '\nMeasTime ', str(datin['0'][0]),
#                                        '\nGPS ', '  '.join((coorout.split(','))[:2])]))
        with open(filename + filenames[0] + '.kml', 'w') as filepos:
            filepos.write(KML.etree.tounicode(curr_pos, pretty_print=True))
        with open(filename + filenames[1] + '.kml', 'w') as filetrck:
            filetrck.write(KML.etree.tounicode(flighttrack, pretty_print=True))
        for filenr, concGAS in enumerate(conckml):
            with open(filename + filenames[filenr+2] + '.kml', 'w') as file:
                file.write(KML.etree.tounicode(concGAS, pretty_print=True))
#            shutil.copy(filename + filenames[filenr+2] + '.kml', filename + filenames[filenr+2]+'/')
#            shutil.make_archive(filename + filenames[filenr+2], 'zip', filename + filenames[filenr+2]+'/')
#            shutil.copy(filename + filenames[filenr+2] + '.zip', filename + filenames[filenr+2] + '.kmz')

    try:
        while True:
            measfiles = os.listdir(remotefolder)
            for file in measfiles[:50]:
                try:
                    shutil.move(remotefolder + file, backupfolder)
                    filepath = pathlib.Path(backupfolder + file)
                except shutil.Error:
                    shutil.move(remotefolder + file, backupfolder + file[:-4] + '_1.dat')
                    filepath = pathlib.Path(backupfolder + file[:-4] + '_1.dat')
                if os.stat(filepath).st_size > 0:
                    try:
                        datin = pandas.read_csv(filepath, delim_whitespace=True,
                                                   header=None,
                                                   usecols=columns,
                                                   names=datanames)
                        data = [np.mean(datin[str(ii+1)]) for ii in range(len(datanames)-1)]
                        coor = ','.join([str(data[4]), str(data[3]), str(data[5])])
                        concentrations = data[0:3]
                        if config['Picarro'].getboolean('coflag'):
                            concentrations.append(data[-1])
    
                        coordsque.appendleft(coor)
                        curr_pos.Document.Placemark.Point.replace(
                                curr_pos.Document.Placemark.Point.coordinates,
                                KML.objectify.fromstring(KML.dicttokmlstring(
                                    {'coordinates': coor}, curr_pos.nsmap)
                                ))
                        KML.add_coordinate(flighttrack, data[3], data[4], data[5])
    
                        if len(coordsque) == timelag:
                            coorout = coordsque.pop()
                            for kml, conc, delval in zip(conckml,
                                                         concentrations,
                                                         vrangearr):
                                stid = KML.selectkmlstyle(conc, delval[0], delval[1])
                                KML.add_Placemark(kml, placemark={'description': ' '.join(['Concentration ', str(conc),
                                                                                           '\n        Altitude ', (coorout.split(','))[2],
                                                                                           '\n        MeasNr ', str(data[6]),
                                                                                           '\n        MeasTime ', str(datin['0']),
                                                                                           '\n        GPS ', '  '.join((coorout.split(','))[:2])]),
                                                                  'styleUrl': stid,
                                                                  'Point': {'coordinates': coorout,
                                                                            'altitudeMode': 'absolute'}
                                                                 })
                                children = kml.Document.getchildren()
                                if len(children) > max_points:
                                    kml.Document.remove(children[18])
                        print('Added Measurement Nr {0:05d}'.format(int(data[6])))
                        print(len(conckml[0].Document.getchildren()))
#                        print(' '.join(['Concentration ', str(conc),
#                                                                                           '\nAltitude ', (coorout.split(','))[2],
#                                                                                           '\nMeasNr ', str(data[6]),
#                                                                                           '\nMeasTime ', str(datin['0']),
#                                                                                           '\nGPS ', '  '.join((coorout.split(','))[:2])]))
                    except ValueError:
                        print('Not valid value in measnr., skipping')
                if data[6] % 1 == 0:
                    with open(filename + filenames[0] + '.kml', 'w') as filepos:
                        filepos.write(KML.etree.tounicode(curr_pos, pretty_print=True))
                    with open(filename + filenames[1] + '.kml', 'w') as filetrck:
                        filetrck.write(KML.etree.tounicode(flighttrack, pretty_print=True))
                    for filenr, concGAS in enumerate(conckml):
                        with open(filename + filenames[filenr+2] + '.kml', 'w') as file:
                            file.write(KML.etree.tounicode(concGAS, pretty_print=True))
#                        try:
#                            os.remove(filename + filenames[filenr+2] + '.kml')
#                        except:
#                            pass
#                        os.rename(filename + filenames[filenr+2] + '.tmp', filename + filenames[filenr+2] + '.kml')
#                        shutil.copy(filename + filenames[filenr+2] + '.kml', filename + filenames[filenr+2]+'/')
#                        subprocess.run(['compact', '/c', filename + filenames[filenr+2] + '/'])
    #                    shutil.make_archive(filename + filenames[filenr+2], 'zip', filename + filenames[filenr+2]+'/')
    #                    shutil.copy(filename + filenames[filenr+2] + '.zip', filename + filenames[filenr+2] + '.kmz')
#                time.sleep(1.0)

            time.sleep(0.2)
    except KeyboardInterrupt:
        print('Keyboard interrupt detected, ending tracking')
        with open(filename + filenames[0] + '.kml', 'w') as filepos:
                    filepos.write(KML.etree.tounicode(curr_pos, pretty_print=True))
        with open(filename + filenames[1] + '.kml', 'w') as filetrck:
            filetrck.write(KML.etree.tounicode(flighttrack, pretty_print=True))
        for filenr, concGAS in enumerate(conckml):
            with open(filename + filenames[filenr+2] + '.kml', 'w') as file:
                file.write(KML.etree.tounicode(concGAS, pretty_print=True))
    finally:
        folder = datetime.datetime.now().strftime(config['Paths']['kmlpath'] + '%y%m%dt%H%M%S_kml_files/')
        os.mkdir(folder)
        for file in glob.glob(filename + '*.kml'):
            shutil.copy2(file, folder)
        shutil.copy2('colorbar.kml', folder)
        os.mkdir(''.join([folder, 'icons/']))
        icons = os.listdir(iconfolder)
        for icon in icons:
            shutil.copy2(''.join([iconfolder, icon]), ''.join([folder, 'icons/']))

    return


def write_msg(msg, fileh):
    """
    Appends a gps message to an already opened file.


    .. note:: Deprecated

        will be removed as track_gps is not used.

    Parameters
    ----------
    msg : `pynmea message`_

    fileh : :term:`file object <python:file object>`
        Must alread be :func:`open <python:open>`. Otherwise, a new file with
        expansion '_err' is opened.

    Returns
    -------
    fileh : :term:`file object <python:file object>`
        If no new file was opened during call to
        :func:`PicarroGPS.code.main.write_msg`, then this is
        the same as the input parameter. Otherwise, the file object of the
        newly opened file will be returned.

    Notes
    -----
    If an `IOError` occurs, a new file with filename expanded by _err will be
    opened and used further.

    """

    try:
        fileh.write(str(msg)+chr(10))
    except IOError:
        if not(fileh.closed):
            fileh.close
        fileh = (fileh.dirname + py.path.local.sep + fileh.purebasename +
                 '_err' + fileh.ext)

    return fileh


def export_colorbar(cmap, bounds, label, cfg):
    """Create colorbar as png for inclusion in kml file.

    Parameters
    ----------
    cmap : `discrete colormap`
        Has to have `cmap.N` attribute

    bounds : 2-element list
        first element: lower bound of values
        second element: upper bound of values

        otherwise negative values will be printed

    label : `str`
        String containing the label for the colorbar.

    cfg : `ConfigParserObj`
        Configuration dictionary read in with configparser. Has to have field
        `Paths` with attribute `iconfolder`

    """

    assert (cmap.N > 0)

    # create figure and norm for colorbar
    fig = plt.figure(figsize=(20,3.5))
    ax = fig.add_axes([0.05,0.51,0.9,0.45])
    bnds = [bounds[0] + (bounds[1]-bounds[0]) / float(cmap.N) * i
             for i in range(cmap.N+1)]
    norm = mpl.colors.BoundaryNorm(bnds, cmap.N)
    ticks = np.arange(bnds[0], bnds[-1]+(bnds[-1]-bnds[0])/5.0, (bnds[-1]-bnds[0])/4.0)
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, boundaries=bnds,
                                   ticks=[round(t,3) for t in ticks], orientation='horizontal')
    cb.ax.tick_params(labelsize=42, pad=20.0)
    cb.set_label(label, size=42, labelpad=20.0)
    # this made error ....
    print(f"cfg cfg['Paths']['iconfolder'] {cfg['Paths']['iconfolder']}")
    fig.savefig(cfg['Paths']['iconfolder'] + 'colorbar_' + label[:4].strip() + '.png')

    return


if __name__ == '__main__':
    pass
