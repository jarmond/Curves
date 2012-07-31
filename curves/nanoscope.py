# Curves force analysis software
# JWA (c)2009
# Import Nanoscope force curve

import re
from numpy import *

class NanoscopeFileReadError(Exception):
    def __init__(self,err,line):
        self.err = err
        self.line = line
    
    def __str__(self):
        return repr(self.err) + ': ' + repr(self.line)

def is_nanoscope_file(filename):
    """Returns True if filename refers to a force curve file. Force curve files
    are identified by the line \*Force file list at the beginning."""
    infile = open(filename,"r") 
    line = infile.readline().strip()
    return line.startswith("\*Force file list")

def read_params(filename):
    # Create dictionary of parameters, indexed by parameter list section
    # and parameter name
    infile = open(filename,"r")
    raw_param = {}
    # TODO: check nanoscope and version
    while infile:
        line = infile.readline().strip()
        # Stop at end marker
        if line.startswith("\*File list end"):
            pass
            #print "Finished scanning parameters."
            break
        # If begins with \* then is a section header
        elif line.startswith("\*"):
            regex_str = r"(\\\*)([\w ]+)"
            match_obj = re.search(regex_str,line)
            if match_obj:
                list_name = cleanKey(match_obj.group(2))
                #print "Section header: " + list_name
                raw_param[list_name] = {}
            else:
                raise NanoscopeFileReadError('Section header',line)
        # Is a parameter, parse, add to dictionary
        elif len(line)>0: # Ignore empty lines
            # Format is: \Parameter name: values
            #print line
            regex_str = r"(\\)((@\d:)?[\w @\-\./\(\),\[\]]*)(\t*:) *(.*)"
            match_obj = re.search(regex_str,line)
            # Grp 1 = \
            # Grp 2 = Param name including @? prefix =
            # Grp 3 = @? prefix
            # Grp 4 = :
            # Grp 5 = Values
            if match_obj:
                #print "Parsing: " + match_obj.group(3)
                values = match_obj.group(5)
                #print "Values: " + match_obj.group(5)
                raw_param[list_name][cleanKey(match_obj.group(2))] = values
            else:
                raise NanoscopeFileReadError('Parameter read',line)
        else:
            # Ignore
            pass
    infile.close()
    return raw_param

def read_curve_data(filename,std_params):
    # Open in binary mode
    infile = open(filename,"rb")
    # Seek to data offset
    data_offset = std_params["DataOffset"]
    data_length = std_params["DataLength"]
    bytes_per_point = std_params["BytesPerPoint"]
    point_count = data_length/bytes_per_point
    points_per_curve = point_count/std_params['NumCurves']
    infile.seek(data_offset)
    # Some curves are truncated, and are preceded with junk data. So we need
    # to work out how many bytes to skip
    approach_skip = (points_per_curve - std_params['ApproachPoints'])*bytes_per_point
    retract_skip = (points_per_curve - std_params['RetractPoints'])*bytes_per_point
    # Read raw data and shape
    # \todo allow for different integer lengths
    if bytes_per_point!=2:
        print "Bytes/point != 2"
    approach = fromfile(file=infile,dtype=int16,count=std_params['ApproachPoints'])
    infile.seek(approach_skip,1) # 1 indicates relative seek
    retract = fromfile(file=infile,dtype=int16,count=std_params['RetractPoints'])
    
    # Convert raw deflection data to Volts
    approach = approach * std_params["HardScale"]
    retract = retract * std_params["HardScale"]

    # Apply Z range
    ramp_offset = std_params["RampOffset"]
    ramp_size = std_params["RampSize"]
    z = linspace(ramp_offset,ramp_offset+ramp_size,max(len(approach),len(retract)))
   
    return (z,approach,retract)

def read_curve_info(raw_params):
    """Return a dictionary of parameters with standardised string keys"""
    # G1: soft scale param, G2: hard scale, G3: hard value
    regex_hard = r"V \[([\w \.]+)\] \(([\d\.]+) V/LSB\) ([-\d\.]+)" # V [param] (value V/LSB)
    # G1: soft value
    regex_soft = r"V ([\d\.]+)" # V soft value
    regex_time = r"([\d\:]+ (PM|AM)) ([\w ]+)" # 1: time, 2: AM or PM, 3: date
    regex_version = r"0x0(\d)(\d\d)" # 1: major version, 2: minor
    regex_samples = r"([\d]+) ([\d]+)"
    std_params = {}

    # Decide which nanoscope version
    m = re.search(regex_version,raw_params['forcefilelist']['version'])
    if m:
        std_params['Version'] = m.group(1) + '.' + m.group(2)
    else:
        raise NanoscopeFileReadError('Version',raw_params['forcefilelist']['version'])
    version = float(std_params['Version'])

    if version >= 7 and version < 8:
        ns = ns7
    elif version >= 6:
        ns = ns6
    else:
        raise NanoscopeFileReadError('Version %f not supported' % version)

    # Some minor versions annoyingly differ in spaces between words in keys, e.g. Z Limit vs ZLimit
    # Clean key strings
    ns = dict(zip(ns.keys(),map(cleanKey, ns.values())))

    # Nanoscope 7 gives AFM for this
    #   if raw_params['Ciao force image list']['Data type'].strip() != 'FORCE':
    #       print('Expected data type = FORCE. Possibly not force curve.')
    
    std_params['Software'] = 'Nanoscope'

    # Data scaling parameters
    m = re.search(regex_hard,raw_params[ns['ImageList']][ns['HardScale']])
    if m:
        soft_param = m.group(1)
        #std_params['HardScale'] = float(m.group(3))/65536 # change made here (divide by 65536 and use group 3), cf manual
        std_params['HardScale'] = float(m.group(2))
    else:
        raise NanoscopeFileReadError('Z scale',raw_params[ns['ImageList']][ns['HardScale']])
    m = re.search(regex_soft,raw_params[ns['ScanList']][cleanKey('@' + soft_param)])
    if m:
        std_params['DeflSens'] = float(m.group(1))
    else:
        raise NanoscopeFileReadError('Z scale',raw_params[ns['ScanList']]['@' + soft_param])
    m = re.search(regex_hard,raw_params[ns['ForceList']][ns['RampSize']]) # identify zsens param
    if m:
        ramp_soft_param = m.group(1)
    else:
        raise NanoscopeFileReadError('Z scan size',raw_params[list_name][ns['RampSize']])
    m = re.search(regex_soft,raw_params[ns['ScannerList']][cleanKey('@' + ramp_soft_param)])
    if m:
        std_params['ZSens'] = float(m.group(1))
    else:
        raise NanoscopeFileReadError('Z sens',raw_params[ns['ScannerList']][cleanKey('@' +
                                                                                     ramp_soft_param)])
    m = re.search(regex_hard,raw_params[ns['ScanList']][ns['ZRange']])
    if m:
        std_params['ZRange'] = float(m.group(3)) * std_params['ZSens']
    else:
        raise NanoscopeFileReadError('Z range',raw_params[ns['ScanList']][ns['ZRange']])


    # Force experimental parameters
    list_name = ns['ForceList']
    std_params['ScanRate'] = float(raw_params[list_name][ns['ScanRate']])
    std_params['ForwardVelocity'] = float(raw_params[list_name][ns['ForwardVelocity']])*std_params['ZSens']
    std_params['ReverseVelocity'] = float(raw_params[list_name][ns['ReverseVelocity']])*std_params['ZSens']
    m = re.search(regex_hard,raw_params[list_name][ns['RampSize']])
    std_params['RampSize'] = float(m.group(3))*std_params['ZSens']
    m = re.search(regex_hard,raw_params[list_name][ns['RampOffset']])
    std_params['RampOffset'] = float(m.group(3))*std_params['ZSens']
    std_params['SurfaceDelay'] = float(raw_params[list_name][ns['SurfaceDelay']])
    std_params['RetractDelay'] = float(raw_params[list_name][ns['RetractDelay']])
    
    # Parameters
    std_params['SpringConst'] = float(raw_params[ns['ImageList']][ns['SpringConst']])
        
    # Information
    list_name = ns['FileList']
    m = re.search(regex_time,raw_params[list_name][ns['Date']])
    if m:
        std_params['Date'] = m.group(3)
        std_params['Time'] = m.group(1)
    else:
        raise NanoscopeFileReadError('Date and time',raw_params[list_name][ns['Date']])

    list_name = ns['EquipmentList']
    std_params['Microscope'] = raw_params[list_name][ns['Microscope']]
    if 'Controller' in ns:
        std_params['Controller'] = 'Nanoscope ' + raw_params[list_name][ns['Controller']]
    else:
        std_params['Controller'] = ''

    list_name = ns['ScannerList']
    std_params['ScannerType'] = raw_params[list_name][ns['ScannerType']]
    std_params['PiezoSize'] = raw_params[list_name][ns['PiezoSize']]

    list_name = ns['ScanList']
    std_params['TipNumber'] = raw_params[list_name][ns['TipNumber']]
    

    list_name = ns['ImageList']
    std_params['DataOffset'] = int(raw_params[list_name][ns['DataOffset']],10)
    std_params['DataLength'] = int(raw_params[list_name][ns['DataLength']],10)
    std_params['BytesPerPoint'] = int(raw_params[list_name][ns['BytesPerPoint']],10)
    m = re.search(regex_samples,raw_params[list_name][ns['SamplesPerLine']])
    if m:
        std_params['RetractPoints'] = int(m.group(1))
        std_params['ApproachPoints'] = int(m.group(2))
        std_params['NumCurves'] = len(m.groups())
    else:
        raise NanoscopeFileReadError('Samps/line',raw_params[list_name][ns['SamplesPerLine']])
    
    
    return std_params


def cleanKey(key):
  try:
    return key.strip().replace(' ','').lower()
  except AttributeError:
    return key

# Nanoscope software keys
# Nanoscope v6.11 - v6.13
ns6 = {
    'Version' : 6,
    'Date' : 'Date',
    'Microscope' : 'Microscope',
    'ForceList' : 'Ciao force list',
    'ImageList' : 'Ciao force image list',
    'FileList' : 'Force file list',
    'EquipmentList' : 'Equipment list',
    'ScannerList' : 'Scanner list',
    'ScanList' : 'Ciao scan list',
    'ScanRate' : 'Scan rate',
    'ForwardVelocity' : 'Forward vel.',
    'ReverseVelocity' : 'Reverse vel.',
#    'RampSize' : '@Z scan size',
    'RampSize' : '@4:Ramp size Zsweep',
#    'RampOffset' : '@Z scan start',
    'RampOffset' : '@4:Ramp offset Zsweep',
    'SpringConst' : 'Spring Constant',
    'Controller' : 'Controller',
    'ScannerType' : 'Scanner type',
    'PiezoSize' : 'Piezo size',
    'TipNumber' : 'Tip serial number',
    'HardScale' : '@4:Z scale',
    'ZRange' : '@1:Z limit',
    'DataOffset' : 'Data offset',
    'DataLength' : 'Data length',
    'BytesPerPoint' : 'Bytes/pixel',
    'SamplesPerLine' : 'Samps/line',
    'DataType' : 'FORCE',
    'SurfaceDelay' : 'Ramp delay',
    'RetractDelay' : 'Reverse delay'
}

# Nanoscope v7.3
ns7 = {
    'Version' : 7,
    'Date' : 'Date',
    'Microscope' : 'Microscope',
    'ForceList' : 'Ciao force list',
    'ImageList' : 'Ciao force image list',
    'FileList' : 'Force file list',
    'EquipmentList' : 'Equipment list',
    'ScannerList' : 'Scanner list',
    'ScanList' : 'Ciao scan list',
    'ScanRate' : 'Scan rate',
    'ForwardVelocity' : 'Forward vel.',
    'ReverseVelocity' : 'Reverse vel.',
    'RampSize' : '@Z scan size',
    'RampOffset' : '@Z scan start',
    'SpringConst' : 'Spring Constant',
    'ScannerType' : 'Scanner type',
    'PiezoSize' : 'Piezo size',
    'TipNumber' : 'Tip Serial Number',
    'HardScale' : '@4:Z scale', # maybe not but have to divide by 65536 (2^16): this might include spring const and 'Z magnify force'
    'ZRange' : '@1:Z Limit',
    'DataOffset' : 'Data offset',
    'DataLength' : 'Data length',
    'BytesPerPoint' : 'Bytes/pixel',
    'SamplesPerLine' : 'Samps/line',
    'DataType' : 'AFM',
    'SurfaceDelay' : 'Ramp delay',
    'RetractDelay' : 'Reverse delay'
}
     
# Units
ns_units = {
    'Version' : '',
    'Date' : '',
    'Time' : '',
    'Microscope' : '',
    'Software' : '',
    'TipNumber' : '',
    'ScanRate' : 'Hz',
    'ForwardVelocity' : 'nm/s',
    'ReverseVelocity' : 'nm/s',
    'RampSize' : 'nm',
    'RampOffset' : 'nm',
    'SpringConst' : 'N/m',
    'Controller' : '',
    'ScannerType' : '',
    'PiezoSize' : '',
    'TipNumber' : '',
    'HardScale' : 'V/LSB',
    'ZRange' : 'nm',
    'ZSens' : 'nm/V',
    'DataOffset' : 'bytes',
    'DataLength' : 'bytes',
    'BytesPerPoint' : 'bytes',
    'SamplesPerLine' : 'Samps/line',
    'SurfaceDelay' : 's',
    'RetractDelay' : 's',
    'DeflSens': 'nm/V',
    'RetractPoints' : '',
    'ApproachPoints' : '',
    'NumCurves': '',
}
