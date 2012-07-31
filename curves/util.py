# Curves force analysis software
# JWA Copyright 19/04/08

import math

from forcecurve import CurveMode
# Utility variables and functions

# Units {Mode: (Long name, Unit, Default Prefix)
unit = {CurveMode.Volts: {'name':'Volts','unit':'V','default_prefix':0},
	CurveMode.Deflection: {'name':'Deflection','unit':'m','default_prefix':-9},
	CurveMode.Force: {'name':'Force','unit':'N','default_prefix':-9},
	CurveMode.Extension: {'name':'Extension','unit':'m','default_prefix':-9},
	CurveMode.Separation: {'name':'Separation','unit':'m','default_prefix':-9}}
prefix = {-15:'f',-12:'p',-9:'n',-6:u'\u03BC',-3:'m',
		0:'',3:'k',6:'M',9:'G',12:'T',15:'P'}

# Generate string for units (bearing in mind raw data scale by adding onto prefix)
def makeUnitStr(mode,pfix=0):
	pfix += unit[mode]['default_prefix']
	return prefix[pfix] + unit[mode]['unit']

def makeTypeStr(mode,pfix=0):
	return unit[mode]['name'] + ' (' + makeUnitStr(mode,pfix) + ')'

def yAxisMode(flags):
	return flags & (CurveMode.Volts|CurveMode.Deflection|CurveMode.Force)

def xAxisMode(flags):
	return flags & (CurveMode.Extension|CurveMode.Separation)

# Determine scale of data (e.g. nN or pN or fN)
def decidePrefix(minX,maxX,mode):
	exp = __calcExp(minX,maxX)
	# Find closest multiple of three 
	pfix = exp - exp % 3
	return pfix

def __calcExp(minX,maxX):
	scale = max(abs(maxX),abs(minX))
	if scale == 0:
		return 0
	# Base 10 exponent
	return int(math.floor(math.log10(scale)))

# Check if ID is valid
def isValid(i):
	return i != -1

# Return list of python strings from QStringList
def PythoniseQStringList(qstrlist):
	y = []
	for x in qstrlist:
		y.append(x.toUtf8().data())
	return y