# Curves force analysis software
# JWA Copyright 19/04/08

# Model fitting dock widget
from PyQt4 import Qt

#import util
import dock_models
import stats
import curveplot
from forcecurve import CurveMode

class ModelsDock(Qt.QDockWidget):
	def __init__(self,dataset,parent=None):
		Qt.QObject.__init__(self,parent) # Must call otherwise C++ object deleted
		self.setAttribute(Qt.Qt.WA_DeleteOnClose)
		# Curve manager
		self._dataset = dataset
		# User interface
		self.ui = dock_models.Ui_dockModels()
		self.ui.setupUi(self)

		self.connectSlots()
		
	def connectSlots(self):
		pass
	
	
	@Qt.pyqtSignature('') # Zero argument version
	def on_btnFitModels_clicked(self):
		if self.ui.rbPeaks.isChecked():
			# List of marked curves
			curve_indices = self._dataset.markedCurves()
			curves = self._dataset.getCurves(curve_indices)
		elif self.ui.rbCursors.isChecked():
			# Find current activated plot window
			current = self.parent().currentWidget()
			if isinstance(current,curveplot.CurvePlotWindow):
				curves = [current.getCurve()]
				return
			else:
				return
		
		# Check selected model
		
		direction = 'retract'
		w=5 # Skip window
		flags = CurveMode.Force | CurveMode.Extension | CurveMode.Baseline | CurveMode.Offset_Z
		# WLC specific - push specifics out to analyse?
		if self.ui.rbCursors.isChecked():
			# Fit within cursors
			#geta
			#getb
			#getdirection
			(Z,A,R) = x.getData(flags)
			if direction == 'retract':
				Y = R
			else:
				Y = A
			wx = Z[a:b]-Z[a] # start from zero
			wy = -Y[a:b] # invert
			fit = stats.fit(wx,wy,model='wormlikechain')
			fits = [(fit,a,b)]
		else:
			# Fit each peak
			for x in curves:
				if x!=None:
					fits = []
					a = x.contactPointIndex() # Start from here
					p = x.getPeaks(direction)
					(Z,A,R) = x.getData(flags)
					if direction == 'retract':
						Y = R
					else:
						Y = A
					for pk in p:
						b = pk[3] # z index of peak centre
						wx = Z[a:b]-Z[a] # start from zero
						wy = -Y[a:b] # invert
						fit = stats.fit(wx,wy,model='wormlikechain')
						fits.append((fit,a,b))
						a=b+w
					x.setModels(fits)
				
	@Qt.pyqtSignature('') # Zero argument version
	def on_btnShowTable_clicked(self):
		# Get curves
		if self.ui.rbPeaks.isChecked():
			# List of marked curves
			curve_indices = self._dataset.markedCurves()
			curves = self._dataset.getCurves(curve_indices)
			if len(curves)==0:
				return
		elif self.ui.rbCursors.isChecked():
			# Find current activated plot window
			current = self.parent().currentWidget()
			if isinstance(current,curveplot.CurvePlotWindow):
				curves = [current.getCurve()]
				if curves[0] == None:
					return
			else:
				return

		# Create tree window
		tree = self.parent().newTreeWindow(5,'Models')
		tree.setHeaders(['Name','Persistence length','Contour length','Offset','Range'])
		self.populateTree(tree.getModel(),curves)
		tree.getView().expandAll()

	def populateTree(self,tree,curves):
		# Add model data for each curve
		for i,x in enumerate(curves):
			parent = tree.insertRow() # Curve name
			tree.setRow(parent,[x.getName()])
			p = x.getModels()
			for i,v in enumerate(p):
				index = tree.insertRow(parent)
				tree.setRow(index,[i]+list(v[0].beta)+[str(v[1])+':'+str(v[2])])
