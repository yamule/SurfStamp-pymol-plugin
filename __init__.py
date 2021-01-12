'''
PyMOL Demo Plugin

The plugin resembles the old "Rendering Plugin" from Michael Lerner, which
was written with Tkinter instead of PyQt.

(c) Schrodinger, Inc.

License: BSD-2-Clause

SurfStamp PyMOL Plugin
(c) yamule
License: Apache License 2.0

'''

from __future__ import absolute_import
from __future__ import print_function
from . import pymol_obj_loader

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

import os;
import re;
import sys;

# entry point to PyMOL's API
from pymol import cmd

# pymol.Qt provides the PyQt5 interface, but may support PyQt4
# and/or PySide as well
from pymol.Qt import QtWidgets
from pymol.Qt import *
from pymol.Qt.utils import loadUi
from pymol.Qt.utils import getSaveFileNameWithExt

surfstamp_jar = os.path.join(os.path.dirname(__file__), 'SurfStamp.jar')

#https://www.finddevguides.com/Pyqt-quick-guide
#https://qiita.com/Nobu12/items/6248c509401b0e666a55
class SurfStampFrame(QtWidgets.QWidget):
	def __init__(self, parent=None):
		super(SurfStampFrame, self).__init__(parent)
		
		self.layout = QtWidgets.QVBoxLayout()
		self.setLayout(self.layout)
		
		
		glayout1 = QtWidgets.QGridLayout();
		self.label_message = QtWidgets.QLabel(self);
		self.label_message.setText("SurfStamp PyMOL plugin");
		self.layout.addWidget(self.label_message);


		self.combo_model = QtWidgets.QComboBox();
		self.combo_model.addItems([]);
		self.layout.addWidget(self.combo_model);
		
		

		self.label_reso = QtWidgets.QLabel(self);
		self.label_reso.setText("Surface Resolution");
		glayout1.addWidget(self.label_reso,1,0);
		self.spin_reso = QtWidgets.QDoubleSpinBox (self)
		self.spin_reso.setRange(0.1,2.0);
		self.spin_reso.setSingleStep(0.05);
		self.spin_reso.setValue(0.7);
		glayout1.addWidget(self.spin_reso,1,1);
		
		self.label_imagesize = QtWidgets.QLabel(self);
		self.label_imagesize.setText("Image Size");
		glayout1.addWidget(self.label_imagesize,2,0);
		self.spin_imagesize = QtWidgets.QSpinBox(self);
		self.spin_imagesize.setRange(500,15000);
		self.spin_imagesize.setSingleStep(10);
		self.spin_imagesize.setValue(4000);
		glayout1.addWidget(self.spin_imagesize,2,1);
		
		
		
		glayout2 = QtWidgets.QGridLayout();
		self.check_outline = QtWidgets.QCheckBox('Outline')
		self.check_outline.setChecked(True);
		glayout2.addWidget(self.check_outline,0,0);
		
		self.check_nowater = QtWidgets.QCheckBox('Remove Waters')
		self.check_nowater.setChecked(True);
		glayout2.addWidget(self.check_nowater,0,1);
		
		self.check_colorall = QtWidgets.QCheckBox('Color All')
		self.check_colorall.setChecked(False);
		glayout2.addWidget(self.check_colorall,1,0);

		self.check_tile = QtWidgets.QCheckBox('Repeating Tile')
		self.check_tile.setChecked(False);
		glayout2.addWidget(self.check_tile,1,1);
		
		self.check_oneletter = QtWidgets.QCheckBox('One Letter');
		glayout2.addWidget(self.check_oneletter,2,0);
		
		self.check_nochainname = QtWidgets.QCheckBox('No Chain Name');
		glayout2.addWidget(self.check_nochainname,2,1);
		
		self.check_ignore_occupancy = QtWidgets.QCheckBox('Ignore Occupancy');
		glayout2.addWidget(self.check_ignore_occupancy,3,0);
		
		
		#MMCIF は AUTH が不完全だ！
		#self.check_label = QtWidgets.QCheckBox('ID Label');
		#self.layout.addWidget(self.check_label);
		
		
		glayout3 = QtWidgets.QGridLayout();
		self.okButton = QtWidgets.QPushButton('Create');
		self.okButton.clicked.connect(self.runSurfStamp);
		glayout3.addWidget(self.okButton,0,0);

		self.closeButton = QtWidgets.QPushButton('Close')
		self.closeButton.clicked.connect(self.hide);
		glayout3.addWidget(self.closeButton,0,1);
		
		self.layout.addLayout(glayout1);
		self.layout.addLayout(glayout2);
		self.layout.addLayout(glayout3);

		screengeom = QtWidgets.qApp.desktop().screenGeometry();

		wwidth = 300;
		hheight = 200;
		self.setGeometry(screengeom.width()/2-wwidth/2,screengeom.height()/2-hheight/2,wwidth,hheight)
		self.setWindowTitle('SurfStamp')
		self.show()
		
		
	def runSurfStamp(self):
		
		import threading
		self.label_message.setText("Please wait...");
		self.okButton.setEnabled(False);
		self.update();
		thread1 = threading.Thread(target=self.runSurfStamp_);
		thread1.start();
		#thread1.join();
		#self.label_message.setText("Finished.");
		#self.update();

	def runSurfStamp_(self):
		import tempfile;
		tmpdir = tempfile.TemporaryDirectory();
		tmp_outfile = tmpdir.name+"/tmpout.obj";
		
		my_view = cmd.get_view()
		
		modelname = self.combo_model.currentText();
		
		usednames = cmd.get_names();
		
		usednames_hs={};
		for mm in usednames:
			usednames_hs[mm] = 100;
		output_modelname = re.sub("[^A-Za-z0-9\\.\\-]","_",modelname)+"_obj";
		cou = 1;
		while output_modelname in usednames_hs:
			output_modelname = modelname+"_obj"+"_"+str(cou);
			#cmd.get_unused_name("XXX"); というのもある？
			cou += 1;

		surf_args = ["java","-jar",surfstamp_jar];
		
		surf_args.extend(["-out",tmp_outfile]);
		
		if False:
			tmp_infile = tmpdir.name+"/tmpin.cif";
			surf_args.extend(["-mmcif_use_label",tmp_infile]);
			cmd.save(tmp_infile,modelname);
		else:
			tmp_infile = tmpdir.name+"/tmpin.pdb";
			cmd.save(tmp_infile,modelname);
			
			surf_args.extend(["-pdb",tmp_infile]);
			
		surf_args.extend(["-surface_resolution",str(self.spin_reso.value())]);
		surf_args.extend(["-image_size",str(self.spin_imagesize.value())]);
		
		
		
		if not self.check_outline.isChecked():
			surf_args.extend(["-nooutline"]);
		if self.check_nowater.isChecked():
			surf_args.extend(["-nowater"]);
		if self.check_oneletter.isChecked():
			surf_args.extend(["-residue_oneletter"]);
		if self.check_nochainname.isChecked():
			surf_args.extend(["-nochainname"]);
		if self.check_tile.isChecked():
			surf_args.extend(["-tile","-no_sep"]);
		if self.check_colorall.isChecked():
			surf_args.extend(["-color_missing","-color_chainbreak"]);
		if self.check_ignore_occupancy.isChecked():
			surf_args.extend(["-ignore_occupancy"]);


		#import subprocess;
		#process = subprocess.run(surf_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT);
		#print(process.stdout.decode("utf-8"));
		os.system(" ".join(surf_args)+" >&2 ");
		
		self.label_message.setText("Finised.");
		self.okButton.setEnabled(True);
		self.update();

		cmd.load_callback(pymol_obj_loader.myOBJCallback(tmp_outfile),output_modelname);
		cmd.set_view(my_view);
		
		

def __init_plugin__(app=None):
	'''
	Add an entry to the PyMOL "Plugin" menu
	'''
	from pymol.plugins import addmenuitemqt
	addmenuitemqt('SurfStamp', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
surf_dialog = None


def run_plugin_gui():
	'''
	Open our custom dialog
	'''
	global surf_dialog

	if surf_dialog is None:
		surf_dialog = make_dialog();
	surf_dialog.combo_model.clear();
	
	surf_dialog.combo_model.addItems(cmd.get_object_list("(all)"));
	
	try :
		sll = cmd.get_names("selections");
		if len(sll) > 0:
			surf_dialog.combo_model.addItems(sll);
	except:
		print("",end="");

	surf_dialog.show();


def make_dialog():
	return SurfStampFrame();
