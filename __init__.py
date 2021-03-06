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
		self.spin_reso.setRange(0.1,1.0);
		self.spin_reso.setSingleStep(0.05);
		self.spin_reso.setValue(0.7);
		glayout1.addWidget(self.spin_reso,1,1);
		
		self.label_imagesize = QtWidgets.QLabel(self);
		self.label_imagesize.setText("Image Size");
		glayout1.addWidget(self.label_imagesize,2,0);
		self.spin_imagesize = QtWidgets.QSpinBox(self);
		self.spin_imagesize.setRange(1000,10000);#PILLOW？は 13000 くらいが Max っぽい
		self.spin_imagesize.setSingleStep(100);
		self.spin_imagesize.setValue(4000);
		glayout1.addWidget(self.spin_imagesize,2,1);

		self.label_fontsize = QtWidgets.QLabel(self);
		self.label_fontsize.setText("Font Size");
		glayout1.addWidget(self.label_fontsize,3,0);
		self.spin_fontsize = QtWidgets.QSpinBox(self);
		self.spin_fontsize.setRange(3,100);
		self.spin_fontsize.setSingleStep(1);
		self.spin_fontsize.setValue(20);
		glayout1.addWidget(self.spin_fontsize,3,1);

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

		self.check_tile = QtWidgets.QCheckBox('Repeating Tile');
		
		self.check_tile.clicked.connect(self.checkTileOn);
		self.check_tile.setChecked(False);
		glayout2.addWidget(self.check_tile,1,1);
		
		self.check_oneletter = QtWidgets.QCheckBox('One Letter');
		glayout2.addWidget(self.check_oneletter,2,0);
		
		self.check_nochainname = QtWidgets.QCheckBox('No Chain Name');
		glayout2.addWidget(self.check_nochainname,2,1);
		
		self.check_ignore_occupancy = QtWidgets.QCheckBox('Ignore Occupancy');
		glayout2.addWidget(self.check_ignore_occupancy,3,0);
		self.check_cartoon = QtWidgets.QCheckBox('Cartoon');
		self.check_cartoon.clicked.connect(self.checkCartoonOn);
		glayout2.addWidget(self.check_cartoon,3,1);

		self.check_mmcif = QtWidgets.QCheckBox('Use MMCIF');
		glayout2.addWidget(self.check_mmcif,4,0);

		self.check_builtin = QtWidgets.QCheckBox('Built-in Surface Generator');
		glayout2.addWidget(self.check_builtin,4,1);
		self.check_builtin.clicked.connect(self.checkBuiltinOn);
		
		#MMCIF は AUTH が不完全だ！
		#self.check_label = QtWidgets.QCheckBox('ID Label');
		#self.layout.addWidget(self.check_label);
		
		# Text field for output file
		glayout4 = QtWidgets.QVBoxLayout();
		
		self.label_outprefix = QtWidgets.QLabel('Output Prefix (Prefix+<something> will be overwritten.)');
		glayout4.addWidget(self.label_outprefix);
		
		glayout4b = QtWidgets.QGridLayout();
		self.text_outprefix = QtWidgets.QLineEdit(self);
		self.text_outprefix.setReadOnly(True);
		glayout4b.addWidget(self.text_outprefix,0,0);
		
		self.button_outprefix = QtWidgets.QPushButton(self);
		self.button_outprefix.setText("Save As");
		
		self.button_outprefix.clicked.connect(self.getFile);
		self.text_outprefix.setStyleSheet("background-color: lightgray;");

		glayout4b.addWidget(self.button_outprefix,0,1);
		glayout4.addLayout(glayout4b);

		glayout3 = QtWidgets.QGridLayout();
		self.button_ok = QtWidgets.QPushButton('Create');
		self.button_ok.clicked.connect(self.runSurfStamp);
		glayout3.addWidget(self.button_ok,0,0);

		self.button_close = QtWidgets.QPushButton('Close')
		self.button_close.clicked.connect(self.hide);
		glayout3.addWidget(self.button_close,0,1);
		
		self.layout.addLayout(glayout1);
		self.layout.addLayout(glayout2);
		self.layout.addLayout(glayout4);
		self.layout.addLayout(glayout3);

		screengeom = QtWidgets.qApp.desktop().screenGeometry();

		wwidth = 300;
		hheight = 200;
		self.setGeometry(screengeom.width()/2-wwidth/2,screengeom.height()/2-hheight/2,wwidth,hheight);
		self.setWindowTitle('SurfStamp');
		self.checkTileOn();
		self.checkBuiltinOn();
		self.checkCartoonOn();
		self.show();
		
	def getFile(self):
		filename = QtWidgets.QFileDialog.getSaveFileName(self,"Output Prefix"
		,""
		,"Wavefront OBJ File (*.obj)");
		if filename:
			if filename[0]:
				self.text_outprefix.setText(filename[0]);

	def checkTileOn(self):
		if self.check_tile.isChecked():
			self.spin_fontsize.setEnabled(True);
		else:
			self.spin_fontsize.setEnabled(False);

	def checkCartoonOn(self):
		if self.check_cartoon.isChecked():
			self.check_builtin.setEnabled(False);
			self.spin_reso.setEnabled(False);
		else:
			self.check_builtin.setEnabled(True);
			self.checkBuiltinOn();

	def checkBuiltinOn(self):
		if self.check_builtin.isChecked():
			self.spin_reso.setEnabled(True);
		else:
			self.spin_reso.setEnabled(False);


	def runSurfStamp(self):
		
		import threading
		self.label_message.setText("Please wait...");
		print("Please wait...")
		self.button_ok.setEnabled(False);
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
		

		if len(re.sub("[\s]*","",self.text_outprefix.text())) > 0:
				tmp_outfile = self.text_outprefix.text();
		
		my_view = cmd.get_view();
		
		modelname = self.combo_model.currentText();
		
		output_modelname = re.sub("[^A-Za-z0-9\\.\\-]","_",modelname)+"_obj";
		output_modelname = cmd.get_unused_name(output_modelname);

		surf_args = ["java","-jar",surfstamp_jar];
		
		surf_args.extend(["-out",tmp_outfile]);
		tmpmodel_created = "";
		
		if self.check_cartoon.isChecked() or (not self.check_builtin.isChecked()):
			#Use .obj generated by PyMOL. If this step is skipped, built-in generator is used.
			tmpmodel_created = cmd.get_unused_name("tmpmodel_");
			cmd.load_model(cmd.get_model(modelname),tmpmodel_created);
			
			cmd.set_view(
			(0.9999905824661255, -0.00367919635027647, -0.002306032460182905
			, 0.003680833615362644, 0.9999929666519165, 0.0007080769282765687
			, 0.0023034177720546722, -0.0007165365968830884, 0.999997079372406
			, 0.0, 0.0, -50.0
			, 0.0, 0.0, 0.0, 40.0, 100.0, -20.0));



			activeobjects =  cmd.get_names('public_objects', enabled_only=1);
			for aa in list(activeobjects):
				cmd.disable(aa);
			cmd.enable(tmpmodel_created);
			modelname = tmpmodel_created;
			cmd.hide("everything",modelname);

			if self.check_cartoon.isChecked():
				cmd.show("cartoon",modelname);
				surf_args.extend(["-use_ca","-force","-sep_block"]);
			else:
				cmd.show("surface",modelname);
				surf_args.extend(["-force"]);

			cmd.reset();
			cmd.origin(position=[0.0,0.0,0.0]);
			cmd.center(origin = 0);
			unusedname = cmd.get_unused_name("pseudo_");
			unused_selectionname = cmd.get_unused_name("pseudo_sel_");
			cmd.pseudoatom(unusedname,pos=[0,0,0]);
			cmd.select(unused_selectionname,"/pseudo_//P/PSD`1/PS1");
			cmd.center(selection=unused_selectionname);
			cmd.save(tmpdir.name+"/tmpin.obj",modelname);
			cmd.delete(unusedname);
			cmd.delete(unused_selectionname);

			for aa in list(activeobjects):
				cmd.enable(aa);

			surf_args.extend(["-obj",tmpdir.name+"/tmpin.obj"]);
			cmd.hide("everything",modelname);
		if self.check_mmcif.isChecked():
			tmp_infile = tmpdir.name+"/tmpin.cif";
			surf_args.extend(["-mmcif_use_label","-mmcif",tmp_infile]);
			cmd.save(tmp_infile,modelname);
		else:
			tmp_infile = tmpdir.name+"/tmpin.pdb";
			cmd.save(tmp_infile,modelname);
			surf_args.extend(["-pdb",tmp_infile]);
		
		if len(tmpmodel_created) > 0:
			cmd.delete(tmpmodel_created);

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
			if self.check_cartoon.isChecked():
				surf_args.extend(["-tile","-font_size",str(self.spin_fontsize.value())]);
			else:
				surf_args.extend(["-tile","-no_sep","-font_size",str(self.spin_fontsize.value())]);
			
		if self.check_colorall.isChecked():
			surf_args.extend(["-color_missing","-color_chainbreak"]);
		if self.check_ignore_occupancy.isChecked():
			surf_args.extend(["-ignore_occupancy"]);

		surf_args.extend(["-quiet","1"]);
		import subprocess;
		process = subprocess.run(surf_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT);
		#process = subprocess.run(surf_args, stdout=None, stderr=subprocess.PIPE);
		print(process.stdout.decode("utf-8"));
		#os.system(" ".join(surf_args)+" >&2 ");
		
		self.label_message.setText("Finished.");
		print("Finished.");
		self.text_outprefix.setText("");
		self.button_ok.setEnabled(True);
		self.update();


		if not re.search("\.obj$",tmp_outfile):
			tmp_outfile = tmp_outfile+".obj";
		cmd.load_callback(pymol_obj_loader.myOBJCallback(tmp_outfile),output_modelname);
		cmd.set_view(my_view);
		
		

def __init_plugin__(app=None):
	'''
	Add an entry to the PyMOL "Plugin" menu
	'''
	from pymol.plugins import addmenuitemqt
	addmenuitemqt('SurfStamp', run_surfstamp_plugin_gui)


# global reference to avoid garbage collection of our dialog
surf_dialog = None


def run_surfstamp_plugin_gui():
	'''
	Open our custom dialog
	'''
	global surf_dialog

	if surf_dialog is None:
		surf_dialog = make_surfstamp_dialog();
		surf_dialog.setWindowFlags( surf_dialog.windowFlags() | QtCore.Qt.WindowStaysOnTopHint );
	surf_dialog.combo_model.clear();
	
	surf_dialog.combo_model.addItems(cmd.get_object_list("(all)"));
	
	try :
		sll = cmd.get_names("selections");
		if len(sll) > 0:
			surf_dialog.combo_model.addItems(sll);
	except:
		print("",end="");
	surf_dialog.text_outprefix.setText("");
	surf_dialog.button_ok.setEnabled(True);
	surf_dialog.show();


def make_surfstamp_dialog():
	return SurfStampFrame();
