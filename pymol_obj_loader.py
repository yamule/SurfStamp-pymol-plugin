
# A simple script to load OBJ file in PyMOL.
# "debug.obj" is the OBJ file made by SurfStamp https://github.com/yamule/SurfStamp-public .
# Because OBJ file is very flexible, the files which are created by other programs may not be loaded.
# Usage:
# 1. Save this script in the same directory with your OBJ.
# 2. Launch PyMOL.
# 3. Load this script from File-> Open (All Files (*))
# 4. Type the following command in the command field (The text field which have PyMOL> label). 
# cmd.load_callback(myOBJCallback("debug.obj"),'unique_name');
# 5. I recommend you to load the structure which is the source of your OBJ file. Because this script does not change the camera position.
# I have run this script with PyMOL open source 2.5.0a0.

#Very new python is needed.
#You cannot run with 2.X, at least.
#I recommend you to use some virtual env like conda.
#
#In addition several uncommon module is needed.
#conda install numpy
#conda install pyopengl
#conda install pyopengl-accelerate
#conda install pillow


'''
   Copyright 2021 yamule

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
'''


from OpenGL.GL import *
from PIL import Image
import re;
import time;
import numpy as np;
import gc;


re_l = re.compile("[\r\n]");
re_pt = re.compile("[\s]+");
re_s = re.compile("/");

class MTLData:
	def __init__(self,pdir):
		self.name = None;
		self.ka = [1.0,1.0,1.0];
		self.kd = [0.8,0.8,0.8];
		self.ks = [0.0,0.0,0.0];
		self.ke = [0.2,0.2,0.2];
		self.ns = 0.0;
		self.ni = 1.0
		self.d = 1.0;
		self.illum = 2;
		self.map_kd = None;
		self.map_ks = None;
		self.map_ka = None;
		self.map_kd_w = None;
		self.map_ks_w = None;
		self.map_ka_w = None;
		self.map_kd_h = None;
		self.map_ks_h = None;
		self.map_ka_h = None;
		self.parent_dir = pdir;
		
	def set_ka(self,flis):
		self.ka = [float(flis[0]),float(flis[1]),float(flis[2])];
		
	def set_ks(self,flis):
		self.ks = [float(flis[0]),float(flis[1]),float(flis[2])];
		
	def set_kd(self,flis):
		self.kd = [float(flis[0]),float(flis[1]),float(flis[2])];
		
	def set_ke(self,flis):
		self.ke = [float(flis[0]),float(flis[1]),float(flis[2])];
		
	def set_ns(self,f):
		self.ns = float(f);
		
	def set_ni(self,f):
		self.ni = float(f);
		
	def set_d(self,f):
		self.d = float(f);
		
	def set_illum(self,i):
		self.d = int(i);
		
	def set_map_kd(self,s):
		imm = Image.open(self.parent_dir+"/"+s);
		w,h = imm.size;
		self.map_kd_w = w;
		self.map_kd_h = h;
		self.map_kd = imm.tobytes();
		
	def set_map_ks(self,s):
		imm = Image.open(self.parent_dir+"/"+s);
		w,h = imm.size;
		self.map_ks_w = w;
		self.map_ks_h = h;
		self.map_ks = imm.tobytes();
		
		
	def set_map_ka(self,s):
		imm = Image.open(self.parent_dir+"/"+s);
		w,h = imm.size;
		self.map_ka_w = w;
		self.map_ka_h = h;
		self.map_ka = imm.tobytes();
		
	def set_name(self,s):
		self.name = s;
	
	@staticmethod
	def from_list(llist,pdir):
		mtll = MTLData(pdir);
		skippednum = 0;
		for ll in llist:
			bb = re_pt.split(ll);
			tag = bb.pop(0).lower();
			if tag == "newmtl":
				mtll.set_name(bb[0]);
			elif tag == "ka":
				mtll.set_ka(bb);	
			elif tag == "kd":
				mtll.set_kd(bb);	
			elif tag == "ks":
				mtll.set_ks(bb);	
			elif tag == "ke":
				mtll.set_ke(bb);
			elif tag == "ns":
				mtll.set_ns(bb[0]);
			elif tag == "ni":
				mtll.set_ni(bb[0]);
			elif tag == "d":
				mtll.set_d(bb[0]);
			elif tag == "illum":
				mtll.set_illum(bb[0]);
			elif tag == "map_kd":
				mtll.set_map_kd(re.sub("^[^\s]+[\s]+","",ll));
			elif tag == "map_ks":
				mtll.set_map_ks(re.sub("^[^\s]+[\s]+","",ll));
			elif tag == "map_ka":
				mtll.set_map_ka(re.sub("^[^\s]+[\s]+","",ll));
			else:
				skippednum += 1;
				if len(tag) == 0:
					continue;
				sys.stderr.write(tag+" is not parsed.\n");
		if len(llist) == skippednum:
			return None;
		return mtll;
		
	@staticmethod
	def load(infile):
		ret = [];
		buff = [];
		pdir = re.sub("[\\/]+[^\\/]+$","",infile)+"/";
		with open(infile) as fin:
			for ll in fin:
				ll = re_l.sub("",ll);
				ptt = re_pt.split(ll);
				if ptt[0] == "newmtl":
					if len(buff) > 0:
						ret.append(MTLData.from_list(buff,pdir));
					buff = [];
				buff.append(ll);
		if len(buff) > 0:
			ret.append(MTLData.from_list(buff,pdir));
		ret_hs = {};
		for rr in ret:
			ret_hs[rr.name] = rr;
		return ret_hs;
class OBJData:
	def __init__(self):
		self.img = None;
		self.vt_tri = [];
		self.v_tri = [];
		self.vn_tri = [];
		self.f_tri = [];
		self.mtl = None;
		self.vi_to_vti = {};
	@staticmethod
	def load(objfile):
		
		v=[];
		vt=[];
		vn=[];

		f=[];
		f_tri=[];
		if not re.search("[\\/]",objfile):
			objfile = "./"+objfile;
		parentdir = re.sub("[\\/]+[^\\/]+$","",objfile)+"/";
		mtllib = None;
		ret = [];
		current_obj = OBJData();
		with open(objfile) as fin:
			for ll in fin:
				ll = re_l.sub("",ll);
				ptt = re_pt.split(ll);
				if ptt[0] == "mtllib":
					mtlfilename = re.sub("mtllib +","",ll);
					mtllib = MTLData.load(parentdir+"/"+mtlfilename);
				elif ptt[0] == "v":
					v.append([float(ptt[1]),float(ptt[2]),float(ptt[3])]);
				elif ptt[0] == "vt":
					vt.append([float(ptt[1]),1.0-float(ptt[2])]);
				elif ptt[0] == "vn":
					vn.append([float(ptt[1]),float(ptt[2]),float(ptt[3])]);
				elif ptt[0] == "usemtl":
					if len(current_obj.f_tri) == 0:
						current_obj.mtl = mtllib[ptt[1]];
					else:
						ret.append(current_obj);
						current_obj = OBJData();
						current_obj.mtl = mtllib[ptt[1]];
				elif ptt[0] == "f":
					if len(ptt) == 4 or len(ptt) == 5:
						pf = [];
						for ii in range(1,len(ptt)):
							pp = re_s.split(ptt[ii]);
							
							vii = int(pp[0])-1;
							vtii = -1;
							if len(pp[1]) > 0:
								vtii = int(pp[1])-1;
							vnii = int(pp[2])-1;
							if not vii in current_obj.vi_to_vti:
								current_obj.vi_to_vti[vii] = {};
							if not vtii in current_obj.vi_to_vti[vii]:
								current_obj.v_tri.append([v[vii][0],v[vii][1],v[vii][2]]);
								
								if len(pp[1]) > 0:
									current_obj.vt_tri.append([vt[vtii][0],vt[vtii][1]]);
									
								current_obj.vn_tri.append([vn[vnii][0],vn[vnii][1],vn[vnii][2]]);
								current_obj.vi_to_vti[vii][vtii] = len(current_obj.v_tri)-1;
							current_obj.f_tri.append(current_obj.vi_to_vti[vii][vtii]);
							
					else:
						print("Face must have 3 or 4 vertices.");
		
		if len(current_obj.f_tri) > 0:
			ret.append(current_obj);
		for rr in ret:
			rr.v_tri = np.array(rr.v_tri);
			if len(rr.vt_tri) > 0:
				rr.vt_tri = np.array(rr.vt_tri);
			rr.vn_tri = np.array(rr.vn_tri);
			rr.f_tri = np.array(rr.f_tri);
		del v;
		del vt;
		del vn;
		gc.collect();
		return ret;
		

	

from pymol.callback import Callback
from pymol import cmd

class myOBJCallback(Callback):
	def __init__(self,infilename):
		super().__init__();
		self.all_obj = OBJData.load(infilename);
		
	def show_obj(self,objj):
		#glClear(GL_COLOR_BUFFER_BIT);#後表示のポリゴンが黒くなる
		glColor3d(1.0, 1.0, 1.0);
		glDisableClientState(GL_COLOR_ARRAY);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnable(GL_TEXTURE_2D);

		glNormalPointer(GL_FLOAT, 0, objj.vn_tri);
		glVertexPointer(3, GL_FLOAT, 0,objj.v_tri);
		glTexCoordPointer(2, GL_FLOAT, 0,objj.vt_tri);
		
		glDrawElements(GL_TRIANGLES, len(objj.f_tri), GL_UNSIGNED_INT, objj.f_tri);
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisable(GL_TEXTURE_2D);
		glEnableClientState(GL_COLOR_ARRAY);
		glFlush();
		
	def __call__(self):
		for aa in self.all_obj:
			glClearColor(0.0, 0.0, 0.0, 1.0)
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, aa.mtl.map_ka_w, aa.mtl.map_ka_h, 0, GL_RGBA, GL_UNSIGNED_BYTE, aa.mtl.map_ka);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			self.show_obj(aa);
		
	def get_extent(self):
		return [[0.0,0.0,0.0],[3.0,3.0,3.0]]
