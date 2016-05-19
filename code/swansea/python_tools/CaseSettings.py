#!/usr/bin/python
import numpy
import scipy
import sys
import re

import ConfigParser
import os
from vecutils import v2string, s2vector

class CaseSettings:
    """Import settings from a case file and provide access to related files"""
    def __init__(self, configfile = None):
        ##  Selects 2D or 3D simulation. Valid values: 2, 3
        self.mode = 2
        ##  Path to and base name component of input files
        self.basename = ""
        ##  Output folder
        self.outputpath = ""
        ##  Name of file defining particle properties and initial state
        self.particlefile = ""
        ##  Run this number of steps for each input timestep. -1 to disable
        self.ptsteps = 1
        ##  Gradient variable prefix - for autodetection of pre-computed gradients
        self.prefix = "GRAD"
        ##  Variable index for Z coordinate data
        self.z = 0
        ##  Variable index for X component of fluid velocity
        self.u = 1
        ##  Variable index for Y component of fluid velocity
        self.v = 2
        ##  Variable index for Z component of fluid velocity
        self.w = 3
        ##  Variable index for local food availability
        self.fish = -1
        ##  Variable index for noise data
        self.noise = -1
        ##  Variable index for water column depth
        self.depth = -1
        ##  Scaling values for food availability
        self.scaling = {-1, -1, -1}
        ##  Threshold depth for behaviour rules
        self.zthreshold = 5
        ##  Noise threshold
        self.nthreshold = 50

        ## Default Speed
        self.defaultspeed = 1.5
        ## Default random mean
        self.defaultspeedrange = 0.5
        ## Y velocity damping
        self.defaultydamp = 0.5
        ## Default food weighting
        self.foodweight = 1

        ## Depth avoidance speed
        self.depthspeed = 2
        ## Depth avoidance random mean
        self.depthspeedrange = 1
        ## Noise avoidance speed 
        self.noisespeed = 2
        ## Noise avoidance random mean
        self.noisespeedrange = 1

        if not configfile is None:
            self.read(configfile)

    def write(self, filename):
        cf = ConfigParser.ConfigParser()

        cf.add_section("paths")
        cf.set("paths", "basename", self.basename)
        cf.set("paths", "particles", self.particles)
        cf.set("paths", "output", self.outputpath)
        
        cf.add_section("settings")
        cf.set("settings", "dimensions", self.mode)
        cf.set("settings", "particle_steps", self.ptsteps)
        cf.set("settings", "prefix", self.prefix)
        cf.set("settings", "threshold", self.zthreshold)
        cf.set("settings", "noise_threshold", self.nthreshold)
        cf.set("settings", "food_weight", self.foodweight)
        
        cf.add_section("indexes")
        cf.set("indexes", "z", self.z)
        cf.set("indexes", "u", self.u)
        cf.set("indexes", "v", self.v)
        cf.set("indexes", "w", self.w)
        cf.set("indexes", "noise", self.noise)
        cf.set("indexes", "fish", self.fish)
        cf.set("indexes", "depth", self.depth)

        cf.add_section("scaling")
        cf.set("scaling", "fish", v2string(self.scaling))
        
        cf.add_section("speeds")
        cf.set("speeds", "default", self.defaultspeed)
        cf.set("speeds", "defaultrange", self.defaultspeedrange)
        cf.set("speeds", "ydamping", self.defaultydamp)
        cf.set("speeds", "noise", self.noisespeed)
        cf.set("speeds", "noiserange", self.noisespeedrange)
        cf.set("speeds", "depth", self.depthspeed)
        cf.set("speeds", "depthrange", self.depthspeedrange)

        with open(filename,'w') as cfd:
            cf.write(cfd)

    def get_optional(self, cfgetter, section, option, default):
        try:
            return cfgetter(section, option)
        except ConfigParser.NoOptionError:
            return default

    def read(self, filename):
        cf = ConfigParser.ConfigParser()
        cf.read(filename)

        get_optional = self.get_optional

        if not cf.has_section("paths"):
            raise Exception("File has no ""paths"" section")

        try:
            self.basename = cf.get("paths", "basename")
            self.particles = cf.get("paths", "particles")
            self.output = cf.get("paths", "output")
        except Exception as e:
            print e
            return False

        if cf.has_section("settings"):
            self.mode = get_optional(cf.getint, "settings", "dimensions", self.mode)
            self.ptsteps = get_optional(cf.getint, "settings", "particle_steps", self.ptsteps)
            self.prefix = get_optional(cf.get, "settings", "prefix", self.prefix)
            self.zthreshold = get_optional(cf.getfloat, "settings", "threshold", self.zthreshold)
            self.nthreshold = get_optional(cf.getfloat, "settings", "noise_threshold", self.nthreshold)
            self.foodweight = get_optional(cf.getfloat, "settings", "food_weight", self.foodweight)

        if cf.has_section("indexes"):
            self.z = get_optional(cf.getint, "settings", "z", self.z)
            self.u = get_optional(cf.getint, "settings", "u", self.u)
            self.v = get_optional(cf.getint, "settings", "v", self.v)
            self.w = get_optional(cf.getint, "settings", "w", self.w)
            self.noise = get_optional(cf.getint, "settings", "noise", self.noise)
            self.depth = get_optional(cf.getint, "settings", "depth", self.depth)
            self.fish = get_optional(cf.getint, "settings", "fish", self.fish)

        if cf.has_section("scaling"):
            self.scaling = s2vector(cf.get("scaling", "fish"))

        if cf.has_section("speeds"):
            self.defaultspeed = get_optional(cf.getfloat, "speeds", "default", self.defaultspeed)
            self.defaultspeedrange = get_optional(cf.getfloat, "speeds", "defaultrange", self.defaultspeedrange)
            self.defaultydamp = get_optional(cf.getfloat, "speeds", "ydamping", self.defaultydamp)
            self.noisespeed = get_optional(cf.getfloat, "speeds", "noise", self.noisespeed)
            self.noisespeedrange = get_optional(cf.getfloat, "speeds", "noiserange", self.noisespeedrange)
            self.depthspeed = get_optional(cf.getfloat, "speeds", "depth", self.depthspeed)
            self.depthspeedrange = get_optional(cf.getfloat, "speeds", "depthrange", self.depthspeedrange)

    def __str__(self):
        rv = []
        for i in vars(self):
            rv.append(".%s = %s" % (i, getattr(self,i)))
        return "\n".join(rv)

    def validate(self):
        """ Validate mode, ptsteps, indices and check outputpath exists """
        if not (self.mode == 2 or self.mode == 3):
            raise Exception("Invalid number of dimensions")

        if self.ptsteps == 0:
            self.ptsteps = 1

        if self.ptsteps < 1 and not self.ptsteps == -1:
            raise Exception("Invalid number of particle timesteps given")

        if self.z < 0:
            raise Exception("Invalid index specified for Z coordinate data")

        if self.u < 0:
            raise Exception("Invalid index given for u velocity data")

        if self.v < 0:
            raise Exception("Invalid index given for v velocity data")

        if self.w < 0:
            raise Exception("Invalid index given for w velocity data")

        if not (os.path.exists(self.outputpath) and os.path.isdir(self.outputpath)):
            raise Exception("Output path does not exist or is not a directory")

        if self.particlefile == "":
            self.particlefile = None

    def __str__(self):
        rv = []
        for i in vars(self):
            rv.append("self.%s = %s" % (i, getattr(self,i)))
        return "\n".join(rv)

    def getVariables(self):
        """Load variables for this case from files"""

        df = open("%s.vars.txt" % self.basename, "r")
        df.next() #Discard header line
        varline=re.compile("^(\d+?)	(.*)$")
        variables=[]
        for var in df:
            v = varline.match(var)
            vn = v.group(2).strip()
            vn = re.sub(' \s+', ' ', vn)
            vi = int(v.group(1))

            variables.append((vi,vn))

        return variables

    def getTimes(self):
        """Load mesh timestep details for this case from file"""

        tf = open("%s.times.txt" % self.basename, "r")
        tf.next()
        timeline = re.compile("^(\d+?)\s([-+]?(\d+(\.\d*)?|\.\d+))$")
        times = []
        for time in tf:
            t = timeline.match(time)
            tv = float(t.group(2).strip())
            ti = int(t.group(1).strip())
            times.append((ti,tv))

        return times
