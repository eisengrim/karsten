#!/usr/bin/python
import ConfigParser
import re
from vecutils import s2vector, v2string

class ParticleFile:
    """Import Particle definitions from file"""

    def __init__(self, pfilepath=None):
        if pfilepath:
            self.read(pfilepath)

    def read(self, pfilepath):
        """Read details from INI style particle file"""
        pfile = ConfigParser.ConfigParser()
        pfile.read(pfilepath)

        if not pfile.has_section("info"):
            raise Exception("Particle file does not include an info section")

        self.np = pfile.getint("info", "np")
        self.comment = ""
        if pfile.has_option("info", "comment"):
            self.comment = pfile.get("info", "comment")

        self.default = Particle()
        if pfile.has_section("defaults"):
            self.default.setPropertiesFromList(pfile.items("defaults"))

        self.particles = []
        for i in range(0, self.np):
            t = Particle(self.default)
            t.id = i
            self.particles.append(t)

        for sect in pfile.sections():
            if re.match("\d+", sect):
                pid = int(sect)
                self.particles[pid].setPropertiesFromList(pfile.items(sect))

    def write(self, pfile):
        """Write details to an INI style particle file"""
        pf = ConfigParser.ConfigParser()

        pf.add_section("info")
        pf.set("info", "comment", self.comment)
        pf.set("info", "np", self.np)

        if (not self.default):
            self.default = Particle()

        pf.add_section("defaults")

        defaultprops = self.default.getProperties()
        for key in iter(defaultprops):
            pf.set("defaults", key, defaultprops[key])

        for i in self.particles:
            pf.add_section(str(i.id))
            p = i.getProperties()
            for key in iter(p):
                if not p[key] == defaultprops[key]:
                    pf.set(str(i.id), key, p[key])
        if isinstance(pfile, str):
            with open(pfile, "w") as pfd:
                pf.write(pfd)
        else:
            pf.write(pfile)


    def __str__(self):
        rv = []
        for i in vars(self):
            rv.append(".%s = %s" % (i, getattr(self,i)))
        return "\n".join(rv)


class Particle:
    """Python equivalent of the particle type defined in particles.h
    At time of writing:
        typedef struct particle {
                int32_t	id;
                double position[3];
                double velocity[3];
                double orientation[3];
                double cent_pressure[3];
                double coeff_drag[3];
                double coeff_lift;
                double area[3];
                double mass;
                int lastelement;
                int exited;
        } particle;
    """
    ## Maps INI file headings to property names and specifies conversion function
    # Format: ini_name: [attribute, decode, encode]
    # decode = function to convert input string to data
    # encode = function to convert data to string
    props = {
            "position": ["position", s2vector, v2string],
            "velocity": ["velocity", s2vector, v2string],
            "orientation": ["orientation", s2vector, v2string],
            "cop": ["cent_pressure", s2vector, v2string],
            "drag": ["coeff_drag", s2vector, v2string],
            "lift": ["coeff_lift", float, str],
            "area": ["area", s2vector, v2string],
            "mass": ["mass", float, str],
            }

    def __init__(self, defaults=None):
        """
        Create a new particle, optionally copying properties from another particle
        Note that the particle ID, lastelement and exited properties are not
        copied from the template particle.
        """

        self.id = -1
        self.position = [0, 0, 0]
        self.velocity = [0, 0, 0]
        self.orientation = [0, 0, 0]
        self.cent_pressure = [0, 0, 0]
        self.coeff_drag = [0, 0, 0]
        self.coeff_lift = 0
        self.area = 0
        self.mass = 0
        self.lastelement = -1
        self.exited = -1

        if defaults:
            if not isinstance(defaults, Particle):
                raise TypeError("Default options provided did not represent a valid Particle")
            self.position = defaults.position[:]
            self.velocity = defaults.velocity[:]
            self.orientation = defaults.orientation[:]
            self.cent_pressure = defaults.cent_pressure[:]
            self.coeff_drag = defaults.coeff_drag
            self.coeff_lift = defaults.coeff_lift
            self.area = defaults.area
            self.mass = defaults.mass

    def setPropertiesFromList(self, proplist):
        """
        Set particle properties from a list of (name, value) pairs

        Expects the list to come from ConfigParser.items(), for example from
        a ParticleFile instance
        """
        props = Particle.props
        for (name, value) in proplist:
            setattr(self, props[name][0], props[name][1](value))

    def getProperties(self):
        """
        Return properties (as defined in Particle.props) as a {name: value} dictionary
        """
        props = Particle.props
        retlist = {}
        for name in props.iterkeys():
            tmp = props[name][2](getattr(self, props[name][0]))
            retlist[name]=tmp
        return retlist

    def __str__(self):
        rv = []
        for i in vars(self):
            rv.append(".%s = %s" % (i, getattr(self,i)))
        return "{%s}" % ", ".join(rv)
