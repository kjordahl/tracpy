"""A simple CherryPy server to set a drifter location and return a track

Kelsey Jordahl
Enthought, Inc.
Time-stamp: <Thu Aug 22 13:48:02 EDT 2013>
"""

import tracmass

import os
import cherrypy
import simplejson

import numpy as np
from shapely.geometry import LineString

from model_run import ModelRun


class WebModelRun(ModelRun):
    @cherrypy.expose
    def index(self):
        return 'CherryPy server here'

    def find_closest_track(self, lon, lat):
        """Find the closest model track to clicked point"""
        r = np.hypot(self.lonp[:, 0] - lon, self.latp[:, 0] - lat)
        return r.argmin()

    @cherrypy.expose
    def drifter(self, lat=None, lng=None):
        """Set a new drifter location and return a track"""
        lat, lng = float(lat), float(lng)
        n = self.find_closest_track(lng, lat)
        print n
        lonp = self.lonp[n, ~np.isnan(self.lonp[n])]
        latp = self.latp[n, ~np.isnan(self.latp[n])]
        # return a shifted track from nearby starting at the point clicked
        lonp = lonp - lonp[0] + float(lng)
        latp = latp - latp[0] + float(lat)
        line = LineString(zip(lonp, latp))
        feature = {"type": "Feature",
                   "geometry": line.__geo_interface__,
                   "properties": {"start_date": self.date.strftime("%d %b %Y %H:%M"),
                                  "ndays": self.ndays,
                                  "nsteps": self.nsteps,
                                  "tseas": self.tseas,
                                  "doturb": self.doturb}}
        return simplejson.dumps(feature)


if __name__ == '__main__':
    lat0 = np.arange(23, 31, 0.5)
    lon0 = np.arange(-98, -87, 0.5)
    print 'creating model...'
    model = WebModelRun(lat0=lat0, lon0=lon0, ndays=5)
    print 'initializing model...'
    model.initialize()
    model.run_steps()
    rootdir = os.path.dirname(os.path.abspath(__file__))
    cherrypy.server.socket_port = 8888
    cherrypy.server.socket_host = '0.0.0.0'
    conf = {'/': {'tools.staticdir.root': rootdir},
            '/tracker': {'tools.staticdir.on': True,
                          'tools.staticdir.dir': 'static'},
            '/favicon.ico': {'tools.staticfile.on': True,
                          'tools.staticfile.filename': os.path.join(rootdir, 'static', 'favicon.ico')}}
    print 'starting web server...'
    cherrypy.quickstart(model, '/', config=conf)
