"""A simple CherryPy server to set a drifter location and return a track

Kelsey Jordahl
Enthought, Inc.
Time-stamp: <Tue Aug 20 16:51:51 EDT 2013>
"""

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

    @cherrypy.expose
    def drifter(self, location=None):
        """Set a new drifter location and return a track"""
        lat, lng = location.split(',')
        r, c = self.lonp.shape
        n = r // 2
        lonp = self.lonp[n, ~np.isnan(self.lonp[n])]
        latp = self.latp[n, ~np.isnan(self.latp[n])]
        # return a dummy track starting at the point clicked
        lonp = lonp - lonp[0] + float(lng)
        latp = latp - latp[0] + float(lat)
        line = LineString(zip(lonp, latp))
        return simplejson.dumps(line.__geo_interface__)


if __name__ == '__main__':
    lat0 = np.arange(23, 31, 0.5)
    lon0 = np.arange(-98, -87, 0.5)
    print 'creating model...'
    model = WebModelRun(lat0=lat0, lon0=lon0, ndays=10)
    print 'initializing model...'
    model.initialize()
    model.run_steps()
    conf = {'/': {'tools.staticdir.root': os.path.dirname(os.path.abspath(__file__))},
            '/tracker': {'tools.staticdir.on': True,
                          'tools.staticdir.dir': 'static'}}
    print 'starting web server...'
    cherrypy.quickstart(model, '/', config=conf)
