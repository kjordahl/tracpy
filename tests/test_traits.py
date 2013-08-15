"""Testing results of new traits model to check that results are
consistent with run.py

Should be run with "nosetests" to ensure relative imports of
tracpy modules work correctly.


Kelsey Jordahl
Enthought
Time-stamp: <Thu Aug 15 19:40:08 EDT 2013>
"""

import unittest
import numpy as np
from datetime import datetime

# tracpy modules
from .. import inout
from .. import tools
from .. import run
from .. import model_run


class TestTraitsModel(unittest.TestCase):
    """Run a tracpy model with run.py and with a Traits model to
    ensure that results are consistent
    """

    def setUp(self):
        self.loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        self.ndays = 5
        self.date = datetime(2009, 11, 25, 0)
        self.tseas = 4 * 3600
        self.nsteps = 5
        self.ff = 1
        self.ah = 0.
        self.av = 0.
        self.doturb = 0
        self.do3d = 0
        self.grid = inout.readgrid(self.loc)
        self.lon0 = np.linspace(-98.5, -87.5, 55)
        self.lat0 = np.linspace(22.5, 31, 49)
        lon0, lat0 = np.meshgrid(self.lon0, self.lat0)
        lon0, lat0 = tools.check_points(lon0, lat0, self.grid)
        self.z0 = 's'
        self.zpar = 29
        self.name = 'Test run.py'
        # TODO: this result could be stored to save time!
        # (but only if run.py is not changing, too)
        self.lonp, self.latp, self.zp, self.t, self.grid = run.run(self.loc, self.nsteps,
                                                                   self.ndays, self.ff, self.date,
                                                                   self.tseas, self.ah, self.av,
                                                                   lon0, lat0,
                                                                   self.z0, self.zpar, self.do3d,
                                                                   self.doturb, self.name)

    def test_model(self):
        model = model_run.ModelRun(lat0=self.lat0, lon0=self.lon0, ndays=self.ndays, loc=self.loc,
                         date=self.date, tseas=self.tseas, nsteps=self.nsteps, ff=self.ff,
                         ah=self.ah, av=self.av, doturb=self.doturb, do3d=self.do3d,
                         z0=self.z0, zpar=self.zpar, name='Test Traits Model')
        model.initialize()
        model.run_steps()
        np.testing.assert_array_almost_equal(self.lonp, model.lonp)
        np.testing.assert_array_almost_equal(self.latp, model.latp)
