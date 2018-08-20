# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:17:10 2016

@author: steven
"""


from .core import *

import pkg_resources

__version__ = pkg_resources.require("genmechanics")[0].version