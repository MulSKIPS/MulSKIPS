"""
pymulskips
"""

__author__ = "Gaetano Calogero, Giuseppe Fisicaro, Ioannis Deretzis, Antonino La Magna"
# __license__ = "MPL-2.0"

#from . import io
#from . import process
#from . import setuprun
#from . import analyze

import glob, os
print("""
  *** Welcome to PyMulSKIPS ***
  Modules:""")
path = os.path.dirname(os.path.realpath(__file__))
for f in glob.glob(path+'/*.py'):
    if f != path+'/__init__.py':
        print('    '+os.path.splitext(os.path.basename(f))[0])
print('')
