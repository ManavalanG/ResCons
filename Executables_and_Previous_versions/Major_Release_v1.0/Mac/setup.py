"""
This is a setup.py script generated by py2applet

Usage:
    python setup.py py2app
"""

from setuptools import setup

APP = ['ResCons_v1.py']
DATA_FILES = ['User_Guide_ResCons.pdf', '128x128_part_no_bg_logo.png']
OPTIONS = {'argv_emulation': True,
			'excludes': ['matplotlib', 'numpy', 'scipy'],
			'includes': ['xmltramp2'],
			'iconfile':'Icon_Mac_ResCons.icns'
			}

setup(
    app=APP,
    name='ResCons 1.0',
    data_files=DATA_FILES,
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],
)
