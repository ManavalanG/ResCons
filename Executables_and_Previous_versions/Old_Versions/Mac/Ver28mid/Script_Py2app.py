"""
This is a setup.py script generated by py2applet

Usage:
    python setup.py py2app
"""

from setuptools import setup

APP = ['ResCon_Ver28mid_beta.py']
DATA_FILES = ['User_Guide_ResCon.pdf', '128x128_part_no_bg_logo.png']
OPTIONS = {'argv_emulation': True,
			'excludes': ['matplotlib', 'numpy', 'scipy'],
			'includes': ['xmltramp2'],
			'iconfile':'icon_ResCon.icns'
			}

setup(
    app=APP,
    name='ResCon (beta v1.28mid)',
    data_files=DATA_FILES,
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],
)