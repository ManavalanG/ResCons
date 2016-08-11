import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
excludes_options = ['matplotlib', 'numpy', 'scipy']
DATA_FILES = ['clustalo_embl_api.py']

# clustal omega API requires module xmltramp2
build_exe_options = {"packages": ["os"], 
					"excludes":excludes_options,
					"includes": ['xmltramp2']
					}


# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "ResCons 1.0",
        version = "1.0",
        description = "ResCons 1.0",
		data_files=DATA_FILES,
        options = {"build_exe": build_exe_options},
        executables = [Executable("ResCons_v1.py",
						base=base, 
						icon= 'ResCons_Files//ResCons_windows_logo.ico'
						)])
		