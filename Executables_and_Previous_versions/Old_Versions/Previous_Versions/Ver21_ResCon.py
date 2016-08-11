#!/usr/bin/env python
'''
Copyright (C) <2015>  Manavalan Gajapathy, Joseph D. Ng

ResCon is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 or version 3 of the License.

ResCon is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with ResCon; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

Author's Email id: manavalan.g@gmail.com
'''

__author__ = "Manavalan Gajapathy, Joseph D. Ng"
# __credits__ = ["Joseph Ng"]
__license__ = "Lesser General Public License"
__version__ = "18 beta"
__email__ = "manavalan.g@gmail.com"


from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename, askdirectory
from tkColorChooser import askcolor
import os, errno, logging, ast, shutil, re
import subprocess
import traceback		# to raise exceptions in gui
import threading
from sys import platform as _platform	# to determine OS under use
Image_in_GUI = False		# if logos needed to be shown in GUI
if Image_in_GUI:
	from PIL import Image, ImageTk		# to add image to the GUI
# note: matplotlib.pyplot may get imported when needed. Scroll down to reading settings from text file to find out.


'''
Ver 21 modifications:
1. Removed the feature of replacing comma with '_$' in tool 'Fasta Description/ID extractor'. Instead standardized the
   way csv file gets written when comma or double quote is part of the string itself.

   As a part of above change, disabled the feature of allowing user to change '_$' with any symbol. Also removed feature
   that looks for '_$' when 'filter by seqs' description' tool is used.

2. Standardized csv file writing in tool 'mismatch analyzer'. If comma or double quote is in the cell value, now
   they are formatted in the standardized manner.

3. Added bar chart to html output showing similarity and identity % (using two methods:
	1. used Matplotlib	- successfully done but final software product size increases to > 150 MB.
	2. Used javascript in html5 - successfully implemented. Final software will remain the same size as now!)

'''

# This is to enable actions specific to particular operating system
# win_os = False		# for both Windows and Ubuntu OS
# linux_os = False	# Keep it as 'False'. This is for threading purposes. But threading is currently causing trouble.
# mac_os	= True		# only for Mac OS
# ubuntu_os = False	# only for Ubuntu OS

# win_os = True		# for both Windows and Ubuntu OS
# linux_os = False	# Keep it as 'False'. This is for threading purposes. But threading is currently causing trouble.
# mac_os	= False		# only for Mac OS
# ubuntu_os = False	# only for Ubuntu OS

# This automatically detects the Operating system under use
if _platform == "linux" or _platform == "linux2":	# for linux OS
	ubuntu_os = True
	mac_os = False
	linux_os = False
	win_os = True
elif _platform == "darwin":		# for Mac
	ubuntu_os = False
	mac_os = True
	linux_os = False
	win_os = False

	try:	# prevents app from opening minimized (at least in OSx 10.9). Uses AppleEvents to give focus to Python app.
		os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
	except Exception as e:
		pass

elif _platform == "win32":		# for Windows
	ubuntu_os = False
	mac_os = False
	linux_os = False
	win_os = True

# For logging. Logging module is utilized but sys.stderr is utilized for writing
# into output log file as it enables to record unexpected errors from the program during its execution.
logger_filename = 'Logs_ResCon.txt'

# This part writes log file to C:/temp folder. However if C: is not
# present, then this part will write into a drive that is present in that computer
if win_os and not ubuntu_os:
	if os.path.exists('C:'):
		logger_filename = 'C://temp//Logs_ResCon.txt'
	else:
		drive_letter = map(chr, range(68, 91))
		drive_letter = ['A', 'B'] + drive_letter		# all alphabets except letter C
		for letter in drive_letter:
			drive = letter + ':'
			if os.path.exists(drive):
				logger_filename = letter + '://temp//Logs_ResCon.txt'
				break

	dir_path = os.path.dirname(logger_filename)	# this part creates temp folder if it is not already present
	try:
		os.makedirs(dir_path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:			# if administrative privilages prevents from creating folder
			tkMessageBox.showerror('Error', "Cannot create file '%s'. Check if you have privileges to create this file "
					"in the path specified." %dir_path)
			raise


sys.stderr = open(logger_filename, 'w')		# Logs and errors are directed to stderr which in turn gets logged.


# Following are dummy variables which will be removed when scripting is completed
html_log = logging.getLogger('value')
genbank_log = logging.getLogger('value')
runscript_log = logging.getLogger('value')
blast_log = logging.getLogger('value')
mismatch_log = logging.getLogger('value')
clustal_log = logging.getLogger('value')
gui_log = logging.getLogger('value')
clades_log = logging.getLogger('value')
descr_filter_log = logging.getLogger('value')
idfilter_log = logging.getLogger('value')
header_extractor_log = logging.getLogger('value')
unexpected_error_log = logging.getLogger('value')
user_error_log = logging.getLogger('value')


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)-19s %(levelname)-8s %(message)s',
					datefmt= '%Y-%m-%d %H:%M:%S')

logging_dict = {'gui_log': 'GUI_Main', 'clustal_log': 'Clustal_Alignment', 'mismatch_log': 'Fetch_Mismatch',
				'html_log': 'html_formatting', 'genbank_log': 'genbank2Fasta', 'clades_log': 'Extract_Clades',
				'blast_log': 'blast_filter', 'runscript_log': 'Run_Script_Main', 'descr_filter_log': 'Filter_by_Descrip',
				'idfilter_log': 'Filter_by_id' , 'header_extractor_log': 'Descrip/ID_Extractor',
				'unexpected_error_log': 'Unexpected_error', 'user_error_log': 'User_error'}

for key, value in logging_dict.iteritems():  # defines multiple loggers using dictionary items
	locals()[key] = logging.getLogger(value)


# This part gets settings from a text file (so as to allow user to modify certain settings)
settings_file_name = "Settings_ResCon.txt"
with open(settings_file_name, 'Ur') as settings_data:
	for line in settings_data:
		line = line .replace('\n', '')
		if line != '' and '##' not in line:			# skips blank and description lines
			for index, char in enumerate(line):		# gets name of variable and also gets rest of the line
				if char == ':':
					var_name = line[0:index].replace(' ', '')
					trimmed_1 = line[index+1 ::]
					break

			for index, char in enumerate(trimmed_1):	# removes leading spaces and gets values provided
				if char != ' ':
					trimmed_2 = trimmed_1[index::]
					break

			for no in range(0, len(trimmed_2)):			# removes spaces at end and gets values provided
				if trimmed_2[-(no+1)] != '^':
					if not no:
						trimmed_3 = trimmed_2
					else:
						trimmed_3 = trimmed_2[:-no]
					break

			# settings for Mismatch Analyzer
			if var_name == 'match_color':		# gets color for matching residues
				# match_color = line_list[1]
				match_color = trimmed_3

			elif var_name == 'similar_color':	# gets color for similar residues
				similar_color = trimmed_3

			elif var_name == 'mismatch_color':	# gets color for mismatching residues
				mismatch_color = trimmed_3

			elif var_name == 'Clustalo_Command':	# gets clustalo command
				line_list = line.split(' {')
				clustalo_command_default = trimmed_3

			elif var_name == 'id_delimiter':
				id_delimiter = trimmed_3

			elif var_name == 'aa_set':			# gets amino acid grouping that user likes to use
				trimmed_2 = trimmed_2.replace(' ', '')
				line_list = trimmed_2.split(',')

				aa_set_upper = [x.upper() for x in line_list]	# makes upper and lower case as a set
				aa_set_lower = [x.lower() for x in line_list]

				aa_set = aa_set_upper + aa_set_lower
				aa_set_str = ', '.join(aa_set_upper)

			elif var_name == 'chart_method':	# gets how to draw bar chart in html
				chart_method = trimmed_3


			# settings for GenPept to fasta converter tool
			elif var_name == 'Fasta_ID_length':		# gets fasta id length
				fasta_id_length_default = trimmed_3

			elif var_name == 'connector_id':		# gets symbol that is used to connect different fields
				connector_id = trimmed_2

			elif var_name == 'newick_sym_replace':	# gets symbol to replace newick sensitive symbols
				newick_sym_replace = trimmed_2


			# settings for Fasta Description/ID Extractor
			# elif var_name == 'symbol_replacing_comma':	# gets symbol for use if user don't want comma in output
			# 	symbol_replacing_comma = trimmed_2

			elif var_name == 'regex_desc_extr':		# regular expression to be used for description extractor
				line_list = line.split(':')
				regex_desc_extr = trimmed_2

			elif var_name == 'regex_id_extr':		# regular expression to be used for ID extractor
				regex_id_extr = trimmed_2


			# Settings for Filter fasta by Sequences' Description
			elif var_name == 'regex_desc_fltr':	# regular expression to be used for filtering by partial description
				regex_desc_fltr = trimmed_2

# chart_method = 'chart.js'		# options: 'pro', 'chart.js' or 'chartnew.js'
if chart_method == 'pro':		# variable's value is obtained from settings text file
	import matplotlib.pyplot as plt


# Function to verify existence of a path and create it if non-existent.
# Also raises warning when folder is already present as it may rewrite files in them.
def verify_filepath_exists(path):
	try:
		os.makedirs(path)
	except OSError as exception:
		runscript_log.error('Error that resulted when creating path "%s":\n  %s' % (path, exception))
		if exception.errno != errno.EEXIST:			# if administrative privilages prevents from creating folder
			temp = ("Cannot create file '%s'. Check if you have privileges to create this file "
					"in the path specified." %path)
			runscript_log.error(temp)
			popup_error(temp)
			# raise
			raise_enabler('stop')

		# Mac & Ubuntu produce OSerror 17 when a folder already exists; for windows (8.1), it is 183
		elif exception.errno == (17 or 183):				# warns if folder pre-exists
			temp = ('Folder "%s" already exists. All or some files in it may get overwritten by ResCons. '
					'Do you wish to continue?' % path)
			runscript_log.warning(temp)
			rewrite_or_not = tkMessageBox.askyesno('Warning', temp, default='yes')
			runscript_log.info("User has selected (True = Yes, False = No):  %s" % rewrite_or_not)
			if not rewrite_or_not:
				runscript_log.info("ResCons will stop as user has requested.\n")

				raise_enabler('User opted not to write in selected output folder!')

			else:
				runscript_log.info("ResCons will continue as user has requested.")

	# Warns user if they intend to write output files into 'Program Files' folder in Windows OS.
	if not ubuntu_os and win_os and 'Program Files' in dir_path:
		temp = "Windows OS may not allow you to write results into 'Program Files' folder unless you have " \
			   "administrator privileges. Click 'Yes' to continue and 'No' to cancel."
		user_error_log.warning(temp)
		ans = tkMessageBox.askyesno('Warning', temp, default = 'no')
		if ans:
			user_error_log.info("User clicked 'Yes'. ResCon continues further.")
		else:
			user_error_log.info("User clicked 'No'. ResCon stops here.")
			raise_enabler('stop')

# Function to pop-up an error in gui
def popup_error(string):
	global top
	try:
		tkMessageBox.showerror('Error', string, parent = top)
	except Exception as e:
		tkMessageBox.showerror('Error', string)


# Function that enables displaying errors in gui when unexpected errors come up during execution of script.
# Thanks to stack-overflow for this function's basic skeleton (@Jochen Ritzel)
def show_error(self, *args):
	err = traceback.format_exception(*args)
	err_string = ''			# turning it a string makes it easier to read.
	for x in err:
		err_string += x
	err_string += ('\n')

	if 'raise_enabler' not in err_string:	# error will pop up only when unexpected error occurs
		temp = "Unexpected error occured. Click 'Yes' for more details on this error."
		unexpected_error_log.error(temp)
		unexpected_error_log.error('Unexpected error that resulted: %s' % err_string)
		try:
			# 'parent' parameter will enable to have pop up window on top of respective windows
			ans = tkMessageBox.askquestion('Error', temp, parent = top, default = 'no')
			if ans == 'yes':
				popup_error(err_string)
		except Exception as e:
			ans = tkMessageBox.askquestion('Error', temp, default = 'no')
			if ans == 'yes':
				popup_error(err_string)

	else:								# for error that can be corrected by user
		user_error_log.error('Error that resulted: %s' % err_string)

	raise_enabler('none')


# Re-enable convert/submit job buttons and hide 'processing' label when 'Raise' function needs to be called
# intentionally to avoid further execution of code when an error is noted.
# Also defines what kind of 'raise' function is called
def raise_enabler(str):
	# to re-enable submit button and hide 'processing job label' of root window during unexpected errors
	try:
		global button_run_script
		global processing, processing_clustal
		button_run_script.configure(state=ACTIVE)
		processing.grid_remove()
		processing_clustal.grid_remove()
	except Exception as e:
		pass

	# to re-enable submit buttons and hide 'processing job label' in child window during unexpected errors
	try:
		global Blast_Submit, genbank_Submit, Extract_Submit, descr_extractor_Submit, fasta_filter_by_name_Submit, fasta_filter_by_id_Submit
		global blast_processing, genbank_processing, subtree_processing, descr_extract_processing, filter_name_processing, filter_id_processing

		submit_buttons = [Blast_Submit, genbank_Submit, Extract_Submit, descr_extractor_Submit,
						fasta_filter_by_name_Submit, fasta_filter_by_id_Submit]
		processing_labels = [blast_processing, genbank_processing, subtree_processing, descr_extract_processing,
						  filter_name_processing, filter_id_processing]

		for no in range(0, len(submit_buttons)):
			submit_buttons[no].configure(state=ACTIVE)
			processing_labels[no].grid_remove()

	except Exception as e:
		pass

	if str == 'stop':
		user_error_log.error("Resulted in an error that can be corrected by User. ResCon is ready for next job now!")
		raise Exception("Error can be rectified by user (you)!. Controlled termination by ResCon!")
	if str == 'OSError':
		raise OSError
	elif str == 'none':		# in cases where submit buttons need to be activated but raise fn not to be called
		pass
	else:
		raise Exception(str)


# Function to index each character (residue) in the Reference sequence
def indexing(string, res):
	index_list = []
	for index, char in enumerate(string):
		if char == res:
			index_list.append(index)
	return index_list


# Works in mac and ubuntu but it is not currently implemented as catching exceptions from thread is not resolved yet.
# Function that runs tools in a new thread and shows a 'process under progress' label
# This function keeps gui stable while a tool is underway.
# this function works in mac and ubuntu but not in windows. In windows, when this function is run, tkmessagebox module
# inside the tool getting called poses trouble as tkinter have trouble being thread-safe.
def gen_thread(proc_label, tool_fn):
	xyz = threading.Thread(target=tool_fn, args= ())
	xyz.start()
	proc_label.grid()
	while xyz.isAlive():
		try:
			app.update()
		except Exception as e:
			top.update()
	xyz.join()
	proc_label.grid_remove()


# function to run a function in a new thread; intended to run clustal omega
def new_thread(tool_fn):
	xyz = threading.Thread(target=tool_fn, args= ())
	xyz.start()
	while xyz.isAlive():
		app.update()
	xyz.join()


# Function to perform clustal alignment, if needed.
# Can accept commands from user
#  (note: user has to supply such commands in dictionary format).
def clustal_alignment():
	clustal_log.debug('Clustal Alignment module initiated.')
	global SeqQuery_file
	global Reference
	global Aligned_Filename
	global Output_Path
	global SeqQuery_FileName
	# global clustal_commandline
	global clustal_user_command_string

	# function to verify if a executable program is available in the system through shell.
	# This function was shamelessly copied from Suraj Barkale's answer in stack exchange
	# source: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
	def which(program):
		def is_exe(fpath):
			return os.path.exists(fpath) and os.access(fpath, os.X_OK)

		def ext_candidates(fpath):
			yield fpath
			for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
				yield fpath + ext

		fpath, fname = os.path.split(program)
		if fpath:
			if is_exe(program):
				return program
		else:
			for path in os.environ["PATH"].split(os.pathsep):
				exe_file = os.path.join(path, program)
				for candidate in ext_candidates(exe_file):
					if is_exe(candidate):
						return candidate
		return None

	is_clustalo_present = which('clustalo')			#tests if clustal omega is available in the system
	if is_clustalo_present == None:
		temp = 'It appears Clustal Omega is not available in your system. Click Cancel if you believe this is not true.'
		clustal_log.error(temp)
		response = tkMessageBox.askokcancel('Clustal Omega not installed', temp, default = 'ok')
		clustal_log.info("User's response (True = OK, False = Cancel) :  %s" %response)
		if response:
			clustal_log.info("User chose not to proceed further.")
			raise_enabler('stop')
		else:
			clustal_log.info("User chose to proceed further despite warning on 'unavailability of Clustal Omega")

	try:
		seq_records = list(SeqIO.parse(SeqQuery_file, "fasta"))
	except Exception as e:
		clustal_log.error('Error that resulted when reading sequences file: %s' % e)
		temp = ("Ran into problem reading file: '%s'" % SeqQuery_file)
		clustal_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	if len(seq_records) == 0:
		temp = ('Sequences are not found in your sequences file. Make sure it is in fasta format.'
				'Try again after fixing this problem')
		clustal_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	clustal_log.info("Number of sequences found in '%s' are '%i'" %(SeqQuery_file, len(seq_records)))

	fasta_IDs_in_Input = [seq_records[SeqNo].id for SeqNo in range(len(seq_records))]

	# Verifies if Reference sequence is present in query sequences and adds to it, if non-existent.
	# And then executes clustal omega.
	OutputFile = None
	if Reference.id not in fasta_IDs_in_Input:
		All_Seqs = list(SeqIO.parse(SeqQuery_file, "fasta"))
		All_Seqs.append(Reference)
		OutputFile_temp = (SeqQuery_FileName).replace('.fasta', '') + "_Reference_SeqAdded.fasta"
		OutputFile = Output_Path + OutputFile_temp
		Output_handle = open(OutputFile, 'w')
		SeqIO.write(All_Seqs, Output_handle, 'fasta')
		Output_handle.close()

		clustal_log.info("Reference seq is not present in file: '%s'. \n\tA new file '%s' is created "
						 "with Reference Seq appended to it." % (SeqQuery_file, OutputFile))
		SeqQuery_file = OutputFile
		Aligned_Filename = Output_Path + "Aligned_Clustalo_" + OutputFile_temp.replace('.fasta', '')
		tree_filename = Output_Path + "Tree_" + OutputFile_temp.replace('.fasta', '') + ".newick"

	else:
		clustal_log.info('Reference seq is present in file: %s' % SeqQuery_file)
		Aligned_Filename = Output_Path + "Aligned_Clustalo_" + (SeqQuery_FileName).replace('.fasta', '')
		tree_filename = Output_Path + "Tree_" + (SeqQuery_FileName).replace('.fasta', '') + ".newick"


	clustal_command_temp_1 = {'infile': SeqQuery_file} 	# Infile and outfile command is provided through script; not dependent on user

	# obtain and verify syntax of user-provided Clustal-O command
	try:
		clustal_command_temp_2 = ast.literal_eval(clustal_user_command_string)	# convert string to a dictionary
	except Exception as e:
		temp = 'Syntax error in Clustal-O command entered. Correct it and Try again!'
		clustal_log.error(temp)
		clustal_log.error('Error that resulted: %s' % e)
		popup_error(temp)
		raise_enabler('stop')

	clustal_commandline = {}									# Merge two dictionaries into one to get complete clustal command
	clustal_commandline.update(clustal_command_temp_1)
	clustal_commandline.update(clustal_command_temp_2)
	clustal_log.info("Clustal commandline provided by user: \n\t%s" %clustal_commandline)

	if 'infile' in clustal_command_temp_2:
		temp = "'infile' parameter is not allowed to be part of Clustal-O command. Remove it and try again!"
		clustal_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	# reads user provided clustal omega output format
	outfmt_allowed = ["a2m", "fa", "fasta", "clu", "clustal", "phy", "phylip", "vie", "vienna"]
	if "outfmt" in clustal_commandline:
		if clustal_commandline['outfmt'] in outfmt_allowed:
			file_extension = clustal_commandline['outfmt']
		else:
			temp = ("'outfmt' parameter '%s' provided in Clustal Omega commandline is not allowed here. " \
				   "Choose among the following parameters:\n a2m, fa, fasta, clu, clustal, phy, phylip, vie or vienna" \
				   % clustal_commandline['outfmt'] )
			clustal_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')
	else:
		file_extension = 'fasta'		# unless specified, clustal output will be in fasta format


	# if user requests to change default filename or filepath of aligned file, following will take care of it.
	if "outfile" in clustal_commandline:
		if clustal_commandline['outfile'] != 'default':
			Aligned_Filename = clustal_commandline['outfile']
			temp_path = os.path.dirname(Aligned_Filename)
			verify_filepath_exists(temp_path)
		else:
			Aligned_Filename += ('.' + file_extension)
			clustal_commandline['outfile'] = Aligned_Filename
	else:
		temp = "Parameter 'outfile' is missing in clustalo command provided. It is required to execute clustalo."
		clustal_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	# this removes clustal aligned file if it is pre-existing in the output folder. Done as a work around to check success of clustal omega alignment
	if os.path.exists(Aligned_Filename):
		os.remove(Aligned_Filename)
		clustal_log.debug("Deleted file: '%s'" % Aligned_Filename)

	# if user requests to change default filename or filepath of phylogenetic tree file, following will take care of it.
	if "guidetree_out" in clustal_commandline:
		if clustal_commandline['guidetree_out'] != 'default.newick':
			tree_filename = clustal_commandline['guidetree_out']
			temp_path = os.path.dirname(tree_filename)
			verify_filepath_exists(temp_path)
		else:
			clustal_commandline['guidetree_out'] = tree_filename

	# Run ClustalO with user provided (partially) commandline
	# clustalomega_cline = ClustalOmegaCommandline(infile=SeqQuery_file, outfile=Aligned_Filename, outfmt='clu',
	# 								guidetree_out=tree_filename, auto=True, verbose=True, force=True)

	try:
		clustalomega_cline = ClustalOmegaCommandline(**clustal_commandline)
	except Exception as e:			# This will catch all the major errors
		temp = 'Error in command directed to Clustalo. Check the command entered!'
		clustal_log.error(temp)
		clustal_log.error('Error that resulted: %s' % e)
		popup_error(temp)
		raise_enabler('stop')

	clustal_log.info("Equivalent command prompt command that will be used to"
					 " execute clustal omega alignment: \n\t'%s'" %str(clustalomega_cline))
	temp = 'Clustal alignment will begin now. Warning: It may take few seconds ' \
		   'to several minutes depending on number of sequences provided and your processor.'
	clustal_log.info(temp)
	tkMessageBox.showinfo('Info', temp)


	# execute clustal omega in a new thread. As catching exceptions when using threading is problematic, ResCon
	# looks for absence of clustal aligned output file as a workaround to verify if clustal omega ran successfully.
	# For the above reason, ResCon deletes if a file with same name as clustal output file was already
	# present in user's output folder.
	processing_clustal.grid()
	if linux_os:
		clustalomega_cline()
	else:
		new_thread(clustalomega_cline)
	processing_clustal.grid_remove()

	# This is if clustal omega results in an error.
	if not os.path.exists(Aligned_Filename):
		processing_clustal.grid_remove()
		temp = 'Error executing clustalo command. Possible problems include: \n  ' \
			   '1. Clustal omega is not installed in your computer  or \n  ' \
			   '2. Error in Clustal omega installation or \n  ' \
			   '3. If already installed, correct path may not have been properly specified or \n  ' \
			   '4. If this is Windows OS, Clustal omega may have trouble if Sequences file has lot of sequences.' \
			   " (This error is caused by Clustal Omega application; not ResCon) \n  " \
			   '5. Error in your Sequences file or in its filename. \n  ' \
			   'Note: This list is not inclusive of all problems though.'
		clustal_log.error(temp)
		popup_error(temp)
		button_run_script.configure(state = ACTIVE)		# Re-enables submit job button while processing
		raise_enabler('Clustal Omega resulted in error! So ultimately ResCon had to terminate the job!')

	# deletes temporarily created fasta file that has reference seq appended to it if wasn't originally present
	if OutputFile and os.path.exists(OutputFile):
		os.remove(OutputFile)
		clustal_log.debug("Deleted reference sequence appended fasta file: '%s'" % OutputFile)


# Function that reads clustal alignment file and extracts details of mismatches, if any, at the sites requested
# and writes detailed info into csv and txt files.
def fetch_mismatch():
	global unique_residues_line
	mismatch_log.debug("Module 'Parsing for Aligned residues' initiated.")
	global resi_positions
	global Aligned_Filename
	global Output_Path
	global query_in_Alignment
	global query_site_residues
	global query_site_actual
	global identity_perc_list
	global similarity_perc_list
	global protein_mode

	mismatch_log.debug("Alignment file '%s' will be read." % Aligned_Filename)
	# data_alignment = AlignIO.read(Aligned_Filename, "clustal")
	# clustal_aligned = False
	format_supported = False
	try:
		# this block identifies if alignment is in any of the supported formats
		format_type_list = ['clustal','emboss','fasta','fasta-m10','ig','maf','nexus','phylip','phylip-sequential','phylip-relaxed','stockholm']
		for format_type in format_type_list:
			try:
				data_alignment = AlignIO.read(Aligned_Filename, format_type)
				# if format_type == "clustal":
					# clustal_aligned = True
				mismatch_log.info('Alignment format is detected as "%s" format' % format_type)
				format_supported = True
				break
			except ValueError:
				pass
	except Exception as e:
		mismatch_log.error("Error that resulted: %s" % e)
		temp = ("Could not open file: '\n\t%s'" % Aligned_Filename)
		mismatch_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	if not format_supported:
		temp = "Your Alignment file does not have alignment format supported by ResCon"
		mismatch_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	mismatch_log.debug("Alignment file '%s' was read successfully." % Aligned_Filename)

	# Checks if Reference sequence ID is present in alignment file.
	Reference_index = None
	for serial_no in range(len(data_alignment)):
		if Reference.id == data_alignment[serial_no].id:
			Reference_index = serial_no

	# Raises error if Reference sequence is not present in alignment file.
	if Reference_index == None:
		temp = ('Error: Alignment file "%s" does not have Sequence ID corresponding to your reference sequence. '
				'Fix it and try again!' % Aligned_Filename)
		mismatch_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	mismatch_log.info("Reference sequence's ID is present in alignment file.")

	# Indexes residues (characters) in reference sequence and also for reference sequence in alignment.
	mismatch_log.debug('Begins indexing residues in reference seq and also of reference sequence in alignment')
	ReferenceSeq_Aligned = str(data_alignment[Reference_index].seq)

	# following line creates list of upper case alphabets.
	# Includes non-amino acid chars as well to support nucleotides other than A, T, G and C
	# works for both upper and lower case characters and considers non-amino acid and non-nucleotide residues as mismatch
	aa_list = map(chr, range(65, 91)) + map(chr, range(97,123))
	Dict_aa_Index_Reference_Seq = {k: [] for k in aa_list}
	Dict_aa_Index_Reference_Seq_inAlignment = {k: [] for k in aa_list}
	for aa in aa_list:
		Dict_aa_Index_Reference_Seq_inAlignment[aa] = indexing(ReferenceSeq_Aligned, aa)
		Dict_aa_Index_Reference_Seq[aa] = indexing(Reference.seq, aa)
	mismatch_log.debug('Done indexing residues in Reference seq and also of Reference sequence in alignment')


	# This section is to get appropriate number of 'name components' in the title row.
	# Also this accounts even if the clustal record ids have excessive number of pipes without any info in them.
	aligned_ids_list = []
	pipes_total_list = []
	justified_list = []
	for record in data_alignment:
		aligned_ids_list.append(record.id)
		pipes = (record.id).split(id_delimiter)

		if len(pipes) == 1:
			justified_list.append(len(pipes[0]))

		pipes_temp = pipes
		for no in range(1, len(pipes)+1):
			if pipes[-no] == '':
				pipes_temp = pipes[:(-no)]
			else:
				break
		pipes_total_list.append(len(pipes_temp))
	pipes_most = max(pipes_total_list)
	mismatch_log.debug('Done determining number of title components in the sequence record IDs.')

	# If sequence IDs are made only of few components, this will enable reasonable alignment in the text file output
	if pipes_most == 1:
		justified = max(justified_list)
	elif pipes_most == 2:
		justified = 45
	elif pipes_most ==3:
		justified = 33
	else:
		justified = 23


	# Writes the title line for csv file and txt file
	Title_components_csv = 'S.No '
	Title_components_txt = '\t\t' + 'S.No'.ljust(7) + '\t'
	for no in range(0, pipes_most):
		Title_components_csv += (',' + 'Title_' + str(no+1))
		Title_components_txt += (('Title_' + str(no + 1)).ljust(justified) + '\t')
	Title_components_csv += ',Sequence Length'
	Title_components_csv += ',No. of Mismatches'
	Title_components_txt += ('Mismatch'.ljust(10)+ '\t' + 'Expected'.ljust(10) + '\n')
	mismatch_log.debug('Determined the title to be used in csv and tect files.')


	if checkboxval_Clustal.get() == True:	# To account for input file used depending on clustal alignment was selected or not
		input_file = SeqQuery_file
	else:
		input_file = Aligned_Filename


	# Extracts residues present at query sites in the reference sequence and also reference sequence in alignment
	mismatch_log.debug('Begins the part to extract residues at the requested query sites.')
	query_site_actual = resi_positions		# query_site_actual refers to actual positioning whereas query_site to python indexes.
	query_site_actual = sorted(query_site_actual)
	query_site = []
	query_site_residues = []
	for n in range(0, len(query_site_actual)):
		query_site.append(int(query_site_actual[n]) - 1)
		query_site_residues.append(Reference.seq[query_site[n]])

	mismatch_log.info('Query sites requested: %s' % (str(query_site_actual)))
	mismatch_log.info('Residues at those query sites: %s' % str(query_site_residues))


	# Finds the position of query site residues in reference's alignment sequence
	query_in_Alignment = []
	for aa in aa_list:
		for query_site_element in query_site:
			if query_site_element in Dict_aa_Index_Reference_Seq[aa]:
				pos = Dict_aa_Index_Reference_Seq[aa].index(query_site_element)
				query_in_Alignment.append(Dict_aa_Index_Reference_Seq_inAlignment[aa][pos])

	query_in_Alignment = sorted(query_in_Alignment)

	#To count for the fact that python indexes from 0 instead of 1
	query_in_Alignment_actual = [(item + 1) for item in query_in_Alignment]

	mismatch_log.info("Corresponding sites in the alignment's Reference sequence : %s" %str(query_in_Alignment_actual))

	# data_alignment = AlignIO.read(Aligned_Filename, "clustal")		# seek function doesn't work here

	# Extracts residues across alignment at the requested query sites
	mismatch_log.info('Begin extracting residues across alignment at the requested query sites')

	dict_clustal_residues = {k: [] for k in aligned_ids_list}
	for record in data_alignment:
		for no in range(0, len(query_in_Alignment)):
			if record.seq[query_in_Alignment[no]] == query_site_residues[no]:
				dict_clustal_residues[record.id].append('')
			else:
				dict_clustal_residues[record.id].append(record.seq[query_in_Alignment[no]])


	# Creates and writes  mismatch details in to a csv file
	mismatch_log.debug('Creates and writes  mismatch details in to a csv file')

	Output_filename_csv = 'Mismatches_Tabulated.csv'
	Output_file_csv = Output_Path + Output_filename_csv

	mismatch_log.info('Mismatch details will be written in to file:  "\n\t%s"' %Output_file_csv )
	try:
		Output_handle_mismatch_csv = open(Output_file_csv, 'w')
	except IOError as exception:
		mismatch_log.error('Error that resulted when opening csv file for writing: %s' % exception)
		if exception.errno != errno.EEXIST:
			temp = ("A file titled '%s' seems to be open in MS Excel. Close that file and try again!" %Output_filename_csv)
			mismatch_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')

	Output_handle_mismatch_csv.write('Data below is obtained from file:    "' + input_file + '"\n')
	Output_handle_mismatch_csv.write('Reference file used:    "' + Reference_file + '"\n')
	Output_handle_mismatch_csv.write('Number of sequences in alignment:    ' +  str(len(data_alignment)) + '\n\n')

	for no in range(0, len(query_site_actual)):
		Title_components_csv += (',' + str(query_site_actual[no]) + ' "' + query_site_residues[no] + '"')
	Title_components_csv += '\n'
	mismatch_log.debug('Created title line for using in csv file')

	Output_handle_mismatch_csv.write("*** Records that have mismatches in at least one of the query sites ***\n\n")
	Output_handle_mismatch_csv.write(Title_components_csv)

	# Gethers residues at query sites and decides if they are mismatch or not.
	mismatch_log.debug('Starting to extract and write mismatch details into csv file')

	sno_1 = 0
	sno_2 = 0
	No_Mismatches_details = ''
	mismatches_none = True			# This var is to get info if there is no mismtaches noted at any of the sites
	for identifier in aligned_ids_list:
		mismatch_count = 0					# to count number of mismatches
		for char in dict_clustal_residues[identifier]:
			if char != '':
				mismatch_count += 1

		pipe_split = identifier.split(id_delimiter)
		if pipe_split[-1] == '':			# Accounts for if pipe is in the end of the id name
			pipe_split = pipe_split[:-1]

		if len(pipe_split) < pipes_most:	# makes sure all IDs have same number of naming components
			for no in range(0, pipes_most - len(pipe_split)):
				pipe_split.append('-na-')

		id_components = ''					# adds name components
		for no in range(0, pipes_most):
			if ',' in pipe_split[no] or '"' in pipe_split[no]:	# standardizes as per csv file format since spreadsheet software will misbehave when a cell value has comma or double quote in it
				pipe_split[no] = pipe_split[no].replace('"','""')
				pipe_split[no] = ('"%s"' % pipe_split[no])

			id_components += (',' + pipe_split[no])

		residues_data = ''					# obtains residue in that position, mismatched or not.
		for no in range(0, len(query_in_Alignment)):
			res = dict_clustal_residues[identifier][no]
			if dict_clustal_residues[identifier][no] == '':		# If residue matches, equal sign is used to improved visibility
				res = '='
			if dict_clustal_residues[identifier][no] == '-':	# Star symbol is used instead of gap symbols to improved visibility
				res = '*'
			residues_data += (',' + res)
		residues_data += '\n'

		len_seq = 'problem_in_script'			# this obtains the length of sequence
		for no in range(0, len(data_alignment)):
				if data_alignment._records[no].id == identifier:
					seq = str(data_alignment._records[no].seq)
					seq = seq.replace('-', '')
					len_seq = len(seq)
					break

		if len_seq == 'problem_in_script':
			mismatch_log.error("Error obtaining sequence length of '%s'." % identifier)

		if mismatch_count > 0:				# Writes into output file if at least one mismatch found in  query positions
			serial_number_1= str (sno_1 + 1)
			Output_handle_mismatch_csv.write(str (serial_number_1) +id_components + ',' + str(len_seq) + ',' +
									str (mismatch_count) + residues_data)
			sno_1 += 1
			mismatches_none = False

		else:								# Gathers details of records without any mismatches in query positions
			serial_number_2= str (sno_2 + 1)
			No_Mismatches_details += (str(serial_number_2) + id_components  + ',' + str(len_seq)  + ',' +
										str(mismatch_count) + residues_data)
			sno_2 += 1

	if mismatches_none:
		Output_handle_mismatch_csv.write(",*** No mismatches at any of the sites requested ***\n")

	# to find unique residue in query sites
	mismatch_log.debug('Unique residues at each requested site will be extracted')

	Dict_Unique_Residues = {k: [] for k in query_site_actual}
	for no in range(0, len(query_site_actual)):
		for res_list in dict_clustal_residues.values():
			Dict_Unique_Residues[query_site_actual[no]].append(res_list[no])

		# this part will enable calculating each unique residue's count and fraction w/o including reference seq's residue.
		# This is so as to offer true calculation as reference seq' residue is always matching anyway.
		if '' in Dict_Unique_Residues[query_site_actual[no]]:
			Dict_Unique_Residues[query_site_actual[no]].remove('')

	Output_handle_mismatch_csv.write("\n\n*** Unique residues seen at the query sites and their count. ***\n"
									 "** (Note: Calculation doesn't include Reference sequences's residue at that position) **\n\n")
	if protein_mode:
		Output_handle_mismatch_csv.write(",Expected Residue, ,Identity_count,% Identity,Similarity_count,% Similarity,"
									 "Unique residues' count and fraction\n")
	else:
		Output_handle_mismatch_csv.write(",Expected Residue, ,Identity_count,% Identity,"
									 "Unique residues' count and fraction\n")

	# aa_set = ['AVFPMILW', 'DE', 'RK', 'STYHCNGQ', 'avfpmilw', 'de', 'rk', 'styhcngq']
	unique_residues_line = ''
	unique_residues_list = []
	identity_count_list = []
	identity_perc_list = []
	similarity_count_list = []
	similarity_perc_list = []
	for no in range(0, len(query_site_actual)):
		for item in range(0, len(Dict_Unique_Residues[query_site_actual[no]])):
			if Dict_Unique_Residues[query_site_actual[no]][item] == '':
				Dict_Unique_Residues[query_site_actual[no]][item] = query_site_residues[no]

		unique_res_detail = ','
		unique_res_detail += (str(query_site_actual[no]) + ' "' + query_site_residues[no] + '",')

		# this writes unique residues for each site in descending order of their count
		Dict_sorting_temp = {}
		for elem in sorted(set(Dict_Unique_Residues[query_site_actual[no]])):
			Dict_sorting_temp[elem] = Dict_Unique_Residues[query_site_actual[no]].count(elem)

		# this obtains identity count and percentage
		if query_site_residues[no] in Dict_sorting_temp:
			identity_count = Dict_sorting_temp[query_site_residues[no]]
			identity_count_list.append(identity_count)
			identity_perc = float(identity_count) / (len(aligned_ids_list) - 1) * 100
			identity_perc_list.append(identity_perc)
		else:
			identity_count = 0
			identity_count_list.append(identity_count)
			identity_perc = 0
			identity_perc_list.append(identity_perc)

		# this obtains similarity count and percentage
		for group in aa_set:		# aa_set is defined in settings text file and allowed to be changed via GUI
			if query_site_residues[no] in group:
				aa_similar = group
				break

		# print Dict_sorting_temp

		similarity_count = 0
		for key in Dict_sorting_temp:
			if key in aa_similar:
				similarity_count += Dict_sorting_temp[key]
		similarity_count_list.append(similarity_count)
		similarity_perc = float(similarity_count) / (len(aligned_ids_list) - 1) * 100
		similarity_perc_list.append(similarity_perc)

		unique_each = ''
		for key, value in sorted(Dict_sorting_temp.iteritems(), key = lambda (k,v): (v,k), reverse = True):
			perc = float (value) / (len(aligned_ids_list) - 1) * 100
			perc = '%3.1f' %perc
			unique_each += (',' + key + ': ' + str(value) + ' (' + str(perc) + '%)')

		# unique_residues_list.append(unique_each)
		if protein_mode:			# residue similarity details only if proteins sequences are used
			unique_res_detail += (',%i,%3.1f,%i,%3.1f' % (identity_count, identity_perc, similarity_count, similarity_perc))
		else:
			unique_res_detail += (',%i,%3.1f' % (identity_count, identity_perc))
		unique_res_detail += unique_each
		unique_residues_line += ( unique_res_detail + '\n')
		unique_residues_list.append(unique_each)		# for text file

	mismatch_log.debug('Done extracting unique residues')

	Output_handle_mismatch_csv.write(unique_residues_line)
	Output_handle_mismatch_csv.write("\n\n*** Records that Do Not have mismatches at any of the query sites ***")
	Output_handle_mismatch_csv.write('\n\n' + Title_components_csv + No_Mismatches_details)
	Output_handle_mismatch_csv.close()
	mismatch_log.info('Done writing mismatch details in to csv file mentioned above.')

	# This part gathers mismatch details to write into a text file.
	mismatch_log.debug('Begin process to write mismatch details in to a "text file"')

	Output_file_Mismatch_txt = Output_Path + 'Mismatches_Detailed.txt'
	Output_handle_mismatch_txt = open(Output_file_Mismatch_txt, 'w')
	Output_handle_mismatch_txt.write('Data below is obtained from file'.ljust(38) + ':  "' + input_file + '"\n' +
				'Reference file used'.ljust(38) + ':  "' + Reference_file + '"\n\n' +
				'Query site(s) for which mismatching are requested'.ljust(55) + ':\t' + str(query_site_actual) + '\n' +
				'Corresponding Query site residue(s) in Reference seq'.ljust(55) + ':\t' +  str(query_site_residues) + '\n' +
				'Their corresponding position(s) in alignment:'.ljust(55) + ':\t' +  str(query_in_Alignment_actual) + '\n\n' +
				"Number of sequences in the provided alignment:".ljust(55) + ':\t' +  str(len(data_alignment)) + '\n\n\n')

	mismatch_log.info('Mismatch details will be written in to text file: \n\t%s' %Output_file_Mismatch_txt)

	# If residues are not matching, this fetches the details about corresponding sequence records
	# for writing in to text file
	mismatch_log.debug('Begins extracting mismatch details to be written in to text file')

	siteno = 0
	for num in range(0, len(query_in_Alignment)):
		aligned_residues = data_alignment[:, query_in_Alignment[num]]		# Obtains residues across alignment in requested position

		element_index = 0
		Mismatch = ''
		increment = 1
		for element in aligned_residues:
			if element != query_site_residues[siteno]:		# Obtains details if residue is a mismatch
				ID_target = data_alignment[element_index].id
				ID_target = str(ID_target)

				pipe_split = ID_target.split(id_delimiter)
				if pipe_split[-1] == '':			# Accounts for if pipe is in the end of the id name
					pipe_split = pipe_split[:-1]

				if len(pipe_split) < pipes_most:	# makes sure all IDs have same number of naming components
					for no in range(0, pipes_most - len(pipe_split)):
						pipe_split.append('-na-')

				id_components = ''					# adds name components
				for no in range(0, pipes_most):
					id_components += (pipe_split[no].ljust(justified) + '\t')

				Details = ('\t\t' + str(increment).ljust(5) + '\t' +
						   id_components + element.ljust(10) + '\t' + query_site_residues[siteno].ljust(10) + '\n')

				Mismatch += Details
				increment += 1

			element_index += 1

		unique_res = unique_residues_list[num].split(',')
		unique_res = unique_res[1:]

		Output_handle_mismatch_txt.write('Site#	' + str(siteno + 1) + '/' + str(len(query_site)) + '\n')
		Output_handle_mismatch_txt.write('Expected Residue'.ljust(45) + ':\t"' +  query_site_residues[siteno]+ '"\n' +
					'Position of this residue in reference seq'.ljust(45) + ':\t' + str(query_site_actual[siteno]) + '\n' +
					'Position of this residue in alignment'.ljust(45) + ':\t' +  str(query_in_Alignment[num] + 1) + '\n' +
					'Residues across alignment at this site'.ljust(45) + ':\t' +  aligned_residues + '\n' +
					'Unique residues present and their count'.ljust(45) + ':\t' +  str(unique_res) + '\n' +
					'% Identity'.ljust(45) + ':\t' +  str(identity_perc_list[num]) +
										 '%%\t(%i residues are identical)' % identity_count_list[num] + '\n')

		if protein_mode:		# residue similarity details only if proteins sequences are used
			Output_handle_mismatch_txt.write('% Similarity'.ljust(45) + ':\t' +  str(similarity_perc_list[num]) +
										 '%%\t(%i residues are similar)' % similarity_count_list[num] + '\n')


		if len(Mismatch) != 0:					# This is to write into output file if mismatches were found
			Output_handle_mismatch_txt.write('\n*** Following are the mismatches noted ***\n\n')
			Output_handle_mismatch_txt.write(Title_components_txt + Mismatch + '\n')
		Output_handle_mismatch_txt.write('\n')
		siteno += 1

		if len(Mismatch) == 0:					# This is if there are no mismatches for a residue
			Output_handle_mismatch_txt.write('*** No mismatches at this site ***\n\n\n')

	Output_handle_mismatch_txt.close()
	mismatch_log.info('Done extracting and writing mismatch details in to text file mentioned above')


# Function that reads clustal alignment( w/ or w/o numbers) and highlights matches and mismatches in color
# at requested sites
def html_formatting():
	html_log.debug("'html_formatting' Module initiated.")
	global query_in_Alignment
	global Output_Path
	# global unique_residues_line
	global identity_perc_list
	global similarity_perc_list
	global query_site_residues
	global query_site_actual
	global Reference
	global protein_mode

	# reads alignment file provided
	html_log.info('Alignment file used for formatting: %s' %Aligned_Filename)
	clustal_aligned = False
	format_supported = False
	try:
		# this block identifies if alignment is in any of the supported formats
		format_type_list = ['clustal','emboss','fasta','fasta-m10','ig','maf','nexus','phylip','phylip-sequential','phylip-relaxed','stockholm']
		for format_type in format_type_list:
			try:
				alignment = AlignIO.read(Aligned_Filename, format_type)
				if format_type == "clustal":
					clustal_aligned = True
				html_log.info('Alignment format is detected as "%s" format' % format_type)
				format_supported = True
				break
			except ValueError:
				pass
	except Exception as e:
		html_log.error("Error that resulted: %s" % e)
		temp = ("Could not open file: '\n\t%s'" % Aligned_Filename)
		html_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	if not format_supported:
		# temp = "Your Alignment file does not have alignment format supported by ResCon"
		temp = 'Format of alignment in your alignment file is not supported by ResCon.'
		html_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	# Extracts alignment data into a dictionary
	id_list = []
	dict_ids = {}
	dict_ids_unedited = {}
	length_id = []
	for no in range(0, len(alignment)):
		temp_id = alignment[no].id
		id_list.append(temp_id)
		dict_ids[temp_id] = list(alignment[no].seq)
		dict_ids_unedited[temp_id] = list(alignment[no].seq)
		length_id.append(len(temp_id))

	if Reference.id not in id_list:
		temp = 'None of the IDs in alignment corresponds to ID of your reference sequence. Fix it and try again!'
		html_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	# determines number of blocks to write output data in clustal format
	max_ID_length = max(length_id)
	len_seq = len(alignment[0].seq)
	if len_seq % 60:
		no_of_blocks = len_seq/60
		no_of_blocks += 1
	else:
		no_of_blocks = len_seq/60

	tag_line_full = [' '] * len(alignment[0].seq)	# enables to have hyperlinked residue positions

	# this part highlights residues depending on they are matching, similar or neither
	# colors:	green - '#42DB33', yellow - '#FFFF00', pink - '#F73D94', orange - '#FFA500', violet_red - '#CC3399'
	# violet - '#9900FF', magenta - '#FF00FF', blue - '#6E9EFF'
	# aa_set = ['AVFPMILW', 'DE', 'RK', 'STYHCNGQ', 'avfpmilw', 'de', 'rk', 'styhcngq']
	# match_color = '#42DB33'
	# similar_color = '#FFFF00'
	# mismatch_color = '#F73D94'
	for key in dict_ids:
		for query_no in range(0, len(query_in_Alignment)):
			residue_ctrl = query_site_residues[query_no]
			similar_aa = ''
			for group in aa_set:		# finds similar amino acid group for user-requested residues
				if residue_ctrl in group:
					similar_aa = group
					break

			residue = dict_ids[key][query_in_Alignment[query_no]]
			if residue == query_site_residues[query_no]:
				dict_ids[key][query_in_Alignment[query_no]] = ('<strong><span style="background-color: %s">%s</span></strong>' % (match_color, residue))
			elif protein_mode and residue in similar_aa:		# similar residues colored only for protein sequences mode
				dict_ids[key][query_in_Alignment[query_no]] = ('<strong><span style="background-color: %s">%s</span></strong>' % (similar_color, residue))
			else:
				dict_ids[key][query_in_Alignment[query_no]] = ('<strong><span style="background-color: %s">%s</span></strong>' % (mismatch_color, residue))

			# for tag line (numbers will be written in vertical format so as to position them directly above the column)
			tag_no_str = str(query_site_actual[query_no])
			tag_no_ver = '<vert>'
			for char_no in range(0, len(tag_no_str)):
				if char_no != len(tag_no_str):
					tag_no_ver += (tag_no_str[char_no] + '<br />')
					# ''.join([tag_no_ver, (tag_no_str[char_no] + '<br />')])
				else:
					tag_no_ver += (tag_no_str[char_no] + '</vert>')
					# ''.join([tag_no_ver, (tag_no_str[char_no] + '</vert>')])

			tag_line_full[query_in_Alignment[query_no]] = '<font color= "red"><a name="%i">%s</a></font>' %(query_site_actual[query_no], tag_no_ver)

	# opens a new html file that will have color formatting
	html_out_name = Output_Path + "Formatted_Alignment.html"
	html_log.info("Formatted aligment will be stored in html file: '%s'" % html_out_name)
	out_html_handle = open(html_out_name, 'w')

	if chart_method == 'chart.js':
		html_template_file = "Template_html5_chartjs.txt"
	elif chart_method == 'chartnew.js':
		html_template_file = "Template_html5_ChartNewjs.txt"

	if chart_method == 'chart.js' or chart_method == 'chartnew.js':
		with open(html_template_file, 'Ur') as html_template_data:
			html_remaining = ''
			check_line = 0
			for line in html_template_data:
				if check_line == 0:
					out_html_handle.write(line)
				else:
					html_remaining += line
				if line.replace('\n', '') == 'check_by_ResCon-->':
					check_line = 1

	if chart_method == 'pro':
		# this part draws bar chart using matplotlib
		plot_filename = 'plot.png'
		plot_filepath = Output_Path + plot_filename

		iden = identity_perc_list
		simi = similarity_perc_list
		width = 0.25
		numb = range(len(query_site_residues))
		numb2 = [x+width for x in numb]

		fig = plt.figure()
		graph = fig.add_subplot(1,1,1)
		plt.bar(numb, iden, width = width, color='blue', label= '% Identity')
		if protein_mode:
			plt.bar(numb2, simi, width = width, color='yellow', label = '% Similarity')
		plt.xticks(numb2, query_site_actual)		# label for bars
		plt.gca().yaxis.grid(True)						# draws horizontal grid lines
		graph.set_yticks([10,30,50,70, 90], minor=True)	# drwas minor horizontal grid lines
		graph.grid(which='minor', alpha=0.5)

		plt.autoscale()		# autoscales for x-axis but y-axis's is manually set later
		plt.ylim([0, 100])
		plt.xlim([numb[0]-0.3, numb[-1]+1])	# if x limits are not set, mac os removes first and last bins if they have zero value
		plt.ylabel('%')
		plt.xlabel('Residue Position')
		if len(query_site_actual) > 5:
			plt.tick_params(labelright=True)	# tick labels on right side as well
		graph.legend(loc='upper center', bbox_to_anchor=(0.5, 1.14),fancybox=True, shadow=True, ncol=2)	# positions legened box on top of chart

		fig = plt.gcf()		# for manipulating image size
		fig.set_size_inches(len(query_site_actual), 5)
		fig.savefig(plot_filepath, bbox_inches='tight')

		html_log.info('Created bar chart using matplotlib successfully!')

		# creates html header and instructs browser to show text using font Lucida Sans Typewriter
		out_html_handle.write(
			'<HTML>\n<HEAD>\n\t<META HTTP-EQUIV="CONTENT-TYPE" CONTENT="text/html; charset=utf-8">\n\t<TITLE>ResCon Formatted</TITLE>\n'
			'\t<STYLE>\n'
			'\t\t* {\n'
			'\t\t\tfont-family:"Lucida Sans Typewriter", "Lucida Console", Monaco, "Bitstream Vera Sans Mono", monospace;\n'
			'\t\t\tfont-size: 13px;\n'
			'\t\t}\n\t\tvert {\n'
			'\t\t\tdisplay: inline-block;\n'
			'\t\t\tvertical-align: middle;\n\t\t}\n'
			'\t</STYLE>\n</HEAD>\n'
			'<BODY LANG="en-US" DIR="LTR">\n<PRE CLASS="western">\n\n')

		# Link bar chart into html file
		out_html_handle.write('    <IMG SRC="%s" ALT="some text", height= 350>\n\n' %plot_filename)


	elif chart_method == 'chart.js':
		out_html_handle.write('<div>\n  <canvas id="canvas" height="350", width = "%i" ></canvas>\n'
							  '<div id="legend"></div>\n</div>\n' %(63*len(query_site_actual)))

		out_html_handle.write('<script>\n\tvar barChartData = {\n'
							  '\t\tlabels : %s,\n'
							  '\t\tdatasets : [\n'
							  '\t\t{	fillColor : "rgba(51,51,255,0.7)",\n'
							  '\t\t\thighlightFill: "rgba(0,0,255,1)",\n'
							  '\t\t\tlabel: "%% Identity",\n'
							  '\t\t\tdata : %s  },\n\n'
							  '\t\t{	fillColor : "rgba(255,140,0, 0.7)",\n'
							  '\t\t\thighlightFill : "rgba(255,140,0, 1)",\n'
							  '\t\t\tlabel: "%% Similarity",\n'
							  '\t\t\tdata : %s  }\n\t\t]\n\t}\n'
							  '\twindow.onload = function(){\n'
							  '\t\tvar bar = new Chart(document.getElementById("canvas").getContext("2d")).Bar(barChartData,\n'
							  '\t\t{\n\t\t\ttooltipTemplate: "<%%if (label){%%><%%=label%%>: <%%}%%><%%= value %%>kb",\n'
							  '\t\t\tanimation: false,\n'
							  '\t\t\tresponsive : false,\n'			# if true, chart gets messed up when browser is resized
							  '\t\t\tbarValueSpacing : 13,\n'
							  '\t\t\tscaleShowVerticalLines: false,\n'
							  '\t\t\tbarShowStroke: false,\n'				# this avoids having misleading color in case the value of a bar is zero
							  '\t\t\tmultiTooltipTemplate: "<%%= datasetLabel %%> - <%%= value %%>",\n'
							  '\t\t});\n\n\tvar legendHolder = document.createElement("div");\n'
							  '\tlegendHolder.innerHTML = bar.generateLegend();\n'
							  '\tdocument.getElementById("legend").appendChild(legendHolder.firstChild);'
							  '\n\t}\n</script>\n' % (query_site_actual, identity_perc_list, similarity_perc_list))

	elif chart_method == 'chartnew.js':
		if ( 63*len(query_site_actual) ) < 325:
			width_canvas = 325
		else:
			width_canvas = 63*len(query_site_actual)

		out_html_handle.write('</HEAD>\n<BODY LANG="en-US" DIR="LTR">\n<PRE CLASS="western">\n\n'
							  '<div>\n  <canvas id="canvas_bar" height="350", width = "%i" ></canvas>\n'
							  '<div id="legend"></div>\n</div>\n\n<script>\n'
							  '\tvar barChartData = \n\t{\n'
							  '\t\tlabels : %s,\n\t\tdatasets : [\n'
							  '\t\t{	fillColor : "rgba(51,51,255,0.8)",\n\t\t\ttitle: "%% Identity",\n'
							  '\t\t\tdata : %s  },\n\n'
							  '\t\t{	fillColor : "rgba(255,140,0, 1.0)",\n\t\t\ttitle: "%% Similarity",\n'
							  '\t\t\tdata : %s  }\n\t\t]\n\t}\n'
							  % (width_canvas,query_site_actual, identity_perc_list, similarity_perc_list))

		out_html_handle.write(html_remaining)


	# Creates a table at top of html page that will show similarity and identity % at requested sites.
	# For that, we borrow query sites, identity and similiarity at those sites from fetch_mismatch module
	row1 = []
	for n in range(0, len(query_site_actual)):
		temp = '<a href="#%i">%i<BR>%s</a>' % (query_site_actual[n], query_site_actual[n], query_site_residues[n])	# hyperlink included
		row1.append(temp)
	row2 = []
	row3 = []
	for n in range(0, len(identity_perc_list)):
		row2.append("%3.1f" % identity_perc_list[n])
		if protein_mode:
			row3.append("%3.1f" % similarity_perc_list[n])

	if protein_mode:		# for protein seqs, similarity data is added
		rows = [row1, row2, row3]
		heading = ['Position', '% Identity', '% Similarity']
	else:
		rows = [row1, row2]
		heading = ['Position', '% Identity']
	limit_col = 20				# this variable determines how many columns are allowed in each table
	no_of_tables = (len(query_site_residues) - 1 + limit_col) / limit_col

	limit_col_temp = limit_col
	for table_no in range(0, no_of_tables):
		table_string = ('<TABLE border="1" cellpadding="0" cellspacing="0" HEIGHT="90px" style="table-layout:fixed">\n')
		if table_no < (no_of_tables - 1):
			end_col = limit_col_temp
		else:
			end_col = len(query_site_residues)

		for no in range(0, len(rows)):
			if not no:
				table_string += ('\t<TR ALIGN="CENTER" BGCOLOR="#E3E3E0">\n')
			else:
				table_string += ('\t<TR ALIGN="CENTER">\n')

			table_string += ('\t\t<TH scope="row" width=130>%s</TH>\n' % heading[no])
			for ele in range((limit_col_temp -limit_col), end_col):
				table_string += ('\t\t<TD width=60>%s</TD>\n' % rows[no][ele])

			table_string += '\t</TR>\n'
		table_string += ('</TABLE>\n')
		out_html_handle.write(table_string)
		limit_col_temp += limit_col


	# this part creates legend for highlighting colors used
	temp = '-' * 90 + '\n\n'
	temp += ' ' * 30 + '<strong>Legend: Colors used for highlighting</strong>\n\n'
	temp += '  <span style="background-color: %s">        </span> Matching residue' % match_color
	if protein_mode:
		temp += (' '* 22 + '<span style="background-color: %s">        </span> Mismatch but a <ins>similar</ins> residue \n\n' % similar_color)
		temp += ('  <span style="background-color: %s">        </span> Mismatch but a <ins>non-similar</ins> residue' % mismatch_color)
	else:
		temp += ('    <span style="background-color: %s">        </span> Mismatching residue' % mismatch_color)

	temp += "    <span style='background-color: #00FFFF'>        </span> Reference Sequence's ID"
	out_html_handle.write(temp + '\n\n' + '-' * 90 + '\n\n')


	# writes formatted data in clustal alignment format (manually done) into a output file
	for block_no in range(1, no_of_blocks+1):
		end_no = block_no * 60
		start_no = end_no - 60

		tag_line = ''.join(tag_line_full[start_no:end_no])
		tag_line = ''.ljust(max_ID_length+6) + tag_line
		out_html_handle.write(tag_line + '\n')

		for identifier in id_list:
			new_line = identifier.ljust(max_ID_length+6)
			if identifier == Reference.id:
				new_line = '<strong><span style="background-color: #00FFFF">' + new_line[0:max_ID_length+1] + '</span></strong>'
				new_line += ''.ljust(5)
			new_line += (''.join(dict_ids[identifier][start_no:end_no]))
			# ''.join([new_line, (''.join(dict_ids[identifier][start_no:end_no]))])

			seq_end_pos = (''.join(dict_ids_unedited[identifier][0:end_no]))
			seq_end_pos = seq_end_pos.replace('-', '')
			seq_end_pos = len(seq_end_pos)

			new_line += ('   ' + str(seq_end_pos))
			# ''.join([new_line, ('   ' + str(seq_end_pos))])
			out_html_handle.write(new_line + '\n')

		if clustal_aligned:			# to write consensus symbol data available in clustal format
			symbol_line = alignment._star_info[start_no:end_no]
			symbol_line = ''.ljust(max_ID_length+6) + symbol_line
			out_html_handle.write(symbol_line + '\n\n')
		else:
			out_html_handle.write('\n')

	temp = '\n<ins>Input files used:</ins>\n\t %s \n\t %s' %(Aligned_Filename, Reference_file)
	out_html_handle.write(temp + "\n</PRE>\n</BODY>\n</HTML>")
	out_html_handle.close()
	html_log.info('Alignment formatted was completed and saved as above mentioned html file.')


# Function that reads blast xml file and then writes fasta files into an output file - all based on
# e-value entered by user
def blast_filter():
	blast_log.info("User selected tool - 'Filtering sequences by Blast e-value'.")
	global Blast_xml, Blast_fasta, E_threshold_entry
	# global Higher_E_val_checkboxval, Lower_E_val_checkboxval
	global Blast_Submit, output_blast_entry, blast_processing
	global lower_or_higher_e_threshold
	global top
	Blast_Submit.configure(state=DISABLED)
	blast_processing.grid()

	# Reads E-value threshold provided by user
	try:
		E_threshold = float(E_threshold_entry.get())
	except ValueError as exception:
		blast_log.error("Error that resulted: %s" % exception)
		temp = 'E-threshold box is empty or number entered is invalid. Enter a valid number'
		blast_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		temp = ('E-threshold entered by user: %s' % str(E_threshold))
		blast_log.info(temp)

	# Determines if user has requested to filter sequences in blast with e-values lesser or higher than threshold
	if lower_or_higher_e_threshold.get() == 1:
		Requested_Lower = True
		Requested_Higher = False
		blast_log.info("User chose 'Lower than E-threshold'")
	else:
		Requested_Lower = False
		Requested_Higher = True
		blast_log.info("User chose 'Higher than and equal to E-threshold'")

		temp = "Note that checking 'Higher than and equal to E-threshold' will also record sequences that are not in blast " \
			   "XML file but are present in provided fasta file that do not have any e-value listed in XML file"
		blast_log.info(temp)
		tkMessageBox.showinfo('Note', temp, parent = top)

	# Reads Blast xml file provided by user
	blast_handle = Blast_xml.get()
	# blast_handle = "C:/Users/Mana/Desktop/ResCons/Anaibogi/seqs/blast/Blast_Against_Family_II_IPPases.xml"
	if not blast_handle:
		temp = "Please enter file path of Blast XML file."
		blast_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		temp = ("Blast XML file provided by user :\n\t%s" %blast_handle)
		blast_log.info(temp)

	# Reads user-provided blast xml file
	try:
		blast_records = NCBIXML.read(open(blast_handle))
	except Exception as e:
		temp = ('Having trouble reading file: \n\t%s. \nMake sure it is in valid xml format' % blast_handle)
		blast_log.error(temp)
		blast_log.error('Error that resulted: %s' % e)
		popup_error(temp)
		raise_enabler('stop')
	else:
		temp = ('Can successfully read Blast xml file provided by user: \n\t%s' % blast_handle)
		blast_log.debug(temp)

	# Reads fasta file provided by user
	if not Blast_fasta.get():
		temp = "Please enter file path of fasta file."
		blast_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		Input_handle_blast_fasta = Blast_fasta.get()
		temp = ("Fasta file provided by user :\n\t%s" %Input_handle_blast_fasta)
		blast_log.info(temp)

	# Reads output folder path provided by user
	if not output_blast_entry.get():
		temp = "Output folder is empty. Pick a folder."
		blast_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		Output_Path_blast = output_blast_entry.get()
		if Output_Path_blast[-1] != '/':
			Output_Path_blast += '/'
		verify_filepath_exists(Output_Path_blast)
		temp = ("Output folder provided by user :\n\t%s" %Output_Path_blast)
		blast_log.info(temp)


	# this part collects IDs of sequences, that have e-values less than requested e-score threshold, into a list
	Lower_than_E_threshold_list = []
	lower = 0
	for index in range(0, len(blast_records.descriptions)):
		E_score = float(blast_records.descriptions[index].e)
		title = str(blast_records.descriptions[index].title)
		title_list = title.split()
		title_id = title_list[1]  # Description id in the seqs used for blast must be unique for each seq.
		if E_score < E_threshold:
			Lower_than_E_threshold_list.append(title_id)
			lower += 1

	Filtered_records = []
	total = 0
	filtered_lower = 0
	filtered_higher = 0
	for seq_record in SeqIO.parse(Input_handle_blast_fasta, "fasta"):
		seq_id = seq_record.id

		# Appends IDs of records that has e-score less than requested threshold
		if Requested_Lower and seq_id in Lower_than_E_threshold_list:
			Filtered_records.append(seq_record)
			filtered_lower += 1

		# Appends IDs of records that has e-score highr than requested threshold
		# This will also record IDs of sequences that are not in blast xml file but are present in fasta file
		if Requested_Higher and seq_id not in Lower_than_E_threshold_list:
			Filtered_records.append(seq_record)
			filtered_higher += 1

		total += 1

	# Output_Path_blast = os.path.dirname(blast_handle)
	# Output_Path_blast += '/Output/'
	# verify_filepath_exists(Output_Path_blast)

	Output_handle_blast = Output_Path_blast + 'Filtered_Sequences_E_' + str(E_threshold) + '.fasta'
	SeqIO.write(Filtered_records, Output_handle_blast, 'fasta')

	blast_log.info('Total number of sequences in input fasta file:	 %s' % str(total))
	blast_log.info("Filtered sequences are stored in file:  \n\t'%s'" % Output_handle_blast)

	# if Requested_Lower and seq_id in Lower_than_E_threshold_list:
	if Requested_Lower:
		temp = ('Total number of sequences in input fasta file:	 %s' % str(total)) \
				+ ('\nTotal number of records with E-value less than entered threshold in blast xml file:  %s' % str(lower)) \
				+ ('\nNumber of sequences written in above output file:  %s' % str(filtered_lower))
		blast_log.info(temp)
		tkMessageBox.showinfo("Job done!", temp, parent = top)

	# if Requested_Higher and seq_id not in Lower_than_E_threshold_list:
	if Requested_Higher:
		temp = ('Total number of sequences in input fasta file:  %s' % str(total)) \
			   	+ ('\nTotal number of records with E-value less than entered threshold in blast xml file:  %s' % str(lower)) \
				+ ('\nNumber of sequences written in above output file:  %s' % str(filtered_higher))
		blast_log.info(temp)
		tkMessageBox.showinfo("Job Done!", temp, parent = top)

	Blast_Submit.configure(state=ACTIVE)
	blast_processing.grid_remove()
	blast_log.info("Blast-Filtering job is done! Ready for next job\n\n")
	blast_log.debug("Exiting 'blast_filter' module.\n")

	log_path = Output_Path_blast + 'log.txt'
	shutil.copyfile(logger_filename, log_path)


# Function that reads GenPept/GenBank files and then outputs those sequences in fasta format - with their
# header having info about their taxanomy
def genpept_converter():
	genbank_log.info("User selected tool - 'GenPept/GenBank to Fasta converter'.")
	global genbank_entry, fasta_id_length_entry, ignore_incomplete_checkboxval, options_dict, genbank_Submit, top, output_genbank_entry, genbank_processing

	genbank_Submit.configure(state=DISABLED)
	genbank_processing.grid()

	# read fasta id length provided
	if fasta_id_length_entry.get():
		fasta_id_length = fasta_id_length_entry.get()
	else:
		temp = 'Fasta ID length box cannot be empty.'
		genbank_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	fasta_id_len_test = re.findall(r'[^0-9 ]', fasta_id_length)
	if fasta_id_len_test:
		temp = "Invalid character(s) found in 'Fasta id length' entered!"
		genbank_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		fasta_id_length = int(fasta_id_length)

	if genbank_entry.get():
		Input_Path_genbank = genbank_entry.get()
		temp = ("GenPept/GenBank file provided by user: \n %s" % Input_Path_genbank)
		genbank_log.info(temp)
	else:
		temp = "Enter valid GenPept/GenBank file path."
		genbank_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	if output_genbank_entry.get():
		Output_Path_genbank = output_genbank_entry.get()
		if Output_Path_genbank[-1] != '/':
			Output_Path_genbank += '/'
		verify_filepath_exists(Output_Path_genbank)
		temp = ("Output folder path provided by user: \n %s" % Output_Path_genbank)
		genbank_log.info(temp)
	else:
		temp = "Output folder path is empty. Pick a folder."
		genbank_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	# This part will extract taxanomy details from GenPept/GenBank records
	genbank_log.info("State of 'Replace newick-hating-symbols with underscores':   %s" % str(bool(Newick_hate_sym_checkboxval.get())))
	genbank_log.info("State of 'Ignore sequence records that are incomplete/partial?':   %s" % str(bool(ignore_incomplete_checkboxval.get())))
	genbank_log.info("Length of fasta ID requested by user':   %s" % str(fasta_id_length))
	genbank_log.debug('Starting to extract taxonomy details for each sequence record.')

	# options_list is used to extract details requested by user for the header id.
	# If this list is changed, change the list in tkinter section as well. Both must have same list elements.
	options_list = ['Locus Name', 'GI', 'Domain', 'Phylum', 'Class', 'Genus', 'Species', 'Seq Length', 'Seq Name']
	taxanomy_options = {k: '' for k in options_list}
	options_log = {k: bool(options_dict[k].get()) for k in options_list}
	genbank_log.info("State of 'Header options':  %s" % options_log)

	Fasta_Output = []
	WithDomainName = 0
	NoDomainName = 0
	hate_128 = 0
	ignore_incomplete = 0
	# Input_file = "NCBI_RefSeq_Data_IPPase.gp"
	for record in SeqIO.parse(Input_Path_genbank, "genbank"):
		seq = record.seq

		taxanomy_options['GI'] = record.annotations['gi']		# gets gi info
		taxanomy_options['Seq Length'] = str(len(seq))			# gets seq length info

		try:
			m = re.match(r"(.+) \[.+", record.description)
			taxanomy_options['Seq Name'] = (m.group(1)).replace(' ', '_')
		except AttributeError:
			taxanomy_options['Seq Name'] = (record.description).replace(' ', '_')

		organism = record.annotations["organism"]
		organism = organism.split()

		if len(organism) == 1:			# Organsims' binomial name in GenPept/GenBank records may vary from 0 to more than 2 words
			taxanomy_options['Genus'] = organism[0]
			taxanomy_options['Species'] = 'NA'

		elif len(organism) > 1:
			taxanomy_options['Genus'] = organism[0]
			taxanomy_options['Species'] = '_'.join(organism[1:(len(organism))])

		else:
			taxanomy_options['Genus'] = 'NA'
			taxanomy_options['Species'] = 'NA'

		if len(record.annotations["taxonomy"]) != 0:			# obtains domain details
			taxanomy_options['Locus Name'] = record.name
			taxanomy_options['Domain'] = record.annotations["taxonomy"][0]
			taxanomy_options['Domain'] = taxanomy_options['Domain'].replace(' ', '_')

			if len(record.annotations["taxonomy"]) > 1:			# obtains phylum details
				taxanomy_options['Phylum'] = record.annotations["taxonomy"][1]
				taxanomy_options['Phylum'] = taxanomy_options['Phylum'].replace(' ', '_')
			else:
				taxanomy_options['Phylum'] = 'NA'

			if len(record.annotations["taxonomy"]) > 2:			# obtains class details
				taxanomy_options['Class'] = record.annotations["taxonomy"][2]
				taxanomy_options['Class'] = taxanomy_options['Class'].replace(' ', '_')
			else:
				taxanomy_options['Class'] = 'NA'

			WithDomainName += 1

		else:													# if taxonomy details are not available
			taxanomy_options['Locus Name'] = record.name
			taxanomy_options['Domain'] = 'NA'
			taxanomy_options['Phylum'] = 'NA'
			taxanomy_options['Class'] = 'NA'

			NoDomainName += 1

		# puts together sequence and its taxanomy details
		# name_modified_temp = name + '|' + domain + '|' + phylum + '|' + class_info + '|' + genus + '|' + species + '|'

		name_modified_temp = ''
		# name_modified_temp = taxanomy_options['Locus_Name'] + '|'
		for option in options_list:
			if options_dict[option].get():
				name_modified_temp += (taxanomy_options[option] + connector_id)		# connector_id is defined in settings and can be edited from GUI settings window

		# depending on user requirement, this part replaces symbols, that can cause trouble when fasta sequences
		# produced are used in making phylogenetic tree, with underscores.

		if Newick_hate_sym_checkboxval.get():
			name_modified = ''
			for char in name_modified_temp:
				if char in ['(', ')', '[', ']', ':', ';', ',', "'"]:
					# name_modified += '_'
					name_modified += newick_sym_replace		# newick_sym_replace is defined in settings text file and can be edited through GUI settings window
				else:
					name_modified += char
		else:
			name_modified = name_modified_temp


		# this part restricts max length of fasta id to the length user requests.
		# Note: Max length of sequence id in clustal omega seems to be 127 (though I haven't found official record
		# suggesting this condition; 127 was rather found to be the limit from personal experience).
		if len(name_modified) > fasta_id_length:
			name_modified = name_modified[0:fasta_id_length]
			hate_128 += 1

		# This part ignores incomplete seqs, if user has requested.
		if ( ( "COMPLETENESS: INCOMPLETE" in (record.annotations["comment"]).upper() or
					"PARTIAL" in (record.description).upper() ) and
					ignore_incomplete_checkboxval.get()  ):

			ignore_incomplete += 1
		else:
			record_modified = SeqRecord(seq, id=name_modified, description=record.description)
			# record_modified = SeqRecord(seq, id=name_modified, description='')
			Fasta_Output.append(record_modified)

	genbank_log.info('Total number of fasta ID names whose names were longer than %i characters and \n hence restricted '
					 'to the maximum of %i characters: %s' % (fasta_id_length, fasta_id_length, str(hate_128)))
	genbank_log.info('Number of sequence records that are incomplete:  %s' % str(ignore_incomplete))
	genbank_log.debug('Done extracting taxonomy details from GenPept/GenBank sequence records')

	# Create output folder, then output fasta file and write fasta sequences with custom header
	# Output_Path_genbank = os.path.dirname(Input_Path_genbank)
	# Output_Path_genbank += '/Output/'
	# verify_filepath_exists(Output_Path_genbank)

	Output_handle_genbank = Output_Path_genbank + "GenPept_to_Fasta_formatted.fasta"
	SeqIO.write(Fasta_Output, Output_handle_genbank, 'fasta')

	# logging details
	temp = (' Number of records present in input GenPept/GenBank file:  ' +  str(NoDomainName + WithDomainName)) \
		   + '\n Number of sequences written into Fasta file:  '+ str(NoDomainName + WithDomainName - ignore_incomplete)\
		   + ('\n GenPept/GenBank to Fasta conversion was successful and stored in file:\n  ' + Output_handle_genbank)
	genbank_log.info(temp)
	tkMessageBox.showinfo('Job Done', temp, parent = top)

	temp = (' Total number of records that have some level of taxonomy details:   ' +  str(WithDomainName)) \
		   + ('\n Total number of records that do not have any level of taxonomy details:   ' + str(NoDomainName))
	genbank_log.info(temp)

	genbank_Submit.configure(state=ACTIVE)
	genbank_processing.grid_remove()
	genbank_log.info("Exiting 'genpept_converter' module. Ready for next job! \n\n")

	# copying log file
	log_path = Output_Path_genbank + 'log.txt'
	shutil.copyfile(logger_filename, log_path)


# Function that reads phylogenetic tree and then outputs a fasta file with sequences from a
# subtree selected by user. Subtree is specified using clade's branch length.
def Extract_Clades():
	clades_log.info("User requested tool 'Extract Subtree Sequences'")
	global Newick_entry, Fasta_Original_entry, BranchLength_entry, Extract_Submit, output_subtree_entry, top, subtree_processing

	Extract_Submit.configure(state = DISABLED)
	subtree_processing.grid()

	# reads user input - Newick file
	if not Newick_entry.get():
		temp = "'Newick file' box is empty. Choose a file."
		clades_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		Extract_Input_Newick = Newick_entry.get()
		clades_log.info('Newick tree file provided by user:\n\t%s' % Extract_Input_Newick)

	# gets fasta file
	if not Fasta_Original_entry.get():
		temp = "'FASTA file' box is empty. Choose a file."
		clades_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		Extract_input_fasta = Fasta_Original_entry.get()
		clades_log.info('Fasta file provided by user:\n\t%s' % Extract_input_fasta)

	# gets output folder
	if not output_subtree_entry.get():
		temp = "'Output folder' box is empty. Choose a folder."
		clades_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		Output_Path_Extract = output_subtree_entry.get()
		if Output_Path_Extract[-1] != '/':
			Output_Path_Extract += '/'
		verify_filepath_exists(Output_Path_Extract)
		clades_log.info('Output folder path provided by user:\n\t%s' % Output_Path_Extract)

	# gets branch length
	Branch_Length_entry = BranchLength_entry.get()
	if Branch_Length_entry == '':
		temp = 'Branch length box cannot be empty.'
		clades_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	Branch_Length_entry = Branch_Length_entry.replace(' ', '')
	Branch_Length_list = Branch_Length_entry.split(",")
	Branch_Length_list = [k for k in Branch_Length_list if k]

	clades_log.info('Branch length(s) provided by user:  %s' % Branch_Length_entry)
	clades_log.info('Branch length(s) recognized by ResCon: %s' % Branch_Length_list)

	format_supported = False
	try:
		formats_tree_list = ['nexml', 'phyloxml', 'newick', 'nexus']
		for format_type in formats_tree_list:
			try:
				tree = Phylo.read(Extract_Input_Newick, format_type)
				clades_log.info('Phylogenetic Tree format is detected as "%s" format' % format_type)
				format_supported = True
				break
			except Exception as e:
				pass
	except Exception as e:
		clades_log.error("Error that resulted: %s" % e)
		temp = ("Could not open file: '\n\t%s'" % Extract_Input_Newick)
		clades_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	if not format_supported:
		temp = "Phylogenetic Tree file provided does not have tree format supported by ResCon.\n " \
			   "'Supported formats: 'newick', 'nexus', 'nexml', 'phyloxml'"
		clades_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	clades_log.debug("Phylogenetic Tree file '%s' was read successfully." % Extract_Input_Newick)

	# this part collects all branch lengths from newick tree
	clades_log.debug('Making list of all branch lengths present in newick tree')
	Module = tree.get_nonterminals()

	Module_Branch_Length = []
	for item in Module:
		Branch_Length_item = item.branch_length
		Module_Branch_Length.append(Branch_Length_item)


	# processes each branch length provided by user
	for number in range(0, len(Branch_Length_list)):
		Branch_Length_str = Branch_Length_list[number]
		clades_log.info('Subtree with Branch length "%s" that will be processed now' % str(Branch_Length_str))

		#determines if user provided branch length has to be integer or float type
		# self-note: Convert this part to pythonic!
		if 'e' in Branch_Length_str or 'E' in Branch_Length_str  or '.' in Branch_Length_str or '^' in Branch_Length_str:
			Branch_Length = float(Branch_Length_str)
		else:
			Branch_Length = int(Branch_Length_str)

		clades_log.info('This software understands the branch length as "%s" and will search for it in the '
						'newick tree provided.' % str(Branch_Length))


		#This part verifies if requested branch length is present in the newick file
		#and it is not in duplicate.
		clades_log.info('Determining if requested branch length is present in the newick file '
						 'and also verifying it is not in duplicate.')
		index_list = []
		for index in range(0, len(Module_Branch_Length)):
			if str(Branch_Length) == str(Module_Branch_Length[index]):
				index_list.append(index)

		if len(index_list) == 1:
			clades_log.info('Branch length found. Sequences extraction in progress....')

			Clade_Of_Interest = tree.get_nonterminals()[index_list[0]].depths().keys()

			Clade_SeqID = []
			for SeqNo in range(0, len(Clade_Of_Interest)):
				id = str(Clade_Of_Interest[SeqNo].name)
				if id != 'None':
					Clade_SeqID.append(id)

			Fasta_id = []
			Extracted_Seqs = []
			for seq_record in SeqIO.parse(Extract_input_fasta, "fasta"):
				title = seq_record.description
				title_list = title.split()
				id = title_list[0]
				Fasta_id.append(id)
				if id in Clade_SeqID:
					Extracted_Seqs.append(seq_record)


			clades_log.info('Total number of sequences in clade with branch length "%s" :		%i' % (
				str(Branch_Length), len(Extracted_Seqs)))

			Output_file_Extract = Output_Path_Extract + 'Clade_' + str(Branch_Length) + '_Specific_sequences.fasta'
			Output_handle_Extract = open(Output_file_Extract, 'w')
			SeqIO.write(Extracted_Seqs, Output_handle_Extract, 'fasta')

			temp = ("Fasta sequences corresponding to clade %s were extracted and saved in to: \n\t '%s'"
					% (str(Branch_Length), Output_file_Extract))
			clades_log.info(temp)
			temp += ("\n\nTotal number of sequences extracted: %i" % len(Extracted_Seqs))
			tkMessageBox.showinfo( ('Subtree Extracted (%i/%i)' %(number + 1, len(Branch_Length_list)) ), temp, parent = top)

			Output_handle_Extract.close()

		# If branch length requested is not found in tree. It has to be exactly same number for this script to work.
		# Rounding off leads to error. This part also provides suggestions by showing coser numbers to branch length.
		elif len(index_list) < 1:
			branches_sorted = sorted(Module_Branch_Length)

			closest = min(range(len(branches_sorted)), key=lambda x: abs(branches_sorted[x] - Branch_Length))
			probables = branches_sorted[closest-2 : closest+3]
			probables_str = ''
			for no in probables:
				probables_str += ('  ' + str(no) + '\n')

			temp = ("Branch Length '%s' is not found in tree. Perhaps your actual branch length may be one among "
					"the following?\n\n%s\nPick correct number and try again!" % (Branch_Length, probables_str))
			clades_log.error(temp)
			tkMessageBox.showerror('Error processing Branch Length (%i/%i)' %(number + 1, len(Branch_Length_list)),temp, parent = top)
			raise_enabler('stop')

		# if branch length requested is in duplicate. Work-around is to edit tree file (preferably a copy) itself
		# so as to have change one of duplicate branch length to some other number
		else:
			temp = 'More than one instances of such branch length (%s) found. ' \
				   'This tool works when there is only one branch length. \nDirty trick: Edit tree file to ' \
				   'change your branch length to a unique number!' % str(Branch_Length_str)
			clades_log.error(temp)
			tkMessageBox.showerror('Error processing Branch Length (%i/%i)' %(number + 1, len(Branch_Length_list)),temp)
			raise_enabler('stop')

		if len(Branch_Length_list) > 1:
			clades_log.info("Time to process next branch length provided in user's list.\n")

		log_path = Output_Path_Extract + 'log.txt'
		shutil.copyfile(logger_filename, log_path)

	subtree_processing.grid_remove()
	Extract_Submit.configure(state = ACTIVE)
	clades_log.info("Exiting 'Extract_Clades' Module. Ready for next job\n\n")


# Function that reads a fasta file and extracts description or identifier from sequences' header
# (their header has to be in specified format).
def description_extractor():
	global descr_extr_file_entry
	global replace_comma_checkboxval
	global descr_extractor_Submit
	global output_descr_extract_entry
	global descr_extract_processing
	global reg_exp_descr_entry
	global partial_descr_checkboxval
	# global descr_option
	global id_option

	header_extractor_log.info("User has requested tool 'Fasta Description/ID Extractor'")
	descr_extractor_Submit.configure(state=DISABLED)
	descr_extract_processing.grid()

	if id_option:
		header_extractor_log.info("User chose to use 'Identifier Extractor'")
		keyword = 'Identifier'
	else:
		header_extractor_log.info("User chose to use 'Description Extractor'")
		keyword = 'Description'

	# reads file provided by user
	if not descr_extr_file_entry.get():
		temp = "'FASTA file' box is empty. Choose a file."
		header_extractor_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		input_fasta_file = descr_extr_file_entry.get()
		header_extractor_log.info("Fasta file provided by user: '%s'" % input_fasta_file)

	# reads output path provided by user
	if not output_descr_extract_entry.get():
		temp = "'Output folder' box is empty. Choose a folder."
		header_extractor_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		output_path = output_descr_extract_entry.get()
		if output_path[-1] != '/':
			output_path += '/'
		verify_filepath_exists(output_path)
		header_extractor_log.info("Output folder path provided by user: '%s'" % output_path)

	# figures if user wants replace comma with user provided symbol or character
	# comma_state = replace_comma_checkboxval.get()
	# if comma_state:
	# 	header_extractor_log.info("User has requested to replace comma in %s with %s" % (keyword, symbol_replacing_comma))
	# else:
	# 	header_extractor_log.info("User has requested NOT to replace comma in %s with %s" % (keyword, symbol_replacing_comma))

	# determines if user want to extract partial description
	partial_descr_state = partial_descr_checkboxval.get()
	header_extractor_log.info("State of checkbox 'Extract Partial %s': %s" %(keyword, str(bool(partial_descr_state))))
	if partial_descr_state:
		reg_exp = reg_exp_descr_entry.get()
		header_extractor_log.info('Regular expression provided by user: %s' % reg_exp)

	protein_descriptions = []
	not_conforming_titles = []
	for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
		title = seq_record.description

		title_temp = ''
		for no in range(0, len(title)):
			if title[no] == ' ':
				title_temp = title[no+1::]		# if description is to be extracted
				break

		if id_option:							# if identifer is to be extracted
			title_temp = title[0:no]

		# this part extracts complete description from fasta header
		if not partial_descr_state:
			protein_descriptions.append(title_temp)

		# this part extracts partial description from fasta header based on user provided regular expression
		# If headers dont match regular expression provided, they will be recorded as 'Header requirement not fuflilled'
		# in output csv file. Also informs user with details if requirements are not met.
		else:
			try:
				m = re.match(r"%s" % reg_exp, title_temp)
				# m = re.match(r"(.+) \[.+", title_temp)
				description = m.group(1)
				protein_descriptions.append(description)
			except AttributeError as exception:
				not_conforming_titles.append(title)
				description = 'Header requirement not fuflilled'
				protein_descriptions.append(description)
				header_extractor_log.error("Requirement error '%s' arose when processing header: '%s'" % (exception, title))

	# warns user if regular expression requirements are not met in any among fasta headers
	if len(not_conforming_titles) and partial_descr_state:
		num = protein_descriptions.count('Header requirement not fuflilled')
		not_conf_string = ''
		for ele in not_conforming_titles:
			not_conf_string += (ele + '\n')

		temp = (' Make sure sequence header in fasta file conforms to requirements provided.'
		'\n Problem arose when reading following fasta header(s): \n %s' %not_conf_string)
		header_extractor_log.error(temp)
		popup_error(temp)

		temp = (' Total number of non-conforming headers: %i' % num +
				"\n Such instances will be saved as 'Header requirement not fuflilled' in output csv file." +
				"\n '%s Extractor' will continue now." % keyword)
		header_extractor_log.info(temp)
		popup_error(temp)

	else:
		header_extractor_log.info('All sequence headers in fasta file conform to the requirements.')

	unique_descriptions = set(protein_descriptions)
	header_extractor_log.debug('Number of sequences found in fasta file: %i' % len(protein_descriptions) +
							'Number of unique %ss found in fasta file:  %i' % (keyword, len(unique_descriptions)))

	# this part gathers count of each description
	Dict_unique_descriptions = {k: [] for k in unique_descriptions}
	for descr in unique_descriptions:
		Dict_unique_descriptions[descr] = protein_descriptions.count(descr)

	# this part enables to get unique descriptions and their count in descending order.
	data_csv = ''
	for key, value in sorted(Dict_unique_descriptions.iteritems(), key = lambda (k,v): (v,k), reverse = True):
		# if comma_state:						# if user wants to replace comma with other symbol or character
		# 	key = key.replace(',', symbol_replacing_comma)		# symbol_replacing_comma is defined in settings text file and can be edited from GUI settings window

		if ',' in key or '"' in key:		# This is to be consistent with csv filetype.
			key = key.replace('"','""')		# double quotes need to be double-double quoted if they need to appear when opening with spreadsheet software
			key = ('"%s"' % key)			# strings within double quote will show up as single celled value

		data_csv += (key + ',' + str(value) + '\n')

	# this part writes protein descriptions in to a csv file
	if len(protein_descriptions):
		output_file = output_path + ('%ss_Extracted.csv' % keyword)
		header_extractor_log.info("%ss will be written into output file: '%s'" % (keyword, output_file))

		# raises error if csv file is open in MS Excel
		try:
			Output_handle = open(output_file, 'w')
		except IOError as exception:
			if exception.errno != errno.EEXIST:
				temp =  ("A file titled '%s' seems to be open in MS Excel. Close that file and try again!" %output_file)
				header_extractor_log.error(temp)
				header_extractor_log.error('Error that resulted: %s' % exception)
				popup_error(temp)
				raise_enabler('stop')

		# opens and writes data
		Output_handle.write('Data below is obtained from fasta file:  "%s"\n' %input_fasta_file +
							'Total number of sequences present in fasta file:  ' + str(len(protein_descriptions)) + '\n' +
							'Number of unique %ss present among them:  ' % keyword + str(len(unique_descriptions)) + '\n\n' +
							'### Following are the unique %ss extracted: ###\n\n' % keyword +
							'%s, Count\n' % keyword + data_csv)
		Output_handle.close()

		header_extractor_log.info('Data writing in to above csv file was completed.')

		temp = (" Total number of sequences present in fasta file:  %i" % len(protein_descriptions) +
				"\n Number of unique %ss present among them:  %i" % (keyword, len(unique_descriptions)))
		header_extractor_log.info(temp)
		tkMessageBox.showinfo('Job Done', temp)

	else:
		temp = ('%s was not found in fasta file provided. ' % keyword +
			   'Make sure title of fasta sequence records are in required format.')
		header_extractor_log.error(temp)
		popup_error(temp)

	descr_extractor_Submit.configure(state=ACTIVE)
	descr_extract_processing.grid_remove()
	header_extractor_log.info('Job done! Ready for next job.\n\n')

	log_path = output_path  + 'log.txt'
	shutil.copyfile(logger_filename, log_path)


# Function that reads fasta sequences and filters (excludes or includes) them based on user-provided Descriptions-list.
# Descriptions must be separated by double commas or provided in a file.
def fasta_filter_by_description():
	global fasta_filename_entry
	global namelist_4filtering_entry
	global namelist_filename_entry
	global checkboxval_fileoption
	global fasta_filter_by_name_Submit
	global output_filter_name_entry
	global filter_name_processing
	global reg_exp_descr_filter_entry
	global partial_dscrp_boxval
	global include_or_exclude_names

	descr_filter_log.info("User has requested tool 'Filter fasta by sequences' Description'.")
	fasta_filter_by_name_Submit.configure(state=DISABLED)
	filter_name_processing.grid()

	if include_or_exclude_names.get() == 1:
		include_names = True
		exclude_names = False
		descr_filter_log.info("User chose 'Include Descriptions in list'")
	else:
		include_names = False
		exclude_names = True
		descr_filter_log.info("User chose 'Exclude Descriptions in list'")

	# reads if user wants to use partial descriptions for filtering
	partial_state = partial_dscrp_boxval.get()
	descr_filter_log.info("State of 'Filter using Partial Descriptions': %s" % str(bool(partial_state)))

	# reads regular expression provided by user
	if partial_state:
		reg_expr = reg_exp_descr_filter_entry.get()
		descr_filter_log.info("Regular expression provided by user: %s" % reg_expr)

	# this part recognises if user has provided Descriptions in a file or as a text in the box and then processes accordingly
	if not checkboxval_fileoption.get():		# if Descriptions are provided in textbox
		names_to_be_filtered = namelist_4filtering_entry.get()
		descr_filter_log.info("Descriptions provided as string by user: '%s'" % names_to_be_filtered)

		# informs if Description list box is empty
		if not names_to_be_filtered:
			temp = "'List of Descriptions' box is empty."
			descr_filter_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')

		## informs if any Description has '_$' (as that might be possble artifact from using Description extractor)
		# elif '_$' in names_to_be_filtered:
		# 	temp = "Characters '_$' are present in Descriptions. Make sure you remove them if you don't want them. " \
		# 		   "Program will continue when you click ok"
		# 	descr_filter_log.error(temp)
		# 	popup_error(temp)

		names_to_be_filtered = names_to_be_filtered.split(',,')		# makes it a list
		total_names = len(names_to_be_filtered)

	else:										# if names are provided in text file
		# reads user-provided Descriptions containing text file
		if not namelist_filename_entry.get():
			temp = "'File with Descriptions' box is empty. Select a file."
			descr_filter_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')
		else:
			list_filename = namelist_filename_entry.get()
			descr_filter_log.info('Text file (containing Descriptions) provided by user:\n\t%s' % list_filename)

		try:
			list_file = open(list_filename, 'Ur')
		except Exception as e:
			temp = "Error arose when reading 'Descriptions-file' provided."
			descr_filter_log.error(temp)
			descr_filter_log.error('Error that resulted: %s' % e)
			popup_error(temp)
			raise_enabler('stop')

		# obtains Descriptions from text file
		names_to_be_filtered = []
		for line in list_file:
			line = line.replace('\n', '')
			names_to_be_filtered.append(line)
		total_names = len(names_to_be_filtered)

		if not names_to_be_filtered:		# raises error if Descriptions were not found
			temp = (" No Description has been found in file:\n %s." % list_filename)
			descr_filter_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')

		# ## informs if a Description has '_$' (as that might be possble artifact from using Description extractor of ResCons)
		# check = 0
		# check = [check+1 for name in names_to_be_filtered if '_$' in name]
		# if check:
		# 	temp = " Warning: Characters '_$' are present in at least one Description in your file. \n " \
		# 		   "Make sure you remove them if you don't want them. \n " \
		# 		   "Program will continue as usual when you click ok"
		# 	descr_filter_log.warning(temp)
		# 	popup_error(temp)

	descr_filter_log.info("Descriptions recognised by software: %s" % names_to_be_filtered)

	# Obtains fasta file provided by user
	if not fasta_filename_entry.get():
		temp = "'FASTA file' box is empty. Choose a file."
		descr_filter_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		input_handle_fasta = fasta_filename_entry.get()
		descr_filter_log.info("FASTA file provided by user: %s" % input_handle_fasta)

	# Obtains output folder provided by user
	if not output_filter_name_entry.get():
		temp = "'Output folder' box is empty. Choose a folder."
		descr_filter_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		output_path_filter = output_filter_name_entry.get()
		if output_path_filter[-1] != '/':
			output_path_filter += '/'
		verify_filepath_exists(output_path_filter)
		descr_filter_log.info("Output folder provided by user: %s" % output_path_filter)


	# this part verifies if provided Descriptions match to fasta seq Descriptions and extract them as user-requested
	extracted_seqs = []
	header_error = []
	header_error_count = 0
	total_seqs = 0
	extracted = 0
	for seq_record in SeqIO.parse(input_handle_fasta, "fasta"):
		title = seq_record.description

		title_temp = ''
		for no in range(0, len(title)):
			if title[no] == ' ':
				title_temp = title[no+1::]
				break

		# if complete description is used for fasta filtering
		if not partial_state:
			name = title_temp

		# if partial description is used for fasta filtering
		else:
			try:											# if open sqaure bracket is present after protein's Description
				title_re = re.match(r"%s" % reg_expr, title_temp)
				# title_re = re.match(r"(.+) \[.+", title_temp)
				name = title_re.group(1)
			except AttributeError as exception:				# if open sqaure bracket is not present after protein's Description
				# name = title_temp
				name = 'Header requirement not fuflilled'
				header_error.append(title)
				header_error_count += 1
				descr_filter_log.error("Requirement error '%s' arose when processing header: '%s'" % (exception, title))

		if include_names :		# if user wants to include sequences corresponding to matching Descriptions
			if name in names_to_be_filtered:
				extracted_seqs.append(seq_record)
				extracted += 1
			total_seqs +=1

		elif exclude_names:	# if user wants to extract sequences of Descriptions that don't match
			if name not in names_to_be_filtered:
				extracted_seqs.append(seq_record)
				extracted += 1
			total_seqs +=1

	# provides details to user if requirements are not met in fasta headers provided
	if header_error_count:
		temp = (' Make sure sequence header in fasta file conforms to requirements provided.'
		'\n Problem arose when reading following fasta header(s): \n "%s".' % header_error)
		descr_filter_log.error(temp)
		popup_error(temp)

		temp = (' Total number of non-conforming headers: %i' % header_error_count +
				"\n 'Filter fasta by sequences name' will continue now.")
		descr_filter_log.info(temp)
		popup_error(temp)
	else:
		descr_filter_log.info('All sequence headers in fasta file seemed to conform to requirements.')


	output_file_filter = output_path_filter +'Filtered_Sequences_by_description.fasta'
	output_handle_filter = open(output_file_filter, 'w')
	SeqIO.write(extracted_seqs, output_handle_filter, 'fasta')
	output_handle_filter.close()

	temp = (' Number of sequences present in fasta file:  %i' % total_seqs +
			'\n Number of names provided by user:  %i' % total_names +
			'\n Number of sequence records extracted and written in output file:  %i' % extracted +
			'\n Output file: %s' %output_path_filter)
	descr_filter_log.info(temp)
	tkMessageBox.showinfo('Job Done!', temp)

	descr_filter_log.info('Job done! Ready for next job.\n\n')

	log_path = output_path_filter + 'log.txt'
	shutil.copyfile(logger_filename, log_path)

	fasta_filter_by_name_Submit.configure(state=ACTIVE)
	filter_name_processing.grid_remove()


# Function that reads fasta sequences and filters (excludes or includes) them based on user-provided sequence ID.
# IDs must be provided in a file.
def fasta_filter_by_id():
	global fasta_filename2_entry
	global idlist_filename_entry
	global include_or_exclude_IDs
	global partial_id_checkval
	global fasta_filter_by_id_Submit
	global output_idlist_entry
	global filter_id_processing

	idfilter_log.info("User has requested tool 'Filter fasta by sequences ID'.")
	fasta_filter_by_id_Submit.configure(state=DISABLED)
	filter_id_processing.grid()

	include_or_exclude= include_or_exclude_IDs.get()
	if include_or_exclude == 1:
		include_IDs = True
		exclude_IDs = False
		idfilter_log.info("User chose 'Include IDs in list'.")
	else:
		include_IDs = False
		exclude_IDs = True
		idfilter_log.info("User chose 'Remove IDs in list'.")


	partial_ids = partial_id_checkval.get()
	idfilter_log.info("State of 'IDs are Partial': %s" % str(bool(partial_ids)))

	# reads user-provided text file with IDs
	if not idlist_filename_entry.get():
		temp = "'File with IDs' box is empty. Select a file."
		idfilter_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		list_filename = idlist_filename_entry.get()
		idfilter_log.info('Text file (containing IDs) provided by user:\n\t%s' % list_filename)

	# reads user-provided output folder
	if not output_idlist_entry.get():
		temp = "'Output folder' box is empty. Select a folder."
		idfilter_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		output_path = output_idlist_entry.get()
		if output_path[-1] != '/':
			output_path += '/'
		verify_filepath_exists(output_path)
		idfilter_log.info('Output folder path provided by user:\n\t%s' % output_path)


	try:
		list_id_file = open(list_filename, 'Ur')
	except Exception as e:
		temp = ("Error arose when reading file: %s" % list_filename)
		idfilter_log.error(temp)
		idfilter_log.error('Error that resulted: %s' % e)
		popup_error(temp)
		raise_enabler('stop')

	ids_to_be_filtered = []
	for line in list_id_file:
		line = line.replace('\n', '')
		line_temp = line
		for no in range(1, len(line)+1):		# strips off pipe symbols in the end of ID.
			if line[-no] == '|':
				line_temp = line[:(-no)]
			else:
				break
		if line_temp:
			ids_to_be_filtered.append(line_temp)
	ids_total = len(ids_to_be_filtered)

	if not ids_to_be_filtered:					# raises error if names were not found
		temp = (" No ID has been found in file:\n %s." % list_filename)
		idfilter_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')

	idfilter_log.info("IDs recognised by software: %s" % ids_to_be_filtered)

	# this dictionary will assist in extracting seqs in order provided in text file by user
	dict_ids = {k: [] for k in ids_to_be_filtered}

	# Obtains fasta file provided by user
	if not fasta_filename2_entry.get():
		temp = "'FASTA file' box is empty. Choose a file."
		idfilter_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		input_file_fasta = fasta_filename2_entry.get()

	# this part verifies if provided IDs match to fasta seq IDs and extract them as user-requested
	extracted_seqs_by_id = []
	total_seqs = 0
	extracted=0
	for seq_record in SeqIO.parse(input_file_fasta, "fasta"):
		title_temp = seq_record.id

		title = title_temp
		for no in range(1, len(title_temp)+1):	# strips off pipe symbol in the end of seq ID.
			if title_temp[-no] =='|':
				title = title_temp[:(-no)]
			else:
				break

		if include_IDs:		# if user wants to include IDs corresponding to matching IDs
			if not partial_ids:
				if title in ids_to_be_filtered:
					# extracted_seqs_by_id.append(seq_record)
					dict_ids[title].append(seq_record)
					extracted += 1
			else:								# if IDs entered by user are partial IDs; not complete
				# if any(partial_id in title for partial_id in ids_to_be_filtered):
				# 	extracted_seqs_by_id.append(seq_record)
				# 	extracted += 1

				for partial_id in ids_to_be_filtered:
					if partial_id in title:
						dict_ids[partial_id].append(seq_record)
						extracted += 1
						break

			total_seqs +=1

		elif exclude_IDs:	# if user wants to extract IDs that don't match
			if not partial_ids:
				if title not in ids_to_be_filtered:
					extracted_seqs_by_id.append(seq_record)
					# dict_ids[title].append(seq_record)
					extracted += 1
			else:								# if IDs entered by user are partial IDs; not complete
				if not any(partial_id in title for partial_id in ids_to_be_filtered):
						extracted_seqs_by_id.append(seq_record)
						extracted += 1
			total_seqs +=1


	# writes extracted data into an output file when user has checked only one checkbox
	# output_path = os.path.dirname(input_file_fasta)
	# output_path += '/Output/'
	# verify_filepath_exists(output_path)

	Output_filename_idfilter = output_path + 'Filtered_Sequences_based_on_IDs.fasta'
	Output_handle_idfilter = open(Output_filename_idfilter, 'w')

	if exclude_IDs:
		SeqIO.write(extracted_seqs_by_id, Output_handle_idfilter, 'fasta')
	else:
		for key in ids_to_be_filtered:
			SeqIO.write(dict_ids[key], Output_handle_idfilter, 'fasta')
	Output_handle_idfilter.close()

	temp = (' Number of sequences present in fasta file:  %i' % total_seqs +
			'\n Number of IDs present in IDs-text file:  %i' % ids_total +
			'\n Number of matching sequence records extracted and written in output file:  %i' % extracted +
			'\n Output file: %s' %Output_filename_idfilter)
	idfilter_log.info(temp)
	tkMessageBox.showinfo('Job Done!', temp)

	idfilter_log.info('Job done! Ready for next job.\n\n')

	log_path = output_path + 'log.txt'
	shutil.copyfile(logger_filename, log_path)

	fasta_filter_by_id_Submit.configure(state=ACTIVE)
	filter_id_processing.grid_remove()



# Dummy function!
def fn_to_pass():
	pass


# Function to hide or show further options based on checkbox status
def hide_expand(checkboxval, frame):
	if not checkboxval.get():
		frame.grid_remove()
	else:
		frame.grid()


# Function to hide or show widgets depending on clustal alignment is requested or not
def hide_expand_clustal(Clade_checkboxval, frame_clade):
	if not Clade_checkboxval.get():
		frame_clade.grid_remove()
		seq_entry.configure(state = DISABLED)
		button_seq.configure(state = DISABLED)
		label_seq.configure(state = DISABLED)
		output_main_entry.delete(0, END)
	else:
		frame_clade.grid()
		seq_entry.configure(state = NORMAL)
		button_seq.configure(state = NORMAL)
		label_seq.configure(state = NORMAL)
		output_main_entry.delete(0, END)


# Function to make a label in tkinter
def label_text(win_name, label_text, height, width, row, column):
	label_t = Label(win_name, text=label_text, height=height, width=width, anchor=W)
	label_t.grid(row=row, column=column, sticky=W)
	return label_t


# Function to make a checkbox in tkinter
def checkbox_fn(appl, string, row_no, col_no, command_name):
	checkboxval = BooleanVar()
	checkboxval.set(False)

	checkbox = Checkbutton(appl, variable=checkboxval, text=string, command=command_name)
	checkbox.grid(row=row_no, column=col_no, sticky=W)
	return checkbox, checkboxval


# Function to make a entry box in tkinter
def entry_box(win_name, var_name, width, row, column):
	templ_entry = Entry(win_name, textvariable=var_name, width=width)
	templ_entry.grid(row=row, column=column, sticky=W)
	return templ_entry


def button_color(win_name, browse_fn, row, column):
	button_1 = Button(win_name, text='Color', width=6, command=browse_fn)
	button_1.grid(row=row, column=column, padx=4, pady=4)
	return button_1


# Function to make a browse button in tkinter
def browse_button(win_name, browse_fn, row, column):
	button_1 = Button(win_name, text='Browse', width=6, command=browse_fn)
	button_1.grid(row=row, column=column, padx=10, pady=10)
	return button_1


# Function that dictates type of input files and inserts selected file path into entry box
def browse_for_file(entry_name, filetype):
	File_path = askopenfilename(filetypes=filetype, title= 'Select a file')
	if File_path:
		entry_name.delete(0, END)
		entry_name.insert(0, File_path)
		entry_name.xview_moveto(1)


# Function that dictates type of input files and inserts selected file path into entry box and output folder
def browse_for_file_change_outpath(entry_name, out_entry, filetype):
	File_path = askopenfilename(filetypes=filetype, title= 'Select a file')
	if File_path:
		entry_name.delete(0, END)
		entry_name.insert(0, File_path)
		entry_name.xview_moveto(1)

		dir_path = os.path.dirname(File_path)
		dir_path += '/Output/'
		out_entry.delete(0, END)
		out_entry.insert(0, dir_path)
		out_entry.xview_moveto(1)


# Function to insert selected folder path into entry box
def browse_for_directory(entry_name):
	File_path = askdirectory(title= 'Select a directory to save output files')
	if File_path:
		entry_name.delete(0, END)
		entry_name.insert(0, File_path)
		entry_name.xview_moveto(1)


# Function to have a label-like text box with copy-able text
def default_text(frame, default, width, row, column):
	def1 = Text(frame, height=1, width = width)
	def1.insert(END, default)
	def1.grid(row=row, column =column, padx= 10)
	return def1


# Function that highlights and shows the frame of corresponding tool selected in more tools window.
def select_tools(frame_name, button):
	# global frame_blast, frame_genbank, frame_clade, frame_descr_extractor, frame_fasta_name_filter
	# global button_blast_filter, button_genbank, button_clades, button_descr_extractor, button_namefilter
	frames_list = [frame_blast, frame_genbank, frame_clade, frame_descr_extractor, frame_fasta_name_filter, frame_fasta_id_filter]
	buttons_list = [button_blast_filter, button_genbank, button_clades, button_descr_extractor, button_namefilter, button_idfilter]

	for frame in frames_list:
		if frame == frame_name:
			frame.grid()
		else:
			frame.grid_remove()

	for button_no in buttons_list:
		if button_no == button:
			if mac_os:
				button_no.configure(width = 40, font=('times', '16', 'bold', 'italic'))	# for MAC
			else:
				button_no.configure(bg = 'steelblue', fg = 'white', font=('times', '11', 'italic'))	# for win and ubunutu

		else:
			if mac_os:
				button_no.configure(width= 32, font=('times', '14'))				# for MAC
			else:
				button_no.configure(bg = '#F8F8F8', fg = 'black', font=('times', '10'))		# for win and ubunutu


# Function that creates right click menu and its contents (thanks to stackoverflow)
def make_menu(w):
	global right_click_menu
	right_click_menu = Menu(w, tearoff=0)
	right_click_menu.add_command(label="Cut")
	right_click_menu.add_command(label="Copy")
	right_click_menu.add_command(label="Paste")


# Function to show right-click menu when needed and what clicking on each would do
def show_menu(e):
	w = e.widget
	right_click_menu.entryconfigure("Cut", command=lambda: w.event_generate("<<Cut>>"))
	right_click_menu.entryconfigure("Copy", command=lambda: w.event_generate("<<Copy>>"))
	right_click_menu.entryconfigure("Paste", command=lambda: w.event_generate("<<Paste>>"))
	right_click_menu.tk.call("tk_popup", right_click_menu, e.x_root, e.y_root)


# Function that reads inputs from user
def read_user_input():
	global resi_positions
	global SeqQuery_file
	global Reference_file
	global Reference
	global SeqQuery_FileName
	global Aligned_Filename
	global Output_Path
	global protein_or_dna_mode
	global protein_mode

	# Obtains file path of Reference file
	# Reference_file = "/Users/Mana/Dropbox/ResCons/Reference.fasta"		# remove this line

	protein_mode = seq_mode()
	if protein_mode:
		runscript_log.info("ResCon will run in 'Protein mode' as user has requested.")
	else:
		runscript_log.info("ResCon will run in 'DNA/RNA mode' as user has requested.")

	if not templ_entry.get():			#uncomment this block, when done
		temp = "'Reference file' box is empty. Choose a file."
		runscript_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		Reference_file = templ_entry.get()
		runscript_log.info("'Reference file provided by user: \n\t%s" % Reference_file)

	if not output_main_entry.get():
		temp = "'Output folder' box is empty. Choose a folder."
		runscript_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		Output_Path = output_main_entry.get()
		if Output_Path[-1] != '/':
			Output_Path += '/'
		verify_filepath_exists(Output_Path)
		runscript_log.info("Output folder path provided by user: %s" % Output_Path)

	if checkboxval_Clustal.get():		# This if if user has requested clustal alignment
		# SeqQuery_file = "/Users/Mana/Dropbox/ResCons/seqs.fasta"			# remove this line when done

		if not seq_entry.get():		#uncomment this block, when done
			temp = "'Sequences file' box is empty. Choose a file."
			runscript_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')
		else:
			SeqQuery_file = seq_entry.get()
			gui_log.info('Sequences File that will be used for alignment:  %s' % SeqQuery_file)

		# this part verifies if seq file path contains single/double quote as clustalo hates them
		if "'" in SeqQuery_file or '"' in SeqQuery_file or '(' in SeqQuery_file or ')' in SeqQuery_file:
			temp = "Single or double quote or parenthesis is present in Sequences file's name or in its file path. " \
				   "Clustal omega will result in error in such cases. Try again!"
			runscript_log.error(temp)
			popup_error(temp)
			raise_enabler('stop')

		SeqQuery_FileName = os.path.basename(SeqQuery_file)

		# log_path = Output_Path + 'log.txt'
		# sys.stderr = open(log_path, 'w')

		gui_log.info("All output files will be stored at: %s" % Output_Path)

	else:								# This if if user has provided  clustal alignment file.
		clustal_log.info('Pre-Aligned file was provided by user.')
		Aligned_Filename = clustal_entry.get()
		clustal_log.info('Pre-Aligned file provided: \n\t%s.' % Aligned_Filename)

		gui_log.info("All output files will be stored at: \n\t%s." % Output_Path)


	# Reads residue positions entered by user
	resi_positions = resi_entry.get()

	if not resi_positions:
		gui_log.error('Residue positions cannot be empty!')
		popup_error('Residue positions cannot be empty!')
		raise_enabler('stop')


	# Verifies it is made only of numbers, commas, spaces and tabs
	resi_pos_test = re.findall(r'[^0-9\, \t\r\n]', resi_positions)
	if resi_pos_test:
		gui_log.error('Invalid character(s) found in residue positions entered!')
		popup_error('Invalid character(s) found in residue positions entered!')
		raise_enabler('stop')

	resi_positions = resi_positions.replace(' ', '')
	resi_positions = resi_positions.split(',')
	resi_positions = [int(x) for x in resi_positions]


	# verifies if zero or duplicate numbers are present in user-provided residue positions list
	if 0 in resi_positions:
		temp = 'Why would you have Zero (0) in residue positions?! Try again!'
		gui_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	elif len(resi_positions) > len(set(resi_positions)):
		gui_log.info('User-provided residue positions: %s' % resi_positions)
		resi_positions = list(set(resi_positions))

		temp = "Duplicate numbers present in residue positions provided. \n  Click 'Yes' if you would like ResCon to " \
			   "remove duplicates and proceed further. \n  Click 'No' to stop further execution."
		gui_log.error(temp)
		ans = tkMessageBox.askquestion('Error', temp, default = 'yes')
		if ans == 'no':
			gui_log.info("User clicked 'No'. ResCon will stop now.")
			raise_enabler('stop')
		else:
			gui_log.info("User clicked 'Yes'. ResCon proceeds further.")

		gui_log.info('Updated residue positions after removing duplicate(s): %s' % resi_positions)

	try:
		Reference = SeqIO.read(Reference_file, "fasta")
	except (ValueError, IOError) as Exception:
		temp = ("Error reading Reference file:  '%s'. "
				"Make sure it is in fasta format and has only one sequence." %Reference_file)
		runscript_log.error(temp)
		runscript_log.error('Error that resulted: %s' % Exception)
		popup_error(temp)
		raise_enabler('stop')

	if '-' in Reference.seq:
		temp = "Gap character '-' is found in your Reference sequence.\n If this is a mistake, click 'Cancel' and " \
			   "ResCon will stop.\n If this is acceptable, click 'OK' and ResCon will proceed"
		runscript_log.warning(temp)
		ans = tkMessageBox.askokcancel('Warning', temp, default= 'cancel')
		if ans:
			runscript_log.info("User clicked 'OK'. ResCon proceeds further.")
		else:
			runscript_log.info("User clicked 'Cancel'. ResCon stops now.")
			raise_enabler('stop')

	if max(resi_positions) > len(Reference.seq):
		temp = ("Error in 'Residue positions' entered. At least one of the residue positions entered "
				"is greater than %i (length of reference sequence)." % len(Reference.seq))
		runscript_log.error(temp)
		popup_error(temp)
		raise_enabler('stop')
	else:
		temp = "All residue positions entered are less than the length of reference sequence."
		runscript_log.info(temp)


# Function to run the main script by calling appropriate functions, based on user input
def run_script():
	global clustal_user_command_string
	global Output_Path

	runscript_log.info('ResCons begins!')
	processing.grid()
	button_run_script.configure(state = DISABLED)		# Disables submit job button while processing
	read_user_input()

	if checkboxval_Clustal.get():		# if clustal alignment is requested
		clustal_log.info('Clustal alignment was requested.')
		clustal_user_command_string = clustal_command.get()
		# button_run_script.configure(state = DISABLED)		# Disables submit job button while processing

		clustal_alignment()
		runscript_log.info('Clustal alignment was completed!')

	runscript_log.debug("Module 'fetch_mismatch' called.")
	# button_run_script.configure(state = DISABLED)		# Disables submit job button while processing
	fetch_mismatch()
	runscript_log.debug("Exited module 'fetch_mismatch' after successful processing.")

	if checkboxval_formatting.get():	# if html formatting of clustal alignment is requested
		gui_log.info('User has requested to format alignment in html format')
		runscript_log.debug("Module 'html_formatting' called.")
		button_run_script.configure(state = DISABLED)		# Disables submit job button while processing

		runscript_log.info('User requested HTML color formatting for clustal alignment')
		html_formatting()
		runscript_log.debug("Exited module 'html_formatting' after successful processing.")

	processing.grid_remove()
	button_run_script.configure(state = ACTIVE)		# Re-enables submit job button while processing

	log_path = Output_Path + 'log.txt'
	shutil.copyfile(logger_filename, log_path)

	tkMessageBox.showinfo('Job Done', "Done. All output files were saved in folder: '%s'." % Output_Path)
	runscript_log.info('All requested jobs are done! Ready for next job!!\n\n')


app = Tk()
app.title("ResCon")
app.geometry("")
app.resizable(width=FALSE, height=FALSE)	# doesn't allow resizing window. This disables maximize button as well
title_bar_logo = '128x128_part_no_bg_logo.ico'
if win_os:			# to show logo in windows task bar and in app border top.
	app.wm_iconbitmap(bitmap=title_bar_logo)


make_menu(app)		# this initiates right-click menu which will pop-up when right-clicked.

if mac_os:
	label_title = Label(app, text='Mismatch analyzer', height= 2, justify=CENTER, font=('Times', '16', 'bold'))
else:
	# label_title = Label(app, text='Mismatch analyzer', height=1, justify=CENTER, font=('Times', '14', 'bold'))
	label_title = Label(app, text='Mismatch analyzer', height=2, justify=CENTER, font=('Times', '14', 'bold'))
label_title.grid(row=0, column=0, sticky= N)  # Didn't use label_text function here as it forces text to be in  left side

frame1 = Frame(app)
frame1.grid(row=2, column=0,sticky='w')

frame_mode = Frame(frame1)		# frame to house protein and dna/rna mode for user's choosing
frame_mode.grid(row=1, column=0, columnspan=3)

if mac_os:
	font_size = 16
else:
	font_size = 12

# Function that determines if user wants ResCon to run in protein mode or dna/rna mode
def seq_mode():
	if protein_or_dna_mode.get() == 1:
		protein_mode = True
		mode1.configure(font=('times', font_size, 'italic'))
		mode2.configure(font=('times', font_size))
	else:
		protein_mode = False
		mode1.configure(font=('times', font_size))
		mode2.configure(font=('times', font_size, 'italic'))
	return protein_mode

protein_or_dna_mode = IntVar()
protein_or_dna_mode.set(1)
mode1 = Radiobutton(frame_mode, text="Protein mode", variable=protein_or_dna_mode, value=1, command = seq_mode, font=('times', font_size, 'italic'))
mode1.grid(row=0, column= 1, padx= 15, sticky = E+W)

mode2 = Radiobutton(frame_mode, text="DNA/RNA mode", variable=protein_or_dna_mode, value=2, command = seq_mode, font=('times', font_size))
mode2.grid(row=0, column= 2, padx=15, sticky = E+W)


label_seq = label_text(frame1, "   Sequences file", 2, 20, 2, 0)
# label_seq.configure(fg='red')

seq_filename = StringVar()
seq_entry = entry_box(frame1, seq_filename, 50, 2, 1)

# this binds the rigt-click menu to entry boxes. This needs to be done only once as bind_class is used.
if mac_os:
	seq_entry.bind_class("Entry", "<Button-2><ButtonRelease-2>", show_menu) 	#button-2 is right click for mac
else:
	seq_entry.bind_class("Entry", "<Button-3><ButtonRelease-3>", show_menu)		#button-3 is right click for windows

# filetype_all = [('All files', '*.*')]
filetype_all = []
filetype_fasta = [('fasta files', '*.fasta'), ('All files', '*.*')]
# filetype_alignment = [('clustal files', '*.aln'), ('clustal files', '*.clustal'), ('All files', '.*')]
filetype_alignment = [('clustal files', '*.aln *.clustal *.clu'), ('fasta files', '*.fa *.fasta *.a2m'),
					  ('vienna files', '*.vie *.vienna'), ('phylip files', '*.phy *.phylip'), ('All files', '.* *.txt')]

# button_seq = browse_button(frame1, lambda: browse_for_file(seq_entry, filetype_fasta), 2, 2)
button_seq = browse_button(frame1, lambda: browse_for_file_change_outpath(seq_entry, output_main_entry, filetype_fasta), 2, 2)

label_templ = label_text(frame1, "   Reference file", 2, 20, 3, 0)
# label_templ.configure(fg='blue')

templ_filename = StringVar()
templ_entry = entry_box(frame1, templ_filename, 50, 3, 1)

button_templ = browse_button(frame1, lambda: browse_for_file(templ_entry, filetype_fasta), 3, 2)
button_templ.grid(padx=25)

label_res_pos = label_text(frame1, "   Residue positions", 2, 20, 6, 0)

resival = StringVar()
resival.set('')
resi_entry = entry_box(frame1, resival, 50, 6, 1)
resi_entry.configure(bg='#FFFF94')

output_main_label = label_text(frame1, "   Output folder", 0, 20, 5, 0)

output_main_val = StringVar()
output_main_val.set('')
output_main_entry = entry_box(frame1, output_main_val, 50, 5, 1)

button_output_main = browse_button(frame1, lambda: browse_for_directory(output_main_entry), 5, 2)


frame2 = Frame(app)
frame2.grid(row=7, column=0, sticky='w')

frame_clustal = Frame(frame2, bd=3, relief=GROOVE)
frame_clustal.configure(pady=10)
frame_clustal.grid(row=9, column=2, sticky='w')
frame_clustal.grid_remove()

checkbox_Clustal, checkboxval_Clustal = checkbox_fn(frame2, 'Clustal Alignment required?', 8, 1,
													lambda: hide_expand_clustal(checkboxval_Clustal, frame_clustal))
checkboxval_Clustal.set(True)
checkbox_Clustal.grid(sticky=E + W)  # to center the label in  window

label_clustal = label_text(frame2, "   Alignment file", 0, 20, 9, 0)

clustal_filename = StringVar()
clustal_entry = entry_box(frame2, clustal_filename, 50, 9, 1)

button_clustal = browse_button(frame2, lambda: browse_for_file_change_outpath(clustal_entry, output_main_entry, filetype_alignment), 9, 2)
button_clustal.grid(padx= 25)

frame_clustal = Frame(frame2, bd=3)
frame_clustal.grid(row=9, column=0, columnspan=10, sticky='w')
# frame_clustal.grid_remove()

clustal_command_label = label_text(frame_clustal, "   Clustal-O Command", 0, 20, 0, 0)

clustal_command_value = StringVar()
clustal_command_value.set(clustalo_command_default)
# clustal_command_value.set("{'outfmt': 'clu', 'iterations': 3, 'residuenumber': 'True', 'force': 'True', 'outfile': 'default', 'guidetree_out': 'default.newick'}")
# clustal_command_value.set("{'outfmt': 'clu', 'auto':'True', 'force': 'True', 'outfile': 'default'}")
clustal_command = entry_box(frame_clustal, clustal_command_value, 50, 0, 1)

dummy = label_text(frame_clustal, '', 0, 20, 0, 2)
dummy.configure(width = 10)
dummy.grid(padx = 15)

checkbox_formatting, checkboxval_formatting = checkbox_fn(frame2, 'HTML Color-coding required?', 12, 1, fn_to_pass)
checkboxval_formatting.set(True)
checkbox_formatting.grid(sticky=E + W)

if linux_os:
	button_run_script = Button(frame2, text='Submit job', width=16, command=lambda: gen_thread(processing, run_script),
							   font=('times', '12', 'bold'), bg='steelblue', fg='white')
else:
	button_run_script = Button(frame2, text='Submit job', width=16, command=run_script, font=('times', '12', 'bold'),
						  bg='steelblue', fg='white')
button_run_script.grid(row=14, column=1, pady=20)
if mac_os:
	button_run_script.configure(font=('times', '14', 'bold'))			# for Mac


if Image_in_GUI:
	# image = Image.open("128x128_no_bg_logo.png")			# adds logo image to the GUI
	image = Image.open("128x128_part_no_bg_logo.png")			# adds logo image to the GUI
	# image = image.resize((170,110), Image.ANTIALIAS)
	photo1 = ImageTk.PhotoImage(image)
	xxx= Label(frame2, image = photo1)
	xxx.grid(row=14, column=0, rowspan = 5, columnspan=2, sticky = W, padx=20)
	# xxx.image(photo1)


# Labels that will show up when ResCon is under execution
global processing, processing_clustal
processing = label_text(frame2, 'Processing...', 0, 20, 15, 1)
processing.configure(fg='red', font= ("times", 16))
processing.grid(sticky= E+W)
processing.grid_remove()

processing_clustal = label_text(frame2, 'Clustal alignment under progress...', 0, 20, 15, 1)
processing_clustal.configure(fg='red', font= ("times", 16),  width=20)
processing_clustal.grid(sticky= E+W)
processing_clustal.grid_remove()


# Function that draws GUI for 'more tools'
def top_win():
	global top
	top = Toplevel(app)
	if win_os:			# to show logo in windows task bar and in app border top.
		top.wm_iconbitmap(bitmap=title_bar_logo)

	# disables root/main window. This causes problem in windows os, and probably in mac, though. When gui is minimized,
	# clicking on it will not open the software again. It can be activated again with task manager meaning the software
	# didn't freeze but it just can't open once minimized with child window open.
	# if mac_os or ubuntu_os:
	# 	top.grab_set()

	top.transient(app)		# keeps toplevel window always on top of root window
	top.geometry("")
	top.resizable(width=FALSE, height=FALSE)	# doesn't allow resizing window. This disables maximize button as well
	top.title("ResCon: More Tools")

	# MoreTools_label = Label(top, text = 'More Tools!!', height = 4)
	MoreTools_label = Label(top, text = '', height = 2)			# for spacing
	MoreTools_label.grid(row = 0, column = 0, columnspan = 2, sticky = EW)

	if Image_in_GUI:
		image = Image.open("128x128_part_no_bg_logo.png")			# adds logo image to the GUI
		image = image.resize((300, 195), Image.ANTIALIAS)
		photo2 = ImageTk.PhotoImage(image)
		yyy= Button(top, image = photo2, bd=0)
		yyy.image = photo2
		yyy.grid(row=0, column=0, rowspan= 5, columnspan= 5)

	# global variables for use in select_tools function
	global frame_blast, frame_genbank, frame_clade, frame_descr_extractor, frame_fasta_name_filter, frame_fasta_id_filter
	global button_blast_filter, button_genbank, button_clades, button_descr_extractor, button_namefilter, button_idfilter
	global Blast_Submit, genbank_Submit, Extract_Submit, descr_extractor_Submit, fasta_filter_by_name_Submit, fasta_filter_by_id_Submit
	global genbank_processing, blast_processing, subtree_processing, filter_id_processing, filter_name_processing, descr_extract_processing
	global output_blast_entry, output_genbank_entry, output_subtree_entry, output_descr_extract_entry, output_filter_name_entry, output_idlist_entry

	'''This part is GUI stuff for E-value based BLAST filter tool'''

	global Blast_xml, Blast_fasta, E_threshold_entry
	# global Higher_E_val_checkboxval, Lower_E_val_checkboxval
	global lower_or_higher_e_threshold
	global genbank_entry
	global Newick_entry, Fasta_Original_entry, BranchLength_entry

	frame_blast = Frame(top, bd=2, relief=GROOVE)
	frame_blast.grid(row = 5, column = 0, columnspan=2, sticky = 'w')
	frame_blast.grid_remove()

	button_blast_filter = Button(top, text = "Filter fasta by Blast's E-value", width =30,
								 command = lambda: select_tools(frame_blast, button_blast_filter))
	button_blast_filter.grid(row = 2, column = 1, padx=10, pady=10)
	if mac_os:
		button_blast_filter.configure(width= 32, font=('times', '14'))		# for Mac
	else:
		button_blast_filter.configure(bg = '#F8F8F8', fg = 'black', font=('times', '10'))

	Blast_xml_label = label_text(frame_blast, "   XML BLAST file   ", 2, 20, 3, 0)

	Blast_xml_filename = StringVar()
	Blast_xml = entry_box(frame_blast, Blast_xml_filename, 50, 3, 1)

	filetype_xml = [('xml files', '*.xml'), ('All files', '*.*')]
	Blast_xml_browse = browse_button(frame_blast, lambda: browse_for_file_change_outpath(Blast_xml, output_blast_entry, filetype_xml), 3, 2)

	Blast_fasta_label = label_text(frame_blast, "   FASTA file   ", 2, 20, 4, 0)

	Blast_fasta_filename = StringVar()
	Blast_fasta = entry_box(frame_blast, Blast_fasta_filename, 50, 4, 1)

	filetype_fasta = [('fasta files', '*.fasta'), ('All files', '*.*')]
	Blast_fasta_browse = browse_button(frame_blast, lambda: browse_for_file(Blast_fasta, filetype_fasta), 4, 2)

	output_blast_label = label_text(frame_blast, "   Output folder   ", 2, 20, 5, 0)

	output_blast_filename = StringVar()
	output_blast_entry = entry_box(frame_blast, output_blast_filename, 50, 5, 1)

	# filetype_all = []
	output_blast_browse = browse_button(frame_blast, lambda: browse_for_directory(output_blast_entry), 5, 2)

	E_threshold_label = label_text(frame_blast, "   E-value Threshold    ", 2, 20, 7, 0)

	E_threshold = StringVar()
	E_threshold_entry = entry_box(frame_blast, E_threshold, 20, 7, 1)
	E_threshold_entry.configure(bg='#FFFF94')

	lower_or_higher_e_threshold = IntVar()
	lower_or_higher_e_threshold.set(1)
	lower_e_radio = Radiobutton(frame_blast, text= 'Lower than E-threshold', variable= lower_or_higher_e_threshold, value=1)
	lower_e_radio.grid(row=8, column=1, sticky='w')

	higher_e_radio = Radiobutton(frame_blast, text= 'Higher than and equal to E-threshold', variable= lower_or_higher_e_threshold, value=2)
	higher_e_radio.grid(row=9, column=1, sticky='w')

	if linux_os:
		Blast_Submit = Button(frame_blast, text='Filter now!', width=12, command= lambda: gen_thread(blast_processing, blast_filter), bg='steelblue2')	# for linux
	else:
		Blast_Submit = Button(frame_blast, text='Filter now!', width=12, command= blast_filter, bg='steelblue2')
	Blast_Submit.grid(row=10, column=2, padx=10, pady=10)

	blast_processing = label_text(frame_blast, 'Processing...', 0, 20, 10, 1)
	blast_processing.configure(fg='red', font= ("times", 16))
	blast_processing.grid(sticky= E+W)
	blast_processing.grid_remove()


	'''This part is GUI stuff for GenPept/GenBank to fasta format converter'''
	global Newick_hate_sym_checkboxval, ignore_incomplete_checkboxval, fasta_id_length_entry, options_dict

	frame_genbank = Frame(top, bd=2, relief=GROOVE)
	frame_genbank.grid(row = 5, column = 0, columnspan=2, sticky = 'w')
	frame_genbank.grid_remove()

	button_genbank = Button(top, text = "GenPept/GenBank To fasta Converter", width =30, command = lambda: select_tools(frame_genbank, button_genbank))
	button_genbank.grid(row = 1, column = 1, padx=10, pady=10)
	if mac_os:
		button_genbank.configure(width= 32, font=('times', '14'))			# for Mac
	else:
		button_genbank.configure(bg = '#F8F8F8', fg = 'black', font=('times', '10'))

	genbank_label = label_text(frame_genbank, "   GenPept/GenBank file", 2, 20, 12, 0)

	genbank_filename = StringVar()
	genbank_entry = entry_box(frame_genbank, genbank_filename, 50, 12, 1)

	filetype_genbank = [('genbank files', '*.gp'),  ('genbank files', '*.gb'), ('All files', '*.*')]
	genbank_browse = browse_button(frame_genbank, lambda: browse_for_file_change_outpath(genbank_entry, output_genbank_entry, filetype_genbank), 12, 2)


	output_genbank_label = label_text(frame_genbank, "   Output folder", 2, 20, 13, 0)

	output_genbank_filename = StringVar()
	output_genbank_entry = entry_box(frame_genbank, output_genbank_filename, 50, 13, 1)

	output_genbank_browse = browse_button(frame_genbank, lambda: browse_for_directory(output_genbank_entry), 13, 2)


	fasta_id_length_label = label_text(frame_genbank, "   Fasta ID length", 2, 20, 15, 0)

	fasta_id_length_entry_val = StringVar()
	fasta_id_length_entry_val.set(fasta_id_length_default)	# value obtained from settings file
	fasta_id_length_entry = entry_box(frame_genbank, fasta_id_length_entry_val, 20, 15, 1)
	fasta_id_length_entry.configure(bg='#FFFF94')

	frame_options = Frame(frame_genbank, bd=2)
	frame_options.grid(row = 17, column = 0, columnspan=4, sticky = 'w')

	options_label = label_text(frame_options, '   Header options', 1, 20, 1, 0)

	# for checkbox values
	options_list = ['Locus Name', 'GI', 'Domain', 'Phylum', 'Class', 'Genus', 'Species', 'Seq Length', 'Seq Name']
	options_dict = {k: 1 for k in options_list}
	row = 1
	column = 1
	for option in options_list:
		temp, options_dict[option] = checkbox_fn(frame_options, option, row, column, fn_to_pass)
		options_dict[option].set(True)
		temp.configure(width = 15, justify = LEFT, anchor = W)
		column += 1
		if column == 4:
			column = 1
			row += 1

	options_dict['GI'].set(False)
	options_dict['Seq Length'].set(False)
	options_dict['Seq Name'].set(False)

	more_options_label = label_text(frame_genbank, '   Retrieval options', 2, 20, 18, 0)

	global Newick_hate_sym_checkbox
	Newick_hate_sym_checkboxval = BooleanVar()
	Newick_hate_sym_checkbox, Newick_hate_sym_checkboxval= checkbox_fn(frame_genbank,
												"Replace Newick-sensitive symbols in header id with '%s'" % newick_sym_replace,18, 1, fn_to_pass)
	# Newick_hate_sym_checkbox = Checkbutton(frame_genbank, variable=Newick_hate_sym_checkboxval, text=string, command=command_name)
	# Newick_hate_sym_checkbox.grid(row=row_no, column=col_no, sticky=W)

	Newick_hate_sym_checkboxval.set(True)
	Newick_hate_sym_checkbox.grid(columnspan= 3, sticky=W)

	ignore_incomplete_checkbox, ignore_incomplete_checkboxval= checkbox_fn(frame_genbank,
													'Ignore incomplete/partial sequence records', 19, 1, fn_to_pass)
	ignore_incomplete_checkboxval.set(True)
	ignore_incomplete_checkbox.grid(columnspan= 3, sticky=W)

	if linux_os:
		genbank_Submit = Button(frame_genbank, text='Convert now!', width=12, command= lambda: gen_thread(genbank_processing, genpept_converter), bg='steelblue2')		# for linux
	else:
		genbank_Submit = Button(frame_genbank, text='Convert now!', width=12, command= genpept_converter, bg='steelblue2')
	genbank_Submit.grid(row=22, column=2, padx=10, pady=10)

	genbank_processing = label_text(frame_genbank, 'Processing...', 0, 20, 22, 1)
	genbank_processing.configure(fg='red', font= ("times", 16))
	genbank_processing.grid(sticky= E+W)
	genbank_processing.grid_remove()


	'''This part is GUI stuff to Extract clades from Phylogenetic tree'''

	frame_clade = Frame(top, bd=2, relief=GROOVE)
	frame_clade.grid(row = 5, column = 0, columnspan=2, sticky = 'w')
	# frame_clade.grid_remove()

	button_clades = Button(top, text = "Subtree Sequences Extractor", width =30, command = lambda: select_tools(frame_clade, button_clades))
	button_clades.grid(row = 1, column = 0, padx=10, pady=10)
	if mac_os:
		button_clades.configure(width = 40, font=('times', '16', 'bold', 'italic'))		# for Mac
	else:
		button_clades.configure(bg = 'steelblue', fg= 'white',font=('times', '11', 'italic'))

	Newick_label = label_text(frame_clade, "   Newick file   ", 2, 20, 18, 0)

	Newick_filename = StringVar()
	Newick_entry = entry_box(frame_clade, Newick_filename, 50, 18, 1)

	filetype_tree = [('Newick files', '*.newick *.nwk'), ('Nexus files', '*.nexus *.nex *.nxs'), ('Text files', '*.txt'),
					   ('NeXML files', '*.xml *.nexml'), ('PhyloXML files', '*.xml *.phyloxml'), ('All files', '*.*')]
	newick_browse = browse_button(frame_clade, lambda: browse_for_file_change_outpath(Newick_entry,
															output_subtree_entry, filetype_tree), 18, 2)

	Fasta_Original_label = label_text(frame_clade, "   FASTA file    ", 2, 20, 19, 0)

	Fasta_Original_filename = StringVar()
	Fasta_Original_entry = entry_box(frame_clade, Fasta_Original_filename, 50, 19, 1)

	Fasta_Original_browse = browse_button(frame_clade, lambda: browse_for_file(Fasta_Original_entry, filetype_fasta), 19, 2)

	output_subtree_label = label_text(frame_clade, "   Output folder    ", 2, 20, 20, 0)

	output_subtree_filename = StringVar()
	output_subtree_entry = entry_box(frame_clade, output_subtree_filename, 50, 20, 1)

	output_subtree_browse = browse_button(frame_clade, lambda: browse_for_directory(output_subtree_entry), 20, 2)

	BranchLength_label = label_text(frame_clade, "   Clade's Branch Length    ", 2, 20, 22, 0)

	BranchLength = StringVar()
	BranchLength_entry = entry_box(frame_clade, BranchLength, 20, 22, 1)
	BranchLength_entry.configure(bg='#FFFF94')

	if linux_os:
		Extract_Submit = Button(frame_clade, text='Extract now!', width=12, command= lambda: gen_thread(subtree_processing, Extract_Clades), bg='steelblue2', fg='grey2')		# for linux
	else:
		Extract_Submit = Button(frame_clade, text='Extract now!', width=12, command= Extract_Clades, bg='steelblue2', fg='grey2')
	Extract_Submit.grid(row=23, column=2, padx=10, pady=10)

	subtree_processing = label_text(frame_clade, 'Processing...', 0, 20, 23, 1)
	subtree_processing.configure(fg='red', font= ("times", 16))
	subtree_processing.grid(sticky= E+W)
	subtree_processing.grid_remove()


	'''This part is GUI stuff for Fasta Description/ID Extractor'''
	global descr_extr_file_entry
	global replace_comma_checkboxval
	global reg_exp_descr_entry
	global partial_descr_checkboxval
	# global descr_option
	global id_option
	global replace_comma_checkbox

	frame_descr_extractor = Frame(top, bd=2, relief=GROOVE)
	frame_descr_extractor.grid(row = 5, column = 0, columnspan=2, sticky = 'w')
	frame_descr_extractor.grid_remove()

	button_descr_extractor = Button(top, text = "Fasta Description / ID Extractor", width =30, command = lambda: select_tools(frame_descr_extractor, button_descr_extractor))
	button_descr_extractor.grid(row = 3, column = 0, padx=10, pady=10)
	if mac_os:
		button_descr_extractor.configure(width= 32, font=('times', '14'))		# for Mac
	else:
		button_descr_extractor.configure(bg = '#F8F8F8', fg = 'black', font=('times', '10'))


	def choose_id_extr():
		global id_option
		id_option = True
		partial_descr_checkbox.configure(text='Extract Partial Identifier')
		# reg_exp_descr.set('(\w+)\|.+')
		reg_exp_descr.set(regex_id_extr)		# regex_id_extr is defined in settings text file
		if mac_os:
			id_choice.configure(text = "=> Identifier Extractor", width = 27, font=('times', '16', 'bold', 'italic'))
			descr_choice.configure(text = "Description Extractor", width= 20, font=('times', '14'))
		else:
			id_choice.configure(text = "=> Identifier Extractor", bg = '#00CC33', fg = 'black', font=('times', '11', 'italic'))
			descr_choice.configure(text = "Description Extractor", bg = '#F8F8F8', fg = 'black', font=('times', '10'))

	def choose_descr_extr():
		global id_option
		id_option = False
		partial_descr_checkbox.configure(text='Extract Partial Description')
		# reg_exp_descr.set('(.+) \[.+')
		reg_exp_descr.set(regex_desc_extr)		# regex_desc_extr is defined in settings text file
		if mac_os:
			descr_choice.configure(text = "=> Description Extractor", width = 27, font=('times', '16', 'bold', 'italic'))
			id_choice.configure(text = "Identifier Extractor", width= 20, font=('times', '14'))
		else:
			descr_choice.configure(text = "=> Description Extractor",bg = '#00CC33', fg = 'black', font=('times', '11', 'italic'))
			id_choice.configure(text = "Identifier Extractor", bg = '#F8F8F8', fg = 'black', font=('times', '10'))

	frame_choice = Frame(frame_descr_extractor, bd=2)
	frame_choice.grid(row = 26, column = 0, columnspan=4)

	id_option = False		# this variable identifies if identifier extraction is chosen by user


	descr_choice = Button(frame_choice, text = "=> Description Extractor", width =20, command = choose_descr_extr)
	descr_choice.grid(row = 0, column = 0, padx=10, pady=10)
	if mac_os:
		descr_choice.configure(width = 27, font=('times', '16', 'bold', 'italic'), relief= GROOVE)	# for MAC
	else:
		descr_choice.configure(bg = '#00CC33', fg = 'black', font=('times', '11', 'italic'), relief= GROOVE)

	id_choice = Button(frame_choice, text = "Identifier Extractor", width =20, command = choose_id_extr)
	id_choice.grid(row = 0, column = 1, padx=10, pady=10)
	if mac_os:
		id_choice.configure(width= 20, font=('times', '14'), relief= GROOVE)				# for MAC
	else:
		id_choice.configure(bg = '#F8F8F8', fg = 'black', font=('times', '10'), relief= GROOVE)

	descr_extractor_label = label_text(frame_descr_extractor, "   FASTA file   ", 2, 20, 27, 0)
	# descr_extractor_label.configure(width = 20, anchor = W)

	descr_extractor_filename = StringVar()
	descr_extr_file_entry = entry_box(frame_descr_extractor, descr_extractor_filename, 50, 27, 1)

	filetype_descr_extractor = [('fasta files', '*.fasta'), ('All files', '*.*')]
	descr_extractor_browse = browse_button(frame_descr_extractor,
				lambda:browse_for_file_change_outpath(descr_extr_file_entry, output_descr_extract_entry, filetype_descr_extractor), 27, 2)

	output_descr_extract_label = label_text(frame_descr_extractor, "   Output folder   ", 2, 20, 28, 0)

	output_descr_extract_filename = StringVar()
	output_descr_extract_entry = entry_box(frame_descr_extractor, output_descr_extract_filename, 50, 28, 1)

	output_descr_extract_browse = browse_button(frame_descr_extractor,
										  lambda:browse_for_directory(output_descr_extract_entry), 28, 2)

	frame_partial_descr = Frame(frame_descr_extractor)
	frame_partial_descr.grid(row=29, column=0, columnspan = 4, sticky='w')
	frame_partial_descr.grid_remove()

	label_text(frame_partial_descr, "   Regular Expression    ", 2, 20, 0, 0)

	reg_exp_descr = StringVar()
	# reg_exp_descr.set('(.+) \[.+')
	reg_exp_descr.set(regex_desc_extr)		# regex_desc_extr is defined in settings text file
	reg_exp_descr_entry = entry_box(frame_partial_descr, reg_exp_descr, 20, 0, 1)
	reg_exp_descr_entry.configure(bg='#FFFF94')

	partial_descr_checkbox, partial_descr_checkboxval = checkbox_fn(frame_descr_extractor, 'Extract Partial Description', 30, 1, lambda: hide_expand(partial_descr_checkboxval, frame_partial_descr))
	partial_descr_checkboxval.set(False)
	if not mac_os:
		partial_descr_checkbox.grid(sticky = E+W)				# to center the label in the window

	# replace_comma_checkbox, replace_comma_checkboxval = checkbox_fn(frame_descr_extractor, "Replace comma with '%s'" % symbol_replacing_comma, 31, 1, fn_to_pass)
	# replace_comma_checkbox.configure(height=2)
	# replace_comma_checkboxval.set(True)
	# if not mac_os:
		# replace_comma_checkbox.grid(sticky = E+W)				# to center the label in the window

	if linux_os:
		descr_extractor_Submit = Button(frame_descr_extractor, text = 'Extract now!', width =12,command= lambda: gen_thread(descr_extract_processing, description_extractor),
									   bg = 'steelblue2')						# for linux
	else:
		descr_extractor_Submit = Button(frame_descr_extractor, text = 'Extract now!', width =12,command= description_extractor,
								   bg = 'steelblue2')
	descr_extractor_Submit.grid(row = 32, column = 2, padx=10, pady=10)

	descr_extract_processing = label_text(frame_descr_extractor, 'Processing...', 0, 20, 32, 1)
	descr_extract_processing.configure(fg='red', font= ("times", 16))
	descr_extract_processing.grid(sticky= E+W)
	descr_extract_processing.grid_remove()


	'''This part is GUI stuff for filtering fasta sequences based on their Description'''
	global fasta_filename_entry
	global namelist_4filtering_entry
	global namelist_filename_entry
	global checkboxval_fileoption
	global reg_exp_descr_filter_entry
	global partial_dscrp_boxval
	global include_or_exclude_names

	frame_fasta_name_filter = Frame(top, bd=2, relief=GROOVE)
	frame_fasta_name_filter.grid(row = 5, column = 0, columnspan = 2, sticky = W)
	frame_fasta_name_filter.grid_remove()

	button_namefilter = Button(top, text = "Filter fasta by Sequences' Description", width =30,
							   command = lambda: select_tools(frame_fasta_name_filter, button_namefilter))
	button_namefilter.grid(row = 3, column = 1, padx=10, pady=10)
	if mac_os:
		button_namefilter.configure(width= 32, font=('times', '14'))			#for Mac
	else:
		button_namefilter.configure(bg = '#F8F8F8', fg = 'black', font=('times', '10'))

	fasta_label = label_text(frame_fasta_name_filter, "   FASTA file   ", 2, 20, 36, 0)

	fasta_filename = StringVar()
	fasta_filename_entry = entry_box(frame_fasta_name_filter, fasta_filename, 50, 36, 1)
	fasta_filename_entry.grid(columnspan = 2)

	filetype_fasta = [('fasta files', '*.fasta'), ('All files', '*.*')]
	fasta_browse = browse_button(frame_fasta_name_filter,
								 lambda:browse_for_file_change_outpath(fasta_filename_entry, output_filter_name_entry, filetype_fasta),36, 2)

	label_text(frame_fasta_name_filter, "   Output folder   ", 2, 20, 37, 0)

	output_filter_name_filename = StringVar()
	output_filter_name_entry = entry_box(frame_fasta_name_filter, output_filter_name_filename, 50, 37, 1)

	browse_button(frame_fasta_name_filter, lambda:browse_for_directory(output_filter_name_entry),37, 2)

	frame_namelist = Frame(frame_fasta_name_filter, bd=3)
	frame_namelist.grid(row = 39, column = 0, columnspan = 10, sticky = W)
	# frame_namelist.grid_remove()

	frame_name_in_file = Frame(frame_fasta_name_filter, bd=3)
	frame_name_in_file.grid(row = 39, column = 0, columnspan = 10, sticky = W)
	frame_name_in_file.grid_remove()

	namelist_4filtering_label = label_text(frame_namelist, "  List of Descriptions  ", 2, 20, 1, 0)

	namelist_4filtering =  StringVar()
	namelist_4filtering_entry = entry_box(frame_namelist, namelist_4filtering, 50, 1, 1)
	namelist_4filtering_entry.configure(bg = '#FFFF94')

	namefile_label = label_text(frame_name_in_file, "  File with Descriptions", 2, 20, 1, 0)
	# namefile_label.configure(width = 19)

	namelist_filename = StringVar()
	namelist_filename_entry = entry_box(frame_name_in_file, namelist_filename , 50, 1, 1)

	filetype_txt = [('txt files', '*.txt'), ('All files', '*.*')]
	filename_browse = browse_button(frame_name_in_file, lambda:browse_for_file(namelist_filename_entry, filetype_txt),1, 2)
	filename_browse.grid(padx = 30)

	frame_reg_ex = Frame(frame_fasta_name_filter, bd=3)
	frame_reg_ex.grid(row = 40, column = 0, columnspan = 3, sticky = W)
	frame_reg_ex.grid_remove()

	label_text(frame_reg_ex, "  Regular Expression", 2, 20, 0, 0)

	reg_exp_descr_filter = StringVar()
	# reg_exp_descr_filter.set('(.+) \[.+')
	reg_exp_descr_filter.set(regex_desc_fltr)	# regex_desc_fltr is defined in settings text file
	reg_exp_descr_filter_entry = entry_box(frame_reg_ex, reg_exp_descr_filter, 20, 0, 1)
	reg_exp_descr_filter_entry.configure(bg='#FEFFCB')

	frame_options = Frame(frame_fasta_name_filter)
	frame_options.grid(row = 41, column = 0, columnspan = 10, sticky = W)

	label_text(frame_options, "   ", 2, 20, 0, 0)

	checkbox_fileoption, checkboxval_fileoption = checkbox_fn(frame_options, 'Get Descriptions from text file',
												0, 2, lambda: hide_expand(checkboxval_fileoption, frame_name_in_file))
	# checkboxval_fileoption.set(True)

	partial_dscrp_box, partial_dscrp_boxval = checkbox_fn(frame_options,
														'Filter using Partial Descriptions', 1, 2, lambda: hide_expand(partial_dscrp_boxval, frame_reg_ex))

	include_or_exclude_names = IntVar()
	include_or_exclude_names.set(1)
	include_names_radio = Radiobutton(frame_options, text= 'Include Descriptions in list', variable= include_or_exclude_names, value=1)
	include_names_radio.grid(row= 0, column= 1)

	remove_names_radio = Radiobutton(frame_options, text= 'Exclude Descriptions in list', variable= include_or_exclude_names, value=2)
	remove_names_radio.grid(row= 1, column= 1)

	if linux_os:
		fasta_filter_by_name_Submit = Button(frame_fasta_name_filter, text = 'Filter now!', width =12,
											 command= lambda: gen_thread(filter_name_processing, fasta_filter_by_description), bg = 'steelblue2')		# for linux
	else:
		fasta_filter_by_name_Submit = Button(frame_fasta_name_filter, text = 'Filter now!', width =12,
											 command= fasta_filter_by_description, bg = 'steelblue2')
	fasta_filter_by_name_Submit.grid(row = 47, column = 2, padx=10, pady=10)

	filter_name_processing = label_text(frame_fasta_name_filter, 'Processing...', 0, 20, 47, 1)
	filter_name_processing.configure(fg='red', font= ("times", 16))
	filter_name_processing.grid(sticky= E+W)
	filter_name_processing.grid_remove()


	'''This part is GUI stuff for filtering fasta sequences based on their ID'''
	global fasta_filename2_entry
	global idlist_filename_entry
	global include_or_exclude_IDs
	global partial_id_checkval

	frame_fasta_id_filter = Frame(top, bd=2, relief=GROOVE)
	frame_fasta_id_filter.grid(row = 5, column = 0, columnspan = 2, sticky = W)
	frame_fasta_id_filter.grid_remove()

	button_idfilter = Button(top, text = "Filter fasta by Sequences' ID", width =30, command = lambda: select_tools(frame_fasta_id_filter, button_idfilter))
	button_idfilter.grid(row = 2, column = 0, padx=10, pady=10)
	if mac_os:
		button_idfilter.configure(width= 32, font=('times', '14'))			# for Mac
	else:
		button_idfilter.configure(bg = '#F8F8F8', fg = 'black', font=('times', '10'))

	fasta_label2 = label_text(frame_fasta_id_filter, "   FASTA file   ", 2, 20, 45, 0)

	fasta_filename2 = StringVar()
	fasta_filename2_entry = entry_box(frame_fasta_id_filter, fasta_filename2, 50, 45, 1)
	fasta_filename2_entry.grid(columnspan=2)

	filetype_fasta = [('fasta files', '*.fasta'), ('All files', '*.*')]
	fasta_browse2 = browse_button(frame_fasta_id_filter,
								  lambda:browse_for_file_change_outpath(fasta_filename2_entry, output_idlist_entry, filetype_fasta),45, 3)

	idfile_labelval = label_text(frame_fasta_id_filter, "   File with IDs", 2, 20, 46, 0)

	idlist_filename = StringVar()
	idlist_filename_entry = entry_box(frame_fasta_id_filter, idlist_filename, 50, 46, 1)
	idlist_filename_entry.grid(columnspan=2)

	filetype_txt = [('txt files', '*.txt'), ('All files', '*.*')]
	filename2_browse = browse_button(frame_fasta_id_filter, lambda:browse_for_file(idlist_filename_entry, filetype_txt),46, 3)
	filename2_browse.grid(padx = 30)

	label_text(frame_fasta_id_filter, "   Output folder", 2, 20, 47, 0)

	output_idlist_filename = StringVar()
	output_idlist_entry = entry_box(frame_fasta_id_filter, output_idlist_filename, 50, 47, 1)
	output_idlist_entry.grid(columnspan=2)

	browse_button(frame_fasta_id_filter, lambda:browse_for_directory(output_idlist_entry), 47, 3)

	include_or_exclude_IDs = IntVar()
	include_or_exclude_IDs.set(1)
	include_IDs_val = Radiobutton(frame_fasta_id_filter, text= 'Include IDs in list', variable=include_or_exclude_IDs, value=1, command=fn_to_pass)
	include_IDs_val.grid(row=49, column=1)

	remove_IDs_val = Radiobutton(frame_fasta_id_filter, text= 'Remove IDs in list', variable=include_or_exclude_IDs, value=2, command=fn_to_pass)
	remove_IDs_val.grid(row=50, column=1)

	partial_id_check, partial_id_checkval = checkbox_fn(frame_fasta_id_filter, 'IDs are Partial', 49, 2, fn_to_pass)

	if linux_os:
		fasta_filter_by_id_Submit = Button(frame_fasta_id_filter, text = 'Filter now!', width =12, command= lambda: gen_thread(filter_id_processing, fasta_filter_by_id), bg = 'steelblue2')		# for linux
	else:
		fasta_filter_by_id_Submit = Button(frame_fasta_id_filter, text = 'Filter now!', width =12, command= fasta_filter_by_id, bg = 'steelblue2')
	fasta_filter_by_id_Submit.grid(row = 51, column = 3, padx=10, pady=10)

	filter_id_processing = label_text(frame_fasta_id_filter, 'Processing...', 0, 20, 51, 1)
	filter_id_processing.configure(fg='red', font= ("times", 16))
	filter_id_processing.grid(sticky= E+W)
	filter_id_processing.grid_remove()


	# Following dummy variables are for spacing before and after the tools' frames
	dummy = Label(top)
	dummy.grid(row = 4, column = 0, sticky = 'w')

	dummy = Label(top)
	dummy.grid(row=9, column=0, sticky='w')

	dummy = Label(top)
	dummy.grid(row = 30, column = 0, sticky = 'w')


button_moretools = Button(frame2, text='More Tools', width=13, command=top_win, font=('times', '12', 'italic'),
						  bg='steelblue1')
button_moretools.grid(row=16, column=1, pady=10)
if mac_os:
	button_moretools.configure(font=('times', '15', 'italic'))			# for Mac

dummy = Label(frame2, height=1, width=16)  # This dummy is just to align other widgets in their proper positions in gui
dummy.grid(row=18, column=0, sticky='w')


# This shows a pop-up when unexpected errors happen
Tk.report_callback_exception = show_error

# Adding file and help menu to GUI
menubar = Menu(app)


# Function that reads settings changed by user through GUI
def edit_settings(settings_win):
	# global match_color, similar_color, mismatch_color, aa_set_str, aa_set, id_delimiter, connector_id, newick_sym_replace, symbol_replacing_comma
	global match_color, similar_color, mismatch_color, aa_set_str, aa_set, id_delimiter, connector_id, newick_sym_replace
	global dict_box_values
	global Newick_hate_sym_checkbox, replace_comma_checkbox

	# settings for mismatch analyzer
	match_color = dict_box_values['match_color_edit'].get()
	similar_color = dict_box_values['similar_color_edit'].get()
	mismatch_color = dict_box_values['mismatch_color_edit'].get()
	id_delimiter = dict_box_values['id_delimiter_edit'].get()

	aa_set_str = dict_box_values['aa_set_str_edit'].get()			# gets similar amino acid group set and parses them
	line_list = aa_set_str.split(',')
	for no in range(0, len(line_list)):
		line_list[no] = line_list[no].replace(' ', '')

	aa_set_upper = [x.upper() for x in line_list]	# makes upper and lower case as a set
	aa_set_lower = [x.lower() for x in line_list]
	aa_set = aa_set_upper + aa_set_lower
	aa_set_str = ', '.join(aa_set_upper)

	# settings for GenPept/GenBank to fasta converter
	connector_id = dict_box_values['connector_id_edit'].get()



	old_newick_sym_replace = newick_sym_replace		# to help with error-checking
	newick_sym_replace = dict_box_values['newick_sym_replace_edit'].get()

	# settings for Fasta Description / ID Extractor
	# old_symbol_replacing_comma = symbol_replacing_comma		# to help with error-checking
	# symbol_replacing_comma = dict_box_values['symbol_replacing_comma_edit'].get()

	# this changes GUI text based on user's changes in settings.
	# works only when 'more tools' windows is open.
	# if newick_sym_replace != old_newick_sym_replace or symbol_replacing_comma != old_symbol_replacing_comma:
	if newick_sym_replace != old_newick_sym_replace:
		try:
			Newick_hate_sym_checkbox.configure(text= "Replace Newick-sensitive symbols in header id with '%s'" %newick_sym_replace)
			# replace_comma_checkbox.configure(text= "Replace comma with '%s'" % symbol_replacing_comma)
			settings_win.destroy()		# closes window if 'more tools' window is open or else throws an error
		except NameError as e:
			newick_sym_replace = old_newick_sym_replace
			# symbol_replacing_comma = old_symbol_replacing_comma
			temp = "Open 'More Tools' windows and then save your changes."
			gui_log.error(temp)
			tkMessageBox.showinfo('Warning', temp)
	else:
		settings_win.destroy()		# closes window


# title_list = ['Color:  Matching', 'Color:  Mismatching but similar', 'Color:  Mismatching and dissimilar', 'Delimiter used in FASTA IDs', 'Amino acids similarity sets', 'Connecting delimiter used in FASTA ID', 'Symbol replacing sensitive symbols', 'Replace comma symbol with']
# default_values = [match_color, similar_color, mismatch_color, id_delimiter, aa_set_str, connector_id, newick_sym_replace, symbol_replacing_comma]
# box_values = ['match_color_edit', 'similar_color_edit', 'mismatch_color_edit', 'id_delimiter_edit', 'aa_set_str_edit', 'connector_id_edit', 'newick_sym_replace_edit', 'symbol_replacing_comma_edit']
title_list = ['Color:  Matching', 'Color:  Mismatching but similar', 'Color:  Mismatching and dissimilar', 'Delimiter used in FASTA IDs', 'Amino acids similarity sets', 'Connecting delimiter used in FASTA ID', 'Symbol replacing sensitive symbols']
default_values = [match_color, similar_color, mismatch_color, id_delimiter, aa_set_str, connector_id, newick_sym_replace]
box_values = ['match_color_edit', 'similar_color_edit', 'mismatch_color_edit', 'id_delimiter_edit', 'aa_set_str_edit', 'connector_id_edit', 'newick_sym_replace_edit']


# Function that makes GUI interface for 'Edit Settings' window available through menubar
def settings_win_fn():
	global dict_box_values, settings_win
	settings_win = Toplevel(app)
	settings_win.transient(app)		# keeps toplevel window always on top of root window
	settings_win.geometry("")
	if win_os:			# to show logo in windows task bar and in app border top.
		settings_win.wm_iconbitmap(bitmap=title_bar_logo)
	settings_win.resizable(width=FALSE, height=FALSE)	# doesn't allow resizing window. This disables maximize button as well
	settings_win.title("ResCon: Edit settings")

	frame1_set = Frame(settings_win)
	frame1_set.grid(row=0, column=0, sticky='w')

	headings_list = ['Settings for Mismatch analyzer', 'Settings for GenPept/GenBank to fasta converter', 'Settings for Fasta Description / ID Extractor']
	# default_values_edited = [match_color, similar_color, mismatch_color, id_delimiter, aa_set_str, connector_id, newick_sym_replace, symbol_replacing_comma]
	default_values_edited = [match_color, similar_color, mismatch_color, id_delimiter, aa_set_str, connector_id, newick_sym_replace]
	dict_box_values = {k: '' for k in box_values}
	dict_entry_box = {k: '' for k in box_values}	# to assist in coloring entry box with respective color

	# function that fills in the respective entry box with color code selected
	def choose_color(key, colour):
		color = askcolor(color= colour)
		try:
			if color[1]:
				dict_box_values[key].set(color[1])
				dict_entry_box[key].configure(bg=color[1])
		except Exception as e:
			temp = 'ResCon has trouble obtaining color code. Enter the code manually'
			gui_log.error(temp)
			popup_error(temp)
			gui_log.error('Error that occured: %s' % e)

	row = 0
	col = 0
	check_no =0
	for num in range(0, len(title_list)):
		if num in [0,5,7]:		# to write title for each section
			title = label_text(frame1_set, headings_list[check_no], 2, 43, row, 0)
			if mac_os:
				title.configure(font=('times', '16', 'bold', 'italic'))
			else:
				title.configure(font=('times', '13', 'bold', 'italic'))
			title.grid(columnspan = 3)
			check_no += 1
			row+=1

		label_text(frame1_set, title_list[num], 0, 37, row, col).grid(padx=10)

		if num == 0:
			button_color(frame1_set, lambda: choose_color('match_color_edit', match_color), row, col).grid(sticky = E)
		elif num == 1:
			button_color(frame1_set, lambda: choose_color('similar_color_edit', similar_color), row, col).grid(sticky = E)
		elif num == 2:
			button_color(frame1_set, lambda: choose_color('mismatch_color_edit', mismatch_color), row, col).grid(sticky = E)

		dict_box_values[box_values[num]] = StringVar()
		dict_box_values[box_values[num]].set(default_values_edited[num])	# value obtained from settings file
		dict_entry_box[box_values[num]] = entry_box(frame1_set, dict_box_values[box_values[num]], 15, row, col+1)
		if num in [0,1,2]:
			dict_entry_box[box_values[num]].configure(bg = default_values_edited[num])	# for respcetive background colors

		label_text(frame1_set, ' Default:', 0, 7, row, col+2)
		text_box = default_text(frame1_set, default_values[num], 15, row, col+3)
		if win_os:
			text_box.configure(bg=settings_win.cget('bg'), relief=FLAT)
			text_box.configure(state="disabled")	# doesn't let selecting text when used in Mac

		row += 1

	label_text(frame1_set, '', 0, 7, 20, 0)		# for spacing

	# this enables to open 'settings file' if changes need to be made permanent
	temp = "Note: To make your changes permanent, change values in 'Settings' file by clicking here!"
	permanent = label_text(frame1_set, temp, 0, 75, 55, 0)		# for spacing
	permanent.configure(fg = 'blue')
	permanent.grid(columnspan= 5)
	permanent.bind("<Button-1>", lambda event: open_file(settings_file_name))	# opens with left mouse button

	if mac_os:
		button_save = Button(frame1_set, text='Save', height=1, command=lambda: edit_settings(settings_win), font=('times', '16', 'bold'))
		button_cancel = Button(frame1_set, text='Cancel', height=1, command=lambda: settings_win.destroy(), font=('times', '16', 'bold'))
	else:
		button_save = Button(frame1_set, text='Save', height=0, width=6, command=lambda: edit_settings(settings_win), font=('times', '12', 'bold'))
		button_cancel = Button(frame1_set, text='Cancel', height=0, width=6, command=lambda: settings_win.destroy(), font=('times', '12', 'bold'))
	button_save.grid(row=50, column=2, columnspan=2, sticky = W)
	button_cancel.grid(row=50, column=3, columnspan=2)


# function that displays text when help -> about button is clicked
def about():
	tkMessageBox.showinfo('About ResCon', 'Author:  Manavalan Gajapathy, Ng Lab \n'
								   'Copyright (c) 2015:  GNU Lesser General Public License')


# function to open a file externally in default application in any OS
def open_file(filename):
	if win_os and not ubuntu_os:
		# os.system('start ' +filename)	# doesn't work if finepath has spaces in it
		os.startfile(filename, 'open')		# works fine even if finepath has spaces in it
	elif ubuntu_os:
		# os.system("xdg-open "+ filename)	# apparently os.system is relatively security riskier
		subprocess.call(('xdg-open', filename))
	elif mac_os:
		# os.system("open "+ filename)		# apparently os.system is relatively security riskier
		subprocess.call(('open', filename))


# File menu with Settings and Quit button
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="Edit Settings", command=settings_win_fn)
filemenu.add_command(label="Open Log File", command=lambda: open_file(logger_filename))
filemenu.add_command(label="Quit ResCon", command=app.quit)
menubar.add_cascade(label="File", menu=filemenu)

# Help menu with About button
helpmenu = Menu(menubar, tearoff=0)
helpmenu.add_command(label="User Guide", command= lambda: open_file('User_Guide_ResCon.pdf'))
helpmenu.add_command(label="About ResCon", command=about)
menubar.add_cascade(label="Help", menu=helpmenu)

# displays the menubar
app.config(menu=menubar)

app.lift()				# opens app in the front of other open windows
app.mainloop()			# Executes gui


