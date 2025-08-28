'''This script takes a list of species .
In the name of each species the genus and specific epithet must be separated by an underscore _.
There must be a directory with the name of the species containing the .fa or .fasta files.
These directories must be located under dirfasta0 (currently dirfasta0='seqs/').

It searches for all .fa or .fasta files and, for each one, creates a command line to perform
the segmentation with siglevel (e.g. =0.95).

It runs it with the alphabet (e.g. -A 4 -bo) and creates the output file in the directory
dirout0 (currently segments/), creating one directory per species (creates the directory
with -p to avoid error if it already exists).

The name of the output file is the same as the fasta file but with the extension .segXX, 
where XX is the significance value (in %).

The list of command lines (commands) is stored in a file cmdlist='cmd_parallel.txt', and with 
a system call it is passed to parallel, then the file is deleted.

From there, the process is monitored with htop from the system.

You can perform the segmentation at different significance levels by including
values in the array siglevel_list'''


import os
import sys
import time
import shutil

######################################################################

#check the presence of the segmentation program scc
if not (shutil.which("scc") or shutil.which("./scc")):
    print('''ERROR: 'scc' required
    
    ''')
    exit(0)
#--------------------------


#from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor,as_completed
#import subprocess as sp

#list_of_species = [
#    "Ornithorhynchus_anatinus",
#    "Pan_paniscus",
#    "Sus_scrofa",
#]
#siglevel_list=('0.91','0.92','0.93','0.94','0.95','0.96','0.97','0.98')
#dirfasta0='seqs/'
#dirout0='segments/'


# Reads parameters from file
from read_input_file import initial_data
siglevel_list,bin_list,dirfasta0,dirout0,histodir,list_of_species,\
  MAX_PROC,dirdistances = initial_data('input_file.txt')
dirfasta0=dirfasta0+'/'
dirout0=dirout0+'/'
# NOT used here
# bin_list, histodir,dirdistances

#right now this is the only option
alfa='-A 4 -bo'

comandos=[]

#-------------- generate a command line for the segmentation of each contig -------

for siglevel in siglevel_list:
	extens='.seg'+siglevel.split('.')[-1]
	############ list of species to segment ######################
	for species in list_of_species: 
		dirfasta=dirfasta0+species+'/'
		dirout=dirout0+species+'/'
		print(f's={siglevel} Alfabeto= {alfa} {dirfasta} --> {dirout}')
		os.system(f'mkdir -p {dirout}')
		for i in [x for x in os.listdir(dirfasta) if (('.fa' in x) or ('.fasta' in x))]:
			fichfasta=dirfasta+i
			fout=dirout+i.replace('.fa',extens).replace('.fasta',extens)
			comandos.append(f"./scc {fichfasta} {siglevel} {siglevel} 1 {alfa} -nsize 10 -s {fout} > /dev/null 2>&1 ")
#-------------------------------------------------------------------------------------------------------			

#-------------- write command lines to file -------------------
cmdlist='cmd_parallel.txt'
f=open(cmdlist,'w',encoding="utf-8")
for i,cmd in enumerate(comandos):
	f.write(cmd+'\n')
f.flush()
f.close()
#--------------------------------------------------------	
	
#--------------- Run the command lines with parallel ------------------
#----------------- (fire and forget) ---------------------

#----------------------
if not shutil.which("parallel"):
    print('''
ERROR: 'parallel' is  not installed
    Please, install it with "sudo apt install parallel"
    
    Alternatively, you can run manually the command lines
    stored in '''+cmdlist+'\n')
    exit(0)
#--------------------------

max_proc = MAX_PROC
os.system(f'cat {cmdlist} | parallel -j {max_proc} &' )
time.sleep(3)
os.system(f'rm {cmdlist}')	
print('******************')
print(f'{len(comandos)} segmentations launched in background')
print('You can also monitor the processes with top or htop')
#-----------------------------------------------------------

#----- control of launched processes ------------------

import psutil
import time
#import keyboard


nproc=-1
nproc0=-1
message='\n\nAll segmentations finished\n'
try:
	while nproc:
		nproc=0
		pids=[]
		for proc in psutil.process_iter(attrs=['cmdline','pid']):
			try:
				cmd=list(proc.info['cmdline'])[0]
			except:
				cmd=''	
			if './scc' in cmd:
				nproc+=1
				#print(nproc)
				pids.append(proc.info['pid'])
				#PID=proc.info['pid']
				#os.system(f'kill {PID}')
		if (nproc!=nproc0):
			print(f'\rProcesses running {nproc} ... (Ctr+C to stop)  ',end='')
		nproc0=nproc		
		time.sleep(10)
except KeyboardInterrupt:
	#print("\n⚡ Interrupción detectada con Ctrl-C")
	message="\n\nAll segmentations stopped\n"
	pass
finally:
	#for PID in pids:
		#os.system(f'kill {PID}')
	print(message)


