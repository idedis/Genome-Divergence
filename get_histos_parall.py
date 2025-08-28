
''' Este obtiene a partir de un directorio donde están todas las segmentaciones de los
trozos de una especie, los histogramas de g+c, g+a (purinas) y g+t (keto)
lo hace tb para una lista de valores de numero de bines
Los escribe cada uno en un fichero en el directorio (ver abajo)

EXT_SEG='.seg95' indica qué segmentación se utiliza

Paraleliza directamente en python: get_histos_parall.py. Está pensado para
ejecutarlo desde
genomesignature
Este toma una lista de especies (se puede tomar por ejemplo del fichero 
lista_mammals.py) y para todas ellas calcula los
histogramas, con los bines que hay en la lista y los guarda en un directorio
con el nombre de la especie pero que puede
estar en otra ubicación (ahora ./histos)

OJO! que si se empieza desde cero, es conveniente hacer un cálculo con un solo
valor de bines para que calcule sin interferencias los ficheros totales
Cuidado cuando se rehacen los cambios, si se quiere empezar de cero hay que borrar
los ficheros totales


This script takes, from a directory containing all the segmentations of the fragments
of a species, the histograms of g+c, g+a (purines), and g+t (keto). 

It also does this for a list of values of number of bins. Each histogram is written
to a file in the directory (see below). The list of species, number of bins and
significance levels are read from  'input_file.txt'

It parallelizes directly in Python.
The calculation of the first value of nbins is a bit slower because it gathers
all the segments an put in a global file. The rest of values of nbins make use of this
file

'''




import numpy as np
import os
import sys
#import matplotlib.pyplot as plt
from scipy.stats import kstest,chisquare
from scipy.spatial.distance import jensenshannon
from scipy.stats import entropy
import subprocess as sp
from multiprocessing import Pool
from multiprocessing import Process, Pipe

#import matplotlib.pyplot as plt

##############################################################
def genome_data(genome_dir):
	#el directorio contiene las segmentaciones de todos
	#los ficheros de una especie (pueden ser los cromosomas
	#o trozos de ellos)
	#Analiza los ficheros que contengan 'seg' en su nombre
	
	directorio=genome_dir #'seg99_Carlito_syrichta'
	lista_fich=os.listdir(directorio)
	tam=[]
	gc=[]
	pur=[]
	ke=[]
	frec=[]
	NTOT=0
	for fname in [x for x in lista_fich if EXT_SEG in x]:   #e.g. EXT_SEG='.seg95'
		#print(fname)
		fseg=open(directorio+'/'+fname,'r')
		nsegs=int(fseg.readline().split()[3])	 #no sirve cuando hay enes
		for i in fseg:
			NTOT=NTOT+1
			l=[int(x) for x in i.replace('\n','').split()]
			if sum(l[3:7])!=0:  #si no son enes
				tam.append(l[2])
				gc.append(float(l[5]+l[6])/float(l[2]))
				pur.append(float(l[3]+l[6])/float(l[2]))	
				ke.append(float(l[4]+l[6])/float(l[2]))	
		fseg.close()
	return gc,pur,ke,tam
###########################################################
def escribe_total(genome_dir,totfile):
	#vuelca a un único fichero todos los gc,pur,ke y tam
	#escribe en un fichero que se llama total_{siglevel}.dat en
	#el mismo directorio

	gc,pur,ke,tam=genome_data(genome_dir)
	
	f=open(totfile,'w')
	for i in range(len(gc)):
		cad=f'{gc[i]}\t{pur[i]}\t{ke[i]}\t{tam[i]}\n'
		f.write(cad)
#############################################################
def calc_histo(genome_dir,tipo=(0,),nbins=50,rango=(0,1),siglevel='95',reset=False):
	#if reset=True eliminates old totfile
	#0=GC, 1=Purinas, 2=Keto  #para las longitudes hago otro
	sbins=f'{int(nbins):03d}' #lo convierte a cadena para añadir al nombre del fichero
	
	#a este hay que pasarle la significacion para que el fichero total
	#sea distinto para cada significacion
	totfile=genome_dir+f'/total_{siglevel}.dat'
	
	if (reset):
		os.system(f'rm -f {totfile}')
		print(f'Deleting old {totfile}')
	
	#exit(0)
	try:
		f=open(totfile,'r')
	except:	
		print(f'Creating {totfile}')
		escribe_total(genome_dir,totfile)
	f=open(totfile,'r')
	histo=[]
	data0=[x for x in f]
	
	f.close()
	for t in tipo:	
		data=[float(x.split('\t')[t]) for x in data0]
		h,bins=np.histogram(data,bins=int(nbins),range=rango,density=True)	
		histo.append(h)
	
	return histo,bins
#############################################################		


########################################################################
#----- esta es la que se llama en paralelo -----------------------------
def write_histo(genome_dir,nbins,out_dir,siglevel='95',reset=False):
	
	if out_dir=='': out_dir=genome_dir  #si es en blanco se manda al mismo dir
	os.system(f'mkdir -p {out_dir}')    #si no existe se crea el directorio de salida
	
	snbins=f'{int(nbins):03d}'  #convertido a cadena de 3 caracteres
	salida=genome_dir.split('/')[-1]+f'\tbins={nbins:3d}\t -->  \t {out_dir}'
	salida0=salida.split('-->')[0] #if less information is required
	
	h,b = calc_histo(genome_dir,(0,1,2),int(nbins),rango=(0,1),siglevel=siglevel,reset=reset)  #0=GC,1=RY,2=KM,3=tam
	
	g=open(out_dir+f'/histo{siglevel}_cg{snbins}.dat','w')
	g.write(f'# Histogram of GC content {genome_dir} segments at 95%\n')
	for i,j in zip(b[1:],h[0]):
		g.write(f'{i:12.4f}\t{j:12.4f}\n')
	g.close()	
	salida=salida+f'\thisto{siglevel}_cg{snbins}.dat'
	#----
	g=open(out_dir+f'/histo{siglevel}_ry{snbins}.dat','w')
	g.write(f'# Histogram of G+A (purine)  content {genome_dir} segments at 95%\n')
	for i,j in zip(b[1:],h[1]):
		g.write(f'{i:12.4f}\t{j:12.4f}\n')
	g.close()	
	salida=salida+f'\thisto{siglevel}_ry{snbins}.dat'
	#----	
	g=open(out_dir+f'/histo{siglevel}_km{snbins}.dat','w')
	g.write(f'# Histogram of G+T (keto)  content {genome_dir} segments at 95%\n')
	for i,j in zip(b[1:],h[2]):
		g.write(f'{i:12.4f}\t{j:12.4f}\n')
	g.close()	
	salida=salida+f'\thisto{siglevel}_km{snbins}.dat'
	
	return salida0
#-----------------------------------------------------------------------
########################################################################	




from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor,as_completed

#list_of_species = [
#    "Ornithorhynchus_anatinus",
#    "Pan_paniscus",
#    "Sus_scrofa",
#]
#
#histodir='histos'
#segdir='segments'
#bin_list=('025','050','075','100','125','150','175','200','300','400','500')	      
#sig_list=('0.91','0.92','0.93','0.94','0.95','0.96','0.97')


# Reads parameters from file
from read_input_file import initial_data
sig_list,bin_list,dirfasta0,segdir,histodir,list_of_species,\
  MAX_PROC,dirdistances = initial_data('input_file.txt')
#fini=open('input_file.txt')
#info=[]
#for line in [x for x in fini if not(x.startswith('#'))]:
#	info.append(line.replace('\n',''))
#sig_list=info[0].split(',')
#bin_list=info[1].split(',')
#dirfasta0=info[2] #not used here
#segdir=info[3]
#histodir=info[4]
#list_of_species=info[5].split(',')
#MAX_PROC=info[6]
#fini.close()

###################################


for siglevel in [x.split('.')[1] for x in sig_list]:
	EXT_SEG='.seg'+siglevel #extension for segmentation files
	
	print('***********************')
	print(f'Starting with {siglevel}% significance level')
	print('')
	params=[] #es la lista de parámetros para cada cálculo de histograma, hay uno para cada
	          #especie y numero de bines porque los tres cg,ry y km se hacen juntos
	
		      
	
	for genome_dir in list_of_species: 
		for nbins in bin_list: 
			par=[f'{segdir}/'+genome_dir,int(nbins),f'{histodir}/{genome_dir}',siglevel,False]
			params.append(par)
		
			
	NTASKS=len(params)
	NTASKS0=NTASKS
	
	
	param_first=[]
	param_rest=[]
	for i in params:
		if int(bin_list[0])==i[1]:   #if the number of bins is the first value
			i[4]=True #deletes the previously existing totfile
			param_first.append(i)
		else:
			param_rest.append(i)
	
			
	#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	
	if __name__ == "__main__":
		
		with ProcessPoolExecutor(max_workers=int(MAX_PROC)) as executor:
			futures = [executor.submit(write_histo, *args) for args in param_first]
	        
			for future in as_completed(futures):
				result = future.result()
				NTASKS-=1
				print(f'{NTASKS0-NTASKS}/{NTASKS0}',result)
		 
		
		
		with ProcessPoolExecutor(max_workers=int(MAX_PROC)) as executor:
			futures = [executor.submit(write_histo, *args) for args in param_rest]
	        
			for future in as_completed(futures):
				result = future.result()
				NTASKS-=1
				print(f'{NTASKS0-NTASKS}/{NTASKS0}',result)
			
		
	print('')	
	print('***********************************')	
	print('')
	print(f'Siglevel {siglevel} %. Histograms saved in directories: ')
	for i in list_of_species:
		print(f'    {histodir}/{i}')
		
	print('')	
	print('***********************************')	
	print('***********************************')	
	print('***********************************')	
	print('')	
		
		
	
