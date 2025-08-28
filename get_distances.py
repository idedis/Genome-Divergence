''' Computes JS distance between the histograms of all posible couples
of species listed in input_file.txt
For each value of significance level, alphabet (cg,ry or km) and number of bins
create a tab separated file with a matrix distance among all species

The files are placed at ./distances and they are named as follows:
matrix_JS_{siglevel}_{alphabet}{nbins}.tsv

e.g. matrix_JS_95_ry_050.tsv

Results are also written to another tsv file with a list of all 
possible couples and the value of JS

e.g. table_JS_95_ry_050.tsv
'''


import numpy as np
import pandas as pd
import os
from scipy.spatial.distance import jensenshannon


#---------------------------------------------------------------------------------------
def calc_js(org1,org2):
	#fname1=f'{root}/{org1}/{histo}'
	#fname2=f'{root}/{org2}/{histo}'
	fname1=org1
	fname2=org2
	p = np.genfromtxt(fname1, comments="#")[:,1]
	q = np.genfromtxt(fname2, comments="#")[:,1]
	p=p/p.sum()
	q=q/q.sum()
	return  jensenshannon(p, q, base=2)  #**2  # base=2 para que JS ∈ [0, 1]
	# sin el cuadrado es la distancia de jensen shannon (la raiz de la divergencia)
#---------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
def calc_matrix(lista,siglevel,alf,nbins,label='todos'):   #label='mammals' 'rodents' 'todos'
	#en principio la idea es ponerlos todos y luego que el programa extraiga las parejas que quiera
	fich_histo=f'histo{siglevel}_{alf}{nbins}.dat'  #solo el nombre del fichero que es común a todos, luego
													#se añade delante el directorio con la especie
	#-----
	#---- writes in matrix form and in table form
	fout=dir_matrix+f'matrix_{label}_{siglevel}_{alf}_{nbins}.tsv'
	fout_table=dir_matrix+f'table_{label}_{siglevel}_{alf}_{nbins}.tsv'
	#----
	f=open(fout,'w')
	ft=open(fout_table,'w')
	cad=''
	for org in lista:
		cad=cad+'\t'+org
	f.write(cad+'\n')	
	x=np.zeros((len(lista),len(lista)))
	for i,org1 in enumerate(lista):
		fich1=dir_histo+org1+'/'+fich_histo
		cad=org1
		for j,org2 in enumerate(lista):
			fich2=dir_histo+org2+'/'+fich_histo
			if i<j: 
				x[i,j]=calc_js(fich1,fich2)
				ft.write(f'{org1}\t{org2}\t{x[i,j]:.8f}\n')
				ft.write(f'{org2}\t{org1}\t{x[i,j]:.8f}\n')
			elif i>j:
				x[i,j]=x[j,i]
			cad=cad+f'\t{x[i,j]:.7f}'
		f.write(cad+'\n')
	f.close()		
#-------------------------------------------------------------------------------------------------------			
	



# Reads parameters from file
fini=open('input_file.txt')
info=[]
for line in [x for x in fini if not(x.startswith('#'))]:
	info.append(line.replace('\n',''))
sig_list=info[0].split(',')
bin_list=info[1].split(',') 
dirfasta0=info[2] #not used here
segdir=info[3] #not used here
dir_histo=info[4]+'/' 
list_of_species=info[5].split(',')
MAX_PROC=info[6] #not used here
dir_matrix=info[7]+'/'
fini.close()

alph_list=('cg','ry','km')
#--------------------------------------------------------------

sig_list=[x.split('.')[1] for x in sig_list] 


os.system(f'mkdir -p {dir_matrix}')
for siglevel in sig_list:
	print(f'Siglevel 0.{siglevel} ...',end='')
	for nbins in bin_list:
		for alph in alph_list:
			#print(f'{nbins}\t{siglevel}\t{alph}')
			calc_matrix(list_of_species,siglevel,alph,nbins,label='JS')
	print('done')






