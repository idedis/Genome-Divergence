'''Hace un plot con 9 histogramas de la lista que hay abajo
se importa de los nonmammals de listados_especies.py , dependiendo de
la que se ponga hace gc,ry o km
La lista es una lista de directorios donde tiene que haber un fichero
llamado histo_XX.dat donde XX es CG,RY o KM. Este fichero tiene
2 columnas, la primera el fin del bin y la segunda el valor del
histograma de densidad normalizado.

Para poner el nbines y el nivel de significación se cambia aqui

snbins='200'
siglevel='95'
file_list=['/home/shared/rick/genomesignature/histos/'+x.replace('histo',f'histo{siglevel}').replace('050',snbins) for x in file_list]
'''



import numpy as np
import matplotlib.pyplot as plt
import sys

# Lista con los nombres de los 6 ficheros

# Reads parameters from file to have the default directories
from read_input_file import initial_data
siglevel_list,bin_list,dirfasta0,dirseg,dirhisto,list_of_species,\
  MAX_PROC,dirdistances = initial_data('input_file.txt')
# Reads parameters
#------------------------------------------------
siglevel=sys.argv[1].replace('0.','')
nbins=f'{int(sys.argv[2]):03d}'
if len(sys.argv)>3:
	list_of_species=[]
	for i in sys.argv[3:]:
		list_of_species.append(i)
		
		

file_list=[x+'/histo_cg050.dat' for x in list_of_species]






snbins=nbins

file_list=[f'{dirhisto}/{x}/histo{siglevel}_cg{nbins}.dat' for x in list_of_species]

letras=('a','b','c','d','e','f','g','h','i')

name_list=[]
for i,(let,fich) in enumerate(zip(letras,list_of_species)):
	name_list.append(f'({let}) '+fich.replace('_','\n     '))


#name_list[8]='Saccharomyces\ncerevisiae'

# Leer todos los datos y calcular el máximo de densidades
hist_data = []
hist_dataRY = []
hist_dataKM = []
y_max = 0

for filename in file_list:
	#esto para el gc
    data = np.loadtxt(filename)
    bin_ends = data[:, 0]
    densities = data[:, 1]
    
    # Insertar 0 al inicio de los bordes de los bins
    bin_edges = np.insert(bin_ends, 0, 0.0)
    bin_edges = bin_edges*100
    
    # Calcular anchos de bin y centros para graficar
    bin_widths = np.diff(bin_edges)
    bin_centers = bin_edges[:-1] + bin_widths / 2
    
    # Guardar para graficar
    hist_data.append((bin_centers, densities, bin_widths))
    y_max = max(y_max, np.max(densities))
    
    ######################################################
    #esto para el RY
    data = np.loadtxt(filename.replace('_cg','_ry'))
    bin_ends = data[:, 0]
    densities = data[:, 1]
    
    # Insertar 0 al inicio de los bordes de los bins
    bin_edges = np.insert(bin_ends, 0, 0.0)
    bin_edges = bin_edges*100
    
    # Calcular anchos de bin y centros para graficar
    bin_widths = np.diff(bin_edges)
    bin_centers = bin_edges[:-1] + bin_widths / 2
    
    # Guardar para graficar
    hist_dataRY.append((bin_centers, densities, bin_widths))
    y_max = max(y_max, np.max(densities))
    
    ######################################################
    #esto para el KM
    data = np.loadtxt(filename.replace('_cg','_km'))
    bin_ends = data[:, 0]
    densities = data[:, 1]
    
    # Insertar 0 al inicio de los bordes de los bins
    bin_edges = np.insert(bin_ends, 0, 0.0)
    bin_edges = bin_edges*100
    
    # Calcular anchos de bin y centros para graficar
    bin_widths = np.diff(bin_edges)
    bin_centers = bin_edges[:-1] + bin_widths / 2
    
    # Guardar para graficar
    hist_dataKM.append((bin_centers, densities, bin_widths))
    y_max = max(y_max, np.max(densities))
    
#------------------------------------------------------------    
    
    

# Crear figura
fig, axs = plt.subplots(3, 3, figsize=(13, 8),sharex=True, sharey=True)
fig.suptitle(f"Histograms for {siglevel}% significance level and {int(nbins)} bins", fontsize=16)

axs = axs.flatten()  # Para iterar fácilmente

#color azul claro #ADD8E6

MAXY=0
# Graficar cada histograma
nfigs=min(len(list_of_species),9)
for i, ax in enumerate(axs[0:nfigs]):
    ax.tick_params(axis='both', labelsize=11,direction='in') 
    #-------
    centers, densities, widths = hist_data[i]
    ax.bar(centers, densities, width=widths, edgecolor='gray', align='center',\
    color='#ACC563',linewidth=0.5,alpha=0.7,label='G+C (strong)',zorder=3)
    MAXY=max(max(densities),MAXY)
    ax.set_ylim(0, y_max * 1.1)  # Misma escala y con margen
    #------
    #centers, densities, widths = hist_dataRY[i]
    #ax.bar(centers, densities, width=widths, edgecolor='gray', align='center',\
    #color='#E9ABEC',linewidth=0.5,label='A+G (purines)',zorder=1)
    #MAXY=max(max(densities),MAXY)
    #ax.set_ylim(0, y_max * 1.1)  # Misma escala y con margen
    #-------
    #centers, densities, widths = hist_dataKM[i]
    #ax.bar(centers, densities, width=widths, edgecolor='gray', align='center',\
    #color='#ADD8E6',linewidth=0.5,alpha=0.7,label='G+T (keto)',zorder=2)
    #MAXY=max(max(densities),MAXY)
    #ax.set_ylim(0, y_max * 1.1)  # Misma escala y con margen
    
    
    ax.legend()
    #ax.set_title(f'{name_list[i]}')
    ax.text(0.05, 0.95, f'{name_list[i]}',
        transform=ax.transAxes,
        ha='left', va='top',
        fontsize=14) , #bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    if i in (7,):
      ax.set_xlabel('Segment composition (%)',fontsize=18,labelpad=10)
    if i in(3,):
      ax.set_ylabel('Normalized histogram',fontsize=18,labelpad=10)

#fijo manualmente la escala Y
# basta con hacerlo en uno porque tenemos sharey=True
axs[0].set_xlim(0,100)
axs[0].set_ylim(0,MAXY*1.1) #(0,8.3)

plt.tight_layout()
plt.subplots_adjust(wspace=0.08, hspace=0.07,left=0.08,bottom=0.08,right=0.98,top=0.92)
plt.savefig('fig_example_various.pdf')
#plt.savefig('/home/rick/trabajo/articulo_genome_signature/version01/fig_example_various.pdf')
plt.show()
