''' reads initial data, it is common for all programs although not all
of them need all parameters'''


def initial_data(file='input_file.txt'):
	fini=open(file)
	info=[]
	for line in [x for x in fini if not(x.startswith('#'))]:
		info.append(line.replace('\n',''))
	sig_list=info[0].split(',')
	bin_list=info[1].split(',')
	dirfasta0=info[2] 
	segdir=info[3]
	histodir=info[4]
	list_of_species=info[5].split(',')
	MAX_PROC=info[6]
	dirdistances=info[7]
	
	
	fini.close()
	
	return sig_list,bin_list,dirfasta0,segdir,histodir,\
			list_of_species,MAX_PROC,dirdistances
	
