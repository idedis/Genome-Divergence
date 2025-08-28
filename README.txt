Step by step instructions ONLY FOR LINUX:

1. Download genome_signature.zip and unpack it. This will create a directory genomedivergence with several files and two subdirectories ./seqs and ./histos
2. Optional: create a Python virtual environment (e.g. gs_env) and activate it by typing:
   $ mkdir gs_env
   $ python3 -m venv gs_env
   $ source gs_env/bin/activate
3. Install the Python required packages: psutil, scipy, numpy, pandas and matplotlib by typing:
   $ pip install psutil scipy numpy pandas matplotlib
4. You will need the utility parallel, if you need install it by typing:
   $ sudo apt install parallel    
5. Find the file input_file.txt and open it. There you can modifify the options to run the different scripts. In this file you specify the genomes you want to segment, the different
   significance levels (if you want to try several), the number of bins for the histograms (also you can do it for serveral values), specify the root directory for the DNA sequences,
   for the segmentation files and for the histograms. You also indicate the list of species being segmented, the maximum number of processes sent in parallel (if you don't know
   what is this, dont modify). Finaly specify the root directory for the files with distances among species.

   All the DNA files (contigs or chromosomes) of each species in fasta format (.fa or .fasta) should be in a subdirectory inside seqs. The name of the subdirectory should be Genus_specific epithet (e.g. Homo_sapiens)

   As an example we have included three "species": species_01,species_02 and species_03 (please check directory seqs)
 
6. Before knowing what we are doing, let's check that everthing works fine, type
   $ python segment_parall.py
   If everthing goes OK, you should obtain something like:

s=0.95 Alfabeto= -A 4 -bo seqs/species_01/ --> segments/species_01/
s=0.95 Alfabeto= -A 4 -bo seqs/species_02/ --> segments/species_02/
s=0.95 Alfabeto= -A 4 -bo seqs/species_03/ --> segments/species_03/
s=0.99 Alfabeto= -A 4 -bo seqs/species_01/ --> segments/species_01/
s=0.99 Alfabeto= -A 4 -bo seqs/species_02/ --> segments/species_02/
s=0.99 Alfabeto= -A 4 -bo seqs/species_03/ --> segments/species_03/
******************
18 segmentations launched in background
You can also monitor the processes with top or htop
Processes running 0 ... (Ctr+C to stop)   

All segmentations finished

    You can try with other sequences, note that a complete genome could take a lot of time. If your computer has many cores you can increase the number of parallel processes from 16 to the desired value. 

7. Now you have the segmentations done. You can explore the directory segments to see the file with the segmentations one for each contig or chromosome.
   They have a header and, for each segment a line with begining of the segment, end of the segment, size, four columns with the number of each nucleotide (ATCG), two columns with composition (R/Y), two
   columns with composition (S/W) and other two with composition (K/M)

8. To create the histograms of R/Y, S/W and K/M composition run:
   $ python get_histos_parall.py
   You should obtain something like:
   
Deleting old segments/species_01/total_95.dat
Deleting old segments/species_02/total_95.dat
Creating segments/species_01/total_95.dat
Creating segments/species_02/total_95.dat
Deleting old segments/species_03/total_95.dat
Creating segments/species_03/total_95.dat
1/6 species_01	bins= 25	 
2/6 species_02	bins= 25	 
3/6 species_03	bins= 25	 
4/6 species_01	bins= 50	 
5/6 species_02	bins= 50	 
6/6 species_03	bins= 50	 

***********************************

Siglevel 95 %. Histograms saved in directories: 
    histos/species_01
    histos/species_02
    histos/species_03

***********************************
***********************************
***********************************
9. To plot the histograms you can type:
   $ python fig_example 0.95 025 species_01 species_02 species_03
  you should obtain something like figure screen_capture.png (in your folder)
   If you type $ python fig_example02 0.95 025 species_01 species_02 species_03
   you plot only RY histograms (figure screen_capture02.png)
    
  
11.
12. (reproduce figures 1 and 2 of "Genome Divergence Based on Entropic Segmentation of DNA" , Bernaola-Galván et al.) run:
$ python 

 
   
   
