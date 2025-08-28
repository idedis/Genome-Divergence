Step-by-step instructions (Linux only):

1. Download `genome_signature.zip` and unpack it. This will create a directory named `genomedivergence` containing several files and two subdirectories: `./seqs` and `./histos`.

2. (Optional) Create a Python virtual environment (e.g., `gs_env`) and activate it:

   $ python3 -m venv gs_env
   $ source gs_env/bin/activate

3. Install the required Python packages (`psutil`, `scipy`, `numpy`, `pandas`, and `matplotlib`):
   
   $ pip install psutil scipy numpy pandas matplotlib

4. You will need the utility **parallel**. Install it if necessary:
   
   $ sudo apt install parallel

5. Locate the file `input_file.txt` and open it. Here you can modify the options for running the scripts. In this file, you specify:
   - The significance levels (you may specify multiple values).
   - The number of bins for the histograms (again, multiple values are possible).
   - The root directories for the DNA sequences, the segmentation files, and the histograms.
   - The list of species being segmented.
   - The maximum number of processes to run in parallel (if unsure, do not change this).
   - The root directory for the distance files among species.

   All DNA files (contigs or chromosomes) in FASTA format (`.fa` or `.fasta`) for each species must be placed in a subdirectory inside `seqs`. The subdirectory name should follow the format Genus_species (e.g., Homo_sapiens).

   As an example, three species are included: `species_01`, `species_02`, and `species_03` (see the `seqs` directory).

6. Before proceeding further, check that everything works correctly:
   
   $ python segment_parall.py
   
   If everything is fine, you should see output similar to:
   s=0.95 Alfabeto= -A 4 -bo seqs/species_01/ --> segments/species_01/
   s=0.95 Alfabeto= -A 4 -bo seqs/species_02/ --> segments/species_02/
   s=0.95 Alfabeto= -A 4 -bo seqs/species_03/ --> segments/species_03/
   ...
   18 segmentations launched in background

   You can monitor the processes with `top` or `htop`.  
   When all segmentations are complete, you will see:
   All segmentations finished

   You can try other sequences, but note that a complete genome may take a long time. If your computer has many cores, you can increase the number of parallel processes (default: 16).

7. Once the segmentations are finished, you can explore the `segments` directory. Each file corresponds to a contig or chromosome and includes:
   - A header.
   - For each segment: start position, end position, size, nucleotide counts (A, T, C, G), composition (R/Y), composition (S/W), and composition (K/M).

8. To create histograms of R/Y, S/W, and K/M composition, run:

   $ python get_histos_parall.py
   
   This will generate output similar to:
   
   Deleting old segments/species_01/total_95.dat
   Creating segments/species_01/total_95.dat
   ...
   Siglevel 95 %. Histograms saved in:
       histos/species_01
       histos/species_02
       histos/species_03

9. To plot the histograms run:

   $ python fig_example.py 0.95 025 species_01 species_02 species_03

   This produces an output similar to `screen_capture.png`.

   Alternatively, to plot only R/Y histograms:

   $ python fig_example02.py 0.95 025 species_01 species_02 species_03
   
   This generates `screen_capture02.png`.

10. The `histos` directory already contains all histograms used in the paper "Genome Divergence Based on Entropic Segmentation of DNA" (Bernaola-Galván et al.).

   To reproduce Figure 1 of the paper run:
   $ python fig_example.py 95 50 Homo_sapiens Xenopus_tropicalis Alligator_mississippiensis Gallus_gallus Danio_rerio Drosophila_melanogaster Asterias_rubens Oryza_sativa Saccharomyces_cerevisiae
   
   This creates `fig_example_various.pdf`.

   To reproduce Figure 2 run:
   
   $ python fig_example02.py 95 50 Homo_sapiens Gorilla_gorilla Pan_troglodytes Felis_catus Canis_lupus Mustela_putorius Rattus_norvegicus Mus_musculus Cricetulus_griseus

11. To compute the distance matrix (i.e., the Jensen-Shannon distances between species), run:

   $ python get_distances.py

   You should see:
   Siglevel 0.95 ...done
   Siglevel 0.99 ...done

   In the `distances` directory you will find:
   - matrix_JS...tsv: square symmetric matrices of JS distances.
   - table_JS...tsv: tables listing all possible species pairs with their JS distances.

   Example: `table_JS_95_cg_025.tsv` is the distance table for segmentations at 95% significance, 25 bins, alphabet SW (cg).

12. Divergence times (from www.timetree.org) are available in the file `time_distance.dat`, which lists divergence times (in Myears) for all species pairs.
