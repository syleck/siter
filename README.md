# How to run **siter.py**
The script siter.py finds binding sites for potassium ions in protein pdb structures by aligning a template of solvated K+ to the eights of protein atoms and calculates RMSD for these alignments. It is a pre-alpha version, so performance is very poor, probably you will also find some bugs. To run the script you will need the following python3 libraries that can be easily installed with pip:
+ numpy (**`pip install numpy`**)
+ scipy (**`pip install scipy`**)
+ Biopython (**`pip install biopython`**)
+ rmsd (**`pip install rmsd`**)
+ multiprocessing (**`pip install multiprocessing`**)

You will also need a pdb file that contains coordinates of potassium ion and water’s oxygens that coordinate the ion. Here we provide the template file named **KOW.pdb**. You will also need files of pdb structures you want to scan.

To run **siter.py** you should make the following steps:
1. Put **siter.py** and **KOW.pdb** in the same directory.
1. Put your protein pdb files to a separate directory which must be located in a folder where you put **siter.py** and **KOW.pdb earlier**.
1. Type the following command in your command line. Be sure that you use python3 by default, otherwise use python3 instead of python in the following command:\
   **`python siter.py <name_of_directories_with_PDBs> KOW.pdb <number_of_threads>`**\
  where **<name_of_directories_with_PDBs>** is the name of a folder where you put your pdb files with protein coordinate and **<number_of_threads>** is a number of threads of your CPU you want to use. 

After the calculations are completed go to the folder where you put your protein pdb files earlier. You will see that pdb files disappeared but instead of them, folders with the names of original pdbs were created. For example, your original dataset contained two pdb structures 1K4C.pdb and 2ITC.pdb, after the calculations you will see the corresponding folders named 1K4C and 2ITC. Go to any of these folders. You will see the following files:
* **<PDB_ID>.pdb** — the original pdb file.
* **<PDB_ID>.ref** — file that contains oxygens and nitrogens from original pdb that were used for scanning.
* **<PDB_ID>_COMBS.txt** — combinations of atoms that were used for calculations.
* **<PDB_ID>_alignment_X.pdb** — original template that was aligned to the protein atoms. X denotes a number of the alignment. 
* **<PDB_ID>_site_X.pdb** — this pdb file contains eight atoms that form the site for K+ and which were used for the corresponding alignment X.
* **<PDB_ID>_RES.txt** — the combinations of protein atoms that form the site are written in square brackets. The RMSD value for the alignment to this site is written to the right of them.
* **<PDB_ID>_RMSD.log** — this file contains RMSD values of the template alignment to the corresponding site. 

Now you can visualize pdb files in any molecular viewer like **PyMol** or **UCSF Chimera**. An example of the visualization is presented at **Fig. 1**. 
Note, that in some cases, when no suitable sites were found, you will not see **<PDB_ID>_alignment_X.pdb**, **<PDB_ID>_site_X.pdb**, **<PDB_ID>_RES.txt** and **<PDB_ID>_RMSD.log** files. 
![](https://github.com/syleck/siter/blob/main/example.png)

**Figure 1. An example of sites** (<PDB_ID>_site_X.pdb) **and alignments** (<PDB_ID>_alignment_X.pdb) **in 1K4C selectivity filter visualized by PyMol**. The alignments are represented by *big purple spheres* and *red crosses*, corresponding to K+ ions and water oxygens, respectively. *Small red spheres* represent protein oxygens that form the sites (in the case of nitrogen atoms the spheres will be deep blue). 

