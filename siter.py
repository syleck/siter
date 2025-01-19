import sys
import os
import glob
import shutil
import numpy as np
import scipy as sp
from scipy.spatial import distance_matrix
from scipy.spatial import distance
from itertools import combinations

from Bio.PDB import *
from Bio.PDB import PDBParser, PDBIO
from rmsd import brute_permutation
from subprocess import check_output, PIPE, STDOUT
from multiprocessing import Pool,TimeoutError

class AtomSelect(Select):
    def accept_atom(self, atom):
        if ((atom.get_parent().get_resname() in ["HIS", "LYS"] and atom.get_name() in ["ND1", "NE2", "NZ"]) or (atom.get_name() in ["O", "OG", "OG1", "OE1", "OE2", "OD1", "OD2", "OH"])) and (atom.altloc == "A" or atom.altloc ==' '):
            return True
        else:
            return (False)

def check (B): 
    Z=np.sort(np.around(np.nanmean(distance_matrix(B, B), axis=0), 1))
    if np.around(Z[-1]-Z[-4], 1)<=0.5 and np.around(Z[-5]-Z[-8], 1)<=0.5 and Z[-1]-Z[0]<1.0:
        condition=True
    else:
        condition=False
    return(condition)

def RMS (A, B):
    import rmsd
    A -= rmsd.centroid(A)
    B -= rmsd.centroid(B)
    C=brute_permutation(B, A)
    A=A[C]
    #print (A)
    U = rmsd.kabsch(A, B)
    A = np.dot(A, U)
    RMSD=rmsd.rmsd(A, B)
    return (RMSD, A, B)

wrdir,WAT,nc=sys.argv[1:4]

p = PDBParser()
io = PDBIO()
template=p.get_structure(WAT[0:-4], WAT)

###getting water coordinates###
water_coordinates=np.array([])
print (template[0])
for model in template:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name()=='OW':
                    print (atom)
                    water_coordinates=np.append(water_coordinates, atom.get_coord())
water_coordinates=water_coordinates.reshape(int(water_coordinates.shape[0]/3), 3)
print (water_coordinates)
###################################

###getting potassium coordinates###
K_coordinates=np.array([])
for model in template:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name()=='K':
                    print (atom)
                    K_coordinates=np.append(K_coordinates, atom.get_coord())
K_coordinates=K_coordinates.reshape(int(K_coordinates.shape[0]/3), 3)
print (K_coordinates)
#####################################

###builiding water distance matrix###
wat_dist_matrix=distance_matrix(water_coordinates, water_coordinates)
E=np.sort (wat_dist_matrix[0, 1:])
print (wat_dist_matrix)
print (E[-1]-E[0])
LNG=E.shape[0]+1

mx=np.nanmax(E)+2.0
mn=np.nanmin(E)-1.0

def aligner(mol):
    import rmsd
        
    print (mol[2:])
    print (mol[2:-4])
    io = PDBIO()
    structure=p.get_structure(mol[2:-4], mol[2:])
    #structure=structure[0]
    io.set_structure(structure)
    io.save(mol[2:-4]+'.ref', AtomSelect(), preserve_atom_numbering = True)
    prot_coordinates=np.array([])
    model=structure[0]
    print (model)
    for chain in model:
        for residue in chain:
            tags = residue.get_full_id()
            if tags[3][0]==" ":
                for atom in residue:
                    if ((atom.get_parent().get_resname() in ["HIS", "LYS"] and atom.get_name() in ["ND1", "NE2", "NZ"]) or (atom.get_name() in ["O", "OG", "OG1", "OE1", "OE2", "OD1", "OD2", "OH"])) and (atom.altloc == "A" or atom.altloc ==' '):
                        print (atom.get_coord)
                        prot_coordinates=np.append(prot_coordinates, atom.get_coord())
    prot_coordinates=prot_coordinates.reshape(int(prot_coordinates.shape[0]/3), 3)
    print (prot_coordinates)
    ALLat=np.array([])
    for chain in model:
        for residue in chain:
            tags = residue.get_full_id()
            if tags[3][0]==" ":
                for atom in residue:
                    print (atom.get_coord)
                    print ("Kostik is a gay")
                    ALLat=np.append(ALLat, atom.get_coord())
                

    ALLat=ALLat.reshape(int(ALLat.shape[0]/3), 3)
    print (ALLat)
    prot_dist_matrix=distance_matrix(prot_coordinates, prot_coordinates)
    DICT={} 
    COMBS=np.array([])
    for row in range(prot_dist_matrix.shape[0]):
        condition=np.logical_and(prot_dist_matrix[row]>=mn, prot_dist_matrix[row]<=mx)
        DICT[row]=np.where(condition==True)[0]
    info=open('atoms.txt', 'w')
    for key in DICT:
        VALS=np.array(DICT[key])
        VALS=np.sort(VALS)
        VARIANTS=np.array (list(combinations(VALS, water_coordinates.shape[0]-1)))
        info.write ("The number of variants for atom "+str(key+1)+" is " +str(len(VARIANTS))+'\n')
        print ("The number of variants fo atom "+str(key+1)+" is " +str(len(VARIANTS)))
        for variant in VARIANTS:
            var=np.concatenate(([key], variant), axis=0)
            var=np.sort(var)
            B=prot_coordinates[var]
            condition=check(B)
            if condition==True:
                print (var)
                COMBS=np.append(COMBS, var)
    info.close()
    COMBS=COMBS.reshape(int(COMBS.shape[0]/water_coordinates.shape[0]), water_coordinates.shape[0])
    COMBS=np.unique(COMBS, axis=0)
    COMBS=COMBS.astype(int)
    if COMBS.shape[0]!=0:
        print ((str(COMBS.shape[0])+' '+'suitable combinations are found!'))
    else:
        print ('No suitable combinations are found!')
    np.savetxt(mol[2:-4]+'_'+'COMBS.txt', COMBS, fmt='%.3f', delimiter=' ', newline='\n', header='')
    ATOMS=[]
    f=open(mol[2:-4]+'.ref', 'r')
    for line in f.readlines():
        line=line.strip()
        line=line.split(' ')
        if ("TER" not in line[0]) and ("END" not in line[0]) and ("HETATM" not in line[0]):
            print (line[0])
            ATOMS.append(line)
    x=1
    RES=[]
    for comb in COMBS:
            A=water_coordinates
            B=prot_coordinates[comb]
            B0=prot_coordinates[comb]
            RMSD, An, Bn = RMS(A, B)
            DIF=Bn-An
            ALG=B0-DIF
            K_new=np.round(ALG.mean(axis=0), 3)
            print (RMSD)
            x+=1
            print (x)
            if RMSD<=1.25:
                cond=True
                for at in ALLat:
                    dist=distance.euclidean (K_new, at)
                    if dist<2.1:
                        print ('The distance is too small: '+str(dist))
                        cond=False
                        print (cond)
                if cond==True:
                    print (" ")
                    print (cond)
                    print (' ')
                    print (' ')
                    print (comb)
                    print ("ШТА?")
                    print (RMSD)
                    b=np.around(np.array([RMSD]), 3)
                    print (comb)
                    k=np.concatenate((comb, b), axis=0)
                    print (k)
                    RES.append(k)
    if len (RES)>0:
        print (RES)
        print (' ')
        RES=np.array(RES)
        RES=RES[RES[:, 8].argsort()]
        print (" ")
        SITES=[]
        print (RES)
        print (" ")
        for i in range(RES.shape[0]):
            i=0
            alt=np.zeros([RES.shape[0]])
            X0=RES[i][0:-1]
            X0=X0.astype(int)
            X=prot_coordinates[X0]
            Xn=prot_coordinates[X0]
            W=water_coordinates
            rms1, X1, Y1 = RMS(W, X)
            A1=Y1-X1
            D1=Xn-A1
            KX=np.array([np.nanmean(D1, axis=0)])
            for j in range (RES.shape[0]):
                Y0=RES[j][0:-1] 
                Y0=Y0.astype(int)
                Y=prot_coordinates[Y0]
                Yn=prot_coordinates[Y0]
                rms2, X2, Y2 = RMS(W, Y)
                A2=Y2-X2
                D2=Yn- A2
                KY=np.array([np.nanmean(D2, axis=0)])
                print (distance.euclidean (KX.flatten(), KY.flatten()))
                dist=distance.euclidean (KX.flatten(), KY.flatten())
                alt[j]+=dist
            print (' ')
            print (alt)
            delete=np.where(alt<2.0)[0]
            SITES.append(RES[delete][0][0:-1])
            RES=np.delete(RES, delete, axis=0)
            print (' ')
            print (RES)
            print (' ')
            print ()
            print (RES.shape[0])
            if RES.shape[0]==0:
                break

        SITES=np.array(SITES, dtype=int)
        print (ATOMS)
        print (SITES)
        num=1
        f5=open(mol[2:-4]+"_RMSDs.log", 'w')
        f6=open(mol[2:-4]+"_RES.txt", 'w')
        for site in SITES:
            A=water_coordinates
            B=prot_coordinates[site]
            B0=prot_coordinates[site]
            RMSD, An, Bn = RMS(A, B)
            DIF=Bn-An
            ALG=B0-DIF
            K_new=np.round(ALG.mean(axis=0), 3)
            f2=open(mol[2:-4]+'_site_'+str(+num)+'.pdb', 'w')
            ALG=np.append(ALG, K_new)
            ALG=ALG.reshape(water_coordinates.shape[0]+1, 3)
            an=0
            f5.write("SITE_"+str(num)+": "+str(np.around(RMSD, 3))+'\n')
            print ("SITE_"+str(num)+": "+str(np.around(RMSD, 3))+'\n')
            f6.write(str(site)+' '+str(np.around(RMSD, 3))+'\n')
            print (str(site)+' '+str(np.around(RMSD, 3))+'\n')
            for atom in site:
                print (atom, site)
                print (ATOMS[atom])
                f2.write((' '.join(ATOMS[atom])+'\n'))
            f2.close()
            for model in template:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atom.set_coord(ALG[an])
                            an+=1
            print (' ')
            io.set_structure(template)
            io.save(mol[2:-4]+"_alignment_"+str(num)+".pdb", preserve_atom_numbering = True)
            print (' ')
            num+=1
        f5.close() 
        f6.close()
    files = glob.glob("./"+mol[2:-4]+"*")
    os.mkdir(mol[2:-4])
    for file in files:
        shutil.move(file, mol[2:-4])
    print (RES)

os.chdir(wrdir)
MOLS=glob.glob('./*pdb')
print (MOLS)
if __name__ == '__main__':
    ncpu=int(nc)
    if ncpu<0:
        ncpu=None
    with Pool(processes=ncpu) as pool:
        pool.map(aligner, MOLS)
        pool.close()
        pool.join()
