# written by Jan H. Jensen
#
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import sys

#from rdkit import rdBase
#print rdBase.rdkitVersion

def TS_conf_search(mol,template,num_confs,name):
    confs = []
    GetFF = lambda x,confId=-1:AllChem.MMFFGetMoleculeForceField(x,AllChem.MMFFGetMoleculeProperties(x),confId=confId)
    for i in range(num_confs):
#       mol = AllChem.ConstrainedEmbed(mol,template,useTethers=True,randomseed=-1,getForceField=GetFF)
        mol = ConstrainedEmbed(mol,template,useTethers=True,randomseed=-1,getForceField=GetFF)
        confs.append(Chem.Mol(mol,False,0))

    w = Chem.SDWriter(name+'=confs.sdf') 
    for i,conf in enumerate(confs):
# write to one file. Convenient for overlays in Pymol
        w.write(conf)
# write to individual files: convenient for QM calculations
        Chem.SDWriter(name+'='+str(i)+'.sdf').write(conf)
        
    
    return

def ConstrainedEmbed(mol, core, useTethers=True, coreConfId=-1, randomseed=2342,getForceField=AllChem.UFFGetMoleculeForceField, **kwargs):
# 
#  Copyright (C) 2006-2017  greg Landrum and Rational Discovery LLC 
# 
#   @@ All Rights Reserved @@ 
#  This file is part of the RDKit. 
#  The contents are covered by the terms of the BSD license 
#  which is included in the file license.txt, found at the root 
#  of the RDKit source tree. 
# 
    force_constant = 1000.
    match = mol.GetSubstructMatch(core)
    if not match:
        raise ValueError("molecule doesn't match the core")
    coordMap = {}
    coreConf = core.GetConformer(coreConfId)
    for i, idxI in enumerate(match):
        corePtI = coreConf.GetAtomPosition(i)
        coordMap[idxI] = corePtI

    if "." in Chem.MolToSmiles(mol):
        ci = AllChem.EmbedMolecule(mol, randomSeed=randomseed, **kwargs)  #jhj
    else:
        ci = AllChem.EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, **kwargs)

    if ci < 0:
        raise ValueError('Could not embed molecule.')

    algMap = [(j, i) for i, j in enumerate(match)]

    if not useTethers:
        # clean up the conformation
        ff = getForceField(mol, confId=0)
        for i, idxI in enumerate(match):
            for j in range(i + 1, len(match)):
                idxJ = match[j]
                d = coordMap[idxI].Distance(coordMap[idxJ])
                ff.AddDistanceConstraint(idxI, idxJ, d, d, force_constant)
        ff.Initialize()
        n = 4
        more = ff.Minimize()
        while more and n:
            more = ff.Minimize()
            n -= 1
        # rotate the embedded conformation onto the core:
        rms = rdMolAlign.AlignMol(mol, core, atomMap=algMap)
    else:
        # rotate the embedded conformation onto the core:
        rms = rdMolAlign.AlignMol(mol, core, atomMap=algMap)
        ff = getForceField(mol, confId=0)
        conf = core.GetConformer()
        for i in range(core.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            q = mol.GetConformer().GetAtomPosition(i)
            pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
            ff.AddDistanceConstraint(pIdx, match[i], 0, 0, force_constant)
        ff.Initialize()
        n = 4
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        while more and n:
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            n -= 1
        # realign
        rms = rdMolAlign.AlignMol(mol, core, atomMap=algMap)
    mol.SetProp('EmbedRMS', str(rms))
    return mol

if __name__ == '__main__':

    mol_filename = sys.argv[1]
    template_filename = sys.argv[2]
    if len(sys.argv) == 4:
       num_confs = int(sys.argv[3])
    else:
       num_confs = 20
    
    if ".sdf" in mol_filename:
        mol = Chem.MolFromMolFile(mol_filename,removeHs=False)
        name = mol_filename.split(".")[0]
    else:
        name, smiles = mol_filename.split(";")
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
    
    if "keepHs" in template_filename:
        template = Chem.MolFromMolFile(template_filename,removeHs=False)
    else:
        template = Chem.MolFromMolFile(template_filename,removeHs=True)
    
    TS_conf_search(mol,template,num_confs,name)
