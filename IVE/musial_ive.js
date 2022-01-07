/*
[1] Definition of biological/bioinformatics related data and variables.
*/
// Matrix of pairwise residue interaction scores respecting their observed frequency in a set of 100 membrane proteins, each with a pairwise identitiy of at most 15% to each other.
var PRISMEM15_SCORES = {
    "CYS-Helix/CYS-Helix": 6.0, "CYS-Helix/ASP-Helix": 2.0, "CYS-Helix/SER-Helix": 2.0, "CYS-Helix/GLN-Helix": 2.0, "CYS-Helix/LYS-Helix": 1.0, "CYS-Helix/ILE-Helix": 4.0, "CYS-Helix/PRO-Helix": 4.0, "CYS-Helix/THR-Helix": 4.0, "CYS-Helix/PHE-Helix": 4.0, "CYS-Helix/ASN-Helix": 2.0, "CYS-Helix/GLY-Helix": 3.0, "CYS-Helix/HIS-Helix": 2.0, "CYS-Helix/LEU-Helix": 4.0, "CYS-Helix/ARG-Helix": 0.0, "CYS-Helix/TRP-Helix": 1.0, "CYS-Helix/ALA-Helix": 4.0, "CYS-Helix/VAL-Helix": 4.0, "CYS-Helix/GLU-Helix": 0.0, "CYS-Helix/TYR-Helix": 3.0, "CYS-Helix/MET-Helix": 4.0, "CYS-Helix/CYS-Coil": 8.0, "CYS-Helix/ASP-Coil": -Infinity, "CYS-Helix/SER-Coil": 1.0, "CYS-Helix/GLN-Coil": -0.0, "CYS-Helix/LYS-Coil": -1.0, "CYS-Helix/ILE-Coil": 2.0, "CYS-Helix/PRO-Coil": 3.0, "CYS-Helix/THR-Coil": -Infinity, "CYS-Helix/PHE-Coil": -Infinity, "CYS-Helix/ASN-Coil": 0.0, "CYS-Helix/GLY-Coil": 3.0, "CYS-Helix/HIS-Coil": 4.0, "CYS-Helix/LEU-Coil": 3.0, "CYS-Helix/ARG-Coil": 1.0, "CYS-Helix/TRP-Coil": -Infinity, "CYS-Helix/ALA-Coil": 2.0, "CYS-Helix/VAL-Coil": 1.0, "CYS-Helix/GLU-Coil": 2.0, "CYS-Helix/TYR-Coil": 1.0, "CYS-Helix/MET-Coil": -Infinity, "CYS-Helix/CYS-Sheet": 5.0, "CYS-Helix/ASP-Sheet": -Infinity, "CYS-Helix/SER-Sheet": 1.0, "CYS-Helix/GLN-Sheet": -Infinity, "CYS-Helix/LYS-Sheet": 3.0, "CYS-Helix/ILE-Sheet": 2.0, "CYS-Helix/PRO-Sheet": -Infinity, "CYS-Helix/THR-Sheet": -Infinity, "CYS-Helix/PHE-Sheet": 1.0, "CYS-Helix/ASN-Sheet": 4.0, "CYS-Helix/GLY-Sheet": 2.0, "CYS-Helix/HIS-Sheet": -Infinity, "CYS-Helix/LEU-Sheet": 1.0, "CYS-Helix/ARG-Sheet": -Infinity, "CYS-Helix/TRP-Sheet": -Infinity, "CYS-Helix/ALA-Sheet": 1.0, "CYS-Helix/VAL-Sheet": 2.0, "CYS-Helix/GLU-Sheet": 1.0, "CYS-Helix/TYR-Sheet": 3.0, "CYS-Helix/MET-Sheet": 3.0, "ASP-Helix/CYS-Helix": 2.0, "ASP-Helix/ASP-Helix": 2.0, "ASP-Helix/SER-Helix": 4.0, "ASP-Helix/GLN-Helix": 4.0, "ASP-Helix/LYS-Helix": 5.0, "ASP-Helix/ILE-Helix": 1.0, "ASP-Helix/PRO-Helix": 3.0, "ASP-Helix/THR-Helix": 3.0, "ASP-Helix/PHE-Helix": 1.0, "ASP-Helix/ASN-Helix": 4.0, "ASP-Helix/GLY-Helix": 2.0, "ASP-Helix/HIS-Helix": 3.0, "ASP-Helix/LEU-Helix": 1.0, "ASP-Helix/ARG-Helix": 5.0, "ASP-Helix/TRP-Helix": 1.0, "ASP-Helix/ALA-Helix": 3.0, "ASP-Helix/VAL-Helix": 1.0, "ASP-Helix/GLU-Helix": 2.0, "ASP-Helix/TYR-Helix": 0.0, "ASP-Helix/MET-Helix": 0.0, "ASP-Helix/CYS-Coil": -1.0, "ASP-Helix/ASP-Coil": 0.0, "ASP-Helix/SER-Coil": 3.0, "ASP-Helix/GLN-Coil": 0.0, "ASP-Helix/LYS-Coil": 2.0, "ASP-Helix/ILE-Coil": 2.0, "ASP-Helix/PRO-Coil": 2.0, "ASP-Helix/THR-Coil": 2.0, "ASP-Helix/PHE-Coil": 2.0, "ASP-Helix/ASN-Coil": 2.0, "ASP-Helix/GLY-Coil": 1.0, "ASP-Helix/HIS-Coil": 2.0, "ASP-Helix/LEU-Coil": 2.0, "ASP-Helix/ARG-Coil": 3.0, "ASP-Helix/TRP-Coil": 1.0, "ASP-Helix/ALA-Coil": 1.0, "ASP-Helix/VAL-Coil": 1.0, "ASP-Helix/GLU-Coil": 1.0, "ASP-Helix/TYR-Coil": -0.0, "ASP-Helix/MET-Coil": -1.0, "ASP-Helix/CYS-Sheet": -1.0, "ASP-Helix/ASP-Sheet": -0.0, "ASP-Helix/SER-Sheet": 0.0, "ASP-Helix/GLN-Sheet": -1.0, "ASP-Helix/LYS-Sheet": 2.0, "ASP-Helix/ILE-Sheet": -2.0, "ASP-Helix/PRO-Sheet": 1.0, "ASP-Helix/THR-Sheet": -1.0, "ASP-Helix/PHE-Sheet": -2.0, "ASP-Helix/ASN-Sheet": -0.0, "ASP-Helix/GLY-Sheet": 0.0, "ASP-Helix/HIS-Sheet": 2.0, "ASP-Helix/LEU-Sheet": -1.0, "ASP-Helix/ARG-Sheet": 0.0, "ASP-Helix/TRP-Sheet": -Infinity, "ASP-Helix/ALA-Sheet": -3.0, "ASP-Helix/VAL-Sheet": -2.0, "ASP-Helix/GLU-Sheet": -Infinity, "ASP-Helix/TYR-Sheet": -3.0, "ASP-Helix/MET-Sheet": -Infinity, "SER-Helix/CYS-Helix": 2.0, "SER-Helix/ASP-Helix": 4.0, "SER-Helix/SER-Helix": 3.0, "SER-Helix/GLN-Helix": 3.0, "SER-Helix/LYS-Helix": 3.0, "SER-Helix/ILE-Helix": 2.0, "SER-Helix/PRO-Helix": 4.0, "SER-Helix/THR-Helix": 4.0, "SER-Helix/PHE-Helix": 3.0, "SER-Helix/ASN-Helix": 3.0, "SER-Helix/GLY-Helix": 4.0, "SER-Helix/HIS-Helix": 3.0, "SER-Helix/LEU-Helix": 3.0, "SER-Helix/ARG-Helix": 3.0, "SER-Helix/TRP-Helix": 1.0, "SER-Helix/ALA-Helix": 3.0, "SER-Helix/VAL-Helix": 3.0, "SER-Helix/GLU-Helix": 3.0, "SER-Helix/TYR-Helix": 2.0, "SER-Helix/MET-Helix": 3.0, "SER-Helix/CYS-Coil": 0.0, "SER-Helix/ASP-Coil": 2.0, "SER-Helix/SER-Coil": 1.0, "SER-Helix/GLN-Coil": 1.0, "SER-Helix/LYS-Coil": 1.0, "SER-Helix/ILE-Coil": 3.0, "SER-Helix/PRO-Coil": 2.0, "SER-Helix/THR-Coil": 2.0, "SER-Helix/PHE-Coil": 1.0, "SER-Helix/ASN-Coil": 1.0, "SER-Helix/GLY-Coil": 2.0, "SER-Helix/HIS-Coil": 1.0, "SER-Helix/LEU-Coil": 2.0, "SER-Helix/ARG-Coil": 0.0, "SER-Helix/TRP-Coil": 1.0, "SER-Helix/ALA-Coil": 2.0, "SER-Helix/VAL-Coil": 1.0, "SER-Helix/GLU-Coil": 2.0, "SER-Helix/TYR-Coil": 2.0, "SER-Helix/MET-Coil": 3.0, "SER-Helix/CYS-Sheet": -1.0, "SER-Helix/ASP-Sheet": 0.0, "SER-Helix/SER-Sheet": -2.0, "SER-Helix/GLN-Sheet": -3.0, "SER-Helix/LYS-Sheet": -1.0, "SER-Helix/ILE-Sheet": -0.0, "SER-Helix/PRO-Sheet": 3.0, "SER-Helix/THR-Sheet": 1.0, "SER-Helix/PHE-Sheet": -2.0, "SER-Helix/ASN-Sheet": -1.0, "SER-Helix/GLY-Sheet": -1.0, "SER-Helix/HIS-Sheet": -0.0, "SER-Helix/LEU-Sheet": -1.0, "SER-Helix/ARG-Sheet": 1.0, "SER-Helix/TRP-Sheet": -Infinity, "SER-Helix/ALA-Sheet": 1.0, "SER-Helix/VAL-Sheet": -1.0, "SER-Helix/GLU-Sheet": -Infinity, "SER-Helix/TYR-Sheet": 0.0, "SER-Helix/MET-Sheet": -Infinity, "GLN-Helix/CYS-Helix": 2.0, "GLN-Helix/ASP-Helix": 4.0, "GLN-Helix/SER-Helix": 3.0, "GLN-Helix/GLN-Helix": 2.0, "GLN-Helix/LYS-Helix": 3.0, "GLN-Helix/ILE-Helix": 3.0, "GLN-Helix/PRO-Helix": 3.0, "GLN-Helix/THR-Helix": 3.0, "GLN-Helix/PHE-Helix": 2.0, "GLN-Helix/ASN-Helix": 4.0, "GLN-Helix/GLY-Helix": 2.0, "GLN-Helix/HIS-Helix": 3.0, "GLN-Helix/LEU-Helix": 3.0, "GLN-Helix/ARG-Helix": 3.0, "GLN-Helix/TRP-Helix": 3.0, "GLN-Helix/ALA-Helix": 2.0, "GLN-Helix/VAL-Helix": 2.0, "GLN-Helix/GLU-Helix": 4.0, "GLN-Helix/TYR-Helix": 4.0, "GLN-Helix/MET-Helix": 2.0, "GLN-Helix/CYS-Coil": 3.0, "GLN-Helix/ASP-Coil": 1.0, "GLN-Helix/SER-Coil": 2.0, "GLN-Helix/GLN-Coil": 2.0, "GLN-Helix/LYS-Coil": 0.0, "GLN-Helix/ILE-Coil": 1.0, "GLN-Helix/PRO-Coil": 1.0, "GLN-Helix/THR-Coil": 1.0, "GLN-Helix/PHE-Coil": 1.0, "GLN-Helix/ASN-Coil": 1.0, "GLN-Helix/GLY-Coil": 1.0, "GLN-Helix/HIS-Coil": 2.0, "GLN-Helix/LEU-Coil": 3.0, "GLN-Helix/ARG-Coil": 1.0, "GLN-Helix/TRP-Coil": 3.0, "GLN-Helix/ALA-Coil": 2.0, "GLN-Helix/VAL-Coil": 0.0, "GLN-Helix/GLU-Coil": 0.0, "GLN-Helix/TYR-Coil": 2.0, "GLN-Helix/MET-Coil": 2.0, "GLN-Helix/CYS-Sheet": 1.0, "GLN-Helix/ASP-Sheet": -1.0, "GLN-Helix/SER-Sheet": -0.0, "GLN-Helix/GLN-Sheet": 2.0, "GLN-Helix/LYS-Sheet": -1.0, "GLN-Helix/ILE-Sheet": -0.0, "GLN-Helix/PRO-Sheet": 3.0, "GLN-Helix/THR-Sheet": -1.0, "GLN-Helix/PHE-Sheet": 0.0, "GLN-Helix/ASN-Sheet": -Infinity, "GLN-Helix/GLY-Sheet": -0.0, "GLN-Helix/HIS-Sheet": -Infinity, "GLN-Helix/LEU-Sheet": 0.0, "GLN-Helix/ARG-Sheet": 1.0, "GLN-Helix/TRP-Sheet": 0.0, "GLN-Helix/ALA-Sheet": 0.0, "GLN-Helix/VAL-Sheet": -0.0, "GLN-Helix/GLU-Sheet": -2.0, "GLN-Helix/TYR-Sheet": -1.0, "GLN-Helix/MET-Sheet": 1.0, "LYS-Helix/CYS-Helix": 1.0, "LYS-Helix/ASP-Helix": 5.0, "LYS-Helix/SER-Helix": 3.0, "LYS-Helix/GLN-Helix": 3.0, "LYS-Helix/LYS-Helix": -1.0, "LYS-Helix/ILE-Helix": 1.0, "LYS-Helix/PRO-Helix": 1.0, "LYS-Helix/THR-Helix": 2.0, "LYS-Helix/PHE-Helix": 2.0, "LYS-Helix/ASN-Helix": 4.0, "LYS-Helix/GLY-Helix": 1.0, "LYS-Helix/HIS-Helix": 4.0, "LYS-Helix/LEU-Helix": 2.0, "LYS-Helix/ARG-Helix": 1.0, "LYS-Helix/TRP-Helix": 4.0, "LYS-Helix/ALA-Helix": 1.0, "LYS-Helix/VAL-Helix": -0.0, "LYS-Helix/GLU-Helix": 5.0, "LYS-Helix/TYR-Helix": 4.0, "LYS-Helix/MET-Helix": 2.0, "LYS-Helix/CYS-Coil": -1.0, "LYS-Helix/ASP-Coil": 3.0, "LYS-Helix/SER-Coil": 1.0, "LYS-Helix/GLN-Coil": 3.0, "LYS-Helix/LYS-Coil": 0.0, "LYS-Helix/ILE-Coil": 1.0, "LYS-Helix/PRO-Coil": -0.0, "LYS-Helix/THR-Coil": 1.0, "LYS-Helix/PHE-Coil": 2.0, "LYS-Helix/ASN-Coil": 2.0, "LYS-Helix/GLY-Coil": 0.0, "LYS-Helix/HIS-Coil": -1.0, "LYS-Helix/LEU-Coil": 1.0, "LYS-Helix/ARG-Coil": -0.0, "LYS-Helix/TRP-Coil": -1.0, "LYS-Helix/ALA-Coil": 0.0, "LYS-Helix/VAL-Coil": 0.0, "LYS-Helix/GLU-Coil": 3.0, "LYS-Helix/TYR-Coil": 2.0, "LYS-Helix/MET-Coil": 2.0, "LYS-Helix/CYS-Sheet": -Infinity, "LYS-Helix/ASP-Sheet": 1.0, "LYS-Helix/SER-Sheet": -0.0, "LYS-Helix/GLN-Sheet": -2.0, "LYS-Helix/LYS-Sheet": -Infinity, "LYS-Helix/ILE-Sheet": -0.0, "LYS-Helix/PRO-Sheet": 2.0, "LYS-Helix/THR-Sheet": 0.0, "LYS-Helix/PHE-Sheet": -2.0, "LYS-Helix/ASN-Sheet": -Infinity, "LYS-Helix/GLY-Sheet": -Infinity, "LYS-Helix/HIS-Sheet": 1.0, "LYS-Helix/LEU-Sheet": -3.0, "LYS-Helix/ARG-Sheet": -2.0, "LYS-Helix/TRP-Sheet": 1.0, "LYS-Helix/ALA-Sheet": -1.0, "LYS-Helix/VAL-Sheet": -1.0, "LYS-Helix/GLU-Sheet": 1.0, "LYS-Helix/TYR-Sheet": -1.0, "LYS-Helix/MET-Sheet": -1.0, "ILE-Helix/CYS-Helix": 4.0, "ILE-Helix/ASP-Helix": 1.0, "ILE-Helix/SER-Helix": 2.0, "ILE-Helix/GLN-Helix": 3.0, "ILE-Helix/LYS-Helix": 1.0, "ILE-Helix/ILE-Helix": 2.0, "ILE-Helix/PRO-Helix": 1.0, "ILE-Helix/THR-Helix": 3.0, "ILE-Helix/PHE-Helix": 3.0, "ILE-Helix/ASN-Helix": 2.0, "ILE-Helix/GLY-Helix": 4.0, "ILE-Helix/HIS-Helix": 2.0, "ILE-Helix/LEU-Helix": 2.0, "ILE-Helix/ARG-Helix": 1.0, "ILE-Helix/TRP-Helix": 3.0, "ILE-Helix/ALA-Helix": 3.0, "ILE-Helix/VAL-Helix": 3.0, "ILE-Helix/GLU-Helix": 1.0, "ILE-Helix/TYR-Helix": 3.0, "ILE-Helix/MET-Helix": 3.0, "ILE-Helix/CYS-Coil": 3.0, "ILE-Helix/ASP-Coil": 1.0, "ILE-Helix/SER-Coil": 0.0, "ILE-Helix/GLN-Coil": 0.0, "ILE-Helix/LYS-Coil": -3.0, "ILE-Helix/ILE-Coil": 2.0, "ILE-Helix/PRO-Coil": 1.0, "ILE-Helix/THR-Coil": -2.0, "ILE-Helix/PHE-Coil": -1.0, "ILE-Helix/ASN-Coil": -1.0, "ILE-Helix/GLY-Coil": -1.0, "ILE-Helix/HIS-Coil": 0.0, "ILE-Helix/LEU-Coil": 1.0, "ILE-Helix/ARG-Coil": -0.0, "ILE-Helix/TRP-Coil": 2.0, "ILE-Helix/ALA-Coil": 1.0, "ILE-Helix/VAL-Coil": 1.0, "ILE-Helix/GLU-Coil": -1.0, "ILE-Helix/TYR-Coil": 1.0, "ILE-Helix/MET-Coil": 2.0, "ILE-Helix/CYS-Sheet": 0.0, "ILE-Helix/ASP-Sheet": 1.0, "ILE-Helix/SER-Sheet": -2.0, "ILE-Helix/GLN-Sheet": -2.0, "ILE-Helix/LYS-Sheet": -Infinity, "ILE-Helix/ILE-Sheet": 1.0, "ILE-Helix/PRO-Sheet": 1.0, "ILE-Helix/THR-Sheet": 0.0, "ILE-Helix/PHE-Sheet": 1.0, "ILE-Helix/ASN-Sheet": 0.0, "ILE-Helix/GLY-Sheet": 0.0, "ILE-Helix/HIS-Sheet": -Infinity, "ILE-Helix/LEU-Sheet": -0.0, "ILE-Helix/ARG-Sheet": -1.0, "ILE-Helix/TRP-Sheet": 1.0, "ILE-Helix/ALA-Sheet": 1.0, "ILE-Helix/VAL-Sheet": 1.0, "ILE-Helix/GLU-Sheet": -1.0, "ILE-Helix/TYR-Sheet": 2.0, "ILE-Helix/MET-Sheet": -0.0, "PRO-Helix/CYS-Helix": 4.0, "PRO-Helix/ASP-Helix": 3.0, "PRO-Helix/SER-Helix": 4.0, "PRO-Helix/GLN-Helix": 3.0, "PRO-Helix/LYS-Helix": 1.0, "PRO-Helix/ILE-Helix": 1.0, "PRO-Helix/PRO-Helix": 2.0, "PRO-Helix/THR-Helix": 4.0, "PRO-Helix/PHE-Helix": 3.0, "PRO-Helix/ASN-Helix": 3.0, "PRO-Helix/GLY-Helix": 4.0, "PRO-Helix/HIS-Helix": 3.0, "PRO-Helix/LEU-Helix": 2.0, "PRO-Helix/ARG-Helix": 2.0, "PRO-Helix/TRP-Helix": 4.0, "PRO-Helix/ALA-Helix": 3.0, "PRO-Helix/VAL-Helix": 3.0, "PRO-Helix/GLU-Helix": 3.0, "PRO-Helix/TYR-Helix": 3.0, "PRO-Helix/MET-Helix": 3.0, "PRO-Helix/CYS-Coil": 2.0, "PRO-Helix/ASP-Coil": 1.0, "PRO-Helix/SER-Coil": 2.0, "PRO-Helix/GLN-Coil": 1.0, "PRO-Helix/LYS-Coil": 1.0, "PRO-Helix/ILE-Coil": 2.0, "PRO-Helix/PRO-Coil": 1.0, "PRO-Helix/THR-Coil": 2.0, "PRO-Helix/PHE-Coil": 2.0, "PRO-Helix/ASN-Coil": -0.0, "PRO-Helix/GLY-Coil": 1.0, "PRO-Helix/HIS-Coil": 3.0, "PRO-Helix/LEU-Coil": 1.0, "PRO-Helix/ARG-Coil": 1.0, "PRO-Helix/TRP-Coil": 4.0, "PRO-Helix/ALA-Coil": 2.0, "PRO-Helix/VAL-Coil": 2.0, "PRO-Helix/GLU-Coil": -0.0, "PRO-Helix/TYR-Coil": 1.0, "PRO-Helix/MET-Coil": 3.0, "PRO-Helix/CYS-Sheet": 2.0, "PRO-Helix/ASP-Sheet": -1.0, "PRO-Helix/SER-Sheet": 1.0, "PRO-Helix/GLN-Sheet": -1.0, "PRO-Helix/LYS-Sheet": -2.0, "PRO-Helix/ILE-Sheet": 1.0, "PRO-Helix/PRO-Sheet": -Infinity, "PRO-Helix/THR-Sheet": -1.0, "PRO-Helix/PHE-Sheet": 1.0, "PRO-Helix/ASN-Sheet": 2.0, "PRO-Helix/GLY-Sheet": 1.0, "PRO-Helix/HIS-Sheet": -Infinity, "PRO-Helix/LEU-Sheet": -3.0, "PRO-Helix/ARG-Sheet": -Infinity, "PRO-Helix/TRP-Sheet": 2.0, "PRO-Helix/ALA-Sheet": 1.0, "PRO-Helix/VAL-Sheet": -0.0, "PRO-Helix/GLU-Sheet": -1.0, "PRO-Helix/TYR-Sheet": 2.0, "PRO-Helix/MET-Sheet": 1.0, "THR-Helix/CYS-Helix": 4.0, "THR-Helix/ASP-Helix": 3.0, "THR-Helix/SER-Helix": 4.0, "THR-Helix/GLN-Helix": 3.0, "THR-Helix/LYS-Helix": 2.0, "THR-Helix/ILE-Helix": 3.0, "THR-Helix/PRO-Helix": 4.0, "THR-Helix/THR-Helix": 1.0, "THR-Helix/PHE-Helix": 3.0, "THR-Helix/ASN-Helix": 3.0, "THR-Helix/GLY-Helix": 4.0, "THR-Helix/HIS-Helix": 3.0, "THR-Helix/LEU-Helix": 2.0, "THR-Helix/ARG-Helix": 2.0, "THR-Helix/TRP-Helix": 2.0, "THR-Helix/ALA-Helix": 3.0, "THR-Helix/VAL-Helix": 3.0, "THR-Helix/GLU-Helix": 3.0, "THR-Helix/TYR-Helix": 3.0, "THR-Helix/MET-Helix": 3.0, "THR-Helix/CYS-Coil": 1.0, "THR-Helix/ASP-Coil": 1.0, "THR-Helix/SER-Coil": -0.0, "THR-Helix/GLN-Coil": -1.0, "THR-Helix/LYS-Coil": 0.0, "THR-Helix/ILE-Coil": 2.0, "THR-Helix/PRO-Coil": 1.0, "THR-Helix/THR-Coil": 2.0, "THR-Helix/PHE-Coil": -0.0, "THR-Helix/ASN-Coil": 0.0, "THR-Helix/GLY-Coil": 2.0, "THR-Helix/HIS-Coil": 2.0, "THR-Helix/LEU-Coil": 1.0, "THR-Helix/ARG-Coil": 2.0, "THR-Helix/TRP-Coil": -2.0, "THR-Helix/ALA-Coil": 2.0, "THR-Helix/VAL-Coil": 2.0, "THR-Helix/GLU-Coil": 1.0, "THR-Helix/TYR-Coil": 1.0, "THR-Helix/MET-Coil": 3.0, "THR-Helix/CYS-Sheet": 2.0, "THR-Helix/ASP-Sheet": -1.0, "THR-Helix/SER-Sheet": -1.0, "THR-Helix/GLN-Sheet": -1.0, "THR-Helix/LYS-Sheet": -3.0, "THR-Helix/ILE-Sheet": -0.0, "THR-Helix/PRO-Sheet": 1.0, "THR-Helix/THR-Sheet": -1.0, "THR-Helix/PHE-Sheet": -Infinity, "THR-Helix/ASN-Sheet": 0.0, "THR-Helix/GLY-Sheet": 0.0, "THR-Helix/HIS-Sheet": -1.0, "THR-Helix/LEU-Sheet": -1.0, "THR-Helix/ARG-Sheet": -1.0, "THR-Helix/TRP-Sheet": -0.0, "THR-Helix/ALA-Sheet": 1.0, "THR-Helix/VAL-Sheet": 0.0, "THR-Helix/GLU-Sheet": -1.0, "THR-Helix/TYR-Sheet": -3.0, "THR-Helix/MET-Sheet": 1.0, "PHE-Helix/CYS-Helix": 4.0, "PHE-Helix/ASP-Helix": 1.0, "PHE-Helix/SER-Helix": 3.0, "PHE-Helix/GLN-Helix": 2.0, "PHE-Helix/LYS-Helix": 2.0, "PHE-Helix/ILE-Helix": 3.0, "PHE-Helix/PRO-Helix": 3.0, "PHE-Helix/THR-Helix": 3.0, "PHE-Helix/PHE-Helix": 2.0, "PHE-Helix/ASN-Helix": 3.0, "PHE-Helix/GLY-Helix": 3.0, "PHE-Helix/HIS-Helix": 3.0, "PHE-Helix/LEU-Helix": 4.0, "PHE-Helix/ARG-Helix": 3.0, "PHE-Helix/TRP-Helix": 4.0, "PHE-Helix/ALA-Helix": 3.0, "PHE-Helix/VAL-Helix": 3.0, "PHE-Helix/GLU-Helix": 1.0, "PHE-Helix/TYR-Helix": 3.0, "PHE-Helix/MET-Helix": 4.0, "PHE-Helix/CYS-Coil": 3.0, "PHE-Helix/ASP-Coil": -4.0, "PHE-Helix/SER-Coil": -1.0, "PHE-Helix/GLN-Coil": -3.0, "PHE-Helix/LYS-Coil": -2.0, "PHE-Helix/ILE-Coil": 2.0, "PHE-Helix/PRO-Coil": 1.0, "PHE-Helix/THR-Coil": -3.0, "PHE-Helix/PHE-Coil": -0.0, "PHE-Helix/ASN-Coil": -4.0, "PHE-Helix/GLY-Coil": -1.0, "PHE-Helix/HIS-Coil": -0.0, "PHE-Helix/LEU-Coil": 1.0, "PHE-Helix/ARG-Coil": -0.0, "PHE-Helix/TRP-Coil": -1.0, "PHE-Helix/ALA-Coil": 1.0, "PHE-Helix/VAL-Coil": 0.0, "PHE-Helix/GLU-Coil": -2.0, "PHE-Helix/TYR-Coil": -0.0, "PHE-Helix/MET-Coil": 3.0, "PHE-Helix/CYS-Sheet": 2.0, "PHE-Helix/ASP-Sheet": -Infinity, "PHE-Helix/SER-Sheet": -1.0, "PHE-Helix/GLN-Sheet": -Infinity, "PHE-Helix/LYS-Sheet": -1.0, "PHE-Helix/ILE-Sheet": -1.0, "PHE-Helix/PRO-Sheet": 1.0, "PHE-Helix/THR-Sheet": -3.0, "PHE-Helix/PHE-Sheet": 2.0, "PHE-Helix/ASN-Sheet": -Infinity, "PHE-Helix/GLY-Sheet": -2.0, "PHE-Helix/HIS-Sheet": -Infinity, "PHE-Helix/LEU-Sheet": 1.0, "PHE-Helix/ARG-Sheet": 2.0, "PHE-Helix/TRP-Sheet": -2.0, "PHE-Helix/ALA-Sheet": -2.0, "PHE-Helix/VAL-Sheet": -1.0, "PHE-Helix/GLU-Sheet": -Infinity, "PHE-Helix/TYR-Sheet": 0.0, "PHE-Helix/MET-Sheet": -Infinity, "ASN-Helix/CYS-Helix": 2.0, "ASN-Helix/ASP-Helix": 4.0, "ASN-Helix/SER-Helix": 3.0, "ASN-Helix/GLN-Helix": 4.0, "ASN-Helix/LYS-Helix": 4.0, "ASN-Helix/ILE-Helix": 2.0, "ASN-Helix/PRO-Helix": 3.0, "ASN-Helix/THR-Helix": 3.0, "ASN-Helix/PHE-Helix": 3.0, "ASN-Helix/ASN-Helix": 3.0, "ASN-Helix/GLY-Helix": 3.0, "ASN-Helix/HIS-Helix": 4.0, "ASN-Helix/LEU-Helix": 2.0, "ASN-Helix/ARG-Helix": 3.0, "ASN-Helix/TRP-Helix": 3.0, "ASN-Helix/ALA-Helix": 2.0, "ASN-Helix/VAL-Helix": 2.0, "ASN-Helix/GLU-Helix": 3.0, "ASN-Helix/TYR-Helix": 3.0, "ASN-Helix/MET-Helix": 3.0, "ASN-Helix/CYS-Coil": -1.0, "ASN-Helix/ASP-Coil": 2.0, "ASN-Helix/SER-Coil": 1.0, "ASN-Helix/GLN-Coil": 2.0, "ASN-Helix/LYS-Coil": 2.0, "ASN-Helix/ILE-Coil": 2.0, "ASN-Helix/PRO-Coil": 0.0, "ASN-Helix/THR-Coil": 1.0, "ASN-Helix/PHE-Coil": 1.0, "ASN-Helix/ASN-Coil": 2.0, "ASN-Helix/GLY-Coil": 2.0, "ASN-Helix/HIS-Coil": 0.0, "ASN-Helix/LEU-Coil": 1.0, "ASN-Helix/ARG-Coil": 1.0, "ASN-Helix/TRP-Coil": 1.0, "ASN-Helix/ALA-Coil": 2.0, "ASN-Helix/VAL-Coil": 2.0, "ASN-Helix/GLU-Coil": 3.0, "ASN-Helix/TYR-Coil": 2.0, "ASN-Helix/MET-Coil": 3.0, "ASN-Helix/CYS-Sheet": -Infinity, "ASN-Helix/ASP-Sheet": 0.0, "ASN-Helix/SER-Sheet": -1.0, "ASN-Helix/GLN-Sheet": -1.0, "ASN-Helix/LYS-Sheet": -Infinity, "ASN-Helix/ILE-Sheet": 0.0, "ASN-Helix/PRO-Sheet": 3.0, "ASN-Helix/THR-Sheet": -1.0, "ASN-Helix/PHE-Sheet": -Infinity, "ASN-Helix/ASN-Sheet": 1.0, "ASN-Helix/GLY-Sheet": -Infinity, "ASN-Helix/HIS-Sheet": 1.0, "ASN-Helix/LEU-Sheet": 0.0, "ASN-Helix/ARG-Sheet": 0.0, "ASN-Helix/TRP-Sheet": -1.0, "ASN-Helix/ALA-Sheet": -1.0, "ASN-Helix/VAL-Sheet": -2.0, "ASN-Helix/GLU-Sheet": -Infinity, "ASN-Helix/TYR-Sheet": -0.0, "ASN-Helix/MET-Sheet": -Infinity, "GLY-Helix/CYS-Helix": 3.0, "GLY-Helix/ASP-Helix": 2.0, "GLY-Helix/SER-Helix": 4.0, "GLY-Helix/GLN-Helix": 2.0, "GLY-Helix/LYS-Helix": 1.0, "GLY-Helix/ILE-Helix": 4.0, "GLY-Helix/PRO-Helix": 4.0, "GLY-Helix/THR-Helix": 4.0, "GLY-Helix/PHE-Helix": 3.0, "GLY-Helix/ASN-Helix": 3.0, "GLY-Helix/GLY-Helix": 3.0, "GLY-Helix/HIS-Helix": 3.0, "GLY-Helix/LEU-Helix": 3.0, "GLY-Helix/ARG-Helix": 2.0, "GLY-Helix/TRP-Helix": 3.0, "GLY-Helix/ALA-Helix": 5.0, "GLY-Helix/VAL-Helix": 4.0, "GLY-Helix/GLU-Helix": 1.0, "GLY-Helix/TYR-Helix": 2.0, "GLY-Helix/MET-Helix": 4.0, "GLY-Helix/CYS-Coil": 3.0, "GLY-Helix/ASP-Coil": 1.0, "GLY-Helix/SER-Coil": 2.0, "GLY-Helix/GLN-Coil": 1.0, "GLY-Helix/LYS-Coil": -1.0, "GLY-Helix/ILE-Coil": 2.0, "GLY-Helix/PRO-Coil": 1.0, "GLY-Helix/THR-Coil": 2.0, "GLY-Helix/PHE-Coil": 1.0, "GLY-Helix/ASN-Coil": 0.0, "GLY-Helix/GLY-Coil": 2.0, "GLY-Helix/HIS-Coil": 0.0, "GLY-Helix/LEU-Coil": 1.0, "GLY-Helix/ARG-Coil": -0.0, "GLY-Helix/TRP-Coil": -1.0, "GLY-Helix/ALA-Coil": 3.0, "GLY-Helix/VAL-Coil": 2.0, "GLY-Helix/GLU-Coil": 1.0, "GLY-Helix/TYR-Coil": -1.0, "GLY-Helix/MET-Coil": 1.0, "GLY-Helix/CYS-Sheet": -0.0, "GLY-Helix/ASP-Sheet": -1.0, "GLY-Helix/SER-Sheet": -2.0, "GLY-Helix/GLN-Sheet": -0.0, "GLY-Helix/LYS-Sheet": -1.0, "GLY-Helix/ILE-Sheet": -0.0, "GLY-Helix/PRO-Sheet": 0.0, "GLY-Helix/THR-Sheet": -0.0, "GLY-Helix/PHE-Sheet": 0.0, "GLY-Helix/ASN-Sheet": -1.0, "GLY-Helix/GLY-Sheet": -0.0, "GLY-Helix/HIS-Sheet": -Infinity, "GLY-Helix/LEU-Sheet": -2.0, "GLY-Helix/ARG-Sheet": 0.0, "GLY-Helix/TRP-Sheet": -1.0, "GLY-Helix/ALA-Sheet": -1.0, "GLY-Helix/VAL-Sheet": -1.0, "GLY-Helix/GLU-Sheet": -2.0, "GLY-Helix/TYR-Sheet": -1.0, "GLY-Helix/MET-Sheet": -Infinity, "HIS-Helix/CYS-Helix": 2.0, "HIS-Helix/ASP-Helix": 3.0, "HIS-Helix/SER-Helix": 3.0, "HIS-Helix/GLN-Helix": 3.0, "HIS-Helix/LYS-Helix": 4.0, "HIS-Helix/ILE-Helix": 2.0, "HIS-Helix/PRO-Helix": 3.0, "HIS-Helix/THR-Helix": 3.0, "HIS-Helix/PHE-Helix": 3.0, "HIS-Helix/ASN-Helix": 4.0, "HIS-Helix/GLY-Helix": 3.0, "HIS-Helix/HIS-Helix": 4.0, "HIS-Helix/LEU-Helix": 2.0, "HIS-Helix/ARG-Helix": 2.0, "HIS-Helix/TRP-Helix": 4.0, "HIS-Helix/ALA-Helix": 2.0, "HIS-Helix/VAL-Helix": 2.0, "HIS-Helix/GLU-Helix": 3.0, "HIS-Helix/TYR-Helix": 3.0, "HIS-Helix/MET-Helix": 3.0, "HIS-Helix/CYS-Coil": 4.0, "HIS-Helix/ASP-Coil": 2.0, "HIS-Helix/SER-Coil": 1.0, "HIS-Helix/GLN-Coil": 2.0, "HIS-Helix/LYS-Coil": 0.0, "HIS-Helix/ILE-Coil": 3.0, "HIS-Helix/PRO-Coil": 2.0, "HIS-Helix/THR-Coil": 2.0, "HIS-Helix/PHE-Coil": 3.0, "HIS-Helix/ASN-Coil": 2.0, "HIS-Helix/GLY-Coil": 1.0, "HIS-Helix/HIS-Coil": 2.0, "HIS-Helix/LEU-Coil": 2.0, "HIS-Helix/ARG-Coil": 2.0, "HIS-Helix/TRP-Coil": -Infinity, "HIS-Helix/ALA-Coil": 2.0, "HIS-Helix/VAL-Coil": 1.0, "HIS-Helix/GLU-Coil": 2.0, "HIS-Helix/TYR-Coil": 2.0, "HIS-Helix/MET-Coil": 2.0, "HIS-Helix/CYS-Sheet": -Infinity, "HIS-Helix/ASP-Sheet": -Infinity, "HIS-Helix/SER-Sheet": -2.0, "HIS-Helix/GLN-Sheet": -1.0, "HIS-Helix/LYS-Sheet": -1.0, "HIS-Helix/ILE-Sheet": -2.0, "HIS-Helix/PRO-Sheet": -Infinity, "HIS-Helix/THR-Sheet": -Infinity, "HIS-Helix/PHE-Sheet": 0.0, "HIS-Helix/ASN-Sheet": -Infinity, "HIS-Helix/GLY-Sheet": -Infinity, "HIS-Helix/HIS-Sheet": -Infinity, "HIS-Helix/LEU-Sheet": 1.0, "HIS-Helix/ARG-Sheet": -1.0, "HIS-Helix/TRP-Sheet": -Infinity, "HIS-Helix/ALA-Sheet": 0.0, "HIS-Helix/VAL-Sheet": 0.0, "HIS-Helix/GLU-Sheet": -Infinity, "HIS-Helix/TYR-Sheet": -Infinity, "HIS-Helix/MET-Sheet": -Infinity, "LEU-Helix/CYS-Helix": 4.0, "LEU-Helix/ASP-Helix": 1.0, "LEU-Helix/SER-Helix": 3.0, "LEU-Helix/GLN-Helix": 3.0, "LEU-Helix/LYS-Helix": 2.0, "LEU-Helix/ILE-Helix": 2.0, "LEU-Helix/PRO-Helix": 2.0, "LEU-Helix/THR-Helix": 2.0, "LEU-Helix/PHE-Helix": 4.0, "LEU-Helix/ASN-Helix": 2.0, "LEU-Helix/GLY-Helix": 3.0, "LEU-Helix/HIS-Helix": 2.0, "LEU-Helix/LEU-Helix": 2.0, "LEU-Helix/ARG-Helix": 2.0, "LEU-Helix/TRP-Helix": 3.0, "LEU-Helix/ALA-Helix": 4.0, "LEU-Helix/VAL-Helix": 3.0, "LEU-Helix/GLU-Helix": 2.0, "LEU-Helix/TYR-Helix": 3.0, "LEU-Helix/MET-Helix": 3.0, "LEU-Helix/CYS-Coil": 1.0, "LEU-Helix/ASP-Coil": -0.0, "LEU-Helix/SER-Coil": 0.0, "LEU-Helix/GLN-Coil": 0.0, "LEU-Helix/LYS-Coil": -3.0, "LEU-Helix/ILE-Coil": 2.0, "LEU-Helix/PRO-Coil": 1.0, "LEU-Helix/THR-Coil": -1.0, "LEU-Helix/PHE-Coil": 2.0, "LEU-Helix/ASN-Coil": -3.0, "LEU-Helix/GLY-Coil": -0.0, "LEU-Helix/HIS-Coil": 0.0, "LEU-Helix/LEU-Coil": 1.0, "LEU-Helix/ARG-Coil": -0.0, "LEU-Helix/TRP-Coil": 1.0, "LEU-Helix/ALA-Coil": 2.0, "LEU-Helix/VAL-Coil": 1.0, "LEU-Helix/GLU-Coil": -1.0, "LEU-Helix/TYR-Coil": 1.0, "LEU-Helix/MET-Coil": 3.0, "LEU-Helix/CYS-Sheet": 2.0, "LEU-Helix/ASP-Sheet": -Infinity, "LEU-Helix/SER-Sheet": -3.0, "LEU-Helix/GLN-Sheet": -1.0, "LEU-Helix/LYS-Sheet": 1.0, "LEU-Helix/ILE-Sheet": 1.0, "LEU-Helix/PRO-Sheet": 1.0, "LEU-Helix/THR-Sheet": -0.0, "LEU-Helix/PHE-Sheet": 1.0, "LEU-Helix/ASN-Sheet": 1.0, "LEU-Helix/GLY-Sheet": 1.0, "LEU-Helix/HIS-Sheet": 0.0, "LEU-Helix/LEU-Sheet": -0.0, "LEU-Helix/ARG-Sheet": -3.0, "LEU-Helix/TRP-Sheet": 2.0, "LEU-Helix/ALA-Sheet": 2.0, "LEU-Helix/VAL-Sheet": -0.0, "LEU-Helix/GLU-Sheet": -0.0, "LEU-Helix/TYR-Sheet": 1.0, "LEU-Helix/MET-Sheet": 2.0, "ARG-Helix/CYS-Helix": 0.0, "ARG-Helix/ASP-Helix": 5.0, "ARG-Helix/SER-Helix": 3.0, "ARG-Helix/GLN-Helix": 3.0, "ARG-Helix/LYS-Helix": 1.0, "ARG-Helix/ILE-Helix": 1.0, "ARG-Helix/PRO-Helix": 2.0, "ARG-Helix/THR-Helix": 2.0, "ARG-Helix/PHE-Helix": 3.0, "ARG-Helix/ASN-Helix": 3.0, "ARG-Helix/GLY-Helix": 2.0, "ARG-Helix/HIS-Helix": 2.0, "ARG-Helix/LEU-Helix": 2.0, "ARG-Helix/ARG-Helix": 1.0, "ARG-Helix/TRP-Helix": 4.0, "ARG-Helix/ALA-Helix": 3.0, "ARG-Helix/VAL-Helix": 2.0, "ARG-Helix/GLU-Helix": 5.0, "ARG-Helix/TYR-Helix": 3.0, "ARG-Helix/MET-Helix": 2.0, "ARG-Helix/CYS-Coil": 2.0, "ARG-Helix/ASP-Coil": 2.0, "ARG-Helix/SER-Coil": 1.0, "ARG-Helix/GLN-Coil": 1.0, "ARG-Helix/LYS-Coil": -1.0, "ARG-Helix/ILE-Coil": 1.0, "ARG-Helix/PRO-Coil": -0.0, "ARG-Helix/THR-Coil": 2.0, "ARG-Helix/PHE-Coil": 2.0, "ARG-Helix/ASN-Coil": 0.0, "ARG-Helix/GLY-Coil": 1.0, "ARG-Helix/HIS-Coil": 2.0, "ARG-Helix/LEU-Coil": 1.0, "ARG-Helix/ARG-Coil": -0.0, "ARG-Helix/TRP-Coil": 2.0, "ARG-Helix/ALA-Coil": -1.0, "ARG-Helix/VAL-Coil": 1.0, "ARG-Helix/GLU-Coil": 3.0, "ARG-Helix/TYR-Coil": 2.0, "ARG-Helix/MET-Coil": 2.0, "ARG-Helix/CYS-Sheet": -Infinity, "ARG-Helix/ASP-Sheet": 3.0, "ARG-Helix/SER-Sheet": -1.0, "ARG-Helix/GLN-Sheet": -1.0, "ARG-Helix/LYS-Sheet": -Infinity, "ARG-Helix/ILE-Sheet": -2.0, "ARG-Helix/PRO-Sheet": 1.0, "ARG-Helix/THR-Sheet": -1.0, "ARG-Helix/PHE-Sheet": -1.0, "ARG-Helix/ASN-Sheet": -Infinity, "ARG-Helix/GLY-Sheet": -1.0, "ARG-Helix/HIS-Sheet": 0.0, "ARG-Helix/LEU-Sheet": -3.0, "ARG-Helix/ARG-Sheet": 0.0, "ARG-Helix/TRP-Sheet": -2.0, "ARG-Helix/ALA-Sheet": -2.0, "ARG-Helix/VAL-Sheet": -4.0, "ARG-Helix/GLU-Sheet": 2.0, "ARG-Helix/TYR-Sheet": -2.0, "ARG-Helix/MET-Sheet": 1.0, "TRP-Helix/CYS-Helix": 1.0, "TRP-Helix/ASP-Helix": 1.0, "TRP-Helix/SER-Helix": 1.0, "TRP-Helix/GLN-Helix": 3.0, "TRP-Helix/LYS-Helix": 4.0, "TRP-Helix/ILE-Helix": 3.0, "TRP-Helix/PRO-Helix": 4.0, "TRP-Helix/THR-Helix": 2.0, "TRP-Helix/PHE-Helix": 4.0, "TRP-Helix/ASN-Helix": 3.0, "TRP-Helix/GLY-Helix": 3.0, "TRP-Helix/HIS-Helix": 4.0, "TRP-Helix/LEU-Helix": 3.0, "TRP-Helix/ARG-Helix": 4.0, "TRP-Helix/TRP-Helix": -2.0, "TRP-Helix/ALA-Helix": 3.0, "TRP-Helix/VAL-Helix": 4.0, "TRP-Helix/GLU-Helix": 1.0, "TRP-Helix/TYR-Helix": 2.0, "TRP-Helix/MET-Helix": 3.0, "TRP-Helix/CYS-Coil": -Infinity, "TRP-Helix/ASP-Coil": -3.0, "TRP-Helix/SER-Coil": -2.0, "TRP-Helix/GLN-Coil": 1.0, "TRP-Helix/LYS-Coil": 0.0, "TRP-Helix/ILE-Coil": 1.0, "TRP-Helix/PRO-Coil": 1.0, "TRP-Helix/THR-Coil": -1.0, "TRP-Helix/PHE-Coil": -0.0, "TRP-Helix/ASN-Coil": 1.0, "TRP-Helix/GLY-Coil": -0.0, "TRP-Helix/HIS-Coil": 2.0, "TRP-Helix/LEU-Coil": 1.0, "TRP-Helix/ARG-Coil": -1.0, "TRP-Helix/TRP-Coil": 1.0, "TRP-Helix/ALA-Coil": 0.0, "TRP-Helix/VAL-Coil": 2.0, "TRP-Helix/GLU-Coil": -2.0, "TRP-Helix/TYR-Coil": 1.0, "TRP-Helix/MET-Coil": 3.0, "TRP-Helix/CYS-Sheet": -Infinity, "TRP-Helix/ASP-Sheet": 1.0, "TRP-Helix/SER-Sheet": -2.0, "TRP-Helix/GLN-Sheet": -Infinity, "TRP-Helix/LYS-Sheet": 0.0, "TRP-Helix/ILE-Sheet": 1.0, "TRP-Helix/PRO-Sheet": 2.0, "TRP-Helix/THR-Sheet": -Infinity, "TRP-Helix/PHE-Sheet": -1.0, "TRP-Helix/ASN-Sheet": -0.0, "TRP-Helix/GLY-Sheet": -0.0, "TRP-Helix/HIS-Sheet": 0.0, "TRP-Helix/LEU-Sheet": -1.0, "TRP-Helix/ARG-Sheet": 1.0, "TRP-Helix/TRP-Sheet": -Infinity, "TRP-Helix/ALA-Sheet": 1.0, "TRP-Helix/VAL-Sheet": -1.0, "TRP-Helix/GLU-Sheet": -Infinity, "TRP-Helix/TYR-Sheet": -Infinity, "TRP-Helix/MET-Sheet": 0.0, "ALA-Helix/CYS-Helix": 4.0, "ALA-Helix/ASP-Helix": 3.0, "ALA-Helix/SER-Helix": 3.0, "ALA-Helix/GLN-Helix": 2.0, "ALA-Helix/LYS-Helix": 1.0, "ALA-Helix/ILE-Helix": 3.0, "ALA-Helix/PRO-Helix": 3.0, "ALA-Helix/THR-Helix": 3.0, "ALA-Helix/PHE-Helix": 3.0, "ALA-Helix/ASN-Helix": 2.0, "ALA-Helix/GLY-Helix": 5.0, "ALA-Helix/HIS-Helix": 2.0, "ALA-Helix/LEU-Helix": 4.0, "ALA-Helix/ARG-Helix": 3.0, "ALA-Helix/TRP-Helix": 3.0, "ALA-Helix/ALA-Helix": 3.0, "ALA-Helix/VAL-Helix": 4.0, "ALA-Helix/GLU-Helix": 2.0, "ALA-Helix/TYR-Helix": 3.0, "ALA-Helix/MET-Helix": 4.0, "ALA-Helix/CYS-Coil": 1.0, "ALA-Helix/ASP-Coil": 1.0, "ALA-Helix/SER-Coil": 1.0, "ALA-Helix/GLN-Coil": 0.0, "ALA-Helix/LYS-Coil": -0.0, "ALA-Helix/ILE-Coil": 2.0, "ALA-Helix/PRO-Coil": 1.0, "ALA-Helix/THR-Coil": 1.0, "ALA-Helix/PHE-Coil": 2.0, "ALA-Helix/ASN-Coil": 0.0, "ALA-Helix/GLY-Coil": 1.0, "ALA-Helix/HIS-Coil": 0.0, "ALA-Helix/LEU-Coil": 2.0, "ALA-Helix/ARG-Coil": -0.0, "ALA-Helix/TRP-Coil": 2.0, "ALA-Helix/ALA-Coil": 2.0, "ALA-Helix/VAL-Coil": 2.0, "ALA-Helix/GLU-Coil": 1.0, "ALA-Helix/TYR-Coil": 2.0, "ALA-Helix/MET-Coil": 2.0, "ALA-Helix/CYS-Sheet": 0.0, "ALA-Helix/ASP-Sheet": -1.0, "ALA-Helix/SER-Sheet": 0.0, "ALA-Helix/GLN-Sheet": -4.0, "ALA-Helix/LYS-Sheet": 0.0, "ALA-Helix/ILE-Sheet": 1.0, "ALA-Helix/PRO-Sheet": 0.0, "ALA-Helix/THR-Sheet": -1.0, "ALA-Helix/PHE-Sheet": 1.0, "ALA-Helix/ASN-Sheet": 0.0, "ALA-Helix/GLY-Sheet": -0.0, "ALA-Helix/HIS-Sheet": 0.0, "ALA-Helix/LEU-Sheet": 1.0, "ALA-Helix/ARG-Sheet": -1.0, "ALA-Helix/TRP-Sheet": 0.0, "ALA-Helix/ALA-Sheet": 1.0, "ALA-Helix/VAL-Sheet": 2.0, "ALA-Helix/GLU-Sheet": -1.0, "ALA-Helix/TYR-Sheet": 0.0, "ALA-Helix/MET-Sheet": 2.0, "VAL-Helix/CYS-Helix": 4.0, "VAL-Helix/ASP-Helix": 1.0, "VAL-Helix/SER-Helix": 3.0, "VAL-Helix/GLN-Helix": 2.0, "VAL-Helix/LYS-Helix": -0.0, "VAL-Helix/ILE-Helix": 3.0, "VAL-Helix/PRO-Helix": 3.0, "VAL-Helix/THR-Helix": 3.0, "VAL-Helix/PHE-Helix": 3.0, "VAL-Helix/ASN-Helix": 2.0, "VAL-Helix/GLY-Helix": 4.0, "VAL-Helix/HIS-Helix": 2.0, "VAL-Helix/LEU-Helix": 3.0, "VAL-Helix/ARG-Helix": 2.0, "VAL-Helix/TRP-Helix": 4.0, "VAL-Helix/ALA-Helix": 4.0, "VAL-Helix/VAL-Helix": 2.0, "VAL-Helix/GLU-Helix": 1.0, "VAL-Helix/TYR-Helix": 3.0, "VAL-Helix/MET-Helix": 3.0, "VAL-Helix/CYS-Coil": 2.0, "VAL-Helix/ASP-Coil": -1.0, "VAL-Helix/SER-Coil": 0.0, "VAL-Helix/GLN-Coil": 0.0, "VAL-Helix/LYS-Coil": -4.0, "VAL-Helix/ILE-Coil": 2.0, "VAL-Helix/PRO-Coil": 1.0, "VAL-Helix/THR-Coil": 1.0, "VAL-Helix/PHE-Coil": 0.0, "VAL-Helix/ASN-Coil": 1.0, "VAL-Helix/GLY-Coil": 1.0, "VAL-Helix/HIS-Coil": -1.0, "VAL-Helix/LEU-Coil": 1.0, "VAL-Helix/ARG-Coil": -2.0, "VAL-Helix/TRP-Coil": 2.0, "VAL-Helix/ALA-Coil": 2.0, "VAL-Helix/VAL-Coil": 1.0, "VAL-Helix/GLU-Coil": -2.0, "VAL-Helix/TYR-Coil": 0.0, "VAL-Helix/MET-Coil": 2.0, "VAL-Helix/CYS-Sheet": 0.0, "VAL-Helix/ASP-Sheet": -1.0, "VAL-Helix/SER-Sheet": 0.0, "VAL-Helix/GLN-Sheet": -Infinity, "VAL-Helix/LYS-Sheet": -Infinity, "VAL-Helix/ILE-Sheet": 2.0, "VAL-Helix/PRO-Sheet": 2.0, "VAL-Helix/THR-Sheet": -0.0, "VAL-Helix/PHE-Sheet": 1.0, "VAL-Helix/ASN-Sheet": -1.0, "VAL-Helix/GLY-Sheet": -0.0, "VAL-Helix/HIS-Sheet": 1.0, "VAL-Helix/LEU-Sheet": -2.0, "VAL-Helix/ARG-Sheet": -2.0, "VAL-Helix/TRP-Sheet": -Infinity, "VAL-Helix/ALA-Sheet": 2.0, "VAL-Helix/VAL-Sheet": 1.0, "VAL-Helix/GLU-Sheet": -Infinity, "VAL-Helix/TYR-Sheet": 1.0, "VAL-Helix/MET-Sheet": 1.0, "GLU-Helix/CYS-Helix": 0.0, "GLU-Helix/ASP-Helix": 2.0, "GLU-Helix/SER-Helix": 3.0, "GLU-Helix/GLN-Helix": 4.0, "GLU-Helix/LYS-Helix": 5.0, "GLU-Helix/ILE-Helix": 1.0, "GLU-Helix/PRO-Helix": 3.0, "GLU-Helix/THR-Helix": 3.0, "GLU-Helix/PHE-Helix": 1.0, "GLU-Helix/ASN-Helix": 3.0, "GLU-Helix/GLY-Helix": 1.0, "GLU-Helix/HIS-Helix": 3.0, "GLU-Helix/LEU-Helix": 2.0, "GLU-Helix/ARG-Helix": 5.0, "GLU-Helix/TRP-Helix": 1.0, "GLU-Helix/ALA-Helix": 2.0, "GLU-Helix/VAL-Helix": 1.0, "GLU-Helix/GLU-Helix": 1.0, "GLU-Helix/TYR-Helix": 2.0, "GLU-Helix/MET-Helix": 2.0, "GLU-Helix/CYS-Coil": -Infinity, "GLU-Helix/ASP-Coil": -3.0, "GLU-Helix/SER-Coil": 2.0, "GLU-Helix/GLN-Coil": 1.0, "GLU-Helix/LYS-Coil": 2.0, "GLU-Helix/ILE-Coil": 1.0, "GLU-Helix/PRO-Coil": -0.0, "GLU-Helix/THR-Coil": 3.0, "GLU-Helix/PHE-Coil": 0.0, "GLU-Helix/ASN-Coil": 2.0, "GLU-Helix/GLY-Coil": 1.0, "GLU-Helix/HIS-Coil": 2.0, "GLU-Helix/LEU-Coil": 2.0, "GLU-Helix/ARG-Coil": 3.0, "GLU-Helix/TRP-Coil": 1.0, "GLU-Helix/ALA-Coil": 1.0, "GLU-Helix/VAL-Coil": 1.0, "GLU-Helix/GLU-Coil": 0.0, "GLU-Helix/TYR-Coil": 0.0, "GLU-Helix/MET-Coil": 1.0, "GLU-Helix/CYS-Sheet": -Infinity, "GLU-Helix/ASP-Sheet": -2.0, "GLU-Helix/SER-Sheet": -1.0, "GLU-Helix/GLN-Sheet": -0.0, "GLU-Helix/LYS-Sheet": 0.0, "GLU-Helix/ILE-Sheet": -1.0, "GLU-Helix/PRO-Sheet": 3.0, "GLU-Helix/THR-Sheet": 0.0, "GLU-Helix/PHE-Sheet": -Infinity, "GLU-Helix/ASN-Sheet": 0.0, "GLU-Helix/GLY-Sheet": -0.0, "GLU-Helix/HIS-Sheet": -1.0, "GLU-Helix/LEU-Sheet": -1.0, "GLU-Helix/ARG-Sheet": 0.0, "GLU-Helix/TRP-Sheet": -Infinity, "GLU-Helix/ALA-Sheet": 0.0, "GLU-Helix/VAL-Sheet": -2.0, "GLU-Helix/GLU-Sheet": -Infinity, "GLU-Helix/TYR-Sheet": -3.0, "GLU-Helix/MET-Sheet": -Infinity, "TYR-Helix/CYS-Helix": 3.0, "TYR-Helix/ASP-Helix": 0.0, "TYR-Helix/SER-Helix": 2.0, "TYR-Helix/GLN-Helix": 4.0, "TYR-Helix/LYS-Helix": 4.0, "TYR-Helix/ILE-Helix": 3.0, "TYR-Helix/PRO-Helix": 3.0, "TYR-Helix/THR-Helix": 3.0, "TYR-Helix/PHE-Helix": 3.0, "TYR-Helix/ASN-Helix": 3.0, "TYR-Helix/GLY-Helix": 2.0, "TYR-Helix/HIS-Helix": 3.0, "TYR-Helix/LEU-Helix": 3.0, "TYR-Helix/ARG-Helix": 3.0, "TYR-Helix/TRP-Helix": 2.0, "TYR-Helix/ALA-Helix": 3.0, "TYR-Helix/VAL-Helix": 3.0, "TYR-Helix/GLU-Helix": 2.0, "TYR-Helix/TYR-Helix": 2.0, "TYR-Helix/MET-Helix": 3.0, "TYR-Helix/CYS-Coil": 1.0, "TYR-Helix/ASP-Coil": 0.0, "TYR-Helix/SER-Coil": 0.0, "TYR-Helix/GLN-Coil": -2.0, "TYR-Helix/LYS-Coil": 2.0, "TYR-Helix/ILE-Coil": 3.0, "TYR-Helix/PRO-Coil": 3.0, "TYR-Helix/THR-Coil": -2.0, "TYR-Helix/PHE-Coil": 1.0, "TYR-Helix/ASN-Coil": 2.0, "TYR-Helix/GLY-Coil": 1.0, "TYR-Helix/HIS-Coil": -1.0, "TYR-Helix/LEU-Coil": 2.0, "TYR-Helix/ARG-Coil": 1.0, "TYR-Helix/TRP-Coil": -Infinity, "TYR-Helix/ALA-Coil": 1.0, "TYR-Helix/VAL-Coil": 2.0, "TYR-Helix/GLU-Coil": -2.0, "TYR-Helix/TYR-Coil": 0.0, "TYR-Helix/MET-Coil": 2.0, "TYR-Helix/CYS-Sheet": 3.0, "TYR-Helix/ASP-Sheet": -0.0, "TYR-Helix/SER-Sheet": 1.0, "TYR-Helix/GLN-Sheet": -Infinity, "TYR-Helix/LYS-Sheet": -2.0, "TYR-Helix/ILE-Sheet": -1.0, "TYR-Helix/PRO-Sheet": 2.0, "TYR-Helix/THR-Sheet": -Infinity, "TYR-Helix/PHE-Sheet": -2.0, "TYR-Helix/ASN-Sheet": -0.0, "TYR-Helix/GLY-Sheet": 0.0, "TYR-Helix/HIS-Sheet": 1.0, "TYR-Helix/LEU-Sheet": 1.0, "TYR-Helix/ARG-Sheet": -1.0, "TYR-Helix/TRP-Sheet": -1.0, "TYR-Helix/ALA-Sheet": 1.0, "TYR-Helix/VAL-Sheet": 0.0, "TYR-Helix/GLU-Sheet": -Infinity, "TYR-Helix/TYR-Sheet": 1.0, "TYR-Helix/MET-Sheet": -1.0, "MET-Helix/CYS-Helix": 4.0, "MET-Helix/ASP-Helix": 0.0, "MET-Helix/SER-Helix": 3.0, "MET-Helix/GLN-Helix": 2.0, "MET-Helix/LYS-Helix": 2.0, "MET-Helix/ILE-Helix": 3.0, "MET-Helix/PRO-Helix": 3.0, "MET-Helix/THR-Helix": 3.0, "MET-Helix/PHE-Helix": 4.0, "MET-Helix/ASN-Helix": 3.0, "MET-Helix/GLY-Helix": 4.0, "MET-Helix/HIS-Helix": 3.0, "MET-Helix/LEU-Helix": 3.0, "MET-Helix/ARG-Helix": 2.0, "MET-Helix/TRP-Helix": 3.0, "MET-Helix/ALA-Helix": 4.0, "MET-Helix/VAL-Helix": 3.0, "MET-Helix/GLU-Helix": 2.0, "MET-Helix/TYR-Helix": 3.0, "MET-Helix/MET-Helix": 2.0, "MET-Helix/CYS-Coil": 3.0, "MET-Helix/ASP-Coil": -2.0, "MET-Helix/SER-Coil": 1.0, "MET-Helix/GLN-Coil": -2.0, "MET-Helix/LYS-Coil": 2.0, "MET-Helix/ILE-Coil": -2.0, "MET-Helix/PRO-Coil": 1.0, "MET-Helix/THR-Coil": -1.0, "MET-Helix/PHE-Coil": 1.0, "MET-Helix/ASN-Coil": -2.0, "MET-Helix/GLY-Coil": 0.0, "MET-Helix/HIS-Coil": -Infinity, "MET-Helix/LEU-Coil": 1.0, "MET-Helix/ARG-Coil": -Infinity, "MET-Helix/TRP-Coil": 3.0, "MET-Helix/ALA-Coil": 3.0, "MET-Helix/VAL-Coil": 3.0, "MET-Helix/GLU-Coil": 0.0, "MET-Helix/TYR-Coil": 2.0, "MET-Helix/MET-Coil": 2.0, "MET-Helix/CYS-Sheet": 1.0, "MET-Helix/ASP-Sheet": -Infinity, "MET-Helix/SER-Sheet": -2.0, "MET-Helix/GLN-Sheet": -Infinity, "MET-Helix/LYS-Sheet": -Infinity, "MET-Helix/ILE-Sheet": 1.0, "MET-Helix/PRO-Sheet": 3.0, "MET-Helix/THR-Sheet": -3.0, "MET-Helix/PHE-Sheet": 1.0, "MET-Helix/ASN-Sheet": 1.0, "MET-Helix/GLY-Sheet": 2.0, "MET-Helix/HIS-Sheet": 2.0, "MET-Helix/LEU-Sheet": 0.0, "MET-Helix/ARG-Sheet": -Infinity, "MET-Helix/TRP-Sheet": 2.0, "MET-Helix/ALA-Sheet": 3.0, "MET-Helix/VAL-Sheet": 1.0, "MET-Helix/GLU-Sheet": -Infinity, "MET-Helix/TYR-Sheet": 2.0, "MET-Helix/MET-Sheet": 1.0, "CYS-Coil/CYS-Helix": 8.0, "CYS-Coil/ASP-Helix": -1.0, "CYS-Coil/SER-Helix": 0.0, "CYS-Coil/GLN-Helix": 3.0, "CYS-Coil/LYS-Helix": -1.0, "CYS-Coil/ILE-Helix": 3.0, "CYS-Coil/PRO-Helix": 2.0, "CYS-Coil/THR-Helix": 1.0, "CYS-Coil/PHE-Helix": 3.0, "CYS-Coil/ASN-Helix": -1.0, "CYS-Coil/GLY-Helix": 3.0, "CYS-Coil/HIS-Helix": 4.0, "CYS-Coil/LEU-Helix": 1.0, "CYS-Coil/ARG-Helix": 2.0, "CYS-Coil/TRP-Helix": -Infinity, "CYS-Coil/ALA-Helix": 1.0, "CYS-Coil/VAL-Helix": 2.0, "CYS-Coil/GLU-Helix": -Infinity, "CYS-Coil/TYR-Helix": 1.0, "CYS-Coil/MET-Helix": 3.0, "CYS-Coil/CYS-Coil": 8.0, "CYS-Coil/ASP-Coil": 1.0, "CYS-Coil/SER-Coil": 4.0, "CYS-Coil/GLN-Coil": 2.0, "CYS-Coil/LYS-Coil": -0.0, "CYS-Coil/ILE-Coil": 3.0, "CYS-Coil/PRO-Coil": 3.0, "CYS-Coil/THR-Coil": 3.0, "CYS-Coil/PHE-Coil": 3.0, "CYS-Coil/ASN-Coil": 3.0, "CYS-Coil/GLY-Coil": 5.0, "CYS-Coil/HIS-Coil": 4.0, "CYS-Coil/LEU-Coil": 3.0, "CYS-Coil/ARG-Coil": 2.0, "CYS-Coil/TRP-Coil": 4.0, "CYS-Coil/ALA-Coil": 3.0, "CYS-Coil/VAL-Coil": 4.0, "CYS-Coil/GLU-Coil": 3.0, "CYS-Coil/TYR-Coil": 1.0, "CYS-Coil/MET-Coil": 4.0, "CYS-Coil/CYS-Sheet": 7.0, "CYS-Coil/ASP-Sheet": -Infinity, "CYS-Coil/SER-Sheet": -0.0, "CYS-Coil/GLN-Sheet": -Infinity, "CYS-Coil/LYS-Sheet": 0.0, "CYS-Coil/ILE-Sheet": 2.0, "CYS-Coil/PRO-Sheet": 4.0, "CYS-Coil/THR-Sheet": 2.0, "CYS-Coil/PHE-Sheet": 3.0, "CYS-Coil/ASN-Sheet": -Infinity, "CYS-Coil/GLY-Sheet": 2.0, "CYS-Coil/HIS-Sheet": 5.0, "CYS-Coil/LEU-Sheet": 3.0, "CYS-Coil/ARG-Sheet": 2.0, "CYS-Coil/TRP-Sheet": -Infinity, "CYS-Coil/ALA-Sheet": 3.0, "CYS-Coil/VAL-Sheet": 2.0, "CYS-Coil/GLU-Sheet": 2.0, "CYS-Coil/TYR-Sheet": 1.0, "CYS-Coil/MET-Sheet": 2.0, "ASP-Coil/CYS-Helix": -Infinity, "ASP-Coil/ASP-Helix": 0.0, "ASP-Coil/SER-Helix": 2.0, "ASP-Coil/GLN-Helix": 1.0, "ASP-Coil/LYS-Helix": 3.0, "ASP-Coil/ILE-Helix": 1.0, "ASP-Coil/PRO-Helix": 1.0, "ASP-Coil/THR-Helix": 1.0, "ASP-Coil/PHE-Helix": -4.0, "ASP-Coil/ASN-Helix": 2.0, "ASP-Coil/GLY-Helix": 1.0, "ASP-Coil/HIS-Helix": 2.0, "ASP-Coil/LEU-Helix": -0.0, "ASP-Coil/ARG-Helix": 2.0, "ASP-Coil/TRP-Helix": -3.0, "ASP-Coil/ALA-Helix": 1.0, "ASP-Coil/VAL-Helix": -1.0, "ASP-Coil/GLU-Helix": -3.0, "ASP-Coil/TYR-Helix": 0.0, "ASP-Coil/MET-Helix": -2.0, "ASP-Coil/CYS-Coil": 1.0, "ASP-Coil/ASP-Coil": 2.0, "ASP-Coil/SER-Coil": 4.0, "ASP-Coil/GLN-Coil": 4.0, "ASP-Coil/LYS-Coil": 5.0, "ASP-Coil/ILE-Coil": 3.0, "ASP-Coil/PRO-Coil": 3.0, "ASP-Coil/THR-Coil": 4.0, "ASP-Coil/PHE-Coil": 2.0, "ASP-Coil/ASN-Coil": 4.0, "ASP-Coil/GLY-Coil": 4.0, "ASP-Coil/HIS-Coil": 3.0, "ASP-Coil/LEU-Coil": 2.0, "ASP-Coil/ARG-Coil": 5.0, "ASP-Coil/TRP-Coil": 2.0, "ASP-Coil/ALA-Coil": 4.0, "ASP-Coil/VAL-Coil": 3.0, "ASP-Coil/GLU-Coil": 2.0, "ASP-Coil/TYR-Coil": 2.0, "ASP-Coil/MET-Coil": 1.0, "ASP-Coil/CYS-Sheet": 2.0, "ASP-Coil/ASP-Sheet": -0.0, "ASP-Coil/SER-Sheet": 1.0, "ASP-Coil/GLN-Sheet": 3.0, "ASP-Coil/LYS-Sheet": 3.0, "ASP-Coil/ILE-Sheet": 1.0, "ASP-Coil/PRO-Sheet": -Infinity, "ASP-Coil/THR-Sheet": 3.0, "ASP-Coil/PHE-Sheet": -1.0, "ASP-Coil/ASN-Sheet": 3.0, "ASP-Coil/GLY-Sheet": 2.0, "ASP-Coil/HIS-Sheet": 4.0, "ASP-Coil/LEU-Sheet": 1.0, "ASP-Coil/ARG-Sheet": 4.0, "ASP-Coil/TRP-Sheet": -Infinity, "ASP-Coil/ALA-Sheet": 1.0, "ASP-Coil/VAL-Sheet": -4.0, "ASP-Coil/GLU-Sheet": -1.0, "ASP-Coil/TYR-Sheet": 1.0, "ASP-Coil/MET-Sheet": 0.0, "SER-Coil/CYS-Helix": 1.0, "SER-Coil/ASP-Helix": 3.0, "SER-Coil/SER-Helix": 1.0, "SER-Coil/GLN-Helix": 2.0, "SER-Coil/LYS-Helix": 1.0, "SER-Coil/ILE-Helix": 0.0, "SER-Coil/PRO-Helix": 2.0, "SER-Coil/THR-Helix": -0.0, "SER-Coil/PHE-Helix": -1.0, "SER-Coil/ASN-Helix": 1.0, "SER-Coil/GLY-Helix": 2.0, "SER-Coil/HIS-Helix": 1.0, "SER-Coil/LEU-Helix": 0.0, "SER-Coil/ARG-Helix": 1.0, "SER-Coil/TRP-Helix": -2.0, "SER-Coil/ALA-Helix": 1.0, "SER-Coil/VAL-Helix": 0.0, "SER-Coil/GLU-Helix": 2.0, "SER-Coil/TYR-Helix": 0.0, "SER-Coil/MET-Helix": 1.0, "SER-Coil/CYS-Coil": 4.0, "SER-Coil/ASP-Coil": 4.0, "SER-Coil/SER-Coil": 3.0, "SER-Coil/GLN-Coil": 4.0, "SER-Coil/LYS-Coil": 4.0, "SER-Coil/ILE-Coil": 3.0, "SER-Coil/PRO-Coil": 3.0, "SER-Coil/THR-Coil": 4.0, "SER-Coil/PHE-Coil": 3.0, "SER-Coil/ASN-Coil": 4.0, "SER-Coil/GLY-Coil": 5.0, "SER-Coil/HIS-Coil": 4.0, "SER-Coil/LEU-Coil": 3.0, "SER-Coil/ARG-Coil": 3.0, "SER-Coil/TRP-Coil": 2.0, "SER-Coil/ALA-Coil": 4.0, "SER-Coil/VAL-Coil": 3.0, "SER-Coil/GLU-Coil": 3.0, "SER-Coil/TYR-Coil": 1.0, "SER-Coil/MET-Coil": 2.0, "SER-Coil/CYS-Sheet": 2.0, "SER-Coil/ASP-Sheet": 4.0, "SER-Coil/SER-Sheet": 2.0, "SER-Coil/GLN-Sheet": 2.0, "SER-Coil/LYS-Sheet": 3.0, "SER-Coil/ILE-Sheet": 2.0, "SER-Coil/PRO-Sheet": 2.0, "SER-Coil/THR-Sheet": 2.0, "SER-Coil/PHE-Sheet": 2.0, "SER-Coil/ASN-Sheet": 4.0, "SER-Coil/GLY-Sheet": 3.0, "SER-Coil/HIS-Sheet": 3.0, "SER-Coil/LEU-Sheet": 1.0, "SER-Coil/ARG-Sheet": 2.0, "SER-Coil/TRP-Sheet": 2.0, "SER-Coil/ALA-Sheet": 2.0, "SER-Coil/VAL-Sheet": 1.0, "SER-Coil/GLU-Sheet": 2.0, "SER-Coil/TYR-Sheet": 1.0, "SER-Coil/MET-Sheet": 2.0, "GLN-Coil/CYS-Helix": -0.0, "GLN-Coil/ASP-Helix": 0.0, "GLN-Coil/SER-Helix": 1.0, "GLN-Coil/GLN-Helix": 2.0, "GLN-Coil/LYS-Helix": 3.0, "GLN-Coil/ILE-Helix": 0.0, "GLN-Coil/PRO-Helix": 1.0, "GLN-Coil/THR-Helix": -1.0, "GLN-Coil/PHE-Helix": -3.0, "GLN-Coil/ASN-Helix": 2.0, "GLN-Coil/GLY-Helix": 1.0, "GLN-Coil/HIS-Helix": 2.0, "GLN-Coil/LEU-Helix": 0.0, "GLN-Coil/ARG-Helix": 1.0, "GLN-Coil/TRP-Helix": 1.0, "GLN-Coil/ALA-Helix": 0.0, "GLN-Coil/VAL-Helix": 0.0, "GLN-Coil/GLU-Helix": 1.0, "GLN-Coil/TYR-Helix": -2.0, "GLN-Coil/MET-Helix": -2.0, "GLN-Coil/CYS-Coil": 2.0, "GLN-Coil/ASP-Coil": 4.0, "GLN-Coil/SER-Coil": 4.0, "GLN-Coil/GLN-Coil": 3.0, "GLN-Coil/LYS-Coil": 3.0, "GLN-Coil/ILE-Coil": 3.0, "GLN-Coil/PRO-Coil": 2.0, "GLN-Coil/THR-Coil": 3.0, "GLN-Coil/PHE-Coil": 2.0, "GLN-Coil/ASN-Coil": 4.0, "GLN-Coil/GLY-Coil": 3.0, "GLN-Coil/HIS-Coil": 4.0, "GLN-Coil/LEU-Coil": 3.0, "GLN-Coil/ARG-Coil": 4.0, "GLN-Coil/TRP-Coil": 4.0, "GLN-Coil/ALA-Coil": 3.0, "GLN-Coil/VAL-Coil": 3.0, "GLN-Coil/GLU-Coil": 4.0, "GLN-Coil/TYR-Coil": 4.0, "GLN-Coil/MET-Coil": 3.0, "GLN-Coil/CYS-Sheet": -Infinity, "GLN-Coil/ASP-Sheet": 1.0, "GLN-Coil/SER-Sheet": 2.0, "GLN-Coil/GLN-Sheet": -0.0, "GLN-Coil/LYS-Sheet": 1.0, "GLN-Coil/ILE-Sheet": 1.0, "GLN-Coil/PRO-Sheet": 1.0, "GLN-Coil/THR-Sheet": 1.0, "GLN-Coil/PHE-Sheet": -1.0, "GLN-Coil/ASN-Sheet": 3.0, "GLN-Coil/GLY-Sheet": -1.0, "GLN-Coil/HIS-Sheet": 3.0, "GLN-Coil/LEU-Sheet": 1.0, "GLN-Coil/ARG-Sheet": 3.0, "GLN-Coil/TRP-Sheet": 2.0, "GLN-Coil/ALA-Sheet": 2.0, "GLN-Coil/VAL-Sheet": 2.0, "GLN-Coil/GLU-Sheet": 3.0, "GLN-Coil/TYR-Sheet": 2.0, "GLN-Coil/MET-Sheet": 1.0, "LYS-Coil/CYS-Helix": -1.0, "LYS-Coil/ASP-Helix": 2.0, "LYS-Coil/SER-Helix": 1.0, "LYS-Coil/GLN-Helix": 0.0, "LYS-Coil/LYS-Helix": 0.0, "LYS-Coil/ILE-Helix": -3.0, "LYS-Coil/PRO-Helix": 1.0, "LYS-Coil/THR-Helix": 0.0, "LYS-Coil/PHE-Helix": -2.0, "LYS-Coil/ASN-Helix": 2.0, "LYS-Coil/GLY-Helix": -1.0, "LYS-Coil/HIS-Helix": 0.0, "LYS-Coil/LEU-Helix": -3.0, "LYS-Coil/ARG-Helix": -1.0, "LYS-Coil/TRP-Helix": 0.0, "LYS-Coil/ALA-Helix": -0.0, "LYS-Coil/VAL-Helix": -4.0, "LYS-Coil/GLU-Helix": 2.0, "LYS-Coil/TYR-Helix": 2.0, "LYS-Coil/MET-Helix": 2.0, "LYS-Coil/CYS-Coil": -0.0, "LYS-Coil/ASP-Coil": 5.0, "LYS-Coil/SER-Coil": 4.0, "LYS-Coil/GLN-Coil": 3.0, "LYS-Coil/LYS-Coil": -1.0, "LYS-Coil/ILE-Coil": 1.0, "LYS-Coil/PRO-Coil": 2.0, "LYS-Coil/THR-Coil": 1.0, "LYS-Coil/PHE-Coil": 2.0, "LYS-Coil/ASN-Coil": 4.0, "LYS-Coil/GLY-Coil": 2.0, "LYS-Coil/HIS-Coil": 2.0, "LYS-Coil/LEU-Coil": 1.0, "LYS-Coil/ARG-Coil": 2.0, "LYS-Coil/TRP-Coil": 4.0, "LYS-Coil/ALA-Coil": 3.0, "LYS-Coil/VAL-Coil": 2.0, "LYS-Coil/GLU-Coil": 5.0, "LYS-Coil/TYR-Coil": 4.0, "LYS-Coil/MET-Coil": 4.0, "LYS-Coil/CYS-Sheet": -Infinity, "LYS-Coil/ASP-Sheet": 4.0, "LYS-Coil/SER-Sheet": 0.0, "LYS-Coil/GLN-Sheet": -Infinity, "LYS-Coil/LYS-Sheet": -1.0, "LYS-Coil/ILE-Sheet": -2.0, "LYS-Coil/PRO-Sheet": 3.0, "LYS-Coil/THR-Sheet": -1.0, "LYS-Coil/PHE-Sheet": -0.0, "LYS-Coil/ASN-Sheet": -0.0, "LYS-Coil/GLY-Sheet": -Infinity, "LYS-Coil/HIS-Sheet": 2.0, "LYS-Coil/LEU-Sheet": -1.0, "LYS-Coil/ARG-Sheet": -1.0, "LYS-Coil/TRP-Sheet": -0.0, "LYS-Coil/ALA-Sheet": -2.0, "LYS-Coil/VAL-Sheet": 0.0, "LYS-Coil/GLU-Sheet": 4.0, "LYS-Coil/TYR-Sheet": 3.0, "LYS-Coil/MET-Sheet": 3.0, "ILE-Coil/CYS-Helix": 2.0, "ILE-Coil/ASP-Helix": 2.0, "ILE-Coil/SER-Helix": 3.0, "ILE-Coil/GLN-Helix": 1.0, "ILE-Coil/LYS-Helix": 1.0, "ILE-Coil/ILE-Helix": 2.0, "ILE-Coil/PRO-Helix": 2.0, "ILE-Coil/THR-Helix": 2.0, "ILE-Coil/PHE-Helix": 2.0, "ILE-Coil/ASN-Helix": 2.0, "ILE-Coil/GLY-Helix": 2.0, "ILE-Coil/HIS-Helix": 3.0, "ILE-Coil/LEU-Helix": 2.0, "ILE-Coil/ARG-Helix": 1.0, "ILE-Coil/TRP-Helix": 1.0, "ILE-Coil/ALA-Helix": 2.0, "ILE-Coil/VAL-Helix": 2.0, "ILE-Coil/GLU-Helix": 1.0, "ILE-Coil/TYR-Helix": 3.0, "ILE-Coil/MET-Helix": -2.0, "ILE-Coil/CYS-Coil": 3.0, "ILE-Coil/ASP-Coil": 3.0, "ILE-Coil/SER-Coil": 3.0, "ILE-Coil/GLN-Coil": 3.0, "ILE-Coil/LYS-Coil": 1.0, "ILE-Coil/ILE-Coil": 0.0, "ILE-Coil/PRO-Coil": 3.0, "ILE-Coil/THR-Coil": 2.0, "ILE-Coil/PHE-Coil": 3.0, "ILE-Coil/ASN-Coil": 3.0, "ILE-Coil/GLY-Coil": 2.0, "ILE-Coil/HIS-Coil": 2.0, "ILE-Coil/LEU-Coil": 2.0, "ILE-Coil/ARG-Coil": 3.0, "ILE-Coil/TRP-Coil": 5.0, "ILE-Coil/ALA-Coil": 3.0, "ILE-Coil/VAL-Coil": 3.0, "ILE-Coil/GLU-Coil": 1.0, "ILE-Coil/TYR-Coil": 3.0, "ILE-Coil/MET-Coil": -Infinity, "ILE-Coil/CYS-Sheet": 3.0, "ILE-Coil/ASP-Sheet": -Infinity, "ILE-Coil/SER-Sheet": 1.0, "ILE-Coil/GLN-Sheet": 1.0, "ILE-Coil/LYS-Sheet": -Infinity, "ILE-Coil/ILE-Sheet": 2.0, "ILE-Coil/PRO-Sheet": 3.0, "ILE-Coil/THR-Sheet": 0.0, "ILE-Coil/PHE-Sheet": 2.0, "ILE-Coil/ASN-Sheet": 3.0, "ILE-Coil/GLY-Sheet": 3.0, "ILE-Coil/HIS-Sheet": 2.0, "ILE-Coil/LEU-Sheet": 3.0, "ILE-Coil/ARG-Sheet": 1.0, "ILE-Coil/TRP-Sheet": -Infinity, "ILE-Coil/ALA-Sheet": 3.0, "ILE-Coil/VAL-Sheet": 1.0, "ILE-Coil/GLU-Sheet": 2.0, "ILE-Coil/TYR-Sheet": 1.0, "ILE-Coil/MET-Sheet": 2.0, "PRO-Coil/CYS-Helix": 3.0, "PRO-Coil/ASP-Helix": 2.0, "PRO-Coil/SER-Helix": 2.0, "PRO-Coil/GLN-Helix": 1.0, "PRO-Coil/LYS-Helix": -0.0, "PRO-Coil/ILE-Helix": 1.0, "PRO-Coil/PRO-Helix": 1.0, "PRO-Coil/THR-Helix": 1.0, "PRO-Coil/PHE-Helix": 1.0, "PRO-Coil/ASN-Helix": 0.0, "PRO-Coil/GLY-Helix": 1.0, "PRO-Coil/HIS-Helix": 2.0, "PRO-Coil/LEU-Helix": 1.0, "PRO-Coil/ARG-Helix": -0.0, "PRO-Coil/TRP-Helix": 1.0, "PRO-Coil/ALA-Helix": 1.0, "PRO-Coil/VAL-Helix": 1.0, "PRO-Coil/GLU-Helix": -0.0, "PRO-Coil/TYR-Helix": 3.0, "PRO-Coil/MET-Helix": 1.0, "PRO-Coil/CYS-Coil": 3.0, "PRO-Coil/ASP-Coil": 3.0, "PRO-Coil/SER-Coil": 3.0, "PRO-Coil/GLN-Coil": 2.0, "PRO-Coil/LYS-Coil": 2.0, "PRO-Coil/ILE-Coil": 3.0, "PRO-Coil/PRO-Coil": 1.0, "PRO-Coil/THR-Coil": 3.0, "PRO-Coil/PHE-Coil": 4.0, "PRO-Coil/ASN-Coil": 4.0, "PRO-Coil/GLY-Coil": 4.0, "PRO-Coil/HIS-Coil": 4.0, "PRO-Coil/LEU-Coil": 2.0, "PRO-Coil/ARG-Coil": 3.0, "PRO-Coil/TRP-Coil": 4.0, "PRO-Coil/ALA-Coil": 4.0, "PRO-Coil/VAL-Coil": 3.0, "PRO-Coil/GLU-Coil": 3.0, "PRO-Coil/TYR-Coil": 5.0, "PRO-Coil/MET-Coil": 3.0, "PRO-Coil/CYS-Sheet": 2.0, "PRO-Coil/ASP-Sheet": 2.0, "PRO-Coil/SER-Sheet": 2.0, "PRO-Coil/GLN-Sheet": 2.0, "PRO-Coil/LYS-Sheet": 1.0, "PRO-Coil/ILE-Sheet": 1.0, "PRO-Coil/PRO-Sheet": 2.0, "PRO-Coil/THR-Sheet": 1.0, "PRO-Coil/PHE-Sheet": 2.0, "PRO-Coil/ASN-Sheet": 2.0, "PRO-Coil/GLY-Sheet": 1.0, "PRO-Coil/HIS-Sheet": 2.0, "PRO-Coil/LEU-Sheet": 2.0, "PRO-Coil/ARG-Sheet": 2.0, "PRO-Coil/TRP-Sheet": 3.0, "PRO-Coil/ALA-Sheet": 3.0, "PRO-Coil/VAL-Sheet": 1.0, "PRO-Coil/GLU-Sheet": 1.0, "PRO-Coil/TYR-Sheet": 4.0, "PRO-Coil/MET-Sheet": 2.0, "THR-Coil/CYS-Helix": -Infinity, "THR-Coil/ASP-Helix": 2.0, "THR-Coil/SER-Helix": 2.0, "THR-Coil/GLN-Helix": 1.0, "THR-Coil/LYS-Helix": 1.0, "THR-Coil/ILE-Helix": -2.0, "THR-Coil/PRO-Helix": 2.0, "THR-Coil/THR-Helix": 2.0, "THR-Coil/PHE-Helix": -3.0, "THR-Coil/ASN-Helix": 1.0, "THR-Coil/GLY-Helix": 2.0, "THR-Coil/HIS-Helix": 2.0, "THR-Coil/LEU-Helix": -1.0, "THR-Coil/ARG-Helix": 2.0, "THR-Coil/TRP-Helix": -1.0, "THR-Coil/ALA-Helix": 1.0, "THR-Coil/VAL-Helix": 1.0, "THR-Coil/GLU-Helix": 3.0, "THR-Coil/TYR-Helix": -2.0, "THR-Coil/MET-Helix": -1.0, "THR-Coil/CYS-Coil": 3.0, "THR-Coil/ASP-Coil": 4.0, "THR-Coil/SER-Coil": 4.0, "THR-Coil/GLN-Coil": 3.0, "THR-Coil/LYS-Coil": 1.0, "THR-Coil/ILE-Coil": 2.0, "THR-Coil/PRO-Coil": 3.0, "THR-Coil/THR-Coil": 2.0, "THR-Coil/PHE-Coil": 3.0, "THR-Coil/ASN-Coil": 4.0, "THR-Coil/GLY-Coil": 4.0, "THR-Coil/HIS-Coil": 3.0, "THR-Coil/LEU-Coil": 3.0, "THR-Coil/ARG-Coil": 3.0, "THR-Coil/TRP-Coil": 1.0, "THR-Coil/ALA-Coil": 4.0, "THR-Coil/VAL-Coil": 2.0, "THR-Coil/GLU-Coil": 4.0, "THR-Coil/TYR-Coil": 2.0, "THR-Coil/MET-Coil": 2.0, "THR-Coil/CYS-Sheet": 3.0, "THR-Coil/ASP-Sheet": 4.0, "THR-Coil/SER-Sheet": 2.0, "THR-Coil/GLN-Sheet": -0.0, "THR-Coil/LYS-Sheet": 3.0, "THR-Coil/ILE-Sheet": 0.0, "THR-Coil/PRO-Sheet": 2.0, "THR-Coil/THR-Sheet": 1.0, "THR-Coil/PHE-Sheet": -0.0, "THR-Coil/ASN-Sheet": 3.0, "THR-Coil/GLY-Sheet": 3.0, "THR-Coil/HIS-Sheet": 4.0, "THR-Coil/LEU-Sheet": -2.0, "THR-Coil/ARG-Sheet": 2.0, "THR-Coil/TRP-Sheet": -Infinity, "THR-Coil/ALA-Sheet": 1.0, "THR-Coil/VAL-Sheet": -1.0, "THR-Coil/GLU-Sheet": 2.0, "THR-Coil/TYR-Sheet": 1.0, "THR-Coil/MET-Sheet": -1.0, "PHE-Coil/CYS-Helix": -Infinity, "PHE-Coil/ASP-Helix": 2.0, "PHE-Coil/SER-Helix": 1.0, "PHE-Coil/GLN-Helix": 1.0, "PHE-Coil/LYS-Helix": 2.0, "PHE-Coil/ILE-Helix": -1.0, "PHE-Coil/PRO-Helix": 2.0, "PHE-Coil/THR-Helix": -0.0, "PHE-Coil/PHE-Helix": -0.0, "PHE-Coil/ASN-Helix": 1.0, "PHE-Coil/GLY-Helix": 1.0, "PHE-Coil/HIS-Helix": 3.0, "PHE-Coil/LEU-Helix": 2.0, "PHE-Coil/ARG-Helix": 2.0, "PHE-Coil/TRP-Helix": -0.0, "PHE-Coil/ALA-Helix": 2.0, "PHE-Coil/VAL-Helix": 0.0, "PHE-Coil/GLU-Helix": 0.0, "PHE-Coil/TYR-Helix": 1.0, "PHE-Coil/MET-Helix": 1.0, "PHE-Coil/CYS-Coil": 3.0, "PHE-Coil/ASP-Coil": 2.0, "PHE-Coil/SER-Coil": 3.0, "PHE-Coil/GLN-Coil": 2.0, "PHE-Coil/LYS-Coil": 2.0, "PHE-Coil/ILE-Coil": 3.0, "PHE-Coil/PRO-Coil": 4.0, "PHE-Coil/THR-Coil": 3.0, "PHE-Coil/PHE-Coil": 2.0, "PHE-Coil/ASN-Coil": 1.0, "PHE-Coil/GLY-Coil": 3.0, "PHE-Coil/HIS-Coil": 2.0, "PHE-Coil/LEU-Coil": 3.0, "PHE-Coil/ARG-Coil": 3.0, "PHE-Coil/TRP-Coil": 3.0, "PHE-Coil/ALA-Coil": 3.0, "PHE-Coil/VAL-Coil": 1.0, "PHE-Coil/GLU-Coil": 3.0, "PHE-Coil/TYR-Coil": 1.0, "PHE-Coil/MET-Coil": 4.0, "PHE-Coil/CYS-Sheet": 3.0, "PHE-Coil/ASP-Sheet": -Infinity, "PHE-Coil/SER-Sheet": 3.0, "PHE-Coil/GLN-Sheet": 1.0, "PHE-Coil/LYS-Sheet": 3.0, "PHE-Coil/ILE-Sheet": -1.0, "PHE-Coil/PRO-Sheet": 4.0, "PHE-Coil/THR-Sheet": 0.0, "PHE-Coil/PHE-Sheet": 2.0, "PHE-Coil/ASN-Sheet": -0.0, "PHE-Coil/GLY-Sheet": 2.0, "PHE-Coil/HIS-Sheet": 0.0, "PHE-Coil/LEU-Sheet": 2.0, "PHE-Coil/ARG-Sheet": 2.0, "PHE-Coil/TRP-Sheet": -Infinity, "PHE-Coil/ALA-Sheet": 3.0, "PHE-Coil/VAL-Sheet": 2.0, "PHE-Coil/GLU-Sheet": 3.0, "PHE-Coil/TYR-Sheet": 4.0, "PHE-Coil/MET-Sheet": 2.0, "ASN-Coil/CYS-Helix": 0.0, "ASN-Coil/ASP-Helix": 2.0, "ASN-Coil/SER-Helix": 1.0, "ASN-Coil/GLN-Helix": 1.0, "ASN-Coil/LYS-Helix": 2.0, "ASN-Coil/ILE-Helix": -1.0, "ASN-Coil/PRO-Helix": -0.0, "ASN-Coil/THR-Helix": 0.0, "ASN-Coil/PHE-Helix": -4.0, "ASN-Coil/ASN-Helix": 2.0, "ASN-Coil/GLY-Helix": 0.0, "ASN-Coil/HIS-Helix": 2.0, "ASN-Coil/LEU-Helix": -3.0, "ASN-Coil/ARG-Helix": 0.0, "ASN-Coil/TRP-Helix": 1.0, "ASN-Coil/ALA-Helix": 0.0, "ASN-Coil/VAL-Helix": 1.0, "ASN-Coil/GLU-Helix": 2.0, "ASN-Coil/TYR-Helix": 2.0, "ASN-Coil/MET-Helix": -2.0, "ASN-Coil/CYS-Coil": 3.0, "ASN-Coil/ASP-Coil": 4.0, "ASN-Coil/SER-Coil": 4.0, "ASN-Coil/GLN-Coil": 4.0, "ASN-Coil/LYS-Coil": 4.0, "ASN-Coil/ILE-Coil": 3.0, "ASN-Coil/PRO-Coil": 4.0, "ASN-Coil/THR-Coil": 4.0, "ASN-Coil/PHE-Coil": 1.0, "ASN-Coil/ASN-Coil": 2.0, "ASN-Coil/GLY-Coil": 4.0, "ASN-Coil/HIS-Coil": 5.0, "ASN-Coil/LEU-Coil": 3.0, "ASN-Coil/ARG-Coil": 4.0, "ASN-Coil/TRP-Coil": 2.0, "ASN-Coil/ALA-Coil": 4.0, "ASN-Coil/VAL-Coil": 3.0, "ASN-Coil/GLU-Coil": 5.0, "ASN-Coil/TYR-Coil": 3.0, "ASN-Coil/MET-Coil": 1.0, "ASN-Coil/CYS-Sheet": -Infinity, "ASN-Coil/ASP-Sheet": 2.0, "ASN-Coil/SER-Sheet": 1.0, "ASN-Coil/GLN-Sheet": 3.0, "ASN-Coil/LYS-Sheet": 2.0, "ASN-Coil/ILE-Sheet": 1.0, "ASN-Coil/PRO-Sheet": -Infinity, "ASN-Coil/THR-Sheet": 1.0, "ASN-Coil/PHE-Sheet": 1.0, "ASN-Coil/ASN-Sheet": 3.0, "ASN-Coil/GLY-Sheet": 1.0, "ASN-Coil/HIS-Sheet": 1.0, "ASN-Coil/LEU-Sheet": -1.0, "ASN-Coil/ARG-Sheet": 2.0, "ASN-Coil/TRP-Sheet": 0.0, "ASN-Coil/ALA-Sheet": 2.0, "ASN-Coil/VAL-Sheet": 2.0, "ASN-Coil/GLU-Sheet": 2.0, "ASN-Coil/TYR-Sheet": 0.0, "ASN-Coil/MET-Sheet": 3.0, "GLY-Coil/CYS-Helix": 3.0, "GLY-Coil/ASP-Helix": 1.0, "GLY-Coil/SER-Helix": 2.0, "GLY-Coil/GLN-Helix": 1.0, "GLY-Coil/LYS-Helix": 0.0, "GLY-Coil/ILE-Helix": -1.0, "GLY-Coil/PRO-Helix": 1.0, "GLY-Coil/THR-Helix": 2.0, "GLY-Coil/PHE-Helix": -1.0, "GLY-Coil/ASN-Helix": 2.0, "GLY-Coil/GLY-Helix": 2.0, "GLY-Coil/HIS-Helix": 1.0, "GLY-Coil/LEU-Helix": -0.0, "GLY-Coil/ARG-Helix": 1.0, "GLY-Coil/TRP-Helix": -0.0, "GLY-Coil/ALA-Helix": 1.0, "GLY-Coil/VAL-Helix": 1.0, "GLY-Coil/GLU-Helix": 1.0, "GLY-Coil/TYR-Helix": 1.0, "GLY-Coil/MET-Helix": 0.0, "GLY-Coil/CYS-Coil": 5.0, "GLY-Coil/ASP-Coil": 4.0, "GLY-Coil/SER-Coil": 5.0, "GLY-Coil/GLN-Coil": 3.0, "GLY-Coil/LYS-Coil": 2.0, "GLY-Coil/ILE-Coil": 2.0, "GLY-Coil/PRO-Coil": 4.0, "GLY-Coil/THR-Coil": 4.0, "GLY-Coil/PHE-Coil": 3.0, "GLY-Coil/ASN-Coil": 4.0, "GLY-Coil/GLY-Coil": 3.0, "GLY-Coil/HIS-Coil": 4.0, "GLY-Coil/LEU-Coil": 3.0, "GLY-Coil/ARG-Coil": 3.0, "GLY-Coil/TRP-Coil": 3.0, "GLY-Coil/ALA-Coil": 5.0, "GLY-Coil/VAL-Coil": 4.0, "GLY-Coil/GLU-Coil": 3.0, "GLY-Coil/TYR-Coil": 3.0, "GLY-Coil/MET-Coil": 2.0, "GLY-Coil/CYS-Sheet": 3.0, "GLY-Coil/ASP-Sheet": 4.0, "GLY-Coil/SER-Sheet": 3.0, "GLY-Coil/GLN-Sheet": 3.0, "GLY-Coil/LYS-Sheet": 2.0, "GLY-Coil/ILE-Sheet": 1.0, "GLY-Coil/PRO-Sheet": 4.0, "GLY-Coil/THR-Sheet": 2.0, "GLY-Coil/PHE-Sheet": 1.0, "GLY-Coil/ASN-Sheet": 3.0, "GLY-Coil/GLY-Sheet": 3.0, "GLY-Coil/HIS-Sheet": 1.0, "GLY-Coil/LEU-Sheet": 2.0, "GLY-Coil/ARG-Sheet": 2.0, "GLY-Coil/TRP-Sheet": 2.0, "GLY-Coil/ALA-Sheet": 2.0, "GLY-Coil/VAL-Sheet": 1.0, "GLY-Coil/GLU-Sheet": 2.0, "GLY-Coil/TYR-Sheet": 2.0, "GLY-Coil/MET-Sheet": 2.0, "HIS-Coil/CYS-Helix": 4.0, "HIS-Coil/ASP-Helix": 2.0, "HIS-Coil/SER-Helix": 1.0, "HIS-Coil/GLN-Helix": 2.0, "HIS-Coil/LYS-Helix": -1.0, "HIS-Coil/ILE-Helix": 0.0, "HIS-Coil/PRO-Helix": 3.0, "HIS-Coil/THR-Helix": 2.0, "HIS-Coil/PHE-Helix": -0.0, "HIS-Coil/ASN-Helix": 0.0, "HIS-Coil/GLY-Helix": 0.0, "HIS-Coil/HIS-Helix": 2.0, "HIS-Coil/LEU-Helix": 0.0, "HIS-Coil/ARG-Helix": 2.0, "HIS-Coil/TRP-Helix": 2.0, "HIS-Coil/ALA-Helix": 0.0, "HIS-Coil/VAL-Helix": -1.0, "HIS-Coil/GLU-Helix": 2.0, "HIS-Coil/TYR-Helix": -1.0, "HIS-Coil/MET-Helix": -Infinity, "HIS-Coil/CYS-Coil": 4.0, "HIS-Coil/ASP-Coil": 3.0, "HIS-Coil/SER-Coil": 4.0, "HIS-Coil/GLN-Coil": 4.0, "HIS-Coil/LYS-Coil": 2.0, "HIS-Coil/ILE-Coil": 2.0, "HIS-Coil/PRO-Coil": 4.0, "HIS-Coil/THR-Coil": 3.0, "HIS-Coil/PHE-Coil": 2.0, "HIS-Coil/ASN-Coil": 5.0, "HIS-Coil/GLY-Coil": 4.0, "HIS-Coil/HIS-Coil": 4.0, "HIS-Coil/LEU-Coil": 4.0, "HIS-Coil/ARG-Coil": 4.0, "HIS-Coil/TRP-Coil": 2.0, "HIS-Coil/ALA-Coil": 4.0, "HIS-Coil/VAL-Coil": 4.0, "HIS-Coil/GLU-Coil": 4.0, "HIS-Coil/TYR-Coil": 2.0, "HIS-Coil/MET-Coil": 4.0, "HIS-Coil/CYS-Sheet": -Infinity, "HIS-Coil/ASP-Sheet": 0.0, "HIS-Coil/SER-Sheet": 2.0, "HIS-Coil/GLN-Sheet": 3.0, "HIS-Coil/LYS-Sheet": 2.0, "HIS-Coil/ILE-Sheet": 1.0, "HIS-Coil/PRO-Sheet": 2.0, "HIS-Coil/THR-Sheet": 0.0, "HIS-Coil/PHE-Sheet": -1.0, "HIS-Coil/ASN-Sheet": -Infinity, "HIS-Coil/GLY-Sheet": 1.0, "HIS-Coil/HIS-Sheet": 1.0, "HIS-Coil/LEU-Sheet": 0.0, "HIS-Coil/ARG-Sheet": 3.0, "HIS-Coil/TRP-Sheet": 4.0, "HIS-Coil/ALA-Sheet": -Infinity, "HIS-Coil/VAL-Sheet": -Infinity, "HIS-Coil/GLU-Sheet": 1.0, "HIS-Coil/TYR-Sheet": 1.0, "HIS-Coil/MET-Sheet": -Infinity, "LEU-Coil/CYS-Helix": 3.0, "LEU-Coil/ASP-Helix": 2.0, "LEU-Coil/SER-Helix": 2.0, "LEU-Coil/GLN-Helix": 3.0, "LEU-Coil/LYS-Helix": 1.0, "LEU-Coil/ILE-Helix": 1.0, "LEU-Coil/PRO-Helix": 1.0, "LEU-Coil/THR-Helix": 1.0, "LEU-Coil/PHE-Helix": 1.0, "LEU-Coil/ASN-Helix": 1.0, "LEU-Coil/GLY-Helix": 1.0, "LEU-Coil/HIS-Helix": 2.0, "LEU-Coil/LEU-Helix": 1.0, "LEU-Coil/ARG-Helix": 1.0, "LEU-Coil/TRP-Helix": 1.0, "LEU-Coil/ALA-Helix": 2.0, "LEU-Coil/VAL-Helix": 1.0, "LEU-Coil/GLU-Helix": 2.0, "LEU-Coil/TYR-Helix": 2.0, "LEU-Coil/MET-Helix": 1.0, "LEU-Coil/CYS-Coil": 3.0, "LEU-Coil/ASP-Coil": 2.0, "LEU-Coil/SER-Coil": 3.0, "LEU-Coil/GLN-Coil": 3.0, "LEU-Coil/LYS-Coil": 1.0, "LEU-Coil/ILE-Coil": 2.0, "LEU-Coil/PRO-Coil": 2.0, "LEU-Coil/THR-Coil": 3.0, "LEU-Coil/PHE-Coil": 3.0, "LEU-Coil/ASN-Coil": 3.0, "LEU-Coil/GLY-Coil": 3.0, "LEU-Coil/HIS-Coil": 4.0, "LEU-Coil/LEU-Coil": 1.0, "LEU-Coil/ARG-Coil": 3.0, "LEU-Coil/TRP-Coil": 3.0, "LEU-Coil/ALA-Coil": 4.0, "LEU-Coil/VAL-Coil": 2.0, "LEU-Coil/GLU-Coil": 2.0, "LEU-Coil/TYR-Coil": 4.0, "LEU-Coil/MET-Coil": 3.0, "LEU-Coil/CYS-Sheet": 3.0, "LEU-Coil/ASP-Sheet": 1.0, "LEU-Coil/SER-Sheet": 2.0, "LEU-Coil/GLN-Sheet": 3.0, "LEU-Coil/LYS-Sheet": -0.0, "LEU-Coil/ILE-Sheet": 2.0, "LEU-Coil/PRO-Sheet": 2.0, "LEU-Coil/THR-Sheet": 0.0, "LEU-Coil/PHE-Sheet": 2.0, "LEU-Coil/ASN-Sheet": 2.0, "LEU-Coil/GLY-Sheet": 3.0, "LEU-Coil/HIS-Sheet": 3.0, "LEU-Coil/LEU-Sheet": 3.0, "LEU-Coil/ARG-Sheet": 2.0, "LEU-Coil/TRP-Sheet": 3.0, "LEU-Coil/ALA-Sheet": 4.0, "LEU-Coil/VAL-Sheet": 3.0, "LEU-Coil/GLU-Sheet": 1.0, "LEU-Coil/TYR-Sheet": 3.0, "LEU-Coil/MET-Sheet": 2.0, "ARG-Coil/CYS-Helix": 1.0, "ARG-Coil/ASP-Helix": 3.0, "ARG-Coil/SER-Helix": 0.0, "ARG-Coil/GLN-Helix": 1.0, "ARG-Coil/LYS-Helix": -0.0, "ARG-Coil/ILE-Helix": -0.0, "ARG-Coil/PRO-Helix": 1.0, "ARG-Coil/THR-Helix": 2.0, "ARG-Coil/PHE-Helix": -0.0, "ARG-Coil/ASN-Helix": 1.0, "ARG-Coil/GLY-Helix": -0.0, "ARG-Coil/HIS-Helix": 2.0, "ARG-Coil/LEU-Helix": -0.0, "ARG-Coil/ARG-Helix": -0.0, "ARG-Coil/TRP-Helix": -1.0, "ARG-Coil/ALA-Helix": -0.0, "ARG-Coil/VAL-Helix": -2.0, "ARG-Coil/GLU-Helix": 3.0, "ARG-Coil/TYR-Helix": 1.0, "ARG-Coil/MET-Helix": -Infinity, "ARG-Coil/CYS-Coil": 2.0, "ARG-Coil/ASP-Coil": 5.0, "ARG-Coil/SER-Coil": 3.0, "ARG-Coil/GLN-Coil": 4.0, "ARG-Coil/LYS-Coil": 2.0, "ARG-Coil/ILE-Coil": 3.0, "ARG-Coil/PRO-Coil": 3.0, "ARG-Coil/THR-Coil": 3.0, "ARG-Coil/PHE-Coil": 3.0, "ARG-Coil/ASN-Coil": 4.0, "ARG-Coil/GLY-Coil": 3.0, "ARG-Coil/HIS-Coil": 4.0, "ARG-Coil/LEU-Coil": 3.0, "ARG-Coil/ARG-Coil": 3.0, "ARG-Coil/TRP-Coil": 4.0, "ARG-Coil/ALA-Coil": 3.0, "ARG-Coil/VAL-Coil": 2.0, "ARG-Coil/GLU-Coil": 4.0, "ARG-Coil/TYR-Coil": 4.0, "ARG-Coil/MET-Coil": 3.0, "ARG-Coil/CYS-Sheet": -Infinity, "ARG-Coil/ASP-Sheet": 2.0, "ARG-Coil/SER-Sheet": -0.0, "ARG-Coil/GLN-Sheet": 2.0, "ARG-Coil/LYS-Sheet": 2.0, "ARG-Coil/ILE-Sheet": 0.0, "ARG-Coil/PRO-Sheet": 2.0, "ARG-Coil/THR-Sheet": -2.0, "ARG-Coil/PHE-Sheet": -2.0, "ARG-Coil/ASN-Sheet": 4.0, "ARG-Coil/GLY-Sheet": -2.0, "ARG-Coil/HIS-Sheet": 1.0, "ARG-Coil/LEU-Sheet": 1.0, "ARG-Coil/ARG-Sheet": 3.0, "ARG-Coil/TRP-Sheet": 4.0, "ARG-Coil/ALA-Sheet": -1.0, "ARG-Coil/VAL-Sheet": 2.0, "ARG-Coil/GLU-Sheet": 4.0, "ARG-Coil/TYR-Sheet": 1.0, "ARG-Coil/MET-Sheet": 2.0, "TRP-Coil/CYS-Helix": -Infinity, "TRP-Coil/ASP-Helix": 1.0, "TRP-Coil/SER-Helix": 1.0, "TRP-Coil/GLN-Helix": 3.0, "TRP-Coil/LYS-Helix": -1.0, "TRP-Coil/ILE-Helix": 2.0, "TRP-Coil/PRO-Helix": 4.0, "TRP-Coil/THR-Helix": -2.0, "TRP-Coil/PHE-Helix": -1.0, "TRP-Coil/ASN-Helix": 1.0, "TRP-Coil/GLY-Helix": -1.0, "TRP-Coil/HIS-Helix": -Infinity, "TRP-Coil/LEU-Helix": 1.0, "TRP-Coil/ARG-Helix": 2.0, "TRP-Coil/TRP-Helix": 1.0, "TRP-Coil/ALA-Helix": 2.0, "TRP-Coil/VAL-Helix": 2.0, "TRP-Coil/GLU-Helix": 1.0, "TRP-Coil/TYR-Helix": -Infinity, "TRP-Coil/MET-Helix": 3.0, "TRP-Coil/CYS-Coil": 4.0, "TRP-Coil/ASP-Coil": 2.0, "TRP-Coil/SER-Coil": 2.0, "TRP-Coil/GLN-Coil": 4.0, "TRP-Coil/LYS-Coil": 4.0, "TRP-Coil/ILE-Coil": 5.0, "TRP-Coil/PRO-Coil": 4.0, "TRP-Coil/THR-Coil": 1.0, "TRP-Coil/PHE-Coil": 3.0, "TRP-Coil/ASN-Coil": 2.0, "TRP-Coil/GLY-Coil": 3.0, "TRP-Coil/HIS-Coil": 2.0, "TRP-Coil/LEU-Coil": 3.0, "TRP-Coil/ARG-Coil": 4.0, "TRP-Coil/TRP-Coil": 5.0, "TRP-Coil/ALA-Coil": -1.0, "TRP-Coil/VAL-Coil": 1.0, "TRP-Coil/GLU-Coil": 1.0, "TRP-Coil/TYR-Coil": 2.0, "TRP-Coil/MET-Coil": -Infinity, "TRP-Coil/CYS-Sheet": -Infinity, "TRP-Coil/ASP-Sheet": -Infinity, "TRP-Coil/SER-Sheet": 2.0, "TRP-Coil/GLN-Sheet": 2.0, "TRP-Coil/LYS-Sheet": 2.0, "TRP-Coil/ILE-Sheet": 1.0, "TRP-Coil/PRO-Sheet": 3.0, "TRP-Coil/THR-Sheet": 2.0, "TRP-Coil/PHE-Sheet": -Infinity, "TRP-Coil/ASN-Sheet": 4.0, "TRP-Coil/GLY-Sheet": -Infinity, "TRP-Coil/HIS-Sheet": 3.0, "TRP-Coil/LEU-Sheet": 3.0, "TRP-Coil/ARG-Sheet": 4.0, "TRP-Coil/TRP-Sheet": -Infinity, "TRP-Coil/ALA-Sheet": -0.0, "TRP-Coil/VAL-Sheet": 2.0, "TRP-Coil/GLU-Sheet": 1.0, "TRP-Coil/TYR-Sheet": -Infinity, "TRP-Coil/MET-Sheet": 3.0, "ALA-Coil/CYS-Helix": 2.0, "ALA-Coil/ASP-Helix": 1.0, "ALA-Coil/SER-Helix": 2.0, "ALA-Coil/GLN-Helix": 2.0, "ALA-Coil/LYS-Helix": 0.0, "ALA-Coil/ILE-Helix": 1.0, "ALA-Coil/PRO-Helix": 2.0, "ALA-Coil/THR-Helix": 2.0, "ALA-Coil/PHE-Helix": 1.0, "ALA-Coil/ASN-Helix": 2.0, "ALA-Coil/GLY-Helix": 3.0, "ALA-Coil/HIS-Helix": 2.0, "ALA-Coil/LEU-Helix": 2.0, "ALA-Coil/ARG-Helix": -1.0, "ALA-Coil/TRP-Helix": 0.0, "ALA-Coil/ALA-Helix": 2.0, "ALA-Coil/VAL-Helix": 2.0, "ALA-Coil/GLU-Helix": 1.0, "ALA-Coil/TYR-Helix": 1.0, "ALA-Coil/MET-Helix": 3.0, "ALA-Coil/CYS-Coil": 3.0, "ALA-Coil/ASP-Coil": 4.0, "ALA-Coil/SER-Coil": 4.0, "ALA-Coil/GLN-Coil": 3.0, "ALA-Coil/LYS-Coil": 3.0, "ALA-Coil/ILE-Coil": 3.0, "ALA-Coil/PRO-Coil": 4.0, "ALA-Coil/THR-Coil": 4.0, "ALA-Coil/PHE-Coil": 3.0, "ALA-Coil/ASN-Coil": 4.0, "ALA-Coil/GLY-Coil": 5.0, "ALA-Coil/HIS-Coil": 4.0, "ALA-Coil/LEU-Coil": 4.0, "ALA-Coil/ARG-Coil": 3.0, "ALA-Coil/TRP-Coil": -1.0, "ALA-Coil/ALA-Coil": 2.0, "ALA-Coil/VAL-Coil": 3.0, "ALA-Coil/GLU-Coil": 2.0, "ALA-Coil/TYR-Coil": 3.0, "ALA-Coil/MET-Coil": 3.0, "ALA-Coil/CYS-Sheet": 3.0, "ALA-Coil/ASP-Sheet": 3.0, "ALA-Coil/SER-Sheet": 1.0, "ALA-Coil/GLN-Sheet": 2.0, "ALA-Coil/LYS-Sheet": 2.0, "ALA-Coil/ILE-Sheet": 2.0, "ALA-Coil/PRO-Sheet": 4.0, "ALA-Coil/THR-Sheet": 2.0, "ALA-Coil/PHE-Sheet": 1.0, "ALA-Coil/ASN-Sheet": 5.0, "ALA-Coil/GLY-Sheet": 2.0, "ALA-Coil/HIS-Sheet": 2.0, "ALA-Coil/LEU-Sheet": 2.0, "ALA-Coil/ARG-Sheet": 1.0, "ALA-Coil/TRP-Sheet": 3.0, "ALA-Coil/ALA-Sheet": 3.0, "ALA-Coil/VAL-Sheet": 2.0, "ALA-Coil/GLU-Sheet": 1.0, "ALA-Coil/TYR-Sheet": 2.0, "ALA-Coil/MET-Sheet": 3.0, "VAL-Coil/CYS-Helix": 1.0, "VAL-Coil/ASP-Helix": 1.0, "VAL-Coil/SER-Helix": 1.0, "VAL-Coil/GLN-Helix": 0.0, "VAL-Coil/LYS-Helix": 0.0, "VAL-Coil/ILE-Helix": 1.0, "VAL-Coil/PRO-Helix": 2.0, "VAL-Coil/THR-Helix": 2.0, "VAL-Coil/PHE-Helix": 0.0, "VAL-Coil/ASN-Helix": 2.0, "VAL-Coil/GLY-Helix": 2.0, "VAL-Coil/HIS-Helix": 1.0, "VAL-Coil/LEU-Helix": 1.0, "VAL-Coil/ARG-Helix": 1.0, "VAL-Coil/TRP-Helix": 2.0, "VAL-Coil/ALA-Helix": 2.0, "VAL-Coil/VAL-Helix": 1.0, "VAL-Coil/GLU-Helix": 1.0, "VAL-Coil/TYR-Helix": 2.0, "VAL-Coil/MET-Helix": 3.0, "VAL-Coil/CYS-Coil": 4.0, "VAL-Coil/ASP-Coil": 3.0, "VAL-Coil/SER-Coil": 3.0, "VAL-Coil/GLN-Coil": 3.0, "VAL-Coil/LYS-Coil": 2.0, "VAL-Coil/ILE-Coil": 3.0, "VAL-Coil/PRO-Coil": 3.0, "VAL-Coil/THR-Coil": 2.0, "VAL-Coil/PHE-Coil": 1.0, "VAL-Coil/ASN-Coil": 3.0, "VAL-Coil/GLY-Coil": 4.0, "VAL-Coil/HIS-Coil": 4.0, "VAL-Coil/LEU-Coil": 2.0, "VAL-Coil/ARG-Coil": 2.0, "VAL-Coil/TRP-Coil": 1.0, "VAL-Coil/ALA-Coil": 3.0, "VAL-Coil/VAL-Coil": 0.0, "VAL-Coil/GLU-Coil": 2.0, "VAL-Coil/TYR-Coil": 2.0, "VAL-Coil/MET-Coil": 2.0, "VAL-Coil/CYS-Sheet": 0.0, "VAL-Coil/ASP-Sheet": -Infinity, "VAL-Coil/SER-Sheet": 3.0, "VAL-Coil/GLN-Sheet": 4.0, "VAL-Coil/LYS-Sheet": 2.0, "VAL-Coil/ILE-Sheet": 2.0, "VAL-Coil/PRO-Sheet": 2.0, "VAL-Coil/THR-Sheet": 1.0, "VAL-Coil/PHE-Sheet": 3.0, "VAL-Coil/ASN-Sheet": 3.0, "VAL-Coil/GLY-Sheet": 4.0, "VAL-Coil/HIS-Sheet": 2.0, "VAL-Coil/LEU-Sheet": 1.0, "VAL-Coil/ARG-Sheet": 3.0, "VAL-Coil/TRP-Sheet": 3.0, "VAL-Coil/ALA-Sheet": 3.0, "VAL-Coil/VAL-Sheet": 2.0, "VAL-Coil/GLU-Sheet": 1.0, "VAL-Coil/TYR-Sheet": 3.0, "VAL-Coil/MET-Sheet": -Infinity, "GLU-Coil/CYS-Helix": 2.0, "GLU-Coil/ASP-Helix": 1.0, "GLU-Coil/SER-Helix": 2.0, "GLU-Coil/GLN-Helix": 0.0, "GLU-Coil/LYS-Helix": 3.0, "GLU-Coil/ILE-Helix": -1.0, "GLU-Coil/PRO-Helix": -0.0, "GLU-Coil/THR-Helix": 1.0, "GLU-Coil/PHE-Helix": -2.0, "GLU-Coil/ASN-Helix": 3.0, "GLU-Coil/GLY-Helix": 1.0, "GLU-Coil/HIS-Helix": 2.0, "GLU-Coil/LEU-Helix": -1.0, "GLU-Coil/ARG-Helix": 3.0, "GLU-Coil/TRP-Helix": -2.0, "GLU-Coil/ALA-Helix": 1.0, "GLU-Coil/VAL-Helix": -2.0, "GLU-Coil/GLU-Helix": 0.0, "GLU-Coil/TYR-Helix": -2.0, "GLU-Coil/MET-Helix": 0.0, "GLU-Coil/CYS-Coil": 3.0, "GLU-Coil/ASP-Coil": 2.0, "GLU-Coil/SER-Coil": 3.0, "GLU-Coil/GLN-Coil": 4.0, "GLU-Coil/LYS-Coil": 5.0, "GLU-Coil/ILE-Coil": 1.0, "GLU-Coil/PRO-Coil": 3.0, "GLU-Coil/THR-Coil": 4.0, "GLU-Coil/PHE-Coil": 3.0, "GLU-Coil/ASN-Coil": 5.0, "GLU-Coil/GLY-Coil": 3.0, "GLU-Coil/HIS-Coil": 4.0, "GLU-Coil/LEU-Coil": 2.0, "GLU-Coil/ARG-Coil": 4.0, "GLU-Coil/TRP-Coil": 1.0, "GLU-Coil/ALA-Coil": 2.0, "GLU-Coil/VAL-Coil": 2.0, "GLU-Coil/GLU-Coil": 2.0, "GLU-Coil/TYR-Coil": 3.0, "GLU-Coil/MET-Coil": 3.0, "GLU-Coil/CYS-Sheet": -Infinity, "GLU-Coil/ASP-Sheet": -1.0, "GLU-Coil/SER-Sheet": 2.0, "GLU-Coil/GLN-Sheet": 1.0, "GLU-Coil/LYS-Sheet": 4.0, "GLU-Coil/ILE-Sheet": -0.0, "GLU-Coil/PRO-Sheet": -Infinity, "GLU-Coil/THR-Sheet": 2.0, "GLU-Coil/PHE-Sheet": 0.0, "GLU-Coil/ASN-Sheet": 3.0, "GLU-Coil/GLY-Sheet": -2.0, "GLU-Coil/HIS-Sheet": 4.0, "GLU-Coil/LEU-Sheet": 1.0, "GLU-Coil/ARG-Sheet": 2.0, "GLU-Coil/TRP-Sheet": -Infinity, "GLU-Coil/ALA-Sheet": -2.0, "GLU-Coil/VAL-Sheet": 1.0, "GLU-Coil/GLU-Sheet": 1.0, "GLU-Coil/TYR-Sheet": 1.0, "GLU-Coil/MET-Sheet": 3.0, "TYR-Coil/CYS-Helix": 1.0, "TYR-Coil/ASP-Helix": -0.0, "TYR-Coil/SER-Helix": 2.0, "TYR-Coil/GLN-Helix": 2.0, "TYR-Coil/LYS-Helix": 2.0, "TYR-Coil/ILE-Helix": 1.0, "TYR-Coil/PRO-Helix": 1.0, "TYR-Coil/THR-Helix": 1.0, "TYR-Coil/PHE-Helix": -0.0, "TYR-Coil/ASN-Helix": 2.0, "TYR-Coil/GLY-Helix": -1.0, "TYR-Coil/HIS-Helix": 2.0, "TYR-Coil/LEU-Helix": 1.0, "TYR-Coil/ARG-Helix": 2.0, "TYR-Coil/TRP-Helix": 1.0, "TYR-Coil/ALA-Helix": 2.0, "TYR-Coil/VAL-Helix": 0.0, "TYR-Coil/GLU-Helix": 0.0, "TYR-Coil/TYR-Helix": 0.0, "TYR-Coil/MET-Helix": 2.0, "TYR-Coil/CYS-Coil": 1.0, "TYR-Coil/ASP-Coil": 2.0, "TYR-Coil/SER-Coil": 1.0, "TYR-Coil/GLN-Coil": 4.0, "TYR-Coil/LYS-Coil": 4.0, "TYR-Coil/ILE-Coil": 3.0, "TYR-Coil/PRO-Coil": 5.0, "TYR-Coil/THR-Coil": 2.0, "TYR-Coil/PHE-Coil": 1.0, "TYR-Coil/ASN-Coil": 3.0, "TYR-Coil/GLY-Coil": 3.0, "TYR-Coil/HIS-Coil": 2.0, "TYR-Coil/LEU-Coil": 4.0, "TYR-Coil/ARG-Coil": 4.0, "TYR-Coil/TRP-Coil": 2.0, "TYR-Coil/ALA-Coil": 3.0, "TYR-Coil/VAL-Coil": 2.0, "TYR-Coil/GLU-Coil": 3.0, "TYR-Coil/TYR-Coil": 2.0, "TYR-Coil/MET-Coil": 1.0, "TYR-Coil/CYS-Sheet": -Infinity, "TYR-Coil/ASP-Sheet": -0.0, "TYR-Coil/SER-Sheet": -2.0, "TYR-Coil/GLN-Sheet": -0.0, "TYR-Coil/LYS-Sheet": 2.0, "TYR-Coil/ILE-Sheet": 0.0, "TYR-Coil/PRO-Sheet": 1.0, "TYR-Coil/THR-Sheet": -0.0, "TYR-Coil/PHE-Sheet": -Infinity, "TYR-Coil/ASN-Sheet": -Infinity, "TYR-Coil/GLY-Sheet": 3.0, "TYR-Coil/HIS-Sheet": -Infinity, "TYR-Coil/LEU-Sheet": 1.0, "TYR-Coil/ARG-Sheet": 2.0, "TYR-Coil/TRP-Sheet": -Infinity, "TYR-Coil/ALA-Sheet": 2.0, "TYR-Coil/VAL-Sheet": 2.0, "TYR-Coil/GLU-Sheet": -1.0, "TYR-Coil/TYR-Sheet": 1.0, "TYR-Coil/MET-Sheet": 0.0, "MET-Coil/CYS-Helix": -Infinity, "MET-Coil/ASP-Helix": -1.0, "MET-Coil/SER-Helix": 3.0, "MET-Coil/GLN-Helix": 2.0, "MET-Coil/LYS-Helix": 2.0, "MET-Coil/ILE-Helix": 2.0, "MET-Coil/PRO-Helix": 3.0, "MET-Coil/THR-Helix": 3.0, "MET-Coil/PHE-Helix": 3.0, "MET-Coil/ASN-Helix": 3.0, "MET-Coil/GLY-Helix": 1.0, "MET-Coil/HIS-Helix": 2.0, "MET-Coil/LEU-Helix": 3.0, "MET-Coil/ARG-Helix": 2.0, "MET-Coil/TRP-Helix": 3.0, "MET-Coil/ALA-Helix": 2.0, "MET-Coil/VAL-Helix": 2.0, "MET-Coil/GLU-Helix": 1.0, "MET-Coil/TYR-Helix": 2.0, "MET-Coil/MET-Helix": 2.0, "MET-Coil/CYS-Coil": 4.0, "MET-Coil/ASP-Coil": 1.0, "MET-Coil/SER-Coil": 2.0, "MET-Coil/GLN-Coil": 3.0, "MET-Coil/LYS-Coil": 4.0, "MET-Coil/ILE-Coil": -Infinity, "MET-Coil/PRO-Coil": 3.0, "MET-Coil/THR-Coil": 2.0, "MET-Coil/PHE-Coil": 4.0, "MET-Coil/ASN-Coil": 1.0, "MET-Coil/GLY-Coil": 2.0, "MET-Coil/HIS-Coil": 4.0, "MET-Coil/LEU-Coil": 3.0, "MET-Coil/ARG-Coil": 3.0, "MET-Coil/TRP-Coil": -Infinity, "MET-Coil/ALA-Coil": 3.0, "MET-Coil/VAL-Coil": 2.0, "MET-Coil/GLU-Coil": 3.0, "MET-Coil/TYR-Coil": 1.0, "MET-Coil/MET-Coil": 3.0, "MET-Coil/CYS-Sheet": 2.0, "MET-Coil/ASP-Sheet": 3.0, "MET-Coil/SER-Sheet": 3.0, "MET-Coil/GLN-Sheet": 3.0, "MET-Coil/LYS-Sheet": 3.0, "MET-Coil/ILE-Sheet": 4.0, "MET-Coil/PRO-Sheet": 5.0, "MET-Coil/THR-Sheet": 2.0, "MET-Coil/PHE-Sheet": -Infinity, "MET-Coil/ASN-Sheet": 3.0, "MET-Coil/GLY-Sheet": 3.0, "MET-Coil/HIS-Sheet": 4.0, "MET-Coil/LEU-Sheet": 3.0, "MET-Coil/ARG-Sheet": 3.0, "MET-Coil/TRP-Sheet": 3.0, "MET-Coil/ALA-Sheet": 2.0, "MET-Coil/VAL-Sheet": -1.0, "MET-Coil/GLU-Sheet": 1.0, "MET-Coil/TYR-Sheet": 4.0, "MET-Coil/MET-Sheet": 2.0, "CYS-Sheet/CYS-Helix": 5.0, "CYS-Sheet/ASP-Helix": -1.0, "CYS-Sheet/SER-Helix": -1.0, "CYS-Sheet/GLN-Helix": 1.0, "CYS-Sheet/LYS-Helix": -Infinity, "CYS-Sheet/ILE-Helix": 0.0, "CYS-Sheet/PRO-Helix": 2.0, "CYS-Sheet/THR-Helix": 2.0, "CYS-Sheet/PHE-Helix": 2.0, "CYS-Sheet/ASN-Helix": -Infinity, "CYS-Sheet/GLY-Helix": -0.0, "CYS-Sheet/HIS-Helix": -Infinity, "CYS-Sheet/LEU-Helix": 2.0, "CYS-Sheet/ARG-Helix": -Infinity, "CYS-Sheet/TRP-Helix": -Infinity, "CYS-Sheet/ALA-Helix": 0.0, "CYS-Sheet/VAL-Helix": 0.0, "CYS-Sheet/GLU-Helix": -Infinity, "CYS-Sheet/TYR-Helix": 3.0, "CYS-Sheet/MET-Helix": 1.0, "CYS-Sheet/CYS-Coil": 7.0, "CYS-Sheet/ASP-Coil": 2.0, "CYS-Sheet/SER-Coil": 2.0, "CYS-Sheet/GLN-Coil": -Infinity, "CYS-Sheet/LYS-Coil": -Infinity, "CYS-Sheet/ILE-Coil": 3.0, "CYS-Sheet/PRO-Coil": 2.0, "CYS-Sheet/THR-Coil": 3.0, "CYS-Sheet/PHE-Coil": 3.0, "CYS-Sheet/ASN-Coil": -Infinity, "CYS-Sheet/GLY-Coil": 3.0, "CYS-Sheet/HIS-Coil": -Infinity, "CYS-Sheet/LEU-Coil": 3.0, "CYS-Sheet/ARG-Coil": -Infinity, "CYS-Sheet/TRP-Coil": -Infinity, "CYS-Sheet/ALA-Coil": 3.0, "CYS-Sheet/VAL-Coil": 0.0, "CYS-Sheet/GLU-Coil": -Infinity, "CYS-Sheet/TYR-Coil": -Infinity, "CYS-Sheet/MET-Coil": 2.0, "CYS-Sheet/CYS-Sheet": 8.0, "CYS-Sheet/ASP-Sheet": -Infinity, "CYS-Sheet/SER-Sheet": 0.0, "CYS-Sheet/GLN-Sheet": 5.0, "CYS-Sheet/LYS-Sheet": 4.0, "CYS-Sheet/ILE-Sheet": 6.0, "CYS-Sheet/PRO-Sheet": 3.0, "CYS-Sheet/THR-Sheet": 4.0, "CYS-Sheet/PHE-Sheet": 4.0, "CYS-Sheet/ASN-Sheet": 3.0, "CYS-Sheet/GLY-Sheet": 6.0, "CYS-Sheet/HIS-Sheet": 5.0, "CYS-Sheet/LEU-Sheet": 6.0, "CYS-Sheet/ARG-Sheet": 1.0, "CYS-Sheet/TRP-Sheet": 4.0, "CYS-Sheet/ALA-Sheet": 5.0, "CYS-Sheet/VAL-Sheet": 6.0, "CYS-Sheet/GLU-Sheet": 4.0, "CYS-Sheet/TYR-Sheet": 5.0, "CYS-Sheet/MET-Sheet": 6.0, "ASP-Sheet/CYS-Helix": -Infinity, "ASP-Sheet/ASP-Helix": -0.0, "ASP-Sheet/SER-Helix": 0.0, "ASP-Sheet/GLN-Helix": -1.0, "ASP-Sheet/LYS-Helix": 1.0, "ASP-Sheet/ILE-Helix": 1.0, "ASP-Sheet/PRO-Helix": -1.0, "ASP-Sheet/THR-Helix": -1.0, "ASP-Sheet/PHE-Helix": -Infinity, "ASP-Sheet/ASN-Helix": 0.0, "ASP-Sheet/GLY-Helix": -1.0, "ASP-Sheet/HIS-Helix": -Infinity, "ASP-Sheet/LEU-Helix": -Infinity, "ASP-Sheet/ARG-Helix": 3.0, "ASP-Sheet/TRP-Helix": 1.0, "ASP-Sheet/ALA-Helix": -1.0, "ASP-Sheet/VAL-Helix": -1.0, "ASP-Sheet/GLU-Helix": -2.0, "ASP-Sheet/TYR-Helix": -0.0, "ASP-Sheet/MET-Helix": -Infinity, "ASP-Sheet/CYS-Coil": -Infinity, "ASP-Sheet/ASP-Coil": -0.0, "ASP-Sheet/SER-Coil": 4.0, "ASP-Sheet/GLN-Coil": 1.0, "ASP-Sheet/LYS-Coil": 4.0, "ASP-Sheet/ILE-Coil": -Infinity, "ASP-Sheet/PRO-Coil": 2.0, "ASP-Sheet/THR-Coil": 4.0, "ASP-Sheet/PHE-Coil": -Infinity, "ASP-Sheet/ASN-Coil": 2.0, "ASP-Sheet/GLY-Coil": 4.0, "ASP-Sheet/HIS-Coil": 0.0, "ASP-Sheet/LEU-Coil": 1.0, "ASP-Sheet/ARG-Coil": 2.0, "ASP-Sheet/TRP-Coil": -Infinity, "ASP-Sheet/ALA-Coil": 3.0, "ASP-Sheet/VAL-Coil": -Infinity, "ASP-Sheet/GLU-Coil": -1.0, "ASP-Sheet/TYR-Coil": -0.0, "ASP-Sheet/MET-Coil": 3.0, "ASP-Sheet/CYS-Sheet": -Infinity, "ASP-Sheet/ASP-Sheet": -Infinity, "ASP-Sheet/SER-Sheet": 6.0, "ASP-Sheet/GLN-Sheet": 5.0, "ASP-Sheet/LYS-Sheet": 5.0, "ASP-Sheet/ILE-Sheet": 2.0, "ASP-Sheet/PRO-Sheet": -Infinity, "ASP-Sheet/THR-Sheet": 5.0, "ASP-Sheet/PHE-Sheet": -0.0, "ASP-Sheet/ASN-Sheet": 6.0, "ASP-Sheet/GLY-Sheet": 4.0, "ASP-Sheet/HIS-Sheet": 5.0, "ASP-Sheet/LEU-Sheet": 3.0, "ASP-Sheet/ARG-Sheet": 7.0, "ASP-Sheet/TRP-Sheet": 3.0, "ASP-Sheet/ALA-Sheet": 4.0, "ASP-Sheet/VAL-Sheet": 2.0, "ASP-Sheet/GLU-Sheet": 4.0, "ASP-Sheet/TYR-Sheet": 4.0, "ASP-Sheet/MET-Sheet": 3.0, "SER-Sheet/CYS-Helix": 1.0, "SER-Sheet/ASP-Helix": 0.0, "SER-Sheet/SER-Helix": -2.0, "SER-Sheet/GLN-Helix": -0.0, "SER-Sheet/LYS-Helix": -0.0, "SER-Sheet/ILE-Helix": -2.0, "SER-Sheet/PRO-Helix": 1.0, "SER-Sheet/THR-Helix": -1.0, "SER-Sheet/PHE-Helix": -1.0, "SER-Sheet/ASN-Helix": -1.0, "SER-Sheet/GLY-Helix": -2.0, "SER-Sheet/HIS-Helix": -2.0, "SER-Sheet/LEU-Helix": -3.0, "SER-Sheet/ARG-Helix": -1.0, "SER-Sheet/TRP-Helix": -2.0, "SER-Sheet/ALA-Helix": 0.0, "SER-Sheet/VAL-Helix": 0.0, "SER-Sheet/GLU-Helix": -1.0, "SER-Sheet/TYR-Helix": 1.0, "SER-Sheet/MET-Helix": -2.0, "SER-Sheet/CYS-Coil": -0.0, "SER-Sheet/ASP-Coil": 1.0, "SER-Sheet/SER-Coil": 2.0, "SER-Sheet/GLN-Coil": 2.0, "SER-Sheet/LYS-Coil": 0.0, "SER-Sheet/ILE-Coil": 1.0, "SER-Sheet/PRO-Coil": 2.0, "SER-Sheet/THR-Coil": 2.0, "SER-Sheet/PHE-Coil": 3.0, "SER-Sheet/ASN-Coil": 1.0, "SER-Sheet/GLY-Coil": 3.0, "SER-Sheet/HIS-Coil": 2.0, "SER-Sheet/LEU-Coil": 2.0, "SER-Sheet/ARG-Coil": -0.0, "SER-Sheet/TRP-Coil": 2.0, "SER-Sheet/ALA-Coil": 1.0, "SER-Sheet/VAL-Coil": 3.0, "SER-Sheet/GLU-Coil": 2.0, "SER-Sheet/TYR-Coil": -2.0, "SER-Sheet/MET-Coil": 3.0, "SER-Sheet/CYS-Sheet": 0.0, "SER-Sheet/ASP-Sheet": 6.0, "SER-Sheet/SER-Sheet": 4.0, "SER-Sheet/GLN-Sheet": 6.0, "SER-Sheet/LYS-Sheet": 5.0, "SER-Sheet/ILE-Sheet": 3.0, "SER-Sheet/PRO-Sheet": 3.0, "SER-Sheet/THR-Sheet": 6.0, "SER-Sheet/PHE-Sheet": 4.0, "SER-Sheet/ASN-Sheet": 5.0, "SER-Sheet/GLY-Sheet": 6.0, "SER-Sheet/HIS-Sheet": 6.0, "SER-Sheet/LEU-Sheet": 3.0, "SER-Sheet/ARG-Sheet": 5.0, "SER-Sheet/TRP-Sheet": 4.0, "SER-Sheet/ALA-Sheet": 4.0, "SER-Sheet/VAL-Sheet": 4.0, "SER-Sheet/GLU-Sheet": 5.0, "SER-Sheet/TYR-Sheet": 3.0, "SER-Sheet/MET-Sheet": 4.0, "GLN-Sheet/CYS-Helix": -Infinity, "GLN-Sheet/ASP-Helix": -1.0, "GLN-Sheet/SER-Helix": -3.0, "GLN-Sheet/GLN-Helix": 2.0, "GLN-Sheet/LYS-Helix": -2.0, "GLN-Sheet/ILE-Helix": -2.0, "GLN-Sheet/PRO-Helix": -1.0, "GLN-Sheet/THR-Helix": -1.0, "GLN-Sheet/PHE-Helix": -Infinity, "GLN-Sheet/ASN-Helix": -1.0, "GLN-Sheet/GLY-Helix": -0.0, "GLN-Sheet/HIS-Helix": -1.0, "GLN-Sheet/LEU-Helix": -1.0, "GLN-Sheet/ARG-Helix": -1.0, "GLN-Sheet/TRP-Helix": -Infinity, "GLN-Sheet/ALA-Helix": -4.0, "GLN-Sheet/VAL-Helix": -Infinity, "GLN-Sheet/GLU-Helix": -0.0, "GLN-Sheet/TYR-Helix": -Infinity, "GLN-Sheet/MET-Helix": -Infinity, "GLN-Sheet/CYS-Coil": -Infinity, "GLN-Sheet/ASP-Coil": 3.0, "GLN-Sheet/SER-Coil": 2.0, "GLN-Sheet/GLN-Coil": -0.0, "GLN-Sheet/LYS-Coil": -Infinity, "GLN-Sheet/ILE-Coil": 1.0, "GLN-Sheet/PRO-Coil": 2.0, "GLN-Sheet/THR-Coil": -0.0, "GLN-Sheet/PHE-Coil": 1.0, "GLN-Sheet/ASN-Coil": 3.0, "GLN-Sheet/GLY-Coil": 3.0, "GLN-Sheet/HIS-Coil": 3.0, "GLN-Sheet/LEU-Coil": 3.0, "GLN-Sheet/ARG-Coil": 2.0, "GLN-Sheet/TRP-Coil": 2.0, "GLN-Sheet/ALA-Coil": 2.0, "GLN-Sheet/VAL-Coil": 4.0, "GLN-Sheet/GLU-Coil": 1.0, "GLN-Sheet/TYR-Coil": -0.0, "GLN-Sheet/MET-Coil": 3.0, "GLN-Sheet/CYS-Sheet": 5.0, "GLN-Sheet/ASP-Sheet": 5.0, "GLN-Sheet/SER-Sheet": 6.0, "GLN-Sheet/GLN-Sheet": 5.0, "GLN-Sheet/LYS-Sheet": 5.0, "GLN-Sheet/ILE-Sheet": 4.0, "GLN-Sheet/PRO-Sheet": 5.0, "GLN-Sheet/THR-Sheet": 5.0, "GLN-Sheet/PHE-Sheet": -0.0, "GLN-Sheet/ASN-Sheet": 6.0, "GLN-Sheet/GLY-Sheet": 4.0, "GLN-Sheet/HIS-Sheet": 4.0, "GLN-Sheet/LEU-Sheet": 4.0, "GLN-Sheet/ARG-Sheet": 5.0, "GLN-Sheet/TRP-Sheet": 5.0, "GLN-Sheet/ALA-Sheet": 4.0, "GLN-Sheet/VAL-Sheet": 4.0, "GLN-Sheet/GLU-Sheet": 4.0, "GLN-Sheet/TYR-Sheet": 3.0, "GLN-Sheet/MET-Sheet": 4.0, "LYS-Sheet/CYS-Helix": 3.0, "LYS-Sheet/ASP-Helix": 2.0, "LYS-Sheet/SER-Helix": -1.0, "LYS-Sheet/GLN-Helix": -1.0, "LYS-Sheet/LYS-Helix": -Infinity, "LYS-Sheet/ILE-Helix": -Infinity, "LYS-Sheet/PRO-Helix": -2.0, "LYS-Sheet/THR-Helix": -3.0, "LYS-Sheet/PHE-Helix": -1.0, "LYS-Sheet/ASN-Helix": -Infinity, "LYS-Sheet/GLY-Helix": -1.0, "LYS-Sheet/HIS-Helix": -1.0, "LYS-Sheet/LEU-Helix": 1.0, "LYS-Sheet/ARG-Helix": -Infinity, "LYS-Sheet/TRP-Helix": 0.0, "LYS-Sheet/ALA-Helix": 0.0, "LYS-Sheet/VAL-Helix": -Infinity, "LYS-Sheet/GLU-Helix": 0.0, "LYS-Sheet/TYR-Helix": -2.0, "LYS-Sheet/MET-Helix": -Infinity, "LYS-Sheet/CYS-Coil": 0.0, "LYS-Sheet/ASP-Coil": 3.0, "LYS-Sheet/SER-Coil": 3.0, "LYS-Sheet/GLN-Coil": 1.0, "LYS-Sheet/LYS-Coil": -1.0, "LYS-Sheet/ILE-Coil": -Infinity, "LYS-Sheet/PRO-Coil": 1.0, "LYS-Sheet/THR-Coil": 3.0, "LYS-Sheet/PHE-Coil": 3.0, "LYS-Sheet/ASN-Coil": 2.0, "LYS-Sheet/GLY-Coil": 2.0, "LYS-Sheet/HIS-Coil": 2.0, "LYS-Sheet/LEU-Coil": -0.0, "LYS-Sheet/ARG-Coil": 2.0, "LYS-Sheet/TRP-Coil": 2.0, "LYS-Sheet/ALA-Coil": 2.0, "LYS-Sheet/VAL-Coil": 2.0, "LYS-Sheet/GLU-Coil": 4.0, "LYS-Sheet/TYR-Coil": 2.0, "LYS-Sheet/MET-Coil": 3.0, "LYS-Sheet/CYS-Sheet": 4.0, "LYS-Sheet/ASP-Sheet": 5.0, "LYS-Sheet/SER-Sheet": 5.0, "LYS-Sheet/GLN-Sheet": 5.0, "LYS-Sheet/LYS-Sheet": 1.0, "LYS-Sheet/ILE-Sheet": 2.0, "LYS-Sheet/PRO-Sheet": 1.0, "LYS-Sheet/THR-Sheet": 5.0, "LYS-Sheet/PHE-Sheet": 3.0, "LYS-Sheet/ASN-Sheet": 5.0, "LYS-Sheet/GLY-Sheet": 4.0, "LYS-Sheet/HIS-Sheet": 4.0, "LYS-Sheet/LEU-Sheet": 0.0, "LYS-Sheet/ARG-Sheet": 3.0, "LYS-Sheet/TRP-Sheet": 6.0, "LYS-Sheet/ALA-Sheet": 2.0, "LYS-Sheet/VAL-Sheet": 3.0, "LYS-Sheet/GLU-Sheet": 7.0, "LYS-Sheet/TYR-Sheet": 4.0, "LYS-Sheet/MET-Sheet": 3.0, "ILE-Sheet/CYS-Helix": 2.0, "ILE-Sheet/ASP-Helix": -2.0, "ILE-Sheet/SER-Helix": -0.0, "ILE-Sheet/GLN-Helix": -0.0, "ILE-Sheet/LYS-Helix": -0.0, "ILE-Sheet/ILE-Helix": 1.0, "ILE-Sheet/PRO-Helix": 1.0, "ILE-Sheet/THR-Helix": -0.0, "ILE-Sheet/PHE-Helix": -1.0, "ILE-Sheet/ASN-Helix": 0.0, "ILE-Sheet/GLY-Helix": -0.0, "ILE-Sheet/HIS-Helix": -2.0, "ILE-Sheet/LEU-Helix": 1.0, "ILE-Sheet/ARG-Helix": -2.0, "ILE-Sheet/TRP-Helix": 1.0, "ILE-Sheet/ALA-Helix": 1.0, "ILE-Sheet/VAL-Helix": 2.0, "ILE-Sheet/GLU-Helix": -1.0, "ILE-Sheet/TYR-Helix": -1.0, "ILE-Sheet/MET-Helix": 1.0, "ILE-Sheet/CYS-Coil": 2.0, "ILE-Sheet/ASP-Coil": 1.0, "ILE-Sheet/SER-Coil": 2.0, "ILE-Sheet/GLN-Coil": 1.0, "ILE-Sheet/LYS-Coil": -2.0, "ILE-Sheet/ILE-Coil": 2.0, "ILE-Sheet/PRO-Coil": 1.0, "ILE-Sheet/THR-Coil": 0.0, "ILE-Sheet/PHE-Coil": -1.0, "ILE-Sheet/ASN-Coil": 1.0, "ILE-Sheet/GLY-Coil": 1.0, "ILE-Sheet/HIS-Coil": 1.0, "ILE-Sheet/LEU-Coil": 2.0, "ILE-Sheet/ARG-Coil": 0.0, "ILE-Sheet/TRP-Coil": 1.0, "ILE-Sheet/ALA-Coil": 2.0, "ILE-Sheet/VAL-Coil": 2.0, "ILE-Sheet/GLU-Coil": -0.0, "ILE-Sheet/TYR-Coil": 0.0, "ILE-Sheet/MET-Coil": 4.0, "ILE-Sheet/CYS-Sheet": 6.0, "ILE-Sheet/ASP-Sheet": 2.0, "ILE-Sheet/SER-Sheet": 3.0, "ILE-Sheet/GLN-Sheet": 4.0, "ILE-Sheet/LYS-Sheet": 2.0, "ILE-Sheet/ILE-Sheet": 4.0, "ILE-Sheet/PRO-Sheet": -Infinity, "ILE-Sheet/THR-Sheet": 3.0, "ILE-Sheet/PHE-Sheet": 5.0, "ILE-Sheet/ASN-Sheet": 2.0, "ILE-Sheet/GLY-Sheet": 4.0, "ILE-Sheet/HIS-Sheet": 2.0, "ILE-Sheet/LEU-Sheet": 5.0, "ILE-Sheet/ARG-Sheet": 3.0, "ILE-Sheet/TRP-Sheet": 3.0, "ILE-Sheet/ALA-Sheet": 6.0, "ILE-Sheet/VAL-Sheet": 5.0, "ILE-Sheet/GLU-Sheet": 3.0, "ILE-Sheet/TYR-Sheet": 3.0, "ILE-Sheet/MET-Sheet": 4.0, "PRO-Sheet/CYS-Helix": -Infinity, "PRO-Sheet/ASP-Helix": 1.0, "PRO-Sheet/SER-Helix": 3.0, "PRO-Sheet/GLN-Helix": 3.0, "PRO-Sheet/LYS-Helix": 2.0, "PRO-Sheet/ILE-Helix": 1.0, "PRO-Sheet/PRO-Helix": -Infinity, "PRO-Sheet/THR-Helix": 1.0, "PRO-Sheet/PHE-Helix": 1.0, "PRO-Sheet/ASN-Helix": 3.0, "PRO-Sheet/GLY-Helix": 0.0, "PRO-Sheet/HIS-Helix": -Infinity, "PRO-Sheet/LEU-Helix": 1.0, "PRO-Sheet/ARG-Helix": 1.0, "PRO-Sheet/TRP-Helix": 2.0, "PRO-Sheet/ALA-Helix": 0.0, "PRO-Sheet/VAL-Helix": 2.0, "PRO-Sheet/GLU-Helix": 3.0, "PRO-Sheet/TYR-Helix": 2.0, "PRO-Sheet/MET-Helix": 3.0, "PRO-Sheet/CYS-Coil": 4.0, "PRO-Sheet/ASP-Coil": -Infinity, "PRO-Sheet/SER-Coil": 2.0, "PRO-Sheet/GLN-Coil": 1.0, "PRO-Sheet/LYS-Coil": 3.0, "PRO-Sheet/ILE-Coil": 3.0, "PRO-Sheet/PRO-Coil": 2.0, "PRO-Sheet/THR-Coil": 2.0, "PRO-Sheet/PHE-Coil": 4.0, "PRO-Sheet/ASN-Coil": -Infinity, "PRO-Sheet/GLY-Coil": 4.0, "PRO-Sheet/HIS-Coil": 2.0, "PRO-Sheet/LEU-Coil": 2.0, "PRO-Sheet/ARG-Coil": 2.0, "PRO-Sheet/TRP-Coil": 3.0, "PRO-Sheet/ALA-Coil": 4.0, "PRO-Sheet/VAL-Coil": 2.0, "PRO-Sheet/GLU-Coil": -Infinity, "PRO-Sheet/TYR-Coil": 1.0, "PRO-Sheet/MET-Coil": 5.0, "PRO-Sheet/CYS-Sheet": 3.0, "PRO-Sheet/ASP-Sheet": -Infinity, "PRO-Sheet/SER-Sheet": 3.0, "PRO-Sheet/GLN-Sheet": 5.0, "PRO-Sheet/LYS-Sheet": 1.0, "PRO-Sheet/ILE-Sheet": -Infinity, "PRO-Sheet/PRO-Sheet": 5.0, "PRO-Sheet/THR-Sheet": 1.0, "PRO-Sheet/PHE-Sheet": 3.0, "PRO-Sheet/ASN-Sheet": -Infinity, "PRO-Sheet/GLY-Sheet": 5.0, "PRO-Sheet/HIS-Sheet": 4.0, "PRO-Sheet/LEU-Sheet": 3.0, "PRO-Sheet/ARG-Sheet": 3.0, "PRO-Sheet/TRP-Sheet": 2.0, "PRO-Sheet/ALA-Sheet": 5.0, "PRO-Sheet/VAL-Sheet": 4.0, "PRO-Sheet/GLU-Sheet": 5.0, "PRO-Sheet/TYR-Sheet": 5.0, "PRO-Sheet/MET-Sheet": 3.0, "THR-Sheet/CYS-Helix": -Infinity, "THR-Sheet/ASP-Helix": -1.0, "THR-Sheet/SER-Helix": 1.0, "THR-Sheet/GLN-Helix": -1.0, "THR-Sheet/LYS-Helix": 0.0, "THR-Sheet/ILE-Helix": 0.0, "THR-Sheet/PRO-Helix": -1.0, "THR-Sheet/THR-Helix": -1.0, "THR-Sheet/PHE-Helix": -3.0, "THR-Sheet/ASN-Helix": -1.0, "THR-Sheet/GLY-Helix": -0.0, "THR-Sheet/HIS-Helix": -Infinity, "THR-Sheet/LEU-Helix": -0.0, "THR-Sheet/ARG-Helix": -1.0, "THR-Sheet/TRP-Helix": -Infinity, "THR-Sheet/ALA-Helix": -1.0, "THR-Sheet/VAL-Helix": -0.0, "THR-Sheet/GLU-Helix": 0.0, "THR-Sheet/TYR-Helix": -Infinity, "THR-Sheet/MET-Helix": -3.0, "THR-Sheet/CYS-Coil": 2.0, "THR-Sheet/ASP-Coil": 3.0, "THR-Sheet/SER-Coil": 2.0, "THR-Sheet/GLN-Coil": 1.0, "THR-Sheet/LYS-Coil": -1.0, "THR-Sheet/ILE-Coil": 0.0, "THR-Sheet/PRO-Coil": 1.0, "THR-Sheet/THR-Coil": 1.0, "THR-Sheet/PHE-Coil": 0.0, "THR-Sheet/ASN-Coil": 1.0, "THR-Sheet/GLY-Coil": 2.0, "THR-Sheet/HIS-Coil": 0.0, "THR-Sheet/LEU-Coil": 0.0, "THR-Sheet/ARG-Coil": -2.0, "THR-Sheet/TRP-Coil": 2.0, "THR-Sheet/ALA-Coil": 2.0, "THR-Sheet/VAL-Coil": 1.0, "THR-Sheet/GLU-Coil": 2.0, "THR-Sheet/TYR-Coil": -0.0, "THR-Sheet/MET-Coil": 2.0, "THR-Sheet/CYS-Sheet": 4.0, "THR-Sheet/ASP-Sheet": 5.0, "THR-Sheet/SER-Sheet": 6.0, "THR-Sheet/GLN-Sheet": 5.0, "THR-Sheet/LYS-Sheet": 5.0, "THR-Sheet/ILE-Sheet": 3.0, "THR-Sheet/PRO-Sheet": 1.0, "THR-Sheet/THR-Sheet": 4.0, "THR-Sheet/PHE-Sheet": 2.0, "THR-Sheet/ASN-Sheet": 4.0, "THR-Sheet/GLY-Sheet": 5.0, "THR-Sheet/HIS-Sheet": 4.0, "THR-Sheet/LEU-Sheet": 4.0, "THR-Sheet/ARG-Sheet": 4.0, "THR-Sheet/TRP-Sheet": 1.0, "THR-Sheet/ALA-Sheet": 4.0, "THR-Sheet/VAL-Sheet": 3.0, "THR-Sheet/GLU-Sheet": 5.0, "THR-Sheet/TYR-Sheet": 3.0, "THR-Sheet/MET-Sheet": 4.0, "PHE-Sheet/CYS-Helix": 1.0, "PHE-Sheet/ASP-Helix": -2.0, "PHE-Sheet/SER-Helix": -2.0, "PHE-Sheet/GLN-Helix": 0.0, "PHE-Sheet/LYS-Helix": -2.0, "PHE-Sheet/ILE-Helix": 1.0, "PHE-Sheet/PRO-Helix": 1.0, "PHE-Sheet/THR-Helix": -Infinity, "PHE-Sheet/PHE-Helix": 2.0, "PHE-Sheet/ASN-Helix": -Infinity, "PHE-Sheet/GLY-Helix": 0.0, "PHE-Sheet/HIS-Helix": 0.0, "PHE-Sheet/LEU-Helix": 1.0, "PHE-Sheet/ARG-Helix": -1.0, "PHE-Sheet/TRP-Helix": -1.0, "PHE-Sheet/ALA-Helix": 1.0, "PHE-Sheet/VAL-Helix": 1.0, "PHE-Sheet/GLU-Helix": -Infinity, "PHE-Sheet/TYR-Helix": -2.0, "PHE-Sheet/MET-Helix": 1.0, "PHE-Sheet/CYS-Coil": 3.0, "PHE-Sheet/ASP-Coil": -1.0, "PHE-Sheet/SER-Coil": 2.0, "PHE-Sheet/GLN-Coil": -1.0, "PHE-Sheet/LYS-Coil": -0.0, "PHE-Sheet/ILE-Coil": 2.0, "PHE-Sheet/PRO-Coil": 2.0, "PHE-Sheet/THR-Coil": -0.0, "PHE-Sheet/PHE-Coil": 2.0, "PHE-Sheet/ASN-Coil": 1.0, "PHE-Sheet/GLY-Coil": 1.0, "PHE-Sheet/HIS-Coil": -1.0, "PHE-Sheet/LEU-Coil": 2.0, "PHE-Sheet/ARG-Coil": -2.0, "PHE-Sheet/TRP-Coil": -Infinity, "PHE-Sheet/ALA-Coil": 1.0, "PHE-Sheet/VAL-Coil": 3.0, "PHE-Sheet/GLU-Coil": 0.0, "PHE-Sheet/TYR-Coil": -Infinity, "PHE-Sheet/MET-Coil": -Infinity, "PHE-Sheet/CYS-Sheet": 4.0, "PHE-Sheet/ASP-Sheet": -0.0, "PHE-Sheet/SER-Sheet": 4.0, "PHE-Sheet/GLN-Sheet": -0.0, "PHE-Sheet/LYS-Sheet": 3.0, "PHE-Sheet/ILE-Sheet": 5.0, "PHE-Sheet/PRO-Sheet": 3.0, "PHE-Sheet/THR-Sheet": 2.0, "PHE-Sheet/PHE-Sheet": 3.0, "PHE-Sheet/ASN-Sheet": -0.0, "PHE-Sheet/GLY-Sheet": 6.0, "PHE-Sheet/HIS-Sheet": 4.0, "PHE-Sheet/LEU-Sheet": 5.0, "PHE-Sheet/ARG-Sheet": 5.0, "PHE-Sheet/TRP-Sheet": 6.0, "PHE-Sheet/ALA-Sheet": 5.0, "PHE-Sheet/VAL-Sheet": 4.0, "PHE-Sheet/GLU-Sheet": 3.0, "PHE-Sheet/TYR-Sheet": 3.0, "PHE-Sheet/MET-Sheet": 4.0, "ASN-Sheet/CYS-Helix": 4.0, "ASN-Sheet/ASP-Helix": -0.0, "ASN-Sheet/SER-Helix": -1.0, "ASN-Sheet/GLN-Helix": -Infinity, "ASN-Sheet/LYS-Helix": -Infinity, "ASN-Sheet/ILE-Helix": 0.0, "ASN-Sheet/PRO-Helix": 2.0, "ASN-Sheet/THR-Helix": 0.0, "ASN-Sheet/PHE-Helix": -Infinity, "ASN-Sheet/ASN-Helix": 1.0, "ASN-Sheet/GLY-Helix": -1.0, "ASN-Sheet/HIS-Helix": -Infinity, "ASN-Sheet/LEU-Helix": 1.0, "ASN-Sheet/ARG-Helix": -Infinity, "ASN-Sheet/TRP-Helix": -0.0, "ASN-Sheet/ALA-Helix": 0.0, "ASN-Sheet/VAL-Helix": -1.0, "ASN-Sheet/GLU-Helix": 0.0, "ASN-Sheet/TYR-Helix": -0.0, "ASN-Sheet/MET-Helix": 1.0, "ASN-Sheet/CYS-Coil": -Infinity, "ASN-Sheet/ASP-Coil": 3.0, "ASN-Sheet/SER-Coil": 4.0, "ASN-Sheet/GLN-Coil": 3.0, "ASN-Sheet/LYS-Coil": -0.0, "ASN-Sheet/ILE-Coil": 3.0, "ASN-Sheet/PRO-Coil": 2.0, "ASN-Sheet/THR-Coil": 3.0, "ASN-Sheet/PHE-Coil": -0.0, "ASN-Sheet/ASN-Coil": 3.0, "ASN-Sheet/GLY-Coil": 3.0, "ASN-Sheet/HIS-Coil": -Infinity, "ASN-Sheet/LEU-Coil": 2.0, "ASN-Sheet/ARG-Coil": 4.0, "ASN-Sheet/TRP-Coil": 4.0, "ASN-Sheet/ALA-Coil": 5.0, "ASN-Sheet/VAL-Coil": 3.0, "ASN-Sheet/GLU-Coil": 3.0, "ASN-Sheet/TYR-Coil": -Infinity, "ASN-Sheet/MET-Coil": 3.0, "ASN-Sheet/CYS-Sheet": 3.0, "ASN-Sheet/ASP-Sheet": 6.0, "ASN-Sheet/SER-Sheet": 5.0, "ASN-Sheet/GLN-Sheet": 6.0, "ASN-Sheet/LYS-Sheet": 5.0, "ASN-Sheet/ILE-Sheet": 2.0, "ASN-Sheet/PRO-Sheet": -Infinity, "ASN-Sheet/THR-Sheet": 4.0, "ASN-Sheet/PHE-Sheet": -0.0, "ASN-Sheet/ASN-Sheet": 6.0, "ASN-Sheet/GLY-Sheet": 6.0, "ASN-Sheet/HIS-Sheet": 5.0, "ASN-Sheet/LEU-Sheet": 1.0, "ASN-Sheet/ARG-Sheet": 4.0, "ASN-Sheet/TRP-Sheet": 5.0, "ASN-Sheet/ALA-Sheet": 5.0, "ASN-Sheet/VAL-Sheet": 1.0, "ASN-Sheet/GLU-Sheet": 5.0, "ASN-Sheet/TYR-Sheet": 5.0, "ASN-Sheet/MET-Sheet": 3.0, "GLY-Sheet/CYS-Helix": 2.0, "GLY-Sheet/ASP-Helix": 0.0, "GLY-Sheet/SER-Helix": -1.0, "GLY-Sheet/GLN-Helix": -0.0, "GLY-Sheet/LYS-Helix": -Infinity, "GLY-Sheet/ILE-Helix": 0.0, "GLY-Sheet/PRO-Helix": 1.0, "GLY-Sheet/THR-Helix": 0.0, "GLY-Sheet/PHE-Helix": -2.0, "GLY-Sheet/ASN-Helix": -Infinity, "GLY-Sheet/GLY-Helix": -0.0, "GLY-Sheet/HIS-Helix": -Infinity, "GLY-Sheet/LEU-Helix": 1.0, "GLY-Sheet/ARG-Helix": -1.0, "GLY-Sheet/TRP-Helix": -0.0, "GLY-Sheet/ALA-Helix": -0.0, "GLY-Sheet/VAL-Helix": -0.0, "GLY-Sheet/GLU-Helix": -0.0, "GLY-Sheet/TYR-Helix": 0.0, "GLY-Sheet/MET-Helix": 2.0, "GLY-Sheet/CYS-Coil": 2.0, "GLY-Sheet/ASP-Coil": 2.0, "GLY-Sheet/SER-Coil": 3.0, "GLY-Sheet/GLN-Coil": -1.0, "GLY-Sheet/LYS-Coil": -Infinity, "GLY-Sheet/ILE-Coil": 3.0, "GLY-Sheet/PRO-Coil": 1.0, "GLY-Sheet/THR-Coil": 3.0, "GLY-Sheet/PHE-Coil": 2.0, "GLY-Sheet/ASN-Coil": 1.0, "GLY-Sheet/GLY-Coil": 3.0, "GLY-Sheet/HIS-Coil": 1.0, "GLY-Sheet/LEU-Coil": 3.0, "GLY-Sheet/ARG-Coil": -2.0, "GLY-Sheet/TRP-Coil": -Infinity, "GLY-Sheet/ALA-Coil": 2.0, "GLY-Sheet/VAL-Coil": 4.0, "GLY-Sheet/GLU-Coil": -2.0, "GLY-Sheet/TYR-Coil": 3.0, "GLY-Sheet/MET-Coil": 3.0, "GLY-Sheet/CYS-Sheet": 6.0, "GLY-Sheet/ASP-Sheet": 4.0, "GLY-Sheet/SER-Sheet": 6.0, "GLY-Sheet/GLN-Sheet": 4.0, "GLY-Sheet/LYS-Sheet": 4.0, "GLY-Sheet/ILE-Sheet": 4.0, "GLY-Sheet/PRO-Sheet": 5.0, "GLY-Sheet/THR-Sheet": 5.0, "GLY-Sheet/PHE-Sheet": 6.0, "GLY-Sheet/ASN-Sheet": 6.0, "GLY-Sheet/GLY-Sheet": 5.0, "GLY-Sheet/HIS-Sheet": 6.0, "GLY-Sheet/LEU-Sheet": 5.0, "GLY-Sheet/ARG-Sheet": 3.0, "GLY-Sheet/TRP-Sheet": 5.0, "GLY-Sheet/ALA-Sheet": 6.0, "GLY-Sheet/VAL-Sheet": 5.0, "GLY-Sheet/GLU-Sheet": 3.0, "GLY-Sheet/TYR-Sheet": 6.0, "GLY-Sheet/MET-Sheet": 4.0, "HIS-Sheet/CYS-Helix": -Infinity, "HIS-Sheet/ASP-Helix": 2.0, "HIS-Sheet/SER-Helix": -0.0, "HIS-Sheet/GLN-Helix": -Infinity, "HIS-Sheet/LYS-Helix": 1.0, "HIS-Sheet/ILE-Helix": -Infinity, "HIS-Sheet/PRO-Helix": -Infinity, "HIS-Sheet/THR-Helix": -1.0, "HIS-Sheet/PHE-Helix": -Infinity, "HIS-Sheet/ASN-Helix": 1.0, "HIS-Sheet/GLY-Helix": -Infinity, "HIS-Sheet/HIS-Helix": -Infinity, "HIS-Sheet/LEU-Helix": 0.0, "HIS-Sheet/ARG-Helix": 0.0, "HIS-Sheet/TRP-Helix": 0.0, "HIS-Sheet/ALA-Helix": 0.0, "HIS-Sheet/VAL-Helix": 1.0, "HIS-Sheet/GLU-Helix": -1.0, "HIS-Sheet/TYR-Helix": 1.0, "HIS-Sheet/MET-Helix": 2.0, "HIS-Sheet/CYS-Coil": 5.0, "HIS-Sheet/ASP-Coil": 4.0, "HIS-Sheet/SER-Coil": 3.0, "HIS-Sheet/GLN-Coil": 3.0, "HIS-Sheet/LYS-Coil": 2.0, "HIS-Sheet/ILE-Coil": 2.0, "HIS-Sheet/PRO-Coil": 2.0, "HIS-Sheet/THR-Coil": 4.0, "HIS-Sheet/PHE-Coil": 0.0, "HIS-Sheet/ASN-Coil": 1.0, "HIS-Sheet/GLY-Coil": 1.0, "HIS-Sheet/HIS-Coil": 1.0, "HIS-Sheet/LEU-Coil": 3.0, "HIS-Sheet/ARG-Coil": 1.0, "HIS-Sheet/TRP-Coil": 3.0, "HIS-Sheet/ALA-Coil": 2.0, "HIS-Sheet/VAL-Coil": 2.0, "HIS-Sheet/GLU-Coil": 4.0, "HIS-Sheet/TYR-Coil": -Infinity, "HIS-Sheet/MET-Coil": 4.0, "HIS-Sheet/CYS-Sheet": 5.0, "HIS-Sheet/ASP-Sheet": 5.0, "HIS-Sheet/SER-Sheet": 6.0, "HIS-Sheet/GLN-Sheet": 4.0, "HIS-Sheet/LYS-Sheet": 4.0, "HIS-Sheet/ILE-Sheet": 2.0, "HIS-Sheet/PRO-Sheet": 4.0, "HIS-Sheet/THR-Sheet": 4.0, "HIS-Sheet/PHE-Sheet": 4.0, "HIS-Sheet/ASN-Sheet": 5.0, "HIS-Sheet/GLY-Sheet": 6.0, "HIS-Sheet/HIS-Sheet": 4.0, "HIS-Sheet/LEU-Sheet": 2.0, "HIS-Sheet/ARG-Sheet": 4.0, "HIS-Sheet/TRP-Sheet": 4.0, "HIS-Sheet/ALA-Sheet": 5.0, "HIS-Sheet/VAL-Sheet": 3.0, "HIS-Sheet/GLU-Sheet": 4.0, "HIS-Sheet/TYR-Sheet": 6.0, "HIS-Sheet/MET-Sheet": 5.0, "LEU-Sheet/CYS-Helix": 1.0, "LEU-Sheet/ASP-Helix": -1.0, "LEU-Sheet/SER-Helix": -1.0, "LEU-Sheet/GLN-Helix": 0.0, "LEU-Sheet/LYS-Helix": -3.0, "LEU-Sheet/ILE-Helix": -0.0, "LEU-Sheet/PRO-Helix": -3.0, "LEU-Sheet/THR-Helix": -1.0, "LEU-Sheet/PHE-Helix": 1.0, "LEU-Sheet/ASN-Helix": 0.0, "LEU-Sheet/GLY-Helix": -2.0, "LEU-Sheet/HIS-Helix": 1.0, "LEU-Sheet/LEU-Helix": -0.0, "LEU-Sheet/ARG-Helix": -3.0, "LEU-Sheet/TRP-Helix": -1.0, "LEU-Sheet/ALA-Helix": 1.0, "LEU-Sheet/VAL-Helix": -2.0, "LEU-Sheet/GLU-Helix": -1.0, "LEU-Sheet/TYR-Helix": 1.0, "LEU-Sheet/MET-Helix": 0.0, "LEU-Sheet/CYS-Coil": 3.0, "LEU-Sheet/ASP-Coil": 1.0, "LEU-Sheet/SER-Coil": 1.0, "LEU-Sheet/GLN-Coil": 1.0, "LEU-Sheet/LYS-Coil": -1.0, "LEU-Sheet/ILE-Coil": 3.0, "LEU-Sheet/PRO-Coil": 2.0, "LEU-Sheet/THR-Coil": -2.0, "LEU-Sheet/PHE-Coil": 2.0, "LEU-Sheet/ASN-Coil": -1.0, "LEU-Sheet/GLY-Coil": 2.0, "LEU-Sheet/HIS-Coil": 0.0, "LEU-Sheet/LEU-Coil": 3.0, "LEU-Sheet/ARG-Coil": 1.0, "LEU-Sheet/TRP-Coil": 3.0, "LEU-Sheet/ALA-Coil": 2.0, "LEU-Sheet/VAL-Coil": 1.0, "LEU-Sheet/GLU-Coil": 1.0, "LEU-Sheet/TYR-Coil": 1.0, "LEU-Sheet/MET-Coil": 3.0, "LEU-Sheet/CYS-Sheet": 6.0, "LEU-Sheet/ASP-Sheet": 3.0, "LEU-Sheet/SER-Sheet": 3.0, "LEU-Sheet/GLN-Sheet": 4.0, "LEU-Sheet/LYS-Sheet": 0.0, "LEU-Sheet/ILE-Sheet": 5.0, "LEU-Sheet/PRO-Sheet": 3.0, "LEU-Sheet/THR-Sheet": 4.0, "LEU-Sheet/PHE-Sheet": 5.0, "LEU-Sheet/ASN-Sheet": 1.0, "LEU-Sheet/GLY-Sheet": 5.0, "LEU-Sheet/HIS-Sheet": 2.0, "LEU-Sheet/LEU-Sheet": 4.0, "LEU-Sheet/ARG-Sheet": 3.0, "LEU-Sheet/TRP-Sheet": 5.0, "LEU-Sheet/ALA-Sheet": 5.0, "LEU-Sheet/VAL-Sheet": 5.0, "LEU-Sheet/GLU-Sheet": 4.0, "LEU-Sheet/TYR-Sheet": 4.0, "LEU-Sheet/MET-Sheet": 6.0, "ARG-Sheet/CYS-Helix": -Infinity, "ARG-Sheet/ASP-Helix": 0.0, "ARG-Sheet/SER-Helix": 1.0, "ARG-Sheet/GLN-Helix": 1.0, "ARG-Sheet/LYS-Helix": -2.0, "ARG-Sheet/ILE-Helix": -1.0, "ARG-Sheet/PRO-Helix": -Infinity, "ARG-Sheet/THR-Helix": -1.0, "ARG-Sheet/PHE-Helix": 2.0, "ARG-Sheet/ASN-Helix": 0.0, "ARG-Sheet/GLY-Helix": 0.0, "ARG-Sheet/HIS-Helix": -1.0, "ARG-Sheet/LEU-Helix": -3.0, "ARG-Sheet/ARG-Helix": 0.0, "ARG-Sheet/TRP-Helix": 1.0, "ARG-Sheet/ALA-Helix": -1.0, "ARG-Sheet/VAL-Helix": -2.0, "ARG-Sheet/GLU-Helix": 0.0, "ARG-Sheet/TYR-Helix": -1.0, "ARG-Sheet/MET-Helix": -Infinity, "ARG-Sheet/CYS-Coil": 2.0, "ARG-Sheet/ASP-Coil": 4.0, "ARG-Sheet/SER-Coil": 2.0, "ARG-Sheet/GLN-Coil": 3.0, "ARG-Sheet/LYS-Coil": -1.0, "ARG-Sheet/ILE-Coil": 1.0, "ARG-Sheet/PRO-Coil": 2.0, "ARG-Sheet/THR-Coil": 2.0, "ARG-Sheet/PHE-Coil": 2.0, "ARG-Sheet/ASN-Coil": 2.0, "ARG-Sheet/GLY-Coil": 2.0, "ARG-Sheet/HIS-Coil": 3.0, "ARG-Sheet/LEU-Coil": 2.0, "ARG-Sheet/ARG-Coil": 3.0, "ARG-Sheet/TRP-Coil": 4.0, "ARG-Sheet/ALA-Coil": 1.0, "ARG-Sheet/VAL-Coil": 3.0, "ARG-Sheet/GLU-Coil": 2.0, "ARG-Sheet/TYR-Coil": 2.0, "ARG-Sheet/MET-Coil": 3.0, "ARG-Sheet/CYS-Sheet": 1.0, "ARG-Sheet/ASP-Sheet": 7.0, "ARG-Sheet/SER-Sheet": 5.0, "ARG-Sheet/GLN-Sheet": 5.0, "ARG-Sheet/LYS-Sheet": 3.0, "ARG-Sheet/ILE-Sheet": 3.0, "ARG-Sheet/PRO-Sheet": 3.0, "ARG-Sheet/THR-Sheet": 4.0, "ARG-Sheet/PHE-Sheet": 5.0, "ARG-Sheet/ASN-Sheet": 4.0, "ARG-Sheet/GLY-Sheet": 3.0, "ARG-Sheet/HIS-Sheet": 4.0, "ARG-Sheet/LEU-Sheet": 3.0, "ARG-Sheet/ARG-Sheet": 3.0, "ARG-Sheet/TRP-Sheet": 5.0, "ARG-Sheet/ALA-Sheet": 5.0, "ARG-Sheet/VAL-Sheet": 4.0, "ARG-Sheet/GLU-Sheet": 7.0, "ARG-Sheet/TYR-Sheet": 4.0, "ARG-Sheet/MET-Sheet": 4.0, "TRP-Sheet/CYS-Helix": -Infinity, "TRP-Sheet/ASP-Helix": -Infinity, "TRP-Sheet/SER-Helix": -Infinity, "TRP-Sheet/GLN-Helix": 0.0, "TRP-Sheet/LYS-Helix": 1.0, "TRP-Sheet/ILE-Helix": 1.0, "TRP-Sheet/PRO-Helix": 2.0, "TRP-Sheet/THR-Helix": -0.0, "TRP-Sheet/PHE-Helix": -2.0, "TRP-Sheet/ASN-Helix": -1.0, "TRP-Sheet/GLY-Helix": -1.0, "TRP-Sheet/HIS-Helix": -Infinity, "TRP-Sheet/LEU-Helix": 2.0, "TRP-Sheet/ARG-Helix": -2.0, "TRP-Sheet/TRP-Helix": -Infinity, "TRP-Sheet/ALA-Helix": 0.0, "TRP-Sheet/VAL-Helix": -Infinity, "TRP-Sheet/GLU-Helix": -Infinity, "TRP-Sheet/TYR-Helix": -1.0, "TRP-Sheet/MET-Helix": 2.0, "TRP-Sheet/CYS-Coil": -Infinity, "TRP-Sheet/ASP-Coil": -Infinity, "TRP-Sheet/SER-Coil": 2.0, "TRP-Sheet/GLN-Coil": 2.0, "TRP-Sheet/LYS-Coil": -0.0, "TRP-Sheet/ILE-Coil": -Infinity, "TRP-Sheet/PRO-Coil": 3.0, "TRP-Sheet/THR-Coil": -Infinity, "TRP-Sheet/PHE-Coil": -Infinity, "TRP-Sheet/ASN-Coil": 0.0, "TRP-Sheet/GLY-Coil": 2.0, "TRP-Sheet/HIS-Coil": 4.0, "TRP-Sheet/LEU-Coil": 3.0, "TRP-Sheet/ARG-Coil": 4.0, "TRP-Sheet/TRP-Coil": -Infinity, "TRP-Sheet/ALA-Coil": 3.0, "TRP-Sheet/VAL-Coil": 3.0, "TRP-Sheet/GLU-Coil": -Infinity, "TRP-Sheet/TYR-Coil": -Infinity, "TRP-Sheet/MET-Coil": 3.0, "TRP-Sheet/CYS-Sheet": 4.0, "TRP-Sheet/ASP-Sheet": 3.0, "TRP-Sheet/SER-Sheet": 4.0, "TRP-Sheet/GLN-Sheet": 5.0, "TRP-Sheet/LYS-Sheet": 6.0, "TRP-Sheet/ILE-Sheet": 3.0, "TRP-Sheet/PRO-Sheet": 2.0, "TRP-Sheet/THR-Sheet": 1.0, "TRP-Sheet/PHE-Sheet": 6.0, "TRP-Sheet/ASN-Sheet": 5.0, "TRP-Sheet/GLY-Sheet": 5.0, "TRP-Sheet/HIS-Sheet": 4.0, "TRP-Sheet/LEU-Sheet": 5.0, "TRP-Sheet/ARG-Sheet": 5.0, "TRP-Sheet/TRP-Sheet": -Infinity, "TRP-Sheet/ALA-Sheet": 6.0, "TRP-Sheet/VAL-Sheet": 3.0, "TRP-Sheet/GLU-Sheet": -Infinity, "TRP-Sheet/TYR-Sheet": -0.0, "TRP-Sheet/MET-Sheet": 5.0, "ALA-Sheet/CYS-Helix": 1.0, "ALA-Sheet/ASP-Helix": -3.0, "ALA-Sheet/SER-Helix": 1.0, "ALA-Sheet/GLN-Helix": 0.0, "ALA-Sheet/LYS-Helix": -1.0, "ALA-Sheet/ILE-Helix": 1.0, "ALA-Sheet/PRO-Helix": 1.0, "ALA-Sheet/THR-Helix": 1.0, "ALA-Sheet/PHE-Helix": -2.0, "ALA-Sheet/ASN-Helix": -1.0, "ALA-Sheet/GLY-Helix": -1.0, "ALA-Sheet/HIS-Helix": 0.0, "ALA-Sheet/LEU-Helix": 2.0, "ALA-Sheet/ARG-Helix": -2.0, "ALA-Sheet/TRP-Helix": 1.0, "ALA-Sheet/ALA-Helix": 1.0, "ALA-Sheet/VAL-Helix": 2.0, "ALA-Sheet/GLU-Helix": 0.0, "ALA-Sheet/TYR-Helix": 1.0, "ALA-Sheet/MET-Helix": 3.0, "ALA-Sheet/CYS-Coil": 3.0, "ALA-Sheet/ASP-Coil": 1.0, "ALA-Sheet/SER-Coil": 2.0, "ALA-Sheet/GLN-Coil": 2.0, "ALA-Sheet/LYS-Coil": -2.0, "ALA-Sheet/ILE-Coil": 3.0, "ALA-Sheet/PRO-Coil": 3.0, "ALA-Sheet/THR-Coil": 1.0, "ALA-Sheet/PHE-Coil": 3.0, "ALA-Sheet/ASN-Coil": 2.0, "ALA-Sheet/GLY-Coil": 2.0, "ALA-Sheet/HIS-Coil": -Infinity, "ALA-Sheet/LEU-Coil": 4.0, "ALA-Sheet/ARG-Coil": -1.0, "ALA-Sheet/TRP-Coil": -0.0, "ALA-Sheet/ALA-Coil": 3.0, "ALA-Sheet/VAL-Coil": 3.0, "ALA-Sheet/GLU-Coil": -2.0, "ALA-Sheet/TYR-Coil": 2.0, "ALA-Sheet/MET-Coil": 2.0, "ALA-Sheet/CYS-Sheet": 5.0, "ALA-Sheet/ASP-Sheet": 4.0, "ALA-Sheet/SER-Sheet": 4.0, "ALA-Sheet/GLN-Sheet": 4.0, "ALA-Sheet/LYS-Sheet": 2.0, "ALA-Sheet/ILE-Sheet": 6.0, "ALA-Sheet/PRO-Sheet": 5.0, "ALA-Sheet/THR-Sheet": 4.0, "ALA-Sheet/PHE-Sheet": 5.0, "ALA-Sheet/ASN-Sheet": 5.0, "ALA-Sheet/GLY-Sheet": 6.0, "ALA-Sheet/HIS-Sheet": 5.0, "ALA-Sheet/LEU-Sheet": 5.0, "ALA-Sheet/ARG-Sheet": 5.0, "ALA-Sheet/TRP-Sheet": 6.0, "ALA-Sheet/ALA-Sheet": 4.0, "ALA-Sheet/VAL-Sheet": 5.0, "ALA-Sheet/GLU-Sheet": 4.0, "ALA-Sheet/TYR-Sheet": 5.0, "ALA-Sheet/MET-Sheet": 5.0, "VAL-Sheet/CYS-Helix": 2.0, "VAL-Sheet/ASP-Helix": -2.0, "VAL-Sheet/SER-Helix": -1.0, "VAL-Sheet/GLN-Helix": -0.0, "VAL-Sheet/LYS-Helix": -1.0, "VAL-Sheet/ILE-Helix": 1.0, "VAL-Sheet/PRO-Helix": -0.0, "VAL-Sheet/THR-Helix": 0.0, "VAL-Sheet/PHE-Helix": -1.0, "VAL-Sheet/ASN-Helix": -2.0, "VAL-Sheet/GLY-Helix": -1.0, "VAL-Sheet/HIS-Helix": 0.0, "VAL-Sheet/LEU-Helix": -0.0, "VAL-Sheet/ARG-Helix": -4.0, "VAL-Sheet/TRP-Helix": -1.0, "VAL-Sheet/ALA-Helix": 2.0, "VAL-Sheet/VAL-Helix": 1.0, "VAL-Sheet/GLU-Helix": -2.0, "VAL-Sheet/TYR-Helix": 0.0, "VAL-Sheet/MET-Helix": 1.0, "VAL-Sheet/CYS-Coil": 2.0, "VAL-Sheet/ASP-Coil": -4.0, "VAL-Sheet/SER-Coil": 1.0, "VAL-Sheet/GLN-Coil": 2.0, "VAL-Sheet/LYS-Coil": 0.0, "VAL-Sheet/ILE-Coil": 1.0, "VAL-Sheet/PRO-Coil": 1.0, "VAL-Sheet/THR-Coil": -1.0, "VAL-Sheet/PHE-Coil": 2.0, "VAL-Sheet/ASN-Coil": 2.0, "VAL-Sheet/GLY-Coil": 1.0, "VAL-Sheet/HIS-Coil": -Infinity, "VAL-Sheet/LEU-Coil": 3.0, "VAL-Sheet/ARG-Coil": 2.0, "VAL-Sheet/TRP-Coil": 2.0, "VAL-Sheet/ALA-Coil": 2.0, "VAL-Sheet/VAL-Coil": 2.0, "VAL-Sheet/GLU-Coil": 1.0, "VAL-Sheet/TYR-Coil": 2.0, "VAL-Sheet/MET-Coil": -1.0, "VAL-Sheet/CYS-Sheet": 6.0, "VAL-Sheet/ASP-Sheet": 2.0, "VAL-Sheet/SER-Sheet": 4.0, "VAL-Sheet/GLN-Sheet": 4.0, "VAL-Sheet/LYS-Sheet": 3.0, "VAL-Sheet/ILE-Sheet": 5.0, "VAL-Sheet/PRO-Sheet": 4.0, "VAL-Sheet/THR-Sheet": 3.0, "VAL-Sheet/PHE-Sheet": 4.0, "VAL-Sheet/ASN-Sheet": 1.0, "VAL-Sheet/GLY-Sheet": 5.0, "VAL-Sheet/HIS-Sheet": 3.0, "VAL-Sheet/LEU-Sheet": 5.0, "VAL-Sheet/ARG-Sheet": 4.0, "VAL-Sheet/TRP-Sheet": 3.0, "VAL-Sheet/ALA-Sheet": 5.0, "VAL-Sheet/VAL-Sheet": 4.0, "VAL-Sheet/GLU-Sheet": 4.0, "VAL-Sheet/TYR-Sheet": 5.0, "VAL-Sheet/MET-Sheet": 5.0, "GLU-Sheet/CYS-Helix": 1.0, "GLU-Sheet/ASP-Helix": -Infinity, "GLU-Sheet/SER-Helix": -Infinity, "GLU-Sheet/GLN-Helix": -2.0, "GLU-Sheet/LYS-Helix": 1.0, "GLU-Sheet/ILE-Helix": -1.0, "GLU-Sheet/PRO-Helix": -1.0, "GLU-Sheet/THR-Helix": -1.0, "GLU-Sheet/PHE-Helix": -Infinity, "GLU-Sheet/ASN-Helix": -Infinity, "GLU-Sheet/GLY-Helix": -2.0, "GLU-Sheet/HIS-Helix": -Infinity, "GLU-Sheet/LEU-Helix": -0.0, "GLU-Sheet/ARG-Helix": 2.0, "GLU-Sheet/TRP-Helix": -Infinity, "GLU-Sheet/ALA-Helix": -1.0, "GLU-Sheet/VAL-Helix": -Infinity, "GLU-Sheet/GLU-Helix": -Infinity, "GLU-Sheet/TYR-Helix": -Infinity, "GLU-Sheet/MET-Helix": -Infinity, "GLU-Sheet/CYS-Coil": 2.0, "GLU-Sheet/ASP-Coil": -1.0, "GLU-Sheet/SER-Coil": 2.0, "GLU-Sheet/GLN-Coil": 3.0, "GLU-Sheet/LYS-Coil": 4.0, "GLU-Sheet/ILE-Coil": 2.0, "GLU-Sheet/PRO-Coil": 1.0, "GLU-Sheet/THR-Coil": 2.0, "GLU-Sheet/PHE-Coil": 3.0, "GLU-Sheet/ASN-Coil": 2.0, "GLU-Sheet/GLY-Coil": 2.0, "GLU-Sheet/HIS-Coil": 1.0, "GLU-Sheet/LEU-Coil": 1.0, "GLU-Sheet/ARG-Coil": 4.0, "GLU-Sheet/TRP-Coil": 1.0, "GLU-Sheet/ALA-Coil": 1.0, "GLU-Sheet/VAL-Coil": 1.0, "GLU-Sheet/GLU-Coil": 1.0, "GLU-Sheet/TYR-Coil": -1.0, "GLU-Sheet/MET-Coil": 1.0, "GLU-Sheet/CYS-Sheet": 4.0, "GLU-Sheet/ASP-Sheet": 4.0, "GLU-Sheet/SER-Sheet": 5.0, "GLU-Sheet/GLN-Sheet": 4.0, "GLU-Sheet/LYS-Sheet": 7.0, "GLU-Sheet/ILE-Sheet": 3.0, "GLU-Sheet/PRO-Sheet": 5.0, "GLU-Sheet/THR-Sheet": 5.0, "GLU-Sheet/PHE-Sheet": 3.0, "GLU-Sheet/ASN-Sheet": 5.0, "GLU-Sheet/GLY-Sheet": 3.0, "GLU-Sheet/HIS-Sheet": 4.0, "GLU-Sheet/LEU-Sheet": 4.0, "GLU-Sheet/ARG-Sheet": 7.0, "GLU-Sheet/TRP-Sheet": -Infinity, "GLU-Sheet/ALA-Sheet": 4.0, "GLU-Sheet/VAL-Sheet": 4.0, "GLU-Sheet/GLU-Sheet": 2.0, "GLU-Sheet/TYR-Sheet": 3.0, "GLU-Sheet/MET-Sheet": 3.0, "TYR-Sheet/CYS-Helix": 3.0, "TYR-Sheet/ASP-Helix": -3.0, "TYR-Sheet/SER-Helix": 0.0, "TYR-Sheet/GLN-Helix": -1.0, "TYR-Sheet/LYS-Helix": -1.0, "TYR-Sheet/ILE-Helix": 2.0, "TYR-Sheet/PRO-Helix": 2.0, "TYR-Sheet/THR-Helix": -3.0, "TYR-Sheet/PHE-Helix": 0.0, "TYR-Sheet/ASN-Helix": -0.0, "TYR-Sheet/GLY-Helix": -1.0, "TYR-Sheet/HIS-Helix": -Infinity, "TYR-Sheet/LEU-Helix": 1.0, "TYR-Sheet/ARG-Helix": -2.0, "TYR-Sheet/TRP-Helix": -Infinity, "TYR-Sheet/ALA-Helix": 0.0, "TYR-Sheet/VAL-Helix": 1.0, "TYR-Sheet/GLU-Helix": -3.0, "TYR-Sheet/TYR-Helix": 1.0, "TYR-Sheet/MET-Helix": 2.0, "TYR-Sheet/CYS-Coil": 1.0, "TYR-Sheet/ASP-Coil": 1.0, "TYR-Sheet/SER-Coil": 1.0, "TYR-Sheet/GLN-Coil": 2.0, "TYR-Sheet/LYS-Coil": 3.0, "TYR-Sheet/ILE-Coil": 1.0, "TYR-Sheet/PRO-Coil": 4.0, "TYR-Sheet/THR-Coil": 1.0, "TYR-Sheet/PHE-Coil": 4.0, "TYR-Sheet/ASN-Coil": 0.0, "TYR-Sheet/GLY-Coil": 2.0, "TYR-Sheet/HIS-Coil": 1.0, "TYR-Sheet/LEU-Coil": 3.0, "TYR-Sheet/ARG-Coil": 1.0, "TYR-Sheet/TRP-Coil": -Infinity, "TYR-Sheet/ALA-Coil": 2.0, "TYR-Sheet/VAL-Coil": 3.0, "TYR-Sheet/GLU-Coil": 1.0, "TYR-Sheet/TYR-Coil": 1.0, "TYR-Sheet/MET-Coil": 4.0, "TYR-Sheet/CYS-Sheet": 5.0, "TYR-Sheet/ASP-Sheet": 4.0, "TYR-Sheet/SER-Sheet": 3.0, "TYR-Sheet/GLN-Sheet": 3.0, "TYR-Sheet/LYS-Sheet": 4.0, "TYR-Sheet/ILE-Sheet": 3.0, "TYR-Sheet/PRO-Sheet": 5.0, "TYR-Sheet/THR-Sheet": 3.0, "TYR-Sheet/PHE-Sheet": 3.0, "TYR-Sheet/ASN-Sheet": 5.0, "TYR-Sheet/GLY-Sheet": 6.0, "TYR-Sheet/HIS-Sheet": 6.0, "TYR-Sheet/LEU-Sheet": 4.0, "TYR-Sheet/ARG-Sheet": 4.0, "TYR-Sheet/TRP-Sheet": -0.0, "TYR-Sheet/ALA-Sheet": 5.0, "TYR-Sheet/VAL-Sheet": 5.0, "TYR-Sheet/GLU-Sheet": 3.0, "TYR-Sheet/TYR-Sheet": 1.0, "TYR-Sheet/MET-Sheet": 4.0, "MET-Sheet/CYS-Helix": 3.0, "MET-Sheet/ASP-Helix": -Infinity, "MET-Sheet/SER-Helix": -Infinity, "MET-Sheet/GLN-Helix": 1.0, "MET-Sheet/LYS-Helix": -1.0, "MET-Sheet/ILE-Helix": -0.0, "MET-Sheet/PRO-Helix": 1.0, "MET-Sheet/THR-Helix": 1.0, "MET-Sheet/PHE-Helix": -Infinity, "MET-Sheet/ASN-Helix": -Infinity, "MET-Sheet/GLY-Helix": -Infinity, "MET-Sheet/HIS-Helix": -Infinity, "MET-Sheet/LEU-Helix": 2.0, "MET-Sheet/ARG-Helix": 1.0, "MET-Sheet/TRP-Helix": 0.0, "MET-Sheet/ALA-Helix": 2.0, "MET-Sheet/VAL-Helix": 1.0, "MET-Sheet/GLU-Helix": -Infinity, "MET-Sheet/TYR-Helix": -1.0, "MET-Sheet/MET-Helix": 1.0, "MET-Sheet/CYS-Coil": 2.0, "MET-Sheet/ASP-Coil": 0.0, "MET-Sheet/SER-Coil": 2.0, "MET-Sheet/GLN-Coil": 1.0, "MET-Sheet/LYS-Coil": 3.0, "MET-Sheet/ILE-Coil": 2.0, "MET-Sheet/PRO-Coil": 2.0, "MET-Sheet/THR-Coil": -1.0, "MET-Sheet/PHE-Coil": 2.0, "MET-Sheet/ASN-Coil": 3.0, "MET-Sheet/GLY-Coil": 2.0, "MET-Sheet/HIS-Coil": -Infinity, "MET-Sheet/LEU-Coil": 2.0, "MET-Sheet/ARG-Coil": 2.0, "MET-Sheet/TRP-Coil": 3.0, "MET-Sheet/ALA-Coil": 3.0, "MET-Sheet/VAL-Coil": -Infinity, "MET-Sheet/GLU-Coil": 3.0, "MET-Sheet/TYR-Coil": 0.0, "MET-Sheet/MET-Coil": 2.0, "MET-Sheet/CYS-Sheet": 6.0, "MET-Sheet/ASP-Sheet": 3.0, "MET-Sheet/SER-Sheet": 4.0, "MET-Sheet/GLN-Sheet": 4.0, "MET-Sheet/LYS-Sheet": 3.0, "MET-Sheet/ILE-Sheet": 4.0, "MET-Sheet/PRO-Sheet": 3.0, "MET-Sheet/THR-Sheet": 4.0, "MET-Sheet/PHE-Sheet": 4.0, "MET-Sheet/ASN-Sheet": 3.0, "MET-Sheet/GLY-Sheet": 4.0, "MET-Sheet/HIS-Sheet": 5.0, "MET-Sheet/LEU-Sheet": 6.0, "MET-Sheet/ARG-Sheet": 4.0, "MET-Sheet/TRP-Sheet": 5.0, "MET-Sheet/ALA-Sheet": 5.0, "MET-Sheet/VAL-Sheet": 5.0, "MET-Sheet/GLU-Sheet": 3.0, "MET-Sheet/TYR-Sheet": 4.0, "MET-Sheet/MET-Sheet": 4.0
};
// Scaled Kyte-Doolittle Hyropathicity scores (DOI: 10.1016/0022-2836(82)90515-0).
var SCALED_KD_HYDROP = {
    "ALA": 0.700,
    "ARG": 0.000,
    "ASN": 0.111,
    "ASP": 0.111, 
    "CYS": 0.778, 
    "GLN": 0.111, 
    "GLU": 0.111, 
    "GLY": 0.456, 
    "HIS": 0.144, 
    "ILE": 1.000, 
    "LEU": 0.922, 
    "LYS": 0.067, 
    "MET": 0.711, 
    "PHE": 0.811, 
    "PRO": 0.322, 
    "SER": 0.411, 
    "THR": 0.422, 
    "TRP": 0.400, 
    "TYR": 0.356, 
    "VAL": 0.967  
};
// Mapping of three letter amino acid codes to one letter amino acid codes. Includes any amino acid (XXX), termination (TER) and ambiguous amino acid (AMB).
var AA_3TO1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "ANY": "X",
    "TER": "Z",
    "INC": "U"
};
// Mapping used to invert nucleotide symbols. Includes lower-case characters and _ (insertion starts), deletion (-) and any nucleotide (N).
var INVERT_BASE = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "-": "-",
    "N": "N"
};
// Style specifications used for 3D protein models, for more information see: https://3dmol.csb.pitt.edu/doc/types.html
var STYLES_3DMOL = {
    "LINE_FADED": {
            colorfunc: (atom) => {
                return COLORS[ atom.resn ];
            },
            radius: 0.01
    },
    "LINE_HIGHLIGHT": {
        colorfunc: (atom) => {
            return COLORS[ atom.resn ];
        },
        radius: 0.2
    },
    "CARTOON_FADED_VARIABILITY": {
        colorfunc: (atom) => {
            let resi = atom.resi;
            let chain = atom.chain;
            return STRCcomputeVariabilityColor( chain, resi );
        },
        opacity: 0.5,
        thickness: 0.2,
        arrows: true
    },
    "CARTOON_HIGHLIGHT_VARIABILITY": {
        colorfunc: (atom) => {
            let resi = atom.resi;
            let chain = atom.chain;
            return STRCcomputeVariabilityColor( chain, resi );
        },
        opacity: 0.9,
        thickness: 0.2,
        arrows: true
    },
    "CARTOON_FADED_HYDROPATHICITY": {
        colorfunc: (atom) => {
            let resi = atom.resi;
            let chain = atom.chain;
            return STRCcomputeScaledKDHydropColor( chain, resi );
        },
        opacity: 0.5,
        thickness: 0.2,
        arrows: true
    },
    "CARTOON_HIGHLIGHT_HYDROPATHICITY": {
        colorfunc: (atom) => {
            let resi = atom.resi;
            let chain = atom.chain;
            return STRCcomputeScaledKDHydropColor( chain, resi );
        },
        opacity: 0.9,
        thickness: 0.2,
        arrows: true
    }
};
// Color specification of all components used for visualizations.
var COLORS = {
    // Nucleotides
    "A": "#DB5461",
    "C": "#0471A6",
    "G": "#E8C547",
    "T": "#3E885B",
    "N": "#75551d",
    "-": "#212121",
    // Secondary Structure
    "coil": "#c9c9c9",
    "helix": "#d13917",
    "sheet": "#1771d1",
    // Amino Acids
    "HIS": "#27269A",   // Polar (positive), Basic
    "LYS": "#2F49B1",
    "ARG": "#3874C8",
    "ASP": "#F75050",   // Polar (negative), Acidic
    "GLU": "#FB6581",
    "SER": "#1EAB13",   // Polar (neutral)
    "THR": "#38B42B",
    "ASN": "#51BC43",
    "GLN": "#6AC55B",
    "CYS": "#DEC443",   // Sulfur bridge forming
    "PHE": "#5826ED",   // Aromatic
    "TRP": "#653CEF",
    "TYR": "#7352F1",
    "ALA": "#239E90",   // Aliphatic
    "VAL": "#34A597",
    "LEU": "#45AD9E",
    "ILE": "#55B4A5",
    "MET": "#66BCAC",
    "PRO": "#77C3B4",
    "GLY": "#88CABC",
    "TER": "#E61090",   // Termination
    "ANY": "#75551D",   // Unknown amino acid
    "INC": "#75551D",   // Incomplete amino acid.
    "DEL": "#212121",   // Deletion
    "NONE": "#E4E5ED"
};
// Encoding of string contents, i.e. nucleotides and amino-acids, used for internal representation in visualizations.
let c = 0;
var CONTENT_ENCODING = {}
for ( const key of Object.keys( COLORS ) ) {
    CONTENT_ENCODING[ key ] = c;
    c += 1;
}
// Decoding of string contents, i.e. nucleotides and amino-acids, used for internal representation in visualizations.
var CONTENT_DECODING = { };
for ( const [key, value] of Object.entries( CONTENT_ENCODING ) ) {
    CONTENT_DECODING[ value ] = key;
}

/*
[2] Definition of runtime variables that hold data or the state of single components.
    The single components will use the following abbreviations for variable names:
    - structureView > STRC
    - variantsView > VARC
    - compositionView > COMC
    - info > INFC
*/
// The main object to hold data.
var DATA = new Object( );
// The view height set for eChart elements.
var VARC_ECHART_VIEW_HEIGHT = 0;
var VARC_ECHART_VIEW_WIDTH = 0;
var COMC_ECHART_VIEW_HEIGHT = 0;
var COMC_ECHART_VIEW_WIDTH = 0;
window.onload = _ => {
    VARC_ECHART_VIEW_HEIGHT = ( ( window.innerHeight * ( 1 - 0.516 ) ) - 60 ).toString( ) + "px";
    VARC_ECHART_VIEW_WIDTH = ( window.innerWidth - 8 ).toString( ) + "px";
    document.getElementById('variantsViewContent').style.height = VARC_ECHART_VIEW_HEIGHT;
    document.getElementById('variantsViewContent').style.width = VARC_ECHART_VIEW_WIDTH;
    new ResizeObserver( ( _ ) => {
        VARC_ECHART_VIEW_HEIGHT = ( ( window.innerHeight * ( 1 - 0.516 ) ) - 60 ).toString( ) + "px";
        VARC_ECHART_VIEW_WIDTH = ( window.innerWidth - 8 ).toString( ) + "px";
        document.getElementById( "variantsViewContent" ).style.height = VARC_ECHART_VIEW_HEIGHT;
        document.getElementById( "variantsViewContent" ).style.width = VARC_ECHART_VIEW_WIDTH;
        if ( typeof( VARC_ECHART ) !== "undefined" ) {
            VARC_ECHART.resize( );
        }
    } ).observe( document.getElementById( "main" ) );

    COMC_ECHART_VIEW_HEIGHT = ( window.innerHeight * 0.513 ).toString( ) + "px";
    COMC_ECHART_VIEW_WIDTH = ( window.innerWidth * 0.3 ).toString( ) + "px";
    document.getElementById('compositionViewContent').style.height = COMC_ECHART_VIEW_HEIGHT;
    document.getElementById('compositionViewContent').style.width = COMC_ECHART_VIEW_WIDTH;
    new ResizeObserver( ( _ ) => {
        COMC_ECHART_VIEW_HEIGHT = ( window.innerHeight * 0.513 ).toString( ) + "px";
        COMC_ECHART_VIEW_WIDTH = ( window.innerWidth * 0.3 ).toString( ) + "px";
        document.getElementById( "compositionViewContent" ).style.height = COMC_ECHART_VIEW_HEIGHT;
        document.getElementById( "compositionViewContent" ).style.width = COMC_ECHART_VIEW_WIDTH;
        if ( typeof( COMC_ECHART ) !== "undefined" ) {
            COMC_ECHART.resize( );
        }
    } ).observe( document.getElementById( "main" ) );
};
// The variantViews eChart element.
var VARC_ECHART;
// The variantViews eChart options object, i.e. the chart specification object, for more information see: https://echarts.apache.org/en/option.html#title
var VARC_ECHART_OPTION;
// The compositionView eChart element.
var COMC_ECHART;
// The compositionViews eChart options object, i.e. the chart specification object, for more information see: https://echarts.apache.org/en/option.html#title
var COMC_ECHART_OPTION;
// If all components finished initialising.
var IS_LOADED = false;
// Timeout in ms for tooltip display delay.
var TIMEOUT = 150;
// Whether the variantView is currently in merged genotypes mode.
var VARC_GENTYPES_MERGED = true;
// The names of present samples.
var SAMPLE_LABELS = [ ];
// List of available chains.
var CHAIN_IDENTIFIERS = [ ];
// The currently selected chain.
var SELECTED_CHAIN = "";
// The .pdb file as File object.
var PDB_FILE;
// The .pdb files content as string.
var PDB_STRING;
// The 3DMol structure viewer.
var STRUCTURE_VIEWER;
// JS object containing per residue information, e.g. secondary-structure and hydropathicity.
var PROTEIN_RESIDUES = new Object( );
// Variables used to compute scaled hydropathicity values as done by Expasy/ProtScale, for more information see: https://web.expasy.org/protscale/protscale-doc.html
var HYDROPATHICITY_WINDOW = 3;
var HYDROPATHICITY_WEIGHT = 0.5;
var HYDROPATHICITY_MIN_VAL = 0;
var HYDROPATHICITY_MAX_VAL = 1;
// Variables to store VARC heatmap data and meta-information.
var VARC_REF_PROTEIN_CONTENT = new Object( );
var VARC_REF_PROTEIN_SS = [ ];
var VARC_REF_PROTEIN_HYDROP = [ ];
var VARC_REF_PROTEIN_VARIABILITY = [ ];
var VARC_REF_TRANSLATED_FEATURE_CONTENT = new Object( );
var VARC_VARIANTS_CONTENT = new Object( );
var VARC_VARIANTS_CONTENT_AMBIGUOUS = new Object( );
var VARC_VARIANTS_LABELS = [ ];
var VARC_PROTEOFORM_ACCESSOR_MAP = new Object( );
var VARC_GENOTYPE_ACCESSOR_MAP = new Object( );
var VARC_PER_SUPERPOSITION_VARIANT_ENTRIES = new Object( );
var VARC_AGGREGATION_MODE = 'PROTEOFORM'; // SAMPLE, GENOTYPE or PROTEOFORM
var VARC_AMBIGUITY_MODE = 'MASK'; // MASK, HIDE, SHOW
// Variables to store STRC meta-information.
var STRC_HIGHLIGHTED_REGION = [ ];
var STRC_COLOR_SCHEME = "VARIABILITY";
var STRC_SELECTED_RESIDUES = [ ];
var STRC_SHOW_CONTACTS = true;
// Variables to store COMC meta-information.
var COMC_DISPLAYED_POSITION = -1;
// Variables to store INFC meta-information.
var INFC_SUBSTITUTION = "None";

var PROTEIN_VIEWER_CONTACT_LABELS = [ ];
var PV_SNV_HIGHLIGHTED = [ ];
var PV_SNV_IS_HIGHLIGHTED = false;
var PV_SNV_HIGHLIGHT_MIN = 0;
var PV_TERMINI_HIGHLIGHTED = [ ];
var PV_TERMINI_IS_HIGHLIGHTED = false;
var PV_BACKBONE_COLOR = "_variants";
var PV_SHOW_CONTACTS = false;

/* Structure View
var VARIANT_POSITIONS = new Object( );
var VARIANTS_ENTROPY = new Object( ); 
// Info Component
var INFOHasContactInfo = false;
var INFOInducedSubstitution = false;
*/

/**
 * Closes the specified element.
 * 
 * @param {Text} elementId 
 */
function closeElement( elementId ) {
    document.getElementById( elementId ) ? document.getElementById( elementId ).style.display = "none" : null;
}

/**
 * Shows the specified element with the specified display style.
 * 
 * @param {Text} elementId 
 * @param {Text} style 
 */
function showElement( elementId, style ) {
    document.getElementById( elementId ) ? document.getElementById( elementId ).style.display = style : null;
}

/**
 * Handles the user selection of a directory.
 * 
 * @param {Event} event: A filedrop event. The target of the event has to be a `FileList` instance. 
 */
 function directoryInputHandler( event ) {
    // Hide input placeholder text and show process loader.
    document.getElementById( "directoryInput" ).style.display = "none";
    document.getElementById( "directoryInputIcon" ).style.display = "none";
    document.getElementById( "directoryInputText" ).style.display = "none";
    document.getElementById( "directoryInputProcessLoader" ).style.display = "flex";
    var files = event.target.files;
    loadData( files, init );
}

/**
 * Replaces the loader with an error message.
 */
function showError( ) {
    document.getElementById( "inputComponent" ).innerHTML = '<i class="fas fa-exclamation-triangle"></i>&nbsp;&nbsp;<b>A critical error occurred during initialization!</b>';
    $('#inputComponent b').css('color','var(--error-color)');
}

/**
 * Parses the MUSIAL output files specified with the passed file list.
 * 
 * @param {FileList} fileList: A `FileList` instance that specifies the output files of a MUSIAL run.
 * @param {_callback} function: A function to call after relevant data has been loaded. The callback function has to accept exactly the fileList parameter.
 */
 function loadData( fileList, _callback ) {
    var file;
    var list = [ ];
    // Parse structure allocation data.
    window.setTimeout(
        function(){
            list = [...fileList].filter( s => s.name.startsWith("structureSuperposition") );
            if ( list.length == 0 ) {
                Swal.fire({
                    icon: 'error',
                    title: 'No structureSuperposition.json file is present in the output directory!',
                    text: 'Please check your MUSIAL output.',
                    background: '#EFF0F8',
                    confirmButtonText: 'Ok',
                    confirmButtonColor: '#607196'
                });
                showError( );
                return;
            } else {
                file = list[ 0 ];
                parseSuperposition( file );
            }
        },
        0
    );
    // Parse structure residues distance map.
    window.setTimeout(
        function(){
            list = [...fileList].filter( s => s.name.startsWith("residueDistanceMap") );
            if ( list.length == 0 ) {
                Swal.fire({
                    icon: 'error',
                    title: 'No residueDistanceMap.tsv file is present in the output directory!',
                    text: 'Please check your MUSIAL output.',
                    background: '#EFF0F8',
                    confirmButtonText: 'Ok',
                    confirmButtonColor: '#607196'
                });
                showError( );
                return;
            } else {
                file = list[ 0 ];
                parseDistanceMap( file );
            }
        },
        100
    );
    // Parse SNV table annotations.
    window.setTimeout(
        function(){
            list = [...fileList].filter( s => s.name.startsWith("variantAnnotationTable") );
            if ( list.length == 0 ) {
                Swal.fire({
                    icon: 'error',
                    title: 'No variantAnnotationTable.tsv file is present in the output directory!',
                    text: 'Please check your MUSIAL output.',
                    background: '#EFF0F8',
                    confirmButtonText: 'Ok',
                    confirmButtonColor: '#607196'
                });
                showError( );
                return;
            } else {
                file = list[ 0 ];
                parseAnnotationTable( file );
            }
        },
        200
    );
    // Run provided callback function to further process parsed data.
    window.setTimeout(
        function(){
            try{
                _callback( fileList );
            } catch (e) {
                Swal.fire({
                    icon: 'error',
                    title: 'An error occurred during the initialization!',
                    text: e.message,
                    background: '#EFF0F8',
                    confirmButtonText: 'Ok',
                    confirmButtonColor: '#607196'
                });
                showError( );
                return;
            }
        },
        300
    );
}

/**
 * Parses structure allocation data from a file.
 * 
 * @param {File} file: File object pointing to the structureAllocation.json file of a MUSIAL output directory.
 */
function parseSuperposition( file ) {
    var fileReader = new FileReader( );
    var fileContent;
    fileReader.onload = function(event) {
        fileContent = event.target.result;
        DATA = JSON.parse( fileContent );
    };
    fileReader.readAsText( file );
}

/**
 * Parses structure residues distance data from a file.
 * 
 * @param {File} file: File object pointing to the residueDistanceMap.tsv file of a MUSIAL output directory.
 */
function parseDistanceMap( file ) {
    var fileReader = new FileReader( );
    var fileContent;
    DATA[ "ProteinResidueDistanceMap" ] = new Object( );
    fileReader.onload = function(event) {
        fileContent = event.target.result;
        var results = fileContent.split( "\r\n" );
        results.shift( );
        var residues = results.shift( ).split( "\t" );
        residues.shift( );
        for ( var residue of residues ) {
            DATA[ "ProteinResidueDistanceMap" ][ residue ] = [ ];
        }
        for ( let r of results ) {
            let row = r.split( "\t" );
            if ( row.length > 1 ) {
                var residue = row[ 0 ];
                row.shift( );
                for ( var index = 0; index < row.length; index++ ) {
                    if ( row[ index ] <= 5.0 && residues[ index ] !== residue ) {
                        DATA[ "ProteinResidueDistanceMap" ][ residues[ index ] ].push( [ residue, row[ index ] ] );
                    }
                }
            }
        }
    };
    fileReader.readAsText( file );
}

/**
 * Parses per sample per variant position annotations from a file.
 * 
 * @param {File} file: File object pointing to the variantAnnotationTable.tsv file of a MUSIAL output directory.
 */
function parseAnnotationTable( file ) {
    var fileReader = new FileReader( );
    var fileContent;
    DATA[ "VariantAnnotationTable" ] = new Object( );
    fileReader.onload = function(event) {
        fileContent = event.target.result;
        var results = fileContent.split( "\r\n" );
        results.shift( );
        var header = results.shift( ).split( "\t" );
        header.shift( );
        header.shift( );
        for ( let r of results ) {
            let row = r.split( "\t" );
            if ( row.length > 1 ) {
                var position = row[ 0 ];
                // Multiply position with -1 if reference feature is located on the anti-sense strand and parse only one annotation for insertions.
                if ( position.includes("+") ) {
                    if ( position.split( "+" )[ 1 ] === '1' ) {
                        if ( DATA.GeneReferenceIsSense === 'false' ) {
                            position = ( - parseInt( position.split( "+" )[ 0 ] ) ).toString( ) + 'I';
                        } else {
                            position = position.split( "+" )[ 0 ] + 'I';
                        }
                    } else {
                        continue;
                    }
                } else if ( DATA.GeneReferenceIsSense === 'false' ) {
                    position = ( - parseInt( position.split( "+" )[ 0 ] ) ).toString( );
                }
                DATA[ "VariantAnnotationTable" ][ position ] = { }; 
                row.shift( );
                row.shift( );
                for ( var index = 0; index < row.length; index++ ) {
                    DATA[ "VariantAnnotationTable" ][ position ][ header[ index ] ] = row[ index ];
                }
            }
        }
    };
    fileReader.readAsText( file );
}

/**
 * Initializes all visualizations from the passed file list and data object.
 * 
 * @param {FileList} files: List of file references to MUSIAL output files.
 */
 function init( files ) {
     // Access and store list of available chains of the protein structure.
     let chainId;
     for ( let key of Object.keys( DATA ) ) {
         if ( key.startsWith( "SuperpositionChain" ) ) {
             chainId = key.replace( "SuperpositionChain", "" );
             CHAIN_IDENTIFIERS.push( chainId );
         }
     }
     if ( CHAIN_IDENTIFIERS.length == 0 ) {
        // If no chains are available, show an error popup.
        Swal.fire({
            icon: 'error',
            title: 'No chain data available!',
            text: 'Please check your MUSIAL output and used .pdb file.',
            background: '#EFF0F8',
            confirmButtonText: 'Ok',
            confirmButtonColor: '#607196'
        });
        showError( );
        return;
     } else {
        // Else the first chain is selected.
        SELECTED_CHAIN = CHAIN_IDENTIFIERS[ 0 ];
     }
    // Access .pdb file from specified files.
    for ( let file of files ) {
        if ( file.name.endsWith(".pdb") ) {
            PDB_FILE = file;
        }
    }
    // Initialize the structure view component.
    window.setTimeout( function(){ initSTRC( PDB_FILE ); }, 0 );
    // Initialize the composition summary component.
    window.setTimeout( function(){ initCOMC( ); }, 250 );
    // Initialize alignment view component.
    window.setTimeout( function(){ initVARC( VARC_AGGREGATION_MODE ); }, 300 );
    // Show all components.
    window.setTimeout( function(){
        closeElement( "inputComponent" );
        showElement( "structureViewComponent", "block" );
        showElement("compositionViewComponent", "block");
        showElement("infoComponent", "block");
        showElement("variantsViewComponent", "block");
        },
        500
    );
    IS_LOADED = true;
    // console.log( DATA );
}

/**
 * Initializes the structure view component.
 * 
 * @param {File} pdbFile: A file object pointing to a .pdb file.
 */
 function initSTRC( pdbFile ) {
     // In order to mark positions VARC variant labels have to be computed here?
     switch ( VARC_AGGREGATION_MODE ) {
        case 'SAMPLE':
            for ( let proteoform of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ] ) ) {
                for ( let genotype of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ][ proteoform ][ "Genotypes" ] ) ) {
                    for ( let sample of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Genotypes" ][ genotype ][ "Samples" ] ) ) {
                        VARC_VARIANTS_LABELS.push( sample );
                    }
                }
            }
            break;
        case 'GENOTYPE':
            for ( let proteoform of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ] ) ) {
                for ( let genotype of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ][ proteoform ][ "Genotypes" ] ) ) {
                    VARC_VARIANTS_LABELS.push( genotype );
                }
            }
            break;
        case 'PROTEOFORM':
            for ( let proteoform of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ] ) ) {
                VARC_VARIANTS_LABELS.push( proteoform );
            }
            break;
    }
    let fileReader = new FileReader( );
    fileReader.readAsText( pdbFile );
    fileReader.onload = function(e) {
        PDB_STRING = fileReader.result;
        let element = $('#structureViewContent');
        let config = { backgroundColor: '#fafafc', id: 'structureViewCanvas', antialias: true, cartoonQuality: 6 };
        STRUCTURE_VIEWER = $3Dmol.createViewer( element, config );
        STRUCTURE_VIEWER.addModels( PDB_STRING, "pdb" );
        STRUCTURE_VIEWER.setClickable(
            {},
            true,
            (sel) => {
                if ( sel.chain != SELECTED_CHAIN ) {
                    VARCSelectChain( sel.chain );
                }
                STRCSelectResidue( sel.chain, sel.resi, sel.resn, true, COLORS[ sel.resn ], COLORS[ sel.resn ], 0.7, 1 );
                VARCMarkPosition( PROTEIN_RESIDUES[ sel.chain ][ sel.resi ][ "superposition" ] );
                STRCDisplayContacts( );
                INFCSetReferenceInformation( PROTEIN_RESIDUES[ sel.chain ][ sel.resi ][ "superposition" ] );
                COMCDisplay( PROTEIN_RESIDUES[ sel.chain ][ sel.resi ][ "superposition" ] );
            }
        );
        STRUCTURE_VIEWER.zoomTo();
        STRUCTURE_VIEWER.render();
        STRUCTURE_VIEWER.zoom(0.95, TIMEOUT);
        // Gather list of residues for extended computations from 3DMol.js GLModel
        let residueIds = new Set( );
        for ( const atom of STRUCTURE_VIEWER.getInternalState( ).models[ 0 ].atoms ) {
            let resi = atom.resi;
            let chain = atom.chain;
            let ss = atom.ss;
            if ( !(chain in PROTEIN_RESIDUES) ) {
                PROTEIN_RESIDUES[ chain ] = new Object( ); 
            }
            if ( residueIds.has( chain + resi ) ) {
                continue;
            } else {
                switch ( ss ) {
                    case 'c':
                        ss = 'coil';
                        break;
                    case 's':
                        ss = 'sheet';
                        break;
                    case 'h':
                        ss = 'helix';
                        break;
                }
                PROTEIN_RESIDUES[ chain ][ resi ] = { "resn": atom.resn, "ss": ss };
                residueIds.add( chain + resi );
            }
        }
        // Add superpositions to residues.
        var proteinPosition;
        for ( const chain of CHAIN_IDENTIFIERS ) {
            for (  const superposition of Object.keys( DATA[ "SuperpositionChain" + chain ] ) ) {
                if ( superposition === "ReferenceProteinLength" || superposition === "Genotypes" || superposition === "Proteoforms" ) {
                    continue;
                } else {
                    proteinPosition = DATA[ "SuperpositionChain" + chain ][ superposition ][ "ProteinPosition" ];
                    if ( proteinPosition !== -1 ) {
                        PROTEIN_RESIDUES[ chain ][ proteinPosition ][ "superposition" ] = superposition;
                    }
                }
            }
        }
        // Compute hydropathicity for collected residues.
        STRCcomputeScaledKDHydrop( );
        // Compute uncertainty for collected residues.
        STRCcomputeVariability( );
    };
}

/**
 * Computes uncertainty scores for each residue stored in `PROTEIN_RESIDUES` as the information entropy using the logarithm to the base 23 (as residues can have 23 types: The 20 naturally occuring residues, termination, any amino acid and ambiguous amino acid).
 * Residues with a low uncertainty, i.e. low variability, will achieve values close to zero and residues with a high uncertainty, i.e. high variability will achieve values close to one.
 */
function STRCcomputeVariability( ) {
    let noSamples = 0;
    let proteinPosition;
    let positionContentCounts = new Object( );
    let proteoformContent;
    let variability;
    for ( let chain of Object.keys( PROTEIN_RESIDUES ) ) {
        // First compute number of samples per proteoform.
        let perProteoformNoSamples = new Object( ); 
        for ( let proteoform of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ] ) ) {
            perProteoformNoSamples[ proteoform ] = 0;
            for ( let genotype of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ][ proteoform ][ "Genotypes" ] ) ) {
                for ( let sample of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Genotypes" ][ genotype ][ "Samples" ] ) ) {
                    perProteoformNoSamples[ proteoform ] += 1;
                }
            }
        }
        // Iterate over each superimposed position to access amino-acid counts.
        for ( let index of Object.keys( DATA[ "SuperpositionChain" + chain ] ) ) {
            noSamples = 0;
            variability = 0;
            positionContentCounts = new Object( );
            if ( index === "ReferenceProteinLength" || index === "Genotypes" || index === "Proteoforms" ) {
                continue;
            }
            proteinPosition = DATA[ "SuperpositionChain" + chain ][ index ][ "ProteinPosition" ];
            if ( proteinPosition !== -1 ) {
                // Count residue types per proteoform and increase sample count wrt. ambiguity mode.
                for ( let proteoform of Object.keys( perProteoformNoSamples ) ) {
                    proteoformContent = DATA[ "SuperpositionChain" + chain ][ index ][ "PerProteoformAminoAcidContent" ][ proteoform ];
                    if ( VARC_AMBIGUITY_MODE != 'SHOW' && proteoformContent === proteoformContent.toLowerCase( ) ) {
                        continue;
                    } else {
                        proteoformContent = proteoformContent.toUpperCase( );
                        if ( proteoformContent === "INC" ) {
                            proteoformContent = "ANY";
                        }
                        if ( proteoformContent in positionContentCounts ) {
                            positionContentCounts[ proteoformContent ] += perProteoformNoSamples[ proteoform ];
                            noSamples += perProteoformNoSamples[ proteoform ];
                        } else {
                            positionContentCounts[ proteoformContent ] = perProteoformNoSamples[ proteoform ];
                            noSamples += perProteoformNoSamples[ proteoform ];
                        }
                    }
                }
                for ( let value of Object.values( positionContentCounts ) ) {
                    variability += ( value / noSamples ) * ( Math.log( value / noSamples ) / Math.log( 23 ) );  
                }
                PROTEIN_RESIDUES[ chain ][ proteinPosition ][ "variability" ] = parseFloat( -1 * variability ).toFixed( 2 );
            }
        }
    }
}

/**
 * Computes a color based on a uncertainty score.
 * The color turns to a blue-green tone the more uncertain the corresponding residue position is.
 */
function STRCcomputeVariabilityColor( chain, resi ) {
    function STRCcomputeColorDiff( to, from, val ) {
        return parseFloat( from + ( ( to - from ) * val ) ).toFixed( 0 );
    }
    let variabilityVal = PROTEIN_RESIDUES[ chain ][ resi ][ "variability" ];
    if ( variabilityVal === 0 ) {
        return "rgb(200,200,200)";
    } else if ( variabilityVal < 0.4  ) {
        return "rgb(" + STRCcomputeColorDiff( 255, 200, variabilityVal / 0.4 ) + ",200," + STRCcomputeColorDiff( 55, 200, variabilityVal / 0.4 ) + ")";
    } else if ( variabilityVal === 0.4  ) {
        return "rgb(255,200,55)";
    } else if (variabilityVal === 1) {
        return "rgb(255,0,55)";
    } else {
        return "rgb(255" + "," + STRCcomputeColorDiff( 0, 200, ( variabilityVal - 0.4 ) / 0.6 ) + ",55)";
    }
}

/**
 * Computes scaled Kyte-Doolittle Hyropathicity scores for each residue stored in `PROTEIN_RESIDUES` using the window and edge weights stored in `HYDROPATHICITY_WINDOW` and `HYDROPATHICITY_WEIGHT` respectively.
 * Hydrophobic residues will achieve a value of up to 1 while hydrophilic values will achieve a value of down to 0.
 */
function STRCcomputeScaledKDHydrop( ) {
    let shift = ( ( HYDROPATHICITY_WINDOW - 1 ) / 2 );
    let normFactor;
    let f;
    let hydropathicity;
    var residueIndex;
    var residueIndexLeft;
    var residueIndexRight;
    var residueIndices;
    for ( let chain of Object.keys( PROTEIN_RESIDUES ) ) {
        residueIndices = Object.keys( PROTEIN_RESIDUES[ chain ] );
        let firstIndex = Math.min( ...residueIndices );
        let lastIndex = Math.max( ...residueIndices );
        for ( residueIndex of residueIndices ) {
            hydropathicity = 0;
            normFactor = 0;
            for ( residueIndexLeft = residueIndex - 1; residueIndexLeft >= Math.max( firstIndex, residueIndex - shift ); residueIndexLeft-- ) {
                f = ( HYDROPATHICITY_WEIGHT + ( ( ( 1 - HYDROPATHICITY_WEIGHT ) / shift ) * ( shift - Math.abs( residueIndex - residueIndexLeft ) ) ) );
                hydropathicity += SCALED_KD_HYDROP[ PROTEIN_RESIDUES[ chain ][ residueIndexLeft ].resn ] * f;
                normFactor += f;
            }
            f = 1;
            hydropathicity += SCALED_KD_HYDROP[ PROTEIN_RESIDUES[ chain ][ residueIndex ].resn ] * f;
            normFactor += 1;
            for ( residueIndexRight = residueIndex + 1; residueIndexRight <= Math.min( lastIndex, residueIndex + shift ); residueIndexRight++ ) {
                f = ( HYDROPATHICITY_WEIGHT + ( ( ( 1 - HYDROPATHICITY_WEIGHT ) / shift ) * ( shift - Math.abs( residueIndex - residueIndexRight ) ) ) );
                hydropathicity += SCALED_KD_HYDROP[ PROTEIN_RESIDUES[ chain ][ residueIndexRight ].resn ] * f;
                normFactor += f;
            }
            PROTEIN_RESIDUES[ chain ][ residueIndex ][ "scKDH" ] = parseFloat( hydropathicity / normFactor ).toFixed( 2 );
        }
    }
}

/**
 * Computes a color based on a scaled Kyte-Doolittle Hyropathicity score.
 * Hydrophobic residues will be assigned a brownish color while hydrophilic residues will be assigned a blueish color.
 */
function STRCcomputeScaledKDHydropColor( chain, resi ) {
    function STRCcomputeColorDiff( to, from, val ) {
        return parseFloat( from + ( ( to - from ) * val ) ).toFixed( 0 );
    }
    let scKDHValue = PROTEIN_RESIDUES[ chain ][ resi ][ "scKDH" ];
    if ( scKDHValue > 0.5 ) {
        return "rgb(" + STRCcomputeColorDiff( 184, 230, ( scKDHValue - 0.5 ) / 0.5 ) + "," + STRCcomputeColorDiff( 119, 230, ( scKDHValue - 0.5 ) / 0.5 ) + "," + STRCcomputeColorDiff( 55, 230, ( scKDHValue - 0.5 ) / 0.5 ) + ")";
    } else if ( scKDHValue < 0.5 ) {
        return "rgb(" + STRCcomputeColorDiff( 69, 230, scKDHValue / 0.5 ) + "," + STRCcomputeColorDiff( 167, 230, scKDHValue / 0.5 ) + "," + STRCcomputeColorDiff( 237, 230, scKDHValue / 0.5 ) + ")";
    } else {
        return "rgb(230, 230, 230)"
    }
}

/**
 * 
 */
function STRCClearSelection( ) {
    // Remove existing selection.
    if ( STRC_SELECTED_RESIDUES.length != 0 ) {
        for ( let selectedResidue of STRC_SELECTED_RESIDUES ) {
            if ( selectedResidue !== undefined ) {
                if ( ( STRC_HIGHLIGHTED_REGION[ 0 ] <= selectedResidue[ 1 ] ) && ( selectedResidue[ 1 ] <= STRC_HIGHLIGHTED_REGION[ 1 ] ) ) {
                    STRUCTURE_VIEWER.setStyle(
                        {
                            chain: selectedResidue[ 0 ],
                            resi: selectedResidue[ 1 ]
                        },
                        {
                            cartoon: STYLES_3DMOL[ "CARTOON_HIGHLIGHT_" + STRC_COLOR_SCHEME ],
                            stick: STYLES_3DMOL.LINE_FADED
                        }
                    );
                } else {
                    STRUCTURE_VIEWER.setStyle(
                        {
                            chain: selectedResidue[ 0 ],
                            resi: selectedResidue[ 1 ]
                        },
                        {
                            cartoon: STYLES_3DMOL[ "CARTOON_FADED_" + STRC_COLOR_SCHEME ],
                            stick: STYLES_3DMOL.LINE_FADED
                        }
                    );
                }
                STRUCTURE_VIEWER.removeLabel( selectedResidue[ 3 ] );
            }
        }
        STRC_SELECTED_RESIDUES = [ ];
    }
}

/**
 * 
 * @param {*} chain 
 * @param {*} resi 
 * @param {*} resn 
 * @param {*} clear 
 * @param {*} labelBackgroundColor 
 * @param {*} labelBorderColor 
 * @param {*} backgroundOpacity
 * @param {*} borderThickness
 * @returns 
 */
function STRCSelectResidue( chain, resi, resn, clear, labelBackgroundColor, labelBorderColor, backgroundOpacity, borderThickness ) {
    // Remove existing selection if specified.
    if ( clear ) {
        STRCClearSelection( );
        STRUCTURE_VIEWER.center( { chain: chain, resi: resi }, 300 );
    }
    STRUCTURE_VIEWER.addStyle(
        {
            chain: chain,
            resi: resi
        },
        {
            stick: STYLES_3DMOL.LINE_HIGHLIGHT
        }
    );
    var label = STRUCTURE_VIEWER.addLabel(
        chain + " " + resi + ": " + resn,
        {
            backgroundColor: labelBackgroundColor,
            backgroundOpacity: backgroundOpacity,
            borderColor: labelBorderColor,
            borderThickness: borderThickness,
            fontColor: "#292926",
            fontSize: 10,
            alignment: "topLeft"
        },
        {
            chain: chain,
            resi: resi
        },
        false
    );
    STRC_SELECTED_RESIDUES.push( [ chain, resi, resn, label ] );
    STRUCTURE_VIEWER.render( );
}

/**
 * Sets the structure backbone coloring scheme to the specified scheme and applies the style.
 * 
 * @param {String} scheme One of `VARIABILITY` and `HYDROP` is currently supported.
 */
function STRCSetColorScheme( scheme ) {
    let colorVariabilityBtn = document.getElementById( 'menuSTRCColorVariabilityBtn' );
    let colorHydropathicityBtn = document.getElementById( 'menuSTRCColorHydropathicityBtn' );
    switch (scheme) {
        case 'VARIABILITY':
            if (STRC_COLOR_SCHEME != 'VARIABILITY'){
                STRC_COLOR_SCHEME = 'VARIABILITY';
                colorVariabilityBtn.classList.add( "menuBtnActive" );
                colorHydropathicityBtn.classList.remove( "menuBtnActive" );
            }
            break;
        case 'HYDROPATHICITY':
            if (STRC_COLOR_SCHEME != 'HYDROPATHICITY'){
                STRC_COLOR_SCHEME = 'HYDROPATHICITY';
                colorHydropathicityBtn.classList.add( "menuBtnActive" );
                colorVariabilityBtn.classList.remove( "menuBtnActive" );
            }
            break;
    }
    STRCApplyStyle( );
}

/**
 * Applies the currently specified structure model style.
 */
function STRCApplyStyle( ) {
    let highlightChain = SELECTED_CHAIN;
    let highlightFrom = STRC_HIGHLIGHTED_REGION[ 0 ];
    let highlightTo = STRC_HIGHLIGHTED_REGION[ 1 ];
    // Apply default style.
    STRUCTURE_VIEWER.setStyle(
        { },
        { cartoon: STYLES_3DMOL[ "CARTOON_FADED_" + STRC_COLOR_SCHEME ], stick: STYLES_3DMOL.LINE_FADED }
    );
    // Apply membrane style.
    STRUCTURE_VIEWER.setStyle(
        {
            chain: "x",
            resn: "DUM"
        },
        {
            cross: {
                radius: 0.4,
                color: '#4d4e4f'
            }
        }
    );
    // Apply backbone highlighting.
    STRUCTURE_VIEWER.setStyle(
        {
            chain: highlightChain,
            resi: highlightFrom.toString( ) + "-" + highlightTo.toString( )
        },
        {
            cartoon: STYLES_3DMOL[ "CARTOON_HIGHLIGHT_" + STRC_COLOR_SCHEME ],
            stick: STYLES_3DMOL.LINE_FADED
        }
    );
    // Apply selection style.
    for (  let selectedResidue of STRC_SELECTED_RESIDUES ) {
        if ( selectedResidue !== undefined ) {
            STRUCTURE_VIEWER.addStyle(
                { chain: selectedResidue[ 0 ], resi: selectedResidue[ 1 ] },
                { stick: STYLES_3DMOL.LINE_HIGHLIGHT }
            );
        }
    }
    STRUCTURE_VIEWER.render( );
}

/**
 * 
 */
function STRCDisplayContacts( ) {
    if ( STRC_SELECTED_RESIDUES.length != 0 && STRC_SHOW_CONTACTS ) {
        let selResChain = STRC_SELECTED_RESIDUES[ 0 ][ 0 ]
        let selResResi = STRC_SELECTED_RESIDUES[ 0 ][ 1 ]
        let selResType = STRC_SELECTED_RESIDUES[ 0 ][ 2 ]
        let selResStruc = PROTEIN_RESIDUES[ selResChain ][ selResResi ][ "ss" ];
        selResStruc = selResStruc.charAt(0).toUpperCase( ) + selResStruc.slice(1);
        let resDistanceMapkey = selResResi + selResChain;
        if ( resDistanceMapkey in DATA[ "ProteinResidueDistanceMap" ] ) {
            let residuesInContact = DATA[ "ProteinResidueDistanceMap" ][ resDistanceMapkey ];
            if ( residuesInContact.length > 0 ) {
                let contactInformationHtmlString = "";
                for ( let residueInContact of residuesInContact ) {
                    let conResChain = residueInContact[ 0 ].slice( -1 );
                    let conResResi = residueInContact[ 0 ].slice( 0, -1 );
                    let conResType = PROTEIN_RESIDUES[ conResChain ][ conResResi ][ "resn" ];
                    let conResStruc = PROTEIN_RESIDUES[ conResChain ][ conResResi ][ "ss" ];
                    let conResColor;
                    let conResStrucColor;
                    let distance = residueInContact[ 1 ];
                    conResColor = COLORS[ conResType ];
                    switch (conResStruc) {
                        case "sheet" :
                            conResStrucColor = COLORS[ "sheet" ];
                            break;
                        case "helix" :
                            conResStrucColor = COLORS[ "helix" ];
                            break;
                        case "coil" :
                            conResStrucColor = COLORS[ "coil" ];
                            break;
                    }
                    conResStruc = conResStruc.charAt(0).toUpperCase( ) + conResStruc.slice(1);
                    STRCSelectResidue( conResChain, conResResi, conResType, false, "#E4E5ED", "", 0.5, 0.0 );
                    contactInformationHtmlString += "<br>&bull; " + conResChain + conResResi + " | " + "<span style='color: " + conResColor + "'>&#9724; </span> " + conResType + " | " + "<span style='color: " + conResStrucColor + "'>&#9724; </span>" + conResStruc + "</i><br>";
                    contactInformationHtmlString += "Distance (SCM): " + parseFloat( distance ).toFixed( 2 ) + " Å<br>";
                    /* Assing PRI15 scores */
                    let PRISMEM15WT = PRISMEM15_SCORES[ selResType + "-" + selResStruc + "/" + conResType + "-" + conResStruc ];
                    if ( PRISMEM15WT == "-Infinity" ) {
                        PRISMEM15WT = "Not observed."
                    } else {
                        PRISMEM15WT = parseInt( PRISMEM15WT );
                    }
                    contactInformationHtmlString += "Score Wildtype (PRISMEM15): " + PRISMEM15WT + "<br>";
                    if ( INFC_SUBSTITUTION != "None" && INFC_SUBSTITUTION != "Termination" && INFC_SUBSTITUTION != "Ambiguous" ) {
                        let PRISMEM15MT = PRISMEM15_SCORES[ INFC_SUBSTITUTION + "-" + selResStruc + "/" + conResType + "-" + conResStruc ];
                        if ( PRISMEM15MT == "-Infinity" ) {
                            PRISMEM15MT = "Not observed."
                        } else {
                            PRISMEM15MT = parseInt( PRISMEM15MT );
                        }
                        contactInformationHtmlString += "Score Mutation (PRISMEM15): " + PRISMEM15MT + "<br>";
                    }
                }
                STRUCTURE_VIEWER.render();
                INFCSetContactInformation( contactInformationHtmlString );
            } else {
                INFCSetContactInformation( "<br><b style='color: #cbd0e0'>No contacts.</b><br>" );
            }
        }
    } else {
        // Remove all contact information.
        for ( let i = 1; i < STRC_SELECTED_RESIDUES.length; i++ ) {
            STRUCTURE_VIEWER.removeLabel( STRC_SELECTED_RESIDUES[ i ][ 3 ] );
        }
        if ( STRC_SELECTED_RESIDUES.length !== 0 ) {
            STRC_SELECTED_RESIDUES = [ STRC_SELECTED_RESIDUES[ 0 ] ];
        } else {
            STRC_SELECTED_RESIDUES = [ ];
        }
        STRCApplyStyle( );
        if ( STRC_SHOW_CONTACTS ) {
            INFCSetContactInformation( "<br><b style='color: #cbd0e0'>No information available.</b><br>" );
        } else {
            INFCSetContactInformation( "<b style='color: #E4E5ED'>Display of contact information disabled.</b>" );
        }
    }
}

/**
 * 
 */
function STRCToggleDisplayContacts( ) {
    let STRCDisplayContactsBtn = document.getElementById( "menuSTRCDisplayContactsBtn" );
    if ( STRC_SHOW_CONTACTS ) {
        STRC_SHOW_CONTACTS = false;
        STRCDisplayContacts( );
        STRCDisplayContactsBtn.classList.remove( "menuBtnActive" );
    } else {
        STRC_SHOW_CONTACTS = true;
        STRCDisplayContacts( );
        STRCDisplayContactsBtn.classList.add( "menuBtnActive" );
    }
}

/**
 * Initializes the position composition component.
 */
 function initCOMC( ) {
    var chartDom = document.getElementById('compositionViewContent');
    COMC_ECHART = echarts.init(chartDom, { "renderer": "canvas" } );
}

/**
 * Initializes the variants view component.
 */
 function initVARC( ) {
    var chartDom = document.getElementById('variantsViewContent');
    VARC_ECHART = echarts.init(chartDom, { "renderer": "canvas" } );
    // Empty lists holding possibly existent VARC data.
    VARC_REF_PROTEIN_CONTENT = {
        "AA_CLUSTER_POLAR": [ ],
        "AA_CLUSTER_UNPOLAR": [ ],
        "AA_CLUSTER_OTHER": [ ]
    };
    VARC_REF_PROTEIN_SS = [ ];
    VARC_REF_PROTEIN_HYDROP = [ ];
    VARC_REF_PROTEIN_VARIABILITY = [ ];
    VARC_REF_TRANSLATED_FEATURE_CONTENT = {
        "AA_CLUSTER_POLAR": [ ],
        "AA_CLUSTER_UNPOLAR": [ ],
        "AA_CLUSTER_OTHER": [ ]
    };
    VARC_VARIANTS_CONTENT = {
        "AA_CLUSTER_POLAR": [ ],
        "AA_CLUSTER_UNPOLAR": [ ],
        "AA_CLUSTER_OTHER": [ ]
    };
    VARC_VARIANTS_CONTENT_AMBIGUOUS = {
        "AA_CLUSTER_POLAR": [ ],
        "AA_CLUSTER_UNPOLAR": [ ],
        "AA_CLUSTER_OTHER": [ ]
    };
    VARC_VARIANTS_LABELS = [ ];
    VARC_PROTEOFORM_ACCESSOR_MAP = new Object( );
    VARC_GENOTYPE_ACCESSOR_MAP = new Object( );
    // Dependent on the specified mode fill variants data labels.
    switch ( VARC_AGGREGATION_MODE ) {
        case 'SAMPLE':
            for ( let proteoform of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ] ) ) {
                for ( let genotype of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ][ proteoform ][ "Genotypes" ] ) ) {
                    for ( let sample of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Genotypes" ][ genotype ][ "Samples" ] ) ) {
                        VARC_VARIANTS_LABELS.push( sample );
                        VARC_PROTEOFORM_ACCESSOR_MAP[ sample ] = proteoform;
                        VARC_GENOTYPE_ACCESSOR_MAP[ sample ] = genotype;
                    }
                }
            }
            break;
        case 'GENOTYPE':
            for ( let proteoform of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ] ) ) {
                for ( let genotype of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ][ proteoform ][ "Genotypes" ] ) ) {
                    VARC_VARIANTS_LABELS.push( genotype );
                    VARC_PROTEOFORM_ACCESSOR_MAP[ genotype ] = proteoform;
                    VARC_GENOTYPE_ACCESSOR_MAP[ genotype ] = genotype;
                }
            }
            break;
        case 'PROTEOFORM':
            for ( let proteoform of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ] ) ) {
                VARC_VARIANTS_LABELS.push( proteoform );
                VARC_PROTEOFORM_ACCESSOR_MAP[ proteoform ] = proteoform;
            }
            break;
    }

    // Iterate over protein positions of selected chain and fill 2D-lists (i.e. heatmap data matrices).
    var c = 0;
    let proteinPosition;
    let refProteinContent;
    let refTransNucContent;
    let variantContent;
    var superpositions = [ ];
    var variantsEditedLabels = [ ];
    var variantYAxisName = 'Samples';
    var referenceCategoricalDataLabels = [ ];
    var referenceCategoricalsCount = 3;
    var referenceNumericalDataLabels = [ ];
    var referenceNumericalsCount = 2;
    var variantPositionsIndicatorData = [ ];
    for( let superposition of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ] ) ) {
        if ( superposition === 'ReferenceProteinLength' || superposition === 'Genotypes' || superposition === 'Proteoforms' ) {
            continue;
        } else {
            proteinPosition = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ superposition ][ "ProteinPosition" ];
        }
        refProteinContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ superposition ][ "ProteinContent" ];
        if ( refProteinContent !== "None" && proteinPosition !== -1 ) {
            refProteinContent = CONTENT_ENCODING[ refProteinContent  ];
            if ( refProteinContent > 28 ) {
                VARC_REF_PROTEIN_CONTENT[ "AA_CLUSTER_OTHER" ].push( [ c, 1, refProteinContent ] );
            } else if ( refProteinContent > 18 ) {
                VARC_REF_PROTEIN_CONTENT[ "AA_CLUSTER_UNPOLAR" ].push( [ c, 1, refProteinContent ] );
            } else if ( refProteinContent > 8 ) {
                VARC_REF_PROTEIN_CONTENT[ "AA_CLUSTER_POLAR" ].push( [ c, 1, refProteinContent ] );
            }
            VARC_REF_PROTEIN_SS.push(
                [ c, 0, CONTENT_ENCODING[ PROTEIN_RESIDUES[ SELECTED_CHAIN ][ proteinPosition ][ "ss" ] ] ]
            );
            VARC_REF_PROTEIN_HYDROP.push(
                [ c, 0, PROTEIN_RESIDUES[ SELECTED_CHAIN ][ proteinPosition ][ "scKDH" ] ]
            );
            VARC_REF_PROTEIN_VARIABILITY.push(
                [ c, 1, PROTEIN_RESIDUES[ SELECTED_CHAIN ][ proteinPosition ][ "variability" ] ]
            );
        }
        refTransNucContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ superposition ][ "ReferenceFeatureTranslatedContent" ];
        if ( refTransNucContent !== "None" ) {
            refTransNucContent = CONTENT_ENCODING[ refTransNucContent  ];
            if ( refTransNucContent > 28 ) {
                VARC_REF_TRANSLATED_FEATURE_CONTENT[ "AA_CLUSTER_OTHER" ].push( [ c, 2, refTransNucContent ] );
            } else if ( refTransNucContent > 18 ) {
                VARC_REF_TRANSLATED_FEATURE_CONTENT[ "AA_CLUSTER_UNPOLAR" ].push( [ c, 2, refTransNucContent ] );
            } else if ( refTransNucContent > 8 ) {
                VARC_REF_TRANSLATED_FEATURE_CONTENT[ "AA_CLUSTER_POLAR" ].push( [ c, 2, refTransNucContent ] );
            }
        }
        let hasUnambiguousVariant = false;
        let labelIndex = 0;
        for ( let label of VARC_VARIANTS_LABELS ) {
            variantContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ superposition ][ "PerProteoformAminoAcidContent" ][ VARC_PROTEOFORM_ACCESSOR_MAP[ label ] ];
            labelIndex = VARC_VARIANTS_LABELS.indexOf( label );
            if ( CONTENT_ENCODING[ variantContent ] !== refTransNucContent ) {
                if ( variantContent == variantContent.toUpperCase( ) ) {
                    // CASE: Unambiguous variant.
                    variantContent = CONTENT_ENCODING[ variantContent ];
                    if ( variantContent > 28 ) {
                        VARC_VARIANTS_CONTENT[ "AA_CLUSTER_OTHER" ].push( [ c, labelIndex, variantContent ] );
                    } else if ( variantContent > 18 ) {
                        VARC_VARIANTS_CONTENT[ "AA_CLUSTER_UNPOLAR" ].push( [ c, labelIndex, variantContent ] );
                    } else if ( variantContent > 8 ) {
                        VARC_VARIANTS_CONTENT[ "AA_CLUSTER_POLAR" ].push( [ c, labelIndex, variantContent ] );
                    }
                    hasUnambiguousVariant = true;
                    variantPositionsIndicatorData.push( [ c, 0, 1 ] );
                } else {
                    // CASE: Ambiguous variant.
                    if ( VARC_AMBIGUITY_MODE === 'SHOW' || VARC_AMBIGUITY_MODE === 'MASK' ) {
                        variantContent = CONTENT_ENCODING[ variantContent.toUpperCase( ) ];
                        if ( variantContent > 28 ) {
                            VARC_VARIANTS_CONTENT_AMBIGUOUS[ "AA_CLUSTER_OTHER" ].push( [ c, labelIndex, variantContent ] );
                        } else if ( variantContent > 18 ) {
                            VARC_VARIANTS_CONTENT_AMBIGUOUS[ "AA_CLUSTER_UNPOLAR" ].push( [ c, labelIndex, variantContent ] );
                        } else if ( variantContent > 8 ) {
                            VARC_VARIANTS_CONTENT_AMBIGUOUS[ "AA_CLUSTER_POLAR" ].push( [ c, labelIndex, variantContent ] );
                        }
                    }
                    if ( VARC_AMBIGUITY_MODE === 'SHOW' ) {
                        variantPositionsIndicatorData.push( [ c, 0, 1 ] );
                    }
                }
            }
        }
        if ( VARC_AMBIGUITY_MODE !== 'HIDE' || refTransNucContent !== "None" || hasUnambiguousVariant ) {
            superpositions.push( superposition );
            c += 1;
        }
    }

    // Edit variant entry labels to display sample counts.
    if ( VARC_AGGREGATION_MODE !== 'SAMPLE' ) {
        let noSamples;
        switch ( VARC_AGGREGATION_MODE ) {
            case 'GENOTYPE':
                variantYAxisName = "Genotypes";
                let perGenotypeNoSamples = new Object( );
                // Compute number of samples per proteoform.
                for ( let genotype of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Genotypes" ] ) ) {
                    perGenotypeNoSamples[ genotype ] = 0;
                    for ( let sample of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Genotypes" ][ genotype ][ "Samples" ] ) ) {
                        perGenotypeNoSamples[ genotype ] += 1;
                    }
                }
                for ( let label of VARC_VARIANTS_LABELS ) {
                    noSamples = perGenotypeNoSamples[ label ];
                    variantsEditedLabels.push( ( VARC_VARIANTS_LABELS.indexOf( label ) + 1 ).toString( ) + " (" + noSamples + ( noSamples == 1 ? " Sample)" : " Samples)" ) );
                }
                break;
            case 'PROTEOFORM':
                variantYAxisName = "Proteoforms";
                let perProteoformNoSamples = new Object( );
                // Compute number of samples per proteoform.
                for ( let proteoform of Object.keys( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ] ) ) {
                    perProteoformNoSamples[ proteoform ] = 0;
                    for ( let genotype of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ][ proteoform ][ "Genotypes" ] ) ) {
                        for ( let sample of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Genotypes" ][ genotype ][ "Samples" ] ) ) {
                            perProteoformNoSamples[ proteoform ] += 1;
                        }
                    }
                }
                for ( let label of VARC_VARIANTS_LABELS ) {
                    noSamples = perProteoformNoSamples[ label ];
                    variantsEditedLabels.push( ( VARC_VARIANTS_LABELS.indexOf( label ) + 1 ).toString( ) + " (" + noSamples + ( noSamples == 1 ? " Sample)" : " Samples)" ) );
                }
                break;
        }
    } else {
        variantsEditedLabels = VARC_VARIANTS_LABELS; 
    }
    referenceNumericalDataLabels = [ "Hydropathicity", "Variability" ];
    referenceCategoricalDataLabels = [ "Secondary Structure", DATA[ "ReferenceProteinName" ] + " (Seq.)", DATA[ "ReferenceFeatureName" ] + " (Transl. Seq.)" ];
    // Add buttons to select chains
    for ( let chainIdentifier of CHAIN_IDENTIFIERS ) {
        let existingBtn = document.getElementById( "menuVARCChainSelect" + chainIdentifier + "Btn" );
        if ( existingBtn == null ) {
            let btn = document.createElement( "button" );
            btn.innerHTML = "Chain " + chainIdentifier;
            btn.id = "menuVARCChainSelect" + chainIdentifier + "Btn";
            if ( chainIdentifier == SELECTED_CHAIN ) {
                btn.classList.add( "menuBtnActive" );
            }
            btn.onclick = () => { VARCSelectChain( chainIdentifier ); };
            document.getElementById( "menuVARCChainSelectContent" ).appendChild( btn );
        }
    }
    VARC_ECHART_OPTION = {
        animation: false,
        title: [
            {
                text: 'Reference Information',
                textStyle: {
                    fontSize: 11
                },
                left: 'center',
                top: '0'
            },
            {
                text: 'Variants Information',
                textStyle: {
                    fontSize: 11
                },
                left: 'center',
                top: (7 + ( referenceNumericalsCount * 4 ) + ( referenceCategoricalsCount * 4 ) ).toString( ) + '%',
            },
            {
                text: 'Color Scheme',
                textStyle: {
                    fontSize: 11
                },
                right: '70px',
                top: '0',
            },
            {
                text: 'Variability',
                textStyle: {
                    fontSize: 9
                },
                right: '70px',
                top: '255px',
            },
            {
                text: 'Hydropathicity',
                textStyle: {
                    fontSize: 9
                },
                right: '60px',
                top: '305px',
            }
        ],
        grid: [
            {   // Grid to display numerical reference information.
                top: '5%',
                bottom: (95 - ( referenceNumericalsCount * 4 ) ).toString( ) + '%',
                left: '9%',
                right: '230px',
                show: true
            },
            {   // Grid to display categorical reference information.
                top: ( 5 + ( referenceNumericalsCount * 4 ) ).toString( ) + '%',
                bottom: (95 - ( referenceNumericalsCount * 4 ) - ( referenceCategoricalsCount * 4 ) ).toString( ) + '%',
                left: '9%',
                right: '230px',
                show: true
            },
            {   // Grid to display variants information.
                top: (12 + ( referenceNumericalsCount * 4 ) + ( referenceCategoricalsCount * 4 ) ).toString( ) + '%',
                bottom: '19%',
                left: '9%',
                right: '230px',
                show: true
            },
            {   // Grid to display position zoom and summary information.
                top: '93.5%',
                bottom: '20px',
                left: '9%',
                right: '230px',
                show: true
            }
        ],
        xAxis: [
            {   // X axis to display numerical reference information.
                type: 'category',
                data: superpositions,
                splitLine: {
                    show: true,
                    interval: 0
                },
                axisTick: {
                    show: false
                },
                axisLabel: {
                    show: false
                },
                axisLine: {
                    show: false
                },
                show: true,
                gridIndex: 0
            },
            {   // X axis to display categorical reference information.
                type: 'category',
                data: superpositions,
                splitLine: {
                    show: true,
                    interval: 0
                },
                axisTick: {
                    show: false
                },
                axisLabel: {
                    show: false
                },
                axisLine: {
                    show: false
                },
                show: true,
                gridIndex: 1
            },
            {   // Primary X axis to display variants information.
                type: 'category',
                data: superpositions,
                name: 'Superimposed Protein Position',
                nameLocation: 'middle',
                nameTextStyle: {
                    fontSize: 11,
                    fontWeight: 'bold'
                },
                axisLabel: {
                    fontSize: 10
                },
                nameGap: 25,
                offset: 0,
                splitLine: {
                    show: true,
                    interval: 0
                },
                show: true,
                gridIndex: 2
            },
            {   // Secondary X axis to display position zoom and summary information (not affected by data zoom).
                type: 'category',
                data: superpositions,
                axisLabel: {
                    fontSize: 10
                },
                position: 'bottom',
                show: true,
                gridIndex: 3
            }
        ],
        yAxis: [
            {   // Y axis to display numerical reference information.
                type: 'category',
                data: referenceNumericalDataLabels,
                splitArea: {
                    show: false
                },
                axisTick: {
                    alignWithLabel: true
                },
                axisLabel: {
                    fontSize: 10,
                    align: 'right'
                },
                gridIndex: 0
            },
            {   // Y axis to display categorical reference information.
                type: 'category',
                data: referenceCategoricalDataLabels,
                splitArea: {
                    show: false
                },
                axisTick: {
                    alignWithLabel: true
                },
                axisLabel: {
                    fontSize: 10,
                    align: 'right'
                },
                gridIndex: 1,
                inverse: true
            },
            {   // Y axis to display variants information.
                type: 'category',
                data: variantsEditedLabels,
                name: variantYAxisName,
                nameLocation: 'middle',
                nameTextStyle: {
                    fontSize: 11,
                    fontWeight: 'bold'
                },
                nameGap: 110,
                splitArea: {
                    show: false
                },
                axisTick: {
                    alignWithLabel: true
                },
                axisLabel: {
                    fontSize: 10,
                    align: 'right'
                },
                gridIndex: 2,
                inverse: true
            },
            {   // Hidden Y axis to display position zoom and summary information.
                type: 'category',
                data: [ "Variant" ],
                show: false,
                gridIndex: 3,
                inverse: true
            }
        ],
        legend: { // Disable default legend.
            show: false
        },
        dataZoom: [
            {   // Data zoom to zoom into position intervals.
                id: "positionZoom",
                type: 'slider',
                xAxisIndex: [ 0, 1, 2 ],
                realtime: false,
                throttle: 100,
                height: "8px",
                bottom: "21px",
                id: "positionZoom",
                showDataShadow: 'auto',
                backgroundColor: "transparent",
                fillerColor: "#60719645",
                borderColor: "#E4E5ED",
                handleStyle: {
                    borderColor: "#607196",
                    borderWidth: 3
                },
                moveHandleSize: 0,
                brushStyle: {
                    color: '#60719685',
                    borderCap: 'round',
                    opacity: 1
                }
            },
            {   // Data zoom to zoom into single samples.
                type: 'slider',
                yAxisIndex: [ 2 ],
                realtime: false,
                throttle: 5,
                width: "6px",
                left: '0',
                labelFormatter: "",
                backgroundColor: "transparent",
                fillerColor: "#60719645",
                borderColor: "#E4E5ED",
                maxValueSpan: 50,
                handleStyle: {
                    borderColor: "#607196",
                    borderWidth: 3
                },
                moveHandleSize: 0,
                brushStyle: {
                    color: '#60719685',
                    borderCap: 'round',
                    opacity: 1
                }
            }
        ],
        visualMap: [
            // Visual maps to color categorical heatmap components.
            { // Cluster one, polar amino acids.
                type: 'piecewise',
                pieces: [
                    // Amino Acids, Polar (positive), Basic
                    { min: 9, max: 9, color: "#27269A" },
                    { min: 10, max: 10, color: "#2F49B1" },
                    { min: 11, max: 11, color: "#3874C8" },
                    // Polar (negative), Acidic
                    { min: 12, max: 12, color: "#F75050" },
                    { min: 13, max: 13, color: "#FB6581" },
                    // Amino Acids, Polar (neutral)
                    { min: 14, max: 14, color: "#1EAB13" },
                    { min: 15, max: 15, color: "#38B42B" },
                    { min: 16, max: 16, color: "#51BC43" },
                    { min: 17, max: 17, color: "#6AC55B" },
                    // Amino Acids, Sulfur bridge forming
                    { min: 18, max: 18, color: "#DEC443" }
                ],
                seriesIndex: [3,6,9,12],
                show: true,
                orient: 'vertical',
                align: 'right',
                top: '25px',
                right: '170px',
                itemHeight: 8,
                itemWidth: 10,
                itemSymbol: 'rect',
                textStyle: {
                    fontSize: 10
                },
                formatter: function (value) {
                    return CONTENT_DECODING[ value ];
                }
            },
            { // Cluster two, unpolar amino acids.
                type: 'piecewise',
                pieces: [
                    // Amino Acids, Aromatic
                    { min: 19, max: 19, color: "#5826ED" },
                    { min: 20, max: 20, color: "#653CEF" },
                    { min: 21, max: 21, color: "#7352F1" },
                    // Amino Acids, Aliphatic
                    { min: 22, max: 22, color: "#239E90" },
                    { min: 23, max: 23, color: "#34A597" },
                    { min: 24, max: 24, color: "#45AD9E" },
                    { min: 25, max: 25, color: "#55B4A5" },
                    { min: 26, max: 26, color: "#66BCAC" },
                    { min: 27, max: 27, color: "#77C3B4" },
                    { min: 28, max: 28, color: "#88CABC" }
                ],
                seriesIndex: [4,7,10,13],
                show: true,
                orient: 'vertical',
                align: 'right',
                top: '25px',
                right: '120px',
                itemHeight: 8,
                itemWidth: 10,
                itemSymbol: 'rect',
                textStyle: {
                    fontSize: 10
                },
                formatter: function (value) {
                    return CONTENT_DECODING[ value ];
                }
            },
            { // Cluster three, other amino acids.
                type: 'piecewise',
                pieces: [
                    // Amino Acids, Other
                    { min: 29, max: 29, color: "#E61090" },
                    { min: 30, max: 30, color: "#75551D" },
                    { min: 31, max: 31, color: "#75551D" },
                    { min: 32, max: 32, color: "#212121" }
                ],
                seriesIndex: [5,8,11,14],
                show: true,
                orient: 'vertical',
                align: 'right',
                top: '25px',
                right: '20px',
                itemHeight: 8,
                itemWidth: 10,
                itemSymbol: 'rect',
                textStyle: {
                    fontSize: 10
                },
                formatter: function (value) {
                    let content = CONTENT_DECODING[ value ];
                    switch (content) {
                        case 'DEL':
                            return 'Deletion';
                        case 'INC':
                            return 'Incomplete';
                        case 'ANY':
                            return 'Any';
                        case 'TER':
                            return 'Termination';
                    }
                }
            },
            { // Secondary structures.
                type: 'piecewise',
                pieces: [
                    { min: 6, max: 6, color: COLORS[ CONTENT_DECODING[ 6 ] ] },
                    { min: 7, max: 7, color: COLORS[ CONTENT_DECODING[ 7 ] ] },
                    { min: 8, max: 8, color: COLORS[ CONTENT_DECODING[ 8 ] ] }
                ],
                seriesIndex: [2],
                show: true,
                orient: 'vertical',
                align: 'right',
                top: '125px',
                right: '20px',
                itemHeight: 8,
                itemWidth: 10,
                itemSymbol: 'rect',
                textStyle: {
                    fontSize: 10
                },
                formatter: function (value) {
                    let content = CONTENT_DECODING[ value ];
                    switch (content) {
                        case 'sheet':
                            return 'Sheet';
                        case 'helix':
                            return 'Helix';
                        case 'coil':
                            return 'Coil';
                    }
                }
            },
            { // Variant position indicator.
                type: 'piecewise',
                pieces: [
                    { min: 1, max: 1, color: '#607196' }
                ],
                seriesIndex: [ 15 ],
                show: false
            },
            {   // Visual map for numerical uncertainty series.
                min: 0 - Math.pow( 10, -6 ),
                max: 1,
                range: [ 0, 1 ],
                calculable: false,
                inRange: {
                    color: [
                      '#c8c8c8',
                      '#e4c880',
                      '#ffc837',
                      '#ff9337',
                      '#ff4a37',
                      '#ff0037'
                    ]
                },
                outOfRange: {
                    color: [
                        '#EFF0F8'
                    ]
                },
                seriesIndex: [ 0 ],
                show: true,
                orient: 'horizontal',
                align: 'right',
                top: '270px',
                right: '25px',
                precision: 2,
                textStyle: {
                    fontSize: 10
                }
            },
            {   // Visual map for numerical hydropathicity series.
                min: HYDROPATHICITY_MIN_VAL - Math.pow( 10, -6 ),
                max: HYDROPATHICITY_MAX_VAL,
                range: [ HYDROPATHICITY_MIN_VAL, HYDROPATHICITY_MAX_VAL ],
                calculable: false,
                inRange: {
                    color: [
                      '#45a7e6',
                      '#e6e6e6',
                      '#b87737'
                    ]
                },
                outOfRange: {
                    color: [
                        '#EFF0F8'
                    ]
                },
                seriesIndex: [ 1 ],
                show: true,
                orient: 'horizontal',
                align: 'right',
                top: '320px',
                right: '25px',
                precision: 2,
                textStyle: {
                    fontSize: 10
                }
            }
        ],
        series: [
            {
                name: 'RefVariability',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_PROTEIN_VARIABILITY,
                xAxisIndex: 0,
                yAxisIndex: 0,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'RefHydropathicity',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_PROTEIN_HYDROP,
                xAxisIndex: 0,
                yAxisIndex: 0,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'RefSecStruc',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_PROTEIN_SS,
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'RefProteinCont_Polar',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_PROTEIN_CONTENT[ "AA_CLUSTER_POLAR" ],
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'RefProteinCont_Unpolar',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_PROTEIN_CONTENT[ "AA_CLUSTER_UNPOLAR" ],
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'RefProteinCont_Other',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_PROTEIN_CONTENT[ "AA_CLUSTER_OTHER" ],
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'RefTransFeatureCont_Polar',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_TRANSLATED_FEATURE_CONTENT[ "AA_CLUSTER_POLAR" ],
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'RefTransFeatureCont_Unpolar',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_TRANSLATED_FEATURE_CONTENT[ "AA_CLUSTER_UNPOLAR" ],
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'RefTransFeatureCont_Other',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_REF_TRANSLATED_FEATURE_CONTENT[ "AA_CLUSTER_OTHER" ],
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'VariantsUnamb_Polar',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_VARIANTS_CONTENT[ "AA_CLUSTER_POLAR" ],
                xAxisIndex: 2,
                yAxisIndex: 2,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'VariantsUnamb_Unpolar',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_VARIANTS_CONTENT[ "AA_CLUSTER_UNPOLAR" ],
                xAxisIndex: 2,
                yAxisIndex: 2,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'VariantsUnamb_Other',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: VARC_VARIANTS_CONTENT[ "AA_CLUSTER_OTHER" ],
                xAxisIndex: 2,
                yAxisIndex: 2,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'VariantsAmb_Polar',
                type: 'heatmap',
                data: VARC_AMBIGUITY_MODE !== 'HIDE' ? VARC_VARIANTS_CONTENT_AMBIGUOUS[ "AA_CLUSTER_POLAR" ] : [ ],
                xAxisIndex: 2,
                yAxisIndex: 2,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1,
                    opacity: VARC_AMBIGUITY_MODE !== 'MASK' ? 1.0 : 0.2
                }
            },
            {
                name: 'VariantsAmb_Unpolar',
                type: 'heatmap',
                data: VARC_AMBIGUITY_MODE !== 'HIDE' ? VARC_VARIANTS_CONTENT_AMBIGUOUS[ "AA_CLUSTER_UNPOLAR" ] : [ ],
                xAxisIndex: 2,
                yAxisIndex: 2,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1,
                    opacity: VARC_AMBIGUITY_MODE !== 'MASK' ? 1.0 : 0.2
                }
            },
            {
                name: 'VariantsAmb_Other',
                type: 'heatmap',
                data: VARC_AMBIGUITY_MODE !== 'HIDE' ? VARC_VARIANTS_CONTENT_AMBIGUOUS[ "AA_CLUSTER_OTHER" ] : [ ],
                xAxisIndex: 2,
                yAxisIndex: 2,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1,
                    opacity: VARC_AMBIGUITY_MODE !== 'MASK' ? 1.0 : 0.2
                }
            },
            {
                name: 'PositionIndicator',
                type: 'heatmap',
                progressive: 500,
                progressiveThreshold: 1000,
                data: variantPositionsIndicatorData,
                xAxisIndex: 3,
                yAxisIndex: 3,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.1
                }
            },
            {
                name: 'PositionLine',
                type: 'line',
                data: [ ],
                xAxisIndex: 2,
                yAxisIndex: 2
            }
        ]
    };
    VARC_ECHART.on('datazoom', function(params) {
        if ( params.dataZoomId == "positionZoom" ) {
            let maxMatchedPosition = superpositions.at( -1 );
            let fromIndex = Math.floor( ( params.start / 100 ) * maxMatchedPosition );
            if ( fromIndex === 0 ) {
                fromIndex = 1;
            }
            let from = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ fromIndex ][ "ProteinPosition" ];
            let toIndex = Math.floor( ( params.end / 100 ) * maxMatchedPosition );
            let to = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ toIndex ][ "ProteinPosition" ];
            while( from === -1 && fromIndex <= toIndex ) {
                fromIndex += 1;
                from = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ fromIndex ][ "ProteinPosition" ];
            }
            while( to === -1 && toIndex >= fromIndex ) {
                toIndex -= 1;
                to = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ toIndex ][ "ProteinPosition" ];
            }
            STRC_HIGHLIGHTED_REGION = [ ];
            STRC_HIGHLIGHTED_REGION.push( from );
            STRC_HIGHLIGHTED_REGION.push( to );
            STRCApplyStyle( );
        }
    });
    VARC_ECHART.on( 'click', function( params ) {
        VARCMarkPosition( params.data[ 0 ] );
        STRCClearSelection( );
        if ( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ params.data[ 0 ] + 1 ][ "ProteinContent" ] !== "None" ) {
            let resi = parseInt( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ params.data[ 0 ] + 1 ][ "ProteinPosition" ] );
            let resn = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ params.data[ 0 ] + 1 ][ "ProteinContent"];
            STRCSelectResidue( SELECTED_CHAIN, resi, resn, true, COLORS[ resn ], COLORS[ resn ], 0.6, 1 );
        }
        COMCDisplay( params.data[ 0 ] + 1 );
        INFCSetReferenceInformation( params.data[ 0 ] + 1  );
        INFCSetVariantInformation(
            params.data[ 0 ] + 1,
            params.seriesName.startsWith( "Variants" ) ? variantsEditedLabels[ params.data[ 1 ] ] : false,
            params.seriesName.startsWith( "Variants" ) ? VARC_VARIANTS_LABELS[ params.data[ 1 ] ] : false,
            VARC_PROTEOFORM_ACCESSOR_MAP[ VARC_VARIANTS_LABELS[ params.data[ 1 ] ] ],
            VARC_GENOTYPE_ACCESSOR_MAP[ VARC_VARIANTS_LABELS[ params.data[ 1 ] ] ]
        );
        STRCDisplayContacts( );
        STRCApplyStyle( );
    } );
    VARC_ECHART.setOption(VARC_ECHART_OPTION);
    STRC_HIGHLIGHTED_REGION.push( 1 );
    STRC_HIGHLIGHTED_REGION.push( superpositions.at( -1 ) );
    STRCApplyStyle( );
}

/**
 *
 */
function VARCSetAmbiguityMode( mode ) {
    if ( VARC_AMBIGUITY_MODE !== mode ) {
        VARC_AMBIGUITY_MODE = mode;
        document.getElementById( "menuVARCAmbiguityModeSelectMaskBtn" ).classList.remove( "menuBtnActive" );
        document.getElementById( "menuVARCAmbiguityModeSelectShowBtn" ).classList.remove( "menuBtnActive" );
        document.getElementById( "menuVARCAmbiguityModeSelectHideBtn" ).classList.remove( "menuBtnActive" );
        switch (mode) {
            case "MASK":
                document.getElementById( "menuVARCAmbiguityModeSelectMaskBtn" ).classList.add( "menuBtnActive" );
                break;
            case "SHOW":
                document.getElementById( "menuVARCAmbiguityModeSelectShowBtn" ).classList.add( "menuBtnActive" );
                break;
            case "HIDE":
                document.getElementById( "menuVARCAmbiguityModeSelectHideBtn" ).classList.add( "menuBtnActive" );
                break;
        }
        STRCcomputeVariability( );
        initVARC( );
    }
}

/**
 *
 */
 function VARCSetAggregationMode( mode ) {
    if ( VARC_AGGREGATION_MODE !== mode ) {
        VARC_AGGREGATION_MODE = mode;
        document.getElementById( "menuVARCAggregationModeSelectSampleBtn" ).classList.remove( "menuBtnActive" );
        document.getElementById( "menuVARCAggregationModeSelectGenotypeBtn" ).classList.remove( "menuBtnActive" );
        document.getElementById( "menuVARCAggregationModeSelectProteoformBtn" ).classList.remove( "menuBtnActive" );
        switch (mode) {
            case "SAMPLE":
                document.getElementById( "menuVARCAggregationModeSelectSampleBtn" ).classList.add( "menuBtnActive" );
                break;
            case "GENOTYPE":
                document.getElementById( "menuVARCAggregationModeSelectGenotypeBtn" ).classList.add( "menuBtnActive" );
                break;
            case "PROTEOFORM":
                document.getElementById( "menuVARCAggregationModeSelectProteoformBtn" ).classList.add( "menuBtnActive" );
                break;
        }
        initVARC( );
    }
}

/**
 * Selects the specified chain as active chain.
 * 
 * @param {String} chainIdentifier 
 */
function VARCSelectChain( chainIdentifier ) {
    if ( chainIdentifier != SELECTED_CHAIN && IS_LOADED ) {
        SELECTED_CHAIN = chainIdentifier;
        initVARC( );
        VARC_ECHART.dispatchAction(
            {
                type: 'dataZoom',
                dataZoomIndex: 0,
                start: 0,
                end: 100
            }
        );
    }
    for ( let cId of CHAIN_IDENTIFIERS ) {
        let btn = document.getElementById( "menuVARCChainSelect" + cId + "Btn" );
        if ( cId == SELECTED_CHAIN ) {
            btn.classList.add( "menuBtnActive" );
        } else {
            btn.classList.remove( "menuBtnActive" );
        }
    }
}

/**
 * Marks a position on the variants view with a dashed line.
 * 
 * @param {Integer} position The position to mark.
 */
function VARCMarkPosition( position ) {
    let d = [ ]
    for ( let label of VARC_VARIANTS_LABELS ) {
        d.push( [ position, VARC_VARIANTS_LABELS.indexOf( label ) ] );
    }
    VARC_ECHART.setOption({
        series: [
            {
                name: 'PositionLine',
                type: 'line',
                data: d,
                xAxisIndex: 2,  
                yAxisIndex: 2,
                symbol: 'none',
                lineStyle: {
                    type: 'dashed',
                    color: COLORS.DEL,
                    width: 0.8,
                    opacity: 0.8
                }
            }
        ]
    });
}

/**
 * TODO
 * @param {*} position 
 */
function COMCDisplay( position ) {
    closeElement( "compositionViewContentDummy" );
    COMC_DISPLAYED_POSITION = position;
    var sunburstData = new Object();
    var referenceAminoAcidContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "ReferenceFeatureTranslatedContent" ].toUpperCase( );
    var referenceNucleotideContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "ReferenceFeatureContent" ];
    for ( var [ proteoform, aminoAcidContent ] of Object.entries( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "PerProteoformAminoAcidContent" ] ) ) {
        if ( VARC_AMBIGUITY_MODE !== "SHOW" && ( aminoAcidContent === aminoAcidContent.toLowerCase( ) ) ) {
            continue;
        } else {
            aminoAcidContent = aminoAcidContent.toUpperCase( );
        }
        if ( !( aminoAcidContent in sunburstData ) ) {
            sunburstData[ aminoAcidContent ] = {
                "Proteoforms": [ ],
                "NucleotideContents": { }
            };
        };
        sunburstData[ aminoAcidContent ][ "Proteoforms" ].push( proteoform );
        for ( let genotype of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Proteoforms" ][ proteoform ][ "Genotypes" ] ) ) {
            let nucleotideContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "PerGenotypeNucleotideContent" ][ genotype ];
            if ( !( nucleotideContent in sunburstData[ aminoAcidContent ][ "NucleotideContents" ] ) ) {
                sunburstData[ aminoAcidContent ][ "NucleotideContents" ][ nucleotideContent ] = {
                    "Genotypes": [ ],
                    "Samples": [ ]
                };
            };
            sunburstData[ aminoAcidContent ][ "NucleotideContents" ][ nucleotideContent ][ "Genotypes" ].push( genotype );
            for ( let sample of Object.values( DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ "Genotypes" ][ genotype ][ "Samples" ] ) ) {
                sunburstData[ aminoAcidContent ][ "NucleotideContents" ][ nucleotideContent ][ "Samples" ].push( sample );
            }
        }
    }
    var sunburstComponents = [ ];
    for ( var [ aminoAcidContent, obj ] of Object.entries( sunburstData ) ) {
        let nucleotideContentObjects = [ ];
        for ( var [ nucleotideContent, obj ] of Object.entries( sunburstData[ aminoAcidContent ][ "NucleotideContents" ] ) ) {
            let sampleContentObjects = [ ];
            for ( var sample of sunburstData[ aminoAcidContent ][ "NucleotideContents" ][ nucleotideContent ][ "Samples" ] ) {
                sampleContentObjects.push( {
                    name: sample,
                    value: 1,
                    itemStyle: {
                        color: '#cbd0e0'
                    }
                } );
            }
            nucleotideContentObjects.push( {
                name: ( nucleotideContent === referenceNucleotideContent ) ? nucleotideContent + "\nReference" : nucleotideContent,
                itemStyle: {
                    color: '#607196'
                },
                children: sampleContentObjects
            } );
        }
        sunburstComponents.push( {
            name: ( aminoAcidContent === referenceAminoAcidContent ) ? aminoAcidContent + "\nReference" : aminoAcidContent,
            itemStyle: {
                color: COLORS[ aminoAcidContent ]
            },
            children: nucleotideContentObjects
        } );
    }
    // Generate and set ECharts option object.
    COMC_ECHART_OPTION = {
        title: [
            {
                text: 'Amino Acid and Nucleotide Composition\nPosition ' + position,
                textStyle: {
                    fontSize: 12
                },
                left: 'center',
                top: '4px'
            }
        ],
        series: [
            {
                type: 'sunburst',
                center: [
                    "50%",
                    "55%"
                ],
                radius: [ 0, '20%' ],
                data: sunburstComponents,
                emphasis: {
                    focus: 'ancestor'
                },
                levels: [
                    {
                        r0: '0%',
                        r: '20%'
                    },
                    {
                        r0: '20%',
                        r: '40%',
                        itemStyle: {
                            borderWidth: 2
                        },
                        label: {
                            rotate: 'tangential',
                            fontSize: 10
                        }
                    },
                    {
                        r0: '40%',
                        r: '60%',
                        label: {
                            rotate: 'tangential',
                            fontSize: 9
                        }
                    },
                    {
                        r0: '60%',
                        r: '63%',
                        label: {
                            position: 'outside',
                            padding: 2,
                            silent: false,
                            fontSize: 8
                        },
                        itemStyle: {
                            borderWidth: 3
                        }
                    }
                ],
                itemStyle: {
                    color: '#607196'
                }
            }
        ]
    };
    COMC_ECHART.setOption( COMC_ECHART_OPTION );
}

/**
 * 
 * @param {*} htmlString 
 */
function INFCSetContactInformation( htmlString ) {
    closeElement( "infoContentDummy" );
    showElement("infoContent", "block");
    document.getElementById( "infoContactContent" ).innerHTML = htmlString;
}

/**
 * 
 * @param {*} htmlString 
 */
function INFCSetReferenceInformation( position ) {
    closeElement( "infoContentDummy" );
    showElement("infoContent", "block");
    let htmlString = "";
    let genomePositions = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "ReferenceFeaturePositions" ];
    let genomeContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "ReferenceFeatureContent" ];
    let translatedGenomeContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "ReferenceFeatureTranslatedContent" ];
    let translatedGenomeContentColor;
    let proteinPosition = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "ProteinPosition" ];
    let proteinContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "ProteinContent" ];
    let proteinContentColor;
    let proteinStruc;
    let proteinStrucColor;
    if ( genomeContent === "None" ) {
        htmlString += "<br><b style='color: #cbd0e0'>No superposition with genome.</b><br>" 
    } else {
        htmlString += "<br>&bull; Genome Positions: " + genomePositions.join( ", " ) + "<br>";
        htmlString += "<br>&bull; Genome Content: " + genomeContent + "<br>";
        if ( translatedGenomeContent === "ANY" ) {
            translatedGenomeContent = "Any (codon contains N)<br>"
            translatedGenomeContentColor = COLORS.ANY;
        } else if ( translatedGenomeContent === "TER" ) {
            translatedGenomeContent = "Termination<br>"
            translatedGenomeContentColor = COLORS.TER;
        } else {
            translatedGenomeContentColor = COLORS[ translatedGenomeContent ];
        }
        htmlString += "<br>&bull; Transl. Genome Content " + "<span style='color: " + translatedGenomeContentColor + "'>&#9724; </span> " + translatedGenomeContent + "<br>";
    }
    if ( proteinPosition === -1 ) {
        htmlString += "<br><b style='color: #cbd0e0'>No superposition with protein structure.</b><br>" 
    } else {
        htmlString += "<br>&bull; Protein Position: " + proteinPosition + "<br>";
        proteinContentColor = COLORS[ proteinContent ];
        htmlString += "<br>&bull; Protein Content: " + "<span style='color: " + proteinContentColor + "'>&#9724; </span> " + proteinContent + "<br>";
        proteinStruc = PROTEIN_RESIDUES[ SELECTED_CHAIN ][ proteinPosition ][ "ss" ];
        proteinStrucColor = COLORS[ proteinStruc ];
        htmlString += "<br>&bull; Protein 2ry Structure: " + "<span style='color: " + proteinStrucColor + "'>&#9724; </span> " + proteinStruc.charAt(0).toUpperCase( ) + proteinStruc.slice(1) + "<br>";
    }
    document.getElementById( "infoReferenceContent" ).innerHTML = htmlString;
}

/**
 * TODO
 * @param {*} htmlString 
 */
 function INFCSetVariantInformation( position, label, variantId, proteoformId, genotypeId ) {
    closeElement( "infoContentDummy" );
    showElement("infoContent", "block");
    let htmlString = "";
    INFC_SUBSTITUTION = "None";
    if ( variantId === false ) {
        htmlString = "<br><b style='color: #cbd0e0'>No variants data selected.</b><br>" ;
        document.getElementById( "infoVariantContent" ).innerHTML = htmlString;
    } else {
        let genomePositions;
        let genomeContent;
        let translatedGenomeContent;
        let translatedGenomeContentColor;
        switch ( VARC_AGGREGATION_MODE ) {
            case 'SAMPLE':
                document.getElementById( "infoVariantHeader" ).innerHTML = "<hr><b>Sample Information</b>";
                htmlString += "<br>&bull; Name: " + label + "<br>";
                genomePositions = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "PerGenotypePositions" ][ genotypeId ];
                if ( JSON.stringify( genomePositions ) == JSON.stringify( [ -1, -1, -1 ] ) ) {
                    genomePositions = [ "None" ];
                }
                genomeContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "PerGenotypeNucleotideContent" ][ genotypeId ];
                htmlString += "<br>&bull; Genome Positions: " + genomePositions.join( ", " ) + "<br>";
                htmlString += "<br>&bull; Genome Content: " + genomeContent + "<br>";
                break;
            case 'GENOTYPE':
                document.getElementById( "infoVariantHeader" ).innerHTML = "<hr><b>Genotype Information</b>";
                htmlString += "<br>&bull; Name: Genotype " + label + "<br>";
                genomePositions = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "PerGenotypePositions" ][ genotypeId ];
                if ( JSON.stringify( genomePositions ) == JSON.stringify( [ -1, -1, -1 ] ) ) {
                    genomePositions = [ "None" ];
                }
                genomeContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "PerGenotypeNucleotideContent" ][ genotypeId ];
                htmlString += "<br>&bull; Genome Positions: " + genomePositions.join( ", " ) + "<br>";
                htmlString += "<br>&bull; Genome Content: " + genomeContent + "<br>";
                break;
            case 'PROTEOFORM':
                document.getElementById( "infoVariantHeader" ).innerHTML = "<hr><b>Proteoform Information</b>";
                htmlString += "<br>&bull; Name: Proteoform " + label + "<br>";
                break;
        }
        translatedGenomeContent = DATA[ "SuperpositionChain" + SELECTED_CHAIN ][ position ][ "PerProteoformAminoAcidContent" ][ proteoformId ];
        if ( translatedGenomeContent === translatedGenomeContent.toLowerCase( ) ) {
            document.getElementById( "infoVariantHeader" ).innerHTML += " <i style='color: #cbd0e0;' class='fas fa-exclamation-circle'></i>";
            translatedGenomeContent = translatedGenomeContent.toUpperCase( );
        }
        translatedGenomeContentColor;
        if ( translatedGenomeContent === "ANY" ) {
            translatedGenomeContent = "Any (codon contains N)<br>"
            translatedGenomeContentColor = COLORS.ANY;
        } else if ( translatedGenomeContent === "TER" ) {
            translatedGenomeContent = "Termination<br>"
            translatedGenomeContentColor = COLORS.TER;
        } else if ( translatedGenomeContent === "INC" ) {
            translatedGenomeContent = "Incomplete (frameshift)<br>"
            translatedGenomeContentColor = COLORS.INC;
        } else if ( translatedGenomeContent === "DEL" ) {
            translatedGenomeContent = "Deletion<br>"
            translatedGenomeContentColor = COLORS.DEL;
        } else {
            translatedGenomeContentColor = COLORS[ translatedGenomeContent ];
            INFC_SUBSTITUTION = translatedGenomeContent;
        }
        htmlString += "<br>&bull; Transl. Genome Content: " + "<span style='color: " + translatedGenomeContentColor + "'>&#9724; </span> " + translatedGenomeContent + "<br>";
        //
        document.getElementById( "infoVariantContent" ).innerHTML = htmlString;
    }
}