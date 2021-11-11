var DATA = new Object( );
var VIEW_HEIGHT = 0;
var ECHART_AV;
var ECHART_AV_OPTION;
var ECHART_PS;
var ECHART_PS_OPTION;
var PRI15Scores = {
    "CYS-helix/CYS-helix": 5.192, "CYS-helix/CYS-coil": 3.787, "CYS-helix/CYS-sheet": -0.14, "CYS-helix/ASP-helix": 0.349, "CYS-helix/ASP-coil": -Infinity, "CYS-helix/ASP-sheet": -Infinity, "CYS-helix/SER-helix": 0.308, "CYS-helix/SER-coil": -2.024, "CYS-helix/SER-sheet": -3.055, "CYS-helix/GLN-helix": -0.857, "CYS-helix/GLN-coil": 0.0, "CYS-helix/GLN-sheet": -Infinity, "CYS-helix/LYS-helix": -2.304, "CYS-helix/LYS-coil": -4.6, "CYS-helix/LYS-sheet": -1.806, "CYS-helix/ILE-helix": 0.844, "CYS-helix/ILE-coil": -2.97, "CYS-helix/ILE-sheet": -2.457, "CYS-helix/PRO-helix": 0.457, "CYS-helix/PRO-coil": -1.587, "CYS-helix/PRO-sheet": -Infinity, "CYS-helix/THR-helix": 0.228, "CYS-helix/THR-coil": -3.875, "CYS-helix/THR-sheet": -Infinity, "CYS-helix/PHE-helix": -0.613, "CYS-helix/PHE-coil": -Infinity, "CYS-helix/PHE-sheet": -3.113, "CYS-helix/ASN-helix": -0.379, "CYS-helix/ASN-coil": -1.747, "CYS-helix/ASN-sheet": -1.581, "CYS-helix/GLY-helix": 1.417, "CYS-helix/GLY-coil": -0.204, "CYS-helix/GLY-sheet": -2.604, "CYS-helix/HIS-helix": 0.193, "CYS-helix/HIS-coil": -1.534, "CYS-helix/HIS-sheet": -Infinity, "CYS-helix/LEU-helix": 0.602, "CYS-helix/LEU-coil": -2.52, "CYS-helix/LEU-sheet": -2.947, "CYS-helix/ARG-helix": -1.378, "CYS-helix/ARG-coil": -0.157, "CYS-helix/ARG-sheet": -Infinity, "CYS-helix/TRP-helix": 0.09, "CYS-helix/TRP-coil": -2.605, "CYS-helix/TRP-sheet": -Infinity, "CYS-helix/ALA-helix": 0.688, "CYS-helix/ALA-coil": -1.994, "CYS-helix/ALA-sheet": -3.119, "CYS-helix/VAL-helix": 0.749, "CYS-helix/VAL-coil": -2.855, "CYS-helix/VAL-sheet": -2.065, "CYS-helix/GLU-helix": -2.346, "CYS-helix/GLU-coil": -2.336, "CYS-helix/GLU-sheet": -2.688, "CYS-helix/TYR-helix": -1.777, "CYS-helix/TYR-coil": -1.238, "CYS-helix/TYR-sheet": -2.043, "CYS-helix/MET-helix": 0.641, "CYS-helix/MET-coil": -Infinity, "CYS-helix/MET-sheet": -1.764, "CYS-coil/CYS-helix": 3.787, "CYS-coil/CYS-coil": 5.338, "CYS-coil/CYS-sheet": 1.218, "CYS-coil/ASP-helix": -1.035, "CYS-coil/ASP-coil": 1.588, "CYS-coil/ASP-sheet": -Infinity, "CYS-coil/SER-helix": -0.661, "CYS-coil/SER-coil": 1.461, "CYS-coil/SER-sheet": -2.843, "CYS-coil/GLN-helix": -1.516, "CYS-coil/GLN-coil": 0.212, "CYS-coil/GLN-sheet": -Infinity, "CYS-coil/LYS-helix": -1.399, "CYS-coil/LYS-coil": -3.694, "CYS-coil/LYS-sheet": -2.693, "CYS-coil/ILE-helix": -0.146, "CYS-coil/ILE-coil": -0.219, "CYS-coil/ILE-sheet": -1.839, "CYS-coil/PRO-helix": -1.634, "CYS-coil/PRO-coil": -0.18, "CYS-coil/PRO-sheet": -0.776, "CYS-coil/THR-helix": 0.127, "CYS-coil/THR-coil": 0.48, "CYS-coil/THR-sheet": -1.766, "CYS-coil/PHE-helix": -1.377, "CYS-coil/PHE-coil": -1.466, "CYS-coil/PHE-sheet": -1.515, "CYS-coil/ASN-helix": -0.051, "CYS-coil/ASN-coil": -1.471, "CYS-coil/ASN-sheet": -Infinity, "CYS-coil/GLY-helix": 0.208, "CYS-coil/GLY-coil": 0.758, "CYS-coil/GLY-sheet": -2.104, "CYS-coil/HIS-helix": -0.782, "CYS-coil/HIS-coil": 0.204, "CYS-coil/HIS-sheet": 0.011, "CYS-coil/LEU-helix": 0.597, "CYS-coil/LEU-coil": -1.733, "CYS-coil/LEU-sheet": -1.349, "CYS-coil/ARG-helix": -2.24, "CYS-coil/ARG-coil": -0.377, "CYS-coil/ARG-sheet": -1.645, "CYS-coil/TRP-helix": -Infinity, "CYS-coil/TRP-coil": -0.601, "CYS-coil/TRP-sheet": -Infinity, "CYS-coil/ALA-helix": -2.458, "CYS-coil/ALA-coil": -1.089, "CYS-coil/ALA-sheet": -0.961, "CYS-coil/VAL-helix": 1.271, "CYS-coil/VAL-coil": -0.301, "CYS-coil/VAL-sheet": -1.748, "CYS-coil/GLU-helix": -4.385, "CYS-coil/GLU-coil": 0.789, "CYS-coil/GLU-sheet": -1.782, "CYS-coil/TYR-helix": -1.816, "CYS-coil/TYR-coil": -2.53, "CYS-coil/TYR-sheet": -2.054, "CYS-coil/MET-helix": -1.484, "CYS-coil/MET-coil": -1.442, "CYS-coil/MET-sheet": -1.552, "CYS-sheet/CYS-helix": -0.14, "CYS-sheet/CYS-coil": 1.218, "CYS-sheet/CYS-sheet": 5.705, "CYS-sheet/ASP-helix": -3.338, "CYS-sheet/ASP-coil": -0.394, "CYS-sheet/ASP-sheet": 2.546, "CYS-sheet/SER-helix": -3.811, "CYS-sheet/SER-coil": -1.876, "CYS-sheet/SER-sheet": 0.741, "CYS-sheet/GLN-helix": -2.615, "CYS-sheet/GLN-coil": -Infinity, "CYS-sheet/GLN-sheet": 3.024, "CYS-sheet/LYS-helix": -Infinity, "CYS-sheet/LYS-coil": -2.19, "CYS-sheet/LYS-sheet": 2.5, "CYS-sheet/ILE-helix": -3.296, "CYS-sheet/ILE-coil": -1.372, "CYS-sheet/ILE-sheet": 1.642, "CYS-sheet/PRO-helix": -1.739, "CYS-sheet/PRO-coil": -1.109, "CYS-sheet/PRO-sheet": -1.064, "CYS-sheet/THR-helix": -2.412, "CYS-sheet/THR-coil": -1.648, "CYS-sheet/THR-sheet": 2.562, "CYS-sheet/PHE-helix": -2.15, "CYS-sheet/PHE-coil": -1.466, "CYS-sheet/PHE-sheet": 0.395, "CYS-sheet/ASN-helix": -Infinity, "CYS-sheet/ASN-coil": -Infinity, "CYS-sheet/ASN-sheet": -0.781, "CYS-sheet/GLY-helix": -3.147, "CYS-sheet/GLY-coil": -1.105, "CYS-sheet/GLY-sheet": 3.072, "CYS-sheet/HIS-helix": -Infinity, "CYS-sheet/HIS-coil": -Infinity, "CYS-sheet/HIS-sheet": 0.129, "CYS-sheet/LEU-helix": -2.126, "CYS-sheet/LEU-coil": -1.209, "CYS-sheet/LEU-sheet": 1.152, "CYS-sheet/ARG-helix": -Infinity, "CYS-sheet/ARG-coil": -Infinity, "CYS-sheet/ARG-sheet": -1.933, "CYS-sheet/TRP-helix": -Infinity, "CYS-sheet/TRP-coil": -Infinity, "CYS-sheet/TRP-sheet": -0.489, "CYS-sheet/ALA-helix": -3.064, "CYS-sheet/ALA-coil": -0.801, "CYS-sheet/ALA-sheet": 2.549, "CYS-sheet/VAL-helix": -3.408, "CYS-sheet/VAL-coil": -2.46, "CYS-sheet/VAL-sheet": 0.632, "CYS-sheet/GLU-helix": -Infinity, "CYS-sheet/GLU-coil": -Infinity, "CYS-sheet/GLU-sheet": 0.127, "CYS-sheet/TYR-helix": -1.277, "CYS-sheet/TYR-coil": -Infinity, "CYS-sheet/TYR-sheet": 1.59, "CYS-sheet/MET-helix": -2.332, "CYS-sheet/MET-coil": 0.503, "CYS-sheet/MET-sheet": 2.786, "ASP-helix/CYS-helix": 0.349, "ASP-helix/CYS-coil": -1.035, "ASP-helix/CYS-sheet": -3.338, "ASP-helix/ASP-helix": 3.779, "ASP-helix/ASP-coil": -0.709, "ASP-helix/ASP-sheet": 1.056, "ASP-helix/SER-helix": 1.126, "ASP-helix/SER-coil": 0.201, "ASP-helix/SER-sheet": -3.209, "ASP-helix/GLN-helix": 0.559, "ASP-helix/GLN-coil": -1.8, "ASP-helix/GLN-sheet": -3.795, "ASP-helix/LYS-helix": 1.165, "ASP-helix/LYS-coil": -2.023, "ASP-helix/LYS-sheet": -2.365, "ASP-helix/ILE-helix": -0.224, "ASP-helix/ILE-coil": -0.568, "ASP-helix/ILE-sheet": -4.402, "ASP-helix/PRO-helix": -0.364, "ASP-helix/PRO-coil": 0.238, "ASP-helix/PRO-sheet": -2.934, "ASP-helix/THR-helix": 0.416, "ASP-helix/THR-coil": -0.2, "ASP-helix/THR-sheet": -4.211, "ASP-helix/PHE-helix": -1.312, "ASP-helix/PHE-coil": 1.591, "ASP-helix/PHE-sheet": -5.059, "ASP-helix/ASN-helix": 1.15, "ASP-helix/ASN-coil": -2.307, "ASP-helix/ASN-sheet": -3.75, "ASP-helix/GLY-helix": -0.831, "ASP-helix/GLY-coil": -2.033, "ASP-helix/GLY-sheet": -3.569, "ASP-helix/HIS-helix": 0.42, "ASP-helix/HIS-coil": -1.592, "ASP-helix/HIS-sheet": -2.147, "ASP-helix/LEU-helix": -0.288, "ASP-helix/LEU-coil": -0.846, "ASP-helix/LEU-sheet": -4.2, "ASP-helix/ARG-helix": 0.889, "ASP-helix/ARG-coil": -0.812, "ASP-helix/ARG-sheet": -3.11, "ASP-helix/TRP-helix": -0.443, "ASP-helix/TRP-coil": -3.164, "ASP-helix/TRP-sheet": -Infinity, "ASP-helix/ALA-helix": -0.362, "ASP-helix/ALA-coil": -1.683, "ASP-helix/ALA-sheet": -5.064, "ASP-helix/VAL-helix": -0.493, "ASP-helix/VAL-coil": -1.36, "ASP-helix/VAL-sheet": -4.416, "ASP-helix/GLU-helix": 1.148, "ASP-helix/GLU-coil": -0.005, "ASP-helix/GLU-sheet": -Infinity, "ASP-helix/TYR-helix": 0.387, "ASP-helix/TYR-coil": -0.653, "ASP-helix/TYR-sheet": -4.905, "ASP-helix/MET-helix": 0.171, "ASP-helix/MET-coil": -2.214, "ASP-helix/MET-sheet": -Infinity, "ASP-coil/CYS-helix": -Infinity, "ASP-coil/CYS-coil": 1.588, "ASP-coil/CYS-sheet": -0.394, "ASP-coil/ASP-helix": -0.709, "ASP-coil/ASP-coil": 3.464, "ASP-coil/ASP-sheet": 0.348, "ASP-coil/SER-helix": -0.323, "ASP-coil/SER-coil": 1.373, "ASP-coil/SER-sheet": 0.055, "ASP-coil/GLN-helix": -0.07, "ASP-coil/GLN-coil": 0.396, "ASP-coil/GLN-sheet": -0.058, "ASP-coil/LYS-helix": -1.291, "ASP-coil/LYS-coil": 0.804, "ASP-coil/LYS-sheet": 0.467, "ASP-coil/ILE-helix": -0.866, "ASP-coil/ILE-coil": 1.104, "ASP-coil/ILE-sheet": -1.662, "ASP-coil/PRO-helix": -3.204, "ASP-coil/PRO-coil": 0.812, "ASP-coil/PRO-sheet": 0.591, "ASP-coil/THR-helix": -2.114, "ASP-coil/THR-coil": 1.35, "ASP-coil/THR-sheet": -0.09, "ASP-coil/PHE-helix": -4.309, "ASP-coil/PHE-coil": 1.394, "ASP-coil/PHE-sheet": 0.213, "ASP-coil/ASN-helix": -1.269, "ASP-coil/ASN-coil": 0.247, "ASP-coil/ASN-sheet": 1.012, "ASP-coil/GLY-helix": -2.794, "ASP-coil/GLY-coil": 1.351, "ASP-coil/GLY-sheet": -0.931, "ASP-coil/HIS-helix": -1.737, "ASP-coil/HIS-coil": -0.072, "ASP-coil/HIS-sheet": 1.081, "ASP-coil/LEU-helix": -0.833, "ASP-coil/LEU-coil": 0.589, "ASP-coil/LEU-sheet": -0.842, "ASP-coil/ARG-helix": -0.678, "ASP-coil/ARG-coil": 0.201, "ASP-coil/ARG-sheet": -0.14, "ASP-coil/TRP-helix": -0.247, "ASP-coil/TRP-coil": -0.619, "ASP-coil/TRP-sheet": 0.053, "ASP-coil/ALA-helix": -1.554, "ASP-coil/ALA-coil": 0.974, "ASP-coil/ALA-sheet": -0.167, "ASP-coil/VAL-helix": -1.372, "ASP-coil/VAL-coil": 1.54, "ASP-coil/VAL-sheet": -0.806, "ASP-coil/GLU-helix": -2.354, "ASP-coil/GLU-coil": 0.88, "ASP-coil/GLU-sheet": 0.111, "ASP-coil/TYR-helix": -2.257, "ASP-coil/TYR-coil": 0.488, "ASP-coil/TYR-sheet": 0.519, "ASP-coil/MET-helix": -3.979, "ASP-coil/MET-coil": 0.311, "ASP-coil/MET-sheet": 1.588, "ASP-sheet/CYS-helix": -Infinity, "ASP-sheet/CYS-coil": -Infinity, "ASP-sheet/CYS-sheet": 2.546, "ASP-sheet/ASP-helix": 1.056, "ASP-sheet/ASP-coil": 0.348, "ASP-sheet/ASP-sheet": 5.478, "ASP-sheet/SER-helix": -3.474, "ASP-sheet/SER-coil": -1.399, "ASP-sheet/SER-sheet": 3.561, "ASP-sheet/GLN-helix": -4.357, "ASP-sheet/GLN-coil": -3.204, "ASP-sheet/GLN-sheet": 3.257, "ASP-sheet/LYS-helix": -3.583, "ASP-sheet/LYS-coil": 1.621, "ASP-sheet/LYS-sheet": 3.095, "ASP-sheet/ILE-helix": -3.652, "ASP-sheet/ILE-coil": -0.341, "ASP-sheet/ILE-sheet": 1.887, "ASP-sheet/PRO-helix": -2.788, "ASP-sheet/PRO-coil": 0.812, "ASP-sheet/PRO-sheet": -Infinity, "ASP-sheet/THR-helix": -4.154, "ASP-sheet/THR-coil": -0.939, "ASP-sheet/THR-sheet": 1.773, "ASP-sheet/PHE-helix": -Infinity, "ASP-sheet/PHE-coil": -Infinity, "ASP-sheet/PHE-sheet": 2.598, "ASP-sheet/ASN-helix": -3.467, "ASP-sheet/ASN-coil": -1.016, "ASP-sheet/ASN-sheet": 3.06, "ASP-sheet/GLY-helix": -4.196, "ASP-sheet/GLY-coil": 1.075, "ASP-sheet/GLY-sheet": 2.908, "ASP-sheet/HIS-helix": -Infinity, "ASP-sheet/HIS-coil": 0.377, "ASP-sheet/HIS-sheet": 0.061, "ASP-sheet/LEU-helix": -Infinity, "ASP-sheet/LEU-coil": 0.15, "ASP-sheet/LEU-sheet": 1.28, "ASP-sheet/ARG-helix": -2.527, "ASP-sheet/ARG-coil": 0.954, "ASP-sheet/ARG-sheet": 2.456, "ASP-sheet/TRP-helix": -3.259, "ASP-sheet/TRP-coil": -Infinity, "ASP-sheet/TRP-sheet": -1.538, "ASP-sheet/ALA-helix": -4.4, "ASP-sheet/ALA-coil": 1.107, "ASP-sheet/ALA-sheet": 2.468, "ASP-sheet/VAL-helix": -4.457, "ASP-sheet/VAL-coil": -Infinity, "ASP-sheet/VAL-sheet": -0.363, "ASP-sheet/GLU-helix": -5.029, "ASP-sheet/GLU-coil": -4.376, "ASP-sheet/GLU-sheet": 2.753, "ASP-sheet/TYR-helix": -3.019, "ASP-sheet/TYR-coil": -0.976, "ASP-sheet/TYR-sheet": 0.481, "ASP-sheet/MET-helix": -Infinity, "ASP-sheet/MET-coil": -2.086, "ASP-sheet/MET-sheet": 2.496, "SER-helix/CYS-helix": 0.308, "SER-helix/CYS-coil": -0.661, "SER-helix/CYS-sheet": -3.811, "SER-helix/ASP-helix": 1.126, "SER-helix/ASP-coil": -0.323, "SER-helix/ASP-sheet": -3.474, "SER-helix/SER-helix": 3.457, "SER-helix/SER-coil": -0.848, "SER-helix/SER-sheet": -4.375, "SER-helix/GLN-helix": 0.079, "SER-helix/GLN-coil": 0.331, "SER-helix/GLN-sheet": -4.961, "SER-helix/LYS-helix": 0.675, "SER-helix/LYS-coil": -0.442, "SER-helix/LYS-sheet": -4.225, "SER-helix/ILE-helix": 0.355, "SER-helix/ILE-coil": 0.167, "SER-helix/ILE-sheet": -3.622, "SER-helix/PRO-helix": 0.628, "SER-helix/PRO-coil": -1.222, "SER-helix/PRO-sheet": -1.702, "SER-helix/THR-helix": 0.308, "SER-helix/THR-coil": -1.323, "SER-helix/THR-sheet": -3.298, "SER-helix/PHE-helix": -0.058, "SER-helix/PHE-coil": 0.273, "SER-helix/PHE-sheet": -4.839, "SER-helix/ASN-helix": 0.677, "SER-helix/ASN-coil": -1.329, "SER-helix/ASN-sheet": -4.223, "SER-helix/GLY-helix": 1.352, "SER-helix/GLY-coil": -0.459, "SER-helix/GLY-sheet": -4.042, "SER-helix/HIS-helix": 0.958, "SER-helix/HIS-coil": 0.237, "SER-helix/HIS-sheet": -3.313, "SER-helix/LEU-helix": 0.392, "SER-helix/LEU-coil": -0.593, "SER-helix/LEU-sheet": -4.113, "SER-helix/ARG-helix": 0.014, "SER-helix/ARG-coil": -0.908, "SER-helix/ARG-sheet": -2.89, "SER-helix/TRP-helix": -0.338, "SER-helix/TRP-coil": -0.546, "SER-helix/TRP-sheet": -Infinity, "SER-helix/ALA-helix": 1.127, "SER-helix/ALA-coil": -1.194, "SER-helix/ALA-sheet": -3.235, "SER-helix/VAL-helix": 1.025, "SER-helix/VAL-coil": -0.814, "SER-helix/VAL-sheet": -4.042, "SER-helix/GLU-helix": 0.822, "SER-helix/GLU-coil": -0.621, "SER-helix/GLU-sheet": -Infinity, "SER-helix/TYR-helix": -1.274, "SER-helix/TYR-coil": -1.126, "SER-helix/TYR-sheet": -3.432, "SER-helix/MET-helix": 0.003, "SER-helix/MET-coil": -0.753, "SER-helix/MET-sheet": -Infinity, "SER-coil/CYS-helix": -2.024, "SER-coil/CYS-coil": 1.461, "SER-coil/CYS-sheet": -1.876, "SER-coil/ASP-helix": 0.201, "SER-coil/ASP-coil": 1.373, "SER-coil/ASP-sheet": -1.399, "SER-coil/SER-helix": -0.848, "SER-coil/SER-coil": 3.367, "SER-coil/SER-sheet": -0.729, "SER-coil/GLN-helix": -1.327, "SER-coil/GLN-coil": 0.92, "SER-coil/GLN-sheet": -0.922, "SER-coil/LYS-helix": 0.271, "SER-coil/LYS-coil": 0.952, "SER-coil/LYS-sheet": -1.028, "SER-coil/ILE-helix": -1.479, "SER-coil/ILE-coil": 1.213, "SER-coil/ILE-sheet": -2.198, "SER-coil/PRO-helix": -2.909, "SER-coil/PRO-coil": 0.309, "SER-coil/PRO-sheet": -2.676, "SER-coil/THR-helix": -2.032, "SER-coil/THR-coil": 0.351, "SER-coil/THR-sheet": 0.122, "SER-coil/PHE-helix": -1.622, "SER-coil/PHE-coil": 0.656, "SER-coil/PHE-sheet": -0.196, "SER-coil/ASN-helix": 0.121, "SER-coil/ASN-coil": 1.311, "SER-coil/ASN-sheet": 0.625, "SER-coil/GLY-helix": -1.364, "SER-coil/GLY-coil": 1.466, "SER-coil/GLY-sheet": -1.844, "SER-coil/HIS-helix": -2.577, "SER-coil/HIS-coil": 0.011, "SER-coil/HIS-sheet": -1.195, "SER-coil/LEU-helix": -0.708, "SER-coil/LEU-coil": 1.039, "SER-coil/LEU-sheet": 0.24, "SER-coil/ARG-helix": -0.962, "SER-coil/ARG-coil": 0.323, "SER-coil/ARG-sheet": -0.353, "SER-coil/TRP-helix": -2.479, "SER-coil/TRP-coil": -0.848, "SER-coil/TRP-sheet": 0.679, "SER-coil/ALA-helix": -0.676, "SER-coil/ALA-coil": 1.11, "SER-coil/ALA-sheet": 0.379, "SER-coil/VAL-helix": -0.937, "SER-coil/VAL-coil": 0.841, "SER-coil/VAL-sheet": -2.025, "SER-coil/GLU-helix": -1.343, "SER-coil/GLU-coil": 1.07, "SER-coil/GLU-sheet": -0.027, "SER-coil/TYR-helix": -1.207, "SER-coil/TYR-coil": 0.59, "SER-coil/TYR-sheet": -0.055, "SER-coil/MET-helix": -2.368, "SER-coil/MET-coil": 0.487, "SER-coil/MET-sheet": -2.065, "SER-sheet/CYS-helix": -3.055, "SER-sheet/CYS-coil": -2.843, "SER-sheet/CYS-sheet": 0.741, "SER-sheet/ASP-helix": -3.209, "SER-sheet/ASP-coil": 0.055, "SER-sheet/ASP-sheet": 3.561, "SER-sheet/SER-helix": -4.375, "SER-sheet/SER-coil": -0.729, "SER-sheet/SER-sheet": 4.982, "SER-sheet/GLN-helix": -3.584, "SER-sheet/GLN-coil": -0.521, "SER-sheet/GLN-sheet": 1.469, "SER-sheet/LYS-helix": -3.909, "SER-sheet/LYS-coil": -1.656, "SER-sheet/LYS-sheet": 2.672, "SER-sheet/ILE-helix": -4.958, "SER-sheet/ILE-coil": -2.118, "SER-sheet/ILE-sheet": 1.652, "SER-sheet/PRO-helix": -2.708, "SER-sheet/PRO-coil": -0.692, "SER-sheet/PRO-sheet": -1.34, "SER-sheet/THR-helix": -4.362, "SER-sheet/THR-coil": 0.091, "SER-sheet/THR-sheet": 2.748, "SER-sheet/PHE-helix": -3.813, "SER-sheet/PHE-coil": -1.875, "SER-sheet/PHE-sheet": 1.722, "SER-sheet/ASN-helix": -4.081, "SER-sheet/ASN-coil": -3.015, "SER-sheet/ASN-sheet": 3.779, "SER-sheet/GLY-helix": -4.404, "SER-sheet/GLY-coil": 0.62, "SER-sheet/GLY-sheet": 3.637, "SER-sheet/HIS-helix": -4.236, "SER-sheet/HIS-coil": -2.579, "SER-sheet/HIS-sheet": 2.338, "SER-sheet/LEU-helix": -5.58, "SER-sheet/LEU-coil": 0.594, "SER-sheet/LEU-sheet": 2.486, "SER-sheet/ARG-helix": -3.246, "SER-sheet/ARG-coil": -2.368, "SER-sheet/ARG-sheet": 2.207, "SER-sheet/TRP-helix": -4.566, "SER-sheet/TRP-coil": -2.551, "SER-sheet/TRP-sheet": 0.118, "SER-sheet/ALA-helix": -3.473, "SER-sheet/ALA-coil": -1.599, "SER-sheet/ALA-sheet": 2.691, "SER-sheet/VAL-helix": -3.818, "SER-sheet/VAL-coil": -0.004, "SER-sheet/VAL-sheet": 0.174, "SER-sheet/GLU-helix": -2.464, "SER-sheet/GLU-coil": -2.592, "SER-sheet/GLU-sheet": 2.029, "SER-sheet/TYR-helix": -2.939, "SER-sheet/TYR-coil": -3.381, "SER-sheet/TYR-sheet": 2.121, "SER-sheet/MET-helix": -4.687, "SER-sheet/MET-coil": -1.852, "SER-sheet/MET-sheet": 2.428, "GLN-helix/CYS-helix": -0.857, "GLN-helix/CYS-coil": -1.516, "GLN-helix/CYS-sheet": -2.615, "GLN-helix/ASP-helix": 0.559, "GLN-helix/ASP-coil": -0.07, "GLN-helix/ASP-sheet": -4.357, "GLN-helix/SER-helix": 0.079, "GLN-helix/SER-coil": -1.327, "GLN-helix/SER-sheet": -3.584, "GLN-helix/GLN-helix": 3.935, "GLN-helix/GLN-coil": 0.935, "GLN-helix/GLN-sheet": -2.666, "GLN-helix/LYS-helix": 0.574, "GLN-helix/LYS-coil": -0.651, "GLN-helix/LYS-sheet": -4.127, "GLN-helix/ILE-helix": -0.795, "GLN-helix/ILE-coil": 0.891, "GLN-helix/ILE-sheet": -3.679, "GLN-helix/PRO-helix": -0.67, "GLN-helix/PRO-coil": -0.833, "GLN-helix/PRO-sheet": -2.211, "GLN-helix/THR-helix": 0.025, "GLN-helix/THR-coil": 0.118, "GLN-helix/THR-sheet": -4.181, "GLN-helix/PHE-helix": 0.096, "GLN-helix/PHE-coil": 0.501, "GLN-helix/PHE-sheet": -3.642, "GLN-helix/ASN-helix": 0.662, "GLN-helix/ASN-coil": -3.039, "GLN-helix/ASN-sheet": -Infinity, "GLN-helix/GLY-helix": -1.453, "GLN-helix/GLY-coil": -1.751, "GLN-helix/GLY-sheet": -3.826, "GLN-helix/HIS-helix": 1.5, "GLN-helix/HIS-coil": -1.763, "GLN-helix/HIS-sheet": -Infinity, "GLN-helix/LEU-helix": 0.479, "GLN-helix/LEU-coil": 0.032, "GLN-helix/LEU-sheet": -0.319, "GLN-helix/ARG-helix": 1.044, "GLN-helix/ARG-coil": -1.406, "GLN-helix/ARG-sheet": -3.079, "GLN-helix/TRP-helix": 0.812, "GLN-helix/TRP-coil": -2.728, "GLN-helix/TRP-sheet": -3.427, "GLN-helix/ALA-helix": 0.524, "GLN-helix/ALA-coil": -0.055, "GLN-helix/ALA-sheet": -3.425, "GLN-helix/VAL-helix": 0.067, "GLN-helix/VAL-coil": -1.452, "GLN-helix/VAL-sheet": -3.693, "GLN-helix/GLU-helix": 0.588, "GLN-helix/GLU-coil": -1.556, "GLN-helix/GLU-sheet": 1.087, "GLN-helix/TYR-helix": 0.13, "GLN-helix/TYR-coil": 0.226, "GLN-helix/TYR-sheet": -4.181, "GLN-helix/MET-helix": 0.429, "GLN-helix/MET-coil": -0.024, "GLN-helix/MET-sheet": -2.58, "GLN-coil/CYS-helix": 0.0, "GLN-coil/CYS-coil": 0.212, "GLN-coil/CYS-sheet": -Infinity, "GLN-coil/ASP-helix": -1.8, "GLN-coil/ASP-coil": 0.396, "GLN-coil/ASP-sheet": -3.204, "GLN-coil/SER-helix": 0.331, "GLN-coil/SER-coil": 0.92, "GLN-coil/SER-sheet": -0.521, "GLN-coil/GLN-helix": 0.935, "GLN-coil/GLN-coil": 4.232, "GLN-coil/GLN-sheet": -1.919, "GLN-coil/LYS-helix": 0.647, "GLN-coil/LYS-coil": 1.015, "GLN-coil/LYS-sheet": -3.262, "GLN-coil/ILE-helix": -0.545, "GLN-coil/ILE-coil": 0.918, "GLN-coil/ILE-sheet": -2.996, "GLN-coil/PRO-helix": -1.039, "GLN-coil/PRO-coil": 0.604, "GLN-coil/PRO-sheet": -3.137, "GLN-coil/THR-helix": -0.609, "GLN-coil/THR-coil": 0.312, "GLN-coil/THR-sheet": -3.433, "GLN-coil/PHE-helix": -1.347, "GLN-coil/PHE-coil": 0.295, "GLN-coil/PHE-sheet": -2.266, "GLN-coil/ASN-helix": -0.194, "GLN-coil/ASN-coil": -0.019, "GLN-coil/ASN-sheet": 1.003, "GLN-coil/GLY-helix": -0.21, "GLN-coil/GLY-coil": 0.828, "GLN-coil/GLY-sheet": -2.386, "GLN-coil/HIS-helix": 0.044, "GLN-coil/HIS-coil": -0.61, "GLN-coil/HIS-sheet": -1.944, "GLN-coil/LEU-helix": -0.548, "GLN-coil/LEU-coil": 0.778, "GLN-coil/LEU-sheet": -2.323, "GLN-coil/ARG-helix": 0.349, "GLN-coil/ARG-coil": 1.048, "GLN-coil/ARG-sheet": 1.063, "GLN-coil/TRP-helix": -1.366, "GLN-coil/TRP-coil": -1.17, "GLN-coil/TRP-sheet": -2.562, "GLN-coil/ALA-helix": -0.211, "GLN-coil/ALA-coil": 1.241, "GLN-coil/ALA-sheet": 1.668, "GLN-coil/VAL-helix": 0.219, "GLN-coil/VAL-coil": 0.487, "GLN-coil/VAL-sheet": -1.934, "GLN-coil/GLU-helix": 0.168, "GLN-coil/GLU-coil": 0.934, "GLN-coil/GLU-sheet": -1.099, "GLN-coil/TYR-helix": -0.116, "GLN-coil/TYR-coil": 1.178, "GLN-coil/TYR-sheet": -3.028, "GLN-coil/MET-helix": -1.314, "GLN-coil/MET-coil": 0.109, "GLN-coil/MET-sheet": -1.833, "GLN-sheet/CYS-helix": -Infinity, "GLN-sheet/CYS-coil": -Infinity, "GLN-sheet/CYS-sheet": 3.024, "GLN-sheet/ASP-helix": -3.795, "GLN-sheet/ASP-coil": -0.058, "GLN-sheet/ASP-sheet": 3.257, "GLN-sheet/SER-helix": -4.961, "GLN-sheet/SER-coil": -0.922, "GLN-sheet/SER-sheet": 1.469, "GLN-sheet/GLN-helix": -2.666, "GLN-sheet/GLN-coil": -1.919, "GLN-sheet/GLN-sheet": 5.214, "GLN-sheet/LYS-helix": -4.783, "GLN-sheet/LYS-coil": 0.592, "GLN-sheet/LYS-sheet": 2.749, "GLN-sheet/ILE-helix": -5.139, "GLN-sheet/ILE-coil": 0.152, "GLN-sheet/ILE-sheet": 1.027, "GLN-sheet/PRO-helix": -3.988, "GLN-sheet/PRO-coil": -0.029, "GLN-sheet/PRO-sheet": -0.605, "GLN-sheet/THR-helix": -4.255, "GLN-sheet/THR-coil": 0.302, "GLN-sheet/THR-sheet": 1.958, "GLN-sheet/PHE-helix": -Infinity, "GLN-sheet/PHE-coil": -3.309, "GLN-sheet/PHE-sheet": 1.748, "GLN-sheet/ASN-helix": -4.262, "GLN-sheet/ASN-coil": -1.117, "GLN-sheet/ASN-sheet": 2.391, "GLN-sheet/GLY-helix": 0.127, "GLN-sheet/GLY-coil": 0.02, "GLN-sheet/GLY-sheet": 2.973, "GLN-sheet/HIS-helix": -3.724, "GLN-sheet/HIS-coil": 0.731, "GLN-sheet/HIS-sheet": 1.569, "GLN-sheet/LEU-helix": -4.375, "GLN-sheet/LEU-coil": 0.028, "GLN-sheet/LEU-sheet": 0.462, "GLN-sheet/ARG-helix": -4.238, "GLN-sheet/ARG-coil": -2.954, "GLN-sheet/ARG-sheet": 2.57, "GLN-sheet/TRP-helix": -Infinity, "GLN-sheet/TRP-coil": -2.444, "GLN-sheet/TRP-sheet": -0.54, "GLN-sheet/ALA-helix": -5.6, "GLN-sheet/ALA-coil": 1.398, "GLN-sheet/ALA-sheet": 2.885, "GLN-sheet/VAL-helix": -0.809, "GLN-sheet/VAL-coil": -1.008, "GLN-sheet/VAL-sheet": -0.364, "GLN-sheet/GLU-helix": -1.492, "GLN-sheet/GLU-coil": -1.387, "GLN-sheet/GLU-sheet": 2.222, "GLN-sheet/TYR-helix": -Infinity, "GLN-sheet/TYR-coil": -3.967, "GLN-sheet/TYR-sheet": 2.574, "GLN-sheet/MET-helix": -Infinity, "GLN-sheet/MET-coil": 0.148, "GLN-sheet/MET-sheet": 1.2, "LYS-helix/CYS-helix": -2.304, "LYS-helix/CYS-coil": -1.399, "LYS-helix/CYS-sheet": -Infinity, "LYS-helix/ASP-helix": 1.165, "LYS-helix/ASP-coil": -1.291, "LYS-helix/ASP-sheet": -3.583, "LYS-helix/SER-helix": 0.675, "LYS-helix/SER-coil": 0.271, "LYS-helix/SER-sheet": -3.909, "LYS-helix/GLN-helix": 0.574, "LYS-helix/GLN-coil": 0.647, "LYS-helix/GLN-sheet": -4.783, "LYS-helix/LYS-helix": 3.649, "LYS-helix/LYS-coil": -0.709, "LYS-helix/LYS-sheet": -Infinity, "LYS-helix/ILE-helix": -0.243, "LYS-helix/ILE-coil": -1.339, "LYS-helix/ILE-sheet": -3.781, "LYS-helix/PRO-helix": -0.377, "LYS-helix/PRO-coil": -1.863, "LYS-helix/PRO-sheet": -2.823, "LYS-helix/THR-helix": -0.538, "LYS-helix/THR-coil": -1.988, "LYS-helix/THR-sheet": -3.813, "LYS-helix/PHE-helix": -0.319, "LYS-helix/PHE-coil": -1.115, "LYS-helix/PHE-sheet": -5.353, "LYS-helix/ASN-helix": 1.5, "LYS-helix/ASN-coil": -0.404, "LYS-helix/ASN-sheet": -Infinity, "LYS-helix/GLY-helix": -0.756, "LYS-helix/GLY-coil": -1.15, "LYS-helix/GLY-sheet": -Infinity, "LYS-helix/HIS-helix": 1.001, "LYS-helix/HIS-coil": -0.799, "LYS-helix/HIS-sheet": -3.135, "LYS-helix/LEU-helix": -0.006, "LYS-helix/LEU-coil": -1.028, "LYS-helix/LEU-sheet": -5.881, "LYS-helix/ARG-helix": 0.488, "LYS-helix/ARG-coil": -0.907, "LYS-helix/ARG-sheet": -4.791, "LYS-helix/TRP-helix": 0.596, "LYS-helix/TRP-coil": -3.053, "LYS-helix/TRP-sheet": -3.346, "LYS-helix/ALA-helix": -0.484, "LYS-helix/ALA-coil": -0.676, "LYS-helix/ALA-sheet": -4.666, "LYS-helix/VAL-helix": -0.564, "LYS-helix/VAL-coil": -1.241, "LYS-helix/VAL-sheet": -2.408, "LYS-helix/GLU-helix": 1.122, "LYS-helix/GLU-coil": -0.348, "LYS-helix/GLU-sheet": -3.136, "LYS-helix/TYR-helix": 0.126, "LYS-helix/TYR-coil": -1.868, "LYS-helix/TYR-sheet": -4.506, "LYS-helix/MET-helix": -0.531, "LYS-helix/MET-coil": -1.042, "LYS-helix/MET-sheet": -0.672, "LYS-coil/CYS-helix": -4.6, "LYS-coil/CYS-coil": -3.694, "LYS-coil/CYS-sheet": -2.19, "LYS-coil/ASP-helix": -2.023, "LYS-coil/ASP-coil": 0.804, "LYS-coil/ASP-sheet": 1.621, "LYS-coil/SER-helix": -0.442, "LYS-coil/SER-coil": 0.952, "LYS-coil/SER-sheet": -1.656, "LYS-coil/GLN-helix": -0.651, "LYS-coil/GLN-coil": 1.015, "LYS-coil/GLN-sheet": 0.592, "LYS-coil/LYS-helix": -0.709, "LYS-coil/LYS-coil": 3.481, "LYS-coil/LYS-sheet": 0.175, "LYS-coil/ILE-helix": -2.026, "LYS-coil/ILE-coil": 0.496, "LYS-coil/ILE-sheet": -0.392, "LYS-coil/PRO-helix": -1.075, "LYS-coil/PRO-coil": 0.336, "LYS-coil/PRO-sheet": -2.192, "LYS-coil/THR-helix": -1.171, "LYS-coil/THR-coil": 0.585, "LYS-coil/THR-sheet": -0.957, "LYS-coil/PHE-helix": -1.694, "LYS-coil/PHE-coil": 0.861, "LYS-coil/PHE-sheet": -0.591, "LYS-coil/ASN-helix": -0.422, "LYS-coil/ASN-coil": 0.793, "LYS-coil/ASN-sheet": 0.526, "LYS-coil/GLY-helix": -2.799, "LYS-coil/GLY-coil": 0.582, "LYS-coil/GLY-sheet": -1.192, "LYS-coil/HIS-helix": -1.156, "LYS-coil/HIS-coil": -0.465, "LYS-coil/HIS-sheet": -1.404, "LYS-coil/LEU-helix": -1.972, "LYS-coil/LEU-coil": 0.686, "LYS-coil/LEU-sheet": -1.274, "LYS-coil/ARG-helix": -1.443, "LYS-coil/ARG-coil": 0.201, "LYS-coil/ARG-sheet": -0.709, "LYS-coil/TRP-helix": -2.199, "LYS-coil/TRP-coil": -2.555, "LYS-coil/TRP-sheet": -0.7, "LYS-coil/ALA-helix": -1.864, "LYS-coil/ALA-coil": 0.608, "LYS-coil/ALA-sheet": -0.826, "LYS-coil/VAL-helix": -3.248, "LYS-coil/VAL-coil": 0.769, "LYS-coil/VAL-sheet": -2.421, "LYS-coil/GLU-helix": -0.358, "LYS-coil/GLU-coil": 0.85, "LYS-coil/GLU-sheet": 0.913, "LYS-coil/TYR-helix": -0.14, "LYS-coil/TYR-coil": 0.214, "LYS-coil/TYR-sheet": -0.068, "LYS-coil/MET-helix": -2.32, "LYS-coil/MET-coil": -0.106, "LYS-coil/MET-sheet": -1.868, "LYS-sheet/CYS-helix": -1.806, "LYS-sheet/CYS-coil": -2.693, "LYS-sheet/CYS-sheet": 2.5, "LYS-sheet/ASP-helix": -2.365, "LYS-sheet/ASP-coil": 0.467, "LYS-sheet/ASP-sheet": 3.095, "LYS-sheet/SER-helix": -4.225, "LYS-sheet/SER-coil": -1.028, "LYS-sheet/SER-sheet": 2.672, "LYS-sheet/GLN-helix": -4.127, "LYS-sheet/GLN-coil": -3.262, "LYS-sheet/GLN-sheet": 2.749, "LYS-sheet/LYS-helix": -Infinity, "LYS-sheet/LYS-coil": 0.175, "LYS-sheet/LYS-sheet": 4.655, "LYS-sheet/ILE-helix": -Infinity, "LYS-sheet/ILE-coil": 0.786, "LYS-sheet/ILE-sheet": 0.83, "LYS-sheet/PRO-helix": -4.35, "LYS-sheet/PRO-coil": -0.28, "LYS-sheet/PRO-sheet": -2.577, "LYS-sheet/THR-helix": -5.311, "LYS-sheet/THR-coil": 0.668, "LYS-sheet/THR-sheet": 2.662, "LYS-sheet/PHE-helix": -4.356, "LYS-sheet/PHE-coil": -0.473, "LYS-sheet/PHE-sheet": 1.228, "LYS-sheet/ASN-helix": -1.985, "LYS-sheet/ASN-coil": -2.578, "LYS-sheet/ASN-sheet": 2.184, "LYS-sheet/GLY-helix": -3.156, "LYS-sheet/GLY-coil": 1.18, "LYS-sheet/GLY-sheet": 2.697, "LYS-sheet/HIS-helix": -4.087, "LYS-sheet/HIS-coil": 0.484, "LYS-sheet/HIS-sheet": 0.226, "LYS-sheet/LEU-helix": -3.558, "LYS-sheet/LEU-coil": 0.131, "LYS-sheet/LEU-sheet": 1.975, "LYS-sheet/ARG-helix": -Infinity, "LYS-sheet/ARG-coil": 0.267, "LYS-sheet/ARG-sheet": -0.187, "LYS-sheet/TRP-helix": -3.723, "LYS-sheet/TRP-coil": -1.554, "LYS-sheet/TRP-sheet": 0.832, "LYS-sheet/ALA-helix": -3.565, "LYS-sheet/ALA-coil": -0.386, "LYS-sheet/ALA-sheet": 2.248, "LYS-sheet/VAL-helix": -1.852, "LYS-sheet/VAL-coil": -0.572, "LYS-sheet/VAL-sheet": 1.301, "LYS-sheet/GLU-helix": -3.883, "LYS-sheet/GLU-coil": -0.186, "LYS-sheet/GLU-sheet": 3.335, "LYS-sheet/TYR-helix": -4.869, "LYS-sheet/TYR-coil": -0.546, "LYS-sheet/TYR-sheet": 1.816, "LYS-sheet/MET-helix": -Infinity, "LYS-sheet/MET-coil": -0.47, "LYS-sheet/MET-sheet": 0.674, "ILE-helix/CYS-helix": 0.844, "ILE-helix/CYS-coil": -0.146, "ILE-helix/CYS-sheet": -3.296, "ILE-helix/ASP-helix": -0.224, "ILE-helix/ASP-coil": -0.866, "ILE-helix/ASP-sheet": -3.652, "ILE-helix/SER-helix": 0.355, "ILE-helix/SER-coil": -1.479, "ILE-helix/SER-sheet": -4.958, "ILE-helix/GLN-helix": -0.795, "ILE-helix/GLN-coil": -0.545, "ILE-helix/GLN-sheet": -5.139, "ILE-helix/LYS-helix": -0.243, "ILE-helix/LYS-coil": -2.026, "ILE-helix/LYS-sheet": -Infinity, "ILE-helix/ILE-helix": 3.359, "ILE-helix/ILE-coil": -2.211, "ILE-helix/ILE-sheet": -3.348, "ILE-helix/PRO-helix": -0.229, "ILE-helix/PRO-coil": -1.279, "ILE-helix/PRO-sheet": -3.18, "ILE-helix/THR-helix": 0.565, "ILE-helix/THR-coil": -1.335, "ILE-helix/THR-sheet": -4.015, "ILE-helix/PHE-helix": -0.777, "ILE-helix/PHE-coil": -2.931, "ILE-helix/PHE-sheet": -3.407, "ILE-helix/ASN-helix": -0.612, "ILE-helix/ASN-coil": -2.804, "ILE-helix/ASN-sheet": -3.995, "ILE-helix/GLY-helix": 1.022, "ILE-helix/GLY-coil": -2.484, "ILE-helix/GLY-sheet": -3.814, "ILE-helix/HIS-helix": -1.694, "ILE-helix/HIS-coil": -1.145, "ILE-helix/HIS-sheet": -Infinity, "ILE-helix/LEU-helix": 0.42, "ILE-helix/LEU-coil": -2.19, "ILE-helix/LEU-sheet": -4.157, "ILE-helix/ARG-helix": -0.598, "ILE-helix/ARG-coil": -1.84, "ILE-helix/ARG-sheet": -4.454, "ILE-helix/TRP-helix": -0.444, "ILE-helix/TRP-coil": -3.004, "ILE-helix/TRP-sheet": 0.728, "ILE-helix/ALA-helix": 0.562, "ILE-helix/ALA-coil": -0.58, "ILE-helix/ALA-sheet": -0.231, "ILE-helix/VAL-helix": 0.859, "ILE-helix/VAL-coil": -0.684, "ILE-helix/VAL-sheet": -3.221, "ILE-helix/GLU-helix": -0.433, "ILE-helix/GLU-coil": -1.199, "ILE-helix/GLU-sheet": -4.591, "ILE-helix/TYR-helix": -0.5, "ILE-helix/TYR-coil": -1.714, "ILE-helix/TYR-sheet": -3.071, "ILE-helix/MET-helix": 0.815, "ILE-helix/MET-coil": -0.031, "ILE-helix/MET-sheet": -3.667, "ILE-coil/CYS-helix": -2.97, "ILE-coil/CYS-coil": -0.219, "ILE-coil/CYS-sheet": -1.372, "ILE-coil/ASP-helix": -0.568, "ILE-coil/ASP-coil": 1.104, "ILE-coil/ASP-sheet": -0.341, "ILE-coil/SER-helix": 0.167, "ILE-coil/SER-coil": 1.213, "ILE-coil/SER-sheet": -2.118, "ILE-coil/GLN-helix": 0.891, "ILE-coil/GLN-coil": 0.918, "ILE-coil/GLN-sheet": 0.152, "ILE-coil/LYS-helix": -1.339, "ILE-coil/LYS-coil": 0.496, "ILE-coil/LYS-sheet": 0.786, "ILE-coil/ILE-helix": -2.211, "ILE-coil/ILE-coil": 4.214, "ILE-coil/ILE-sheet": 0.212, "ILE-coil/PRO-helix": -1.453, "ILE-coil/PRO-coil": 1.118, "ILE-coil/PRO-sheet": -1.949, "ILE-coil/THR-helix": -1.519, "ILE-coil/THR-coil": 0.483, "ILE-coil/THR-sheet": -0.838, "ILE-coil/PHE-helix": -0.998, "ILE-coil/PHE-coil": 0.417, "ILE-coil/PHE-sheet": 0.742, "ILE-coil/ASN-helix": -2.53, "ILE-coil/ASN-coil": 0.569, "ILE-coil/ASN-sheet": 1.141, "ILE-coil/GLY-helix": 1.04, "ILE-coil/GLY-coil": 1.342, "ILE-coil/GLY-sheet": -0.662, "ILE-coil/HIS-helix": -1.613, "ILE-coil/HIS-coil": 0.461, "ILE-coil/HIS-sheet": -2.26, "ILE-coil/LEU-helix": -0.753, "ILE-coil/LEU-coil": 0.159, "ILE-coil/LEU-sheet": 0.913, "ILE-coil/ARG-helix": -1.045, "ILE-coil/ARG-coil": 0.29, "ILE-coil/ARG-sheet": -0.004, "ILE-coil/TRP-helix": 0.182, "ILE-coil/TRP-coil": 0.306, "ILE-coil/TRP-sheet": -Infinity, "ILE-coil/ALA-helix": 0.204, "ILE-coil/ALA-coil": 0.473, "ILE-coil/ALA-sheet": -1.651, "ILE-coil/VAL-helix": 0.14, "ILE-coil/VAL-coil": 1.894, "ILE-coil/VAL-sheet": -1.351, "ILE-coil/GLU-helix": -2.827, "ILE-coil/GLU-coil": 0.621, "ILE-coil/GLU-sheet": -0.416, "ILE-coil/TYR-helix": -0.996, "ILE-coil/TYR-coil": -1.087, "ILE-coil/TYR-sheet": 0.781, "ILE-coil/MET-helix": -3.062, "ILE-coil/MET-coil": 0.582, "ILE-coil/MET-sheet": -1.743, "ILE-sheet/CYS-helix": -2.457, "ILE-sheet/CYS-coil": -1.839, "ILE-sheet/CYS-sheet": 1.642, "ILE-sheet/ASP-helix": -4.402, "ILE-sheet/ASP-coil": -1.662, "ILE-sheet/ASP-sheet": 1.887, "ILE-sheet/SER-helix": -3.622, "ILE-sheet/SER-coil": -2.198, "ILE-sheet/SER-sheet": 1.652, "ILE-sheet/GLN-helix": -3.679, "ILE-sheet/GLN-coil": -2.996, "ILE-sheet/GLN-sheet": 1.027, "ILE-sheet/LYS-helix": -3.781, "ILE-sheet/LYS-coil": -0.392, "ILE-sheet/LYS-sheet": 0.83, "ILE-sheet/ILE-helix": -3.348, "ILE-sheet/ILE-coil": 0.212, "ILE-sheet/ILE-sheet": 4.35, "ILE-sheet/PRO-helix": -2.986, "ILE-sheet/PRO-coil": -1.144, "ILE-sheet/PRO-sheet": -Infinity, "ILE-sheet/THR-helix": 0.009, "ILE-sheet/THR-coil": -0.584, "ILE-sheet/THR-sheet": 1.529, "ILE-sheet/PHE-helix": -3.685, "ILE-sheet/PHE-coil": -2.212, "ILE-sheet/PHE-sheet": 2.093, "ILE-sheet/ASN-helix": -3.077, "ILE-sheet/ASN-coil": -0.97, "ILE-sheet/ASN-sheet": 0.895, "ILE-sheet/GLY-helix": -3.519, "ILE-sheet/GLY-coil": -0.608, "ILE-sheet/GLY-sheet": 2.825, "ILE-sheet/HIS-helix": -4.331, "ILE-sheet/HIS-coil": -3.185, "ILE-sheet/HIS-sheet": 0.269, "ILE-sheet/LEU-helix": -0.252, "ILE-sheet/LEU-coil": -0.064, "ILE-sheet/LEU-sheet": 1.775, "ILE-sheet/ARG-helix": -2.594, "ILE-sheet/ARG-coil": -3.562, "ILE-sheet/ARG-sheet": 1.474, "ILE-sheet/TRP-helix": -3.052, "ILE-sheet/TRP-coil": -3.051, "ILE-sheet/TRP-sheet": -0.049, "ILE-sheet/ALA-helix": -0.249, "ILE-sheet/ALA-coil": -0.079, "ILE-sheet/ALA-sheet": 1.877, "ILE-sheet/VAL-helix": -0.231, "ILE-sheet/VAL-coil": -1.967, "ILE-sheet/VAL-sheet": 2.077, "ILE-sheet/GLU-helix": -4.128, "ILE-sheet/GLU-coil": -2.687, "ILE-sheet/GLU-sheet": 0.316, "ILE-sheet/TYR-helix": -3.727, "ILE-sheet/TYR-coil": -3.476, "ILE-sheet/TYR-sheet": 1.037, "ILE-sheet/MET-helix": -2.703, "ILE-sheet/MET-coil": 0.476, "ILE-sheet/MET-sheet": 2.354, "PRO-helix/CYS-helix": 0.457, "PRO-helix/CYS-coil": -1.634, "PRO-helix/CYS-sheet": -1.739, "PRO-helix/ASP-helix": -0.364, "PRO-helix/ASP-coil": -3.204, "PRO-helix/ASP-sheet": -2.788, "PRO-helix/SER-helix": 0.628, "PRO-helix/SER-coil": -2.909, "PRO-helix/SER-sheet": -2.708, "PRO-helix/GLN-helix": -0.67, "PRO-helix/GLN-coil": -1.039, "PRO-helix/GLN-sheet": -3.988, "PRO-helix/LYS-helix": -0.377, "PRO-helix/LYS-coil": -1.075, "PRO-helix/LYS-sheet": -4.35, "PRO-helix/ILE-helix": -0.229, "PRO-helix/ILE-coil": -1.453, "PRO-helix/ILE-sheet": -2.986, "PRO-helix/PRO-helix": 4.238, "PRO-helix/PRO-coil": -3.076, "PRO-helix/PRO-sheet": -Infinity, "PRO-helix/THR-helix": 0.58, "PRO-helix/THR-coil": 0.442, "PRO-helix/THR-sheet": -4.116, "PRO-helix/PHE-helix": 0.07, "PRO-helix/PHE-coil": -1.043, "PRO-helix/PHE-sheet": -2.949, "PRO-helix/ASN-helix": -0.233, "PRO-helix/ASN-coil": -3.821, "PRO-helix/ASN-sheet": -2.556, "PRO-helix/GLY-helix": 0.721, "PRO-helix/GLY-coil": -0.186, "PRO-helix/GLY-sheet": -2.845, "PRO-helix/HIS-helix": -0.23, "PRO-helix/HIS-coil": -2.509, "PRO-helix/HIS-sheet": -Infinity, "PRO-helix/LEU-helix": 0.275, "PRO-helix/LEU-coil": -1.392, "PRO-helix/LEU-sheet": -1.589, "PRO-helix/ARG-helix": -0.414, "PRO-helix/ARG-coil": -1.228, "PRO-helix/ARG-sheet": -Infinity, "PRO-helix/TRP-helix": -0.676, "PRO-helix/TRP-coil": 0.405, "PRO-helix/TRP-sheet": -2.264, "PRO-helix/ALA-helix": 0.291, "PRO-helix/ALA-coil": -1.999, "PRO-helix/ALA-sheet": -3.178, "PRO-helix/VAL-helix": -0.082, "PRO-helix/VAL-coil": -1.69, "PRO-helix/VAL-sheet": -3.405, "PRO-helix/GLU-helix": -0.967, "PRO-helix/GLU-coil": -1.313, "PRO-helix/GLU-sheet": 0.4, "PRO-helix/TYR-helix": 0.251, "PRO-helix/TYR-coil": -3.271, "PRO-helix/TYR-sheet": -2.458, "PRO-helix/MET-helix": -0.513, "PRO-helix/MET-coil": -2.589, "PRO-helix/MET-sheet": -2.516, "PRO-coil/CYS-helix": -1.587, "PRO-coil/CYS-coil": -0.18, "PRO-coil/CYS-sheet": -1.109, "PRO-coil/ASP-helix": 0.238, "PRO-coil/ASP-coil": 0.812, "PRO-coil/ASP-sheet": 0.812, "PRO-coil/SER-helix": -1.222, "PRO-coil/SER-coil": 0.309, "PRO-coil/SER-sheet": -0.692, "PRO-coil/GLN-helix": -0.833, "PRO-coil/GLN-coil": 0.604, "PRO-coil/GLN-sheet": -0.029, "PRO-coil/LYS-helix": -1.863, "PRO-coil/LYS-coil": 0.336, "PRO-coil/LYS-sheet": -0.28, "PRO-coil/ILE-helix": -1.279, "PRO-coil/ILE-coil": 1.118, "PRO-coil/ILE-sheet": -1.144, "PRO-coil/PRO-helix": -3.076, "PRO-coil/PRO-coil": 3.429, "PRO-coil/PRO-sheet": -2.092, "PRO-coil/THR-helix": -2.503, "PRO-coil/THR-coil": 0.75, "PRO-coil/THR-sheet": 0.572, "PRO-coil/PHE-helix": -1.197, "PRO-coil/PHE-coil": 0.819, "PRO-coil/PHE-sheet": -0.008, "PRO-coil/ASN-helix": -2.138, "PRO-coil/ASN-coil": 0.14, "PRO-coil/ASN-sheet": -1.683, "PRO-coil/GLY-helix": -1.823, "PRO-coil/GLY-coil": 1.231, "PRO-coil/GLY-sheet": -1.391, "PRO-coil/HIS-helix": -2.646, "PRO-coil/HIS-coil": -0.063, "PRO-coil/HIS-sheet": -1.815, "PRO-coil/LEU-helix": -1.246, "PRO-coil/LEU-coil": 1.241, "PRO-coil/LEU-sheet": -0.74, "PRO-coil/ARG-helix": -1.717, "PRO-coil/ARG-coil": -0.554, "PRO-coil/ARG-sheet": 0.421, "PRO-coil/TRP-helix": 0.526, "PRO-coil/TRP-coil": 0.42, "PRO-coil/TRP-sheet": 0.503, "PRO-coil/ALA-helix": -1.463, "PRO-coil/ALA-coil": 1.56, "PRO-coil/ALA-sheet": -0.523, "PRO-coil/VAL-helix": -2.381, "PRO-coil/VAL-coil": 0.839, "PRO-coil/VAL-sheet": -0.495, "PRO-coil/GLU-helix": -2.953, "PRO-coil/GLU-coil": 0.738, "PRO-coil/GLU-sheet": -2.142, "PRO-coil/TYR-helix": -0.203, "PRO-coil/TYR-coil": 0.481, "PRO-coil/TYR-sheet": -0.234, "PRO-coil/MET-helix": -2.404, "PRO-coil/MET-coil": -0.416, "PRO-coil/MET-sheet": -1.01, "PRO-sheet/CYS-helix": -Infinity, "PRO-sheet/CYS-coil": -0.776, "PRO-sheet/CYS-sheet": -1.064, "PRO-sheet/ASP-helix": -2.934, "PRO-sheet/ASP-coil": 0.591, "PRO-sheet/ASP-sheet": -Infinity, "PRO-sheet/SER-helix": -1.702, "PRO-sheet/SER-coil": -2.676, "PRO-sheet/SER-sheet": -1.34, "PRO-sheet/GLN-helix": -2.211, "PRO-sheet/GLN-coil": -3.137, "PRO-sheet/GLN-sheet": -0.605, "PRO-sheet/LYS-helix": -2.823, "PRO-sheet/LYS-coil": -2.192, "PRO-sheet/LYS-sheet": -2.577, "PRO-sheet/ILE-helix": -3.18, "PRO-sheet/ILE-coil": -1.949, "PRO-sheet/ILE-sheet": -Infinity, "PRO-sheet/PRO-helix": -Infinity, "PRO-sheet/PRO-coil": -2.092, "PRO-sheet/PRO-sheet": 5.976, "PRO-sheet/THR-helix": -2.989, "PRO-sheet/THR-coil": -2.448, "PRO-sheet/THR-sheet": -3.036, "PRO-sheet/PHE-helix": -3.133, "PRO-sheet/PHE-coil": -0.944, "PRO-sheet/PHE-sheet": -1.686, "PRO-sheet/ASN-helix": -2.015, "PRO-sheet/ASN-coil": -Infinity, "PRO-sheet/ASN-sheet": -Infinity, "PRO-sheet/GLY-helix": -3.437, "PRO-sheet/GLY-coil": -1.33, "PRO-sheet/GLY-sheet": -0.602, "PRO-sheet/HIS-helix": -Infinity, "PRO-sheet/HIS-coil": 0.849, "PRO-sheet/HIS-sheet": -0.566, "PRO-sheet/LEU-helix": -3.514, "PRO-sheet/LEU-coil": -0.4, "PRO-sheet/LEU-sheet": 0.782, "PRO-sheet/ARG-helix": -3.377, "PRO-sheet/ARG-coil": -2.787, "PRO-sheet/ARG-sheet": -1.529, "PRO-sheet/TRP-helix": -2.5, "PRO-sheet/TRP-coil": -2.276, "PRO-sheet/TRP-sheet": 4.127, "PRO-sheet/ALA-helix": -3.641, "PRO-sheet/ALA-coil": -1.666, "PRO-sheet/ALA-sheet": 3.03, "PRO-sheet/VAL-helix": -2.599, "PRO-sheet/VAL-coil": -2.75, "PRO-sheet/VAL-sheet": 0.63, "PRO-sheet/GLU-helix": -2.477, "PRO-sheet/GLU-coil": -Infinity, "PRO-sheet/GLU-sheet": 2.888, "PRO-sheet/TYR-helix": -2.259, "PRO-sheet/TYR-coil": -2.008, "PRO-sheet/TYR-sheet": 1.339, "PRO-sheet/MET-helix": -1.928, "PRO-sheet/MET-coil": -1.039, "PRO-sheet/MET-sheet": -1.435, "THR-helix/CYS-helix": 0.228, "THR-helix/CYS-coil": 0.127, "THR-helix/CYS-sheet": -2.412, "THR-helix/ASP-helix": 0.416, "THR-helix/ASP-coil": -2.114, "THR-helix/ASP-sheet": -4.154, "THR-helix/SER-helix": 0.308, "THR-helix/SER-coil": -2.032, "THR-helix/SER-sheet": -4.362, "THR-helix/GLN-helix": 0.025, "THR-helix/GLN-coil": -0.609, "THR-helix/GLN-sheet": -4.255, "THR-helix/LYS-helix": -0.538, "THR-helix/LYS-coil": -1.171, "THR-helix/LYS-sheet": -5.311, "THR-helix/ILE-helix": 0.565, "THR-helix/ILE-coil": -1.519, "THR-helix/ILE-sheet": 0.009, "THR-helix/PRO-helix": 0.58, "THR-helix/PRO-coil": -2.503, "THR-helix/PRO-sheet": -2.989, "THR-helix/THR-helix": 3.509, "THR-helix/THR-coil": -2.567, "THR-helix/THR-sheet": -3.372, "THR-helix/PHE-helix": 1.005, "THR-helix/PHE-coil": -0.42, "THR-helix/PHE-sheet": -Infinity, "THR-helix/ASN-helix": 0.143, "THR-helix/ASN-coil": -1.151, "THR-helix/ASN-sheet": -3.517, "THR-helix/GLY-helix": 1.014, "THR-helix/GLY-coil": -1.509, "THR-helix/GLY-sheet": -3.806, "THR-helix/HIS-helix": -0.01, "THR-helix/HIS-coil": -2.984, "THR-helix/HIS-sheet": -3.993, "THR-helix/LEU-helix": 0.502, "THR-helix/LEU-coil": -0.979, "THR-helix/LEU-sheet": -0.723, "THR-helix/ARG-helix": 0.129, "THR-helix/ARG-coil": -1.76, "THR-helix/ARG-sheet": -2.391, "THR-helix/TRP-helix": -0.104, "THR-helix/TRP-coil": -1.139, "THR-helix/TRP-sheet": -3.918, "THR-helix/ALA-helix": 0.319, "THR-helix/ALA-coil": -0.005, "THR-helix/ALA-sheet": -3.222, "THR-helix/VAL-helix": 1.041, "THR-helix/VAL-coil": -0.018, "THR-helix/VAL-sheet": -3.41, "THR-helix/GLU-helix": 0.21, "THR-helix/GLU-coil": -1.908, "THR-helix/GLU-sheet": -4.401, "THR-helix/TYR-helix": 1.401, "THR-helix/TYR-coil": -0.81, "THR-helix/TYR-sheet": -5.365, "THR-helix/MET-helix": -0.399, "THR-helix/MET-coil": -1.398, "THR-helix/MET-sheet": -3.071, "THR-coil/CYS-helix": -3.875, "THR-coil/CYS-coil": 0.48, "THR-coil/CYS-sheet": -1.648, "THR-coil/ASP-helix": -0.2, "THR-coil/ASP-coil": 1.35, "THR-coil/ASP-sheet": -0.939, "THR-coil/SER-helix": -1.323, "THR-coil/SER-coil": 0.351, "THR-coil/SER-sheet": 0.091, "THR-coil/GLN-helix": 0.118, "THR-coil/GLN-coil": 0.312, "THR-coil/GLN-sheet": 0.302, "THR-coil/LYS-helix": -1.988, "THR-coil/LYS-coil": 0.585, "THR-coil/LYS-sheet": 0.668, "THR-coil/ILE-helix": -1.335, "THR-coil/ILE-coil": 0.483, "THR-coil/ILE-sheet": -0.584, "THR-coil/PRO-helix": 0.442, "THR-coil/PRO-coil": 0.75, "THR-coil/PRO-sheet": -2.448, "THR-coil/THR-helix": -2.567, "THR-coil/THR-coil": 3.417, "THR-coil/THR-sheet": -0.036, "THR-coil/PHE-helix": -0.925, "THR-coil/PHE-coil": 1.205, "THR-coil/PHE-sheet": 0.657, "THR-coil/ASN-helix": -1.723, "THR-coil/ASN-coil": 0.521, "THR-coil/ASN-sheet": 1.734, "THR-coil/GLY-helix": -1.632, "THR-coil/GLY-coil": 1.383, "THR-coil/GLY-sheet": 0.373, "THR-coil/HIS-helix": 0.364, "THR-coil/HIS-coil": -0.948, "THR-coil/HIS-sheet": 0.355, "THR-coil/LEU-helix": -1.423, "THR-coil/LEU-coil": 0.255, "THR-coil/LEU-sheet": 0.024, "THR-coil/ARG-helix": -2.152, "THR-coil/ARG-coil": 0.238, "THR-coil/ARG-sheet": -0.033, "THR-coil/TRP-helix": -2.988, "THR-coil/TRP-coil": 0.326, "THR-coil/TRP-sheet": -0.891, "THR-coil/ALA-helix": -0.61, "THR-coil/ALA-coil": 1.443, "THR-coil/ALA-sheet": -1.725, "THR-coil/VAL-helix": -1.436, "THR-coil/VAL-coil": 0.782, "THR-coil/VAL-sheet": -2.726, "THR-coil/GLU-helix": -0.509, "THR-coil/GLU-coil": 0.534, "THR-coil/GLU-sheet": 1.128, "THR-coil/TYR-helix": -1.166, "THR-coil/TYR-coil": 0.887, "THR-coil/TYR-sheet": -0.079, "THR-coil/MET-helix": -3.715, "THR-coil/MET-coil": 0.846, "THR-coil/MET-sheet": -2.935, "THR-sheet/CYS-helix": -Infinity, "THR-sheet/CYS-coil": -1.766, "THR-sheet/CYS-sheet": 2.562, "THR-sheet/ASP-helix": -4.211, "THR-sheet/ASP-coil": -0.09, "THR-sheet/ASP-sheet": 1.773, "THR-sheet/SER-helix": -3.298, "THR-sheet/SER-coil": 0.122, "THR-sheet/SER-sheet": 2.748, "THR-sheet/GLN-helix": -4.181, "THR-sheet/GLN-coil": -3.433, "THR-sheet/GLN-sheet": 1.958, "THR-sheet/LYS-helix": -3.813, "THR-sheet/LYS-coil": -0.957, "THR-sheet/LYS-sheet": 2.662, "THR-sheet/ILE-helix": -4.015, "THR-sheet/ILE-coil": -0.838, "THR-sheet/ILE-sheet": 1.529, "THR-sheet/PRO-helix": -4.116, "THR-sheet/PRO-coil": 0.572, "THR-sheet/PRO-sheet": -3.036, "THR-sheet/THR-helix": -3.372, "THR-sheet/THR-coil": -0.036, "THR-sheet/THR-sheet": 4.686, "THR-sheet/PHE-helix": -5.509, "THR-sheet/PHE-coil": -2.185, "THR-sheet/PHE-sheet": 2.387, "THR-sheet/ASN-helix": -3.985, "THR-sheet/ASN-coil": 0.478, "THR-sheet/ASN-sheet": 1.938, "THR-sheet/GLY-helix": -0.24, "THR-sheet/GLY-coil": 0.756, "THR-sheet/GLY-sheet": 2.49, "THR-sheet/HIS-helix": -Infinity, "THR-sheet/HIS-coil": -3.804, "THR-sheet/HIS-sheet": 0.95, "THR-sheet/LEU-helix": -1.136, "THR-sheet/LEU-coil": 0.909, "THR-sheet/LEU-sheet": 2.435, "THR-sheet/ARG-helix": -0.345, "THR-sheet/ARG-coil": -0.577, "THR-sheet/ARG-sheet": 1.404, "THR-sheet/TRP-helix": -Infinity, "THR-sheet/TRP-coil": 0.318, "THR-sheet/TRP-sheet": 1.734, "THR-sheet/ALA-helix": -3.203, "THR-sheet/ALA-coil": 0.216, "THR-sheet/ALA-sheet": 2.205, "THR-sheet/VAL-helix": -4.127, "THR-sheet/VAL-coil": 0.575, "THR-sheet/VAL-sheet": 1.071, "THR-sheet/GLU-helix": -4.006, "THR-sheet/GLU-coil": -2.814, "THR-sheet/GLU-sheet": 3.187, "THR-sheet/TYR-helix": -1.185, "THR-sheet/TYR-coil": -0.877, "THR-sheet/TYR-sheet": 1.611, "THR-sheet/MET-helix": -4.997, "THR-sheet/MET-coil": -2.498, "THR-sheet/MET-sheet": -0.173, "PHE-helix/CYS-helix": -0.613, "PHE-helix/CYS-coil": -1.377, "PHE-helix/CYS-sheet": -2.15, "PHE-helix/ASP-helix": -1.312, "PHE-helix/ASP-coil": -4.309, "PHE-helix/ASP-sheet": -Infinity, "PHE-helix/SER-helix": -0.058, "PHE-helix/SER-coil": -1.622, "PHE-helix/SER-sheet": -3.813, "PHE-helix/GLN-helix": 0.096, "PHE-helix/GLN-coil": -1.347, "PHE-helix/GLN-sheet": -Infinity, "PHE-helix/LYS-helix": -0.319, "PHE-helix/LYS-coil": -1.694, "PHE-helix/LYS-sheet": -4.356, "PHE-helix/ILE-helix": -0.777, "PHE-helix/ILE-coil": -0.998, "PHE-helix/ILE-sheet": -3.685, "PHE-helix/PRO-helix": 0.07, "PHE-helix/PRO-coil": -1.197, "PHE-helix/PRO-sheet": -3.133, "PHE-helix/THR-helix": 1.005, "PHE-helix/THR-coil": -0.925, "PHE-helix/THR-sheet": -5.509, "PHE-helix/PHE-helix": 3.534, "PHE-helix/PHE-coil": -1.503, "PHE-helix/PHE-sheet": -2.955, "PHE-helix/ASN-helix": -0.681, "PHE-helix/ASN-coil": -3.199, "PHE-helix/ASN-sheet": -Infinity, "PHE-helix/GLY-helix": 0.33, "PHE-helix/GLY-coil": -1.349, "PHE-helix/GLY-sheet": -4.461, "PHE-helix/HIS-helix": -0.574, "PHE-helix/HIS-coil": -2.03, "PHE-helix/HIS-sheet": -Infinity, "PHE-helix/LEU-helix": 0.063, "PHE-helix/LEU-coil": -1.45, "PHE-helix/LEU-sheet": -2.607, "PHE-helix/ARG-helix": -0.449, "PHE-helix/ARG-coil": -2.061, "PHE-helix/ARG-sheet": -2.297, "PHE-helix/TRP-helix": -0.845, "PHE-helix/TRP-coil": 0.351, "PHE-helix/TRP-sheet": -4.349, "PHE-helix/ALA-helix": 0.166, "PHE-helix/ALA-coil": -0.979, "PHE-helix/ALA-sheet": -4.57, "PHE-helix/VAL-helix": 0.146, "PHE-helix/VAL-coil": -2.131, "PHE-helix/VAL-sheet": -3.922, "PHE-helix/GLU-helix": -1.282, "PHE-helix/GLU-coil": -3.787, "PHE-helix/GLU-sheet": -Infinity, "PHE-helix/TYR-helix": 0.345, "PHE-helix/TYR-coil": -2.401, "PHE-helix/TYR-sheet": -3.494, "PHE-helix/MET-helix": 0.679, "PHE-helix/MET-coil": -2.064, "PHE-helix/MET-sheet": -Infinity, "PHE-coil/CYS-helix": -Infinity, "PHE-coil/CYS-coil": -1.466, "PHE-coil/CYS-sheet": -1.466, "PHE-coil/ASP-helix": 1.591, "PHE-coil/ASP-coil": 1.394, "PHE-coil/ASP-sheet": -Infinity, "PHE-coil/SER-helix": 0.273, "PHE-coil/SER-coil": 0.656, "PHE-coil/SER-sheet": -1.875, "PHE-coil/GLN-helix": 0.501, "PHE-coil/GLN-coil": 0.295, "PHE-coil/GLN-sheet": -3.309, "PHE-coil/LYS-helix": -1.115, "PHE-coil/LYS-coil": 0.861, "PHE-coil/LYS-sheet": -0.473, "PHE-coil/ILE-helix": -2.931, "PHE-coil/ILE-coil": 0.417, "PHE-coil/ILE-sheet": -2.212, "PHE-coil/PRO-helix": -1.043, "PHE-coil/PRO-coil": 0.819, "PHE-coil/PRO-sheet": -0.944, "PHE-coil/THR-helix": -0.42, "PHE-coil/THR-coil": 1.205, "PHE-coil/THR-sheet": -2.185, "PHE-coil/PHE-helix": -1.503, "PHE-coil/PHE-coil": 4.02, "PHE-coil/PHE-sheet": -2.781, "PHE-coil/ASN-helix": -1.339, "PHE-coil/ASN-coil": 0.221, "PHE-coil/ASN-sheet": -0.912, "PHE-coil/GLY-helix": -0.671, "PHE-coil/GLY-coil": 1.375, "PHE-coil/GLY-sheet": -2.39, "PHE-coil/HIS-helix": -2.166, "PHE-coil/HIS-coil": -0.137, "PHE-coil/HIS-sheet": -3.047, "PHE-coil/LEU-helix": -0.52, "PHE-coil/LEU-coil": 0.24, "PHE-coil/LEU-sheet": -1.922, "PHE-coil/ARG-helix": -1.022, "PHE-coil/ARG-coil": 0.538, "PHE-coil/ARG-sheet": -1.813, "PHE-coil/TRP-helix": -0.955, "PHE-coil/TRP-coil": -1.12, "PHE-coil/TRP-sheet": -Infinity, "PHE-coil/ALA-helix": -2.743, "PHE-coil/ALA-coil": 0.614, "PHE-coil/ALA-sheet": -1.87, "PHE-coil/VAL-helix": -0.309, "PHE-coil/VAL-coil": 0.833, "PHE-coil/VAL-sheet": -2.39, "PHE-coil/GLU-helix": -1.979, "PHE-coil/GLU-coil": 0.492, "PHE-coil/GLU-sheet": -2.068, "PHE-coil/TYR-helix": 0.076, "PHE-coil/TYR-coil": 1.401, "PHE-coil/TYR-sheet": -1.087, "PHE-coil/MET-helix": -1.77, "PHE-coil/MET-coil": 0.253, "PHE-coil/MET-sheet": -2.53, "PHE-sheet/CYS-helix": -3.113, "PHE-sheet/CYS-coil": -1.515, "PHE-sheet/CYS-sheet": 0.395, "PHE-sheet/ASP-helix": -5.059, "PHE-sheet/ASP-coil": 0.213, "PHE-sheet/ASP-sheet": 2.598, "PHE-sheet/SER-helix": -4.839, "PHE-sheet/SER-coil": -0.196, "PHE-sheet/SER-sheet": 1.722, "PHE-sheet/GLN-helix": -3.642, "PHE-sheet/GLN-coil": -2.266, "PHE-sheet/GLN-sheet": 1.748, "PHE-sheet/LYS-helix": -5.353, "PHE-sheet/LYS-coil": -0.591, "PHE-sheet/LYS-sheet": 1.228, "PHE-sheet/ILE-helix": -3.407, "PHE-sheet/ILE-coil": 0.742, "PHE-sheet/ILE-sheet": 2.093, "PHE-sheet/PRO-helix": -2.949, "PHE-sheet/PRO-coil": -0.008, "PHE-sheet/PRO-sheet": -1.686, "PHE-sheet/THR-helix": -Infinity, "PHE-sheet/THR-coil": 0.657, "PHE-sheet/THR-sheet": 2.387, "PHE-sheet/PHE-helix": -2.955, "PHE-sheet/PHE-coil": -2.781, "PHE-sheet/PHE-sheet": 4.842, "PHE-sheet/ASN-helix": -Infinity, "PHE-sheet/ASN-coil": -0.301, "PHE-sheet/ASN-sheet": 2.052, "PHE-sheet/GLY-helix": -2.617, "PHE-sheet/GLY-coil": 0.253, "PHE-sheet/GLY-sheet": 1.976, "PHE-sheet/HIS-helix": -3.602, "PHE-sheet/HIS-coil": -4.246, "PHE-sheet/HIS-sheet": -0.745, "PHE-sheet/LEU-helix": -0.975, "PHE-sheet/LEU-coil": 0.374, "PHE-sheet/LEU-sheet": 1.937, "PHE-sheet/ARG-helix": -4.808, "PHE-sheet/ARG-coil": -1.019, "PHE-sheet/ARG-sheet": 1.494, "PHE-sheet/TRP-helix": -1.58, "PHE-sheet/TRP-coil": -0.663, "PHE-sheet/TRP-sheet": 1.225, "PHE-sheet/ALA-helix": -2.044, "PHE-sheet/ALA-coil": -0.964, "PHE-sheet/ALA-sheet": 3.08, "PHE-sheet/VAL-helix": -3.424, "PHE-sheet/VAL-coil": -0.156, "PHE-sheet/VAL-sheet": 1.468, "PHE-sheet/GLU-helix": -Infinity, "PHE-sheet/GLU-coil": -0.433, "PHE-sheet/GLU-sheet": 1.542, "PHE-sheet/TYR-helix": -5.077, "PHE-sheet/TYR-coil": -1.973, "PHE-sheet/TYR-sheet": 2.304, "PHE-sheet/MET-helix": 0.538, "PHE-sheet/MET-coil": -3.856, "PHE-sheet/MET-sheet": 1.843, "ASN-helix/CYS-helix": -0.379, "ASN-helix/CYS-coil": -0.051, "ASN-helix/CYS-sheet": -Infinity, "ASN-helix/ASP-helix": 1.15, "ASN-helix/ASP-coil": -1.269, "ASN-helix/ASP-sheet": -3.467, "ASN-helix/SER-helix": 0.677, "ASN-helix/SER-coil": 0.121, "ASN-helix/SER-sheet": -4.081, "ASN-helix/GLN-helix": 0.662, "ASN-helix/GLN-coil": -0.194, "ASN-helix/GLN-sheet": -4.262, "ASN-helix/LYS-helix": 1.5, "ASN-helix/LYS-coil": -0.422, "ASN-helix/LYS-sheet": -1.985, "ASN-helix/ILE-helix": -0.612, "ASN-helix/ILE-coil": -2.53, "ASN-helix/ILE-sheet": -3.077, "ASN-helix/PRO-helix": -0.233, "ASN-helix/PRO-coil": -2.138, "ASN-helix/PRO-sheet": -2.015, "ASN-helix/THR-helix": 0.143, "ASN-helix/THR-coil": -1.723, "ASN-helix/THR-sheet": -3.985, "ASN-helix/PHE-helix": -0.681, "ASN-helix/PHE-coil": -1.339, "ASN-helix/PHE-sheet": -Infinity, "ASN-helix/ASN-helix": 4.036, "ASN-helix/ASN-coil": -1.016, "ASN-helix/ASN-sheet": -2.83, "ASN-helix/GLY-helix": 0.646, "ASN-helix/GLY-coil": -0.245, "ASN-helix/GLY-sheet": -Infinity, "ASN-helix/HIS-helix": -1.445, "ASN-helix/HIS-coil": 0.267, "ASN-helix/HIS-sheet": -2.613, "ASN-helix/LEU-helix": -0.266, "ASN-helix/LEU-coil": -2.853, "ASN-helix/LEU-sheet": -0.917, "ASN-helix/ARG-helix": 0.969, "ASN-helix/ARG-coil": -0.922, "ASN-helix/ARG-sheet": 0.173, "ASN-helix/TRP-helix": -0.871, "ASN-helix/TRP-coil": -1.926, "ASN-helix/TRP-sheet": -3.924, "ASN-helix/ALA-helix": -0.552, "ASN-helix/ALA-coil": -1.999, "ASN-helix/ALA-sheet": -3.046, "ASN-helix/VAL-helix": 0.031, "ASN-helix/VAL-coil": -1.994, "ASN-helix/VAL-sheet": -4.19, "ASN-helix/GLU-helix": 1.001, "ASN-helix/GLU-coil": 0.511, "ASN-helix/GLU-sheet": -Infinity, "ASN-helix/TYR-helix": -0.246, "ASN-helix/TYR-coil": -1.858, "ASN-helix/TYR-sheet": -3.58, "ASN-helix/MET-helix": 0.123, "ASN-helix/MET-coil": -1.833, "ASN-helix/MET-sheet": -Infinity, "ASN-coil/CYS-helix": -1.747, "ASN-coil/CYS-coil": -1.471, "ASN-coil/CYS-sheet": -Infinity, "ASN-coil/ASP-helix": -2.307, "ASN-coil/ASP-coil": 0.247, "ASN-coil/ASP-sheet": -1.016, "ASN-coil/SER-helix": -1.329, "ASN-coil/SER-coil": 1.311, "ASN-coil/SER-sheet": -3.015, "ASN-coil/GLN-helix": -3.039, "ASN-coil/GLN-coil": -0.019, "ASN-coil/GLN-sheet": -1.117, "ASN-coil/LYS-helix": -0.404, "ASN-coil/LYS-coil": 0.793, "ASN-coil/LYS-sheet": -2.578, "ASN-coil/ILE-helix": -2.804, "ASN-coil/ILE-coil": 0.569, "ASN-coil/ILE-sheet": -0.97, "ASN-coil/PRO-helix": -3.821, "ASN-coil/PRO-coil": 0.14, "ASN-coil/PRO-sheet": -Infinity, "ASN-coil/THR-helix": -1.151, "ASN-coil/THR-coil": 0.521, "ASN-coil/THR-sheet": 0.478, "ASN-coil/PHE-helix": -3.199, "ASN-coil/PHE-coil": 0.221, "ASN-coil/PHE-sheet": -0.301, "ASN-coil/ASN-helix": -1.016, "ASN-coil/ASN-coil": 3.66, "ASN-coil/ASN-sheet": 0.512, "ASN-coil/GLY-helix": -3.266, "ASN-coil/GLY-coil": 0.767, "ASN-coil/GLY-sheet": -0.619, "ASN-coil/HIS-helix": 1.429, "ASN-coil/HIS-coil": -0.301, "ASN-coil/HIS-sheet": -2.647, "ASN-coil/LEU-helix": -2.704, "ASN-coil/LEU-coil": 0.578, "ASN-coil/LEU-sheet": -1.837, "ASN-coil/ARG-helix": -0.947, "ASN-coil/ARG-coil": -0.162, "ASN-coil/ARG-sheet": 0.485, "ASN-coil/TRP-helix": -1.402, "ASN-coil/TRP-coil": -1.138, "ASN-coil/TRP-sheet": 0.796, "ASN-coil/ALA-helix": -1.055, "ASN-coil/ALA-coil": 0.219, "ASN-coil/ALA-sheet": -2.099, "ASN-coil/VAL-helix": -2.664, "ASN-coil/VAL-coil": 0.764, "ASN-coil/VAL-sheet": -2.613, "ASN-coil/GLU-helix": 0.313, "ASN-coil/GLU-coil": 0.808, "ASN-coil/GLU-sheet": -1.496, "ASN-coil/TYR-helix": -1.492, "ASN-coil/TYR-coil": -0.267, "ASN-coil/TYR-sheet": -0.927, "ASN-coil/MET-helix": -0.377, "ASN-coil/MET-coil": 0.265, "ASN-coil/MET-sheet": -1.319, "ASN-sheet/CYS-helix": -1.581, "ASN-sheet/CYS-coil": -Infinity, "ASN-sheet/CYS-sheet": -0.781, "ASN-sheet/ASP-helix": -3.75, "ASN-sheet/ASP-coil": 1.012, "ASN-sheet/ASP-sheet": 3.06, "ASN-sheet/SER-helix": -4.223, "ASN-sheet/SER-coil": 0.625, "ASN-sheet/SER-sheet": 3.779, "ASN-sheet/GLN-helix": -Infinity, "ASN-sheet/GLN-coil": 1.003, "ASN-sheet/GLN-sheet": 2.391, "ASN-sheet/LYS-helix": -Infinity, "ASN-sheet/LYS-coil": 0.526, "ASN-sheet/LYS-sheet": 2.184, "ASN-sheet/ILE-helix": -3.995, "ASN-sheet/ILE-coil": 1.141, "ASN-sheet/ILE-sheet": 0.895, "ASN-sheet/PRO-helix": -2.556, "ASN-sheet/PRO-coil": -1.683, "ASN-sheet/PRO-sheet": -Infinity, "ASN-sheet/THR-helix": -3.517, "ASN-sheet/THR-coil": 1.734, "ASN-sheet/THR-sheet": 1.938, "ASN-sheet/PHE-helix": -Infinity, "ASN-sheet/PHE-coil": -0.912, "ASN-sheet/PHE-sheet": 2.052, "ASN-sheet/ASN-helix": -2.83, "ASN-sheet/ASN-coil": 0.512, "ASN-sheet/ASN-sheet": 5.551, "ASN-sheet/GLY-helix": 0.065, "ASN-sheet/GLY-coil": 0.676, "ASN-sheet/GLY-sheet": 3.493, "ASN-sheet/HIS-helix": -Infinity, "ASN-sheet/HIS-coil": -Infinity, "ASN-sheet/HIS-sheet": 2.746, "ASN-sheet/LEU-helix": -3.636, "ASN-sheet/LEU-coil": 0.486, "ASN-sheet/LEU-sheet": 1.293, "ASN-sheet/ARG-helix": -Infinity, "ASN-sheet/ARG-coil": -0.74, "ASN-sheet/ARG-sheet": 2.369, "ASN-sheet/TRP-helix": -4.009, "ASN-sheet/TRP-coil": -1.483, "ASN-sheet/TRP-sheet": 1.612, "ASN-sheet/ALA-helix": -3.609, "ASN-sheet/ALA-coil": 1.29, "ASN-sheet/ALA-sheet": 2.439, "ASN-sheet/VAL-helix": -4.108, "ASN-sheet/VAL-coil": 0.316, "ASN-sheet/VAL-sheet": 2.168, "ASN-sheet/GLU-helix": -3.986, "ASN-sheet/GLU-coil": -0.695, "ASN-sheet/GLU-sheet": 2.108, "ASN-sheet/TYR-helix": -3.768, "ASN-sheet/TYR-coil": 0.976, "ASN-sheet/TYR-sheet": 1.008, "ASN-sheet/MET-helix": -3.031, "ASN-sheet/MET-coil": -2.547, "ASN-sheet/MET-sheet": 1.7, "GLY-helix/CYS-helix": 1.417, "GLY-helix/CYS-coil": 0.208, "GLY-helix/CYS-sheet": -3.147, "GLY-helix/ASP-helix": -0.831, "GLY-helix/ASP-coil": -2.794, "GLY-helix/ASP-sheet": -4.196, "GLY-helix/SER-helix": 1.352, "GLY-helix/SER-coil": -1.364, "GLY-helix/SER-sheet": -4.404, "GLY-helix/GLN-helix": -1.453, "GLY-helix/GLN-coil": -0.21, "GLY-helix/GLN-sheet": 0.127, "GLY-helix/LYS-helix": -0.756, "GLY-helix/LYS-coil": -2.799, "GLY-helix/LYS-sheet": -3.156, "GLY-helix/ILE-helix": 1.022, "GLY-helix/ILE-coil": 1.04, "GLY-helix/ILE-sheet": -3.519, "GLY-helix/PRO-helix": 0.721, "GLY-helix/PRO-coil": -1.823, "GLY-helix/PRO-sheet": -3.437, "GLY-helix/THR-helix": 1.014, "GLY-helix/THR-coil": -1.632, "GLY-helix/THR-sheet": -0.24, "GLY-helix/PHE-helix": 0.33, "GLY-helix/PHE-coil": -0.671, "GLY-helix/PHE-sheet": -2.617, "GLY-helix/ASN-helix": 0.646, "GLY-helix/ASN-coil": -3.266, "GLY-helix/ASN-sheet": 0.065, "GLY-helix/GLY-helix": 3.531, "GLY-helix/GLY-coil": -0.48, "GLY-helix/GLY-sheet": -3.848, "GLY-helix/HIS-helix": -0.577, "GLY-helix/HIS-coil": -2.19, "GLY-helix/HIS-sheet": -Infinity, "GLY-helix/LEU-helix": 0.928, "GLY-helix/LEU-coil": 0.266, "GLY-helix/LEU-sheet": -4.479, "GLY-helix/ARG-helix": 0.194, "GLY-helix/ARG-coil": -0.822, "GLY-helix/ARG-sheet": -3.053, "GLY-helix/TRP-helix": -0.416, "GLY-helix/TRP-coil": -0.224, "GLY-helix/TRP-sheet": -3.96, "GLY-helix/ALA-helix": 1.817, "GLY-helix/ALA-coil": -0.205, "GLY-helix/ALA-sheet": -4.181, "GLY-helix/VAL-helix": 1.349, "GLY-helix/VAL-coil": -1.163, "GLY-helix/VAL-sheet": -3.619, "GLY-helix/GLU-helix": -0.152, "GLY-helix/GLU-coil": -1.839, "GLY-helix/GLU-sheet": -4.443, "GLY-helix/TYR-helix": 0.069, "GLY-helix/TYR-coil": -0.866, "GLY-helix/TYR-sheet": -3.798, "GLY-helix/MET-helix": -0.268, "GLY-helix/MET-coil": -1.982, "GLY-helix/MET-sheet": -Infinity, "GLY-coil/CYS-helix": -0.204, "GLY-coil/CYS-coil": 0.758, "GLY-coil/CYS-sheet": -1.105, "GLY-coil/ASP-helix": -2.033, "GLY-coil/ASP-coil": 1.351, "GLY-coil/ASP-sheet": 1.075, "GLY-coil/SER-helix": -0.459, "GLY-coil/SER-coil": 1.466, "GLY-coil/SER-sheet": 0.62, "GLY-coil/GLN-helix": -1.751, "GLY-coil/GLN-coil": 0.828, "GLY-coil/GLN-sheet": 0.02, "GLY-coil/LYS-helix": -1.15, "GLY-coil/LYS-coil": 0.582, "GLY-coil/LYS-sheet": 1.18, "GLY-coil/ILE-helix": -2.484, "GLY-coil/ILE-coil": 1.342, "GLY-coil/ILE-sheet": -0.608, "GLY-coil/PRO-helix": -0.186, "GLY-coil/PRO-coil": 1.231, "GLY-coil/PRO-sheet": -1.33, "GLY-coil/THR-helix": -1.509, "GLY-coil/THR-coil": 1.383, "GLY-coil/THR-sheet": 0.756, "GLY-coil/PHE-helix": -1.349, "GLY-coil/PHE-coil": 1.375, "GLY-coil/PHE-sheet": 0.253, "GLY-coil/ASN-helix": -0.245, "GLY-coil/ASN-coil": 0.767, "GLY-coil/ASN-sheet": 0.676, "GLY-coil/GLY-helix": -0.48, "GLY-coil/GLY-coil": 3.058, "GLY-coil/GLY-sheet": -0.226, "GLY-coil/HIS-helix": -1.502, "GLY-coil/HIS-coil": 0.363, "GLY-coil/HIS-sheet": 0.017, "GLY-coil/LEU-helix": -0.501, "GLY-coil/LEU-coil": 0.753, "GLY-coil/LEU-sheet": 0.781, "GLY-coil/ARG-helix": -1.089, "GLY-coil/ARG-coil": -0.049, "GLY-coil/ARG-sheet": -2.081, "GLY-coil/TRP-helix": -1.03, "GLY-coil/TRP-coil": -0.12, "GLY-coil/TRP-sheet": 0.962, "GLY-coil/ALA-helix": -0.855, "GLY-coil/ALA-coil": 1.156, "GLY-coil/ALA-sheet": -0.106, "GLY-coil/VAL-helix": -0.931, "GLY-coil/VAL-coil": 1.426, "GLY-coil/VAL-sheet": 0.412, "GLY-coil/GLU-helix": -1.593, "GLY-coil/GLU-coil": 0.488, "GLY-coil/GLU-sheet": 0.112, "GLY-coil/TYR-helix": -1.418, "GLY-coil/TYR-coil": 0.451, "GLY-coil/TYR-sheet": 0.255, "GLY-coil/MET-helix": -1.438, "GLY-coil/MET-coil": -0.403, "GLY-coil/MET-sheet": -1.7, "GLY-sheet/CYS-helix": -2.604, "GLY-sheet/CYS-coil": -2.104, "GLY-sheet/CYS-sheet": 3.072, "GLY-sheet/ASP-helix": -3.569, "GLY-sheet/ASP-coil": -0.931, "GLY-sheet/ASP-sheet": 2.908, "GLY-sheet/SER-helix": -4.042, "GLY-sheet/SER-coil": -1.844, "GLY-sheet/SER-sheet": 3.637, "GLY-sheet/GLN-helix": -3.826, "GLY-sheet/GLN-coil": -2.386, "GLY-sheet/GLN-sheet": 2.973, "GLY-sheet/LYS-helix": -Infinity, "GLY-sheet/LYS-coil": -1.192, "GLY-sheet/LYS-sheet": 2.697, "GLY-sheet/ILE-helix": -3.814, "GLY-sheet/ILE-coil": -0.662, "GLY-sheet/ILE-sheet": 2.825, "GLY-sheet/PRO-helix": -2.845, "GLY-sheet/PRO-coil": -1.391, "GLY-sheet/PRO-sheet": -0.602, "GLY-sheet/THR-helix": -3.806, "GLY-sheet/THR-coil": 0.373, "GLY-sheet/THR-sheet": 2.49, "GLY-sheet/PHE-helix": -4.461, "GLY-sheet/PHE-coil": -2.39, "GLY-sheet/PHE-sheet": 1.976, "GLY-sheet/ASN-helix": -Infinity, "GLY-sheet/ASN-coil": -0.619, "GLY-sheet/ASN-sheet": 3.493, "GLY-sheet/GLY-helix": -3.848, "GLY-sheet/GLY-coil": -0.226, "GLY-sheet/GLY-sheet": 5.133, "GLY-sheet/HIS-helix": -Infinity, "GLY-sheet/HIS-coil": -2.756, "GLY-sheet/HIS-sheet": 2.469, "GLY-sheet/LEU-helix": -3.589, "GLY-sheet/LEU-coil": 0.802, "GLY-sheet/LEU-sheet": 3.442, "GLY-sheet/ARG-helix": -4.705, "GLY-sheet/ARG-coil": -1.17, "GLY-sheet/ARG-sheet": 0.699, "GLY-sheet/TRP-helix": -3.828, "GLY-sheet/TRP-coil": -1.995, "GLY-sheet/TRP-sheet": 1.558, "GLY-sheet/ALA-helix": -3.988, "GLY-sheet/ALA-coil": -0.748, "GLY-sheet/ALA-sheet": 3.444, "GLY-sheet/VAL-helix": -4.109, "GLY-sheet/VAL-coil": -0.566, "GLY-sheet/VAL-sheet": 3.152, "GLY-sheet/GLU-helix": -4.21, "GLY-sheet/GLU-coil": -0.937, "GLY-sheet/GLU-sheet": 1.733, "GLY-sheet/TYR-helix": -3.587, "GLY-sheet/TYR-coil": 1.758, "GLY-sheet/TYR-sheet": 2.809, "GLY-sheet/MET-helix": -2.445, "GLY-sheet/MET-coil": -2.143, "GLY-sheet/MET-sheet": 1.91, "HIS-helix/CYS-helix": 0.193, "HIS-helix/CYS-coil": -0.782, "HIS-helix/CYS-sheet": -Infinity, "HIS-helix/ASP-helix": 0.42, "HIS-helix/ASP-coil": -1.737, "HIS-helix/ASP-sheet": -Infinity, "HIS-helix/SER-helix": 0.958, "HIS-helix/SER-coil": -2.577, "HIS-helix/SER-sheet": -4.236, "HIS-helix/GLN-helix": 1.5, "HIS-helix/GLN-coil": 0.044, "HIS-helix/GLN-sheet": -3.724, "HIS-helix/LYS-helix": 1.001, "HIS-helix/LYS-coil": -1.156, "HIS-helix/LYS-sheet": -4.087, "HIS-helix/ILE-helix": -1.694, "HIS-helix/ILE-coil": -1.613, "HIS-helix/ILE-sheet": -4.331, "HIS-helix/PRO-helix": -0.23, "HIS-helix/PRO-coil": -2.646, "HIS-helix/PRO-sheet": -Infinity, "HIS-helix/THR-helix": -0.01, "HIS-helix/THR-coil": 0.364, "HIS-helix/THR-sheet": -Infinity, "HIS-helix/PHE-helix": -0.574, "HIS-helix/PHE-coil": -2.166, "HIS-helix/PHE-sheet": -3.602, "HIS-helix/ASN-helix": -1.445, "HIS-helix/ASN-coil": 1.429, "HIS-helix/ASN-sheet": -Infinity, "HIS-helix/GLY-helix": -0.577, "HIS-helix/GLY-coil": -1.502, "HIS-helix/GLY-sheet": -Infinity, "HIS-helix/HIS-helix": 4.632, "HIS-helix/HIS-coil": -0.636, "HIS-helix/HIS-sheet": -Infinity, "HIS-helix/LEU-helix": 0.752, "HIS-helix/LEU-coil": -0.076, "HIS-helix/LEU-sheet": -3.03, "HIS-helix/ARG-helix": 0.368, "HIS-helix/ARG-coil": -1.078, "HIS-helix/ARG-sheet": -3.732, "HIS-helix/TRP-helix": -0.469, "HIS-helix/TRP-coil": -Infinity, "HIS-helix/TRP-sheet": -Infinity, "HIS-helix/ALA-helix": -0.438, "HIS-helix/ALA-coil": -0.169, "HIS-helix/ALA-sheet": -3.202, "HIS-helix/VAL-helix": -0.313, "HIS-helix/VAL-coil": -1.315, "HIS-helix/VAL-sheet": 1.236, "HIS-helix/GLU-helix": 1.112, "HIS-helix/GLU-coil": -0.968, "HIS-helix/GLU-sheet": -Infinity, "HIS-helix/TYR-helix": -0.165, "HIS-helix/TYR-coil": -0.81, "HIS-helix/TYR-sheet": -Infinity, "HIS-helix/MET-helix": 0.294, "HIS-helix/MET-coil": -1.989, "HIS-helix/MET-sheet": -Infinity, "HIS-coil/CYS-helix": -1.534, "HIS-coil/CYS-coil": 0.204, "HIS-coil/CYS-sheet": -Infinity, "HIS-coil/ASP-helix": -1.592, "HIS-coil/ASP-coil": -0.072, "HIS-coil/ASP-sheet": 0.377, "HIS-coil/SER-helix": 0.237, "HIS-coil/SER-coil": 0.011, "HIS-coil/SER-sheet": -2.579, "HIS-coil/GLN-helix": -1.763, "HIS-coil/GLN-coil": -0.61, "HIS-coil/GLN-sheet": 0.731, "HIS-coil/LYS-helix": -0.799, "HIS-coil/LYS-coil": -0.465, "HIS-coil/LYS-sheet": 0.484, "HIS-coil/ILE-helix": -1.145, "HIS-coil/ILE-coil": 0.461, "HIS-coil/ILE-sheet": -3.185, "HIS-coil/PRO-helix": -2.509, "HIS-coil/PRO-coil": -0.063, "HIS-coil/PRO-sheet": 0.849, "HIS-coil/THR-helix": -2.984, "HIS-coil/THR-coil": -0.948, "HIS-coil/THR-sheet": -3.804, "HIS-coil/PHE-helix": -2.03, "HIS-coil/PHE-coil": -0.137, "HIS-coil/PHE-sheet": -4.246, "HIS-coil/ASN-helix": 0.267, "HIS-coil/ASN-coil": -0.301, "HIS-coil/ASN-sheet": -Infinity, "HIS-coil/GLY-helix": -2.19, "HIS-coil/GLY-coil": 0.363, "HIS-coil/GLY-sheet": -2.756, "HIS-coil/HIS-helix": -0.636, "HIS-coil/HIS-coil": 3.917, "HIS-coil/HIS-sheet": -2.721, "HIS-coil/LEU-helix": -0.429, "HIS-coil/LEU-coil": -0.16, "HIS-coil/LEU-sheet": 0.124, "HIS-coil/ARG-helix": -1.561, "HIS-coil/ARG-coil": 0.121, "HIS-coil/ARG-sheet": 0.961, "HIS-coil/TRP-helix": -1.71, "HIS-coil/TRP-coil": -0.183, "HIS-coil/TRP-sheet": -1.729, "HIS-coil/ALA-helix": -1.606, "HIS-coil/ALA-coil": 0.564, "HIS-coil/ALA-sheet": -1.256, "HIS-coil/VAL-helix": -1.94, "HIS-coil/VAL-coil": -1.378, "HIS-coil/VAL-sheet": -Infinity, "HIS-coil/GLU-helix": -0.458, "HIS-coil/GLU-coil": 1.154, "HIS-coil/GLU-sheet": -3.128, "HIS-coil/TYR-helix": -0.564, "HIS-coil/TYR-coil": 0.086, "HIS-coil/TYR-sheet": 0.711, "HIS-coil/MET-helix": -3.166, "HIS-coil/MET-coil": -0.591, "HIS-coil/MET-sheet": -Infinity, "HIS-sheet/CYS-helix": -Infinity, "HIS-sheet/CYS-coil": 0.011, "HIS-sheet/CYS-sheet": 0.129, "HIS-sheet/ASP-helix": -2.147, "HIS-sheet/ASP-coil": 1.081, "HIS-sheet/ASP-sheet": 0.061, "HIS-sheet/SER-helix": -3.313, "HIS-sheet/SER-coil": -1.195, "HIS-sheet/SER-sheet": 2.338, "HIS-sheet/GLN-helix": -Infinity, "HIS-sheet/GLN-coil": -1.944, "HIS-sheet/GLN-sheet": 1.569, "HIS-sheet/LYS-helix": -3.135, "HIS-sheet/LYS-coil": -1.404, "HIS-sheet/LYS-sheet": 0.226, "HIS-sheet/ILE-helix": -Infinity, "HIS-sheet/ILE-coil": -2.26, "HIS-sheet/ILE-sheet": 0.269, "HIS-sheet/PRO-helix": -Infinity, "HIS-sheet/PRO-coil": -1.815, "HIS-sheet/PRO-sheet": -0.566, "HIS-sheet/THR-helix": -3.993, "HIS-sheet/THR-coil": 0.355, "HIS-sheet/THR-sheet": 0.95, "HIS-sheet/PHE-helix": -Infinity, "HIS-sheet/PHE-coil": -3.047, "HIS-sheet/PHE-sheet": -0.745, "HIS-sheet/ASN-helix": -2.613, "HIS-sheet/ASN-coil": -2.647, "HIS-sheet/ASN-sheet": 2.746, "HIS-sheet/GLY-helix": -Infinity, "HIS-sheet/GLY-coil": 0.017, "HIS-sheet/GLY-sheet": 2.469, "HIS-sheet/HIS-helix": -Infinity, "HIS-sheet/HIS-coil": -2.721, "HIS-sheet/HIS-sheet": 5.378, "HIS-sheet/LEU-helix": -3.419, "HIS-sheet/LEU-coil": 0.876, "HIS-sheet/LEU-sheet": 1.213, "HIS-sheet/ARG-helix": -3.282, "HIS-sheet/ARG-coil": -2.692, "HIS-sheet/ARG-sheet": -0.182, "HIS-sheet/TRP-helix": -3.099, "HIS-sheet/TRP-coil": -1.489, "HIS-sheet/TRP-sheet": -0.684, "HIS-sheet/ALA-helix": -3.036, "HIS-sheet/ALA-coil": -2.383, "HIS-sheet/ALA-sheet": 1.347, "HIS-sheet/VAL-helix": -3.198, "HIS-sheet/VAL-coil": -2.25, "HIS-sheet/VAL-sheet": -0.581, "HIS-sheet/GLU-helix": -4.175, "HIS-sheet/GLU-coil": -1.577, "HIS-sheet/GLU-sheet": -0.473, "HIS-sheet/TYR-helix": -2.453, "HIS-sheet/TYR-coil": -Infinity, "HIS-sheet/TYR-sheet": 1.452, "HIS-sheet/MET-helix": -2.121, "HIS-sheet/MET-coil": -1.232, "HIS-sheet/MET-sheet": 0.045, "LEU-helix/CYS-helix": 0.602, "LEU-helix/CYS-coil": 0.597, "LEU-helix/CYS-sheet": -2.126, "LEU-helix/ASP-helix": -0.288, "LEU-helix/ASP-coil": -0.833, "LEU-helix/ASP-sheet": -Infinity, "LEU-helix/SER-helix": 0.392, "LEU-helix/SER-coil": -0.708, "LEU-helix/SER-sheet": -5.58, "LEU-helix/GLN-helix": 0.479, "LEU-helix/GLN-coil": -0.548, "LEU-helix/GLN-sheet": -4.375, "LEU-helix/LYS-helix": -0.006, "LEU-helix/LYS-coil": -1.972, "LEU-helix/LYS-sheet": -3.558, "LEU-helix/ILE-helix": 0.42, "LEU-helix/ILE-coil": -0.753, "LEU-helix/ILE-sheet": -0.252, "LEU-helix/PRO-helix": 0.275, "LEU-helix/PRO-coil": -1.246, "LEU-helix/PRO-sheet": -3.514, "LEU-helix/THR-helix": 0.502, "LEU-helix/THR-coil": -1.423, "LEU-helix/THR-sheet": -1.136, "LEU-helix/PHE-helix": 0.063, "LEU-helix/PHE-coil": -0.52, "LEU-helix/PHE-sheet": -0.975, "LEU-helix/ASN-helix": -0.266, "LEU-helix/ASN-coil": -2.704, "LEU-helix/ASN-sheet": -3.636, "LEU-helix/GLY-helix": 0.928, "LEU-helix/GLY-coil": -0.501, "LEU-helix/GLY-sheet": -3.589, "LEU-helix/HIS-helix": 0.752, "LEU-helix/HIS-coil": -0.429, "LEU-helix/HIS-sheet": -3.419, "LEU-helix/LEU-helix": 2.803, "LEU-helix/LEU-coil": -0.828, "LEU-helix/LEU-sheet": -0.638, "LEU-helix/ARG-helix": 0.178, "LEU-helix/ARG-coil": -1.071, "LEU-helix/ARG-sheet": -5.076, "LEU-helix/TRP-helix": -0.184, "LEU-helix/TRP-coil": -0.29, "LEU-helix/TRP-sheet": 0.291, "LEU-helix/ALA-helix": 0.68, "LEU-helix/ALA-coil": -0.214, "LEU-helix/ALA-sheet": -0.244, "LEU-helix/VAL-helix": 0.738, "LEU-helix/VAL-coil": -0.555, "LEU-helix/VAL-sheet": -0.882, "LEU-helix/GLU-helix": 0.294, "LEU-helix/GLU-coil": -2.058, "LEU-helix/GLU-sheet": -2.014, "LEU-helix/TYR-helix": 0.276, "LEU-helix/TYR-coil": -3.475, "LEU-helix/TYR-sheet": -2.622, "LEU-helix/MET-helix": 0.291, "LEU-helix/MET-coil": -0.573, "LEU-helix/MET-sheet": -2.209, "LEU-coil/CYS-helix": -2.52, "LEU-coil/CYS-coil": -1.733, "LEU-coil/CYS-sheet": -1.209, "LEU-coil/ASP-helix": -0.846, "LEU-coil/ASP-coil": 0.589, "LEU-coil/ASP-sheet": 0.15, "LEU-coil/SER-helix": -0.593, "LEU-coil/SER-coil": 1.039, "LEU-coil/SER-sheet": 0.594, "LEU-coil/GLN-helix": 0.032, "LEU-coil/GLN-coil": 0.778, "LEU-coil/GLN-sheet": 0.028, "LEU-coil/LYS-helix": -1.028, "LEU-coil/LYS-coil": 0.686, "LEU-coil/LYS-sheet": 0.131, "LEU-coil/ILE-helix": -2.19, "LEU-coil/ILE-coil": 0.159, "LEU-coil/ILE-sheet": -0.064, "LEU-coil/PRO-helix": -1.392, "LEU-coil/PRO-coil": 1.241, "LEU-coil/PRO-sheet": -0.4, "LEU-coil/THR-helix": -0.979, "LEU-coil/THR-coil": 0.255, "LEU-coil/THR-sheet": 0.909, "LEU-coil/PHE-helix": -1.45, "LEU-coil/PHE-coil": 0.24, "LEU-coil/PHE-sheet": 0.374, "LEU-coil/ASN-helix": -2.853, "LEU-coil/ASN-coil": 0.578, "LEU-coil/ASN-sheet": 0.486, "LEU-coil/GLY-helix": 0.266, "LEU-coil/GLY-coil": 0.753, "LEU-coil/GLY-sheet": 0.802, "LEU-coil/HIS-helix": -0.076, "LEU-coil/HIS-coil": -0.16, "LEU-coil/HIS-sheet": 0.876, "LEU-coil/LEU-helix": -0.828, "LEU-coil/LEU-coil": 3.45, "LEU-coil/LEU-sheet": 0.675, "LEU-coil/ARG-helix": -0.502, "LEU-coil/ARG-coil": -0.309, "LEU-coil/ARG-sheet": -0.694, "LEU-coil/TRP-helix": -0.337, "LEU-coil/TRP-coil": -1.4, "LEU-coil/TRP-sheet": 0.601, "LEU-coil/ALA-helix": -0.287, "LEU-coil/ALA-coil": 1.293, "LEU-coil/ALA-sheet": 0.388, "LEU-coil/VAL-helix": -1.551, "LEU-coil/VAL-coil": 0.649, "LEU-coil/VAL-sheet": 0.942, "LEU-coil/GLU-helix": -1.958, "LEU-coil/GLU-coil": 0.831, "LEU-coil/GLU-sheet": 0.367, "LEU-coil/TYR-helix": -1.911, "LEU-coil/TYR-coil": 0.275, "LEU-coil/TYR-sheet": 1.402, "LEU-coil/MET-helix": -0.394, "LEU-coil/MET-coil": 0.16, "LEU-coil/MET-sheet": -2.168, "LEU-sheet/CYS-helix": -2.947, "LEU-sheet/CYS-coil": -1.349, "LEU-sheet/CYS-sheet": 1.152, "LEU-sheet/ASP-helix": -4.2, "LEU-sheet/ASP-coil": -0.842, "LEU-sheet/ASP-sheet": 1.28, "LEU-sheet/SER-helix": -4.113, "LEU-sheet/SER-coil": 0.24, "LEU-sheet/SER-sheet": 2.486, "LEU-sheet/GLN-helix": -0.319, "LEU-sheet/GLN-coil": -2.323, "LEU-sheet/GLN-sheet": 0.462, "LEU-sheet/LYS-helix": -5.881, "LEU-sheet/LYS-coil": -1.274, "LEU-sheet/LYS-sheet": 1.975, "LEU-sheet/ILE-helix": -4.157, "LEU-sheet/ILE-coil": 0.913, "LEU-sheet/ILE-sheet": 1.775, "LEU-sheet/PRO-helix": -1.589, "LEU-sheet/PRO-coil": -0.74, "LEU-sheet/PRO-sheet": 0.782, "LEU-sheet/THR-helix": -0.723, "LEU-sheet/THR-coil": 0.024, "LEU-sheet/THR-sheet": 2.435, "LEU-sheet/PHE-helix": -2.607, "LEU-sheet/PHE-coil": -1.922, "LEU-sheet/PHE-sheet": 1.937, "LEU-sheet/ASN-helix": -0.917, "LEU-sheet/ASN-coil": -1.837, "LEU-sheet/ASN-sheet": 1.293, "LEU-sheet/GLY-helix": -4.479, "LEU-sheet/GLY-coil": 0.781, "LEU-sheet/GLY-sheet": 3.442, "LEU-sheet/HIS-helix": -3.03, "LEU-sheet/HIS-coil": 0.124, "LEU-sheet/HIS-sheet": 1.213, "LEU-sheet/LEU-helix": -0.638, "LEU-sheet/LEU-coil": 0.675, "LEU-sheet/LEU-sheet": 4.335, "LEU-sheet/ARG-helix": -5.336, "LEU-sheet/ARG-coil": -1.219, "LEU-sheet/ARG-sheet": 1.112, "LEU-sheet/TRP-helix": 0.341, "LEU-sheet/TRP-coil": 0.172, "LEU-sheet/TRP-sheet": 1.427, "LEU-sheet/ALA-helix": -1.971, "LEU-sheet/ALA-coil": 0.919, "LEU-sheet/ALA-sheet": 2.447, "LEU-sheet/VAL-helix": -1.459, "LEU-sheet/VAL-coil": -1.125, "LEU-sheet/VAL-sheet": 2.189, "LEU-sheet/GLU-helix": -4.282, "LEU-sheet/GLU-coil": -1.965, "LEU-sheet/GLU-sheet": 1.962, "LEU-sheet/TYR-helix": -1.356, "LEU-sheet/TYR-coil": -0.775, "LEU-sheet/TYR-sheet": 1.978, "LEU-sheet/MET-helix": -3.481, "LEU-sheet/MET-coil": -2.081, "LEU-sheet/MET-sheet": 0.072, "ARG-helix/CYS-helix": -1.378, "ARG-helix/CYS-coil": -2.24, "ARG-helix/CYS-sheet": -Infinity, "ARG-helix/ASP-helix": 0.889, "ARG-helix/ASP-coil": -0.678, "ARG-helix/ASP-sheet": -2.527, "ARG-helix/SER-helix": 0.014, "ARG-helix/SER-coil": -0.962, "ARG-helix/SER-sheet": -3.246, "ARG-helix/GLN-helix": 1.044, "ARG-helix/GLN-coil": 0.349, "ARG-helix/GLN-sheet": -4.238, "ARG-helix/LYS-helix": 0.488, "ARG-helix/LYS-coil": -1.443, "ARG-helix/LYS-sheet": -Infinity, "ARG-helix/ILE-helix": -0.598, "ARG-helix/ILE-coil": -1.045, "ARG-helix/ILE-sheet": -2.594, "ARG-helix/PRO-helix": -0.414, "ARG-helix/PRO-coil": -1.717, "ARG-helix/PRO-sheet": -3.377, "ARG-helix/THR-helix": 0.129, "ARG-helix/THR-coil": -2.152, "ARG-helix/THR-sheet": -0.345, "ARG-helix/PHE-helix": -0.449, "ARG-helix/PHE-coil": -1.022, "ARG-helix/PHE-sheet": -4.808, "ARG-helix/ASN-helix": 0.969, "ARG-helix/ASN-coil": -0.947, "ARG-helix/ASN-sheet": -Infinity, "ARG-helix/GLY-helix": 0.194, "ARG-helix/GLY-coil": -1.089, "ARG-helix/GLY-sheet": -4.705, "ARG-helix/HIS-helix": 0.368, "ARG-helix/HIS-coil": -1.561, "ARG-helix/HIS-sheet": -3.282, "ARG-helix/LEU-helix": 0.178, "ARG-helix/LEU-coil": -0.502, "ARG-helix/LEU-sheet": -5.336, "ARG-helix/ARG-helix": 3.521, "ARG-helix/ARG-coil": -1.562, "ARG-helix/ARG-sheet": -3.329, "ARG-helix/TRP-helix": 0.023, "ARG-helix/TRP-coil": 0.349, "ARG-helix/TRP-sheet": -4.593, "ARG-helix/ALA-helix": 0.059, "ARG-helix/ALA-coil": -1.333, "ARG-helix/ALA-sheet": -4.814, "ARG-helix/VAL-helix": -0.229, "ARG-helix/VAL-coil": -0.626, "ARG-helix/VAL-sheet": -5.957, "ARG-helix/GLU-helix": 1.409, "ARG-helix/GLU-coil": -0.792, "ARG-helix/GLU-sheet": -2.678, "ARG-helix/TYR-helix": 0.238, "ARG-helix/TYR-coil": -1.664, "ARG-helix/TYR-sheet": -2.708, "ARG-helix/MET-helix": -0.354, "ARG-helix/MET-coil": -0.516, "ARG-helix/MET-sheet": -3.053, "ARG-coil/CYS-helix": -0.157, "ARG-coil/CYS-coil": -0.377, "ARG-coil/CYS-sheet": -Infinity, "ARG-coil/ASP-helix": -0.812, "ARG-coil/ASP-coil": 0.201, "ARG-coil/ASP-sheet": 0.954, "ARG-coil/SER-helix": -0.908, "ARG-coil/SER-coil": 0.323, "ARG-coil/SER-sheet": -2.368, "ARG-coil/GLN-helix": -1.406, "ARG-coil/GLN-coil": 1.048, "ARG-coil/GLN-sheet": -2.954, "ARG-coil/LYS-helix": -0.907, "ARG-coil/LYS-coil": 0.201, "ARG-coil/LYS-sheet": 0.267, "ARG-coil/ILE-helix": -1.84, "ARG-coil/ILE-coil": 0.29, "ARG-coil/ILE-sheet": -3.562, "ARG-coil/PRO-helix": -1.228, "ARG-coil/PRO-coil": -0.554, "ARG-coil/PRO-sheet": -2.787, "ARG-coil/THR-helix": -1.76, "ARG-coil/THR-coil": 0.238, "ARG-coil/THR-sheet": -0.577, "ARG-coil/PHE-helix": -2.061, "ARG-coil/PHE-coil": 0.538, "ARG-coil/PHE-sheet": -1.019, "ARG-coil/ASN-helix": -0.922, "ARG-coil/ASN-coil": -0.162, "ARG-coil/ASN-sheet": -0.74, "ARG-coil/GLY-helix": -0.822, "ARG-coil/GLY-coil": -0.049, "ARG-coil/GLY-sheet": -1.17, "ARG-coil/HIS-helix": -1.078, "ARG-coil/HIS-coil": 0.121, "ARG-coil/HIS-sheet": -2.692, "ARG-coil/LEU-helix": -1.071, "ARG-coil/LEU-coil": -0.309, "ARG-coil/LEU-sheet": -1.219, "ARG-coil/ARG-helix": -1.562, "ARG-coil/ARG-coil": 3.62, "ARG-coil/ARG-sheet": -1.213, "ARG-coil/TRP-helix": -1.368, "ARG-coil/TRP-coil": 0.242, "ARG-coil/TRP-sheet": -1.806, "ARG-coil/ALA-helix": -1.701, "ARG-coil/ALA-coil": 0.761, "ARG-coil/ALA-sheet": -4.224, "ARG-coil/VAL-helix": -1.335, "ARG-coil/VAL-coil": 0.336, "ARG-coil/VAL-sheet": -2.071, "ARG-coil/GLU-helix": -0.471, "ARG-coil/GLU-coil": 0.268, "ARG-coil/GLU-sheet": 0.745, "ARG-coil/TYR-helix": -0.768, "ARG-coil/TYR-coil": 0.384, "ARG-coil/TYR-sheet": 0.618, "ARG-coil/MET-helix": -2.607, "ARG-coil/MET-coil": -1.293, "ARG-coil/MET-sheet": -1.952, "ARG-sheet/CYS-helix": -Infinity, "ARG-sheet/CYS-coil": -1.645, "ARG-sheet/CYS-sheet": -1.933, "ARG-sheet/ASP-helix": -3.11, "ARG-sheet/ASP-coil": -0.14, "ARG-sheet/ASP-sheet": 2.456, "ARG-sheet/SER-helix": -2.89, "ARG-sheet/SER-coil": -0.353, "ARG-sheet/SER-sheet": 2.207, "ARG-sheet/GLN-helix": -3.079, "ARG-sheet/GLN-coil": 1.063, "ARG-sheet/GLN-sheet": 2.57, "ARG-sheet/LYS-helix": -4.791, "ARG-sheet/LYS-coil": -0.709, "ARG-sheet/LYS-sheet": -0.187, "ARG-sheet/ILE-helix": -4.454, "ARG-sheet/ILE-coil": -0.004, "ARG-sheet/ILE-sheet": 1.474, "ARG-sheet/PRO-helix": -Infinity, "ARG-sheet/PRO-coil": 0.421, "ARG-sheet/PRO-sheet": -1.529, "ARG-sheet/THR-helix": -2.391, "ARG-sheet/THR-coil": -0.033, "ARG-sheet/THR-sheet": 1.404, "ARG-sheet/PHE-helix": -2.297, "ARG-sheet/PHE-coil": -1.813, "ARG-sheet/PHE-sheet": 1.494, "ARG-sheet/ASN-helix": 0.173, "ARG-sheet/ASN-coil": 0.485, "ARG-sheet/ASN-sheet": 2.369, "ARG-sheet/GLY-helix": -3.053, "ARG-sheet/GLY-coil": -2.081, "ARG-sheet/GLY-sheet": 0.699, "ARG-sheet/HIS-helix": -3.732, "ARG-sheet/HIS-coil": 0.961, "ARG-sheet/HIS-sheet": -0.182, "ARG-sheet/LEU-helix": -5.076, "ARG-sheet/LEU-coil": -0.694, "ARG-sheet/LEU-sheet": 1.112, "ARG-sheet/ARG-helix": -3.329, "ARG-sheet/ARG-coil": -1.213, "ARG-sheet/ARG-sheet": 4.73, "ARG-sheet/TRP-helix": -2.963, "ARG-sheet/TRP-coil": -1.353, "ARG-sheet/TRP-sheet": 1.889, "ARG-sheet/ALA-helix": -3.816, "ARG-sheet/ALA-coil": -0.313, "ARG-sheet/ALA-sheet": 1.649, "ARG-sheet/VAL-helix": -4.566, "ARG-sheet/VAL-coil": -0.323, "ARG-sheet/VAL-sheet": -0.777, "ARG-sheet/GLU-helix": -2.941, "ARG-sheet/GLU-coil": 0.079, "ARG-sheet/GLU-sheet": 1.249, "ARG-sheet/TYR-helix": -3.128, "ARG-sheet/TYR-coil": 0.915, "ARG-sheet/TYR-sheet": 1.391, "ARG-sheet/MET-helix": -Infinity, "ARG-sheet/MET-coil": -1.907, "ARG-sheet/MET-sheet": 2.229, "TRP-helix/CYS-helix": 0.09, "TRP-helix/CYS-coil": -Infinity, "TRP-helix/CYS-sheet": -Infinity, "TRP-helix/ASP-helix": -0.443, "TRP-helix/ASP-coil": -0.247, "TRP-helix/ASP-sheet": -3.259, "TRP-helix/SER-helix": -0.338, "TRP-helix/SER-coil": -2.479, "TRP-helix/SER-sheet": -4.566, "TRP-helix/GLN-helix": 0.812, "TRP-helix/GLN-coil": -1.366, "TRP-helix/GLN-sheet": -Infinity, "TRP-helix/LYS-helix": 0.596, "TRP-helix/LYS-coil": -2.199, "TRP-helix/LYS-sheet": -3.723, "TRP-helix/ILE-helix": -0.444, "TRP-helix/ILE-coil": 0.182, "TRP-helix/ILE-sheet": -3.052, "TRP-helix/PRO-helix": -0.676, "TRP-helix/PRO-coil": 0.526, "TRP-helix/PRO-sheet": -2.5, "TRP-helix/THR-helix": -0.104, "TRP-helix/THR-coil": -2.988, "TRP-helix/THR-sheet": -Infinity, "TRP-helix/PHE-helix": -0.845, "TRP-helix/PHE-coil": -0.955, "TRP-helix/PHE-sheet": -1.58, "TRP-helix/ASN-helix": -0.871, "TRP-helix/ASN-coil": -1.402, "TRP-helix/ASN-sheet": -4.009, "TRP-helix/GLY-helix": -0.416, "TRP-helix/GLY-coil": -1.03, "TRP-helix/GLY-sheet": -3.828, "TRP-helix/HIS-helix": -0.469, "TRP-helix/HIS-coil": -1.71, "TRP-helix/HIS-sheet": -3.099, "TRP-helix/LEU-helix": -0.184, "TRP-helix/LEU-coil": -0.337, "TRP-helix/LEU-sheet": 0.341, "TRP-helix/ARG-helix": 0.023, "TRP-helix/ARG-coil": -1.368, "TRP-helix/ARG-sheet": -2.963, "TRP-helix/TRP-helix": 4.274, "TRP-helix/TRP-coil": 1.147, "TRP-helix/TRP-sheet": -Infinity, "TRP-helix/ALA-helix": 1.29, "TRP-helix/ALA-coil": -2.206, "TRP-helix/ALA-sheet": -3.021, "TRP-helix/VAL-helix": 0.301, "TRP-helix/VAL-coil": -2.104, "TRP-helix/VAL-sheet": -3.982, "TRP-helix/GLU-helix": -0.766, "TRP-helix/GLU-coil": -1.139, "TRP-helix/GLU-sheet": -Infinity, "TRP-helix/TYR-helix": 0.264, "TRP-helix/TYR-coil": 1.38, "TRP-helix/TYR-sheet": -Infinity, "TRP-helix/MET-helix": -0.085, "TRP-helix/MET-coil": -1.866, "TRP-helix/MET-sheet": -3.275, "TRP-coil/CYS-helix": -2.605, "TRP-coil/CYS-coil": -0.601, "TRP-coil/CYS-sheet": -Infinity, "TRP-coil/ASP-helix": -3.164, "TRP-coil/ASP-coil": -0.619, "TRP-coil/ASP-sheet": -Infinity, "TRP-coil/SER-helix": -0.546, "TRP-coil/SER-coil": -0.848, "TRP-coil/SER-sheet": -2.551, "TRP-coil/GLN-helix": -2.728, "TRP-coil/GLN-coil": -1.17, "TRP-coil/GLN-sheet": -2.444, "TRP-coil/LYS-helix": -3.053, "TRP-coil/LYS-coil": -2.555, "TRP-coil/LYS-sheet": -1.554, "TRP-coil/ILE-helix": -3.004, "TRP-coil/ILE-coil": 0.306, "TRP-coil/ILE-sheet": -3.051, "TRP-coil/PRO-helix": 0.405, "TRP-coil/PRO-coil": 0.42, "TRP-coil/PRO-sheet": -2.276, "TRP-coil/THR-helix": -1.139, "TRP-coil/THR-coil": 0.326, "TRP-coil/THR-sheet": 0.318, "TRP-coil/PHE-helix": 0.351, "TRP-coil/PHE-coil": -1.12, "TRP-coil/PHE-sheet": -0.663, "TRP-coil/ASN-helix": -1.926, "TRP-coil/ASN-coil": -1.138, "TRP-coil/ASN-sheet": -1.483, "TRP-coil/GLY-helix": -0.224, "TRP-coil/GLY-coil": -0.12, "TRP-coil/GLY-sheet": -1.995, "TRP-coil/HIS-helix": -Infinity, "TRP-coil/HIS-coil": -0.183, "TRP-coil/HIS-sheet": -1.489, "TRP-coil/LEU-helix": -0.29, "TRP-coil/LEU-coil": -1.4, "TRP-coil/LEU-sheet": 0.172, "TRP-coil/ARG-helix": 0.349, "TRP-coil/ARG-coil": 0.242, "TRP-coil/ARG-sheet": -1.353, "TRP-coil/TRP-helix": 1.147, "TRP-coil/TRP-coil": 5.425, "TRP-coil/TRP-sheet": -Infinity, "TRP-coil/ALA-helix": 0.59, "TRP-coil/ALA-coil": -0.482, "TRP-coil/ALA-sheet": -3.713, "TRP-coil/VAL-helix": -2.915, "TRP-coil/VAL-coil": 0.454, "TRP-coil/VAL-sheet": 0.461, "TRP-coil/GLU-helix": -3.4, "TRP-coil/GLU-coil": -0.222, "TRP-coil/GLU-sheet": -1.896, "TRP-coil/TYR-helix": -0.762, "TRP-coil/TYR-coil": 1.684, "TRP-coil/TYR-sheet": -1.762, "TRP-coil/MET-helix": 0.132, "TRP-coil/MET-coil": -Infinity, "TRP-coil/MET-sheet": -1.665, "TRP-sheet/CYS-helix": -Infinity, "TRP-sheet/CYS-coil": -Infinity, "TRP-sheet/CYS-sheet": -0.489, "TRP-sheet/ASP-helix": -Infinity, "TRP-sheet/ASP-coil": 0.053, "TRP-sheet/ASP-sheet": -1.538, "TRP-sheet/SER-helix": -Infinity, "TRP-sheet/SER-coil": 0.679, "TRP-sheet/SER-sheet": 0.118, "TRP-sheet/GLN-helix": -3.427, "TRP-sheet/GLN-coil": -2.562, "TRP-sheet/GLN-sheet": -0.54, "TRP-sheet/LYS-helix": -3.346, "TRP-sheet/LYS-coil": -0.7, "TRP-sheet/LYS-sheet": 0.832, "TRP-sheet/ILE-helix": 0.728, "TRP-sheet/ILE-coil": -Infinity, "TRP-sheet/ILE-sheet": -0.049, "TRP-sheet/PRO-helix": -2.264, "TRP-sheet/PRO-coil": 0.503, "TRP-sheet/PRO-sheet": 4.127, "TRP-sheet/THR-helix": -3.918, "TRP-sheet/THR-coil": -0.891, "TRP-sheet/THR-sheet": 1.734, "TRP-sheet/PHE-helix": -4.349, "TRP-sheet/PHE-coil": -Infinity, "TRP-sheet/PHE-sheet": 1.225, "TRP-sheet/ASN-helix": -3.924, "TRP-sheet/ASN-coil": 0.796, "TRP-sheet/ASN-sheet": 1.612, "TRP-sheet/GLY-helix": -3.96, "TRP-sheet/GLY-coil": 0.962, "TRP-sheet/GLY-sheet": 1.558, "TRP-sheet/HIS-helix": -Infinity, "TRP-sheet/HIS-coil": -1.729, "TRP-sheet/HIS-sheet": -0.684, "TRP-sheet/LEU-helix": 0.291, "TRP-sheet/LEU-coil": 0.601, "TRP-sheet/LEU-sheet": 1.427, "TRP-sheet/ARG-helix": -4.593, "TRP-sheet/ARG-coil": -1.806, "TRP-sheet/ARG-sheet": 1.889, "TRP-sheet/TRP-helix": -Infinity, "TRP-sheet/TRP-coil": -Infinity, "TRP-sheet/TRP-sheet": 5.528, "TRP-sheet/ALA-helix": -3.471, "TRP-sheet/ALA-coil": 0.79, "TRP-sheet/ALA-sheet": 1.816, "TRP-sheet/VAL-helix": -Infinity, "TRP-sheet/VAL-coil": -2.174, "TRP-sheet/VAL-sheet": 2.202, "TRP-sheet/GLU-helix": -Infinity, "TRP-sheet/GLU-coil": -0.674, "TRP-sheet/GLU-sheet": 0.879, "TRP-sheet/TYR-helix": -4.169, "TRP-sheet/TYR-coil": -Infinity, "TRP-sheet/TYR-sheet": 2.169, "TRP-sheet/MET-helix": 0.467, "TRP-sheet/MET-coil": -2.255, "TRP-sheet/MET-sheet": -0.012, "ALA-helix/CYS-helix": 0.688, "ALA-helix/CYS-coil": -2.458, "ALA-helix/CYS-sheet": -3.064, "ALA-helix/ASP-helix": -0.362, "ALA-helix/ASP-coil": -1.554, "ALA-helix/ASP-sheet": -4.4, "ALA-helix/SER-helix": 1.127, "ALA-helix/SER-coil": -0.676, "ALA-helix/SER-sheet": -3.473, "ALA-helix/GLN-helix": 0.524, "ALA-helix/GLN-coil": -0.211, "ALA-helix/GLN-sheet": -5.6, "ALA-helix/LYS-helix": -0.484, "ALA-helix/LYS-coil": -1.864, "ALA-helix/LYS-sheet": -3.565, "ALA-helix/ILE-helix": 0.562, "ALA-helix/ILE-coil": 0.204, "ALA-helix/ILE-sheet": -0.249, "ALA-helix/PRO-helix": 0.291, "ALA-helix/PRO-coil": -1.463, "ALA-helix/PRO-sheet": -3.641, "ALA-helix/THR-helix": 0.319, "ALA-helix/THR-coil": -0.61, "ALA-helix/THR-sheet": -3.203, "ALA-helix/PHE-helix": 0.166, "ALA-helix/PHE-coil": -2.743, "ALA-helix/PHE-sheet": -2.044, "ALA-helix/ASN-helix": -0.552, "ALA-helix/ASN-coil": -1.055, "ALA-helix/ASN-sheet": -3.609, "ALA-helix/GLY-helix": 1.817, "ALA-helix/GLY-coil": -0.855, "ALA-helix/GLY-sheet": -3.988, "ALA-helix/HIS-helix": -0.438, "ALA-helix/HIS-coil": -1.606, "ALA-helix/HIS-sheet": -3.036, "ALA-helix/LEU-helix": 0.68, "ALA-helix/LEU-coil": -0.287, "ALA-helix/LEU-sheet": -1.971, "ALA-helix/ARG-helix": 0.059, "ALA-helix/ARG-coil": -1.701, "ALA-helix/ARG-sheet": -3.816, "ALA-helix/TRP-helix": 1.29, "ALA-helix/TRP-coil": 0.59, "ALA-helix/TRP-sheet": -3.471, "ALA-helix/ALA-helix": 3.025, "ALA-helix/ALA-coil": -0.591, "ALA-helix/ALA-sheet": -0.076, "ALA-helix/VAL-helix": 0.622, "ALA-helix/VAL-coil": -1.779, "ALA-helix/VAL-sheet": -0.798, "ALA-helix/GLU-helix": 0.185, "ALA-helix/GLU-coil": -0.077, "ALA-helix/GLU-sheet": -4.359, "ALA-helix/TYR-helix": 0.67, "ALA-helix/TYR-coil": -0.39, "ALA-helix/TYR-sheet": -3.532, "ALA-helix/MET-helix": 0.732, "ALA-helix/MET-coil": -1.023, "ALA-helix/MET-sheet": -2.256, "ALA-coil/CYS-helix": -1.994, "ALA-coil/CYS-coil": -1.089, "ALA-coil/CYS-sheet": -0.801, "ALA-coil/ASP-helix": -1.683, "ALA-coil/ASP-coil": 0.974, "ALA-coil/ASP-sheet": 1.107, "ALA-coil/SER-helix": -1.194, "ALA-coil/SER-coil": 1.11, "ALA-coil/SER-sheet": -1.599, "ALA-coil/GLN-helix": -0.055, "ALA-coil/GLN-coil": 1.241, "ALA-coil/GLN-sheet": 1.398, "ALA-coil/LYS-helix": -0.676, "ALA-coil/LYS-coil": 0.608, "ALA-coil/LYS-sheet": -0.386, "ALA-coil/ILE-helix": -0.58, "ALA-coil/ILE-coil": 0.473, "ALA-coil/ILE-sheet": -0.079, "ALA-coil/PRO-helix": -1.999, "ALA-coil/PRO-coil": 1.56, "ALA-coil/PRO-sheet": -1.666, "ALA-coil/THR-helix": -0.005, "ALA-coil/THR-coil": 1.443, "ALA-coil/THR-sheet": 0.216, "ALA-coil/PHE-helix": -0.979, "ALA-coil/PHE-coil": 0.614, "ALA-coil/PHE-sheet": -0.964, "ALA-coil/ASN-helix": -1.999, "ALA-coil/ASN-coil": 0.219, "ALA-coil/ASN-sheet": 1.29, "ALA-coil/GLY-helix": -0.205, "ALA-coil/GLY-coil": 1.156, "ALA-coil/GLY-sheet": -0.748, "ALA-coil/HIS-helix": -0.169, "ALA-coil/HIS-coil": 0.564, "ALA-coil/HIS-sheet": -2.383, "ALA-coil/LEU-helix": -0.214, "ALA-coil/LEU-coil": 1.293, "ALA-coil/LEU-sheet": 0.919, "ALA-coil/ARG-helix": -1.333, "ALA-coil/ARG-coil": 0.761, "ALA-coil/ARG-sheet": -0.313, "ALA-coil/TRP-helix": -2.206, "ALA-coil/TRP-coil": -0.482, "ALA-coil/TRP-sheet": 0.79, "ALA-coil/ALA-helix": -0.591, "ALA-coil/ALA-coil": 3.519, "ALA-coil/ALA-sheet": 0.699, "ALA-coil/VAL-helix": -0.23, "ALA-coil/VAL-coil": 0.843, "ALA-coil/VAL-sheet": 0.246, "ALA-coil/GLU-helix": -1.39, "ALA-coil/GLU-coil": -0.077, "ALA-coil/GLU-sheet": 1.178, "ALA-coil/TYR-helix": -0.403, "ALA-coil/TYR-coil": 0.795, "ALA-coil/TYR-sheet": -1.644, "ALA-coil/MET-helix": -0.572, "ALA-coil/MET-coil": 0.818, "ALA-coil/MET-sheet": -1.999, "ALA-sheet/CYS-helix": -3.119, "ALA-sheet/CYS-coil": -0.961, "ALA-sheet/CYS-sheet": 2.549, "ALA-sheet/ASP-helix": -5.064, "ALA-sheet/ASP-coil": -0.167, "ALA-sheet/ASP-sheet": 2.468, "ALA-sheet/SER-helix": -3.235, "ALA-sheet/SER-coil": 0.379, "ALA-sheet/SER-sheet": 2.691, "ALA-sheet/GLN-helix": -3.425, "ALA-sheet/GLN-coil": 1.668, "ALA-sheet/GLN-sheet": 2.885, "ALA-sheet/LYS-helix": -4.666, "ALA-sheet/LYS-coil": -0.826, "ALA-sheet/LYS-sheet": 2.248, "ALA-sheet/ILE-helix": -0.231, "ALA-sheet/ILE-coil": -1.651, "ALA-sheet/ILE-sheet": 1.877, "ALA-sheet/PRO-helix": -3.178, "ALA-sheet/PRO-coil": -0.523, "ALA-sheet/PRO-sheet": 3.03, "ALA-sheet/THR-helix": -3.222, "ALA-sheet/THR-coil": -1.725, "ALA-sheet/THR-sheet": 2.205, "ALA-sheet/PHE-helix": -4.57, "ALA-sheet/PHE-coil": -1.87, "ALA-sheet/PHE-sheet": 3.08, "ALA-sheet/ASN-helix": -3.046, "ALA-sheet/ASN-coil": -2.099, "ALA-sheet/ASN-sheet": 2.439, "ALA-sheet/GLY-helix": -4.181, "ALA-sheet/GLY-coil": -0.106, "ALA-sheet/GLY-sheet": 3.444, "ALA-sheet/HIS-helix": -3.202, "ALA-sheet/HIS-coil": -1.256, "ALA-sheet/HIS-sheet": 1.347, "ALA-sheet/LEU-helix": -0.244, "ALA-sheet/LEU-coil": 0.388, "ALA-sheet/LEU-sheet": 2.447, "ALA-sheet/ARG-helix": -4.814, "ALA-sheet/ARG-coil": -4.224, "ALA-sheet/ARG-sheet": 1.649, "ALA-sheet/TRP-helix": -3.021, "ALA-sheet/TRP-coil": -3.713, "ALA-sheet/TRP-sheet": 1.816, "ALA-sheet/ALA-helix": -0.076, "ALA-sheet/ALA-coil": 0.699, "ALA-sheet/ALA-sheet": 5.019, "ALA-sheet/VAL-helix": -3.12, "ALA-sheet/VAL-coil": 1.349, "ALA-sheet/VAL-sheet": 2.133, "ALA-sheet/GLU-helix": -3.914, "ALA-sheet/GLU-coil": -1.029, "ALA-sheet/GLU-sheet": 1.959, "ALA-sheet/TYR-helix": -3.291, "ALA-sheet/TYR-coil": -1.653, "ALA-sheet/TYR-sheet": 2.336, "ALA-sheet/MET-helix": -2.112, "ALA-sheet/MET-coil": -2.476, "ALA-sheet/MET-sheet": 0.495, "VAL-helix/CYS-helix": 0.749, "VAL-helix/CYS-coil": 1.271, "VAL-helix/CYS-sheet": -3.408, "VAL-helix/ASP-helix": -0.493, "VAL-helix/ASP-coil": -1.372, "VAL-helix/ASP-sheet": -4.457, "VAL-helix/SER-helix": 1.025, "VAL-helix/SER-coil": -0.937, "VAL-helix/SER-sheet": -3.818, "VAL-helix/GLN-helix": 0.067, "VAL-helix/GLN-coil": 0.219, "VAL-helix/GLN-sheet": -0.809, "VAL-helix/LYS-helix": -0.564, "VAL-helix/LYS-coil": -3.248, "VAL-helix/LYS-sheet": -1.852, "VAL-helix/ILE-helix": 0.859, "VAL-helix/ILE-coil": 0.14, "VAL-helix/ILE-sheet": -0.231, "VAL-helix/PRO-helix": -0.082, "VAL-helix/PRO-coil": -2.381, "VAL-helix/PRO-sheet": -2.599, "VAL-helix/THR-helix": 1.041, "VAL-helix/THR-coil": -1.436, "VAL-helix/THR-sheet": -4.127, "VAL-helix/PHE-helix": 0.146, "VAL-helix/PHE-coil": -0.309, "VAL-helix/PHE-sheet": -3.424, "VAL-helix/ASN-helix": 0.031, "VAL-helix/ASN-coil": -2.664, "VAL-helix/ASN-sheet": -4.108, "VAL-helix/GLY-helix": 1.349, "VAL-helix/GLY-coil": -0.931, "VAL-helix/GLY-sheet": -4.109, "VAL-helix/HIS-helix": -0.313, "VAL-helix/HIS-coil": -1.94, "VAL-helix/HIS-sheet": -3.198, "VAL-helix/LEU-helix": 0.738, "VAL-helix/LEU-coil": -1.551, "VAL-helix/LEU-sheet": -1.459, "VAL-helix/ARG-helix": -0.229, "VAL-helix/ARG-coil": -1.335, "VAL-helix/ARG-sheet": -4.566, "VAL-helix/TRP-helix": 0.301, "VAL-helix/TRP-coil": -2.915, "VAL-helix/TRP-sheet": -Infinity, "VAL-helix/ALA-helix": 0.622, "VAL-helix/ALA-coil": -0.23, "VAL-helix/ALA-sheet": -3.12, "VAL-helix/VAL-helix": 3.418, "VAL-helix/VAL-coil": -0.058, "VAL-helix/VAL-sheet": -3.233, "VAL-helix/GLU-helix": 0.087, "VAL-helix/GLU-coil": -1.842, "VAL-helix/GLU-sheet": -Infinity, "VAL-helix/TYR-helix": 0.491, "VAL-helix/TYR-coil": -3.371, "VAL-helix/TYR-sheet": -3.589, "VAL-helix/MET-helix": 0.768, "VAL-helix/MET-coil": -0.951, "VAL-helix/MET-sheet": -0.196, "VAL-coil/CYS-helix": -2.855, "VAL-coil/CYS-coil": -0.301, "VAL-coil/CYS-sheet": -2.46, "VAL-coil/ASP-helix": -1.36, "VAL-coil/ASP-coil": 1.54, "VAL-coil/ASP-sheet": -Infinity, "VAL-coil/SER-helix": -0.814, "VAL-coil/SER-coil": 0.841, "VAL-coil/SER-sheet": -0.004, "VAL-coil/GLN-helix": -1.452, "VAL-coil/GLN-coil": 0.487, "VAL-coil/GLN-sheet": -1.008, "VAL-coil/LYS-helix": -1.241, "VAL-coil/LYS-coil": 0.769, "VAL-coil/LYS-sheet": -0.572, "VAL-coil/ILE-helix": -0.684, "VAL-coil/ILE-coil": 1.894, "VAL-coil/ILE-sheet": -1.967, "VAL-coil/PRO-helix": -1.69, "VAL-coil/PRO-coil": 0.839, "VAL-coil/PRO-sheet": -2.75, "VAL-coil/THR-helix": -0.018, "VAL-coil/THR-coil": 0.782, "VAL-coil/THR-sheet": 0.575, "VAL-coil/PHE-helix": -2.131, "VAL-coil/PHE-coil": 0.833, "VAL-coil/PHE-sheet": -0.156, "VAL-coil/ASN-helix": -1.994, "VAL-coil/ASN-coil": 0.764, "VAL-coil/ASN-sheet": 0.316, "VAL-coil/GLY-helix": -1.163, "VAL-coil/GLY-coil": 1.426, "VAL-coil/GLY-sheet": -0.566, "VAL-coil/HIS-helix": -1.315, "VAL-coil/HIS-coil": -1.378, "VAL-coil/HIS-sheet": -2.25, "VAL-coil/LEU-helix": -0.555, "VAL-coil/LEU-coil": 0.649, "VAL-coil/LEU-sheet": -1.125, "VAL-coil/ARG-helix": -0.626, "VAL-coil/ARG-coil": 0.336, "VAL-coil/ARG-sheet": -0.323, "VAL-coil/TRP-helix": -2.104, "VAL-coil/TRP-coil": 0.454, "VAL-coil/TRP-sheet": -2.174, "VAL-coil/ALA-helix": -1.779, "VAL-coil/ALA-coil": 0.843, "VAL-coil/ALA-sheet": 1.349, "VAL-coil/VAL-helix": -0.058, "VAL-coil/VAL-coil": 3.828, "VAL-coil/VAL-sheet": 0.012, "VAL-coil/GLU-helix": -1.075, "VAL-coil/GLU-coil": 0.88, "VAL-coil/GLU-sheet": 0.922, "VAL-coil/TYR-helix": -2.516, "VAL-coil/TYR-coil": 0.353, "VAL-coil/TYR-sheet": -0.458, "VAL-coil/MET-helix": -1.343, "VAL-coil/MET-coil": 1.092, "VAL-coil/MET-sheet": 0.941, "VAL-sheet/CYS-helix": -2.065, "VAL-sheet/CYS-coil": -1.748, "VAL-sheet/CYS-sheet": 0.632, "VAL-sheet/ASP-helix": -4.416, "VAL-sheet/ASP-coil": -0.806, "VAL-sheet/ASP-sheet": -0.363, "VAL-sheet/SER-helix": -4.042, "VAL-sheet/SER-coil": -2.025, "VAL-sheet/SER-sheet": 0.174, "VAL-sheet/GLN-helix": -3.693, "VAL-sheet/GLN-coil": -1.934, "VAL-sheet/GLN-sheet": -0.364, "VAL-sheet/LYS-helix": -2.408, "VAL-sheet/LYS-coil": -2.421, "VAL-sheet/LYS-sheet": 1.301, "VAL-sheet/ILE-helix": -3.221, "VAL-sheet/ILE-coil": -1.351, "VAL-sheet/ILE-sheet": 2.077, "VAL-sheet/PRO-helix": -3.405, "VAL-sheet/PRO-coil": -0.495, "VAL-sheet/PRO-sheet": 0.63, "VAL-sheet/THR-helix": -3.41, "VAL-sheet/THR-coil": -2.726, "VAL-sheet/THR-sheet": 1.071, "VAL-sheet/PHE-helix": -3.922, "VAL-sheet/PHE-coil": -2.39, "VAL-sheet/PHE-sheet": 1.468, "VAL-sheet/ASN-helix": -4.19, "VAL-sheet/ASN-coil": -2.613, "VAL-sheet/ASN-sheet": 2.168, "VAL-sheet/GLY-helix": -3.619, "VAL-sheet/GLY-coil": 0.412, "VAL-sheet/GLY-sheet": 3.152, "VAL-sheet/HIS-helix": 1.236, "VAL-sheet/HIS-coil": -Infinity, "VAL-sheet/HIS-sheet": -0.581, "VAL-sheet/LEU-helix": -0.882, "VAL-sheet/LEU-coil": 0.942, "VAL-sheet/LEU-sheet": 2.189, "VAL-sheet/ARG-helix": -5.957, "VAL-sheet/ARG-coil": -2.071, "VAL-sheet/ARG-sheet": -0.777, "VAL-sheet/TRP-helix": -3.982, "VAL-sheet/TRP-coil": 0.461, "VAL-sheet/TRP-sheet": 2.202, "VAL-sheet/ALA-helix": -0.798, "VAL-sheet/ALA-coil": 0.246, "VAL-sheet/ALA-sheet": 2.133, "VAL-sheet/VAL-helix": -3.233, "VAL-sheet/VAL-coil": 0.012, "VAL-sheet/VAL-sheet": 4.189, "VAL-sheet/GLU-helix": -0.727, "VAL-sheet/GLU-coil": -2.285, "VAL-sheet/GLU-sheet": 1.98, "VAL-sheet/TYR-helix": -2.588, "VAL-sheet/TYR-coil": -0.369, "VAL-sheet/TYR-sheet": 1.334, "VAL-sheet/MET-helix": -0.096, "VAL-sheet/MET-coil": -4.312, "VAL-sheet/MET-sheet": 1.87, "GLU-helix/CYS-helix": -2.346, "GLU-helix/CYS-coil": -4.385, "GLU-helix/CYS-sheet": -Infinity, "GLU-helix/ASP-helix": 1.148, "GLU-helix/ASP-coil": -2.354, "GLU-helix/ASP-sheet": -5.029, "GLU-helix/SER-helix": 0.822, "GLU-helix/SER-coil": -1.343, "GLU-helix/SER-sheet": -2.464, "GLU-helix/GLN-helix": 0.588, "GLU-helix/GLN-coil": 0.168, "GLU-helix/GLN-sheet": -1.492, "GLU-helix/LYS-helix": 1.122, "GLU-helix/LYS-coil": -0.358, "GLU-helix/LYS-sheet": -3.883, "GLU-helix/ILE-helix": -0.433, "GLU-helix/ILE-coil": -2.827, "GLU-helix/ILE-sheet": -4.128, "GLU-helix/PRO-helix": -0.967, "GLU-helix/PRO-coil": -2.953, "GLU-helix/PRO-sheet": -2.477, "GLU-helix/THR-helix": 0.21, "GLU-helix/THR-coil": -0.509, "GLU-helix/THR-sheet": -4.006, "GLU-helix/PHE-helix": -1.282, "GLU-helix/PHE-coil": -1.979, "GLU-helix/PHE-sheet": -Infinity, "GLU-helix/ASN-helix": 1.001, "GLU-helix/ASN-coil": 0.313, "GLU-helix/ASN-sheet": -3.986, "GLU-helix/GLY-helix": -0.152, "GLU-helix/GLY-coil": -1.593, "GLU-helix/GLY-sheet": -4.21, "GLU-helix/HIS-helix": 1.112, "GLU-helix/HIS-coil": -0.458, "GLU-helix/HIS-sheet": -4.175, "GLU-helix/LEU-helix": 0.294, "GLU-helix/LEU-coil": -1.958, "GLU-helix/LEU-sheet": -4.282, "GLU-helix/ARG-helix": 1.409, "GLU-helix/ARG-coil": -0.471, "GLU-helix/ARG-sheet": -2.941, "GLU-helix/TRP-helix": -0.766, "GLU-helix/TRP-coil": -3.4, "GLU-helix/TRP-sheet": -Infinity, "GLU-helix/ALA-helix": 0.185, "GLU-helix/ALA-coil": -1.39, "GLU-helix/ALA-sheet": -3.914, "GLU-helix/VAL-helix": 0.087, "GLU-helix/VAL-coil": -1.075, "GLU-helix/VAL-sheet": -0.727, "GLU-helix/GLU-helix": 3.498, "GLU-helix/GLU-coil": 0.403, "GLU-helix/GLU-sheet": -Infinity, "GLU-helix/TYR-helix": 1.066, "GLU-helix/TYR-coil": -1.568, "GLU-helix/TYR-sheet": 0.102, "GLU-helix/MET-helix": -0.009, "GLU-helix/MET-coil": -1.197, "GLU-helix/MET-sheet": -Infinity, "GLU-coil/CYS-helix": -2.336, "GLU-coil/CYS-coil": 0.789, "GLU-coil/CYS-sheet": -Infinity, "GLU-coil/ASP-helix": -0.005, "GLU-coil/ASP-coil": 0.88, "GLU-coil/ASP-sheet": -4.376, "GLU-coil/SER-helix": -0.621, "GLU-coil/SER-coil": 1.07, "GLU-coil/SER-sheet": -2.592, "GLU-coil/GLN-helix": -1.556, "GLU-coil/GLN-coil": 0.934, "GLU-coil/GLN-sheet": -1.387, "GLU-coil/LYS-helix": -0.348, "GLU-coil/LYS-coil": 0.85, "GLU-coil/LYS-sheet": -0.186, "GLU-coil/ILE-helix": -1.199, "GLU-coil/ILE-coil": 0.621, "GLU-coil/ILE-sheet": -2.687, "GLU-coil/PRO-helix": -1.313, "GLU-coil/PRO-coil": 0.738, "GLU-coil/PRO-sheet": -Infinity, "GLU-coil/THR-helix": -1.908, "GLU-coil/THR-coil": 0.534, "GLU-coil/THR-sheet": -2.814, "GLU-coil/PHE-helix": -3.787, "GLU-coil/PHE-coil": 0.492, "GLU-coil/PHE-sheet": -0.433, "GLU-coil/ASN-helix": 0.511, "GLU-coil/ASN-coil": 0.808, "GLU-coil/ASN-sheet": -0.695, "GLU-coil/GLY-helix": -1.839, "GLU-coil/GLY-coil": 0.488, "GLU-coil/GLY-sheet": -0.937, "GLU-coil/HIS-helix": -0.968, "GLU-coil/HIS-coil": 1.154, "GLU-coil/HIS-sheet": -1.577, "GLU-coil/LEU-helix": -2.058, "GLU-coil/LEU-coil": 0.831, "GLU-coil/LEU-sheet": -1.965, "GLU-coil/ARG-helix": -0.792, "GLU-coil/ARG-coil": 0.268, "GLU-coil/ARG-sheet": 0.079, "GLU-coil/TRP-helix": -1.139, "GLU-coil/TRP-coil": -0.222, "GLU-coil/TRP-sheet": -0.674, "GLU-coil/ALA-helix": -0.077, "GLU-coil/ALA-coil": -0.077, "GLU-coil/ALA-sheet": -1.029, "GLU-coil/VAL-helix": -1.842, "GLU-coil/VAL-coil": 0.88, "GLU-coil/VAL-sheet": -2.285, "GLU-coil/GLU-helix": 0.403, "GLU-coil/GLU-coil": 3.598, "GLU-coil/GLU-sheet": -1.578, "GLU-coil/TYR-helix": -0.329, "GLU-coil/TYR-coil": 1.187, "GLU-coil/TYR-sheet": -0.564, "GLU-coil/MET-helix": -1.147, "GLU-coil/MET-coil": 0.187, "GLU-coil/MET-sheet": -0.991, "GLU-sheet/CYS-helix": -2.688, "GLU-sheet/CYS-coil": -1.782, "GLU-sheet/CYS-sheet": 0.127, "GLU-sheet/ASP-helix": -Infinity, "GLU-sheet/ASP-coil": 0.111, "GLU-sheet/ASP-sheet": 2.753, "GLU-sheet/SER-helix": -Infinity, "GLU-sheet/SER-coil": -0.027, "GLU-sheet/SER-sheet": 2.029, "GLU-sheet/GLN-helix": 1.087, "GLU-sheet/GLN-coil": -1.099, "GLU-sheet/GLN-sheet": 2.222, "GLU-sheet/LYS-helix": -3.136, "GLU-sheet/LYS-coil": 0.913, "GLU-sheet/LYS-sheet": 3.335, "GLU-sheet/ILE-helix": -4.591, "GLU-sheet/ILE-coil": -0.416, "GLU-sheet/ILE-sheet": 0.316, "GLU-sheet/PRO-helix": 0.4, "GLU-sheet/PRO-coil": -2.142, "GLU-sheet/PRO-sheet": 2.888, "GLU-sheet/THR-helix": -4.401, "GLU-sheet/THR-coil": 1.128, "GLU-sheet/THR-sheet": 3.187, "GLU-sheet/PHE-helix": -Infinity, "GLU-sheet/PHE-coil": -2.068, "GLU-sheet/PHE-sheet": 1.542, "GLU-sheet/ASN-helix": -Infinity, "GLU-sheet/ASN-coil": -1.496, "GLU-sheet/ASN-sheet": 2.108, "GLU-sheet/GLY-helix": -4.443, "GLU-sheet/GLY-coil": 0.112, "GLU-sheet/GLY-sheet": 1.733, "GLU-sheet/HIS-helix": -Infinity, "GLU-sheet/HIS-coil": -3.128, "GLU-sheet/HIS-sheet": -0.473, "GLU-sheet/LEU-helix": -2.014, "GLU-sheet/LEU-coil": 0.367, "GLU-sheet/LEU-sheet": 1.962, "GLU-sheet/ARG-helix": -2.678, "GLU-sheet/ARG-coil": 0.745, "GLU-sheet/ARG-sheet": 1.249, "GLU-sheet/TRP-helix": -Infinity, "GLU-sheet/TRP-coil": -1.896, "GLU-sheet/TRP-sheet": 0.879, "GLU-sheet/ALA-helix": -4.359, "GLU-sheet/ALA-coil": 1.178, "GLU-sheet/ALA-sheet": 1.959, "GLU-sheet/VAL-helix": -Infinity, "GLU-sheet/VAL-coil": 0.922, "GLU-sheet/VAL-sheet": 1.98, "GLU-sheet/GLU-helix": -Infinity, "GLU-sheet/GLU-coil": -1.578, "GLU-sheet/GLU-sheet": 4.954, "GLU-sheet/TYR-helix": -Infinity, "GLU-sheet/TYR-coil": -2.321, "GLU-sheet/TYR-sheet": 1.191, "GLU-sheet/MET-helix": -Infinity, "GLU-sheet/MET-coil": -1.033, "GLU-sheet/MET-sheet": -1.055, "TYR-helix/CYS-helix": -1.777, "TYR-helix/CYS-coil": -1.816, "TYR-helix/CYS-sheet": -1.277, "TYR-helix/ASP-helix": 0.387, "TYR-helix/ASP-coil": -2.257, "TYR-helix/ASP-sheet": -3.019, "TYR-helix/SER-helix": -1.274, "TYR-helix/SER-coil": -1.207, "TYR-helix/SER-sheet": -2.939, "TYR-helix/GLN-helix": 0.13, "TYR-helix/GLN-coil": -0.116, "TYR-helix/GLN-sheet": -Infinity, "TYR-helix/LYS-helix": 0.126, "TYR-helix/LYS-coil": -0.14, "TYR-helix/LYS-sheet": -4.869, "TYR-helix/ILE-helix": -0.5, "TYR-helix/ILE-coil": -0.996, "TYR-helix/ILE-sheet": -3.727, "TYR-helix/PRO-helix": 0.251, "TYR-helix/PRO-coil": -0.203, "TYR-helix/PRO-sheet": -2.259, "TYR-helix/THR-helix": 1.401, "TYR-helix/THR-coil": -1.166, "TYR-helix/THR-sheet": -1.185, "TYR-helix/PHE-helix": 0.345, "TYR-helix/PHE-coil": 0.076, "TYR-helix/PHE-sheet": -5.077, "TYR-helix/ASN-helix": -0.246, "TYR-helix/ASN-coil": -1.492, "TYR-helix/ASN-sheet": -3.768, "TYR-helix/GLY-helix": 0.069, "TYR-helix/GLY-coil": -1.418, "TYR-helix/GLY-sheet": -3.587, "TYR-helix/HIS-helix": -0.165, "TYR-helix/HIS-coil": -0.564, "TYR-helix/HIS-sheet": -2.453, "TYR-helix/LEU-helix": 0.276, "TYR-helix/LEU-coil": -1.911, "TYR-helix/LEU-sheet": -1.356, "TYR-helix/ARG-helix": 0.238, "TYR-helix/ARG-coil": -0.768, "TYR-helix/ARG-sheet": -3.128, "TYR-helix/TRP-helix": 0.264, "TYR-helix/TRP-coil": -0.762, "TYR-helix/TRP-sheet": -4.169, "TYR-helix/ALA-helix": 0.67, "TYR-helix/ALA-coil": -0.403, "TYR-helix/ALA-sheet": -3.291, "TYR-helix/VAL-helix": 0.491, "TYR-helix/VAL-coil": -2.516, "TYR-helix/VAL-sheet": -2.588, "TYR-helix/GLU-helix": 1.066, "TYR-helix/GLU-coil": -0.329, "TYR-helix/GLU-sheet": -Infinity, "TYR-helix/TYR-helix": 3.837, "TYR-helix/TYR-coil": 0.933, "TYR-helix/TYR-sheet": -2.844, "TYR-helix/MET-helix": 0.257, "TYR-helix/MET-coil": -2.232, "TYR-helix/MET-sheet": -3.727, "TYR-coil/CYS-helix": -1.238, "TYR-coil/CYS-coil": -2.53, "TYR-coil/CYS-sheet": -Infinity, "TYR-coil/ASP-helix": -0.653, "TYR-coil/ASP-coil": 0.488, "TYR-coil/ASP-sheet": -0.976, "TYR-coil/SER-helix": -1.126, "TYR-coil/SER-coil": 0.59, "TYR-coil/SER-sheet": -3.381, "TYR-coil/GLN-helix": 0.226, "TYR-coil/GLN-coil": 1.178, "TYR-coil/GLN-sheet": -3.967, "TYR-coil/LYS-helix": -1.868, "TYR-coil/LYS-coil": 0.214, "TYR-coil/LYS-sheet": -0.546, "TYR-coil/ILE-helix": -1.714, "TYR-coil/ILE-coil": -1.087, "TYR-coil/ILE-sheet": -3.476, "TYR-coil/PRO-helix": -3.271, "TYR-coil/PRO-coil": 0.481, "TYR-coil/PRO-sheet": -2.008, "TYR-coil/THR-helix": -0.81, "TYR-coil/THR-coil": 0.887, "TYR-coil/THR-sheet": -0.877, "TYR-coil/PHE-helix": -2.401, "TYR-coil/PHE-coil": 1.401, "TYR-coil/PHE-sheet": -1.973, "TYR-coil/ASN-helix": -1.858, "TYR-coil/ASN-coil": -0.267, "TYR-coil/ASN-sheet": 0.976, "TYR-coil/GLY-helix": -0.866, "TYR-coil/GLY-coil": 0.451, "TYR-coil/GLY-sheet": 1.758, "TYR-coil/HIS-helix": -0.81, "TYR-coil/HIS-coil": 0.086, "TYR-coil/HIS-sheet": -Infinity, "TYR-coil/LEU-helix": -3.475, "TYR-coil/LEU-coil": 0.275, "TYR-coil/LEU-sheet": -0.775, "TYR-coil/ARG-helix": -1.664, "TYR-coil/ARG-coil": 0.384, "TYR-coil/ARG-sheet": 0.915, "TYR-coil/TRP-helix": 1.38, "TYR-coil/TRP-coil": 1.684, "TYR-coil/TRP-sheet": -Infinity, "TYR-coil/ALA-helix": -0.39, "TYR-coil/ALA-coil": 0.795, "TYR-coil/ALA-sheet": -1.653, "TYR-coil/VAL-helix": -3.371, "TYR-coil/VAL-coil": 0.353, "TYR-coil/VAL-sheet": -0.369, "TYR-coil/GLU-helix": -1.568, "TYR-coil/GLU-coil": 1.187, "TYR-coil/GLU-sheet": -2.321, "TYR-coil/TYR-helix": 0.933, "TYR-coil/TYR-coil": 4.259, "TYR-coil/TYR-sheet": -0.918, "TYR-coil/MET-helix": -2.988, "TYR-coil/MET-coil": -0.959, "TYR-coil/MET-sheet": -3.188, "TYR-sheet/CYS-helix": -2.043, "TYR-sheet/CYS-coil": -2.054, "TYR-sheet/CYS-sheet": 1.59, "TYR-sheet/ASP-helix": -4.905, "TYR-sheet/ASP-coil": 0.519, "TYR-sheet/ASP-sheet": 0.481, "TYR-sheet/SER-helix": -3.432, "TYR-sheet/SER-coil": -0.055, "TYR-sheet/SER-sheet": 2.121, "TYR-sheet/GLN-helix": -4.181, "TYR-sheet/GLN-coil": -3.028, "TYR-sheet/GLN-sheet": 2.574, "TYR-sheet/LYS-helix": -4.506, "TYR-sheet/LYS-coil": -0.068, "TYR-sheet/LYS-sheet": 1.816, "TYR-sheet/ILE-helix": -3.071, "TYR-sheet/ILE-coil": 0.781, "TYR-sheet/ILE-sheet": 1.037, "TYR-sheet/PRO-helix": -2.458, "TYR-sheet/PRO-coil": -0.234, "TYR-sheet/PRO-sheet": 1.339, "TYR-sheet/THR-helix": -5.365, "TYR-sheet/THR-coil": -0.079, "TYR-sheet/THR-sheet": 1.611, "TYR-sheet/PHE-helix": -3.494, "TYR-sheet/PHE-coil": -1.087, "TYR-sheet/PHE-sheet": 2.304, "TYR-sheet/ASN-helix": -3.58, "TYR-sheet/ASN-coil": -0.927, "TYR-sheet/ASN-sheet": 1.008, "TYR-sheet/GLY-helix": -3.798, "TYR-sheet/GLY-coil": 0.255, "TYR-sheet/GLY-sheet": 2.809, "TYR-sheet/HIS-helix": -Infinity, "TYR-sheet/HIS-coil": 0.711, "TYR-sheet/HIS-sheet": 1.452, "TYR-sheet/LEU-helix": -2.622, "TYR-sheet/LEU-coil": 1.402, "TYR-sheet/LEU-sheet": 1.978, "TYR-sheet/ARG-helix": -2.708, "TYR-sheet/ARG-coil": 0.618, "TYR-sheet/ARG-sheet": 1.391, "TYR-sheet/TRP-helix": -Infinity, "TYR-sheet/TRP-coil": -1.762, "TYR-sheet/TRP-sheet": 2.169, "TYR-sheet/ALA-helix": -3.532, "TYR-sheet/ALA-coil": -1.644, "TYR-sheet/ALA-sheet": 2.336, "TYR-sheet/VAL-helix": -3.589, "TYR-sheet/VAL-coil": -0.458, "TYR-sheet/VAL-sheet": 1.334, "TYR-sheet/GLU-helix": 0.102, "TYR-sheet/GLU-coil": -0.564, "TYR-sheet/GLU-sheet": 1.191, "TYR-sheet/TYR-helix": -2.844, "TYR-sheet/TYR-coil": -0.918, "TYR-sheet/TYR-sheet": 4.661, "TYR-sheet/MET-helix": -2.394, "TYR-sheet/MET-coil": -0.658, "TYR-sheet/MET-sheet": 1.892, "MET-helix/CYS-helix": 0.641, "MET-helix/CYS-coil": -1.484, "MET-helix/CYS-sheet": -2.332, "MET-helix/ASP-helix": 0.171, "MET-helix/ASP-coil": -3.979, "MET-helix/ASP-sheet": -Infinity, "MET-helix/SER-helix": 0.003, "MET-helix/SER-coil": -2.368, "MET-helix/SER-sheet": -4.687, "MET-helix/GLN-helix": 0.429, "MET-helix/GLN-coil": -1.314, "MET-helix/GLN-sheet": -Infinity, "MET-helix/LYS-helix": -0.531, "MET-helix/LYS-coil": -2.32, "MET-helix/LYS-sheet": -Infinity, "MET-helix/ILE-helix": 0.815, "MET-helix/ILE-coil": -3.062, "MET-helix/ILE-sheet": -2.703, "MET-helix/PRO-helix": -0.513, "MET-helix/PRO-coil": -2.404, "MET-helix/PRO-sheet": -1.928, "MET-helix/THR-helix": -0.399, "MET-helix/THR-coil": -3.715, "MET-helix/THR-sheet": -4.997, "MET-helix/PHE-helix": 0.679, "MET-helix/PHE-coil": -1.77, "MET-helix/PHE-sheet": 0.538, "MET-helix/ASN-helix": 0.123, "MET-helix/ASN-coil": -0.377, "MET-helix/ASN-sheet": -3.031, "MET-helix/GLY-helix": -0.268, "MET-helix/GLY-coil": -1.438, "MET-helix/GLY-sheet": -2.445, "MET-helix/HIS-helix": 0.294, "MET-helix/HIS-coil": -3.166, "MET-helix/HIS-sheet": -2.121, "MET-helix/LEU-helix": 0.291, "MET-helix/LEU-coil": -0.394, "MET-helix/LEU-sheet": -3.481, "MET-helix/ARG-helix": -0.354, "MET-helix/ARG-coil": -2.607, "MET-helix/ARG-sheet": -Infinity, "MET-helix/TRP-helix": -0.085, "MET-helix/TRP-coil": 0.132, "MET-helix/TRP-sheet": 0.467, "MET-helix/ALA-helix": 0.732, "MET-helix/ALA-coil": -0.572, "MET-helix/ALA-sheet": -2.112, "MET-helix/VAL-helix": 0.768, "MET-helix/VAL-coil": -1.343, "MET-helix/VAL-sheet": -0.096, "MET-helix/GLU-helix": -0.009, "MET-helix/GLU-coil": -1.147, "MET-helix/GLU-sheet": -Infinity, "MET-helix/TYR-helix": 0.257, "MET-helix/TYR-coil": -2.988, "MET-helix/TYR-sheet": -2.394, "MET-helix/MET-helix": 4.326, "MET-helix/MET-coil": -2.776, "MET-helix/MET-sheet": -2.703, "MET-coil/CYS-helix": -Infinity, "MET-coil/CYS-coil": -1.442, "MET-coil/CYS-sheet": 0.503, "MET-coil/ASP-helix": -2.214, "MET-coil/ASP-coil": 0.311, "MET-coil/ASP-sheet": -2.086, "MET-coil/SER-helix": -0.753, "MET-coil/SER-coil": 0.487, "MET-coil/SER-sheet": -1.852, "MET-coil/GLN-helix": -0.024, "MET-coil/GLN-coil": 0.109, "MET-coil/GLN-sheet": 0.148, "MET-coil/LYS-helix": -1.042, "MET-coil/LYS-coil": -0.106, "MET-coil/LYS-sheet": -0.47, "MET-coil/ILE-helix": -0.031, "MET-coil/ILE-coil": 0.582, "MET-coil/ILE-sheet": 0.476, "MET-coil/PRO-helix": -2.589, "MET-coil/PRO-coil": -0.416, "MET-coil/PRO-sheet": -1.039, "MET-coil/THR-helix": -1.398, "MET-coil/THR-coil": 0.846, "MET-coil/THR-sheet": -2.498, "MET-coil/PHE-helix": -2.064, "MET-coil/PHE-coil": 0.253, "MET-coil/PHE-sheet": -3.856, "MET-coil/ASN-helix": -1.833, "MET-coil/ASN-coil": 0.265, "MET-coil/ASN-sheet": -2.547, "MET-coil/GLY-helix": -1.982, "MET-coil/GLY-coil": -0.403, "MET-coil/GLY-sheet": -2.143, "MET-coil/HIS-helix": -1.989, "MET-coil/HIS-coil": -0.591, "MET-coil/HIS-sheet": -1.232, "MET-coil/LEU-helix": -0.573, "MET-coil/LEU-coil": 0.16, "MET-coil/LEU-sheet": -2.081, "MET-coil/ARG-helix": -0.516, "MET-coil/ARG-coil": -1.293, "MET-coil/ARG-sheet": -1.907, "MET-coil/TRP-helix": -1.866, "MET-coil/TRP-coil": -Infinity, "MET-coil/TRP-sheet": -2.255, "MET-coil/ALA-helix": -1.023, "MET-coil/ALA-coil": 0.818, "MET-coil/ALA-sheet": -2.476, "MET-coil/VAL-helix": -0.951, "MET-coil/VAL-coil": 1.092, "MET-coil/VAL-sheet": -4.312, "MET-coil/GLU-helix": -1.197, "MET-coil/GLU-coil": 0.187, "MET-coil/GLU-sheet": -1.033, "MET-coil/TYR-helix": -2.232, "MET-coil/TYR-coil": -0.959, "MET-coil/TYR-sheet": -0.658, "MET-coil/MET-helix": -2.776, "MET-coil/MET-coil": 4.37, "MET-coil/MET-sheet": -2.507, "MET-sheet/CYS-helix": -1.764, "MET-sheet/CYS-coil": -1.552, "MET-sheet/CYS-sheet": 2.786, "MET-sheet/ASP-helix": -Infinity, "MET-sheet/ASP-coil": 1.588, "MET-sheet/ASP-sheet": 2.496, "MET-sheet/SER-helix": -Infinity, "MET-sheet/SER-coil": -2.065, "MET-sheet/SER-sheet": 2.428, "MET-sheet/GLN-helix": -2.58, "MET-sheet/GLN-coil": -1.833, "MET-sheet/GLN-sheet": 1.2, "MET-sheet/LYS-helix": -0.672, "MET-sheet/LYS-coil": -1.868, "MET-sheet/LYS-sheet": 0.674, "MET-sheet/ILE-helix": -3.667, "MET-sheet/ILE-coil": -1.743, "MET-sheet/ILE-sheet": 2.354, "MET-sheet/PRO-helix": -2.516, "MET-sheet/PRO-coil": -1.01, "MET-sheet/PRO-sheet": -1.435, "MET-sheet/THR-helix": -3.071, "MET-sheet/THR-coil": -2.935, "MET-sheet/THR-sheet": -0.173, "MET-sheet/PHE-helix": -Infinity, "MET-sheet/PHE-coil": -2.53, "MET-sheet/PHE-sheet": 1.843, "MET-sheet/ASN-helix": -Infinity, "MET-sheet/ASN-coil": -1.319, "MET-sheet/ASN-sheet": 1.7, "MET-sheet/GLY-helix": -Infinity, "MET-sheet/GLY-coil": -1.7, "MET-sheet/GLY-sheet": 1.91, "MET-sheet/HIS-helix": -Infinity, "MET-sheet/HIS-coil": -Infinity, "MET-sheet/HIS-sheet": 0.045, "MET-sheet/LEU-helix": -2.209, "MET-sheet/LEU-coil": -2.168, "MET-sheet/LEU-sheet": 0.072, "MET-sheet/ARG-helix": -3.053, "MET-sheet/ARG-coil": -1.952, "MET-sheet/ARG-sheet": 2.229, "MET-sheet/TRP-helix": -3.275, "MET-sheet/TRP-coil": -1.665, "MET-sheet/TRP-sheet": -0.012, "MET-sheet/ALA-helix": -2.256, "MET-sheet/ALA-coil": -1.999, "MET-sheet/ALA-sheet": 0.495, "MET-sheet/VAL-helix": -0.196, "MET-sheet/VAL-coil": 0.941, "MET-sheet/VAL-sheet": 1.87, "MET-sheet/GLU-helix": -Infinity, "MET-sheet/GLU-coil": -0.991, "MET-sheet/GLU-sheet": -1.055, "MET-sheet/TYR-helix": -3.727, "MET-sheet/TYR-coil": -3.188, "MET-sheet/TYR-sheet": 1.892, "MET-sheet/MET-helix": -2.703, "MET-sheet/MET-coil": -2.507, "MET-sheet/MET-sheet": 5.581
};

/* Runtime variables and options */
var IS_LOADED = false;
var TIMEOUT = 150;
/* Alignemnt View */
var AV_IS_MERGED = true;
var PROTEIN_LABEL = "";
var REFERENCE_LABEL = "";
var REFERENCE_TRANS_LABEL = "";
var SAMPLE_LABELS = [ ];
var HEATMAP_REFERENCE_DATA = [ ];
var HEATMAP_MERGEDSAMPLE_DATA = [ ];
var HEATMAP_SAMPLE_DATA = [ ];
var MERGED_SAMPLES = { };
var ALIGNMENT_POSITION_MAP = new Object( );
var PROTEIN_POSITION_MAP = new Object( );
var IS_SENSE = true;
var CHAIN_IDENTIFIERS = [ ];
var SELECTED_CHAIN = "";

/* Structure View */
var PROTEIN_VIEWER;
var PROTEIN_VIEWER_SELECTED_RESIDUES = [ ];
var PROTEIN_VIEWER_CONTACT_LABELS = [ ];
var PROTEIN_VIEWER_HIGHLIGHTED_REGION = new Set( );
var PROTEIN_SS_DATA = new Object( );
var PROTEIN_RESIDUE_CONTACTS = new Object( );

/* Info Component */
var INFOHasContactInfo = false;
var INFOInducedSubstitution = false;

/* Mappings */
var INVERT_BASE = {
    "A": "T",
    "a": "t",
    "C": "G",
    "c": "g",
    "G": "C",
    "g": "c",
    "T": "A",
    "t": "a",
    ".": ".",
    ":": ":",
    "-": "-",
    "N": "N"
};
var MOL_STYLES = {
    "default_line": {
            colorfunc: (atom) => {
                return COLOR_SCHEMES.aminoAcids3[ atom.resn ];
            } 
    },
    "default_cartoon": {
        color: "#b8b8b8",
        opacity: 0.4,
        thickness: 0.2,
        arrows: true
    },
    "highlight_cartoon": {
        color: "#ff8400",
        opacity: 0.7,
        thickness: 0.2,
        arrows: true
    }
};
var COLOR_SCHEMES = {
    "nucleotides": {
        "A": "#3E885B",
        "a": "#93CDAA", 
        "C": "#0471A6",
        "c": "#79D0FC",
        "G": "#E8C547",
        "g": "#F3E1A0",
        "T": "#DB5461",
        "t": "#ECA7AE"
    },
    "other": {
        "-": "292926",
        ".": "#b0b0b0",
        ":": "#dbdbdb",
        "___": "#292926",
        "N": "#75551d",
        "STOP": "#8c00ff",
        "coil": "#c9c9c9",
        "helix": "#d13917",
        "sheet": "#1771d1"
    },
    "aminoAcids1": {
        "H": "#87ABFF", // Polar (positive), Basic
        "K": "#4e72c7",
        "R": "#aec3f2",
        "D": "#F75050", // Polar (negative), Acidic
        "E": "#e69797",
        "S": "#AEDEAA", // Polar (neutral)
        "T": "#43c238",
        "N": "#208018",
        "Q": "#61ff54",
        "C": "#1eab13",
        "F": "#FFEAB1", // Aromatic
        "W": "#e6c46a",
        "Y": "#f7ee99",
        "A": "#38ffff", // Aliphatic
        "V": "#a1eded",
        "L": "#679999",
        "I": "#1fadad",
        "M": "#16fac9",
        "P": "#9be8dd",
        "G": "#94fffb"
    },
    "aminoAcids3": {
        "HIS": "#87ABFF", // Polar (positive), Basic
        "LYS": "#4e72c7",
        "ARG": "#aec3f2",
        "ASP": "#F75050", // Polar (negative), Acidic
        "GLU": "#e69797",
        "SER": "#AEDEAA", // Polar (neutral)
        "THR": "#43c238",
        "ASN": "#208018",
        "GLN": "#61ff54",
        "CYS": "#1eab13",
        "PHE": "#FFEAB1", // Aromatic
        "TRP": "#e6c46a",
        "TYR": "#f7ee99",
        "ALA": "#38ffff", // Aliphatic
        "VAL": "#a1eded",
        "LEU": "#679999",
        "ILE": "#1fadad",
        "MET": "#16fac9",
        "PRO": "#9be8dd",
        "GLY": "#94fffb"
    }
};
var STRING_ENCODING = {
    "A": 0,
    "a": 1,
    "C": 2,
    "c": 3,
    "G": 4,
    "g": 5,
    "T": 6,
    "t": 7,
    ".": 8,
    ":": 9,
    "-": 10,
    "___": 11,
    "STOP": 12,
    "N": 13,
    "HIS": 14,
    "LYS": 15,
    "ARG": 16,
    "ASP": 17,
    "GLU": 18,
    "SER": 19,
    "THR": 20,
    "ASN": 21,
    "GLN": 22,
    "CYS": 23,
    "PHE": 24,
    "TRP": 25,
    "TYR": 26,
    "ALA": 27,
    "VAL": 28,
    "LEU": 29,
    "ILE": 30,
    "MET": 31,
    "PRO": 32,
    "GLY": 33,
    "coil": 34,
    "helix": 35,
    "sheet": 36
}
var STRING_DECODING = { };
for ( const [key, value] of Object.keys( STRING_ENCODING ) ) {
    STRING_DECODING[ value ] = key;
}

window.onload = _ => {
    new ResizeObserver( ( entries ) => {
        document.getElementById( "alignmentViewContent" ).style.height = ( VIEW_HEIGHT - entries[ 0 ].contentRect.height ).toString( ) + "px";
        ECHART_AV.resize( );
    } ).observe( document.getElementById( "alignmentViewContent" ) );
    new ResizeObserver( ( entries ) => {
        document.getElementById( "positionSummaryViewContent" ).style.height = ( VIEW_HEIGHT - entries[ 0 ].contentRect.height ).toString( ) + "px";
        ECHART_PS.resize( );
    } ).observe( document.getElementById( "positionSummaryViewContent" ) );
};

/**
 * Initializes drag functionality for the specified element.
 * <p>
 * Source: https://www.w3schools.com/howto/howto_js_draggable.asp
 * 
 * @param {Element} element: The element to enable dragging for.
 */
function dragElement(element) {
  var pos1 = 0, pos2 = 0, pos3 = 0, pos4 = 0;
  if (document.getElementById(element.id + "Header")) {
    // If present, the header is where you move the DIV from:
    document.getElementById(element.id + "Header").onmousedown = dragMouseDown;
  } else {
    // Otherwise, move the DIV from anywhere inside the DIV:
    element.onmousedown = dragMouseDown;
  }

  function dragMouseDown(e) {
    e = e || window.event;
    e.preventDefault();
    // Get the mouse cursor position at startup:
    pos3 = e.clientX;
    pos4 = e.clientY;
    document.onmouseup = closeDragElement;
    // Call a function whenever the cursor moves:
    document.onmousemove = elementDrag;
  }

  function elementDrag(e) {
    e = e || window.event;
    e.preventDefault();
    // Calculate the new cursor position:
    pos1 = pos3 - e.clientX;
    pos2 = pos4 - e.clientY;
    pos3 = e.clientX;
    pos4 = e.clientY;
    // Set the element's new position:
    element.style.top = (element.offsetTop - pos2) + "px";
    element.style.left = (element.offsetLeft - pos1) + "px";
  }

  function closeDragElement() {
    // Stop moving when mouse button is released:
    document.onmouseup = null;
    document.onmousemove = null;
  }
}

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
 * Initializes the protein viewer component.
 * 
 * @param {File} pdbFile: A file object referencing to a .pdb file. 
 */
function initStructureView( pdbFile ) {
    let pdbString;
    let fileReader = new FileReader( );
    fileReader.readAsText( pdbFile );
    fileReader.onload = function(e) {
        pdbString = fileReader.result;
        let element = $('#structureViewContent');
        let config = { backgroundColor: '#FFFFFF', id: 'structureViewCanvas' };
        PROTEIN_VIEWER = $3Dmol.createViewer( element, config );
        PROTEIN_VIEWER.addModel( pdbString, "pdb" );
        PVSetDefaultStyle( );
        PROTEIN_VIEWER.setClickable(
            {},
            true,
            (sel) => {
                PVSelectResidue( sel.resi, sel.chain, sel.resn, true );
                if( ! PROTEIN_VIEWER_HIGHLIGHTED_REGION.has( sel.chain + sel.resi ) ) {
                    ECHART_AV.dispatchAction(
                        {
                            type: 'dataZoom',
                            dataZoomIndex: 0,
                            startValue: Math.max( PROTEIN_POSITION_MAP[ sel.resi ] - 17, 0  ),
                            endValue: Math.min( PROTEIN_POSITION_MAP[ sel.resi ] + 15, Object.keys( DATA[ SELECTED_CHAIN ] ).length )
                        }
                    );
                    PVSetDefaultStyle( );
                    PROTEIN_VIEWER_HIGHLIGHTED_REGION = new Set( );
                    PVHighlightRegion(  Math.max( sel.resi - 5, 0  ), Math.min( sel.resi + 5, Object.keys( DATA[ SELECTED_CHAIN ] ).length ), SELECTED_CHAIN.split( "." )[ 1 ] );
                }
                PSSet( PROTEIN_POSITION_MAP[ sel.resi ] - 1 );
                AVTooltip( PROTEIN_POSITION_MAP[ sel.resi ] - 1, 0, null );
            }
        );
        PROTEIN_VIEWER.zoomTo();
        PROTEIN_VIEWER.render();
        PROTEIN_VIEWER.zoom(0.95, TIMEOUT);
        /* Gather secondary structure information from 3DMol.js GLModel */
        internalState = PROTEIN_VIEWER.getInternalState( );
        for ( const atom of internalState.models[ 0 ].atoms ) {
            switch ( atom.ss ) {
                case 'c':
                    PROTEIN_SS_DATA[ atom.resi.toString( ) + atom.chain ] = 'coil';
                    break;
                case 's':
                    PROTEIN_SS_DATA[ atom.resi.toString( ) + atom.chain ] = 'sheet';
                    break;
                case 'h':
                    PROTEIN_SS_DATA[ atom.resi.toString( ) + atom.chain ] = 'helix';
                    break;
            }
        }
        dragElement(document.getElementById("structureView"));
        showElement( "structureView", "block" );
        document.getElementById( "structureViewCanvas" ).style.top = "24px";
    };
}

/**
 * Initializes the legend component.
 */
function initLegend( ) {
    dragElement(document.getElementById("legend"));
    showElement( "legend", "block" );
}

/**
 * Initializes the position summary component.
 */
function initPositionSummary( ) {
    dragElement(document.getElementById("positionSummaryView"));
    showElement("positionSummaryView", "block");
    var chartDom = document.getElementById('positionSummaryViewContent');
    ECHART_PS = echarts.init(chartDom, { "renderer": "canvas" } );
}

/**
 * Initializes the tooltip component.
 */
 function initTooltip( ) {
    dragElement(document.getElementById("tooltip"));
    showElement("tooltip", "block");
}

/**
 * Initializes alignment visualization.
 * <p>
 * Currently only first chain is used.
 */
function initAlignmentView( ) {
    dragElement( document.getElementById("alignmentView") );
    showElement( "alignmentView", "block" );

    var chartDom = document.getElementById('alignmentViewContent');
    ECHART_AV = echarts.init(chartDom, { "renderer": "canvas" } );
    /* Initialize reference sequence variables. */
    var referenceData = [ ];
    var referenceTransData = [ ];
    /* Initialize protein sequence variables. */
    var proteinData = [ ];
    var ssData = [ ];
    /* Initialize sample variables. */
    var sampleData = { };

    /* Fill sample labels. */
    for (let key of Object.keys( DATA[ "sampleSummary" ] ) ) {
        SAMPLE_LABELS.push( key );
        sampleData[ key ] = [ ];
    }

    /* Parse information from data: */
    var positionCounter = 0;
    for (const [key, value] of Object.entries( DATA[ SELECTED_CHAIN ] ) ) {

        ALIGNMENT_POSITION_MAP[ positionCounter ] = [ key, value[ "proteinPosition" ] ];
        if ( value[ "proteinPosition" ] != "_" ) {
            PROTEIN_POSITION_MAP[ value[ "proteinPosition"] ] = positionCounter;
        }
        if ( positionCounter % 3 == 1 ) {
            proteinData.push( [ positionCounter - 1, 0, value[ "proteinContent" ] ] );
            proteinData.push( [ positionCounter, 0, value[ "proteinContent" ] ] );
            proteinData.push( [ positionCounter + 1, 0, value[ "proteinContent" ] ] );
            if ( value[ "proteinContent" ] != "_" ) {
                ssData.push( [ positionCounter - 1, 1, PROTEIN_SS_DATA[ value[ "proteinPosition" ] + SELECTED_CHAIN.split( "." )[ 1 ] ] ] );
                ssData.push( [ positionCounter, 1, PROTEIN_SS_DATA[ value[ "proteinPosition" ]  + SELECTED_CHAIN.split( "." )[ 1 ] ] ] );
                ssData.push( [ positionCounter + 1, 1, PROTEIN_SS_DATA[ value[ "proteinPosition" ]  + SELECTED_CHAIN.split( "." )[ 1 ] ] ] );
            }
            referenceTransData.push( [ positionCounter - 1, 2, value[ "genomeTranslatedContent" ] ] );
            referenceTransData.push( [ positionCounter, 2, value[ "genomeTranslatedContent" ] ] );
            referenceTransData.push( [ positionCounter + 1, 2, value[ "genomeTranslatedContent" ] ] );
        }
        referenceData.push( [ positionCounter, 3, value[ "genomeContent" ] ] )
        var hasMissense = Object.keys( value[ "missenseMutations" ] ).length != 0;
        if ( hasMissense ){
            for( const sample of SAMPLE_LABELS ) {
                if ( sample in value[ "missenseMutations" ] ) {
                    sampleData[ sample ].push( [ positionCounter, SAMPLE_LABELS.indexOf( sample ), value[ "missenseMutations" ][ sample ].split( ">" )[ 0 ] ] );
                } else {
                    if ( IS_SENSE ) {
                        sampleData[ sample ].push(
                            [ positionCounter, SAMPLE_LABELS.indexOf( sample ), DATA[ "snvs" ][ key ][ sample ] ]
                        );
                    } else {
                        sampleData[ sample ].push(
                            [ positionCounter, SAMPLE_LABELS.indexOf( sample ), INVERT_BASE[ DATA[ "snvs" ][ ( - parseInt( key ) ).toString( ) ][ sample ] ] ]
                        );
                    }
                }
            }
        }
        positionCounter += 1;
    }

    /* Convert gathered data into numeric data for heatmap */
    proteinData.map( entry => {
        HEATMAP_REFERENCE_DATA.push( [ entry[0], entry[1], STRING_ENCODING[ entry[ 2 ] ] ] );
    } );
    ssData.map( entry => {
        HEATMAP_REFERENCE_DATA.push( [ entry[0], entry[1], STRING_ENCODING[ entry[ 2 ] ] ] );
    } );
    referenceTransData.map( entry => {
        HEATMAP_REFERENCE_DATA.push( [ entry[0], entry[1], STRING_ENCODING[ entry[ 2 ] ] ] );
    } );
    referenceData.map( entry => {
        HEATMAP_REFERENCE_DATA.push( [ entry[0], entry[1], STRING_ENCODING[ entry[ 2 ] ] ] );
    } );
    for ( const [key, value] of Object.entries( sampleData ) ) {
        value.map( entry => {
            HEATMAP_SAMPLE_DATA.push( [ entry[0], entry[1], STRING_ENCODING[ entry[ 2 ] ] ] );
        } );
    }

    /* Construct merged sample data */
    for ( const [key, value] of Object.entries( sampleData ) ) {
        let sampleFingerprint = "";
        for ( let v of value ) {
            sampleFingerprint += v[ 2 ];
        }
        if ( sampleFingerprint in MERGED_SAMPLES ) {
            MERGED_SAMPLES[ sampleFingerprint ].push( key );
        } else {
            MERGED_SAMPLES[ sampleFingerprint ] = [ key ];
            value.map( entry => {
                HEATMAP_MERGEDSAMPLE_DATA.push( [ entry[0], Object.keys( MERGED_SAMPLES ).length - 1, STRING_ENCODING[ entry[ 2 ] ] ] );
            } );
        }
    }
    let counter = 1;
    for ( key of Object.keys( MERGED_SAMPLES ) ) {
        let label = "Genotype " + counter + "\n(" + MERGED_SAMPLES[ key ].length;
        if ( MERGED_SAMPLES[ key ].length == 1 ) {
            label += " sample)";
        } else {
            label += " samples)";
        }
        MERGED_SAMPLES[ label ] = MERGED_SAMPLES[ key ];
        delete MERGED_SAMPLES[ key ];
        counter += 1;
    }
    
    /* Add buttons to select chains */
    for ( let chainIdentifier of CHAIN_IDENTIFIERS ) {
        let existingBtn = document.getElementById( "AVSelChainBtn" + chainIdentifier.split( "." )[ 1 ] );
        if ( existingBtn == null ) {
            let btn = document.createElement( "button" );
            btn.innerHTML = chainIdentifier.split( "." )[ 1 ];
            btn.id = "AVSelChainBtn" + chainIdentifier.split( "." )[ 1 ];
            if ( chainIdentifier == SELECTED_CHAIN ) {
                btn.classList.add( "compHeaderBtnActive" );
            } else {
                btn.classList.add( "compHeaderBtn" );
            }
            btn.classList.add( "compHeaderBtn" );
            btn.onclick = function(){
                if ( chainIdentifier != SELECTED_CHAIN ) {
                    SELECTED_CHAIN = chainIdentifier;
                    initAlignmentView( );
                    ECHART_AV.dispatchAction(
                        {
                            type: 'dataZoom',
                            dataZoomIndex: 0,
                            start: 0,
                            end: 100
                        }
                    );
                }
            };
            document.getElementById( "alignmentViewHeader" ).appendChild( btn );
        } else {
            if ( chainIdentifier == SELECTED_CHAIN ) {
                existingBtn.classList.remove( "compHeaderBtn" );
                existingBtn.classList.add( "compHeaderBtnActive" );
            } else {
                existingBtn.classList.remove( "compHeaderBtnActive" );
                existingBtn.classList.add( "compHeaderBtn" );
            }
        }
    }

    ECHART_AV_OPTION = {
        animation: false,
        grid: [
            {   // Grid to display reference alignment information.
                top: '2%',
                bottom: '80%',
                left: '7%',
                right: '4%',
                show: true,
                backgroundColor: "#292926"
            },
            {   // Grid to display sample alignment information.
                top: '21%',
                bottom: '10%',
                left: '7%',
                right: '4%'
            }
        ],
        xAxis: [
            {   // X axis to display reference alignment information.
                type: 'category',
                data: Object.keys( ALIGNMENT_POSITION_MAP ),
                splitArea: {
                    show: false
                },
                show: false,
                gridIndex: 0
            },
            {   // X axis to display sample alignment information.
                type: 'category',
                data: Object.keys( ALIGNMENT_POSITION_MAP ),
                splitArea: {
                    show: false
                },
                show: true,
                gridIndex: 1
            }
        ],
        yAxis: [
            {   // Y axis to display reference alignment information.
                type: 'category',
                data: [ PROTEIN_LABEL, "Sec. Struc.  ", REFERENCE_TRANS_LABEL, REFERENCE_LABEL ],
                splitArea: {
                    show: false
                },
                gridIndex: 0,
                inverse: true
            },
            {   // Y axis to display sample alignment information.
                type: 'category',
                data: Object.keys( MERGED_SAMPLES ),
                splitArea: {
                    show: false
                },
                gridIndex: 1,
                inverse: true
            }
        ],
        legend: { // Disable default legend.
            show: false
        },
        tooltip: { // Disable default tooltip.
            show: true,
            trigger: 'item',
            triggerOn: 'click',
            formatter: function (params, ticket, callback) {
                /* Collect information about position: */
                AVTooltip( params.data[ 0 ], params.seriesIndex, params.data[ 1 ] );
            }
        },
        dataZoom: [
            {   // Data zoom to zoom into position intervals.
                id: "positionZoom",
                type: 'slider',
                xAxisIndex: [ 0, 1 ],
                realtime: false,
                throttle: 100,
                height: "10px",
                bottom: "2%",
                id: "positionZoom",
                showDataShadow: 'auto',
                backgroundColor: "transparent",
                fillerColor: "transparent",
                borderColor: "rgba(228, 229, 237, 0.733)",
                handleStyle: {
                    borderColor: "rgba(96, 113, 150, 1.0)"
                },
                moveHandleStyle: {
                    color: "rgba(96, 113, 150, 0.533)"
                },
                emphasis: {
                    moveHandleStyle: {
                        color: "rgba(96, 113, 150, 1.0)"
                    }
                }
            },
            {   // Data zoom to zoom into single samples.
                type: 'slider',
                yAxisIndex: [ 1 ],
                throttle: 5,
                width: "10px",
                right: '2%',
                backgroundColor: "transparent",
                fillerColor: "transparent",
                borderColor: "rgba(228, 229, 237, 0.733)",
                handleStyle: {
                    borderColor: "rgba(96, 113, 150, 1.0)"
                },
                moveHandleStyle: {
                    color: "rgba(96, 113, 150, 0.533)"
                },
                emphasis: {
                    moveHandleStyle: {
                        color: "rgba(96, 113, 150, 1.0)"
                    }
                }
            }
        ],
        visualMap: [
            {   // Visual map to color heatmap.
                type: 'piecewise',
                pieces: [
                    // ACGT$.:;-aqygfhtrzcvd
                    { min: STRING_ENCODING.A, max: STRING_ENCODING.A, color: COLOR_SCHEMES.nucleotides.A }, // A
                    { min: STRING_ENCODING.C, max: STRING_ENCODING.C, color: COLOR_SCHEMES.nucleotides.C }, // C
                    { min: STRING_ENCODING.G, max: STRING_ENCODING.G, color: COLOR_SCHEMES.nucleotides.G }, // G
                    { min: STRING_ENCODING.T, max: STRING_ENCODING.T, color: COLOR_SCHEMES.nucleotides.T }, // T
                    { min: STRING_ENCODING.a, max: STRING_ENCODING.a, color: COLOR_SCHEMES.nucleotides.a }, // A, declined
                    { min: STRING_ENCODING.c, max: STRING_ENCODING.c, color: COLOR_SCHEMES.nucleotides.c }, // C, declined
                    { min: STRING_ENCODING.g, max: STRING_ENCODING.g, color: COLOR_SCHEMES.nucleotides.g }, // G, declined
                    { min: STRING_ENCODING.t, max: STRING_ENCODING.t, color: COLOR_SCHEMES.nucleotides.t }, // T, declined
                    { min: STRING_ENCODING.___, max: STRING_ENCODING.___, color: COLOR_SCHEMES.other.___ }, // No match with protein sequence
                    { min: STRING_ENCODING["."], max: STRING_ENCODING["."], color: COLOR_SCHEMES.other["."] }, // Reference
                    { min: STRING_ENCODING[":"], max: STRING_ENCODING[":"], color: COLOR_SCHEMES.other[":"] }, // Reference, declined
                    { min: STRING_ENCODING.STOP, max: STRING_ENCODING.STOP, color: COLOR_SCHEMES.other.STOP }, // Translated stop codon
                    { min: STRING_ENCODING["-"], max: STRING_ENCODING["-"], color: COLOR_SCHEMES.other["-"] }, // Deletion
                    { min: STRING_ENCODING.N, max: STRING_ENCODING.N, color: COLOR_SCHEMES.other.N }, // No call
                    { min: STRING_ENCODING.ALA, max: STRING_ENCODING.ALA, color: COLOR_SCHEMES.aminoAcids3.ALA }, // Amino acids
                    { min: STRING_ENCODING.ARG, max: STRING_ENCODING.ARG, color: COLOR_SCHEMES.aminoAcids3.ARG },
                    { min: STRING_ENCODING.ASP, max: STRING_ENCODING.ASP, color: COLOR_SCHEMES.aminoAcids3.ASP },
                    { min: STRING_ENCODING.ASN, max: STRING_ENCODING.ASN, color: COLOR_SCHEMES.aminoAcids3.ASN },
                    { min: STRING_ENCODING.GLN, max: STRING_ENCODING.GLN, color: COLOR_SCHEMES.aminoAcids3.GLN },
                    { min: STRING_ENCODING.GLU, max: STRING_ENCODING.GLU, color: COLOR_SCHEMES.aminoAcids3.GLU },
                    { min: STRING_ENCODING.GLY, max: STRING_ENCODING.GLY, color: COLOR_SCHEMES.aminoAcids3.GLY },
                    { min: STRING_ENCODING.ILE, max: STRING_ENCODING.ILE, color: COLOR_SCHEMES.aminoAcids3.ILE },
                    { min: STRING_ENCODING.LEU, max: STRING_ENCODING.LEU, color: COLOR_SCHEMES.aminoAcids3.LEU },
                    { min: STRING_ENCODING.THR, max: STRING_ENCODING.THR, color: COLOR_SCHEMES.aminoAcids3.THR },
                    { min: STRING_ENCODING.TRP, max: STRING_ENCODING.TRP, color: COLOR_SCHEMES.aminoAcids3.TRP },
                    { min: STRING_ENCODING.TYR, max: STRING_ENCODING.TYR, color: COLOR_SCHEMES.aminoAcids3.TYR },
                    { min: STRING_ENCODING.PHE, max: STRING_ENCODING.PHE, color: COLOR_SCHEMES.aminoAcids3.PHE },
                    { min: STRING_ENCODING.PRO, max: STRING_ENCODING.PRO, color: COLOR_SCHEMES.aminoAcids3.PRO },
                    { min: STRING_ENCODING.CYS, max: STRING_ENCODING.CYS, color: COLOR_SCHEMES.aminoAcids3.CYS },
                    { min: STRING_ENCODING.SER, max: STRING_ENCODING.SER, color: COLOR_SCHEMES.aminoAcids3.SER },
                    { min: STRING_ENCODING.MET, max: STRING_ENCODING.MET, color: COLOR_SCHEMES.aminoAcids3.MET },
                    { min: STRING_ENCODING.VAL, max: STRING_ENCODING.VAL, color: COLOR_SCHEMES.aminoAcids3.VAL },
                    { min: STRING_ENCODING.HIS, max: STRING_ENCODING.HIS, color: COLOR_SCHEMES.aminoAcids3.HIS },
                    { min: STRING_ENCODING.LYS, max: STRING_ENCODING.LYS, color: COLOR_SCHEMES.aminoAcids3.LYS },
                    { min: STRING_ENCODING.coil, max: STRING_ENCODING.coil, color: COLOR_SCHEMES.other.coil },
                    { min: STRING_ENCODING.sheet, max: STRING_ENCODING.sheet, color: COLOR_SCHEMES.other.sheet },
                    { min: STRING_ENCODING.helix, max: STRING_ENCODING.helix, color: COLOR_SCHEMES.other.helix }
                ],
                seriesIndex: [ 0, 1 ],
                show: false
            }
        ],
        series: [
            {
                name: 'ReferenceHeatmap',
                type: 'heatmap',
                progressive: 500,
                data: HEATMAP_REFERENCE_DATA,
                xAxisIndex: 0,
                yAxisIndex: 0,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.2
                }
            },
            {
                name: 'SampleHeatmap',
                type: 'heatmap',
                progressive: 500,
                data: HEATMAP_MERGEDSAMPLE_DATA,
                xAxisIndex: 1,
                yAxisIndex: 1,
                itemStyle: {
                    borderColor: "#fafafc",
                    borderWidth: 0.2
                }
            }
        ]
    };

    ECHART_AV.on('datazoom', function(params) {
        let proteinFirstIndex = Math.min( ... Object.keys( PROTEIN_POSITION_MAP ).map( e => parseInt( e ) ) );
        let proteinLastIndex = Math.max( ... Object.keys( PROTEIN_POSITION_MAP ).map( e => parseInt( e ) ) );
        if ( params.dataZoomId == "positionZoom" ) {
            let fromIsStart = false;
            let toIsEnd = false;
            let length = Object.keys( ALIGNMENT_POSITION_MAP ).length;
            let fromIndex = Math.max( 0, Math.floor( ( params.start / 100 ) * length ) );
            let from = ALIGNMENT_POSITION_MAP[ fromIndex ][ 1 ];
            if ( from == "_" || from <= proteinFirstIndex ) {
                fromIsStart = true;
                from = proteinFirstIndex;
            }
            let toIndex = Math.min( length - 1, Math.floor( ( params.end / 100 ) * length ) );
            let to = ALIGNMENT_POSITION_MAP[ toIndex ][ 1 ];
            if ( to == "_" || to >= proteinLastIndex ) {
                toIsEnd = true;
                to = proteinLastIndex;
            }
            /*
            if ( to == -1 ) {
                while ( to == -1 ) {
                    toIndex -= 1;
                    to = ALIGNMENT_POSITION_MAP[ toIndex ][ 1 ];
                }
            }
            */
            PVSetDefaultStyle( );
            PROTEIN_VIEWER_HIGHLIGHTED_REGION = new Set( );
            if ( ! fromIsStart || ! toIsEnd ) {
                PVHighlightRegion( from, to, SELECTED_CHAIN.split( "." )[ 1 ] );
            }
        }
    });

    ECHART_AV.setOption(ECHART_AV_OPTION);
}

/**
 * Initializes all visualizations from the passed file list and data object.
 * 
 * @param {FileList} files: List of file references to MUSIAL output files.
 */
function init( files ) {
    // Close input component.
    closeElement( "input" );
    // Set initial chain to be the first one parsed.
    SELECTED_CHAIN = CHAIN_IDENTIFIERS[ 0 ];
    // Initialize structure view component.
    var pdbFile;
    for ( let file of files ) {
        if ( file.name.endsWith(".pdb") ) {
            pdbFile = file;
        }
    }
    window.setTimeout( function(){ initStructureView( pdbFile ); }, 100 );
    // Initialize alignment view component.
    window.setTimeout( function(){
        initAlignmentView( );
    }, 150 );
    // Initialize legend component.
    window.setTimeout( function(){ initLegend( ); }, 180 );
    // Initialize the position summary component.
    window.setTimeout( function(){ initPositionSummary( ); }, 210 );
    // Initialize the tooltip component.
    window.setTimeout( function(){ initTooltip( ); }, 240 );
    IS_LOADED = true;
}

/**
 * Handles the user selection of a directory.
 * 
 * @param {Event} event: A filedrop event. The target of the event has to be a `FileList` instance. 
 */
function directoryInputHandler( event ) {
    document.getElementById( "directoryInput" ).style.display = "none";
    document.getElementById( "directoryInputIcon" ).style.display = "none";
    document.getElementById( "directoryInputText" ).style.display = "none";
    document.getElementById( "directoryInputProcessLoader" ).style.display = "flex";
    var files = event.target.files;
    parseResults( files );
    window.setTimeout( function(){ init( files ); }, 1000 );
}

/**
 * Parses the MUSIAL output files specified with the passed file list.
 * 
 * @param {FileList} fileList: A `FileList` instance that specifies the output files of a MUSIAL run. 
 */
function parseResults( fileList ) {
    // Parse data from protein alignment.
    var file;
    for ( file of fileList ) {
        if ( file.name.startsWith( "proteinAlignment" ) ) {
            var chainId = "chain." + file.name.split( "_" )[ 1 ].split( "." )[ 0 ];
            CHAIN_IDENTIFIERS.push( chainId );
            DATA[ chainId ] = new Object( );
            parseProteinAlignment( file, chainId );
        }
    }

    // Parse data from protein residue contact map.
    file = [...fileList].filter( s => s.name.startsWith("residueContactMap") )[ 0 ];
    parseContactMap( file );
    
    // Parse per position summary.
    file = [...fileList].filter( s => s.name.startsWith("positionStatistics") )[ 0 ];
    DATA[ "positionSummary" ] = new Object( );
    parsePositionSummary( file );
    
    // Parse per sample summary.
    file = [...fileList].filter( s => s.name.startsWith("sampleStatistics") )[ 0 ];
    DATA[ "sampleSummary" ] = new Object( );
    parseSampleSummary( file );

    // Parse SNV table annotations.
    file = [...fileList].filter( s => s.name.startsWith("snvAnnotations") )[ 0 ];
    DATA[ "annotations" ] = new Object( );
    parseAnnotations( file );

    // Parse SNV table annotations.
    file = [...fileList].filter( s => s.name.startsWith("snvTable") )[ 0 ];
    DATA[ "snvs" ] = new Object( );
    parseSNVs( file );
}

function parseProteinAlignment( file, chainId ) {
    var reader = new FileReader( );
    var fileString;
    reader.onload = function(event) {
        fileString = event.target.result;
        var results = fileString.split( "\r\n" );
        for ( let row of results ) {
            var rowData = row.split( "\t" );
            if ( rowData[ 0 ] == "TRC" ) {
                var genomePositions = rowData[ 1 ].split( "," );
                var genomeContent = rowData[ 2 ].split( "" );
                var proteinPosition = "-1";
                var proteinContent = "___";
                for (var index = 0; index < genomePositions.length; index++) {
                    DATA[ chainId ][ genomePositions[ index ] ] = {
                        "genomeContent": genomeContent[ index ],
                        "genomeTranslatedContent": "_",
                        "missenseMutations": { },
                        "proteinPosition": proteinPosition,
                        "proteinContent": proteinContent
                    }
                }
            } else if ( rowData[ 0 ] == "ALN" ) {
                var genomePositions = rowData[ 1 ].split( "," );
                var genomeContent = rowData[ 2 ].split( "" );
                var genomeTranslatedContent = rowData[ 3 ];
                var proteinPositions = rowData[ 4 ];
                var proteinContent = rowData[ 5 ];
                var missenseMutations = [ { }, { }, { } ];
                // Check if feature is located on anti-sense strand:
                if ( genomePositions[ 0 ] < 0 ) {
                    IS_SENSE = false;
                }
                if ( rowData[ 6 ] != "" ) {
                    var missenseMutationsList = rowData[ 6 ].split( "," );
                    for ( var missenseMutation of missenseMutationsList ) {
                        var missenseMutationFields = missenseMutation.split( "|" );
                        var sampleName = missenseMutationFields[ 0 ];
                        var sampleContent = missenseMutationFields[ 1 ].split( "" );
                        var sampleTranslatedContent = missenseMutationFields[ 2 ];
                        for ( var index = 0; index < sampleContent.length; index++ ) {
                            if ( sampleContent[ index ] != genomeContent[ index ] ) {
                                missenseMutations[ index ][ sampleName ] = sampleContent[ index ] + ">" + sampleTranslatedContent;
                            }
                        }
                    }
                }
                for (var index = 0; index < genomePositions.length; index++) {
                    DATA[ chainId ][ genomePositions[ index ] ] = {
                        "genomeContent": genomeContent[ index ],
                        "genomeTranslatedContent": genomeTranslatedContent,
                        "missenseMutations": missenseMutations[ index ],
                        "proteinPosition": ( proteinContent == "@" ) ? "-1" : proteinPositions,
                        "proteinContent": ( proteinContent == "@" ) ? "___" : proteinContent
                    }
                }
            } else if ( rowData[ 0 ] == "#input" ) {
                PROTEIN_LABEL = rowData[ 1 ].split( "=" )[ 1 ];
                REFERENCE_LABEL = rowData[ 2 ].split( "=" )[ 1 ]; 
                REFERENCE_TRANS_LABEL = REFERENCE_LABEL + ", trans.";
            }
        }
    };
    reader.readAsText( file );
}

function parseContactMap( file ) {
    var reader = new FileReader( );
    var fileString;
    reader.onload = function(event) {
        fileString = event.target.result;
        var results = fileString.split( "\r\n" );
        var residues = results.shift( ).split( "\t" );
        residues.shift( );
        for ( var residue of residues ) {
            PROTEIN_RESIDUE_CONTACTS[ residue ] = [ ];
        }
        for ( let r of results ) {
            let row = r.split( "\t" );
            if ( row.length > 1 ) {
                var residue = row[ 0 ];
                row.shift( );
                for ( var index = 0; index < row.length; index++ ) {
                    if ( row[ index ] <= 5.0 && residues[ index ] !== residue ) {
                        PROTEIN_RESIDUE_CONTACTS[ residues[ index ] ].push( [ residue, row[ index ] ] );
                    }
                }
            }
        }
    };
    reader.readAsText( file );
}

function parsePositionSummary( file ) {
    var reader = new FileReader( );
    var fileString;
    reader.onload = function(event) {
        fileString = event.target.result;
        var results = fileString.split( "\r\n" );
        var header = results.shift( ).split( "\t" );
        header.shift( );
        for ( let r of results ) {
            let row = r.split( "\t" );
            if ( row.length > 1 ) {
                var position = row[ 0 ];
                DATA[ "positionSummary" ][ position ] = { };
                row.shift( );
                for ( var index = 0; index < row.length; index++ ) {
                    DATA[ "positionSummary" ][ position ][ header[ index ] ] = row[ index ];
                }
            }
        }
    };
    reader.readAsText( file );
}

function parseSampleSummary( file ) {
    var reader = new FileReader( );
    var fileString;
    reader.onload = function(event) {
        fileString = event.target.result;
        var results = fileString.split( "\r\n" );
        var header = results.shift( ).split( "\t" );
        header.shift( );
        for ( let r of results ) {
            let row = r.split( "\t" );
            if ( row.length > 1 ) {
                var sample = row[ 0 ];
                DATA[ "sampleSummary" ][ sample ] = { };
                row.shift( );
                for ( var index = 0; index < row.length; index++ ) {
                    DATA[ "sampleSummary" ][ sample ][ header[ index ] ] = row[ index ];
                }
            }
        }
    };
    reader.readAsText( file );
}

function parseAnnotations( file ) {
    var reader = new FileReader( );
    var fileString;
    reader.onload = function(event) {
        fileString = event.target.result;
        var results = fileString.split( "\r\n" );
        var header = results.shift( ).split( "\t" );
        header.shift( );
        header.shift( );
        for ( let r of results ) {
            let row = r.split( "\t" );
            if ( row.length > 1 ) {
                var position = row[ 0 ];
                DATA[ "annotations" ][ position ] = { }; 
                row.shift( );
                row.shift( );
                for ( var index = 0; index < row.length; index++ ) {
                    DATA[ "annotations" ][ position ][ header[ index ] ] = row[ index ];
                }
            }
        }
    };
    reader.readAsText( file );
}

function parseSNVs( file ) {
    var reader = new FileReader( );
    var fileString;
    reader.onload = function(event) {
        fileString = event.target.result;
        var results = fileString.split( "\r\n" );
        var header = results.shift( ).split( "\t" );
        header.shift( );
        header.shift( );
        for ( let r of results ) {
            let row = r.split( "\t" );
            if ( row.length > 1 ) {
                var position = row[ 0 ];
                DATA[ "snvs" ][ position ] = { };
                row.shift( );
                row.shift( );
                for ( var index = 0; index < row.length; index++ ) {
                    DATA[ "snvs" ][ position ][ header[ index ] ] = row[ index ];
                }
            }
        }
    };
    reader.readAsText( file );
}

function AVUngroup( ) {
    if ( IS_LOADED ) {
        ECHART_AV_OPTION = {
            yAxis: [
                {   // Y axis to display reference alignment information.
                    type: 'category',
                    data: [ PROTEIN_LABEL, "", REFERENCE_TRANS_LABEL, REFERENCE_LABEL ],
                    splitArea: {
                        show: false
                    },
                    gridIndex: 0,
                    inverse: true
                },
                {   // Y axis to display sample alignment information.
                    type: 'category',
                    data: SAMPLE_LABELS,
                    splitArea: {
                        show: false
                    },
                    gridIndex: 1,
                    inverse: true
                }
            ],
            series: [
                {
                    name: 'ReferenceHeatmap',
                    type: 'heatmap',
                    data: HEATMAP_REFERENCE_DATA,
                    /*
                    label: {
                        show: false,
                        formatter: function ( params ) {
                            return dnaIntToChar( params.value[ 2 ] );
                        },
                        color: "black"
                    },
                    */
                    emphasis: {
                        itemStyle: {
                            shadowBlur: 10,
                            shadowColor: 'rgba(0, 0, 0, 0.5)'
                        }
                    },
                    xAxisIndex: 0,
                    yAxisIndex: 0,
                    itemStyle: {
                        borderColor: "#E8E9ED",
                        borderWidth: 0.1
                    }
                },
                {
                    name: 'SampleHeatmap',
                    type: 'heatmap',
                    data: HEATMAP_SAMPLE_DATA,
                    emphasis: {
                        itemStyle: {
                            shadowBlur: 10,
                            shadowColor: 'rgba(0, 0, 0, 0.5)'
                        }
                    },
                    xAxisIndex: 1,
                    yAxisIndex: 1,
                    itemStyle: {
                        borderColor: "#E8E9ED",
                        borderWidth: 0.1
                    }
                }
            ]
        };
        ECHART_AV.setOption(ECHART_AV_OPTION);
        AV_IS_MERGED = false;
        document.getElementById("AVGroupingBtn").innerHTML = '<i class="fas fa-object-group"></i>';
        document.getElementById("AVGroupingBtn").onclick = AVGroup;
    }
}

function AVGroup( ) {
    if ( IS_LOADED ) {
        ECHART_AV_OPTION = {
            yAxis: [
                {   // Y axis to display reference alignment information.
                    type: 'category',
                    data: [ PROTEIN_LABEL, "", REFERENCE_TRANS_LABEL, REFERENCE_LABEL ],
                    splitArea: {
                        show: false
                    },
                    gridIndex: 0,
                    inverse: true
                },
                {   // Y axis to display sample alignment information.
                    type: 'category',
                    data: Object.keys( MERGED_SAMPLES ),
                    splitArea: {
                        show: false
                    },
                    gridIndex: 1,
                    inverse: true
                }
            ],
            series: [
                {
                    name: 'ReferenceHeatmap',
                    type: 'heatmap',
                    data: HEATMAP_REFERENCE_DATA,
                    /*
                    label: {
                        show: false,
                        formatter: function ( params ) {
                            return dnaIntToChar( params.value[ 2 ] );
                        },
                        color: "black"
                    },
                    */
                    emphasis: {
                        itemStyle: {
                            shadowBlur: 10,
                            shadowColor: 'rgba(0, 0, 0, 0.5)'
                        }
                    },
                    xAxisIndex: 0,
                    yAxisIndex: 0,
                    itemStyle: {
                        borderColor: "#E8E9ED",
                        borderWidth: 0.1
                    }
                },
                {
                    name: 'SampleHeatmap',
                    type: 'heatmap',
                    data: HEATMAP_MERGEDSAMPLE_DATA,
                    emphasis: {
                        itemStyle: {
                            shadowBlur: 10,
                            shadowColor: 'rgba(0, 0, 0, 0.5)'
                        }
                    },
                    xAxisIndex: 1,
                    yAxisIndex: 1,
                    itemStyle: {
                        borderColor: "#E8E9ED",
                        borderWidth: 0.1
                    }
                }
            ]
        };
        ECHART_AV.setOption(ECHART_AV_OPTION);
        AV_IS_MERGED = true;
        document.getElementById("AVGroupingBtn").innerHTML = '<i class="fas fa-object-ungroup"></i>';
        document.getElementById("AVGroupingBtn").onclick = AVUngroup;
    }
}

function AVTooltip( position, seriesIndex, yValue ) {
    let tooltipHtmlString = "";
    let genomePosition = ALIGNMENT_POSITION_MAP[ position ][ 0 ];
    /* Based on position collect information about reference and sample: */
    let referenceContent = "";
    switch ( position % 3 ) {
        case 0:
            referenceContent = "<u>" + DATA[ SELECTED_CHAIN ][ genomePosition ].genomeContent + "</u>"
            + DATA[ SELECTED_CHAIN ][ parseInt( genomePosition ) + 1 ].genomeContent
            + DATA[ SELECTED_CHAIN ][ parseInt( genomePosition ) + 2 ].genomeContent;
            break;
        case 1:
            referenceContent = DATA[ SELECTED_CHAIN ][ String( parseInt( genomePosition ) - 1 )].genomeContent
            + "<u>" + DATA[ SELECTED_CHAIN ][ genomePosition ].genomeContent + "</u>"
            + DATA[ SELECTED_CHAIN ][ String( parseInt( genomePosition ) + 1 ) ].genomeContent;
            break;
        case 2:
            referenceContent = DATA[ SELECTED_CHAIN ][ String( parseInt( genomePosition ) - 2 ) ].genomeContent
            + DATA[ SELECTED_CHAIN ][ String( parseInt( genomePosition ) - 1 ) ].genomeContent
            + "<u>" + DATA[ SELECTED_CHAIN ][ genomePosition ].genomeContent + "</u>";
            break;
    }
    if ( referenceContent.includes( "___" ) ) {
        referenceContent = "None"
    }
    let referenceTransContent = DATA[ SELECTED_CHAIN ][ genomePosition ].genomeTranslatedContent;
    if ( referenceTransContent == "___" || referenceTransContent == "_" ) {
        referenceTransContent = "None"
    } else if ( referenceTransContent == "STOP" ) {
        referenceTransContent = "None (Termination)"
    }
    let proteinSequenceContent = DATA[ SELECTED_CHAIN ][ genomePosition ].proteinContent;
    if ( proteinSequenceContent == "___" || proteinSequenceContent == "_" ) {
        proteinSequenceContent = "None";
    } else if ( proteinSequenceContent == "STOP" ) {
        proteinSequenceContent = "None (Termination)";
    }
    let sampleName = false;
    let sampleContent = false;
    let codonPositions = [ ];
    let inducedSubstitution = false;
    let sampleAnnotation = false;
    /* If series index is not zero, sample information is available: */
    if ( seriesIndex != 0 ) {
        let assignContent = function( pos, sId ) {
            let content;
            if ( IS_SENSE ) {
                if ( pos in DATA[ "snvs" ] ) {
                    content = DATA[ "snvs" ][ pos ][ sId ];
                } else {
                    content = ".";
                }
            } else {
                if ( ( - parseInt( pos ) ) in DATA[ "snvs" ] ) {
                    content = INVERT_BASE[ DATA[ "snvs" ][ ( - parseInt( pos ) ) ][ sId ] ];
                } else {
                    content = ".";
                }
            }
            if ( content.includes( "." ) || content.includes( ":" ) ) {
                return DATA[ SELECTED_CHAIN ][ pos ].genomeContent;
            } else {
                return content;
            }
        }
        let sampleAccessor;
        inducedSubstitution = "None";
        if ( AV_IS_MERGED ) {
            sampleName = Object.keys( MERGED_SAMPLES )[ yValue ];
            sampleAccessor = MERGED_SAMPLES[ sampleName ][ 0 ];
        } else {
            sampleName = SAMPLE_LABELS[ yValue ];
            sampleAccessor = sampleName;
            if ( IS_SENSE ) {
                sampleAnnotation = DATA[ "annotations" ][ genomePosition ][ sampleName ];
            } else {
                sampleAnnotation = DATA[ "annotations" ][ ( - parseInt( genomePosition ) ).toString( ) ][ sampleName ];
            }
        }
        switch ( position % 3 ) {
            case 0:
                sampleContent = "<u>" + assignContent( genomePosition, sampleAccessor ) + "</u>"
                + assignContent( parseInt( genomePosition ) + 1, sampleAccessor )
                + assignContent( parseInt( genomePosition ) + 2, sampleAccessor );
                codonPositions = [ genomePosition, parseInt( genomePosition ) + 1, parseInt( genomePosition ) + 2 ];
                break;
            case 1:
                sampleContent = assignContent( parseInt( genomePosition ) - 1, sampleAccessor )
                + "<u>" + assignContent( genomePosition, sampleAccessor ) + "</u>"
                + assignContent( parseInt( genomePosition ) + 1, sampleAccessor );
                codonPositions = [ parseInt( genomePosition ) - 1, genomePosition, parseInt( genomePosition ) + 1 ];
                break;
            case 2:
                sampleContent = assignContent( parseInt( genomePosition ) - 2, sampleAccessor )
                + assignContent( parseInt( genomePosition ) - 1, sampleAccessor )
                + "<u>" + assignContent( genomePosition, sampleAccessor ) + "</u>";
                codonPositions = [ parseInt( genomePosition ) - 2, parseInt( genomePosition ) - 1, genomePosition ];
                break;
        }
        // Check for induced substitutions.
        for ( let codonPosition of codonPositions ) {
            if ( sampleAccessor in DATA[ SELECTED_CHAIN ][ codonPosition ][ "missenseMutations" ] && inducedSubstitution == "None" ) {
                inducedSubstitution = DATA[ SELECTED_CHAIN ][ codonPosition ][ "missenseMutations" ][ sampleAccessor ].split( ">" )[ 1 ];
            }
        }
    }
    tooltipHtmlString += "<b>Reference Information</b><br>"
    tooltipHtmlString += "Position: " + genomePosition + "<br>";
    tooltipHtmlString += "Nucleotide codon: " + referenceContent + "<br>";
    tooltipHtmlString += "Translated amino acid: " + referenceTransContent + "<br>";
    tooltipHtmlString += "Aligned amino acid: " + proteinSequenceContent + "<br>";
    if ( sampleName ) {
        tooltipHtmlString += "<hr>";
        tooltipHtmlString += "<b>Sample Information</b><br>"
        tooltipHtmlString += "Name: " + sampleName + "<br>";
        tooltipHtmlString += "Variant codon: " + sampleContent + "<br>";
        tooltipHtmlString += "Amino acid substitution: " + inducedSubstitution + "<br>";
        if ( inducedSubstitution != "None" ) {
            INFOInducedSubstitution = inducedSubstitution;
        } else {
            INFOInducedSubstitution = false;
        }
        if ( sampleAnnotation ) {
            tooltipHtmlString += "<hr>"
            tooltipHtmlString += "<b>Sample Annotation</b><br>"
            let coverage = "";
            let quality = "";
            let frequency = "";
            for ( let annotationField of sampleAnnotation.split( "@" ) ) {
                if ( annotationField == "LOW_COVERAGE" ) {
                    coverage = '<i style="color: #D56062;" class="fas fa-angle-double-down"></i> ';
                } else if ( annotationField == "LOW_QUALITY" ) {
                    quality = '<i style="color: #D56062;" class="fas fa-angle-double-down"></i> ';
                } else if ( annotationField == "LOW_FREQUENCY" ) {
                    frequency = '<i style="color: #D56062;" class="fas fa-angle-double-down"></i> ';
                } else {
                    let field = annotationField.split( "=" )[ 0 ];
                    let value = annotationField.split( "=" )[ 1 ];
                    if ( field == "COVERAGE" ) {
                        coverage += "Coverage: " + parseFloat( value ).toFixed( 2 ) + "<br>";
                        tooltipHtmlString += coverage;
                    } else if ( field == "QUALITY" ) {
                        quality += "Quality: " + parseFloat( value ).toFixed( 2 ) + "<br>";
                        tooltipHtmlString += quality;
                    } else if ( field == "FREQUENCY" ) {
                        frequency += "Frequency: " + parseFloat( value ).toFixed( 2 ) + "<br>";
                        tooltipHtmlString += frequency;
                    }
                }
            }
        }
    }
    /* Invoke position summary tooltip with timeout (i.e. to implement delay) */
    PSSet( position );
    PVSelectResidue( ALIGNMENT_POSITION_MAP[ position ][ 1 ], SELECTED_CHAIN.split( "." )[ 1 ], proteinSequenceContent, true );
    INFOSetContent( tooltipHtmlString );
    PVShowContacts( );
    return;
}

function INFOSetContent( htmlString ) {
    INFOHasContactInfo = false;
    document.getElementById( "tooltipContent" ).innerHTML = htmlString;
}

function INFOAddContactInfo( htmlString ) {
    if ( ! INFOHasContactInfo ) {
        document.getElementById( "tooltipContent" ).innerHTML += htmlString;
    }
    INFOHasContactInfo = true;
}

function INFOGetContent( ) {
    return document.getElementById( "tooltipContent" ).innerHTML;
}

function PSSet( position ) {
    /* Collect count summary wrt. nucleotides of position */
    var outerPieData = [ ];
    let genomePosition;
    if ( IS_SENSE ) {
        genomePosition = ALIGNMENT_POSITION_MAP[ position ][ 0 ];
    } else {
        genomePosition = - parseInt( ALIGNMENT_POSITION_MAP[ position ][ 0 ] );
    }
    if ( genomePosition in DATA[ "snvs" ] ) {
        let counts = {
            "A": 0,
            "A, dec.": 0,
            "C": 0,
            "C, dec.": 0,
            "G": 0,
            "G, dec.": 0,
            "T": 0,
            "T, dec.": 0,
            "Ref.": 0,
            "Ref., dec.": 0,
            "No Call": 0,
            "Del.": 0,
        }
        for ( const sample of SAMPLE_LABELS ) {
            let sampleContent;
            if ( IS_SENSE ) {
                sampleContent = DATA[ "snvs" ][ genomePosition ][ sample ];
            } else {
                sampleContent = INVERT_BASE[ DATA[ "snvs" ][ genomePosition ][ sample ] ];
            }
            switch (sampleContent) {
                case "A":
                    counts[ "A" ] += 1;
                    break;
                case "a":
                    counts[ "A, dec." ] += 1;
                    break;
                case "C":
                    counts[ "C" ] += 1;
                    break;
                case "c":
                    counts[ "C, dec." ] += 1;
                    break;
                case "G":
                    counts[ "G" ] += 1;
                    break;
                case "g":
                    counts[ "G, dec." ] += 1;
                    break;
                case "T":
                    counts[ "T" ] += 1;
                    break;
                case "t":
                    counts[ "T, dec." ] += 1;
                    break;
                case ".":
                    counts[ "Ref." ] += 1;
                    break;
                case ":":
                    counts[ "Ref., dec." ] += 1;
                    break;
                case "N":
                    counts[ "No Call" ] += 1;
                    break;
                case "-":
                    counts[ "Del." ] += 1;
                    break;
            }
        }
        for ( const [key, value] of Object.entries( counts ) ) {
            if ( value > 0 ) {
                let color;
                switch (key) {
                    case "A":
                    color = COLOR_SCHEMES.nucleotides.A;
                    break;
                case "A, dec.":
                    color = COLOR_SCHEMES.nucleotides.a;
                    break;
                case "C":
                    color = COLOR_SCHEMES.nucleotides.C;
                    break;
                case "C, dec.":
                    color = COLOR_SCHEMES.nucleotides.c;
                    break;
                case "G":
                    color = COLOR_SCHEMES.nucleotides.G;
                    break;
                case "G, dec.":
                    color = COLOR_SCHEMES.nucleotides.g;
                    break;
                case "T":
                    color = COLOR_SCHEMES.nucleotides.T;
                    break;
                case "T, dec.":
                    color = COLOR_SCHEMES.nucleotides.t;
                    break;
                case "Ref.":
                    color = COLOR_SCHEMES.other[ "." ];
                    break;
                case "Ref., dec.":
                    color = COLOR_SCHEMES.other[ ":" ];
                    break;
                case "No Call":
                    color = COLOR_SCHEMES.other[ "N" ];
                    break;
                case "Del.":
                    color = COLOR_SCHEMES.other[ "-" ];
                    break;
                }
                outerPieData.push( {
                    name: key,
                    value: value,
                    itemStyle: {
                        color: color
                    }
                } );
            }
        }
    } else {
        outerPieData.push( {
            name: "Ref.",
            value: SAMPLE_LABELS.length,
            itemStyle: {
                color: COLOR_SCHEMES.other[ "." ]
            }
        } );
    }
    /* Collect count summary wrt. induced amino-acid substitutions of position */
    var innerPieData = [ ];
    let aaCounts = {
        "HIS": 0,
        "LYS": 0,
        "ARG": 0,
        "ASP": 0,
        "GLU": 0,
        "SER": 0,
        "THR": 0,
        "ASN": 0,
        "GLN": 0,
        "CYS": 0,
        "PHE": 0,
        "TRP": 0,
        "TYR": 0,
        "ALA": 0,
        "VAL": 0,
        "LEU": 0,
        "ILE": 0,
        "MET": 0,
        "PRO": 0,
        "GLY": 0,
        "STOP": 0
    };
    let totalAlt = 0;
    for ( const [_, value] of Object.entries( DATA[ SELECTED_CHAIN ][ ALIGNMENT_POSITION_MAP[ position ][ 0 ] ][ "missenseMutations" ] ) ) {
        let aaAlt = value.split( ">" )[ 1 ];
        totalAlt += 1;
        aaCounts[ aaAlt ] += 1;
    }
    let nones = SAMPLE_LABELS.length - totalAlt; 
    if ( nones > 0 ) {
        innerPieData.push( {
            name: "None",
            value: nones,
            itemStyle: {
                color: COLOR_SCHEMES.other[ "." ]
            }
        } );
    }
    for ( const [key, value] of Object.entries( aaCounts ) ) {
        if ( value > 0 ) {
            if ( key == "STOP" ) {
                innerPieData.push( {
                    name: "Termination",
                    value: value,
                    itemStyle: {
                        color: "#8c00ff"
                    }
                } );
            } else {
                innerPieData.push( {
                    name: key,
                    value: value,
                    itemStyle: {
                        color: COLOR_SCHEMES.aminoAcids3[ key ]
                    }
                } );
            }
        }
    }
    var titleSubText = "Alignment: " + position + "; Genome: " + ALIGNMENT_POSITION_MAP[ position ][ 0 ] + "; Protein: " + ( ALIGNMENT_POSITION_MAP[ position ][ 1 ] !== "_" ? ALIGNMENT_POSITION_MAP[ position ][ 1 ] : "None" );
    document.getElementById( "positionSummaryViewHeader" ).innerHTML = "Position Summary | " + titleSubText;
    ECHART_PS_OPTION = {
        series: [
            {
                name: 'Induced Substitutions',
                type: 'pie',
                radius: [ '35%', '25%' ],
                emphasis: {
                    scale: false
                },
                labelLine: {
                    length: 20,
                    length2: 10
                },
                label: {
                    formatter: '{b|{b}}\n{num|{d}% ({c})}',
                    rich: {
                    b: {
                        color: '#292926',
                        fontSize: 11,
                        fontWeight: 'bold',
                        align: 'center'
                    },
                    num: {
                        color: '#292926',
                        fontSize: 11,
                        backgroundColor: '#E4E5ED',
                        padding: [3, 4],
                        borderRadius: 4,
                        align: 'left'
                    }
                    }
                },
                data: innerPieData
            },
            {
                name: 'Nucleotide Composition',
                type: 'pie',
                radius: ['70%', '80%'],
                emphasis: {
                    scale: false
                },
                labelLine: {
                    length: 15,
                    length2: 5
                },
                label: {
                    formatter: '{b|{b}}\n{num|{d}% ({c})}',
                    rich: {
                    b: {
                        color: '#292926',
                        fontSize: 11,
                        fontWeight: 'bold',
                        align: 'center'
                    },
                    num: {
                        color: '#292926',
                        fontSize: 11,
                        backgroundColor: '#E4E5ED',
                        padding: [3, 4],
                        borderRadius: 4,
                        align: 'left'
                    }
                    }
                },
                data: outerPieData
            }
          ]
    };
    ECHART_PS.setOption( ECHART_PS_OPTION );
}

function PVSetDefaultStyle( ) {
    PROTEIN_VIEWER.setStyle(
        { },
        {
            cartoon: MOL_STYLES.default_cartoon,
            line: MOL_STYLES.default_line
        }
    );
    // Re-apply selection style.
    for (  let selectedResidue of PROTEIN_VIEWER_SELECTED_RESIDUES ) {
        PROTEIN_VIEWER.addStyle(
            {
                chain: selectedResidue[ 1 ],
                resi: selectedResidue[ 0 ]
            },
            {
                stick: {
                    radius: 0.2,
                    color: COLOR_SCHEMES.aminoAcids3[ selectedResidue[ 2 ] ]
                }
            }
        );
    }
    PROTEIN_VIEWER.render( );
};

function PVHighlightRegion( from, to, chain ) {
    PROTEIN_VIEWER.setStyle(
        {
            chain: chain,
            resi: from.toString( ) + "-" + to.toString( )
        },
        {
            cartoon: MOL_STYLES.highlight_cartoon,
            line: MOL_STYLES.default_line
        }
    );
    PROTEIN_VIEWER_HIGHLIGHTED_REGION = new Set( );
    for ( let resi = parseInt( from ); resi <= parseInt( to ); resi++ ) {
        PROTEIN_VIEWER_HIGHLIGHTED_REGION.add( chain + resi.toString( ) );
    }
    // Re-apply selection style.
    for (  let selectedResidue of PROTEIN_VIEWER_SELECTED_RESIDUES ) {
        PROTEIN_VIEWER.addStyle(
            {
                chain: selectedResidue[ 1 ],
                resi: selectedResidue[ 0 ]
            },
            {
                stick: {
                    radius: 0.2,
                    color: COLOR_SCHEMES.aminoAcids3[ selectedResidue[ 2 ] ]
                }
            }
        );
    }
    PROTEIN_VIEWER.render( );
}

function PVClearSelection( ) {
    // Remove existing selection.
    if ( PROTEIN_VIEWER_SELECTED_RESIDUES.length != 0 ) {
        for ( let selectedResidue of PROTEIN_VIEWER_SELECTED_RESIDUES ) {
            if ( PROTEIN_VIEWER_HIGHLIGHTED_REGION.has( selectedResidue[ 1 ] + selectedResidue[ 0 ] ) ) {
                PROTEIN_VIEWER.setStyle(
                    {
                        chain: selectedResidue[ 1 ],
                        resi: selectedResidue[ 0 ]
                    },
                    {
                        cartoon: MOL_STYLES.highlight_cartoon,
                        line: MOL_STYLES.default_line
                    }
                );
            } else {
                PROTEIN_VIEWER.setStyle(
                    {
                        chain: selectedResidue[ 1 ],
                        resi: selectedResidue[ 0 ]
                    },
                    {
                        cartoon: MOL_STYLES.default_cartoon,
                        line: MOL_STYLES.default_line
                    }
                );
            }
            PROTEIN_VIEWER.removeLabel( selectedResidue[ 3 ] );
        }
        PROTEIN_VIEWER_SELECTED_RESIDUES = [ ];
    }
}

function PVSelectResidue( resi, chain, resn, clear ) {
    // Remove existing selection.
    if ( clear ) {
        PVClearSelection( );
    }
    // Visualize new selection.
    if ( resi == "_" ) {
        return;
    }
    PROTEIN_VIEWER.addStyle(
        {
            chain: chain,
            resi: resi
        },
        {
            stick: {
                radius: 0.2,
                color: COLOR_SCHEMES.aminoAcids3[ resn ]
            }
        }
    );
    var label = PROTEIN_VIEWER.addLabel(
        chain + resi + ": " + resn,
        {
            backgroundColor: "#E4E5ED",
            backgroundOpacity: "0.8",
            fontColor: "#292926",
            fontSize: 11,
            borderColor: "#292926",
            borderThickness: 0.2
        },
        {
            chain: chain,
            resi: resi
        },
        false
    );
    PROTEIN_VIEWER_SELECTED_RESIDUES.push( [ resi, chain, resn, label ] );
    PROTEIN_VIEWER.render( );
}

function PVShowContacts( ) {
    if ( PROTEIN_VIEWER_SELECTED_RESIDUES.length != 0 ) {
        let internalState = PROTEIN_VIEWER.getInternalState( );
        let selectedChainIdentifier = PROTEIN_VIEWER_SELECTED_RESIDUES[ 0 ][ 1 ]
        let selectedResidueNumber = PROTEIN_VIEWER_SELECTED_RESIDUES[ 0 ][ 0 ]
        let selectedResidueType = PROTEIN_VIEWER_SELECTED_RESIDUES[ 0 ][ 2 ]
        let selectedSecondaryStructure = PROTEIN_SS_DATA[ selectedResidueNumber + selectedChainIdentifier ];
        let key = selectedResidueNumber + selectedChainIdentifier;
        if ( key in PROTEIN_RESIDUE_CONTACTS ) {
            let contactInformationHtmlString = "<hr><b>Contact Information</b><br>";
            let residuesInContact = PROTEIN_RESIDUE_CONTACTS[ key ];
            if ( residuesInContact.length > 0 ) {
                for ( let residueInContact of residuesInContact ) {
                    let contactResidueIdentifier = residueInContact[ 0 ];
                    let contactChainIdentifier = contactResidueIdentifier.slice( -1 );
                    let contactResidueNumber = contactResidueIdentifier.slice( 0, -1 );
                    let contactResidueType = "?"
                    let contactDistance = residueInContact[ 1 ];
                    /* Gather residue information from 3DMol.js GLModel */
                    for ( const atom of internalState.models[ 0 ].atoms ) {
                        if ( atom.resi == contactResidueNumber && atom.chain == contactChainIdentifier && atom.atom == "CA" ) {
                            contactResidueType = atom.resn;
                            break;
                        }
                    }
                    let contactSecondaryStructure = PROTEIN_SS_DATA[ contactResidueNumber + contactChainIdentifier ];
                    PVSelectResidue( contactResidueNumber, contactChainIdentifier, contactResidueType, false );
                    contactInformationHtmlString += "<i>" + contactChainIdentifier + contactResidueNumber + ": " + contactResidueType + ", " + contactSecondaryStructure + "</i><br>";
                    contactInformationHtmlString += "Distance (SCM): " + parseFloat( contactDistance ).toFixed( 2 ) + " Å<br>";
                    /* Assing PRI15 scores */
                    let pri15ScoreWildtype = PRI15Scores[ selectedResidueType + "-" + selectedSecondaryStructure + "/" + contactResidueType + "-" + contactSecondaryStructure ];
                    if ( pri15ScoreWildtype == "-Infinity" ) {
                        pri15ScoreWildtype = "Not observed."
                    } else {
                        pri15ScoreWildtype = parseFloat( pri15ScoreWildtype ).toFixed( 2 );
                    }
                    contactInformationHtmlString += "Interaction Score (Wildtype): " + pri15ScoreWildtype + "<br>";
                    if ( INFOInducedSubstitution ) {
                        let pri15ScoreInducedSubstitution = PRI15Scores[ INFOInducedSubstitution + "-" + selectedSecondaryStructure + "/" + contactResidueType + "-" + contactSecondaryStructure ];
                        if ( pri15ScoreInducedSubstitution == "-Infinity" ) {
                            pri15ScoreInducedSubstitution = "Not observed."
                        } else {
                            pri15ScoreInducedSubstitution = parseFloat( pri15ScoreInducedSubstitution ).toFixed( 2 );
                        }
                        contactInformationHtmlString += "Interaction Score (Mut.): " + pri15ScoreInducedSubstitution + "<br>";
                    }
                }
                var label = PROTEIN_VIEWER.addLabel(
                    selectedChainIdentifier + selectedResidueNumber + ": " + selectedResidueType,
                    {
                        backgroundColor: "#fceaca",
                        backgroundOpacity: "0.8",
                        fontColor: "#292926",
                        fontSize: 11,
                        borderColor: "#ff8400",
                        borderThickness: 1
                    },
                    {
                        chain: selectedChainIdentifier,
                        resi: selectedResidueNumber
                    },
                    false
                );
                PROTEIN_VIEWER.removeLabel( PROTEIN_VIEWER_SELECTED_RESIDUES[ 0 ][ 3 ] );
                PROTEIN_VIEWER_SELECTED_RESIDUES[ 0 ][ 3 ] = label;
                PROTEIN_VIEWER.render();
            } else {
                contactInformationHtmlString += "No contacts"
            }
            /* Add contact information to Info component */
            INFOAddContactInfo( contactInformationHtmlString );
        }
    }
}