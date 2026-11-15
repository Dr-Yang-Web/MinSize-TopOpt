# MinSize-TopOpt
Self-Emergent Fillet Feature and Explicit Concurrent Minimum Size Control in Density-Based Topology Optimization
1. Overview
This repository provides the open-source MATLAB implementation accompanying the paper:

Bo Yang, Shuntian Jia, Liang Gao, Hao Li. Self-emergent fillet feature and explicit concurrent minimum size control in density-based topology optimization. Computers & Structures, 2026, Article 108546.
DOI: 10.1016/j.compstruc.2026.108546 (inferred from the manuscript identifier “CAS-108546”; please verify against the final published version).

The code implements an explicit constraint framework based on the interface saturation criterion for density-based topology optimization. It enables the concurrent control of minimum solid thickness, minimum void width, and fillet radii through a single geometric constraint. When only the solid-phase constraint is imposed, fillet features emerge automatically at re-entrant corners while the minimum solid and void sizes are satisfied simultaneously. An optional void-phase constraint provides active control for larger voids or fillets. The method is validated on benchmark problems (e.g., the MBB beam) and produces rounded, manufacturable geometries with minimal stiffness loss.

2. Repository Structure
File / Folder	Description
Run_Top_MinSize.m	Main driver script. Defines the problem, sets parameters, runs the optimization loop, performs FE analysis, applies the size constraints, calls the MMA optimizer, and post-processes the result.
filter_2d.m	Standard density filter (cone filter).
size_filter_2d.m	Circular size filter used for the interface saturation criterion (solid/void phase).
size_sen.m	Computes the size-constraint functions and their sensitivities for the solid and void phases.
E2N.m	Element-to-node density interpolation for smooth geometry representation.
compute_threshold.m	Computes the density cutoff threshold for level-set extraction.
mmasub.m	NOT included in this repository. The MMA optimizer subroutine by Prof. Krister Svanberg must be obtained separately (see Section 4).
subsolv.m	NOT included in this repository. The sub-solver called by mmasub.m, also from Prof. Svanberg’s MMA package.
MBB_top_data_*/	Output folder created by Run_Top_MinSize.m for saving .mat history files.
Note on the base code: The FE analysis and filter structures are adapted from the well-known top88 MATLAB code developed by the Ole Sigmund group at DTU (www.topopt.dtu.dk, 2011).

3. Educational and Non-Commercial Use
This code is distributed for educational and academic research purposes only. It is not intended for commercial use, commercial product design, or integration into proprietary software without explicit permission from the authors. Users are encouraged to study, modify, and extend the code for learning and non-commercial research, provided that the original paper is properly cited.

4. Important: Obtaining the MMA Optimizer
The optimization loop calls the Method of Moving Asymptotes (MMA) solver developed by Prof. Krister Svanberg (KTH Royal Institute of Technology, Stockholm, Sweden). The files mmasub.m and subsolv.m are NOT distributed with this repository because the MMA MATLAB routines are subject to the author’s distribution policy.

To run the code, you must obtain the MMA routines directly from Prof. Svanberg:

Contact / request page: https://people.kth.se/~krille/Welcome.html

Email: krille@math.kth.se

Reference version used in this code: MMA Version September 2007 (with a small change in August 2008). The package contains mmasub.m and subsolv.m.

Please follow Prof. Svanberg’s instructions for academic use. Do not redistribute the MMA files without his consent.

5. Citation
If you use this code or build upon the proposed method in your work, please cite the following paper:

bibtex
@article{Yang2026SelfEmergent,
  author  = {Yang, Bo and Jia, Shuntian and Gao, Liang and Li, Hao},
  title   = {Self-emergent fillet feature and explicit concurrent minimum size control in density-based topology optimization},
  journal = {Computers \& Structures},
  year    = {2026},
  pages   = {108546},
  doi     = {10.1016/j.compstruc.2026.108546},
  note    = {Article number: 108546}
}
Plain-text format:
Bo Yang, Shuntian Jia, Liang Gao, Hao Li. Self-emergent fillet feature and explicit concurrent minimum size control in density-based topology optimization. Computers & Structures, 2026, Article 108546. https://doi.org/10.1016/j.compstruc.2026.108546

(The DOI and article number are based on the manuscript identifier; please update the volume/page fields once the final version is assigned.)

6. Quick Start
Obtain mmasub.m and subsolv.m from Prof. Svanberg (Section 4) and place them in the same directory as Run_Top_MinSize.m.

Ensure that filter_2d.m, size_filter_2d.m, size_sen.m, E2N.m, and compute_threshold.m are also in the MATLAB path.

Open Run_Top_MinSize.m and adjust the parameters at the top of the script:

rs, rv: desired minimum radii for the solid and void phases.

flag: 0 (no size control), 1 (solid-phase constraint only → self-emergent fillets), 2 (bi-phase constraint).

nx, ny, Vsolid: mesh resolution and target volume fraction.

Run the script. The optimization history and final design will be saved in the specified folder.

7. License and Disclaimer
This code is provided as-is, without warranty of any kind, express or implied. The authors are not responsible for any consequences arising from the use of this code. By downloading or using this code, you agree to use it only for lawful, educational, and non-commercial purposes.
