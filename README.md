# Sequential-Deposition-MatlabToolbox
This repository contains several Matlab implementations of sequential deposition algorithms, used for the generation of particle packs composed of poly- and mono-dispersed disks. Sequential deposition algorithms are a class of particle-packing methods used to model the formation of materials by adding particles one at a time into a domain, often under simplified physical constraints. These algorithms aim to generate realistic spatial arrangements of particles—typically spheres or disks—by simulating a step-by-step addition process that mimics natural deposition or sedimentation.

A foundational example is Random Sequential Deposition (RSD), where particles are placed randomly in the container, then fall under a quasi-gravitational force, rolling around other stable particles until they settle into a stable position. The process uses closed-form analytical solutions for particle movement, collision, and immobilization, avoiding the need for Newtonian physics (Figure 1: e.g., Visscher and Bolsterli, 1972). This continues until no more particles can be added without violation of the non-overspill condition, resulting in a jammed / saturated structure. RSD is computationally efficient and widely used for studying disordered packings and porous media.

In addition to basic RSD, this toolbox contains several variants, namely (1) _Zoned Sequential Deposition_ (ZSD), for the generation of Patterned Porous Media (see Pavuluri et al. in-review), (2) _Bimodal Random Sequential Deposition_, and (3) _Layered Random Sequential Deposition_.

![collision_fig](https://github.com/user-attachments/assets/9c757651-021c-40c9-a492-08d0980a1573)
**Figure 1:** Collision and stablization conditions for RSD.
