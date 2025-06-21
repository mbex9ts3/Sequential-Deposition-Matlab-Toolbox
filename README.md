# Sequential-Deposition-MatlabToolbox
This repository contains several Matlab implementations of sequential deposition algorithms, used for the generation of particle packs composed of poly- and mono-dispersed disks. Sequential deposition algorithms are a class of particle-packing methods used to model the formation of materials by adding particles one at a time into a domain, often under simplified physical constraints. These algorithms aim to generate realistic spatial arrangements of particles—typically spheres or disks—by simulating a step-by-step addition process that mimics natural deposition or sedimentation.

A foundational example is Random Sequential Deposition (RSD), where particles are placed randomly in the container, then fall under a quasi-gravitational force, rolling around other stable particles until they settle into a stable position. The process uses closed-form analytical solutions for particle movement, collision, and immobilization, avoiding the need for Newtonian physics (Figure 1: e.g., Visscher and Bolsterli, 1972). This continues until no more particles can be added without violation of the non-overspill condition, resulting in a jammed / saturated structure. RSD is computationally efficient and widely used for studying disordered packings and porous media.

In addition to basic RSD, this toolbox contains several variants, namely (1) _Zoned Sequential Deposition_ (ZSD), for the generation of Patterned Porous Media (see Pavuluri et al. in-review / Figure 2), (2) _Bimodal Random Sequential Deposition_ (Figure 3), and (3) _Layered Random Sequential Deposition_. If you use this toolbox in your published work please cite the following article: 

_Controlling capillary fingering morphology in patterned porous1 media. Saideep Pavuluri, Thomas Daniel Seers, Ali Saeibehrouzi, Ran
Holtzman, Soroush Abolfathi, Petr Denissenko, and Harris Sajjad Rabbani._

For questions about its use please contact: thomas.seers@qatar.tamu.edu

![collision_fig](https://github.com/user-attachments/assets/9c757651-021c-40c9-a492-08d0980a1573)
**Figure 1:** Collision and stablization conditions for RSD.

![shark](https://github.com/user-attachments/assets/0de6bb0d-93ba-45ce-8d6e-4068b387e304)
**Figure 2:** Shark patterned porous media.

![layers LQ](https://github.com/user-attachments/assets/6a80a2de-dbc4-4427-b265-c92bd6ca62fa)
**Figure 3:** Layered porous media.

