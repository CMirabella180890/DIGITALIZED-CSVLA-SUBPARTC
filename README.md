# DIGITALIZED-CSVLA-SUBPARTC
## Introduction 
This Python repository aims to digitalize Certification Specifications for Very Light Airplane, Subpart C through a Model Based approach. ***Model-based-systems engineering*** is a formalized methodology that is used to support the requirements, design, analysis, verification, and validation associated with the development of complex systems. This methodology, when applied to aircraft certification procedures, can lead to a lean, agile and smooth process, faster than the classic approach. Engineers can then focus their effort in exploring the design envelope and further improve the product since early phases of the developing program.

## Main features 
Subpart C is, in general, structured as follow: 
- **General**
- **Flight loads**
- **Control surface and system loads**
- **Horizontal tail surfaces** 
- **Vertical tail surfaces**
- **Supplementary condition for tail surfaces**
- **Ailerons, wing, flaps and special devices**
- **Ground loads**
- **Water loads**
- **Emergency landing conditions**
- **Fatigue evaluation**

Some of these Paragraphs (for instance, Emergency landing conditions and Fatigue evaluation) have to be examined through actual testing. The code inside this repository is focused on Flight envelope diagrams; these diagrams are then used to assess the sizing aerodynamic loads acting on aerostructures (main wing, tailplane, but also command lines, flaps, engine support structures, ...). The final structure of the code is the following: 

### Structure of the code 
1. All the flight envelope diagrams are evaluated, following the airworthiness regulations selected by the user (currently, only CS-VLA airworthiness rules are available).
2. Balancing loads are then calculated. 
3. Shear, Bending, Torsion distributions along the wing semi span are assessed. These are the actual forces and moments that engineers will use to size aerostructures. 
4. Loads on command lines, flaps and engine support structures. 

## Disclaimers
For the sake of clarity, the following disclaimers must be well understood by the interested reader: 
1. The examined, example aircraft is the Tecnam P92. This aircraft is a good example of Very Light Airplane, with a relatively simple wing-body aerodynamics. The author of this repository, unfortunately, had no access to actual wing-body aerodynamics data relative to this aircraft, so **do not rely on the numerical results of the code**! Hopefully, results are inside the "right ballapark", but it must be understood that several, simplifying hypothesis had been made.
2. Graphs and diagrams requires a full LaTeX installation.
3. The repository is not already complete, but all the functionalities of the code will be implemented very soon. 
