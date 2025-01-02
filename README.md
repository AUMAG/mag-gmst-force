# Force and torque calculations between cylindrical arc magnets and coils
This page features a [Mathematica notebook](https://github.com/AUMAG/mag-gmst-force/blob/8b237dd854b60771678d33ff9b432291a644afba/Force%20and%20Torque%20Between%20Cylindrical%20Arc%20Coils%20and%20Magnets.nb) ([PDF version](https://github.com/AUMAG/mag-gmst-force/blob/8b237dd854b60771678d33ff9b432291a644afba/Force%20and%20Torque%20Between%20Cylindrical%20Arc%20Coils%20and%20Magnets.pdf)) that is supplementary material for the research article ‘A Generalised Maxwell Stress Tensor for Semi-Analytic Force and Torque Between Permanent Magnets, Coils, and Soft Iron’. The notebook contains an example case study implementing the Generalised and Classical Maxwell Stress Tensors (GMST, MST) from the article, and compares the results to a Finite Element Analysis result (from ANSYS Maxwell 3D). Magnetic field results can be generated using prior work, with a package [MagCylField](https://github.com/AUMAG/mag-gmst-force/blob/8b237dd854b60771678d33ff9b432291a644afba/MagCylField.wl) containing all analytic solutions and an [example Mathematica document](https://github.com/AUMAG/mag-gmst-force/blob/8b237dd854b60771678d33ff9b432291a644afba/MagCylField-package-test/MagCylField%20example.nb) for creating field plots.

A summary for implementing this elemental modelling method for linear/rotary motors is given below, with a detailed outlined in future work (preprint to come) and past work ([10.1109/TMAG.2020.3010566](https://doi.org/10.1109/TMAG.2020.3010566)).

1. For a known analytic solution, assign the coil and permanent magnet: geometry, location, and magnetisation/current. If there is no permeable matter, then skip to step 6.
1. Calculate the total magnetic flux density from all coils and permanent magnets with unity relative permeability across the regions with permeable matter (initially assumed to be a vacuum).
1. Using the prior magnetic flux density and the B-H curve for the permeable material, assign an initial uniform vector magnetisation to each volume element.
1. Calculate the total magnetic flux density from all new volume elements at all centroids defined in step 2.
1. Follow the same process of step 3, assigning new magnetisations. If convergence criteria are met, then proceed to step 6.  Otherwise, repeat step 4.
1. Form a closed-surface in the airgap around either the stator or rotor/mover volume. Mesh the surface and calculate the discretised magnetic field at each area centroid from all sources. Calculate the force and torque using GMST, with the local surface normal, area, radii, and field.

## Citation

If you use this work, please cite the appropriate paper using the following references:

    @Article {Forbes2024,
        author = {Forbes, M. and Robertson, W.S.P. and Zander, A.C. and Paulides, J.J.H},
        journal = {Advanced Physics Research},
        title = {The Magnetic Field from Cylindrical Arc Coils and Magnets: {A} Compendium with New Analytic Solutions for Radial Magnetisation and Azimuthal Current},
        doi = {10.1002/apxr.202300136},
        publisher = {Wiley},
        number = {7},
        pages = {2300136},
        volume = {3},
    }

    @Article {Forbes2025,
        author = {Forbes, M. and Robertson, W.S.P. and Zander, A.C. and Vidler, J. and Paulides, J.J.H},
        journal = {Preprint},
        title = {A Generalised Maxwell Stress Tensor for Semi-Analytic Force and Torque Between Permanent Magnets, Coils, and Soft Iron},
    }
