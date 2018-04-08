# Lost-circulation
streamline modeling
pressure distribution and geomechanics modeling.
roughness induced pressure drop calculation


1. wellbore model
2. formation model   (Euler forward method to numerical solve mud invasion radius)
3. well trajectory
4. geomechanics
5. coupling model
6. cutting transport and mud density
7. coding

The drilling can induce the pore pressure change of fault. and the fracture pressure should be change also well.

so we can use the green function and boundary integration to calculate the pressure change of fault.

See the pressure change code of Dr. Younis!

Non-isotherm hydraulic model, when comparision with the field data, the error only for 0.63% !!!!!!!!!!



!!! caution, the numerical solution is not stable, the perturbation equation is unstable.

check the second high-oder, and the Taylor series!
