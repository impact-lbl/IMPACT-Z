Author: Biaobin Li

Email: biaobin@ustc.edu.cn

# Introduction

This python script could transform string-type input file into `ImpactZ.in` file. Users have experiences with `ELEGANT` will enjoy this tool. Right now, not all `IMPACT-Z` elements are added in the code. Users could go to the python scipts: `impactz_parser.py/__default_lattice()` and `impactz_parser.py/impactzin_lattice()` to add the new elements you want to use. 



# How to run it

Add `genimpactzin` into your `PATH`, for my  case, add the following line to your `.bashrc`:

```bash
export PATH=/mnt/d/githubProj/IMPACT-Z/utilities/lattice_parser:$PATH
```

Given the `lte.impz` input file:

```bash
genimpactzin lte.impz line
```

which will generate the `ImpactZ.in` file. Now you can run the `ImpactZexe` in parallel version as:

```
mpirun -np 4 ImpactZexe
```

4 processes are used as $core\_num\_T\times core\_num\_L=4$.

If you are using the single process version, `core_num_T=1,core_num_L=1` should be given. Then just type:

```
ImpactZ.exe
```

`ImpactZ.in` is automatically read.

The user is encouraged to have a look in `utilities/lattice_parser/examples`, one example is given to show how this work.



# A simple example

For the convenience of illustration, `lte.impz` refers to the python level read-in file, `ImpactZ.in` refers to `ImpactZexe` read-in file. You can rename `lte.impz`  any other names you like.

The `lte.impz` file mainly consists of three sections, `control, beam and lattice` sections. The detailed mapping relationships between`ImpactZ.in` and `lte.impz` are listed in the following section. Here we give a simple example for the  usage of `lte.impz` :

```python
!control section
!===============
&control

core_num_T = 2;
core_num_L = 2;
meshx = 32;
meshy = 32;
meshz = 64;
kinetic_energy = 300e6;
freq_scale = 1.3e9;

&end

!beam section
!==============
&beam
mass = 0.511001e6;
charge = -1.0;

distribution_type = 2;
Np = 5000;
total_charge = 1e-9;

emit_nx=0.176e-6, beta_x=12.73, alpha_x=-0.85;
emit_ny=0.176e-6, beta_y=12.73, alpha_y=-0.85;

sigz=1e-3, sigE=5e3;

&end

!lattice section
!=====================
&lattice

!rpn expression is supported,
!only a few mathematical operator are added, please see 
!lattice_parser.py/rpn_cal() for more details.
!------------------------------------------------------
% 0.2 sto LB1
% -4.410 pi * 180 / sto AB1  ! Bend angle

BCX11: BEND,L= LB1,ANGLE=AB1,       E2=AB1,       steps=1, pipe_radius=2.1640E-02, fint=0.3893
BCX12: BEND,L= LB1,ANGLE= "0 AB1 -",E1= "0 AB1 -",steps=1, pipe_radius=2.1640E-02, fint=0.3893
BCX13: BEND,L= LB1,ANGLE= "0 AB1 -",E2= "0 AB1 -",steps=1, pipe_radius=2.1640E-02, fint=0.3893
BCX14: BEND,L= LB1,ANGLE=AB1,       E1=AB1,       steps=1, pipe_radius=2.1640E-02, fint=0.3893

D1  : DRIF, L=5.0
Dm  : DRIF, L=0.5

BC1 : LINE=(BCX11,D1,BCX12,Dm,BCX13,D1,BCX14)

W0:   watch, filename_ID=1000
W1:   watch, filename_ID=1001

Line : LINE=(W0,BC1,W1)

&end
```



# Control and beam section

The mapping relationship between `ImpactZ.in` and `lte.impz` in `control` sections are listed as following:

```bash
line1: 
core_num_T core_num_L

line2:
6 Np integrator error output_ratio

line3：
meshx meshy meshz flagbc x_pipe_width y_pipe_width period_len

line4:
distribution_type restart sub_cycle 1

line5:
Np

line6:
current  #Q*f_scale

line7:
# value defined automatically by charge/mass in the code; 
# the definition of charge and mass, see beam section.

line8-line10:
# defined by beam section
alpha_x beta_x emit_x mismatchx mismatchpx offsetX     offsetPx
alpha_y beta_y emit_y mismatchy mismatchpy offsetY     offsetPy
alpha_z beta_z emit_z mismatchz mismatchE  offsetPhase offsetEnergy

line11:
current kinetic_energy mass charge freq_scale ini_phase 

```



## Control parameters

All control parameters in `lte.impz` are listed:

| Parameter Name | Units | Type   | Default | Description                                                  |
| -------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| core_num_T     |       | int    | 1       | processor number for the transverse direction.               |
| core_num_L     |       | int    | 1       | processor number for the longitudinal direction.             |
| integrator     |       | int    | 1       | 1 for linear map integrator, 2 for nonlinear Lorentz integrator. |
| error          |       | int    | 0       | 0 for OFF, 1 for ON error studies.                           |
| output_ratio   |       | int    | 1       | 1 for standard output.                                       |
| meshx          |       | int    | 64      | space charge grid for x-direction.                           |
| meshy          |       | int    | 64      | space charge grid for y-direction.                           |
| meshz          |       | int    | 64      | space charge grid for z-direction.                           |
| flagbc         |       | int    | 1       | space charge boundary situation, 1 for 3D open.  See Ji’s manual for more details. |
| x_pipe_width   | m     | double | 0.014   | x pipe width.                                                |
| y_pipe_width   | m     | double | 0.014   | y pipe width.                                                |
| period_len     | m     | double | 0.10    | period length.                                               |
| restart        |       | int    | 0       | 0 for OFF restart, 1 for ON. Restart from some point after stop. |
| sub_cycle      |       | int    | 0       | 0 for no sub-cycle, 1 for ON.                                |
| kinetic_energy | eV    | double | 0       | kinetic energy of the beam.                                  |
| freq_rf_scale  | Hz    | double | 2856e6  | scale frequency $f_{scal}$,  $Scxl=c/(2\pi f_{scal}) $.      |
| ini_phase      |       | double | 0.0     | initial phase of the reference particle.                     |
| steps_permeter |       | int    | 0       | how many sc kicks per meter in a single element. `steps_permeter=0`, element `steps=1`. |
| maps_permeter  |       | int    | 0       | how many maps per meter in a single element. `maps_permeter=0`, element `maps=1`. |
| sample_out     |       | int    | 1e5     | how many particles to sample out in `watch` elements.        |
| slice_bin      |       | int    | 128     | how many slice bins are used in `watch` elements.            |
| pipe_radius    | m     | double | 0.014   | pipe radius for all elements, which are not given values to `pipe_radius` by users in lattice section (default is `pipe_radius=0.0` m). |



## Beam parameters

All beam section parameters in `lte.impz` are listed:

| Parameter Name    | Units   | Type   | Default | Description                                                |
| ----------------- | ------- | ------ | ------- | ---------------------------------------------------------- |
| mass              | eV      | double | 0.511e6 | mass of the particle.                                      |
| charge            |         | double | -1.0    | -1 for electron.                                           |
| distribution_type |         | int    | 2       | 6D gaussian distribution. See more options in Ji’s manual. |
| Np                |         | int    | 1e3     | particle number.                                           |
| total_charge      | C       | double | 1e-9    | charge of the beam.                                        |
| emit_x            | m rad   | double | 0.0     | emitance.                                                  |
| emit_nx           | m rad   | double | 0.0     | normalized emittance.                                      |
| beta_x            | m       | double | 1.0     | twiss para.                                                |
| alpha_x           |         | double | 0.0     | twiss para.                                                |
| sigx              | m       | double | 0.0     | rms bunch size.                                            |
| sigpx             |         | double | 0.0     | rms value of $\gamma\beta_x/\gamma_0\beta_0$               |
| emit_y            | m rad   | double | 0.0     | emittance.                                                 |
| emit_ny           | m rad   | double | 0.0     | normalized emittance.                                      |
| beta_y            | m       | double | 1.0     | twiss para.                                                |
| alpha_y           |         | double | 0.0     | twiss para.                                                |
| sigy              | m       | double | 0.0     | rms bunch size.                                            |
| sigpy             |         | double | 0.0     | rms value of $\gamma\beta_y/\gamma_0\beta_0$               |
| emit_z            | deg MeV | double | 0.0     | twiss para.                                                |
| beta_z            | deg/MeV | double | 1.0     | twiss para.                                                |
| alpha_z           |         | double | 0.0     | twiss para.                                                |
| sigz              | m       | double | 0.0     | rms bunch length.                                          |
| sigE              | eV      | double | 0.0     | rms energy spread.                                         |

Users could either use twiss parameters to define initial beam distribution, or use rms values. For $\sigma_{ij}\neq0$ cases, please use twiss-para. 



# Lattice section

Right now, only a few frequently used elements in `ImpactZ.in` are added into the python parser.



## Elements

### DRIFT

0 element.

| Parameter Name | Units | Type   | Default | Description                                                  |
| -------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| L              | m     | double | 0.0     | length of drift                                              |
| steps          |       | int    | 0       | how many sc kicks for this element.                          |
| maps           |       | int    | 0       | each half-drift involves computing a map for that half-element, computed by numerical integration with 1 maps |
| pipe_radius    | m     | double | 0.0     | pipe radius                                                  |



### QUAD

1 element.

| Parameter Name | Units       | Type   | Default | Description                                                  |
| -------------- | ----------- | ------ | ------- | ------------------------------------------------------------ |
| L              | m           | double | 0.0     | length                                                       |
| steps          |             | int    | 0       | how many sc kicks for this element.                          |
| maps           |             | int    | 0       | map steps.                                                   |
| $K_1$          | $\rm{/m^2}$ | double | 0.0     | quadrupole strength, $K_1=\frac{1}{(B\rho)_0}\frac{\partial B_y}{\partial x}$ |
|                |             |        |         |                                                              |
| pipe_radius    | m           | double | 0.0     | pipe radius                                                  |
| Dx             | m           | double | 0.0     | x misalignment error                                         |
| Dy             | m           | double | 0.0     | y misalignment error                                         |
| rotate_x       | rad         | double | 0.0     | rotation error in x direction                                |
| rotate_y       | rad         | double | 0.0     | rotation error in y direction                                |
| ratate_z       | rad         | double | 0.0     | rotation error in y direction                                |



### BEND

4 element. A magnetic dipole implemented as a matrix, up to 2nd order. See K. Brown paper for more information.

| Parameter Name | Units        | Type   | Default | Description                                                  |
| -------------- | ------------ | ------ | ------- | ------------------------------------------------------------ |
| L              | m            | double | 0.0     | arc length                                                   |
| steps          |              | int    | 0       | how many SC/CSR kicks for this element.                      |
| maps           |              | int    | 0       | map steps                                                    |
| angle          | rad          | double | 0.0     | bend angle                                                   |
| E1             | rad          | double | 0.0     | entrance edge angle                                          |
| E2             | rad          | double | 0.0     | exit edge angle                                              |
| $K_1$          | $1\rm{/m^2}$ | double | 0.0     | quadrupole strength, $K_1=\frac{1}{(B\rho)_0}\frac{\partial B_y}{\partial x}$, not added yet in V2.1 version. |
| PIPE_RADIUS    | m            | double | 0.0     | pipe radius                                                  |
| h1             | 1/m          | double | 0.0     | entrance pole-face curvature                                 |
| h2             | 1/m          | double | 0.0     | exit pole-face curvature                                     |
| fint           |              | double | 0.5     | integrated fringe field, set to 0.5 to keep same with ELEGANT. |
| Dx             | m            | double | 0.0     | x misalignment error                                         |
| Dy             | m            | double | 0.0     | y misalignment error                                         |
| rotate_x       | rad          | double | 0.0     | rotation error in x direction                                |
| rotate_y       | rad          | double | 0.0     | rotation error in y direction                                |
| ratate_z       | rad          | double | 0.0     | rotation error in y direction. This parameter could be used as TILT from ELEGANT. BUT REMEMBER TO SET ERROR=1 IN CONTROL SECTION. |
| CSR            |              | int    | 0       | 0/1, whether to include 1D-CSR effects or not.               |



### RFCW

Ideal sinusoidal RF model, combined with -41 wake element.

| Parameter Name | Units  | Type   | Default | Description                                                  |
| -------------- | ------ | ------ | ------- | ------------------------------------------------------------ |
| L              | m      | double | 0.0     | length                                                       |
| steps          |        | int    | 0       | how many SC kicks for this element.                          |
| maps           |        | int    | 0       | map steps                                                    |
| volt           | V      | double | 0.0     | peak voltage                                                 |
| gradient       | V/m    | double | 0.0     | peak acceleration gradient.                                  |
| phase          | degree | double | 0.0     | driven phase,  sin() function is used (same as ELEGANT, different with IMPACT-Z), $E_z=A\cdot \rm{sin}(kz+\phi)$, phase=90 is the crest for acceleration |
| freq           | Hz     | double | 2.856e9 | RF frequency                                                 |
| pipe_radius    | m      | double | 0.0     | pipe radius                                                  |
| Dx             | m      | double | 0.0     | x misalignment error                                         |
| Dy             | m      | double | 0.0     | y misalignment error                                         |
| rotate_x       | rad    | double | 0.0     | rotation error in x direction                                |
| rotate_y       | rad    | double | 0.0     | rotation error in y direction                                |
| ratate_z       | rad    | double | 0.0     | rotation error in z direction                                |
|                |        |        |         |                                                              |
| wakeflag       |        | int    | -1      | `wakeflag=-1` turn off RF wakefield, `wakeflag=5`  only turn on longitudinal wake, `wakeflag=15` include both longitudinal and transverse wake. |
| wakefile_ID    |        | int    | None    | If WAKEFIEL_ID=41, it refers to `rfdata41.in` , which contains RF structure wakefield, 1st column is s [m],  2nd column is longitudinal wakefield $w_L$ [V/C/m], 3rd and 4th columns are transverse wakefield $w_x, w_y$ [$\rm{V/C/m^2}$]. |



### DTL

`101` element.

| Parameter Name | Units  | Type   | Default | Description                                                  |
| -------------- | ------ | ------ | ------- | ------------------------------------------------------------ |
| L              | m      | double | 0.0     | length                                                       |
| steps          |        | int    | 0       | how many segments  for the element. DIFFERENT from steps in control section, not nseg/m. |
| maps           |        | int    | 0       | map steps                                                    |
| scale          |        | double | 1.0     |                                                              |
| freq           | Hz     | double | 324e6   | RF frequency                                                 |
| phase          | degree | double | 0.0     | driven phase,  sin() function is used (same as ELEGANT, different with IMPACT-Z), $E_z=A\cdot \rm{sin}(kz+\phi)$, phase=90 is the crest for acceleration |
| ID             |        | int    | 100     | file ID for the external field                               |
| pipe_radius    | m      | double | 0.0     | pipe radius                                                  |
| Lq1            | m      | double | 0.0     | quad 1 length                                                |
| grad1          | T/m    | double | 0.0     | quad 1 gradient                                              |
| Lq2            | m      | double | 0.0     | quad 2 length                                                |
| grad2          | T/m    | double | 0.0     | quad 2 gradient                                              |
|                |        |        |         |                                                              |
| Dx_q           | m      | double | 0.0     | x misalignment error for quad                                |
| Dy_q           | m      | double | 0.0     | y misalignment error                                         |
| rotate_x_q     | rad    | double | 0.0     | rotation error in x direction                                |
| rotate_y_q     | rad    | double | 0.0     | rotation error in y direction                                |
| ratate_z_q     | rad    | double | 0.0     | rotation error in z direction                                |
|                |        |        |         |                                                              |
| Dx_RF          | m      | double | 0.0     | x misalignment error for quad                                |
| Dy_RF          | m      | double | 0.0     | y misalignment error                                         |
| rotate_x_RF    | rad    | double | 0.0     | rotation error in x direction                                |
| rotate_y_RF    | rad    | double | 0.0     | rotation error in y direction                                |
| ratate_z_RF    | rad    | double | 0.0     | rotation error in z direction                                |



### SC

104 element.

| Parameter Name | Units  | Type   | Default | Description                                                  |
| -------------- | ------ | ------ | ------- | ------------------------------------------------------------ |
| L              | m      | double | 0.0     | length                                                       |
| steps          |        | int    | 0       | how many segments  for the element. DIFFERENT from steps in control section, not nseg/m. |
| maps           |        | int    | 0       | map steps                                                    |
| scale          |        | double | 1.0     | field scale factor                                           |
| ID             |        | int    | 100     | file ID for the external field                               |
| phase          | degree | double | 0.0     | driven phase,  sin() function is used (same as ELEGANT, different with IMPACT-Z), $E_z=A\cdot \rm{sin}(kz+\phi)$, phase=90 is the crest for acceleration |
| freq           | Hz     | double | 324e6   | RF frequency                                                 |
| pipe_radius    | m      | double | 0.0     | pipe radius                                                  |
| Dx             | m      | double | 0.0     | x misalignment error                                         |
| Dy             | m      | double | 0.0     | y misalignment error                                         |
| rotate_x       | rad    | double | 0.0     | rotation error in x direction                                |
| rotate_y       | rad    | double | 0.0     | rotation error in y direction                                |
| ratate_z       | rad    | double | 0.0     | rotation error in z direction                                |



### FIELDMAP

110 element.

| Parameter Name | Units  | Type   | Default | Description                                                  |
| -------------- | ------ | ------ | ------- | ------------------------------------------------------------ |
| L              | m      | double | 0.0     | length                                                       |
| steps          |        | int    | 0       | how many segments  for the element. DIFFERENT from steps in control section, not nseg/m. |
| maps           |        | int    | 0       | map steps                                                    |
| scale          |        | double | 1.0     |                                                              |
| ID             |        | int    | 100     | file ID for the external field                               |
| phase          | degree | double | 0.0     | driven phase,  sin() function is used (same as ELEGANT, different with IMPACT-Z), $E_z=A\cdot \rm{sin}(kz+\phi)$, phase=90 is the crest for acceleration |
| freq           | Hz     | double | 324e6   | RF frequency                                                 |
| Xradius        | m      | double | 0.0     | pipe radius                                                  |
| Yradius        | m      | double | 0.0     | pipe radius                                                  |
| Dx             | m      | double | 0.0     | x misalignment error                                         |
| Dy             | m      | double | 0.0     | y misalignment error                                         |
| rotate_x       | rad    | double | 0.0     | rotation error in x direction                                |
| rotate_y       | rad    | double | 0.0     | rotation error in y direction                                |
| ratate_z       | rad    | double | 0.0     | rotation error in z direction                                |
| datatype       |        | int    | 1       | 1 using discrete data, 2 using both discrete data and analytical function, other using analytical function only |
| coordinate     |        | int    | 2       | 1 in cylindrical coordinate, 2 in Cartesian                  |



### SHIFTCENTER

-1 and -19 element.

| Parameter Name | Units | Type   | Default | Description                                                  |
| -------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| option         |       | string | “zdE”   | (1). `option=”zdE”`, then shift the beam longitudinally to the bunch centroid so that $<z>=<\Delta E>=0$. (2). `option=“xy”`, then shift the beam so that $<x>=<y>=0$. |

usage:

```
elem1: shiftcenter, option="zdE";
```



### WATCH

-2 and -8 elements.

| Parameter Name | Units | Type | Default | Description                                                  |
| -------------- | ----- | ---- | ------- | ------------------------------------------------------------ |
| filename_ID    |       | int  | 9999    | number larger than 1000 is recommended.                      |
| sample_freq    |       | int  | 0       | particles sample out frequency. If not set, `sample_out` in `control` section will take effects. |
| slice_bin      |       | int  | 0       | bins number for getting histogram slice information. If not set, `slice_bin` in `control` section will take effects. |
| coord_info     |       | int  | 1       | `0/1`, whether to output particles phase space, i.e. whether to add -2 element. |
| slice_info     |       | int  | 1       | `0/1`, whether to output slice information, i.e. whether to add -8 element. |

Output particle distribution and beam slice information into `fort.N` and `fort.(N+1e4)` files repectively, where N is the filename_ID. 

If `filename_ID = 1001`, then the output file would be `fort.1001` and `fort.11001`. `fort.11001` refers to Impact-Z `-8` element output file (+10000), which outputs slice information. The columns in this file are as following:

| Column number | Units | Description                            |
| ------------- | ----- | -------------------------------------- |
| 1             | m     | bunch length                           |
| 2             |       | particle number per slice              |
| 3             | A     | current per slice                      |
| 4             | mrad  | x-direction normalized slice emittance |
| 5             | mrad  | y-direction normalized slice emittance |
| 6             |       | relative slice energy spread, dE/E     |
| 7             | eV    | uncorrelated energy spread per slice   |
| 8             | m     | $<x>$                                  |
| 9             | m     | $<y>$                                  |



### RCOL

-13 element, collimate the beam with transverse rectangular aperture sizes. Name same with Elegant.

| Parameter Name | Units | Type   | Default | Description                                           |
| -------------- | ----- | ------ | ------- | ----------------------------------------------------- |
| x_max          | m     | double | 0.04    | xmax for x direction.                                 |
| y_max          | m     | double | 0.04    | ymax for y direction.                                 |
| x_min          | m     | double | None    | By default, x_min=None, then x_min=-x_max is applied. |
| y_min          | m     | double | None    | By default, y_min=None, then y_min=-y_max is applied. |

The None default values are from the reason that ELEGANT only has `x_max` and `y_max` settings. 



### ECOL

-14 element, collimate the beam with transverse elliptical aperture size.

| Parameter Name | Units | Type   | Default | Description               |
| -------------- | ----- | ------ | ------- | ------------------------- |
| x_max          | m     | double | 0.04    | half-axis in x direction. |
| y_max          | m     | double | 0.04    | half-axis in y direction. |



### ROTATE

-18 element, rotate the beam with respect to the longitudinal axis.

| Parameter Name | Units | Type   | Default | Description                          |
| -------------- | ----- | ------ | ------- | ------------------------------------ |
| angle          | rad   | double | 0       | Both (x,y), and (px,py) are rotated. |

This element refers to `-18` element, the souce code is:

```fortran
do i = 1, innp
    tmpx = Pts1(1,i)*cos(phi)+Pts1(3,i)*sin(phi)
    tmpy = -Pts1(1,i)*sin(phi)+Pts1(3,i)*cos(phi)
    tmppx = Pts1(2,i)*cos(phi)+Pts1(4,i)*sin(phi)
    tmppy = -Pts1(2,i)*sin(phi)+Pts1(4,i)*cos(phi)
    Pts1(1,i) = tmpx
    Pts1(2,i) = tmppx
    Pts1(3,i) = tmpy
    Pts1(4,i) = tmppy
enddo
```



### SCATTER

-20 element.

| Parameter Name | Units | Type   | Default | Description                                                  |
| -------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| dE             | eV    | double | 0       | rms scattering for $\Delta E$ [eV]. dE=1000, then increase energy spread 1keV. |

