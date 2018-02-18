# xStellTools, v0.1a

A collection of scripts, codes and documents to assist with stellarators.

See also: http://vmecwiki.pppl.wikispaces.net/

Codes here are categorized by directory. 
This is a loose collection of codes, not a one-code-fits-all. 

Code : Description
------------/-----------
Biot_Savart/calc_b/

calc_B_BiotSavart.m:
    Input: Cartesian 'observation point' at which components B (Bx, By, Bz) are desired. Coil: a 'struct' that describes the coil geometry, and the coil current
    Output: Bx, By, Bz
    The method of <ref goes here> is used to rapidly and efficiently calculate the magnetic field components due to currents in the coilset. Check are included to see if the distance between the observxeation point and a coil is less than some specified tolerance.

calc_B_HSX.m:
    Input: Cartesian observation point (P_x, P_y, P_z), the current in the main field coilset, and the taper array describing the current in teh auxilliary coils.
    Output: Bx, By, Bz, Br, Bphi

calc_B_HSX_RPhiZ.m:
    Input: Cylindrical observation point (P_R, P_Phi, P_Z), the current in the main field coilset, and the taper array describing the current in teh auxilliary coils.
    Output: Bx, By, Bz, Br, Bphi

calc_B_LTX.m:
    Input: Cartesian observation point (P_x, P_y, P_z), and an array specfying the current in each of teh coils ([ToroidalField OhmicHeating Red_Upper Red_Lower Orange_U Orange_L Yellow_U Yellow_L Green_U Green_L Blue_U Blue_L Internal_U Internal_L]
    Output: Bx, By, Bz, Br, Bphi

load_HSX_coils.m
    Input: None
    Output: A structure describing the coilset for HSX. This structure is used internally by calc_B_HSX. This file does not need to be executed by the user, except for diagnostic purposes.

load_LTX_coils.m
    Input: None
    Output: A structure describing the coilset for LTX. This structure is used internally by calc_B_LTX. This file does not need to be executed by the user, except for diagnostic purposes.

Bugs happen. Please report them.
~xStellTools team
