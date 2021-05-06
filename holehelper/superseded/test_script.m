clc
clear all
close all

wpath=what;
currentlocation=wpath.path;
if ispc %pc
    addpath([currentlocation]);
    addpath([currentlocation,'\analysis']);
    addpath([currentlocation,'\analysis\cFSM']);
    addpath([currentlocation,'\helpers']);
    addpath([currentlocation,'\interface']);
    addpath([currentlocation,'\plotters']);
    addpath([currentlocation,'\holehelper']);
    addpath(genpath(pwd));
else %mac! or unix
    addpath([currentlocation,'/analysis']);
    addpath([currentlocation,'/analysis/cFSM']);
    addpath([currentlocation,'/helpers']);
    addpath([currentlocation,'/interface']);
    addpath([currentlocation,'/plotters']);
     addpath([currentlocation,'/holehelper']);
    addpath(genpath(pwd));
end

var_length=300;
node=[1 5.0000 1.0000 1 1 1 1 33.333  
2 5.0000 0.0000 1 1 1 1 50.000  
3 2.5000 0.0000 1 1 1 1 50.000  
4 0.0000 0.0000 1 1 1 1 50.000  
5 0.0000 3.0000 1 1 1 1 16.667  
6 0.0000 6.0000 1 1 1 1 -16.667 
7 0.0000 9.0000 1 1 1 1 -50.000 
8 2.5000 9.0000 1 1 1 1 -50.000 
9 5.0000 9.0000 1 1 1 1 -50.000 
10 5.0000 8.0000 1 1 1 1 -33.333
                                ];
elem=[1 1 2 0.100000 100 
2 2 3 0.100000 100 
3 3 4 0.100000 100 
4 4 5 0.100000 100 
5 5 6 0.100000 100 
6 6 7 0.100000 100 
7 7 8 0.100000 100 
8 8 9 0.100000 100 
9 9 10 0.100000 100];

prop=[100 29500.00 29500.00 0.30 0.30 11346.15];

ed_location=[1 7.5 9.5 40   ];
ed_location_x=[80];
P=1;
Mxx=0;
Mzz=0;
M11=0;
M22=0;

distor_len=[50 50];

%node and elem for general information,g,d,l
 [ cross_section_range,node_all,elem_all ] =hole_fun_general( var_length,ed_location,ed_location_x,elem,node );
 [ node4local,elem4local ] =hole_fun_local( node_all,elem_all,P,Mxx,Mzz,M11,M22);
 [ props4global,netProp ] =hole_fun_global( var_length, node_all,elem_all,cross_section_range);
 [ elem_distor ] =hole_fun_distor( cross_section_range,elem,distor_len );
 
 %runner for g,d,l
 [curve_local{2},shapes_local{2},length_local{2}] = run_strip_local( prop,node4local{2},elem4local{2},cross_section_range{2} );
 [ curve_distor{2},shapes_distor{2}] = run_strip_distortional( prop,node,elem_distor{2},distor_len(2));
 [Pe] = cutwp_prop_4hole(node,elem,var_length,prop,props4global,P,M11,M22);
 

