function [node_all4plot, cross_section_range,node_all,elem_all,binary_exist ] =...
    hole_fun_general( var_length,ed_location,ed_location_x,elem,node )

%Created by Junle Cai in 2015
%Correcting for variable thickness BWS in 2016

%%%%%%%%General
%cross_section_range{i}: ith cross-section: [start end start end....]

%nodel_all{i}: node information of the ith cross-section, including the
%original nodes and hole edges   node_cross_section  node_holes

%elem_all: elem, the width at hole elements are 0's  elem_holes

%netProp{i}: cross-section properties of net cross-section

%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%local buckling
%also a function
%node4local: node renumbered and unneccesary nodes deleted, such that
%"strip" knows there are several segments for net cross-section

%elem4local: elem for local buckling, reference stress calculated based on
%net cross-section property

%Note for local buckling analysis: need to input half-waves lengths for
%"strip", which can be easily done using "cross_section_range". The first
%cross-section is gross cross-section.
%%%%%%%%%%%%%%%%%

%%%%%%%%%global buckling

%input node elem ,, out prop net gross smeared
%props4global
%props4global.J_avg=temp_J/var_length;
%props4global.Ix_avg=temp_x/var_length;
%props4global.Iy_avg=temp_y/var_length;
%props4global.I2_avg=temp_2/var_length;
%props4global.Cw_min=smallest Cw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%distortional
%need input of distor_len 
%elem_distor for elem, just use node for node

%%%%%seperate funs for local distor and global , input matrices in the
%%%%%interface, len, locations, local,d,g--spit node elem, addtional
%%%%%functions, to call l,d,g. for one c.s.


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%"calculate cross-sections" button
%determining all cross-sectional types: for cross_section{i}, outputs a 
%vector indicating the ranges where it occurs
%-------------------------------------------------------------------
%get length from "length" edit box  

% var_length= str2num(get(ed_length,'string'));
% var_elem= str2num(get(ed_elem,'string'));
var_elem=elem;
thickness=var_elem(1,4); %this is the problem as it assumes thickness is a single value BWS 11/29/2016
% 
% prop= str2num(get(ed_prop,'string'));
% assignin('base', 'var_length', var_length);
%creates node|x|z|s matrix("node_location")
%retrieves x and z coordinates of all nodes
%  node_prop= str2num(get(ed_node,'string')); 
node_prop=node;

%inserts x and z information into node_location matrix
node_location= node_prop(:,1:3);
node_location(1,4)=0;
[maxNode,~]=size(node_prop);
%calculate distance along the cross section to each node("s" coordinate)
for i=2:maxNode
    dx{i}=node_location(i,2)- node_location((i-1),2);
    dz{i}=node_location(i,3)- node_location((i-1),3);
    ds{i}=sqrt(dx{i}^2 + dz{i}^2);
    node_location(i,4)=ds{i};
end 
    
%retrieve info from "location" edit boxes    
% s = str2num(get(ed_location_x,'string')); 
% 
% 
% [row,column]=size(s);
hole_dimension=ed_location;
%%
ilocationx_char=ed_location_x;

[row,~]=size(ed_location);

for i=1:row
    hole_location{i}=str2num(ilocationx_char(i,:));
    
end
%%
hole_location
% hole_location=num2cell(s,2);%convert to cell array

%create start/end vectors for each hole ("hole_range")
for i=1:row
    column(i)=length(hole_location{i});
    matrix{i}=zeros(1,column(i)*2);
end
    hole_range=matrix.';
for i=1:row
    hole_range{i}(1,1:2:end)=hole_location{i}; %filling in "start" info
    hole_range{i}(1,2:2:end)=hole_location{i}(1,:)+hole_dimension(i,4);
    %add hole length to hole start to find hole end
end

%create "remainder" matrix
for i=1:row
    matrix{i}=zeros(1,column(i)*2+2);
    matrix{i}(1,2:end-1)=hole_range{i};
    matrix{i}(1,end)=var_length(1,1);
end
    remainder=matrix.';

%delete extra columns
for i=1:row
    min1=hole_range{i}(1);
    hole_range{i}(hole_range{i}<min1)=[];
end

for i=1:row
    min1=remainder{i}(2);
    remainder{i}(remainder{i}<min1)=[];
    remainder{i}(2:end+1)=remainder{i};
    remainder{i}(1)=0;
end
  
%call "group_intersection" function
%loop outputs all possible combinations of holes, in binary
%then implements cross_section function
n=row; %number of types of holes
count=1;
for i=0:2^n-1
    j=dec2bin(i);
    nn=length(j);
    binary{i+1}=zeros(n,1);
    binary{i+1}(n-nn+1:n,1)=bitget(i,nn:-1:1);
    temp=group_intersection(n,binary{i+1},hole_range,remainder);
    if isempty(temp)~=1
        cross_section_range{count}=temp;
        binary_exist{count}=binary{i+1};
        count=count+1;
    end
end

total_exist_cs=count-1;

temp_n=n;

% assignin('base', 'binary_exist', binary_exist);
% assignin('base', 'cross_section_range', cross_section_range);




% ilocation=str2num(get(ed_location,'String'));
ilocation=ed_location;

[hole_xy_coordinate_x,hole_xy_coordinate_z,s_coordinate ] = scoordinate2xy( ilocation,node );

for i=total_exist_cs:-1:1;
[ node_all4plot{i},node_all{i},net_elem{i},small_piece{i},large_piece{i}] =Net_cross_section_plotter_CUFSM4fun( elem,binary_exist{i},ilocation,node,hole_xy_coordinate_x,hole_xy_coordinate_z,s_coordinate,...
    cross_section_range{i},var_length);
elem_all{i}=net_elem{i};
% [node_all{i},elem_all{i}]=renumbernodes(node_all{i},elem_all{i});
end





end

