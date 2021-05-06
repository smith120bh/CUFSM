function [ node_gross,elem_gross ] = hole_node_elem_4_local( all_coord,net_elem,P,Mxx,Mzz,M11,M22)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

node_gross=all_coord;
% node_gross=[1 5.00 1.00 1 1 1 1 33.33
%     2 5.00 0.00 1 1 1 1 50.00
%     3 2.50 0.00 1 1 1 1 50.00
%     4 0.00 0.00 1 1 1 1 50.00
%     5 0.00 3.00 1 1 1 1 16.67
%     6 0.00 6.00 1 1 1 1 -16.67
%     7 0.00 9.00 1 1 1 1 -50.00
%     8 2.50 9.00 1 1 1 1 -50.00
%     9 5.00 9.00 1 1 1 1 -50.00
%     10 5.00 8.00 1 1 1 1 -33.33];
%
%Elements

elem_gross=net_elem;
% elem_gross=[1 1 2 0.100000 100
%     2 2 3 0.100000 100
%     3 3 4 0.100000 100
%     4 4 5 0.100000 100
%     5 5 6 0.000000 100
%     6 6 7 0.100000 100
%     7 7 8 0.100000 100
%     8 8 9 0.100000 100
%     9 9 10 0.100000 100];




%-----------------------------------------------------------------
%------------TWEAKING MODEL USING OTHER CUFSM FEATURES------------
%-----------------------------------------------------------------
%Features available in the GUI may also be used in these batch programs,
%for instance the mesh in this default file is rather course, let's double
%the number of elements
% [node_gross,elem_gross]=doubler(node_gross,elem_gross);
%
%In this example the stresses (last column of nodes) are already defined,
%but it is common to use the properties page to define these values instead
%of entering in the nodal stresses. Right now this problem applies a
%reference bending moment, let's apply a reference compressive load using
%the subroutines normally used on the properties page of CUFSM
%
%first calculate the global properties
[A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,J,Xs,Ys,Cw,B1,B2,w] = cutwp_prop2(node_gross(:,2:3),elem_gross(:,2:4));
thetap=thetap*180/pi; %degrees...
Bx=NaN; By=NaN;
%
%second set the refernce stress
fy=50;
%
%third calculate the P and M associated with the reference stress
unsymmetric=0; %i.e. do a restrained bending calculation
% [P,Mxx,Mzz,M11,M22]=yieldMP(node_gross,fy,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymmetric);
%
%fourth apply just the P to the model
node_gross=stresgen(node_gross,P*1,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymmetric);
%

[rows,~]=size(node_gross);
% number=1:rows-1;
% for i=1:rows-1
%     if elem_gross(i,4)==0
%         number(i)=0;
%         
%     end
%     
% end
clear number
number=[];
count=1;
for i=1:rows-1
    if elem_gross(i,4)==0
        number(count)=i;
        count=count+1;
    end
    
end

if isempty(number)~=1
elem_gross(number,:)=[];

end

[len,~]=size(elem_gross);
elem_gross(:,1)=1:len;
number2 = intersect(node_gross(:,1),elem_gross(:,2:3));
node_gross=node_gross(number2,:);
[node_gross,elem_gross]=renumbernodes(node_gross,elem_gross);
% [len,~]=size(node_gross);
% node_gross(:,8)=zeros(len,1);

end

