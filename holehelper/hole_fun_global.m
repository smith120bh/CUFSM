function [ props4global,netProp ] =...
    hole_fun_global( var_length, node_all,elem_all,cross_section_range)


%%%%%%%%%global buckling

%input node elem ,, out prop net gross smeared
%props4global
%props4global.J_avg=temp_J/var_length;
%props4global.Ix_avg=temp_x/var_length;
%props4global.Iy_avg=temp_y/var_length;
%props4global.I2_avg=temp_2/var_length;
%props4global.Cw_min=smallest Cw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,total_exist_cs]=size(elem_all);

%calculate avg I
temp_A=0;
temp_x=0;
temp_y=0;
temp_1=0;
temp_2=0;
temp_J=0;
temp_xz=0;
temp_Cw=0;
for i=1:total_exist_cs
    
    cs_seg_length=cross_section_range{i}(2:2:end)-cross_section_range{i}(1:2:end);
    cs_total_length=sum(cs_seg_length);
    
    if i==1
node_gross=node_all{i};
elem_gross=elem_all{i};


    else
node_gross=node_all{i};
elem_gross=elem_all{i}; 
%         
% node_gross=node4local{i};
%  elem_gross=elem4local{i}; 
    end
    [nrelem,~]=size(elem_all{i});

a=(elem_gross(:,4)==0);
ct=0;
for m=2:nrelem
    if a(m)==1 && a(m-1)==1
    ct=ct+1;
    mark(ct)=m
    end
end

if ct~=0
    elem_gross(mark,:)=[]; 
    [nrelem,~]=size(elem_gross);
    elem_gross(:,2)=1:nrelem;
    elem_gross(:,3)=2:nrelem+1;
    node_gross(mark,:)=[];
end





% assignin('base', 'node_gross', node_gross );
% assignin('base', 'elem_gross', elem_gross );

% [total,~]=sizenode_gross
% for n=1:total
%     if node_gross(n,2)==node_gross(n+1,2) && node_gross(n,2)==node_gross(n+1,2)
% end

% [A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,J,Xs,Ys,Cw,B1,B2,w] = cutwp_prop2_for_hole(node_gross(:,2:3),elem_gross(:,2:4));

[A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,J,Xs,Ys,Cw,B1,B2,w] = cutwp_prop2(node_gross(:,2:3),elem_gross(:,2:4));

% if i==1
%    length_gross=cs_total_length;
% end
% if i~=1
%    temp_Cw=temp_Cw+cs_total_length*Cw; 
%     
% end
if i==1
    Cw_temp=Cw;
           netProp{i}.A=A;
    netProp{i}.Ixx=Ixx;
    netProp{i}.Izz=Izz;
    netProp{i}.I11=I11;
    netProp{i}.I22=I22;
    netProp{i}.J=J;
    netProp{i}.xcg=xcg; 
    netProp{i}.zcg=zcg; 
    netProp{i}.Ixz=Ixz; 
    netProp{i}.thetap=thetap; 
    netProp{i}.Xs=Xs; 
    netProp{i}.Ys=Ys; 
    netProp{i}.B1=B1; 
    netProp{i}.B2=B2; 
    netProp{i}.w=w; 
    netProp{i}.Cw=Cw;
end

if i~=1

    if Cw<Cw_temp
        Cw_temp=Cw;
    end
           netProp{i}.A=A;
    netProp{i}.Ixx=Ixx;
    netProp{i}.Izz=Izz;
    netProp{i}.I11=I11;
    netProp{i}.I22=I22;
    netProp{i}.J=J;
    netProp{i}.xcg=xcg; 
    netProp{i}.zcg=zcg; 
    netProp{i}.Ixz=Ixz; 
    netProp{i}.thetap=thetap; 
    netProp{i}.Xs=Xs; 
    netProp{i}.Ys=Ys; 
    netProp{i}.B1=B1; 
    netProp{i}.B2=B2; 
    netProp{i}.w=w; 
    netProp{i}.Cw=Cw; 
end
temp_A=temp_A+cs_total_length*A;
temp_x=temp_x+cs_total_length*Ixx;
temp_y=temp_y+cs_total_length*Izz;
temp_xz=temp_xz+cs_total_length*Ixz;
temp_1=temp_1+cs_total_length*I11;
temp_2=temp_2+cs_total_length*I22;
temp_J=temp_J+cs_total_length*J;

end
props4global.A=temp_A/var_length;
props4global.J=temp_J/var_length;
props4global.Ixx=temp_x/var_length;
props4global.Izz=temp_y/var_length;
props4global.Ixz=temp_xz/var_length;
props4global.I11=temp_1/var_length;
props4global.I22=temp_2/var_length;
props4global.Cw=Cw_temp;





end

