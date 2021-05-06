function [fib_c,nodef_c,node_c,element_c,node_c90,element_c90,A,xcg,zcg,Ixx,Izz,Ixz,thetapf,I11,I22]=...
    PMM_fiber_node_element_Z_rot_general(section,Dim,n_elem_cufsm,n_elem_plastic,NL,CorZ,Round_corner)
%%% In the name of Allah
%%% Fiber element Model
%%% March 2013
%%% Adding Z; November 2013
%%% Shahabeddin Torabian
%%% Generalizing the shape; May 2015
%%% Shahabeddin Torabian


% Reading SSMA Section Properties
% [ndata,text,alldata]=xlsread('SSMA1.xlsx');

% % Number of Elements for CUFSM stability analysis
% n_elem_cufsm=[4 4 4 4 8 4 4 4 4];
% 
% % Number of Elements and fiber layers for Plastic Surface analysis
% n_elem_plastic=3*[4 4 4 4 16 4 4 4 4];
% NL_plastic=7;
% NL=7;
% CorZ=2;


% Assigning dimensions

% Dim=[6.272	1.747	1.865	0.587	0.651	0.314	0.251	0.364	0.285	88.716	88.759	42.210	42.157 0.05873];

h=Dim(1);
b1=Dim(2);
b2=Dim(3);
d1=Dim(4);
d2=Dim(5);
r1=Dim(6);
r3=Dim(7);
r2=Dim(8);
r4=Dim(9);
qf1=Dim(10);
qf2=Dim(11);
ql1=Dim(12);
ql2=Dim(13);
t=Dim(14);


% for 90  deg.
        h_90=h+r1*tan((qf1/2)/180*pi)+r3*tan((qf2/2)/180*pi);
        b1_90=b1+r1*tan((qf1/2)/180*pi)+r2*tan((ql1/2)/180*pi);
        b2_90=b2+r3*tan((qf2/2)/180*pi)+r4*tan((ql2/2)/180*pi);
        d1_90=d1+r2*tan((ql1/2)/180*pi);
        d2_90=d2+r4*tan((ql2/2)/180*pi);
        qf1_90=qf1;
        qf2_90=qf2;
        ql1_90=ql1;
        ql2_90=ql2;


kipin=1;

% if CorZ==1
%     theta=90;
% else 
%     theta=ndata(I,21);
% end



% n_elem_plastic=3*[4 4 4 4 16 4 4 4 4];

% Number of layers (Even Numbers 1,3,5,7,9,11,...)

% NL=5;
Dt=t/NL;

% Preparing templatecalc.m input data
if Round_corner==0;
    r1=0;r2=0;r3=0;r4=0;
end

if r1==0&r2==0&r3==0&r4==0;
    
    if CorZ==1
        for i=1:NL      
            hh(i)=h+(i-(NL+1)/2)*2*Dt;
            bb1(i)=b1+(i-(NL+1)/2)*2*Dt;
            bb2(i)=b2+(i-(NL+1)/2)*2*Dt;
            dd1(i)=d1+(i-(NL+1)/2)*Dt;
            dd2(i)=d2+(i-(NL+1)/2)*Dt;
            rr1(i)=0;
            rr2(i)=0;
            rr3(i)=0;
            rr4(i)=0;
            qqf1(i)=qf1;
            qqf2(i)=qf2;
            qql1(i)=ql1;
            qql2(i)=ql2;           
            tt(i)=Dt;    
        end
    else
        for i=1:NL 
            hh(i)=h-t;
            bb1(i)=b1_90-(i-(NL+1)/2)*Dt-(i-(NL+1)/2)*Dt*tan((theta/2)/180*pi);
            bb2(i)=b2_90+(i-(NL+1)/2)*Dt+(i-(NL+1)/2)*Dt*tan((theta/2)/180*pi);
            dd1(i)=d1_90-(i-(NL+1)/2)*Dt*tan((theta/2)/180*pi);
            dd2(i)=d2_90+(i-(NL+1)/2)*Dt*tan((theta/2)/180*pi);           
            rr1(i)=0;
            rr2(i)=0;
            rr3(i)=0;
            rr4(i)=0;
            qqf1(i)=qf1;
            qqf2(i)=qf2;
            qql1(i)=ql1;
            qql2(i)=ql2;
            tt(i)=Dt;
        end
    end
        
    
else
   
    for i=1:NL        
        hh(i)=h;  
        if CorZ==1
        rr1(i)=(i-(NL+1)/2)*Dt+r1;
        rr2(i)=(i-(NL+1)/2)*Dt+r2;
        else
        rr1(i)=((NL+1)/2-i)*Dt+r1;
        rr2(i)=((NL+1)/2-i)*Dt+r2;
        end
        rr3(i)=(i-(NL+1)/2)*Dt+r3;
        rr4(i)=(i-(NL+1)/2)*Dt+r4;
        bb1(i)=b1;
        bb2(i)=b2;
        dd1(i)=d1;
        dd2(i)=d2;
        qqf1(i)=qf1;
        qqf2(i)=qf2;
        qql1(i)=ql1;
        qql2(i)=ql2;
        tt(i)=Dt;

    end

end

[prop1,node_c,element_c,lengths1,springs1,constraints1,geom1,cz1]=...
    templatecalc2_g(CorZ,h,b1,b2,d1,d2,r1,r2,r3,r4,qf1,qf2,ql1,ql2,t,kipin,n_elem_cufsm);

[A,xcg,zcg,Ixx,Izz,Ixz,thetap_c,I11,I22]=grosprop(node_c,element_c);

node_c(:,2)=node_c(:,2)-xcg;
node_c(:,3)=node_c(:,3)-zcg;

R = [cos(-thetap_c/180*pi)  -sin(-thetap_c/180*pi) ; sin(-thetap_c/180*pi)  cos(-thetap_c/180*pi)];
node_c_rot=node_c;
Dim_r=R*node_c(:,2:3)';
node_c_rot(:,2:3)=Dim_r';
node_c=node_c_rot;
           
[prop1,node_c90,element_c90,lengths1,springs1,constraints1,geom1,cz1]=...
    templatecalc2_g(CorZ,h_90,b1_90,b2_90,d1_90,d2_90,0,0,0,0,qf1_90,qf2_90,ql1_90,ql2_90,t,kipin,n_elem_cufsm);

[A,xcg,zcg,Ixx,Izz,Ixz,thetap_c90,I11,I22]=grosprop(node_c90,element_c90);

node_c90(:,2)=node_c90(:,2)-xcg;
node_c90(:,3)=node_c90(:,3)-zcg;

R_90 = [cos(-thetap_c90/180*pi)  -sin(-thetap_c90/180*pi) ; sin(-thetap_c90/180*pi)  cos(-thetap_c90/180*pi)];
node_c90_rot=node_c90;
Dim_r90=R*node_c90(:,2:3)';
node_c90_rot(:,2:3)=Dim_r90';
node_c90=node_c90_rot;

for i=1:NL

[prop1,node1,elem1,lengths1,springs1,constraints1,geom1,cz1]=...
    templatecalc2_g(CorZ,hh(i),bb1(i),bb2(i),dd1(i),dd2(i),rr1(i),rr2(i),rr3(i),rr4(i),qqf1(i),qqf2(i),qql1(i),qql2(i),tt(i),kipin,n_elem_plastic);
NN=length(node1);
NE=length(elem1);
nodeR(:,1)=node1(:,1)+NN*(i-1);
if CorZ==1
corrx=node_c(1,2)+(i-(NL+1)/2)*Dt*sin(ql1/180*pi()+qf1/180*pi()-pi()/2)-node1(1,2);
corry=node_c(1,3)-(i-(NL+1)/2)*Dt*cos(ql1/180*pi()+qf1/180*pi()-pi()/2)-node1(1,3);
else
corrx=node_c(1,2)-(i-(NL+1)/2)*Dt*sin(ql1/180*pi()+qf1/180*pi()-pi()/2)-node1(1,2);
corry=node_c(1,3)+(i-(NL+1)/2)*Dt*cos(ql1/180*pi()+qf1/180*pi()-pi()/2)-node1(1,3);
end
nodeR(:,2)=node1(:,2)+corrx;
nodeR(:,3)=node1(:,3)+corry;
nodeR(:,4:8)=node1(:,4:8);

elemR(:,1)=elem1(:,1)+NE*(i-1);
elemR(:,2:3)=elem1(:,2:3)+NN*(i-1);
elemR(:,4)=Dt;
elemR(:,5)=elem1(:,5);

if i==1
node=nodeR;
elem=elemR;
springs=springs1;
constraints=constraints1;
else
node=[node;nodeR];
elem=[elem;elemR]; 
springs=[springs;springs1];
constraints=[constraints;constraints1];
end

end


%%% Section Properties and Fiber generation
[A,xcg,zcg,Ixx,Izz,Ixz,thetapf,I11,I22,fib]=grosprop1(node,elem);

axestemp=1;
cz=1;

flags=[0,0,0,0,0,0,0,0,0];
axesnum=1;


NF=length(fib);

AX=0;
Az=0;

for i=1:NF
 
  AX=AX+fib(i,2)*fib(i,4);
  Az=Az+fib(i,3)*fib(i,4);
  
       
end

Af=sum(fib(:,4));
xgf=AX/Af;
zgf=Az/Af;


%%% Moving center to Centroid
%fib:[fiber# x_f z_f A_f]

fib_c=fib;
fib_c(:,2)=fib(:,2)-xcg;
fib_c(:,3)=fib(:,3)-zcg;

R = [cos(-thetapf/180*pi)  -sin(-thetapf/180*pi) ; sin(-thetapf/180*pi)  cos(-thetapf/180*pi)];

fib_c_rot=fib_c;
Dim_r=R*fib_c(:,2:3)';
fib_c_rot(:,2:3)=Dim_r';
fib_c=fib_c_rot;

nodef_c=node;
nodef_c(:,2)=node(:,2)-xcg;
nodef_c(:,3)=node(:,3)-zcg;

nodef_c_rot=nodef_c;
Dim_r=R*nodef_c(:,2:3)';
nodef_c_rot(:,2:3)=Dim_r';
nodef_c=nodef_c_rot;

% for rotated Z
FIx=Ixx;
FIz=Izz;
Ixx=I11;
Izz=I22;
I11=FIx;
I22=FIz;
Ixz=0;


end