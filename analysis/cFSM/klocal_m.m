function [k]=klocal_m(Ex,Ey,vx,vy,G,t,a,b,m,BC)

%assemble local elastic stiffness matrix for a single longitudinal term m
%
%created on Jul 10, 2009 by Z. Li

%Generate element stiffness matrix (k) in local coordinates
E1=Ex/(1-vx*vy);
E2=Ey/(1-vx*vy);
Dx=Ex*t^3/(12*(1-vx*vy));
Dy=Ey*t^3/(12*(1-vx*vy));
D1=vx*Ey*t^3/(12*(1-vx*vy));
Dxy=G*t^3/12;
%
% k=sparse(zeros(8*m,8*m));
z0=zeros(4,4);
kk=m;
nn=m;
km_mp=zeros(4,4);
kf_mp=zeros(4,4);
um=kk*pi;
up=nn*pi;
c1=um/a;
c2=up/a;
%
[I1,I2,I3,I4,I5] = BC_I1_5(BC,kk,nn,a);
%
%asemble the matrix of Km_mp
km_mp(1,1)=E1*I1/b+G*b*I5/3;
km_mp(1,2)=E2*vx*(-1/2/c2)*I3-G*I5/2/c2;
km_mp(1,3)=-E1*I1/b+G*b*I5/6;
km_mp(1,4)=E2*vx*(-1/2/c2)*I3+G*I5/2/c2;

km_mp(2,1)=E2*vx*(-1/2/c1)*I2-G*I5/2/c1;
km_mp(2,2)=E2*b*I4/3/c1/c2+G*I5/b/c1/c2;
km_mp(2,3)=E2*vx*(1/2/c1)*I2-G*I5/2/c1;
km_mp(2,4)=E2*b*I4/6/c1/c2-G*I5/b/c1/c2;

km_mp(3,1)=-E1*I1/b+G*b*I5/6;
km_mp(3,2)=E2*vx*(1/2/c2)*I3-G*I5/2/c2;
km_mp(3,3)=E1*I1/b+G*b*I5/3;
km_mp(3,4)=E2*vx*(1/2/c2)*I3+G*I5/2/c2;

km_mp(4,1)=E2*vx*(-1/2/c1)*I2+G*I5/2/c1;
km_mp(4,2)=E2*b*I4/6/c1/c2-G*I5/b/c1/c2;
km_mp(4,3)=E2*vx*(1/2/c1)*I2+G*I5/2/c1;
km_mp(4,4)=E2*b*I4/3/c1/c2+G*I5/b/c1/c2;
km_mp=km_mp*t;
%
%
%asemble the matrix of Kf_mp
kf_mp(1,1)=(5040*Dx*I1-504*b^2*D1*I2-504*b^2*D1*I3+156*b^4*Dy*I4+2016*b^2*Dxy*I5)/420/b^3;
kf_mp(1,2)=(2520*b*Dx*I1-462*b^3*D1*I2-42*b^3*D1*I3+22*b^5*Dy*I4+168*b^3*Dxy*I5)/420/b^3;
kf_mp(1,3)=(-5040*Dx*I1+504*b^2*D1*I2+504*b^2*D1*I3+54*b^4*Dy*I4-2016*b^2*Dxy*I5)/420/b^3;
kf_mp(1,4)=(2520*b*Dx*I1-42*b^3*D1*I2-42*b^3*D1*I3-13*b^5*Dy*I4+168*b^3*Dxy*I5)/420/b^3;

kf_mp(2,1)=(2520*b*Dx*I1-462*b^3*D1*I3-42*b^3*D1*I2+22*b^5*Dy*I4+168*b^3*Dxy*I5)/420/b^3;
kf_mp(2,2)=(1680*b^2*Dx*I1-56*b^4*D1*I2-56*b^4*D1*I3+4*b^6*Dy*I4+224*b^4*Dxy*I5)/420/b^3;
kf_mp(2,3)=(-2520*b*Dx*I1+42*b^3*D1*I2+42*b^3*D1*I3+13*b^5*Dy*I4-168*b^3*Dxy*I5)/420/b^3;
kf_mp(2,4)=(840*b^2*Dx*I1+14*b^4*D1*I2+14*b^4*D1*I3-3*b^6*Dy*I4-56*b^4*Dxy*I5)/420/b^3;

kf_mp(3,1)=kf_mp(1,3);
kf_mp(3,2)=kf_mp(2,3);
kf_mp(3,3)=(5040*Dx*I1-504*b^2*D1*I2-504*b^2*D1*I3+156*b^4*Dy*I4+2016*b^2*Dxy*I5)/420/b^3;
kf_mp(3,4)=(-2520*b*Dx*I1+462*b^3*D1*I2+42*b^3*D1*I3-22*b^5*Dy*I4-168*b^3*Dxy*I5)/420/b^3;

kf_mp(4,1)=kf_mp(1,4);
kf_mp(4,2)=kf_mp(2,4);
kf_mp(4,3)=(-2520*b*Dx*I1+462*b^3*D1*I3+42*b^3*D1*I2-22*b^5*Dy*I4-168*b^3*Dxy*I5)/420/b^3;%not symmetric
kf_mp(4,4)=(1680*b^2*Dx*I1-56*b^4*D1*I2-56*b^4*D1*I3+4*b^6*Dy*I4+224*b^4*Dxy*I5)/420/b^3;

%assemble the membrane and flexural stiffness matrices
kmp=[km_mp  z0
    z0  kf_mp];
%add it into local element stiffness matrix by corresponding to m
k=kmp;


       
        
        