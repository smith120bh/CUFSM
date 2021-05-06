load temptest2
[A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,J,Xs,Ys,Cw,B1,B2,w] = cutwp_prop2(node(:,2:3),elem(:,2:4));
[P,M11,M22,B,err] = stress_to_action(node,xcg,zcg,thetap,A,I11,I22,w,Cw);
P
M11
M22
B
err