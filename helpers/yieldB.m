function [By]=yieldB(fy,Cw,w)
%BWS
%October 2015
%
B=1;
stress=B*w/Cw;
factor=fy/max(abs(stress));
By=factor*B;

