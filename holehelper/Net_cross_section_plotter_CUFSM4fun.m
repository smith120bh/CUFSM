function [ all_coord_for_plot,all_coord,net_elem,small_piece,large_piece] = Net_cross_section_plotter_CUFSM4fun( elem,binary,ilocation,user_coordinate,hole_xy_coordinate_x,hole_xy_coordinate_z,s_coordinate,...
    cross_section_range,member_len)

%Junle Cai 2015
%BWS edits 2016

% )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes her

t=elem(1,4); %here again is the constant thickness assumption BWS 11/29/2016
matnum=elem(1,5);

gross_coordinate=[user_coordinate(:,1:3) s_coordinate];



%4 being s

[nrholes,~]=size(ilocation);
for i=1:nrholes
    
hole_point(i*2-1,1)=0;
hole_point(i*2-1,2)=hole_xy_coordinate_x{i}(1);
hole_point(i*2-1,3)=hole_xy_coordinate_z{i}(1);
hole_point(i*2-1,4)=ilocation(i,2);


hole_point(i*2,1)=0;
hole_point(i*2,2)=hole_xy_coordinate_x{i}(2);
hole_point(i*2,3)=hole_xy_coordinate_z{i}(2);
hole_point(i*2,4)=ilocation(i,3);

end

row=[];
count=1;

for i=1:nrholes
if binary(i)==0
 row(count)=i;
 count=count+1;
end
end

count=count-1;



if isempty(row)~=1
 rowrow(1:2:2*count-1)=row*2-1;
 rowrow(2:2:2*count)=row*2;
 hole_point(rowrow,:)=[];   
%ilocation(row,:)=[];
%[nrholes,~]=size(ilocation);

end

% hole_point=[];
% for i=1:nrholes
%     1
% hole_point(i*2-1,1)=0;
% hole_point(i*2-1,2)=hole_xy_coordinate_x{i}(1);
% hole_point(i*2-1,3)=hole_xy_coordinate_z{i}(1);
% hole_point(i*2-1,4)=ilocation(i,2);
% 
% 
% hole_point(i*2,1)=0;
% hole_point(i*2,2)=hole_xy_coordinate_x{i}(2);
% hole_point(i*2,3)=hole_xy_coordinate_z{i}(2);
% hole_point(i*2,4)=ilocation(i,3);
% 
% end
% 
% hole_point

if isempty(hole_point)
all_coord=gross_coordinate;
else
    all_coord=[gross_coordinate
           hole_point];
end


all_coord = sortrows(all_coord,4) ;      
 
%
% axes(axesnum)
% cla
%
%
% figure
%
flag=0;

[nr_seg,~]=size(all_coord);
for i=1:nr_seg-1
    
if all_coord(i,1)==0
    flag=flag+1;
end

if mod(flag,2)==1 || (all_coord(i,2)==all_coord(i+1,2) && all_coord(i,3)==all_coord(i+1,3))
%     plot([all_coord(i,2) all_coord(i+1,2)],[all_coord(i,3) all_coord(i+1,3)],'color','k','LineWidth',1,'marker','*')
    net_elem(i,:)=[i i i+1 0 matnum]; %material number also hardcoded BWS 
%       hold on
else
    net_elem(i,:)=[i i i+1 t matnum];
  
    
end


end
net_elem;


% xmax=max(DefPlot)
% axis auto normal
% axis equal
% axis off
% 
% axes(axes_longd)
% 
% cla
%  plot([0 member_len],[0 0],'color','k','LineWidth',0.5)
 hold on
 
 len=length(cross_section_range);
%  for i=1:len/2
%   plot([cross_section_range(2*i-1) cross_section_range(2*i)],[0 0],'color','r','LineWidth',1,'marker','square')
%  
%  
%  end
 
 hole_occupy=cross_section_range(2:2:len)-cross_section_range(1:2:len-1);
 small_piece=min(hole_occupy);
 large_piece=max(hole_occupy);
  
 all_coord_for_plot=all_coord;
 all_coord_for_plot(:,8)=zeros(nr_seg,1);
 
all_coord(:,1)=1:nr_seg;
all_coord(:,4:7)=ones(nr_seg,4);
all_coord(:,8)=zeros(nr_seg,1);





% row=[];
% count=1;
% 
% [nrNode,~]=size(all_coord);
% for i=1:nrNode-1
%     if all_coord(i,2:3)==all_coord(i+1,2:3)
%          row(count)=i;
%            count=count+1;
%     end
% end
% if isempty(row)~=1
%     all_coord(row+1,:)=[];
%     net_elem(row)=[];
% end

end