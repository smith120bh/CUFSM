function [ ] = Net_cross_section_plotter4500(...
    axesnum,axes_longd,cross_section_range,member_len,all_coord)

% )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes her

% t=elem(1,4);
% 
% gross_coordinate=[user_coordinate(:,1:3) s_coordinate];
% 
% 
% 
% %4 being s
% 
% [nrholes,~]=size(ilocation);
% for i=1:nrholes
%     
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
% row=[];
% count=1;
% 
% for i=1:nrholes
% if binary(i)==0
%  row(count)=i;
%  count=count+1;
% end
% end
% 
% count=count-1;
% 
% 
% 
% if isempty(row)~=1
%  rowrow(1:2:2*count-1)=row*2-1;
%  rowrow(2:2:2*count)=row*2;
%  hole_point(rowrow,:)=[];   
% end
% 
% if isempty(hole_point)
% all_coord=gross_coordinate;
% else
%     all_coord=[gross_coordinate
%            hole_point];
% end
% 
% 
% all_coord = sortrows(all_coord,4) ;      
 
%
axes(axesnum)
cla
%
%
%  figure
%
flag=0;

[nr_seg,~]=size(all_coord);

for i=1:nr_seg-1
    
if all_coord(i,1)==0
    flag=flag+1;
end

if mod(flag,2)==0
    plot([all_coord(i,2) all_coord(i+1,2)],[all_coord(i,3) all_coord(i+1,3)],'color','k','LineWidth',1,'marker','*')

      hold on
else

  
    
end


end



% xmax=max(DefPlot)
axis auto normal
axis equal
axis off

%%
axes(axes_longd)

cla
 plot([0 member_len],[0 0],'color','k','LineWidth',3)
 hold on
 
 len=length(cross_section_range);
 for i=1:len/2
  plot([cross_section_range(2*i-1) cross_section_range(2*i)],[0 0],'color','r','LineWidth',3,'marker','square')
 
 
 end
  axis on
  set(gca,'ytick',[])
%  
%  hole_occupy=cross_section_range(2:2:len)-cross_section_range(1:2:len-1);
%  small_piece=min(hole_occupy);
%  large_piece=max(hole_occupy);
%  
% all_coord(:,1)=1:nr_seg;
% all_coord(:,4:7)=ones(nr_seg,4);


end

