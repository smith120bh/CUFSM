function [all_coord,sub_startend]=Net_cross_section_finder4fun( elem,binary_exist,ilocation,user_coordinate,hole_xy_coordinate_x,hole_xy_coordinate_z,s_coordinate,...
    total_exist_cs)

% )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes her
% global elem
t=elem(1,4);

gross_coordinate=[user_coordinate(:,1:3) s_coordinate];
[nrnode,~]=size(gross_coordinate);

gross_coordinate(1:nrnode,5)=1;  %5th row, not hole


%4 being s

[nrholes,~]=size(ilocation);
for i=1:nrholes
    
hole_point(i*2-1,1)=0;
hole_point(i*2-1,2)=hole_xy_coordinate_x{i}(1);
hole_point(i*2-1,3)=hole_xy_coordinate_z{i}(1);
hole_point(i*2-1,4)=ilocation(i,2);
hole_point(i*2-1,5)=0;

hole_point(i*2,1)=0;
hole_point(i*2,2)=hole_xy_coordinate_x{i}(2);
hole_point(i*2,3)=hole_xy_coordinate_z{i}(2);
hole_point(i*2,4)=ilocation(i,3);
hole_point(i*2,5)=0;

end

ct=[];
count=1;
for i=1:nrnode

    if any(hole_point(:,4)==gross_coordinate(i,4))
        ct(count)=i;
        count=count+1;
    end
end

if isempty(ct)~=1
    gross_coordinate(ct,:)=[];
end

assignin('base', 'gross_coordinate',gross_coordinate);


hole_point
all_coord=[gross_coordinate
           hole_point];
       
[all_coord,hole_row] = sortrows(all_coord,[4 5]) ;


[all_node,~]=size(all_coord);

% all_coord(1:end,1)=1:all_node;

% hole_row
% hole_row=hole_row(all_node-2*nrholes+1:all_node)
hole_row=find(all_coord(:,1)==0);
 count=1;
for n=1:total_exist_cs
    total=1:all_node;
    temp=total(binary_exist{n}==1);
    temp1=2*temp-1;
    temp2=2*temp;
    temp=sort([temp1 temp2]);
 hole_number=hole_row(temp);
 total(hole_number)=-1
 
 flag=0;
count=1;
 if total(1)~=-1
 sub_startend{n}(count,1)=1;
 flag=1;
 for j=2:all_node
     if (total(j)==-1 || j==all_node) && flag==1
         sub_startend{n}(count,2)=j;
         flag=0;
         count=count+1;
     elseif total(j)==-1 && j~=all_node && flag==0
        sub_startend{n}(count,1)=j;
         flag=1;
     end
     
 end
 end
 
 if total(1)==-1
     for j=2:all_node
      if (total(j)==-1) && flag==0
         sub_startend{n}(count,1)=j;
         flag=1;
      elseif (total(j)==-1 || j==all_node) && flag==1
        sub_startend{n}(count,2)=j;
         flag=0;
         count=count+1;
      end

      
      
      
     end
     
 end
 
%  for j=1:all_node
%      
%     
%   if total(j)~=-1 && flag==0
%  sub_startend{n}(count,1)=j;
%  flag=1;
%   end
%   if total(j)>0 && total(j+1)==-1 && flag==1  %end of one segment
%   flag=0;
%   sub_startend{n}(count,2)=j;
%   count=count+1;
%   end
% %   sub_startend{n,:}=SplitVec(total,'cons');
%  end
%  
 end
 
end