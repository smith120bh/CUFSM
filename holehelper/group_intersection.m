function [ cross_section_range ] = group_intersection(n, binary, hole_range, remainder)
%   for a given vector of binary numbers, ex. [1 0 0 1 0] with "1" 
%   representing the incorporation of that hole, finds the regions where
%   those holes intersect to form a particular cross-sectional shape, using
%   a range_intersection function.

count=1;
while count<=n;
if count==1;
        if binary(1)==0;
            temp=remainder{1};
        else if binary(1)==1
                temp=hole_range{1};
            end
        end
    else
        if binary(count)==0;
            temp=range_intersection(temp, remainder{count});
        else if binary(count)==1;
                temp=range_intersection(temp, hole_range{count});
            end
        end  
          
end
  count=count+1;
end
  cross_section_range=temp;
  
  if isempty(cross_section_range)~=1
  len=length(cross_section_range);
count=1;
number=[];
  for i=1:len/2
      if cross_section_range(2*i-1)==cross_section_range(2*i)
          number(count*2-1)=2*i-1;
          number(count*2)=2*i;
          count=count+1;
      end
  end
if isempty(number)~=1
      cross_section_range(number)=[];
end

%   if cross_section_range(1:2:len-1)==cross_section_range(2:2:len)
%       cross_section_range=[];
%   end
  
  
  end
