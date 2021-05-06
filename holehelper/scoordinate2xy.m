function [hole_xy_coordinate_x,hole_xy_coordinate_z,s_coordinate ] = scoordinate2xy( ilocation,user_coordinate )

startpoint=ilocation(:,2);
endpoint = ilocation(:,3);

[MaxS,~]=size(startpoint);
s_coordinate  = 0;
[MaxCoordinate,~]=size(user_coordinate);
for i=2:MaxCoordinate;
    dist_x=user_coordinate(i,2)-user_coordinate((i-1),2);
    dist_z  =user_coordinate(i,3)-user_coordinate((i-1),3); 
    dist_s =  sqrt(dist_x^2+dist_z^2);
    s_coordinate(i,1)= s_coordinate(i-1) + dist_s;
end
   
for i = 1:MaxS;
    counting_start=ilocation(i,2);
    counting_end  =ilocation(i,3); 
    
  for j=1:MaxCoordinate-1  
    if counting_start >= s_coordinate(j) && counting_start <= s_coordinate(j+1);
        As=s_coordinate(j);
        Bs=s_coordinate(j+1);
        
        d_s=counting_start-As;
        before_x=user_coordinate(j,2);
        after_x=user_coordinate(j+1,2);
        
         before_z=user_coordinate(j,3);
        after_z=user_coordinate(j+1,3);
        
        x_start=before_x+d_s/(Bs-As)*(after_x-before_x);
        z_start=before_z+d_s/(Bs-As)*(after_z-before_z);
    end
 

    
    
    if counting_end >= s_coordinate(j) && counting_end <= s_coordinate(j+1);
        Ae=s_coordinate(j);
        Be=s_coordinate(j+1);

        d_s=counting_end-Ae;
        before_x=user_coordinate(j,2);
        after_x=user_coordinate(j+1,2);
        
         before_z=user_coordinate(j,3);
        after_z=user_coordinate(j+1,3);
        
        x_end=before_x+d_s/(Be-Ae)*(after_x-before_x);
        z_end=before_z+d_s/(Be-Ae)*(after_z-before_z);
        break
        
    end
   
   
  end
  
  
     hole_xy_coordinate_x{i} = [x_start x_end];
     hole_xy_coordinate_z{i} = [z_start z_end];
     
    
end
   
end

