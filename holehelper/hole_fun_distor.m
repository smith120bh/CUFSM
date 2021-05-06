function [ elem_distor ] =...
    hole_fun_distor( cross_section_range,elem,distor_len,plate_range,s_coordinate,binary )


% [nrElem,~]=size(elem_local);
% 
% for i=2:nrElem
%     if elem_local(i,1)~=elem_local(i-1,2) %hole
%         start=elem_local(i-1,2);
%         last=elem_local(i,1);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%distortional
%need input of distor_len 
%elem_distor for elem, just use node for node

[~,total_exist_cs]=size(cross_section_range);

for i=2:total_exist_cs  
    [nrnodes]=length(s_coordinate{i});
    %%%%%%%  
    ct=1;
    nrhole=length(binary{i});
    for n=1:nrhole
        if binary{i}(n)==1
            plate_range_exist(ct,:)=plate_range(n,:);
            ct=ct+1; %bws addition
        end
    end
    %plate_range_exist=unique(plate_range,'rows') %commented out by BWS
    %%%%%%%
    [nr_range,~]=size(plate_range_exist);
    for m=1:nr_range 
    %     for j=1:nrnodes
    %         if s_coordinate{i}(j)<=plate_range_exist(m,1) && s_coordinate{i}(j)>=plate_range_exist(m,1)
    %             start{i}(m)=j;
    %             break
    %         end
    %     end
        tmp = abs(s_coordinate{1}-plate_range_exist(m,1)) %should always be to s_coordinate of 1 since that is the mesh used in the DB modeling
        [idx idx] = min(tmp) %index of closest value
        start{i}(m)=idx;
        %     for j=start{i}(ct):nrnodes
    %         if s_coordinate{i}(j)<=plate_range_exist(m,2) && s_coordinate{i}(j)>=plate_range_exist(m,2)
    %             last{i}(m)=j;
    %             break
    %         end 
    %     end
         tmp = abs(s_coordinate{1}-plate_range_exist(m,2)); %again changed from i to 1 since the gross section mesh not the mesh with holes is used
        [idx idx] = min(tmp); %index of closest value
        last{i}(m)=idx-1;    
    end

    for v=1:nr_range
        L_hole=max(cross_section_range{i}(2:2:end)-cross_section_range{i}(1:2:end));
        t=elem(start{i}(v),4);
        teq=(1-L_hole/distor_len)^(1/3)*t;
        elem_distor{i}=elem;
        elem_distor{i}(start{i}(v):last{i}(v),4)=teq; %had to do -1 on last in one case, not in another     
    end
    
end
elem_distor{1}=elem;
end


