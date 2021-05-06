function [ node4local,elem4local ] =...
    hole_fun_local( node_all,elem_all,P,Mxx,Mzz,M11,M22)



%%%%%%%%%%%%local buckling
%also a function
%node4local: node renumbered and unneccesary nodes deleted, such that
%"strip" knows there are several segments for net cross-section

%elem4local: elem for local buckling, reference stress calculated based on
%net cross-section property

%Note for local buckling analysis: need to input half-waves lengths for
%"strip", which can be easily done using "cross_section_range". The first
%cross-section is gross cross-section.
%%%%%%%%%%%%%%%%%

[~,total_exist_cs]=size(elem_all);


for i=1:total_exist_cs
    [ node4local{i},elem4local{i} ] = hole_node_elem_4_local( node_all{i},elem_all{i},P,Mxx,Mzz,M11,M22);
end
% 
% node4local{1}=node_all{1};
% elem4local{1}=elem_all{1};

end

