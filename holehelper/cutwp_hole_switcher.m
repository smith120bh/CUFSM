function cutwp_hole_switcher(command_str)
global Pedit_hole M11edit_hole M22edit_hole
global rad_Pe rad_Me1 rad_Me2 rad_Me12 text_Pe text_Me1 text_Me2 statictext_Pe statictext_Me1 statictext_Me2 statictext_KL1 statictext_KL2 statictext_KL3 statictext_exy1
global statictext_exy2 ed_KL1 ed_KL2 ed_KL3 ed_exy ed_c text_maxmode slider rad_Me text_mode props4global rad_shear
global rad_origin rad_centroid rad_axisxy rad_axis12 rad_axial rad_deform rad_node
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if strcmp(command_str,'slider')  
    set(text_mode,'string',num2str(round(get(slider,'value'))));
    holehelper_cb(105);
elseif strcmp(command_str,'Pe')
    set(rad_Pe,'value',1);
    set(rad_Me,'value',0);
    set(rad_Me1,'value',0);
    set(rad_Me2,'value',0);
    set(rad_Me12,'value',0);
    set(text_mode,'string','1');
    set(slider,'value',1);
     holehelper_cb(105);
elseif strcmp(command_str,'Me')
    set(rad_Pe,'value',0);
    set(rad_Me,'value',1);
    set(rad_Me1,'value',1);
    set(rad_Me2,'value',0);
    set(rad_Me12,'value',0);
    set(text_mode,'string','1');
    set(slider,'value',1);
     holehelper_cb(105);
elseif strcmp(command_str,'Me1')
    set(rad_Pe,'value',0);
    set(rad_Me,'value',1);
    set(rad_Me1,'value',1);
    set(rad_Me2,'value',0);
    set(rad_Me12,'value',0);
    set(text_mode,'string','1');
    set(slider,'value',1);
     holehelper_cb(105);
elseif strcmp(command_str,'Me2')
    set(rad_Pe,'value',0);
    set(rad_Me,'value',1);
    set(rad_Me1,'value',0);
    set(rad_Me2,'value',1);
    set(rad_Me12,'value',0);
    set(text_mode,'string','1');
    set(slider,'value',1);
     holehelper_cb(105);
elseif strcmp(command_str,'Me12')
    set(rad_Pe,'value',0);
    set(rad_Me,'value',1);
    set(rad_Me1,'value',0);
    set(rad_Me2,'value',0);
    set(rad_Me12,'value',1);
    set(text_mode,'string','1');
    set(slider,'value',1);
     holehelper_cb(105);
end
end

