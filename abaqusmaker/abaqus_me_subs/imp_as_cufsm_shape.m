function [impshape,m_a]=imp_as_cufsm_shape(filename,mag)
%imp_as_cufsm_shape
%Ben Schafer
%December 2005
%
%input
%filename of a completed analysis .mat
%mag vector of imperfection magnitudes
%
%output
%impshape in the same form as shape, but holding the imperfect geometry can
%be used for plotting the midspan summation of imperfection geomtries, this
%approximates the imperfection that will be created by cufsm_to_abaqus
%
%
load(filename)
%
if iscell(shapes)
    %     shape=shapes{j}(:,1);
    impshape=zeros(length(shapes{1}(1,1)),1);
    if sum(mag>0)>1
        msgbox(['Warning! You can only select the imperfection at you physical length'],'imperfection','help');
    else
            i=find(mag~=0);
            shape=shapes{i}(:,1); %only do first mode, shapes for now
            %norm
            shape=shape/max(abs(shape));
            %imp shape
            impshape=impshape+mag(i)*shape;
            m_a=m_all{i};
    end
    
else
    impshape=zeros(length(shapes(:,1,1)),1);
    for i=1:length(lengths)
        magi=mag(i);
        shape=shapes(:,i,1); %only do first mode, shapes for now
        %norm
        shape=shape/max(abs(shape));
        %imp shape
        impshape=impshape+magi*shape;
    end
    m_a=[1];
end
