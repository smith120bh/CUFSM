function []=abaqus_me_cb(num)
%
global pathname filename prop node elem lengths curve shapes clas springs constraints GBTcon BC m_all elemtype elemtp nelbsp nelbct
global axessect flags filename_txt ed_filename ed_impmag ed_length ed_numele txt_msg buckling static ed_elemtype txt_nelbsp ed_nelbsp txt_nelbct ed_nelbct
%
switch num
%
case 1
%-------------------------------------------------------------------
%code for the load button
%standard load
[pathname,filename,prop,node,elem,lengths,curve,shapes,springs,constraints,GBTcon,clas,BC,m_all]=loader; %need to change to add BC, m_all
crossect(node,elem,axessect,springs,constraints,flags);
set(filename_txt,'String',[pathname,filename]);
savefile=filename;
savefile(length(filename)-3:length(filename))=['.inp'];
set(ed_filename,'String',[pathname,savefile]);
mag=zeros(length(lengths),1);

if iscell(curve)
    nlength=max(size(curve));    
    if nlength<length(lengths)
        set(txt_msg,'String',['Warning! CUFSM file does not seem to have results, can''t generate imperfections...']);
        return
    end
    for i=1:nlength
        curve1(i,:)=curve{i}(1,:);
    end
    set(ed_impmag,'String',sprintf('%f %f %f\n',[curve1(:,1:2) mag]'));
else
    if length(curve(:,1))<length(lengths)
        set(txt_msg,'String',['Warning! CUFSM file does not seem to have results, can''t generate imperfections...']);
        return
    end
    set(ed_impmag,'String',sprintf('%f %f %f\n',[curve(:,1:2) mag]'));
end

%-------------------------------------------------------------------
%
case 2
%-------------------------------------------------------------------
%code for the generate ABAQUS input file button
L=str2num(get(ed_length,'String'));
nele=str2num(get(ed_numele,'String'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define element type
elemtp=str2num(get(ed_elemtype,'String'));
nelbsp=str2num(get(ed_nelbsp,'String'));
nelbct=str2num(get(ed_nelbct,'String'));
if elemtp==1
    elemtype='S9R5';
elseif elemtp==2
    elemtype='S4';
elseif elemtp==3
    elemtype='S4R';
else
    ['Warning! Your element type must be one of S9R5, S4, and S4R!']
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curveandmag=str2num(get(ed_impmag,'String'));
mag=curveandmag(:,3);
if length(mag)~=length(lengths)
    set(txt_msg,'String',['Warning! Imperfection magnitudes not read in properly, reload CUFSM file and try once more']);
    return
end
if rem(nele,2)>0
    set(txt_msg,'String',['Warning! Number of elements should be an even number for boundary conditions to be correct']);
end
if get(buckling,'Value')==1
    analysis=1;
else
    analysis=2;
end

loadfile=filename(1:length(filename)-4);
savefile=get(ed_filename,'String');
savefile=savefile(1:length(savefile)-4);
cufsm_to_abaqus([pathname,loadfile],savefile,L,nele,mag,analysis,elemtype,nelbsp,nelbct)
set(txt_msg,'String',[savefile,'.inp Generated']);
%
%
case 3
%-------------------------------------------------------------------
%set analysis to buckling, i.e., 1
set(buckling,'Value',1);
set(static,'Value',0);
%-------------------------------------------------------------------
%
case 4
%-------------------------------------------------------------------
%set analysis type to static, i.e., 2
set(buckling,'Value',0);
set(static,'Value',1);
%-------------------------------------------------------------------
%
case 5
%-------------------------------------------------------------------
%plot the imperfection selected
loadfile=filename(1:length(filename)-4);
curveandmag=str2num(get(ed_impmag,'String'));
mag=curveandmag(:,3);
if length(mag)~=length(lengths)
    set(txt_msg,'String',['Warning! Imperfection magnitudes not read in properly, reload CUFSM file and try once more']);
    return
end
[impshape,m_a]=imp_as_cufsm_shape([pathname,loadfile],mag);
xborder=(max(node(:,2))-min(node(:,2)))*0.2;
yborder=(max(node(:,3))-min(node(:,3)))*0.2;
SurfPos=1/2;
dispshap_absolutescale(1,node,elem,impshape,axessect,1.0,0,[min(node(:,2))-xborder max(node(:,2))+xborder min(node(:,3))-yborder max(node(:,3))+yborder],m_a,BC,SurfPos);
%[min(node(:,2)) max(node(:,2)) min(node(:,3)) max(node(:,3))]
%-------------------------------------------------------------------

        
end