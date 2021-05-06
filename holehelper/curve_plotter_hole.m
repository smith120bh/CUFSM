function [] = curve_plotter_hole( prop,node,elem,curve,shapes,lengths,axescurve,axes2dshape_hole)

for i=1:length(lengths)
    m_all{i}=[1];
end
BC='S-S';
%cFSM features to utilize constrained finite strip method. We turned off here.
%Use of the cFSM variables are outside the scope of this batchcode example (at least for now)
%however, the following template for cFSM is ready to use with modification.
nnodes = length(node(:,1));
ndof_m= 4*nnodes;
GBTcon.ospace=1;GBTcon.couple=1;GBTcon.orth=2;GBTcon.norm=1; %see strip for possible solutions

[elprop,m_node,m_elem,node_prop,nmno,ncno,nsno,ndm,nlm,DOFperm]=base_properties(node,elem);

ngm=4;nom=2*(length(node(:,1))-1);
GBTcon.local=zeros(1,nlm);
GBTcon.dist=zeros(1,ndm);
GBTcon.glob=zeros(1,ngm); %can be less than 4 for special sections
GBTcon.other=zeros(1,nom);
%to turn on any cFSM classification use 1's on the GBTcon vectors defined
%above. Leave off unless you want it.
%for example say you wanted distortional only->GBTcon.dist=ones(1,ndm);
%Set default initial values
springs=0;

%In this case we have no constraint
% if any, constraint should be defined in this format:
% constraints=[node#e DOFe coeff. noder#k DOFk
%               ...   ...   ...    ...    ... ];
constraints=0;

%set the eigenmode you want to output, optional
neigs=10; %GUI default is 20
clas=0;
i=1;
filenamecell{i}=['CUFSM results'];
propcell{i}=prop;
nodecell{i}=node;
elemcell{i}=elem;
lengthscell{i}=lengths;
curvecell{i}=curve;
shapescell{i}=shapes;
assignin('base','shapesl',shapes);
clascell{i}=clas;
springscell{i}=springs;
constraintscell{i}=constraints;
GBTconcell{i}=GBTcon;
BCcell{i}=BC;
m_allcell{i}=m_all;

fileindex=1;
lengthindex=1;
modeindex=1;
clasopt=0;
logopt=1;
minopt=1;
files=1;
modes=(1:1:length(curve{lengthindex}(:,2)));
modedisplay=1;
filedisplay=files;

lengthindex=ceil(length(lengths)/2);
modeindex=1;
clasopt=0;
logopt=1;
minopt=1;
xmin=min(lengths)*10/11;
ymin=0;
xmax=max(lengths)*11/10;
for j=1:max(size(curve));
    curve_sign(j,1)=curve{j}(modeindex,1);
    curve_sign(j,2)=curve{j}(modeindex,2);
end
ymax=min([max(curve_sign(:,2)),3*median(curve_sign(:,2))]);

modes=(1:1:length(curve{lengthindex}(:,2)));
modedisplay=1;
% filedisplay=files;

threed=0;
undef=1;
scale=1;
SurfPos=1/2;
curveoption=1;
ifcheck3d=0;
ifpatch=0;
% m_a=m_all{lengthindex};

if exist('springs')
else
    springs=0;
end
%
%Plot the first mode shape as a start
% mode=shapes{lengthindex}(:,modeindex);
% if threed==1
%     dispshap(undef,node,elem,mode,axes2dshape,scale,springs,m_a,BC,SurfPos);
%     dispshp2(undef,lengths(lengthindex),node,elem,mode,axes3dshape,scale,m_a,BC,ifpatch);
% else
%     dispshap(undef,node,elem,mode,axes2dshapelarge,scale,springs,m_a,BC,SurfPos);
% end
%Plot the buckling curve as a start
if length(lengths)==1
    curveoption=2;
end
if curveoption==1
    xmin=min(lengths)*10/11;
    ymin=0;
    xmax=max(lengths)*11/10;
    for j=1:max(size(curve));
        curve_sign(j,1)=curve{j}(modeindex,1);
        curve_sign(j,2)=curve{j}(modeindex,2);
    end
    ymax=min([max(curve_sign(:,2)),3*median(curve_sign(:,2))]);
elseif curveoption==2
    xmin=1;
    ymin=0;
    xmax=length(curve{lengthindex}(:,2));
    ymax=min([max(curve{lengthindex}(:,2)),3*median(curve{lengthindex}(:,2))]);
end

if curveoption==1
    picpoint=[curve{lengthindex}(modeindex,1) curve{lengthindex}(modeindex,2)];
    thecurve3(curvecell,filenamecell,clascell,filedisplay,minopt,logopt,clasopt,axescurve,xmin,xmax,ymin,ymax,modedisplay,fileindex,modeindex,picpoint)
else
    axes(axescurve);
    cla reset, axis off
    %
%     set(axescurvemode,'visible','on');
    
    %plot the load factor vs mode number curve
    picpoint=[modes(modeindex) curve{lengthindex}(modeindex,2)];
    thecurve3mode(curvecell,filenamecell,clascell,fileindex,minopt,logopt,clasopt,axescurvemode,xmin,xmax,ymin,ymax,fileindex,lengthindex,picpoint);
     axes(axes2dshape_hole);


  %plot buckling mode shape
% change this to plot buckling mode shape at specified length and specified
% mode
lengthindex=1;
modeindex=1;

%additional initial input
threed=0;%whether plot 3d mode shape
undef=1;
scale=-1;
SurfPos=1/2;
curveoption=1;
ifcheck3d=0;
ifpatch=0;
m_a=m_all{lengthindex};

%

%Plot the first mode shape as a start
mode=shapes{lengthindex}(:,modeindex);
dispshap(undef,node,elem,mode,axes2dshape_hole,scale,springs,m_a,BC,SurfPos);
if threed==1
%     figure(3)
%     axes3dshape=axes('Units','normalized','Position',[0.1 0.1 0.8 0.8],'visible','off');
  cla axes3dshape
dispshp2(undef,lengths(lengthindex),node,elem,mode,axes3dshape,scale,m_a,BC,ifpatch);
end
    
    %plot the participation vs longitudinal terms
%     mode=shapes{lengthindex}(:,modeindex);
%     nnodes=length(node(:,1));
%     m_a=m_all{lengthindex};
%     [d_part]=longtermpart(nnodes,mode,m_a);
%     axes(axesparticipation)
%     cla
%     bar(m_a,d_part);hold on
%     xlabel('m, longitudinal term');hold on;
%     ylabel('Participation');hold on;
%     legendstring=[filenamecell{fileindex},', length = ',num2str(lengths(lengthindex)), ', mode = ', num2str(modeindex)];
%     hlegend=legend(legendstring);hold on
%     set(hlegend,'Location','best');
%     hold off    
end   