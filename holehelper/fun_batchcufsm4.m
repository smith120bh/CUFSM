function [curve,curvecell,filenamecell,clascell,filedisplay,minopt,logopt,clasopt,axescurve,xmin,xmax,ymin,ymax,modedisplay,fileindex,modeindex,picpoint...
    ,lengths,lengthindex,shapes,shapescell...
    undef,node,elem,mode,axes2dshape,scale,springs,m_a,BC,SurfPos] = fun_batchcufsm4( prop,node,elem,number,small_piece,large_piece,axes_handle,axes2dshape )
%-----------------------------------------------------------------
%----------------additional input definitions---------------------
%-----------------------------------------------------------------

%Lengths:
%for signiture curve, the length is interpreted as half-wave length the same as older cufsm version (3.x and previous)
lengths=linspace(large_piece/20,large_piece+small_piece/2,20);
% %for general boundary conditions, the length is interpreted as physical member length
% lengths=[92 192];

%Boundary conditions: a string specifying boundary conditions to be analyzed
%'S-S' simply-pimply supported boundary condition at loaded edges
%'C-C' clamped-clamped boundary condition at loaded edges
%'S-C' simply-clamped supported boundary condition at loaded edges
%'C-F' clamped-free supported boundary condition at loaded edges
%'C-G' clamped-guided supported boundary condition at loaded edges
BC='S-S';
%note, the signiture curve soultion is only meaningful when BC is S-S. Be sure to set your BC as S-S when performing signiture curve analysis.

%Longitudinal terms in cell: associated with each length
%for signiture curve: longitudinal term is defaultly 1 for each length
for i=1:length(lengths)
    m_all{i}=[1];
end
%for general boundary conditions
% multiple longitudinal terms are defined
% the longitudinal terms for each length can be different
% for each length, longitudinal terms are not necessarily consecutive as shown here
% (always recommend) the recommended longitudinal terms can be defined by calling:
% [m_all]=m_recommend(prop,node,elem,lengths);
% careful here, you likely want m's around L/Lcre, L/Lcrd, L/Lcrl.. maybe be
% quite high in terms of m, and usually much less terms than brute force..
% for i=1:length(lengths)
%     m_all{i}=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
% end

%In this case we have no springs
%if any, springs should be defined in this format (kflag: 1 or 0):
%springs=[node# dof stiffness kflag
%          ...  ...    ...     ... ];
springs=0;

%In this case we have no constraint
% if any, constraint should be defined in this format:
% constraints=[node#e DOFe coeff. noder#k DOFk
%               ...   ...   ...    ...    ... ];
constraints=0;

%set the eigenmode you want to output, optional
neigs=10; %GUI default is 20


%cFSM features to utilize constrained finite strip method. We turned off here.
%Use of the cFSM variables are outside the scope of this batchcode example (at least for now)
%however, the following template for cFSM is ready to use with modification.
nnodes = length(node(:,1));
ndof_m= 4*nnodes;
GBTcon.ospace=1;GBTcon.couple=1;GBTcon.orth=2;GBTcon.norm=1; %see strip for possible solutions

node
elem
[elprop,m_node,m_elem,node_prop,nmno,ncno,nsno,ndm,nlm,DOFperm]=base_properties(node,elem);

ngm=4;nom=2*(length(node(:,1))-1);
GBTcon.local=zeros(1,nlm);
GBTcon.dist=zeros(1,ndm);
GBTcon.glob=zeros(1,ngm); %can be less than 4 for special sections
GBTcon.other=zeros(1,nom);
%to turn on any cFSM classification use 1's on the GBTcon vectors defined
%above. Leave off unless you want it.
%for example say you wanted distortional only->GBTcon.dist=ones(1,ndm);


%-----------------------------------------------------------------
%---------------RUN THE ANALYSIS----------------------------------
%-----------------------------------------------------------------
%
[curve,shapes]=strip(prop,node,elem,lengths,springs,constraints,GBTcon,BC,m_all,neigs);
%

%-----------------------------------------------------------------
%---------------Modal identification------------------------------
%-----------------------------------------------------------------
% although cFSM feature is not explained in this batch code, the interface
% of how to call them is provided. I provide clas (set as 0) here for the
% later use in plot.
%
% GBTcon.norm=1;
% %norm:
% %   method=1=vector norm
% %   method=2=strain energy norm
% %   method=3=work norm
% clas=classify(prop,node,elem,lengths,shapes,GBTcon,BC,m_all);
clas=0;


%for example you could save the filename and look at the results in the GUI
fname=['batchcufsm4results' num2str(number)];
save (fname);

%Hopefully it is obvious that all of the above could be put in a loop. One
%of the dimensions could be a variable and then all of the preceding steps
%could go inside a looop - in the preceding only the use of CUFSM in batch
%mode is demonstrated.

%
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%
%
%-----------------------------------------------------------------
%-------------POST PROCESSING EXAMPLES----------------------------
%-----------------------------------------------------------------
%Traditional post-processing can be done within the GUI, but perhaps in
%some circumstances it is desirable to do some post-processing here. If you
%were really doing a parameter study you would likely want to make your own
%plots and not used canned stuff from CUFSM but here I at least show how
%you would used the canned CUFSM stuff in case you would want to.
%note, the major thing is you need to understand how the results are
%stored, particularly curve and shapes.
%curve: buckling curve (load factor) for each length
%curve{l} = [ length mode#1
%             length mode#2
%             ...    ...
%             length mode#]
%for curve, inside each cell that corresponds a length, the length in matrix is the
%same
%shapes = mode shapes for each length
%shapes{l} = mode, mode is a matrix, each column corresponds to a mode.

%Figure 1 lets plot the usual buckling curve
% figure
% axes(axes_handle)
% 
% cla
% axescurve=axes('Units','normalized','Position',[0.1 0.1 0.8 0.8],'Box','on','XTickLabel','','YTickLabel','');
axescurve=axes_handle;
set(axescurve,'FontName','times','FontSize',14)
% figure
% axes2dshape=axes('Units','normalized','Position',[0.1 0.1 0.8 0.8],'visible','off');


%Set default initial values
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
%whether the solution is a signature curve solution or general boundary
%condition solution
for j=1:max(size(m_all))
    if length(m_all{j})==1&m_all{j}==m_all{1}
        solutiontype=1;
    else
        solutiontype=2;
        break
    end
end
solutiontypecell{i}=solutiontype;

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

%Plot the buckling curve as a start
curveoption=1; %plot the buckling curve vs length
if length(lengths)==1
    curveoption=2;%plot the buckling curve vs mode#
end
%set the xmax, ymax
if curveoption==1
    xmin=min(lengths)*10/11;
    ymin=0;
    xmax=max(lengths)*11/10;
    for j=1:max(size(curve));
        curve_sign(j,1)=curve{j}(modeindex,1);
        curve_sign(j,2)=curve{j}(modeindex,2);
    end
    ymax=min([max(curve_sign(:,2)),3*median(curve_sign(:,2))]);
   
    %plot the buckling curve vs length
    picpoint=[curve{lengthindex}(modeindex,1) curve{lengthindex}(modeindex,2)];
%    
    thecurve3(curvecell,filenamecell,clascell,filedisplay,minopt,logopt,clasopt, axescurve,xmin,xmax,ymin,ymax,modedisplay,fileindex,modeindex,picpoint)
    
elseif curveoption==2
    
    figure(11)
    axescurvemode=axes('Units','normalized','Position',[0.1 0.1 0.8 0.8],'Box','on','XTickLabel','','YTickLabel','');
    set(axescurvemode,'FontName','times','FontSize',14)
    xmin=1;
    ymin=0;
    xmax=length(curve{lengthindex}(:,2));
    ymax=min([max(curve{lengthindex}(:,2)),3*median(curve{lengthindex}(:,2))]);
    %plot the buckling load vs mode#
    picpoint=[modes(modeindex) curve{lengthindex}(modeindex,2)];
    thecurve3mode(curvecell,filenamecell,clascell,fileindex,minopt,logopt,clasopt,axescurvemode,xmin,xmax,ymin,ymax,fileindex,lengthindex,picpoint);
   
    %plot longitudinal term participation for modeindex
    figure(12)
    set(gca,'FontName','times','FontSize',14)
    
    mode=shapes{lengthindex}(:,modeindex);
    nnodes=length(node(:,1));
    m_a=m_all{lengthindex};
    [d_part]=longtermpart(nnodes,mode,m_a);
    bar(m_a,d_part);hold on
    xlabel('m, longitudinal term');hold on;
    ylabel('Participation');hold on;
    legendstring=[filenamecell{fileindex},', length = ',num2str(lengths(lengthindex)), ', mode = ', num2str(modeindex)];
        hlegend=legend(legendstring);hold on
    set(hlegend,'Location','best');
    hold off  
end
hold on
ylim=get(gca,'Ylim');
plot([large_piece large_piece],ylim,'color','r');

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
dispshap(undef,node,elem,mode,axes2dshape,scale,springs,m_a,BC,SurfPos);
if threed==1
%     figure(3)
%     axes3dshape=axes('Units','normalized','Position',[0.1 0.1 0.8 0.8],'visible','off');
  cla axes3dshape
dispshp2(undef,lengths(lengthindex),node,elem,mode,axes3dshape,scale,m_a,BC,ifpatch);
end


end

