function []=cufsm_to_abaqus(filename,savename,L,nele,mag,analysis,elemtype,nelbsp,nelbct)
%
%cufsm_to_abaqus.m
%Ben Schafer
%December 2005
%
%Objective
%Creation of an ABAQUS model from a CUFSM model
%
%Assumptions
%(1)ABAQUS S9R5 element is used, which means 1 FE element for every two from
%CUFSM, since S9R5 has mid side nodes (it is a parabolic element). It is
%assumed that the user has set up the CUFSM input file so that this will
%work properly, i.e., an even number of elements fits in any location, etc.
%the code proceeds through the model by FSM element pairs, so it can work for
%multi-branched and disconnected models if the user is careful
%(2)Linear elastic material assumed, Ex,vx are read from CUFSM prop
%(3)Applied stresses vary linearly
%(4)Should use an even number of FE elements along the length so that the
%midle section actually falls in the midle.
%
%CUFSM INPUT
%CUFSM input file
%%%filename=['templateC'];
load(filename);
%from a CUFSM file with a completed analysis in it, we should be given
% CUFSM INPUTS
% prop: [matnum Ex Ey vx vy G] 6 x nmats
% node: [node# x z dofx dofz dofy dofrot stress] nnodes x 8;
% elem: [elem# nodei nodej t matnum] nelems x 5;
% lengths: lengths to be analyzed
% springs: [node# d.o.f. kspring kflag] where 1=x dir 2= z dir 3 = y dir 4 = q dir (twist) flag says if k is a foundation stiffness or a total stiffness
% constraints: [node#e dofe coeff node#k dofk] e=dof to be eliminated k=kept dof dofe_node = coeff*dofk_nodek 
% GBTcon: GBTcon.glob,GBTcon.dist, GBTcon.local, GBTcon.other vectors of 1's
%   and 0's referring to inclusion (1) or exclusion (0) of a given mode 
% CUFSM OUTPUTS
% curve: 3D array: [lengths LF] x mode#
% shapes: 3D array: [nodal_displacement lengths] x mode#
% clas: 3D array [L D G O lengths] x mode# L,D,G,O are est. participations
%
%ADDITIONAL INPUT
%L=physical member length in ABAQUS model
%%%L=100;
%nele=number of S9R5 elements along the length
%%%nele=40;
%mag=magnitude of imperfection, based off CUFSM modes, size=lengths
%%%mag=zeros(length(lengths),1);
%analysis=1 for buckling, 2 for static analysis
%%%analysis=1;
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Added for selection of element types besides of S9R5, but also S4 and S4R.
%Moreover, it can also take account for springs and constraints by
%defining springs and constraints along the longitudinal length. 
%(Zhanjie 2008)
%
%%%elemtype: a string of 'S9R5', 'S4', 'S4R' which is the element types in
%%%ABAQUS. (Zhanjie 2008)
%nelbsp: Number of elements between springs, used to discretize springs
%nelbct: Number of elements between constraints, used to discretize const.
%(Luiz 2008)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%WARNINGS
%Has to be an even number of FSM elements for this to work
if strcmp(elemtype,'S9R5')
    if rem(length(elem(:,1)),2)>0
        ['Warning! Your CUFSM model has an odd number of elements this will not convert to ABAQUS S9R5 elements. Please modify your model so that the number of elements is an even number']
    end
end
   
%

%PRELIMINARIES
%Count FSM nodes and modes
nnodes=length(node(:,1)); %Number of FSM cross-section nodes
if iscell(curve)
    nmodes=max(size(curve)); %Number of FSM mode shapes for first mode, same as number of lengths
else
    nmodes=length(curve(:,1)); %Number of FSM mode shapes for first mode, same as number of lengths
end
%Determine FE number of nodes and increment (Zhanjie 2007)
if strcmp(elemtype,'S9R5')
    nL=2*nele+1; %Number of FE nodes along the length
elseif strcmp(elemtype,'S4') || strcmp(elemtype, 'S4R')
    nL=nele+1;   %Number of FE nodes along the length
end
%Determine the node numbering increment along the length
if nnodes<100
    FEsection_increment=100; %so along the length the numbering goes up by 100's
else
    FEsection_increment=nnodes+1;
end
%
%
%NODAL COORDINATES
%Using mag, build up the imperfect nodal coordiantes
%Begin with the undeformed geometry
undefx=zeros(nnodes,nL);
undefz=zeros(nnodes,nL);
for i=1:nL
    undefx(:,i)=node(:,2);
    undefz(:,i)=node(:,3);
    undefy(:,i)=ones(nnodes,1)*(i-1)/(nL-1)*L;
end
%Define variable for the deformed/imperfect shape
shapex=undefx;
shapez=undefz;
shapey=undefy;
if analysis==2
    %Loop through all the possible modes to build the shape
    multa=find(mag~=0);
    for j=1:length(multa)
        if iscell(shapes)
            shape=shapes{multa(j)}(:,1);
            m_a=m_all{multa(j)};
        else
            shape=shapes(:,multa(j));
            m_a=[1];
        end
        if iscell(curve)
            wvl=curve{multa(j)}(1,1);
        else
            wvl=curve(multa(j),1);
        end
        %         shape=shapes(:,j);
        %         wvl=curve(j,1);
        mult=mag(multa(j)); %this is the variable which determines magnitude
        defx=zeros(nnodes,nL);
        defy=zeros(nnodes,nL);
        defz=zeros(nnodes,nL);
        
        %get the basic x,y,z displacements of the mode
        if length(m_a)>1
            %general boundary condition results
            totalm=length(m_a);
            for mn=1:totalm
                longiterm=m_a(mn);
                km=longiterm*pi/L;
                for i=1:nnodes
                    modexm(i,1)=shape(4*nnodes*(mn-1)+2*i-1,1);
                    modezm(i,1)=shape(4*nnodes*(mn-1)+2*nnodes+2*i-1,1);
                    modeym(i,1)=shape(4*nnodes*(mn-1)+2*i,1);
                end
                for dd=1:nL
                    yloc=(dd-1)/(nL-1)*L;
                    if strcmp(BC,'S-S')   %consider longit shape function
                        s=sin(longiterm*pi*yloc/L);
                        c=cos(longiterm*pi*yloc/L)*longiterm*pi/L;
                    elseif strcmp(BC,'C-C')   %consider longit shape function
                        s=sin(longiterm*pi*yloc/L)*sin(pi*yloc/L);
                        c=cos(longiterm*pi*yloc/L)*longiterm*pi/L*sin(pi*yloc/L)+sin(longiterm*pi*yloc/L)*cos(pi*yloc/L)*pi/L;                        
                    elseif strcmp(BC,'S-C')   %consider longit shape function
                        s=sin((longiterm+1)*pi*yloc/L)+(longiterm+1)*sin(longiterm*pi*yloc/L)/longiterm;
                        c=cos((longiterm+1)*pi*yloc/L)*(longiterm+1)*pi/L+(longiterm+1)*cos(longiterm*pi*yloc/L)*pi/L;
                    elseif strcmp(BC,'C-F')   %consider longit shape function
                        s=1-cos((longiterm-1/2)*pi*yloc/L);
                        c=sin((longiterm-1/2)*pi*yloc/L)*(longiterm-1/2)*pi/L;
                    elseif strcmp(BC,'C-G')   %consider longit shape function
                        s=sin((longiterm-1/2)*pi*yloc/L)*sin(pi*yloc/2/L);
                        c=cos((longiterm-1/2)*pi*yloc/L)*(longiterm-1/2)*pi/L*sin(pi*yloc/L/2)+sin((longiterm-1/2)*pi*yloc/L)*cos(pi*yloc/2/L)*pi/2/L;
                    end
                    
                    defx(:,dd)=modexm(:,1)*s+defx(:,dd);
                    defz(:,dd)=modezm(:,1)*s+defz(:,dd);
                    defy(:,dd)=modeym(:,1)*c/km+defy(:,dd);
                end
            end
            maxis=max([max(max(abs(defx))) max(max(abs(defz))) max(max(abs(defy)))]);
            modex=defx/maxis*mult;
            modey=defy/maxis*mult;
            modez=defz/maxis*mult;
            %add perturbed shape to shape
            shapex=shapex+modex;
            shapez=shapez+modez;
            shapey=shapey+modey;
        else
            % signature curve results
            for i=1:nnodes
                modex(i,1)=shape(2*i-1,1);
                modez(i,1)=shape(2*nnodes+2*i-1,1);
                modey(i,1)=shape(2*i,1);
            end
            %normalize the basic displacements to 1
            maxis=max([max(abs(modex)) max(abs(modez)) max(abs(modey))]);
            modex=modex/maxis;
            modey=modey/maxis;
            modez=modez/maxis;
            %multiply the displacement to the appropriate maximum magnifier
            modex=modex*mult;
            modey=modey*mult;
            modez=modez*mult;
            %generate along the length
            for i=1:nL
                yloc=(i-1)/(nL-1)*L;
                defx(:,i)=modex(:,1)*sin(pi*yloc/wvl);
                defz(:,i)=modez(:,1)*sin(pi*yloc/wvl);
                defy(:,i)=modey(:,1)*sin(pi*yloc/wvl);
            end
            %add perturbed shape to shape
            shapex=shapex+defx;
            shapez=shapez+defz;
            shapey=shapey+defy;
        end
    end
end
%You can make a quick dirty plot of the nodal coordinates in the perturbed
%geometry if you want to see them. This was used in debugging.
%plot3(shapex,shapey,-shapez)
%
%NODAL COORDINATES IN FE FORM
%FEnode=[node# x y z] x is along the length
for i=1:nL
	FEnode((i-1)*nnodes+1:(i-1)*nnodes+nnodes,1)=node(:,1)+FEsection_increment*(i-1);
	%x in finite elements coordinate system
	FEnode((i-1)*nnodes+1:(i-1)*nnodes+nnodes,2)=shapey(:,i);
	%y in finite elements coordinate system
	FEnode((i-1)*nnodes+1:(i-1)*nnodes+nnodes,3)=shapex(:,i);
	%z in finite elements coordinate system
	FEnode((i-1)*nnodes+1:(i-1)*nnodes+nnodes,4)=shapez(:,i);
end
%
if strcmp(elemtype,'S9R5')
    %ELEMENT DEFINITIONS
    %Note, elements are done in strips, to reflect thickness or material
    %changes which could exist in the FSM model. BUT! the S9R5 element used
    %grabs the elements in pairs, so thickness or material should be the same
    %for pairs of elements.
    %Note
    %Element Connectivity
    %ABAQUS S9R5 for example
    %       face3
    %      4--7--3
    %      |     |
    %face4 8  9  6 face2
    %      |     |
    %      1--5--2
    %       face1
    %Go through the elements 2 at a time
    k=1; %counter for strip, i.e., strip along the length (group of 2 FSM strips)
    elemnum=1;
    for i=1:2:length(elem(:,1))-1
        for j=1:(nL-1)/2
            n1=elem(i,2);
            n2=elem(i+1,3);
            n3=elem(i+1,3) + 2*FEsection_increment;
            n4=elem(i,2)   + 2*FEsection_increment;
            n5=elem(i,3);
            n6=elem(i+1,3) + FEsection_increment;
            n7=elem(i,3)   + 2*FEsection_increment;
            n8=elem(i,2)   + FEsection_increment;
            n9=elem(i,3) + FEsection_increment;
            nnum=[n1 n2 n3 n4 n5 n6 n7 n8 n9]+(j-1)*2*FEsection_increment;
            enum=elemnum;%nnum(5);
            FEelem(j,:,k)=[enum nnum];
            elemnum=elemnum+1;            
        end
        k=k+1;
    end
    %ELEMENT THICKNESS AND MATERIAL NUMBER
    k=1;
    for i=1:2:length(elem(:,1))-1
        t(k)=elem(i,4);
        matnum(k)=elem(i,5);
        k=k+1;
    end
elseif strcmp(elemtype,'S4') || strcmp(elemtype, 'S4R')
    %For S4 or S4R   (Zhanjie 2007)
    %Similar procedure as explained for S9R5
    %Element Connectivity
    %ABAQUS S4 or S4R for example
    %       face3
    %      4----3
    %face4 |    |face2
    %      1----2
    %      face1
    %
    k=1; %counter for strip, i.e., strip along the length
    elemnum=1;
    for i=1:1:length(elem(:,1))
        for j=1:(nL-1)
            n1=elem(i,2);
            n2=elem(i,3);
            n3=elem(i,3) + FEsection_increment;
            n4=elem(i,2) + FEsection_increment;
            nnum=[n1 n2 n3 n4]+(j-1)*FEsection_increment;
            enum=elemnum;%nnum(1)
            FEelem(j,:,k)=[enum nnum];
            elemnum=elemnum+1;
        end
        k=k+1;
    end
    %ELEMENT THICKNESS AND MATERIAL NUMBER
    k=1;
    for i=1:1:length(elem(:,1))
        t(k)=elem(i,4);
        matnum(k)=elem(i,5);
        k=k+1;
    end
end



%MATERIAL DEFINITION (ELASTIC ONLY)
for i=1:length(prop(:,1))
    mat(i,:)=[prop(i,1) prop(i,2) prop(i,4)];
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add Springs (Luiz 2007)
firstelem=max(max(max(FEelem)))+1;
if springs==0
else
    %Create the element spring and specific the node
    %%%%Warning:User has to change nelbsp, default=1%%%%
    %nelbsp=12; %Number of elements between springs, used to discretize springs,
    %next step it should attached to the GUI.
    springelem2=0;
    springnode2=0;
    for j=1:length(springs(:,1)) %number of springs defined on CUFSM
        for i=1:((nele/nelbsp)+1) %(nele+1) is the number of springs along the length.  %Zhanjie ?
            springelem(i,1)=firstelem+(i-1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Modified to be applicable for S4 and S4R  (Zhanjie 2008)
            if strcmp(elemtype,'S9R5')
                springnode(i,1)=springs(j,1)+2*nelbsp*FEsection_increment*(i-1);
            elseif strcmp(elemtype,'S4') || strcmp(elemtype, 'S4R')
                springnode(i,1)=springs(j,1)+nelbsp*FEsection_increment*(i-1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %x in finite elements coordinate system, select the nodes to apply the
            %springs
        end
        %
        %Match the direction of springs between CUFSM and Abaqus
        if springs(j,2)==1
            springs(j,2)=2;
        elseif springs(j,2)==2
            springs(j,2)=3;
        elseif springs(j,2)==3
            springs(j,2)=1;
        else
            springs(j,2)=4;
        end
        %
        %Spring stiffness: 0: total stiffness and 1: foundation stiffness
        if springs(j,4)==0
            springs(j,3)=springs(j,3)/((nele/nelbsp)+1);
        else
            springs(j,3)=springs(j,3)*(L/((nele/nelbsp)+1));
        end
        %Matrix springelem
        firstelem=max(springelem(:,1))+1;
        if springelem2==0
            springelem2=springelem;
        else
            springelem2=[springelem2' springelem'];
            springelem2=springelem2';
        end
        %Matrix springnode
        if springnode2==0
            springnode2=springnode;
        else
            springnode2=[springnode2' springnode'];
            springnode2=springnode2';
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add General Constraints (Luiz 2007), improved 2008 to add spring between
%nodes
%%%%Warning:User has to change nelbct (default=1) and kspring (default=0))%%%%
%%%%%%%%%%%%%%Input%%%%%%%%%%%%%
% nelbct=1; %Number of elements between contraints, used to discretize const.
%Add General Constraints to consider springs between nodes (Luiz 2008)
kspring1=0; %direction 3 in CUFSM or 1 in Abaqus
kspring2=0; %direction 1 in CUFSM or 2 in Abaqus
kspring3=0; %direction 2 in CUFSM or 3 in Abaqus
kspring4=0; %direction 4 in CUFSM or 4 in Abaqus
%Spring stiffness, if equal zero, general constraint are 
%considered, if not spring between nodes are considered
%%%%%%%%%%%%End Input%%%%%%%%%%%
gct=0; %Final general constraint matrix
if constraints==0
else
    for i=1:length(constraints(:,1))
        for j=2:3:5
            %Match the direction of springs between CUFSM and Abaqus
            if constraints(i,j)==1
                constraints(i,j)=2;
            elseif constraints(i,j)==2
                constraints(i,j)=3;
            elseif constraints(i,j)==3
                constraints(i,j)=1;
            else
                constraints(i,j)=4;
            end
        end
        % Create all the nodes to be restrained and relations
        for j=1:((nele/nelbct)+1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Modified to be applicable for S4 and S4R  (Zhanjie 2008)
            if strcmp(elemtype,'S9R5')
                gc(j,1)=constraints(i,1)+2*nelbct*FEsection_increment*(j-1);
                gc(j,2)=constraints(i,2);
                gc(j,3)=constraints(i,3);
                gc(j,4)=constraints(i,4)+2*nelbct*FEsection_increment*(j-1);
                gc(j,5)=constraints(i,5);
            elseif strcmp(elemtype,'S4') || strcmp(elemtype, 'S4R')
                gc(j,1)=constraints(i,1)+nelbct*FEsection_increment*(j-1);
                gc(j,2)=constraints(i,2);
                gc(j,3)=constraints(i,3);
                gc(j,4)=constraints(i,4)+nelbct*FEsection_increment*(j-1);
                gc(j,5)=constraints(i,5);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if gct==0
            gct=gc;
        else
            gct=[gct' gc'];
            gct=gct';
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%BOUNDARY CONDITIONS
%The default boundary conditions are set to mimic FSM shape functions to
%the extent possible. So, pinned ends, but longitudinal pinned at the mid.
%Nodes on the boundaries and the middle
end1=[1:1:nnodes];
end2=end1 + (nL-1)*FEsection_increment;
middle=end1 + (nL-1)/2*FEsection_increment;
%
%LOAD
%Generate consistent nodal loads from the CUFSM stresses
if strcmp(elemtype,'S9R5')
    %Consistent nodal loads for a parabolic element with linear variation in
    %stress across the element for S9R5
    %f1-----.n1        1/6*f1*h*t-->.n1
    %   ----|                       |
    %    ---.n2   1/3*(f1+f2)*h*t-->.n2
    %     --|                       |
    %    f2-.n3        1/6*f2*h*t-->.n3
    %set up cload vector
    for i=1:nnodes
        cload(i,:)=[i 0];
    end
    for i=1:2:length(elem(:,1))-1
        n1=elem(i,2);
        n2=elem(i,3);
        n3=elem(i+1,3);
        f1=node(n1,8);
        f2=node(n3,8);
        xn1=node(n1,2);
        yn1=node(n1,3);
        xn3=node(n3,2);
        yn3=node(n3,3);
        h=sqrt((xn3-xn1)^2+(yn3-yn1)^2);
        ti=elem(i,4);
        Fn1=1/6*f1*h*ti;
        Fn2=1/3*(f1+f2)*h*ti;
        Fn3=1/6*f2*h*ti;
        cload(n1,2)=cload(n1,2)+Fn1;
        cload(n2,2)=cload(n2,2)+Fn2;
        cload(n3,2)=cload(n3,2)+Fn3;
    end
    end1cload=cload;
    end2cload=[end2' -cload(:,2)];
elseif strcmp(elemtype,'S4') || strcmp(elemtype, 'S4R')
    %for S4 or S4R
    %Consistent nodal loads for an element with linear variation in
    %stress across the element based on the shape function
    %f1-----.n1        1/3*f1*h*t+1/6*f2*h*t-->.n1
    %   ----|                                   |
    %     --|                                   |
    %    f2-.n2        1/3*f2*h*t+1/6*f1*h*t-->.n2
    %set up cload vector
    for i=1:nnodes
        cload(i,:)=[i 0];
    end
    for i=1:1:length(elem(:,1))
        n1=elem(i,2);
        n2=elem(i,3);
        f1=node(n1,8);
        f2=node(n2,8);
        xn1=node(n1,2);
        yn1=node(n1,3);
        xn2=node(n2,2);
        yn2=node(n2,3);
        h=sqrt((xn2-xn1)^2+(yn2-yn1)^2);
        ti=elem(i,4);
        Fn1=1/3*f1*h*ti+1/6*f2*h*ti;
        Fn2=1/3*f2*h*ti+1/6*f1*h*ti;
        cload(n1,2)=cload(n1,2)+Fn1;
        cload(n2,2)=cload(n2,2)+Fn2;
    end
    end1cload=cload;
    end2cload=[end2' -cload(:,2)];
end
%
%
%
%Generate ABAQUS file
fid=fopen([savename,'.inp'],'w');
fprintf(fid,'*HEADING\n');
fprintf(fid,'%s\n',[savename,' ',datestr(now)]);
fprintf(fid,'*NODE,NSET=ALL\n');
fprintf(fid,'%d,%f,%f,%f \n',FEnode');
for k=1:size(FEelem,3)
    fprintf(fid,'*ELEMENT,TYPE=%s,ELSET=STRIP%d\n',elemtype,k);
    if strcmp(elemtype,'S9R5')
        fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d \n',FEelem(:,:,k)');
    elseif strcmp(elemtype,'S4') || strcmp(elemtype, 'S4R') %added for S4 and S4R (Zhanjie 2007)
        fprintf(fid,'%d,%d,%d,%d,%d \n',FEelem(:,:,k)');
    end
    fprintf(fid,'*SHELL SECTION, ELSET=STRIP%d, MATERIAL=MAT%d\n',[k matnum(k)]);
    fprintf(fid,'%f \n',t(k));
end
for i=1:size(mat,1)
    fprintf(fid,'*MATERIAL, NAME=MAT%d\n',mat(i,1));
    fprintf(fid,'*ELASTIC, TYPE=ISOTROPIC\n');
    fprintf(fid,'%f,%f \n',[mat(i,2:3)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spring (Luiz 2007)
if springs==0
else
fprintf(fid,'*Element, type=spring1\n');
for i=1:length(springelem2)
    fprintf(fid,'%g, %g\n',springelem2(i),springnode2(i));
end
%
for i=1:length(springs(:,1))
fprintf(fid,'*Elset, elset=spring%d, generate\n',i);
fprintf(fid,'%g,%g,1\n',springelem2(1+(length(springelem)*(i-1))),springelem2((length(springelem)*(i-1))+length(springelem)));
end
%
for i=1:length(springs(:,1))
fprintf(fid,'*Spring, elset=spring%d\n',i);
fprintf(fid,'%g, %g\n',springs(i,2),springs(i,2));
fprintf(fid,'%g\n',springs(i,3));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Constraints (Luiz 2007) improved (2008) to considere spring
%between nodes
if constraints==0
else
    if kspring1==0 & kspring2==0 & kspring3==0 & kspring4==0
        fprintf(fid,'*EQUATION\n');
        for i=1:length(gct(:,1))
        fprintf(fid,'%g\n',2);
        fprintf(fid,'%g, %g, %g, %g, %g, %g\n',gct(i,1),gct(i,2),1,gct(i,4),gct(i,5),-gct(i,3));
        end
    else
        fprintf(fid,'*Element, type=spring2\n');
        for i=1:length(gct(:,1))
        fprintf(fid,'%g, %g, %g\n',firstelem+i,gct(i,1),gct(i,4));
        end
%   
        for j=1:length(constraints(:,1))
        fprintf(fid,'*Elset, elset=dispspring%d, generate\n',j);
        fprintf(fid,'%g,%g,1\n',firstelem+1+((j-1)*((nele/nelbct)+1)),firstelem+1+(nele/nelbct)+((j-1)*((nele/nelbct)+1)));
        end
%
        for j=1:length(constraints(:,1))
        fprintf(fid,'*Spring, elset=dispspring%d\n',j);
        fprintf(fid,'%g, %g\n',gct(1+((j-1)*((nele/nelbct)+1)),2),gct(1+((j-1)*((nele/nelbct)+1)),2));
            if gct(1+((j-1)*((nele/nelbct)+1)),2)==1
                fprintf(fid,'%g,\n',kspring1);
            elseif gct(1+((j-1)*((nele/nelbct)+1)),2)==2
                fprintf(fid,'%g,\n',kspring2);             
            elseif gct(1+((j-1)*((nele/nelbct)+1)),2)==3
                fprintf(fid,'%g,\n',kspring3);
            else
                fprintf(fid,'%g,\n',kspring4);
            end
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'*NSET, NSET=END1\n');
fprintf(fid,'%d\n',end1');
fprintf(fid,'*NSET, NSET=END2\n');
fprintf(fid,'%d\n',end2');
fprintf(fid,'*NSET, NSET=MIDDLE\n');
fprintf(fid,'%d\n',middle');
fprintf(fid,'*BOUNDARY\n');
fprintf(fid,'END1,2,3\n');
fprintf(fid,'END2,2,3\n');
fprintf(fid,'MIDDLE,1,1\n');
fprintf(fid,'*STEP\n');
if analysis==1
    fprintf(fid,'*BUCKLE\n');
    fprintf(fid,'10,9,100\n');
elseif analysis==2
    fprintf(fid,'*STATIC\n');
    fprintf(fid,'0.1\n');  
end
fprintf(fid,'*CLOAD\n');
fprintf(fid,'%d,1,%f\n',end1cload');
fprintf(fid,'%d,1,%f\n',end2cload');
fprintf(fid,'*OUTPUT, FIELD\n');
fprintf(fid,'*ELEMENT OUTPUT\n');
fprintf(fid,'1,3,5\n');
fprintf(fid,'S,E\n');
fprintf(fid,'*NODE OUTPUT\n');
fprintf(fid,'U\n');
fprintf(fid,'*END STEP\n');
status=fclose(fid);
