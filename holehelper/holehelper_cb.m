function []=holehelper_cb(num);
%bws
%March 2016
%J. Cai 2016
%
%general
global fig screen prop node elem lengths curve shapes clas springs constraints GBTcon BC m_all neigs version screen
%output from pre2
global subfig ed_prop ed_node ed_elem ed_lengths axestop screen flags modeflag ed_springs ed_constraints
%output from load file
global filename
%output from holehelper
global popanel flags axessect input1edit input2edit input3edit label_cs node_all cross_section_range var_length node_all4plot axeslong
global b bl1 bl2 b13  node4local elem4local A 
global curve_local shapes_local length_local axescurve_hole axes2dshape_hole
global halfwave hint waveedit label_title len_cur_local uplengthl downlengthl

global undefl nodel eleml model axes2dshapel scalel springsl m_al BCl SurfPosl shapesl len_cur_distorl idx_hole_l lengths_local  
global curvecelll filenamecelll clascelll filedisplayl minoptl logoptl clasoptl axescurvel xminl xmaxl yminl ymaxl modedisplayl  fileindexl modeindexl picpointl
global screen Ix_avg Iy_avg I1_avg I2_avg len_cur_local modeindex lengthindex
global xmin xmax ymin ymax modeindex_local picpoint curvel shapescelll J_avg Cw_avg record_hole table_hole table popanel_hole s_coordinate large_piece flags_hole
global xcg zcg thetap I11 I22 w Cw B S2edit input_Lcrd_edit elem_distor distor_len lgd_flag curve_distor shapes_distor distor_len_sig elem_all upperbox
global label_title_load_ratio
global   out2A out2J out2xcg out2zcg out2Ixx out2Izz out2Ixz out2thetap out2I11 out2I22 out2C_w netProp label_title lowerbox currentlocation
global KL1edit KL2edit KL3edit
global Pedit_hole M11edit_hole M22edit_hole
global rad_Pe rad_Me1 rad_Me2 rad_Me12 text_Pe text_Me1 text_Me2 statictext_Pe statictext_Me1 statictext_Me2 statictext_KL1 statictext_KL2 statictext_KL3 statictext_exy1
global statictext_exy2 ed_KL1 ed_KL2 ed_KL3 ed_exy ed_c text_maxmode slider rad_Me text_mode props4global rad_shear
global rad_origin rad_centroid rad_axisxy rad_axis12 rad_axial rad_deform rad_node t_axes CS_prop ed_c_hole range_slider bl3

%     if get(rad_Pe,'Value');
%         force = 'Pe';
%         set(rad_Me1,'Visible','off');
%         set(rad_Me2,'Visible','off');
%         set(rad_Me12,'Visible','off');   
%         set(text_Pe,'Visible','on');
%         set(text_Me1,'Visible','off');
%         set(text_Me2,'Visible','off');
%         set(statictext_Pe,'Visible','on','string','Pe:');
%         set(statictext_Me1,'Visible','off');
%         set(statictext_Me2,'Visible','off');                        
%         set(statictext_KL1,'Visible','on');
%         set(statictext_KL2,'Visible','on');
%         set(statictext_exy1,'Visible','on');
%         set(statictext_exy2,'Visible','on');
%         set(ed_KL1,'Visible','on');
%         set(ed_KL2,'Visible','on');
%         set(ed_exy,'Visible','on');
%         set(ed_c,'Visible','off');
%         set(text_maxmode,'String','/3');
%         set(slider,'Max',3,'sliderstep',[0.5 0.5]);

switch num 

    case 1
    %---------------------
    %plot the cross-section properly
     crossNo=str2num(get(label_cs,'string'));
     if crossNo==1
          crossect_hole(node,elem,axessect,springs,constraints,flags_hole)
     else
     
%    crossect_hole(node4local{crossNo},elem4local{crossNo},axessect,springs,constraints,flags_hole)  
    crossect_hole(node4local{crossNo},elem4local{crossNo},axessect,springs,constraints,flags_hole)
     end
    %---------------------
    
    case 3 % help
        name=['Definiton of Variable for Hole Effect Analysis Tool'];
        fhh=figure('Name',name,...
            'NumberTitle','off');
        % set(fig,'Units','normalized')%
        set(fhh,'MenuBar','none');
        apic=imread('hole_help_pic.jpg','jpg');
        image(apic);
        axis('equal')
        axis('off')
    
    case 5
    %%
    %Previous Section Button    
    %-----------------------------
         crossNo=str2num(get(label_cs,'string')); 
        if crossNo~=1
            crossNo=crossNo-1;
             set(label_cs,'string',num2str(crossNo)); 
      Net_cross_section_plotter4500( axessect,axeslong,cross_section_range{crossNo},var_length,node_all4plot{crossNo});
      set(label_cs,'string',num2str(crossNo)); 
      cla(axessect)
crossect(node4local{crossNo},elem4local{crossNo},axessect,springs,constraints,flags_hole)
 isEmpytAxes = isempty(get(axes2dshape_hole, 'Children'));
 if isEmpytAxes~=1
     if lgd_flag==1
  [curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4( curve_local{crossNo},shapes_local{crossNo},prop,node4local{crossNo},elem4local{crossNo},crossNo,length_local{crossNo}, axescurve_hole,axes2dshape_hole);
     elseif lgd_flag==2
         [curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4(curve_distor{crossNo},shapes_distor{crossNo},prop,node,elem_distor{crossNo},crossNo,distor_len_sig, axescurve_hole,axes2dshape_hole);
    
%        for m=1:length(distor_len_sig)
% 
%             if distor_len_sig(m)==distor_len
%                 idx_hole_l=m;
%                 m
% break
%                 
%             end
%         end
%         
% set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
% picpointl=[curvel{idx_hole_l}(modeindexl,1) curvel{idx_hole_l}(modeindexl,2)];
%   thecurve3(curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl)
     end
     
%     set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
 label=[' load factor=',num2str(curvel{idx_hole_l}(modeindexl,2))];
set(label_title,'String',label);
set(label_title,'visible','on');
     
     
set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
 end
 
 [maxNode,~]=size(node4local{crossNo});
%calculate distance along the cross section to each node("s" coordinate)
node_location= node4local{crossNo}(:,1:3);
node_location(1,4)=0;
[maxNode,~]=size(node4local{crossNo});
ss_coordinate(1)=0;
for i=2:maxNode
    dx{i}=node_location(i,2)- node_location((i-1),2);
    dz{i}=node_location(i,3)- node_location((i-1),3);
    ds(i)=sqrt(dx{i}^2 + dz{i}^2);
    
    ss_coordinate(i,1)=ss_coordinate(i-1)+ds(i);
end 
s_display(:,1)=1:maxNode;
s_display(:,2)=ss_coordinate;
set(S2edit,'String',num2str(s_display));
if lgd_flag==1
set(hint,'String',['Hint: the longest piece of this CS is ' num2str(large_piece(crossNo))]);

end
 if lgd_flag==3
   set(CS_prop,...
        'String',['Cross-Section ' num2str(crossNo) ' Properties']);  
     
set(out2A,...
        'String',num2str(netProp{crossNo}.A));
  
    set(out2J,...
        'String',num2str(netProp{crossNo}.J));
    set(out2xcg,...
        'String',num2str(netProp{crossNo}.xcg));
    set(out2zcg,...
        'String',num2str(netProp{crossNo}.zcg));
    set(out2Ixx,...
        'String',num2str(netProp{crossNo}.Ixx));
    set(out2Izz,...
        'String',num2str(netProp{crossNo}.Izz));
    set(out2Ixz,...
        'String',num2str(netProp{crossNo}.Ixz));
    set(out2thetap,...
        'String',num2str(netProp{crossNo}.thetap));
    set(out2I11,...
        'String',num2str(netProp{crossNo}.I11));
    set(out2I22,...
        'String',num2str(netProp{crossNo}.I22));
    set(out2C_w,...
        'String',num2str(netProp{crossNo}.Cw));
          set(label_title_load_ratio,'Visible',' off');
      set(label_title,'Visible',' off');
 end

        end
    %-----------------------------

    case 6
    %%
    %Next Section Button
        %-----------------------------
        maxNo=length(node_all)
            crossNo=str2num(get(label_cs,'string')); 
        if crossNo~=maxNo
            crossNo=crossNo+1;
              set(label_cs,'string',num2str(crossNo)); 
      Net_cross_section_plotter4500( axessect,axeslong,cross_section_range{crossNo},var_length,node_all4plot{crossNo});
      set(label_cs,'string',num2str(crossNo)); 
      cla(axessect)
      cla(axescurve_hole)
crossect_hole(node4local{crossNo},elem4local{crossNo},axessect,springs,constraints,flags_hole)
 isEmpytAxes = isempty(get(axes2dshape_hole, 'Children'));
  if isEmpytAxes~=1
     if lgd_flag==1
  [curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4( curve_local{crossNo},shapes_local{crossNo},prop,node4local{crossNo},elem4local{crossNo},crossNo,length_local{crossNo}, axescurve_hole,axes2dshape_hole);
     elseif lgd_flag==2
         [curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4(curve_distor{crossNo},shapes_distor{crossNo},prop,node,elem_distor{crossNo},crossNo,distor_len_sig, axescurve_hole,axes2dshape_hole);
     elseif lgd_flag==3
            set(CS_prop,...
        'String',['Cross-Section ' num2str(crossNo) ' Properties']);  
        set(out2A,...
        'String',num2str(netProp{crossNo}.A));
  
    set(out2J,...
        'String',num2str(netProp{crossNo}.J));
    set(out2xcg,...
        'String',num2str(netProp{crossNo}.xcg));
    set(out2zcg,...
        'String',num2str(netProp{crossNo}.zcg));
    set(out2Ixx,...
        'String',num2str(netProp{crossNo}.Ixx));
    set(out2Izz,...
        'String',num2str(netProp{crossNo}.Izz));
    set(out2Ixz,...
        'String',num2str(netProp{crossNo}.Ixz));
    set(out2thetap,...
        'String',num2str(netProp{crossNo}.thetap));
    set(out2I11,...
        'String',num2str(netProp{crossNo}.I11));
    set(out2I22,...
        'String',num2str(netProp{crossNo}.I22));
    set(out2C_w,...
        'String',num2str(netProp{crossNo}.Cw));
%        for m=1:length(distor_len_sig)
% 
%             if distor_len_sig(m)==distor_len
%                 idx_hole_l=m;
%                 m
% break
%                 
%             end
%         end
%         
% set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
% picpointl=[curvel{idx_hole_l}(modeindexl,1) curvel{idx_hole_l}(modeindexl,2)];
%   thecurve3(curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl)
     end
%      set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
 label=[' load factor=',num2str(curvel{idx_hole_l}(modeindexl,2))];
set(label_title,'String',label);
set(label_title,'visible','on');
set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
 end
 [maxNode,~]=size(node4local{crossNo});
%calculate distance along the cross section to each node("s" coordinate)
node_location= node4local{crossNo}(:,1:3);
node_location(1,4)=0;
[maxNode,~]=size(node4local{crossNo});
ss_coordinate(1)=0;
for i=2:maxNode
    dx{i}=node_location(i,2)- node_location((i-1),2);
    dz{i}=node_location(i,3)- node_location((i-1),3);
    ds(i)=sqrt(dx{i}^2 + dz{i}^2);
    
    ss_coordinate(i,1)=ss_coordinate(i-1)+ds(i);
end 
s_display(:,1)=1:maxNode;
s_display(:,2)=ss_coordinate;
set(S2edit,'String',num2str(s_display));
if lgd_flag==1
set(hint,'String',['Hint: the longest piece of this CS is ' num2str(large_piece(crossNo))]);
end
if lgd_flag==3
    
                set(CS_prop,...
        'String',['Cross-Section ' num2str(crossNo) ' Properties']);  
set(out2A,...
        'String',num2str(netProp{crossNo}.A));
  
    set(out2J,...
        'String',num2str(netProp{crossNo}.J));
    set(out2xcg,...
        'String',num2str(netProp{crossNo}.xcg));
    set(out2zcg,...
        'String',num2str(netProp{crossNo}.zcg));
    set(out2Ixx,...
        'String',num2str(netProp{crossNo}.Ixx));
    set(out2Izz,...
        'String',num2str(netProp{crossNo}.Izz));
    set(out2Ixz,...
        'String',num2str(netProp{crossNo}.Ixz));
    set(out2thetap,...
        'String',num2str(netProp{crossNo}.thetap));
    set(out2I11,...
        'String',num2str(netProp{crossNo}.I11));
    set(out2I22,...
        'String',num2str(netProp{crossNo}.I22));
    %
    set(out2C_w,...
        'String',num2str(netProp{crossNo}.Cw));
          set(label_title_load_ratio,'Visible',' off');
      set(label_title,'Visible',' off');
 end
        end 
    %-----------------------------

    case 7
    %-----------------------------
    %Coarse button for gridded surface discretization
    Dtheta=str2num(get(Dthetaedit,'String'));
    Dphi=str2num(get(Dphiedit,'String'));
    Dtheta=ceil(Dtheta*2);
    Dphi=ceil(Dphi*2);
    set(Dthetaedit,'String',num2str(Dtheta));
    set(Dphiedit,'String',num2str(Dphi));    
    %-----------------------------

    case 8
    %-----------------------------
    %Fine button for gridded surface discretization
    Dtheta=str2num(get(Dthetaedit,'String'));
    Dphi=str2num(get(Dphiedit,'String'));
    Dtheta=ceil(Dtheta/2);
    Dphi=ceil(Dphi/2);
    set(Dthetaedit,'String',num2str(Dtheta));
    set(Dphiedit,'String',num2str(Dphi));    
    %-----------------------------

    
    case 10
    %-----------------------------------------------
    %Generate cross-sections with holes..
    %Remove RHS analysis from view, since new section generated.
        set(upperbox,'Visible',' Off');
        set(lowerbox,'Visible',' Off');
        %
        set(halfwave,'visible','off')
        set(hint,'visible','off')
        set(waveedit,'visible','off')
        set(axescurve_hole,'visible','off');
        set(len_cur_local,'visible','off')
        set(uplengthl,'visible','off')
        set(downlengthl,'visible','off')
        set(axes2dshape_hole,'visible','off')
        cla(axes2dshape_hole)
        cla(axescurve_hole)
        set(label_title,'visible','off')
  
    var_length= str2num(get(input1edit,'string'));
    ed_location=str2num(get(input2edit,'string'));
   
 
%     len=length(ed_location);
%     for i=1:len
%     ed_location{i}(:,1)=[];
%     end
%      ed_location_x=str2num(get(input3edit,'string'));
% ed_location=(get(input2edit,'string'));

  ed_location_x=get(input3edit,'string');              
 ed_location_x(:,1)=[];
   
% axessect=axes('Units','normalized','Position',[0.27 0.2 0.22 0.78],'visible','off');
%plot current applied stress
set(label_cs,'string',2); 
      
     [  node_all4plot,cross_section_range,node_all,elem_all,binary_exist ] =hole_fun_general( var_length,ed_location, ed_location_x,elem,node );
     Net_cross_section_plotter4500( axessect,axeslong,cross_section_range{2},var_length,node_all{2}); 
assignin('base', 'binary_exist', binary_exist);

     %      P=1;
%             binary_exist
%        assignin('base', 'A', A);
 [A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,J,Xs,Ys,Cw,B1,B2,w] = cutwp_prop2(node(:,2:3),elem(:,2:4));
[P,M11,M22,B,err] = stress_to_action(node,xcg,zcg,thetap,A,I11,I22,w,Cw);
  maxP=max([P M11 M22 B]);
        if abs(P/maxP)<0.001, P=0;, end
        if abs(B/maxP)<0.001, B=0;, end
        if abs(M11/maxP)<0.001, M11=0;, end
        if abs(M22/maxP)<0.001, M22=0;, end  
%         fy=max(abs(node(:,8))); %guess at fy assuming max stress in the model = fy
%         report the error in the conversion so the user knows and can
%         decide whether or not this is going to be good enough
%         text={['Member actions determined from applied stresses already in model.'];
%             ['Sum squared error in stress across nodes, determined from basing stress'];
%             ['on these actions instead of original stress read from pre, err=',num2str(err)]};
%         choice = questdlg(text,'Accept reference applied actions','OK, Continue','Don''t Use Stress','OK, Continue');
%         switch choice
%             case 'OK, Continue'
%                 %ok...
%             case 'Don''t Use Stress'
%                 return
%         end
%         If a generated action is really small we should set it to zero
%         lets say if the action is less than 0.1% of the yield for that
%         action then we set it to zero, so we need the yield then we can
%         make the corrections as needed
%         find the yield values
%         [Py,Mxxy,Mzzy,M11y,M22y]=yieldMP(node,fy,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,1);
%         [By]=yieldB(fy,Cw,w);
%         correct as needed
%         if abs(P/Py)<0.001, P=0;, end
%         if abs(B/By)<0.001, B=0;, end
%         if abs(M11/M11y)<0.001, M11=0;, end
%         if abs(M22/M22y)<0.001, M22=0;, end  



      [ node4local,elem4local ] =hole_fun_local( node_all,elem_all,P,0,0,M11,M22);
      
     assignin('base', 'elem4local', elem4local);  
     assignin('base', 'node4local', node4local); 
%calculate distance along the cross section to each node("s" coordinate)

      noCS=length(node4local);
      for j=1:noCS
          large_piece=0;
          [maxNode,~]=size(node4local{j});
          [NoElem,~]=size(elem4local{j});
          s_coordinate{j}(1)=0;
for i=2:maxNode
    dx{i}=node4local{j}(i,2)- node4local{j}((i-1),2);
    dz{i}=node4local{j}(i,3)- node4local{j}((i-1),3);
    ds(i)=sqrt(dx{i}^2 + dz{i}^2);
    
    s_coordinate{j}(i,1)=s_coordinate{j}(i-1)+ds(i);
end   
      end
      
                    temp_location=str2num(get(input2edit,'string'));
ed_location=temp_location(:,1:4);
plate_range=temp_location(:,5:6);
      distor_len=str2num(get(input_Lcrd_edit,'string'));

[ elem_distor ] =...
    hole_fun_distor( cross_section_range,elem,distor_len,plate_range,s_coordinate,binary_exist );
flag=0;

      
      
      
 inode=node4local{2};
%        inode(:,8)=0*node4local{1}(:,8);
%strespic(inode,elem,axesstres,scale)
cla(axessect)
scale=1;
 crossect(inode,elem4local{2},axessect,springs,constraints,flags)
 strespic(inode,elem4local{2},axessect,scale)
       assignin('base', 'nnode_all', node_all);
       assignin('base', 'elem_all', elem_all);
     cross_section_range{1};
crossect_hole(inode,elem4local{2},axessect,springs,constraints,flags_hole)
%      set(bl2,'String',num2str(var_length));
%      set(b,'min',0,'max',var_length);
 noCS=length(elem4local);
for j=1:noCS
    start{j}(1)=1;
    [NoElem,~]=size(elem4local{j});
    ct=1;
for i=1:NoElem-1
    if elem4local{j}(i+1,2)~=elem4local{j}(i,3)
        ct=ct+1;
        start{j}(ct)=elem4local{j}(i+1,2);
        last{j}(ct-1)=elem4local{j}(i,3);
    end
    last{j}(ct)=elem4local{j}(NoElem,3);
large_piece(j)=max(s_coordinate{j}(last{j})-s_coordinate{j}(start{j}));
end

crossNo=str2num(get(label_cs,'string')); 
[maxNode,~]=size(node4local{crossNo});
%calculate distance along the cross section to each node("s" coordinate)
node_location= node4local{crossNo}(:,1:3);
node_location(1,4)=0;
[maxNode,~]=size(node4local{crossNo});
ss_coordinate(1)=0;
for i=2:maxNode
    dx{i}=node_location(i,2)- node_location((i-1),2);
    dz{i}=node_location(i,3)- node_location((i-1),3);
    ds(i)=sqrt(dx{i}^2 + dz{i}^2);
    
    ss_coordinate(i,1)=ss_coordinate(i-1)+ds(i);
end 
s_display(:,1)=1:maxNode;
s_display(:,2)=ss_coordinate;
set(S2edit,'String',num2str(s_display));


     set(bl2,'String',num2str(var_length));
     set(b,'min',0,'max',var_length);

end
     
    %-----------------------------------------------
    
    case 11
    %-----------------------------------------------
     crossNo=str2num(get(label_cs,'string')); 
%  hole_occupy=cross_section_range{crossNo}(2:2:end)-cross_section_range{crossNo}(1:2:end-1);



%  large_piece=max(hole_occupy);


% set(waveedit,'String',num2str(lengths_hole{crossNo}));   
% set(hint,'String',['Hint: the longest piece of this CS is ' num2str(large_piece)]);
    
       
noCS=length(node4local);
 for i=1:noCS
     lengths_hole{i}=linspace(large_piece(i)/20,large_piece(i),20);
       
 end
BC_hole='S-S';
%  [filename,pathname]=uiputfile('*','Save as');

  saver_hole(prop,node4local,elem4local,lengths_hole,curve,shapes,springs,constraints,GBTcon,clas,BC_hole,m_all,noCS,elem_distor,node,distor_len);

%     %Grid the raw plastic surface
%     if M11pr==0
%         %let user know they must do the raw construction first
%         errordlg('You must select and build raw surface before interpolation', 'Error - Need raw surface first');
%     else
%         Dtheta=str2num(get(Dthetaedit,'String'));
%         Dphi=str2num(get(Dphiedit,'String'));
%         Plasticsurface=1;
%         [M11p,M22p,Pp,beta_P]=Plastic_Int(M11pr, M22pr, Ppr,Dtheta,Dphi,Plasticsurface);
%         %plot the results
%         rawpts=0;, gridit=1; gridedge=0; gridpts=1;
%         set(rawpointsui,'Value',rawpts);
%         set(gridui,'Value',gridit);
%         set(gridedgeui,'Value',gridedge);
%         set(gridptsui,'Value',gridpts);    
%         %Grid the raw plastic surface
%         betaray=1;
%         thetabetap=str2num(get(Thetaedit,'String'));
%         phibetap=str2num(get(Phiedit,'String'));
%         PMMplasticplotter(axesPMM,M11pr,M22pr,Ppr,rawpts,M11p,M22p,Pp,gridit,gridedge,gridpts,betaray,thetabetap,phibetap)
%     end
    %-----------------------------------------------

    case 12
    %-----------------------------------------------
    %Calculate the Plastic Anchors
    thetaMM=0  ;,PhiPM=90;,[M11ponM11y,~,~,~]=Plastic_Int(M11pr, M22pr, Ppr,thetaMM,PhiPM,0);
    thetaMM=90 ;,PhiPM=90;,[~,M22ponM22y,~,~]=Plastic_Int(M11pr, M22pr, Ppr,thetaMM,PhiPM,0);
    %Get Mpp values
    fy=str2num(get(fyedit,'String'));
    unsymm=1;
    %calculate new values for yield moments
    [A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,J,Xs,Ys,Cw,B1,B2,w] = cutwp_prop2(node(:,2:3),elem(:,2:4));
    [Py,Mxxy,Mzzy,M11y,M22y]=yieldMP(node,fy,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm);
    M11panchor=M11ponM11y*M11y;
    M22panchor=M22ponM22y*M22y;
    Ppanchor=Py;
    %feed the results back to the interface
    set(M11ponM11yedit,'String',num2str(M11ponM11y));
    set(M22ponM22yedit,'String',num2str(M22ponM22y));
    set(M11pedit,'String',num2str(M11panchor));
    set(M22pedit,'String',num2str(M22panchor));
    set(Ppanchoredit,'String',num2str(Ppanchor));    
    %-----------------------------------------------
    
    
    case 13
    %------------------------------------------------
    %Calculate betap
    thetaMM=str2num(get(Thetaedit,'String'));
    phiPM=str2num(get(Phiedit,'String'));
    [~,~,~,betap]=Plastic_Int(M11pr, M22pr, Ppr,thetaMM,phiPM,0);
    set(Betapedit,'String',num2str(betap));
    %and plot the updated beta ray vector
    holehelper_cb(50)
    
    case 50
    %-----------------------------------------------
    rawpts=get(rawpointsui,'Value');
    gridit=get(gridui,'Value');
    gridedge=get(gridedgeui,'Value');
    gridpts=get(gridptsui,'Value');    
    %Grid the raw plastic surface
    betaray=1;
    thetabetap=str2num(get(Thetaedit,'String'));
    phibetap=str2num(get(Phiedit,'String'));
    PMMplasticplotter(axesPMM,M11pr,M22pr,Ppr,rawpts,M11p,M22p,Pp,gridit,gridedge,gridpts,betaray,thetabetap,phibetap)
    
    %-----------------------------------------------
    %----------------------------------------------------
    %viewpoint and rotate button for 3D PMM visualization
    %----------------------------------------------------
   case 51
        axes(axesPMM(1));
        view(0,90)
   case 52
        axes(axesPMM(1));
        view(0,0)
   case 53
        axes(axesPMM(1));
        view(90,0)
   case 54
        axes(axesPMM(1));
        view(37.5,30)
   case 55
        axes(axesPMM(1));
        rotate3d
        %----------------------------------------------------
  
        %generated stress plot options
        %Set the various plotting option flags
        %flags:[node# element# mat# stress# stresspic coord constraints springs origin] 1 means show
        case 204 
            if flags_hole(1)==1, flags_hole(1)=0;, else flags_hole(1)=1;, end
            holehelper_cb(1);
        case 205 
            if flags_hole(2)==1, flags_hole(2)=0;, else flags_hole(2)=1;, end
            holehelper_cb(1);
        case 206 
            if flags_hole(3)==1, flags_hole(3)=0;, else flags_hole(3)=1;, end
            holehelper_cb(1);
        case 207 
            if flags_hole(4)==1, flags_hole(4)=0;, else flags_hole(4)=1;, end
            holehelper_cb(1);
        case 208 
            if flags_hole(5)==1, flags_hole(5)=0;, else flags_hole(5)=1;, end
            holehelper_cb(1);
        case 209 
            if flags_hole(6)==1, flags_hole(6)=0;, else flags_hole(6)=1;, end
            holehelper_cb(1);
        case 210 
            if flags_hole(7)==1, flags_hole(7)=0;, else flags_hole(7)=1;, end
            holehelper_cb(1);
        case 211 
            if flags_hole(8)==1, flags_hole(8)=0;, else flags_hole(8)=1;, end
            holehelper_cb(1);
        case 212 
            if flags_hole(9)==1, flags_hole(9)=0;, else flags_hole(9)=1;, end
            holehelper_cb(1);
        case 213 
            if flags_hole(10)==1, flags_hole(10)=0;, else flags_hole(10)=1;, end
            holehelper_cb(1);
        case 214 
            if flags_hole(11)==1, flags_hole(11)=0;, else flags_hole(11)=1;, end
            holehelper_cb(1);
        case 220
            if strcmp(popanel_hole.Visible,'off')
                popanel_hole.Visible = 'on'
            else
                popanel_hole.Visible = 'off'
            end 
 
    case 900
        helpdlg('The plastic surface is constructing by evaluating a large number of different plastic neutral axis (PNA) locations and determining the P,M11,M22 that results at the selected PNA. All the angles in the angles box are investigated as PNAs and the given angle is investigated at separate locations in the section - grades determines the number of different PNAs that are investigated at a given angle.', 'Help on raw plastic surface construction');
    case 901
        helpdlg('The raw plastic surface is interpolated to an evenly spaced grid of points for viewing and exporting. Select the increments in theta and phi to construct the gridded surface', 'Help on gridding/interpolating plastic surface');

    case 101
        %%
        %Set default initial values
        set(upperbox,'Visible',' Off');
        set(lowerbox,'Visible',' Off');
        %
        set(halfwave,'visible','off')
        set(hint,'visible','off')
        set(waveedit,'visible','off')
        set(axescurve_hole,'visible','off');
        set(len_cur_local,'visible','off')
        set(uplengthl,'visible','off')
        set(downlengthl,'visible','off')
        set(axes2dshape_hole,'visible','off')
        cla(axes2dshape_hole)
        cla(axescurve_hole)
        set(label_title,'visible','off')
        

        
        crossNo=str2num(get(label_cs,'string')); 
%  hole_occupy=cross_section_range{crossNo}(2:2:end)-cross_section_range{crossNo}(1:2:end-1);
%  large_piece=max(hole_occupy);

% lengths_hole{crossNo}=linspace(large_piece/20,large_piece,19);



    maxCS=length(node4local);
    for i=1:maxCS
        lengths_hole{i}=linspace(large_piece(i)/20,large_piece(i),20);
         [curve_local{i},shapes_local{i},length_local{i}] = run_strip_local( prop,node4local{i},elem4local{i},lengths_hole{i});
    end
    
    set(waveedit,'String',num2str(lengths_hole{crossNo}));   
set(hint,'String',['Hint: the longest piece of this CS is ' num2str(large_piece(crossNo))]);


%        [curve_local{crossNo},shapes_local{crossNo},length_local{crossNo}] = run_strip_local( prop,node4local{crossNo},elem4local{crossNo},lengths_hole{crossNo});
% assignin('base', 'curve_local', curve_local);
% assignin('base', 'shapes_local', shapes_local);
% assignin('base', 'length_local', length_local);
%      curve_plotter_hole( prop,node4local{crossNo},elem4local{crossNo},curve_local{crossNo},shapes_local{crossNo},length_local{crossNo},axescurve_hole,axes2dshape_hole);
%      axes2dshape_hole
%      idx_hole_l=1;
[curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4( curve_local{crossNo},shapes_local{crossNo},prop,node4local{crossNo},elem4local{crossNo},crossNo,length_local{crossNo}, axescurve_hole,axes2dshape_hole);
%  cr=0;
%     
%         for m=1:length(curvel(:,1))-2
%             load1=curve_sign(m,2);
%             load2=curve_sign(m+1,2);
%             load3=curve_sign(m+2,2);
%             if (load2<load1)&(load2<=load3)
%                 cr=cr+1;
%                 hold on
%                 idx_hole_l=m+1;
%             end
%         end
set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
label=[' load factor =',num2str(curvel{idx_hole_l}(modeindexl,2))];
set(label_title,'String',label);
set(label_title,'visible','on');
lgd_flag=1;
%
%calcs done show everything
        set(halfwave,'visible','on')
        set(hint,'visible','on')
        set(waveedit,'visible','on')
        set(axescurve_hole,'visible','on');
        set(len_cur_local,'visible','on')
        set(uplengthl,'visible','on')
        set(downlengthl,'visible','on')
        set(label_title,'visible','on')

    case 102
          %%
          %DISTORTIONAL ANALYSIS
        %Set default initial values
        crossNo=str2num(get(label_cs,'string')); 
%  hole_occupy=cross_section_range{crossNo}(2:2:end)-cross_section_range{crossNo}(1:2:end-1);
%  large_piece=max(hole_occupy);
        %Set default initial values
        set(upperbox,'Visible',' Off');
        set(lowerbox,'Visible',' Off');
        %
        set(halfwave,'visible','off')
        set(hint,'visible','off')
        set(waveedit,'visible','off')
        set(axescurve_hole,'visible','off');
        set(len_cur_local,'visible','off')
        set(uplengthl,'visible','off')
        set(downlengthl,'visible','off')
        set(axes2dshape_hole,'visible','off')
        cla(axes2dshape_hole)
        cla(axescurve_hole)
        set(label_title,'visible','off')
  
% lengths_hole{crossNo}=linspace(large_piece/20,large_piece,19);


        temp_location=str2num(get(input2edit,'string'));
% plate_range=temp_location(:,5:6);
%  [ elem_distor ] =hole_fun_distor( cross_section_range,elem,distor_len );
    maxCS=length(node4local);
%     distor_len_sig=linspace(distor_len/5*4,distor_len/5*6,9);
    for i=1:maxCS
       
%          lengths_hole{i}=linspace(large_piece(i)/20,large_piece(i),20);
distor_len_sig=distor_len;
 [ curve_distor{i},shapes_distor{i}] = run_strip_distortional( prop,node,elem_distor{i},distor_len_sig);
    end
   assignin('base', 'curve_distor', curve_distor) 
    
    set(waveedit,'String',num2str(distor_len_sig));   
    
%        [curve_local{crossNo},shapes_local{crossNo},length_local{crossNo}] = run_strip_local( prop,node4local{crossNo},elem4local{crossNo},lengths_hole{crossNo});
% assignin('base', 'curve_local', curve_local);
% assignin('base', 'shapes_local', shapes_local);
% assignin('base', 'length_local', length_local);
%      curve_plotter_hole( prop,node4local{crossNo},elem4local{crossNo},curve_local{crossNo},shapes_local{crossNo},length_local{crossNo},axescurve_hole,axes2dshape_hole);
%      axes2dshape_hole
%      idx_hole_l=1;
[curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4( curve_distor{crossNo},shapes_distor{crossNo},prop,node,elem_distor{crossNo},crossNo,distor_len_sig, axescurve_hole,axes2dshape_hole);
% 
%  cr=0;
%     
%         for m=1:length(distor_len_sig)
% 
%             if distor_len_sig(m)==distor_len
%                 idx_hole_l=m;
%                 m
% break
%                 
%             end
%         end
%         
% set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
% picpointl=[curvel{idx_hole_l}(modeindexl,1) curvel{idx_hole_l}(modeindexl,2)];
%   thecurve3(curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl)
%     
label=[' load factor =',num2str(curvel{idx_hole_l}(modeindexl,2))];
set(label_title,'String',label);
set(label_title,'visible','on');
set(hint,'String',['Distortional results only provided at Lcrd (as defined by user)']);
set(hint,'visible','on')
lgd_flag=2;
%  set(hint,'String',[num2str(idx_hole_l)]);       
%calcs done show everything
        set(hint,'visible','on')
        
%%
case 104
        %%
        %Set default initial values
        crossNo=str2num(get(label_cs,'string')); 
 hole_occupy=cross_section_range{crossNo}(2:2:end)-cross_section_range{crossNo}(1:2:end-1);
 large_piece=max(hole_occupy);

lengths_hole{crossNo}=str2num(get(waveedit,'String'));
set(waveedit,'String',num2str(lengths_hole{crossNo}));   
set(hint,'String',['Hint: the longest piece of this CS is ' num2str(large_piece)]);
    
       [curve_local{crossNo},shapes_local{crossNo},length_local{crossNo}] = run_strip_local( prop,node4local{crossNo},elem4local{crossNo},lengths_hole{crossNo});
% assignin('base', 'curve_local', curve_local);
% assignin('base', 'shapes_local', shapes_local);
% assignin('base', 'length_local', length_local);
%      curve_plotter_hole( prop,node4local{crossNo},elem4local{crossNo},curve_local{crossNo},shapes_local{crossNo},length_local{crossNo},axescurve_hole,axes2dshape_hole);
%      axes2dshape_hole
%      idx_hole_l=1;
[curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4( curve_local{crossNo},shapes_local{crossNo},prop,node4local{crossNo},elem4local{crossNo},crossNo,length_local{crossNo}, axescurve_hole,axes2dshape_hole);


set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
lgd_flag=1;

    case 1044
     pickpoint=get(axescurve_hole,'CurrentPoint'); 
     
     for j=1:max(size(curvel));
    curve_sign(j,1)=curvel{j}(modedisplayl(1),1);
    curve_sign(j,2)=curvel{j}(modedisplayl(1),2);
    if length(modedisplayl)>1
        for mn=2:length(modedisplayl)
            templ(j,modedisplayl(mn))=curvel{j}(modedisplayl(mn),1);
            templf(j,modedisplayl(mn))=curvel{j}(modedisplayl(mn),2);
        end
    end
end
if length(modedisplayl)>1
    [reldiff(1),pickminindex(1)]=min(sqrt((curve_sign(:,1)-pickpoint(1,1)).^2+(curve_sign(:,2)-pickpoint(1,2)).^2));
    for i=2:length(modedisplayl)
        [reldiff(i),pickminindex(i)]=min(sqrt((templ(:,i)-pickpoint(1,1)).^2+(templf(:,i)-pickpoint(1,2)).^2));
    end
    [minreldiff,mindiffindex]=min(reldiff);
    modeindexl=mindiffindex;
    idx_hole_l=pickminindex(modeindexl);
else
    [reldiff,pickminindex]=min(sqrt((curve_sign(:,1)-pickpoint(1,1)).^2+(curve_sign(:,2)-pickpoint(1,2)).^2));
    idx_hole_l=pickminindex;
    % modeindex=1;
end
set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
picpointl=[curvel{idx_hole_l}(modeindexl,1) curvel{idx_hole_l}(modeindexl,2)];
    thecurve3_hole(curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl)
     label=[' load factor=',num2str(curvel{idx_hole_l}(modeindexl,2))];
set(label_title,'String',label);
set(label_title,'visible','on');
% set(len_cur,'String',['length = ',num2str(lengths(lengthindex))]);
% modes=(1:1:length(curvel{idx_hole_l}(:,2)));
% set(mode_cur,'String',num2str(modes(modeindexl)));
        
%%
case 302
      modeindexl=1;
idx_hole_l=idx_hole_l+1;
% modell=shapesl{idx_hole_l}(:,modeindexl);
axes2dshapel=axes2dshape_hole;
% dispshap(undefl,nodel,eleml,modell,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl);
%  set(len_cur_distorl,'String',num2str(lengths_distor(idx_hole_l)));
 
 set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
 
 picpointl=[curvel{idx_hole_l}(modeindexl,1) curvel{idx_hole_l}(modeindexl,2)];
    thecurve3_hole(curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl)
    
    model=shapesl{idx_hole_l}(:,modeindexl);
dispshap(undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl);
 label=[' load factor=',num2str(curvel{idx_hole_l}(modeindexl,2))];
set(label_title,'String',label);
set(label_title,'visible','on');
    case 301


 modeindexl=1;
idx_hole_l=idx_hole_l-1;
% modell=shapesl{idx_hole_l}(:,modeindexl);
axes2dshapel=axes2dshape_hole;
% dispshap(undefl,nodel,eleml,modell,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl);
%  set(len_cur_distorl,'String',num2str(lengths_distor(idx_hole_l)));
 
 set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
 
 picpointl=[curvel{idx_hole_l}(modeindexl,1) curvel{idx_hole_l}(modeindexl,2)];
    thecurve3_hole(curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl)
 model=shapesl{idx_hole_l}(:,modeindexl);
dispshap(undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl);
 label=[' load factor=',num2str(curvel{idx_hole_l}(modeindexl,2))];
label_title_load_ratio=uicontrol(subfig,...
    'Style','text','units','normalized',...
    'Position',[0.580 0.52 0.3 0.03],...
    'String',label);


    case 401
        crossNo=str2num(get(label_cs,'string')); 
        table(crossNo,1)=crossNo;
        table(crossNo,2)=curvel{idx_hole_l}(modeindexl,2);
        
        set(table_hole,'String',num2str(table))

    case 103
        set(axescurve_hole,'visible','off');
        set(halfwave,'visible','off')
        set(hint,'visible','off')
        set(waveedit,'visible','off')
        set(len_cur_local,'visible','off')
        set(uplengthl,'visible','off')
        set(downlengthl,'visible','off')

        set(upperbox,'Visible',' on');
        %uistack(upperbox,'top');
        
        lgd_flag=3;
        set(waveedit,'String',''); 
        [ props4global,netProp ] =hole_fun_global( var_length, node_all,elem_all,cross_section_range);
        assignin('base', 'props4global', props4global);
       assignin('base', 'netProp', netProp );
       assignin('base', 'node_all', node_all );
       assignin('base', 'elem_all', elem_all );
       %%
     crossNo=str2num(get(label_cs,'string'));    
      
        thetap=thetap*180/pi; %degrees...
    Bx=NaN;, By=NaN;
    %initial plot
%     propplot(node,elem,xcg,zcg,thetap,axesprop)
    %
    %
%     if exist(label_title_load_ratio)
%           set(label_title_load_ratio,'Visible',' off');
%     end
%     if
%        set(label_title,'Visible',' off');
%     end


    
     
    %SECTION PROPERTIES
%     upperbox=uicontrol(subfig,...
%         'Style','frame','units','normalized',...
%         'HorizontalAlignment','Left',...
%         'Position',[0.5 0.0 0.5 0.850]);
  
%     upperlabel=uicontrol(subfig,...
%         'Style','text','units','normalized',...
%         'HorizontalAlignment','Center',...
%         'FontName','Arial','FontSize',13,...
%         'Position',[0.01 0.93 .98 0.05],...
%         'String','Calculated Section Properties (Fully Composite Always Assumed)');
    
    %basic sectional properties
%     upperframe=uicontrol(subfig,...
%         'Style','frame','units','normalized',...
%         'Position',[0.0 0.45 0.45 0.45]);
    CS_prop=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Center','fontweight','bold',...
        'Position',[0.01 0.95 0.38 0.05],...
        'String',['Cross-Section ' num2str(crossNo) ' Properties']);
    outA=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.02 0.90 0.05 0.05],...
        'String','A = ');
    outxcg=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.02 0.85 0.05 0.05],...
        'String','xcg = ');
    outzcg=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.25 0.85 0.05 0.05],...
        'String','zcg = ');
    outIxx=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.02 0.80 0.05 0.05],...
        'String','Ixx = ');
    outIzz=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.25 0.80 0.05 0.05],...
        'String','Izz = ');
    outIxz=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.02 0.75 0.05 0.05],...
        'String','Ixz = ');
    outthetap=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.25 0.75 0.05 0.05],...
        'fontname','symbol',...
        'String',('q = '));
    outI11=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.02 0.70 0.05 0.05],...
        'String','I11 = ');
    outI22=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.25 0.70 0.05 0.05],...
        'String','I22 = ');
    outJ=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.01 0.65 0.05 0.05],...
        'String','J = ');
    outC_w=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.25 0.65 0.05 0.05],...
        'String','Cw = ');
    set(out2A,...
        'String',num2str(netProp{crossNo}.A));
  
    set(out2J,...
        'String',num2str(netProp{crossNo}.J));
    set(out2xcg,...
        'String',num2str(netProp{crossNo}.xcg));
    set(out2zcg,...
        'String',num2str(netProp{crossNo}.zcg));
    set(out2Ixx,...
        'String',num2str(netProp{crossNo}.Ixx));
    set(out2Izz,...
        'String',num2str(netProp{crossNo}.Izz));
    set(out2Ixz,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.Ixz));
    set(out2thetap,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.thetap));
    set(out2I11,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.I11));
    set(out2I22,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.I22));
    %
    set(out2C_w,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.Cw));
%     out2X_s=uicontrol(upperbox,...
%         'Style','text','units','normalized',...
%         'HorizontalAlignment','Left',...
%         'Position',[0.07 0.3 0.1 0.05],...
%         'String',num2str(Xs));
%     out2Y_s=uicontrol(upperbox,...
%         'Style','text','units','normalized',...
%         'HorizontalAlignment','Left',...
%         'Position',[0.3 0.3 0.1 0.05],...
%         'String',num2str(Ys));
    
    % Beta properties from CUFSM 3, using CUTWP engine
    % errors corrected in CUTWP engine December 2006 - BWS
    %----------------------------------------------------
%     outB1=uicontrol(upperbox,...
%         'Style','text','units','normalized',...
%         'HorizontalAlignment','Right',...
%         'Position',[0.02 0.15 0.05 0.05],...
%         'fontname','symbol',...
%         'String','b1 = ');
%     outB2=uicontrol(upperbox,...
%         'Style','text','units','normalized',...
%         'HorizontalAlignment','Right',...
%         'Position',[0.02 0.05 0.05 0.05],...
%         'fontname','symbol',...
%         'String','b2 = ');
%     out2B1=uicontrol(upperbox,...
%         'Style','text','units','normalized',...
%         'HorizontalAlignment','Left',...
%         'Position',[0.07 0.15 0.08 0.05],...
%         'String',num2str(B1));
%     out2B2=uicontrol(upperbox,...
%         'Style','text','units','normalized',...
%         'HorizontalAlignment','Left',...
%         'Position',[0.07 0.05 0.08 0.05],...
%         'String',num2str(B2));
    %----------------------------------------------------
    
%     
%     Bas_Adv=uicontrol(upperbox,...
%         'Style','popupmenu','units','normalized',...
%         'Position',[0.15 0.15 0.15 0.05],...
%         'String',['Basic Plot   ';...
%         'Warping Plot '],...
%         'Value',1,...
%         'Callback',[...
%         'propout_cb(6);']);
    
%     scaletitle_w=uicontrol(upperbox,...
%         'Style','text','units','normalized',...
%         'HorizontalAlignment','Right',...
%         'Position',[0.3 0.15 0.03 0.05],...
%         'fontname','symbol',...
%         'String','w');

out_open_prop2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Center','fontweight','bold',...
        'Position',[0.51 0.95 0.38 0.05],...
        'String','Smeared Section Properties');
    outA2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.52 0.90 0.05 0.05],...
        'String','A = ');
    outIxx2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.52 0.80 0.05 0.05],...
        'String','Ixx = ');
    outIzz2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.75 0.80 0.05 0.05],...
        'String','Izz = ');
    outIxz2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.52 0.75 0.05 0.05],...
        'String','Ixz = ');
    outI112=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.52 0.70 0.05 0.05],...
        'String','I11 = ');
    outI222=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.75 0.70 0.05 0.05],...
        'String','I22 = ');
    outJ2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.51 0.65 0.05 0.05],...
        'String','J = ');
    outC_w2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Right',...
        'Position',[0.75 0.65 0.05 0.05],...
        'String','Cw = ');
    
    
    out2A2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'Position',[0.58 0.90 0.12 0.05],...
        'String',num2str(props4global.A));  
    out2Ixx2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'Position',[0.57 0.80 0.12 0.05],...
        'String',num2str(props4global.Ixx));
    out2Izz2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'Position',[0.8 0.80 0.12 0.05],...
        'String',num2str(props4global.Izz));
    out2Ixz2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'Position',[0.57 0.75 0.12 0.05],...
        'String',num2str(props4global.Ixz));
    out2I112=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'Position',[0.57 0.70 0.12 0.05],...
        'String',num2str(props4global.I11));
    out2I222=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'Position',[0.8 0.70 0.12 0.05],...
        'String',num2str(props4global.I22));
    out2J2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'Position',[0.58 0.65 0.12 0.05],...
        'String',num2str(props4global.J));
    out2C_w2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'Position',[0.8 0.65 0.12 0.05],...
        'String',num2str(props4global.Cw));
    %
    %
    %Additional Section properties
    %
%     lowerframe=uicontrol(upperbox,...
%         'Style','frame','units','normalized',...
%         'Position',[0.0 0.0 0.45 0.45]);
    

    %
          set(lowerbox,'Visible',' on');
           KL1=str2num(get(KL1edit,'string')); 
          KL2=str2num(get(KL2edit,'string')); 
             KL3=str2num(get(KL3edit,'string'));
             Pe=get(Pedit_hole,'String');  
M11=get(M11edit_hole,'String');  
M22=get(M22edit_hole,'String');  

         
%%
out_open_prop2=uicontrol(upperbox,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Center','fontweight','bold',...
        'Position',[0.01 0.55 0.98 0.05],...
        'String','CUTWP Driven Buckling Analysis with Smeared Section Properties');

set(lowerbox,'Visible','on');
     uistack(lowerbox,'top');
  W = 963.5; H = 423.5;  
    AW = 384; AH = 240;
rad_origin=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);', ...
        'units','normalized',...
        'Position',[0.01 390/H 0.22 20/H],... 
        'Value',0, ...
        'String','Origin', ...
        'Style','checkbox', ...
        'Tag','rad_origin');
    
    rad_centroid=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);', ...
        'units','normalized',...
        'Position',[0.01 365/H 0.22 20/H],... 
        'Value',0, ...
        'String','Centroid', ...
        'Style','checkbox', ...
        'Tag','rad_centroid');
    
    rad_axisxy=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);', ...
        'units','normalized',...
        'Position',[0.01 340/H 0.22 20/H],... 
        'Value',0, ...
        'String','Axis (x,y)', ...
        'Style','checkbox', ...
        'Tag','rad_axisxy');
    
    rad_axis12=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);', ...
        'units','normalized',...
        'Position',[0.01 315/H 0.22 20/H],... 
        'Value',1, ...
        'String','Axis (1,2)', ...
        'Style','checkbox', ...
        'Tag','rad_axis12');
    
    rad_shear=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);', ...
        'units','normalized',...
        'Position',[0.25 390/H 0.22 20/H],... 
        'Value',1, ...
        'String','Shear Center', ...
        'Style','checkbox', ...
        'Tag','rad_shear');
    
    rad_axial=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);', ...
        'units','normalized',...
        'Position',[0.25 365/H 0.22 20/H],... 
        'Value',1, ...
        'String','Axial Force', ...
        'Style','checkbox', ...
        'Tag','rad_axial');
    
    rad_deform=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);', ...
        'units','normalized',...
        'Position',[0.25 340/H 0.22 20/H],... 
        'Value',1, ...
        'String','Deformed Shape', ...
        'Style','checkbox', ...
        'Tag','rad_deform');
    
    rad_node=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);', ...
        'units','normalized',...
        'Position',[0.25 315/H 0.22 20/H],... 
        'Value',0, ...
        'String','Node & Segment', ...
        'Style','checkbox', ...
        'Tag','rad_node');
    
      rad_Pe=uicontrol(lowerbox, ...
        'CallBack','cutwp_hole_switcher(''Pe'');', ...
        'units','normalized',...
        'Position',[15/W 99.5/H 0.45 20/H],... 
        'Value',1, ...
        'String','Elastic Critical Axial Force, Pe', ...
        'Style','radio', ...
        'Tag','rad_Pe');
    set(rad_Pe,'Value',1);
    rad_Me=uicontrol(lowerbox, ...
        'CallBack','cutwp_hole_switcher(''Me'');', ...
        'units','normalized',...
        'Position',[15/W 79.5/H 0.45 20/H],... 
        'Value',0, ...
        'String','Elastic Critical Moment, Me', ...
        'Style','radio', ...
        'Tag','rad_Me');
    set(rad_Me,'Value',0);
     rad_Me1=uicontrol(lowerbox, ...
        'CallBack','cutwp_hole_switcher(''Me1'');', ...
        'units','normalized',...
        'Position',[15/W 50/H 0.45 20/H],... 
        'Value',0, ...
        'String','Bending about the 1-axis, Me1', ...
        'Style','radio', ...
        'Tag','rad_Me1');
    set(rad_Me1,'Value',0);
    rad_Me2=uicontrol(lowerbox, ...
        'CallBack','cutwp_hole_switcher(''Me2'');', ...
        'units','normalized',...
        'Position',[15/W 30/H 0.45 20/H],... 
        'Value',0, ...
        'String','Bending about the 2-axis, Me2', ...
        'Style','radio', ...
        'Tag','rad_Me2');
    set(rad_Me2,'Value',0);
    rad_Me12=uicontrol(lowerbox, ...
        'CallBack','cutwp_hole_switcher(''Me12'');', ...
        'units','normalized',...
        'Position',[15/W 10/H 0.45 20/H],... 
        'Value',0, ...
        'String','Biaxial Bending, Me1/Me2:', ...
        'Style','radio', ...
        'Tag','rad_Me12');
   set(rad_Me12,'Value',0); 
    ed_c=uicontrol(lowerbox, ...
        'CallBack','holehelper_cb(105);',... 
        'units','normalized',...
        'Position',[0.4 10/H 90/W 20/H],... 
        'String','0', ...
        'Style','edit');
    
   
        slider=uicontrol(lowerbox, ...
        'CallBack','cutwp_hole_switcher(''slider'');', ...
        'units','normalized',...
        'Position',[285/W+0.3 10/H 159/W 17/H],... 
        'Style','slider', ...
        'sliderstep',[0.5 0.5], ...
        'value',1, ...
        'Min',1,'Max',3, ...
        'Tag','slider');
    
        statictext_Pe=uicontrol(lowerbox, ...
            'units','normalized',...
        'Position',[464/W+0.3 20/H 35/W 17/H],... 
        'String','Pe:', ...
        'Style','text', ...
        'Tag','statictext_Pe');
    
            statictext_Me1=uicontrol(lowerbox, ...
            'units','normalized',...
        'Position',[464/W+0.3 30/H 35/W 17/H],... 
        'String','Me1:', ...
        'Style','text');
    
            statictext_Me2=uicontrol(lowerbox, ...
            'units','normalized',...
        'Position',[464/W+0.3 10/H 35/W 17/H],... 
        'String','Me2:', ...
        'Style','text');
    
        text_Pe=uicontrol(lowerbox, ...
                'units','normalized',...
        'Position',[495/W+0.3 20/H 133.5/W 20/H],... 
        'String','15.509', ...
        'Style','edit',...
        'Tag','text_Pe');  
            text_Me1=uicontrol(lowerbox, ...
                'units','normalized',...
        'Position',[495/W+0.3 30/H 133.5/W 20/H],... 
        'String','15.509', ...
        'Style','edit'); 
            text_Me2=uicontrol(lowerbox, ...
                'units','normalized',...
        'Position',[495/W+0.3 10/H 133.5/W 20/H],... 
        'String','15.509', ...
        'Style','edit'); 
     
    
        text_mode=uicontrol(lowerbox, ...
            'units','normalized',...
        'Position',[385/W+0.3 30/H 20/W 17/H], ... 
        'HorizontalAlignment','right', ...
        'String','1', ...
        'Style','text', ...
        'Tag','text_mode');
    
    text_maxmode=uicontrol(lowerbox, ...
        'units','normalized',...
        'Position',[405/W+.3 30/H 30/W 17/H], ... 
        'String','/3', ...
        'Style','text', ...
        'Tag','text_maxmode');
    
%         ed_c=uicontrol(...
%         'CallBack','cutwp_hole_switcher(''Me12'');',... 
%         'Position',[175/W+0.3 10/H 90/W 20/H],... 
%         'BackgroundColor','w',...
%         'String','0', ...
%         'Style','edit');
%     
    t_axes = axes('Parent',lowerbox,'Position',[0.5 0.15 0.49 0.8]); 
    holehelper_cb(105);
    %%
     case 105
           W = 963.5; H = 423.5;  
    AW = 384; AH = 240;
coord=node(:,2:3);
   ends=elem(:,2:4);
    KL1 = str2num(get(KL1edit,'String'));
    KL2 = str2num(get(KL2edit,'String'));
    KL3 = str2num(get(KL3edit,'String'));
   
%     AllHandleList = get(lowerbox,'Children');
%      findobj('Tag', 'rad_Pe')
    if get(rad_Pe,'Value');
        forcee = 'Pe';
        set(rad_Me1,'Visible','off');
        set(rad_Me2,'Visible','off');
        set(rad_Me12,'Visible','off');   
        set(text_Pe,'Visible','on');
        set(text_Me1,'Visible','off');
        set(text_Me2,'Visible','off');
        set(statictext_Pe,'Visible','on','string','Pe:');
        set(statictext_Me1,'Visible','off');
        set(statictext_Me2,'Visible','off');                        
        set(statictext_KL1,'Visible','on');
        set(statictext_KL2,'Visible','on');
        set(statictext_exy1,'Visible','on');
        set(statictext_exy2,'Visible','on');
        set(ed_KL1,'Visible','on');
        set(ed_KL2,'Visible','on');
        set(ed_exy,'Visible','on');
        set(ed_c,'Visible','off');
        set(text_maxmode,'String','/3');
        set(slider,'Max',3,'sliderstep',[0.5 0.5]);
        set(rad_Me1,'Value',0);
        set(rad_Me2,'Value',0);
        set(rad_Me12,'Value',0);
         set(rad_Me,'Value',0);
    elseif get(rad_Me,'Value');
        set(rad_Me1,'Visible','on');
        set(rad_Me2,'Visible','on');
        set(rad_Me12,'Visible','on'); 
        set(statictext_exy1,'Visible','off');
        set(statictext_exy2,'Visible','off');
        set(ed_exy,'Visible','off');
        set(ed_c,'Visible','on');
        set(text_maxmode,'String','/2');
        set(slider,'Max',2,'sliderstep',[1 1]);
        if get(rad_Me1,'Value');
            forcee = 'Me1';
            set(text_Pe,'Visible','on');
            set(text_Me1,'Visible','off');
            set(text_Me2,'Visible','off');
            set(statictext_Pe,'Visible','on','string','Me1:');
            set(statictext_Me1,'Visible','off');
            set(statictext_Me2,'Visible','off');             
            set(statictext_KL1,'Visible','off');
            set(statictext_KL2,'Visible','on');
            set(ed_KL1,'Visible','off');
            set(ed_KL2,'Visible','on');
        elseif get(rad_Me2,'Value');
            forcee = 'Me2';
            set(text_Pe,'Visible','on');
            set(text_Me1,'Visible','off');
            set(text_Me2,'Visible','off');
            set(statictext_Pe,'Visible','on','string','Me2:');
            set(statictext_Me1,'Visible','off');
            set(statictext_Me2,'Visible','off');     
            set(statictext_KL1,'Visible','on');
            set(statictext_KL2,'Visible','off');
            set(ed_KL1,'Visible','on');
            set(ed_KL2,'Visible','off');
        elseif get(rad_Me12,'Value');
            forcee = 'Me12';
            set(text_Pe,'Visible','off');
            set(text_Me1,'Visible','on');
            set(text_Me2,'Visible','on');
            set(statictext_Pe,'Visible','off');
            set(statictext_Me1,'Visible','on');
            set(statictext_Me2,'Visible','on');              
            set(statictext_KL1,'Visible','on');
            set(statictext_KL2,'Visible','on');
            set(ed_KL1,'Visible','on');
            set(ed_KL2,'Visible','on');
            
        end
    end
 c = str2num(get(ed_c,'String'));   
    % check input data   
%     coord = str2num(get(ed_coord,'String'));
%     ends = str2num(get(ed_ends,'String'));
%     KL1 = str2num(get(ed_KL1,'String'));
%     KL2 = str2num(get(ed_KL2,'String'));
%     KL3 = str2num(get(ed_KL3,'String'));
%     exy = str2num(get(ed_exy,'String'));

%     c = str2num(get(ed_c,'String'));
%     E = str2num(get(ed_E,'String'));
%     v = str2num(get(ed_v,'String'));
%     dist = str2num(get(ed_dist,'String'));
    mode = str2num(get(text_mode,'string'));
%     
%     if (~strcmp(class(coord),'double'))|isempty(coord)|(size(coord,2)~=2)
%         msgbox('Please redefine Node Data','Error Message','error')
%         return
%     elseif (~strcmp(class(ends),'double'))|isempty(ends)|(size(ends,2)~=3)|(size(ends,1)==1)|(size(coord,1)<max(max(ends(:,1:2))))|(1>min(min(ends(:,1:2))))|any(ends(:,3)<0)
%         msgbox('Please redefine Element Data','Error Message','error')
%         return
%     elseif (~strcmp(class(KL1),'double'))|isempty(KL1)|(KL1 <= 0)|(size(KL1,2)~=1)|(size(KL1,1)~=1)
%         msgbox('Please redefine KL1','Error Message','error')
%         return
%     elseif (~strcmp(class(KL2),'double'))|isempty(KL2)|(KL2 <= 0)|(size(KL2,2)~=1)|(size(KL2,1)~=1)
%         msgbox('Please redefine KL2','Error Message','error')
%         return
%     elseif (~strcmp(class(KL3),'double'))|isempty(KL3)|(KL3 <= 0)|(size(KL3,2)~=1)|(size(KL3,1)~=1)
%         msgbox('Please redefine KL3','Error Message','error')
%         return
%     elseif (~strcmp(class(exy),'double'))|isempty(exy)|(size(exy,2)~=2)|(size(exy,1)~=1)
%         msgbox('Please redefine Eccentricity coordinate','Error Message','error')
%         return
%     elseif (~strcmp(class(c),'double'))|isempty(c)|(size(c,2)~=1)|(size(c,1)~=1)
%         msgbox('Please redefine Me1/Me2','Error Message','error')
%         return
%     elseif (~strcmp(class(E),'double'))|isempty(E)|(E <= 0)|(size(E,2)~=1)|(size(E,1)~=1)
%         msgbox('Please redefine Elastic Modulus','Error Message','error')
%         return
%     elseif (~strcmp(class(v),'double'))|isempty(v)|(size(v,2)~=1)|(size(v,1)~=1)
%         msgbox('Please redefine Poisson''s ratio','Error Message','error')
%         return
%     elseif (~strcmp(class(dist),'double'))|isempty(dist)|(size(dist,2)~=1)|(size(dist,1)~=1)
%         msgbox('Please redefine Displacement Factor','Error Message','error')
%         return
%     end       
    
%     if get(rad_Me12,'Value');
%         [A,xc,yc,Ix,Iy,Ixy,theta,I1,I2,J,xs,ys,Cw,B1,B2,Pe,dcoord] = cutwp_prop(coord,ends,KL1,KL2,KL3,forcee,c,E,v,dist);
%     else
%         [A,xc,yc,Ix,Iy,Ixy,theta,I1,I2,J,xs,ys,Cw,B1,B2,Pe,dcoord] = cutwp_prop(coord,ends,KL1,KL2,KL3,forcee,exy,E,v,dist);
%     end

   [xc,yc,Ix,Iy,Ixy,theta,I1,I2,J,xs,ys,~,B1,B2,Pe,dcoord] = cutwp_prop_4hole(node, elem,KL1,KL2,KL3,prop,props4global,forcee,c);
    assignin('base', 'Pe', Pe); 
%     set(text_A,'string',sprintf('%-13.6G',A));
%     set(text_Ix,'string',sprintf('%-13.6G',Ix));
%     set(text_Iy,'string',sprintf('%-13.6G',Iy));
%     set(text_Ixy,'string',sprintf('%-13.6G',Ixy));
%     set(text_I1,'string',sprintf('%-13.6G',I1));
%     set(text_I2,'string',sprintf('%-13.6G',I2));
%     set(text_theta,'string',sprintf('%-13.6G',theta));
%     set(text_J,'string',sprintf('%-13.6G',J));
%     set(text_centroid,'string',sprintf('( %6.6G, %6.6G )',xc,yc)); 
    
%     if isnan(Cw)
%         set(statictext_J,'string','Torsion Constant-Closed Section, J:');
%         set(statictext_shear,'Visible','off'); 
%         set(statictext_Cw,'Visible','off'); 
%         set(statictext_B1,'Visible','off'); 
%         set(statictext_B2,'Visible','off'); 
%         set(statictext_Pe,'Visible','off'); 
%         set(statictext_Me1,'Visible','off');
%         set(statictext_Me2,'Visible','off');
%         set(text_shear,'Visible','off'); 
%         set(text_Cw,'Visible','off');
%         set(text_B1,'Visible','off');
%         set(text_B2,'Visible','off');
%         set(text_Pe,'Visible','off');
%         set(text_Me1,'Visible','off');
%         set(text_Me2,'Visible','off');
%         set(rad_shear,'Visible','off','Value',0);
%         set(rad_axial,'Visible','off','Value',0);
%         set(rad_deform,'Visible','off','Value',0);
%         set(slider,'Visible','off');
%         set(statictext_dist,'Visible','off');
%         set(statictext_buck,'Visible','off');
%         set(text_mode,'Visible','off');
%         set(text_maxmode,'Visible','off');
%         set(ed_dist,'Visible','off');
%         
%     else
%         
%         set(statictext_J,'string','Torsion Constant-Open Section, J:');
%         set(statictext_shear,'Visible','on'); 
%         set(statictext_Cw,'Visible','on'); 
%         set(statictext_B1,'Visible','on'); 
%         set(statictext_B2,'Visible','on'); 
%         
%         set(text_shear,'Visible','on','string',sprintf('( %6.6G, %6.6G )',xs,ys));
%         set(text_Cw,'Visible','on','string',sprintf('%-13.6G',Cw));
%         set(text_B1,'Visible','on','string',sprintf('%-13.6G',B1));
%         set(text_B2,'Visible','on','string',sprintf('%-13.6G',B2));
        
        set(rad_shear,'Visible','on');
        if get(rad_Pe,'Value')|get(rad_Me1,'Value')|get(rad_Me2,'Value')
            set(statictext_Pe,'Visible','on'); 
            set(statictext_Me1,'Visible','off');
            set(statictext_Me2,'Visible','off');
            set(text_Pe,'Visible','on','string',sprintf('%-13.6G',Pe(mode)));
            set(text_Me1,'Visible','off');
            set(text_Me2,'Visible','off');
        elseif get(rad_Me12,'Value')
            set(statictext_Pe,'Visible','off'); 
            
            set(statictext_Me1,'Visible','on');
            set(statictext_Me2,'Visible','on');
            set(text_Pe,'Visible','off');
            set(text_Me1,'Visible','on','string',sprintf('%-13.6G',Pe(mode)*c));
            set(text_Me2,'Visible','on','string',sprintf('%-13.6G',Pe(mode)));
        end
        if get(rad_Pe,'Value')
            set(rad_axial,'Visible','on');
        else
            set(rad_axial,'Visible','off','Value',0);
        end
        
        set(rad_deform,'Visible','on');
        set(slider,'Visible','on');
%         set(statictext_dist,'Visible','on');
%         set(statictext_buck,'Visible','on');
        set(text_mode,'Visible','on');
        set(text_maxmode,'Visible','on');
%         set(ed_dist,'Visible','on');
        
            origin = get(rad_origin,'Value');
    centroid = get(rad_centroid,'Value');
    axisxy = get(rad_axisxy,'Value');
    axis12 = get(rad_axis12,'Value');
    shear = get(rad_shear,'Value');
    axial = get(rad_axial,'Value');
    deform = get(rad_deform,'Value');
    nnode = get(rad_node,'Value');
    
        AW = 384; AH = 240;
%     im = load('imcutwp');
%      
%     image(im.im); 
%     axis off;
%     clear im;
    cla(t_axes); 
   axes(t_axes);
    exy=[0 0];
    cutwp_draw(coord,ends,exy,xc,yc,theta,xs,ys,dcoord,origin,centroid,axisxy,axis12,shear,axial,deform,nnode,mode);    
%       set(label_title_load_ratio,'Visible',' off');
% if exist('label_title')
%        set(label_title,'Visible',' off');
% end
  
    case 106
        totalCS=length(cross_section_range);
        L=(var_length*(get(range_slider,'value')));
        set(bl3,'String',num2str(L))
        assignin('base', 'cross_section_range', cross_section_range);
        assignin('base', 'L', L);
        assignin('base', 'totalCS', totalCS);
for i=1:totalCS
    for j=1:length(cross_section_range{i})/2
    if (L>=cross_section_range{i}(2*j-1) && L<=cross_section_range{i}(2*j))
              crossNo=i;  
        break  
    end
    end
end
        
       
              set(label_cs,'string',num2str(crossNo)); 
              crossNo=str2num(get(label_cs,'string')); 
      Net_cross_section_plotter4500( axessect,axeslong,cross_section_range{crossNo},var_length,node_all4plot{crossNo});
%       set(label_cs,'string',num2str(crossNo)); 
      cla(axessect)
      cla(axescurve_hole)
crossect_hole(node4local{crossNo},elem4local{crossNo},axessect,springs,constraints,flags_hole)
 isEmpytAxes = isempty(get(axes2dshape_hole, 'Children'));
  if isEmpytAxes~=1
     if lgd_flag==1
  [curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4( curve_local{crossNo},shapes_local{crossNo},prop,node4local{crossNo},elem4local{crossNo},crossNo,length_local{crossNo}, axescurve_hole,axes2dshape_hole);
     elseif lgd_flag==2
         [curvel,curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl...
   ,lengths_local,idx_hole_l,shapesl,shapescelll...
   ,undefl,nodel,eleml,model,axes2dshapel,scalel,springsl,m_al,BCl,SurfPosl]= plotter_fun_batchcufsm4(curve_distor{crossNo},shapes_distor{crossNo},prop,node,elem_distor{crossNo},crossNo,distor_len_sig, axescurve_hole,axes2dshape_hole);
     elseif lgd_flag==3
            set(CS_prop,...
        'String',['Cross-Section ' num2str(crossNo) ' Properties']);  
        set(out2A,...
        'String',num2str(netProp{crossNo}.A));
  
    set(out2J,...
        'String',num2str(netProp{crossNo}.J));
    set(out2xcg,...
        'String',num2str(netProp{crossNo}.xcg));
    set(out2zcg,...
        'String',num2str(netProp{crossNo}.zcg));
    set(out2Ixx,...
        'String',num2str(netProp{crossNo}.Ixx));
    set(out2Izz,...
        'String',num2str(netProp{crossNo}.Izz));
    set(out2Ixz,...
        'String',num2str(netProp{crossNo}.Ixz));
    set(out2thetap,...
        'String',num2str(netProp{crossNo}.thetap));
    set(out2I11,...
        'String',num2str(netProp{crossNo}.I11));
    set(out2I22,...
        'String',num2str(netProp{crossNo}.I22));
    set(out2C_w,...
        'String',num2str(netProp{crossNo}.Cw));
%        for m=1:length(distor_len_sig)
% 
%             if distor_len_sig(m)==distor_len
%                 idx_hole_l=m;
%                 m
% break
%                 
%             end
%         end
%         
% set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
% picpointl=[curvel{idx_hole_l}(modeindexl,1) curvel{idx_hole_l}(modeindexl,2)];
%   thecurve3_hole(curvecelll,filenamecelll,clascelll,filedisplayl,minoptl,logoptl,clasoptl,axescurvel,xminl,xmaxl,yminl,ymaxl,modedisplayl,fileindexl,modeindexl,picpointl)
     end
%      set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
 label=[' load factor=',num2str(curvel{idx_hole_l}(modeindexl,2))];
set(label_title,'String',label);
set(label_title,'visible','on');
set(len_cur_local,'String',num2str(lengths_local(idx_hole_l)));
 end
 [maxNode,~]=size(node4local{crossNo});
%calculate distance along the cross section to each node("s" coordinate)
node_location= node4local{crossNo}(:,1:3);
node_location(1,4)=0;
[maxNode,~]=size(node4local{crossNo});
ss_coordinate(1)=0;
for i=2:maxNode
    dx{i}=node_location(i,2)- node_location((i-1),2);
    dz{i}=node_location(i,3)- node_location((i-1),3);
    ds(i)=sqrt(dx{i}^2 + dz{i}^2);
    
    ss_coordinate(i,1)=ss_coordinate(i-1)+ds(i);
end 
s_display(:,1)=1:maxNode;
s_display(:,2)=ss_coordinate;
set(S2edit,'String',num2str(s_display));
if lgd_flag==1
set(hint,'String',['Hint: the longest piece of this CS is ' num2str(large_piece(crossNo))]);
end
if lgd_flag==3
    
                set(CS_prop,...
        'String',['Cross-Section ' num2str(crossNo) ' Properties']);  
set(out2A,...
        'String',num2str(netProp{crossNo}.A));
  
    set(out2J,...
        'String',num2str(netProp{crossNo}.J));
    set(out2xcg,...
        'String',num2str(netProp{crossNo}.xcg));
    set(out2zcg,...
        'String',num2str(netProp{crossNo}.zcg));
    set(out2Ixx,...
        'String',num2str(netProp{crossNo}.Ixx));
    set(out2Izz,...
        'String',num2str(netProp{crossNo}.Izz));
    set(out2Ixz,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.Ixz));
    set(out2thetap,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.thetap));
    set(out2I11,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.I11));
    set(out2I22,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.I22));
    %
    set(out2C_w,...
        'Style','text','units','normalized',...
        'HorizontalAlignment','Left',...
        'String',num2str(netProp{crossNo}.Cw));
          set(label_title_load_ratio,'Visible',' off');
      set(label_title,'Visible',' off');
 end
       


        end
    
    

% if screen==0
%     pre2;
% end
% %if lengths empty (since it is done in a separate step)
% if isempty(lengths)
%     lengths(1)=1;
% end
% %save key variable to temp file in cutwp directory
% if ispc %pc
%     save([currentlocation,'\cutwp\fromcufsm'],'prop','node','elem','lengths');
% else %mac! or unix
%     save([currentlocation,'/cutwp/fromcufsm'],'prop','node','elem','lengths');
% end
% % cutwp_hole;
% cutwp_hole('FromCUFSM');
% %     fig1 = get(gcf,'children');
% %     copyobj(fig1,lowerbox);
%        
end
        
        
