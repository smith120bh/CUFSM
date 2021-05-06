function [pathname,filename]=saver_hole(prop,node,elem,lengths,curve,shapes,springs,constraints,GBTcon,clas,BC,m_all,noCS,elem_distor,node_distor,distor_len)
%BWS
%August 2000 (last modified)
%J Cai. 2016
%
 [filename,pathname]=uiputfile('*','Save as');
if filename==0
	return
else
    %if isempty(curve)|curve(1,1)==0|isempty(shapes)|shapes(1,1)==0
%     if exist('curve')<1 & exist('shapes')<1
%         save([pathname,filename],'prop','node','elem','lengths','springs','constraints','GBTcon','BC','m_all');
%     else
for i=1:noCS
        if exist('clas')<1
            for j=1:length(lengths{i})
    m_all{j}=[1];
            end
             S.prop=prop;
            S.node=node{i};
            S.elem=elem{i};
            S.lengths=lengths{i};
            S.springs=springs;S.constraints=constraints;S.GBTcon=GBTcon;S.curve=curve;S.shapes=shapes;S.BC=BC;S.m_all=m_all;
             Sd.prop=prop;
            Sd.node=node_distor;
            Sd.elem=elem_distor{i};
            Sd.lengths=distor_len;
            Sd.springs=springs;Sd.constraints=constraints;Sd.GBTcon=GBTcon;Sd.curve=curve;Sd.shapes=shapes;Sd.BC=BC;Sd.m_all=m_all;           
            save([pathname,[ filename '_L' num2str(i)] '.mat'],'-struct', 'S');
            save([pathname,[filename '_D' num2str(i) '.mat']],'-struct', 'Sd');
%             save('Data.mat', '-struct', 'S')
%             save([pathname,[ filename '_L' num2str(i)] '.mat'],'prop','node','elem{i}','lengths{i}','springs','constraints','GBTcon','curve','shapes','BC','m_all');
%             save([pathname,[filename '_D' num2str(i) '.mat']],'prop','node_distor','elem_distor{i}','distor_len','springs','constraints','GBTcon','curve','shapes','BC','m_all');
        elseif exist('clas')==1
            for j=1:length(lengths{i})
    m_all{j}=[1];
            end
            S.prop=prop;
              S.node=node{i};
            S.elem=elem{i};
            S.lengths=lengths{i};
            S.springs=springs;S.constraints=constraints;S.GBTcon=GBTcon;S.curve=curve;S.shapes=shapes;S.BC=BC;S.m_all=m_all;
             S.clas=clas;
              Sd.prop=prop;
            Sd.node=node_distor;
            Sd.elem=elem_distor{i};
            Sd.lengths=distor_len;
            Sd.clas=clas;
            Sd.springs=springs;Sd.constraints=constraints;Sd.GBTcon=GBTcon;Sd.curve=curve;Sd.shapes=shapes;Sd.BC=BC;Sd.m_all=m_all;           
            save([pathname,[ filename '_L' num2str(i)] '.mat'],'-struct', 'S');
            save([pathname,[filename '_D' num2str(i) '.mat']],'-struct', 'Sd');
%             save([pathname,[filename '_L' num2str(i) '.mat']],'prop','node{i}','elem{i}','lengths{i}','springs','constraints','GBTcon','curve','shapes','clas','BC','m_all');
%             save([pathname,[ filename '_D' num2str(i) '.mat']],'prop','node_distor','elem_distor{i}','distor_len','springs','constraints','GBTcon','curve','shapes','BC','m_all');
            
        end
end
end


