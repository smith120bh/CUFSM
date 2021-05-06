function [pathname,filename,prop,node,elem,lengths,curve,shapes,springs,constraints,GBTcon,clas,BC,m_all]=loader();
%BWS
%August 2000 (last modified)
%
%
[filename,pathname]=uigetfile('*.mat','Select File for use in CUFSM');
if filename==0
	return
else
	prop=[];
	node=[];
	elem=[];
	lengths=[];
	curve=[];
	shapes=[];
	springs=[];
	constraints=[];
    GBTcon.glob=[];
    GBTcon.dist=[];
    GBTcon.local=[];
    GBTcon.other=[];
    clas=[];

	load([pathname,filename]);
    if exist('m_all')<1
        for i=1:length(lengths)
            m_all{i}=[1];
        end
    end
    if exist('BC')<1
        BC='S-S';
    end
%     if exist('GBTcon.glob')<1
%         GBTcon.glob=0;
%     end
%     if exist('GBTcon.dist')<1
%         GBTcon.dist=0;
%     end
%     if exist('GBTcon.local')<1
%         GBTcon.local=0;
%     end
%     if exist('GBTcon.other')<1
%         GBTcon.other=0;
%     end
    if exist('GBTcon.ospace')<1
        GBTcon.ospace=1;
    end
    if exist('GBTcon.orth')<1
        GBTcon.orth=2;
    end
    if exist('GBTcon.couple')<1
        GBTcon.couple=1;
    end
    if exist('GBTcon.norm')<1
        GBTcon.norm=1;
    end
end

