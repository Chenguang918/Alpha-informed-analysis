%% glm fit
function [Data_allblock]=step5_glmfit_SGO_func(duration2,log,stat,Data_allblock,nsub,sub)
% path_out=['D:\Zhaochenguang\Paper2\Alpha-informed analysis\GLMfit\'];
% path_in=['D:\Zhaochenguang\Paper2\Alpha-informed analysis\Remmotiondata\'];

% for nsub=1:length(sub)
    subnum=['sub',num2str(sub(nsub))];
%     load([path_in,subnum,'_remmotiondata.mat']);
    nStim=length(Data_allblock.type);
    rCount=0;
    durations1=0.2;   %stimulus event
    durations2=duration2;   %alpha-informed event
    Data_allblock.sub=subnum;
    Data_allblock.M1={};
    Data_allblock.M2={};
    Data_allblock.D1={};
    Data_allblock.D2={};
    Data_allblock.onsets1 ={};
    Data_allblock.onsets2 ={};
%% SGO 
    for rIn = 1 : nStim
            rCount = rCount + 1;
            type{rIn}=['cond',num2str(Data_allblock.type(rIn))];
            typelog{rIn}=['condlog',num2str(Data_allblock.type(rIn))];
            eval([type{rIn},'=nonzeros(Data_allblock.regressor1(:,',num2str(rIn),'));']);
            eval([typelog{rIn},'=nonzeros(Data_allblock.regressor2(:,',num2str(rIn),'))+log;']);
            Data_allblock.M1(1,rCount)={ones(1,eval(['length(',type{rIn},')']))} ;                   % standard design matrix
            Data_allblock.M2(1,rCount)={eval(['Data_allblock.M_data.Cue',num2str(Data_allblock.type(rIn))])} ;  
            Data_allblock.D1(1,rCount)={ones(1,eval(['length(',type{rIn},')']))*durations1} ;
            Data_allblock.D2(1,rCount)={ones(1,eval(['length(',type{rIn},')']))*durations2} ;
            Data_allblock.onsets1(1,rCount) ={eval([type{rIn},'/Data_allblock.fs'])};
            Data_allblock.onsets2(1,rCount) ={eval([typelog{rIn},'/Data_allblock.fs'])};
            Data_allblock.P(1,rCount)={pinv(cell2mat(Data_allblock.onsets1(1,rCount)))};
            Data_allblock.y=Data_allblock.oriHbO_after;  
            % Schmidt-Gram orthogonalization
            Xn=Data_allblock.M2;
            Xg=Data_allblock.onsets1;
            XG=Data_allblock.P;
            Xnn=((-1)^rIn)*cell2mat(Xn(1,rCount));
            % clear
            Xnc=(Xnn-mean(Xnn))/std(Xnn);
            outlier=find(abs(Xnc)>2.5);%The extreme value returns to zero
            Xnc(outlier)=0;
            Xgg=cell2mat(Xg(1,rCount));
            XGG=cell2mat(XG(1,rCount));
            XNN=Xnc-Xgg*XGG*Xnc;    %Xnn
            Data_allblock.XN(1,rCount)={XNN};
    end
%%    
for xi = 1 : nStim   
[beta,a1,a2] =glm_singletrial(Data_allblock.y,Data_allblock.onsets1(1,xi),Data_allblock.onsets1(1,xi),Data_allblock.D1(1,xi),Data_allblock.D2(1,xi),Data_allblock.M1(1,xi),Data_allblock.XN(1,xi), 0, 0, Data_allblock.fs);        
r=corr(a1,a2);
Data_allblock.beta1(:,:,xi) = beta;%beta value
Data_allblock.r1(:,xi) = r;%beta value
end  
%%
for xi = 1 : 3   
 onsets1=[cell2mat(Data_allblock.onsets1(1,2*xi-1)); cell2mat(Data_allblock.onsets1(1,2*xi))];
 DD1=[cell2mat(Data_allblock.D1(1,2*xi-1)) cell2mat(Data_allblock.D1(1,2*xi))];
 DD2=[cell2mat(Data_allblock.D2(1,2*xi-1)) cell2mat(Data_allblock.D2(1,2*xi))];
 MM1=[cell2mat(Data_allblock.M1(1,2*xi-1)) cell2mat(Data_allblock.M1(1,2*xi))];
 MM2=[cell2mat(Data_allblock.XN(1,2*xi-1)); cell2mat(Data_allblock.XN(1,2*xi))];
[beta,a1,a2] =glm_singletrial(Data_allblock.y,{onsets1},{onsets1},{DD1'},{DD2'},{MM1'},{MM2}, 0, stat, Data_allblock.fs);        
r=corr(a1,a2);
Data_allblock.beta3(:,:,xi) = beta;%beta value
Data_allblock.r3(:,xi) = r;%beta value
end  
%     save([path_out,subnum,'_glmfit.mat'],'Data_allblock');
%     disp(subnum)
