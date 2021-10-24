function [Data_allblock]=step4_remove_trial_select_func(ERP_MEASURES,nsub,sub)
std_delt = 2;%%%%motion artifact criteria, outline std
    subnum=['sub',num2str(sub(nsub))];
    in_path = 'G:\Paper2\Alpha-informed analysis\HbO_interpodata\';
    load('G:\Paper2\Alpha-informed analysis\42sub_triinfo.mat')
    in_name = [subnum,'_Interpodata.mat'];
    load([in_path,in_name]);
    %%
    M_data=ERP_MEASURES(:,nsub);
    M_type=MI.M_type(:,nsub);
    M_type(find(M_type==0))=[];
    %%
        datastand=[1:9,13,28:39,43,44];%
        art_data=mean(Data_allblock.HbO(:,datastand),2);%select data always good &using average to do artifact
        art_marker=Data_allblock.markerTime(1:2:end);
        trial_keep = s_artifact(art_data,art_marker,std_delt,40);%%%%remove motion. 40,前后20个point内的marker
        tem_nTrial = length(trial_keep);
        Data_allblock.markerTimeafter(1:2:tem_nTrial*2)=Data_allblock.markerTime(trial_keep*2-1); %%%%%before remove motion, only kept trials marker
        Data_allblock.markerTimeafter(2:2:tem_nTrial*2)=Data_allblock.markerTime(trial_keep*2);%%%%%before remove motion, only kept trials marker
        Data_allblock.markerTxt(1:2:tem_nTrial*2) = Data_allblock.markertype( trial_keep*2-1);
        Data_allblock.markerTxt(2:2:tem_nTrial*2) = Data_allblock.markertype( trial_keep*2);
        Data_allblock.trial_keep = trial_keep';
        Data_allblock.markerTimeafter = Data_allblock.markerTimeafter';
        Data_allblock.markerTxt = Data_allblock.markerTxt';
        Data_allblock.type=[11;21;12;22;13;23];
        ori_marker=Data_allblock.markertype;
        clear trial_keep
        trial_keep = Data_allblock.trial_keep;
        oricue_marker=Data_allblock.markertype(1:2:end);
        trial_all = [1:length(Data_allblock.markertype)./2];%%%each trial have 2 events;
%%     
        if length(M_type)~=length(oricue_marker)
        step=10; 
        [match_MI,match_M_data]=Timeline_match(oricue_marker,M_type,M_data,step);
        M_data=match_M_data;
        indx_miss_trial=find(match_MI==0);
        [~,ai,~]=intersect(trial_keep,indx_miss_trial);
        trial_keep(ai)=[];
        end
%% 
             j = 1; i = 1;
             str_rmv = [];%%%%%the bad marker
             stp_rmv = [];%%%%%the first good marker next to bad marker
             trial_out = setdiff(trial_all,trial_keep);%get trial to remove
             str_rmv(1) = trial_out(1);
             while i < length(trial_out)
                 if j > length(str_rmv) %%%% the pre trial is not continuous
                     str_rmv(j) = trial_out(i);%%%%j is number of trials not continuous
                 end
                 if trial_out(i+1) - trial_out(i) > 1%%%%remove trials are not continuous
                      stp_rmv(j) = trial_out(i) + 1;%%%the trial next to the bad trial
                      j = j + 1;
                 end
                 i = i + 1;
             end
             
             if j == length(str_rmv) %%%the last remove trial is next to the last second trial
                 stp_rmv(j) = trial_out( i ) + 1;
             else %%%%the last trial is not continous 
                 str_rmv(j) = trial_out( i );
                 if trial_out( i )==length(trial_all) 
                     stp_rmv(j) = trial_out( i );
                 else
                     stp_rmv(j) = trial_out( i )+1;
                 end 
             end
              Data_allblock.nRmv = j;%the num of remove time
              Data_allblock.strRmvTrial = str_rmv;% Remove a section of bad trial from the first trial  
              Data_allblock.stpRmvTrial = stp_rmv;%The trial after the last bad trial is the first good trial  
              Data_allblock.strRmvTime = Data_allblock.markerTime(str_rmv*2-1);%%%each trial have 2 events
              Data_allblock.stpRmvTime = zeros(size(Data_allblock.strRmvTime));
              Data_allblock.trial_out=trial_out';
              
              if stp_rmv(end) > length(Data_allblock.markerTime)./2%%%%the last trial is removed
                  Data_allblock.stpRmvTime(1:j-1) = Data_allblock.markerTime(stp_rmv(1:j-1)*2-1)-1;
                  Data_allblock.stpRmvTime(j) = length(Data_allblock.HbO);
              else
                  Data_allblock.stpRmvTime = Data_allblock.markerTime(stp_rmv*2-1)-1;%%%%the point before the good trial
              end
               Data_allblock.nRmvPoint = Data_allblock.stpRmvTime - Data_allblock.strRmvTime + 1;
               points_keep = [1:30];
               Data_allblock.markerTimeRem=Data_allblock.markerTimeafter;
               if trial_keep(1)==1%%%%the first trial is removed
                   points_keep = [1 : Data_allblock.strRmvTime(1)-1];
               end            
               ii = find(Data_allblock.markerTimeafter>=Data_allblock.strRmvTime(1));
               Data_allblock.markerTimeRem(ii) = Data_allblock.markerTimeRem(ii) - Data_allblock.nRmvPoint(1);
               for i = 2 : j
                   points_keep = [points_keep,Data_allblock.stpRmvTime(i-1)+1:Data_allblock.strRmvTime(i)-1];
                   ii = find(Data_allblock.markerTimeafter>=Data_allblock.strRmvTime(i));
                   Data_allblock.markerTimeRem(ii) = Data_allblock.markerTimeRem(ii) - Data_allblock.nRmvPoint(i);
               end  
               if Data_allblock.stpRmvTime(end) < length(Data_allblock.HbO)
                   points_keep = [points_keep,Data_allblock.stpRmvTime(end)+1:length(Data_allblock.HbO)];
               end
               Data_allblock.HbO_aftRem = Data_allblock.HbO(points_keep,:);
               Data_allblock.points_keep = points_keep';
               nStim=length(Data_allblock.type);

%% HbO epoch
                a11=find(ori_marker==11);
                Cue.Cue11=(a11+1)*0.5;
                RevCue.Cue11=intersect(Cue.Cue11,trial_keep);
                M.Cue11=M_data(RevCue.Cue11);
                Timecue.cue11=art_marker(RevCue.Cue11);
                Timecuelog.cue11=art_marker(RevCue.Cue11)+1;    
                
                
                a21=find(ori_marker==21);
                Cue.Cue21=(a21+1)*0.5;
                RevCue.Cue21=intersect(Cue.Cue21,trial_keep);
                M.Cue21=M_data(RevCue.Cue21);
                Timecue.cue21=art_marker(RevCue.Cue21);
                Timecuelog.cue21=art_marker(RevCue.Cue21)+1;   
                
                a12=find(ori_marker==12);
                Cue.Cue12=(a12+1)*0.5;
                RevCue.Cue12=intersect(Cue.Cue12,trial_keep);
                M.Cue12=M_data(RevCue.Cue12);
                Timecue.cue12=art_marker(RevCue.Cue12);
                Timecuelog.cue12=art_marker(RevCue.Cue12)+1;   
                
                a22=find(ori_marker==22);
                Cue.Cue22=(a22+1)*0.5;
                RevCue.Cue22=intersect(Cue.Cue22,trial_keep);
                M.Cue22=M_data(RevCue.Cue22);
                Timecue.cue22=art_marker(RevCue.Cue22);
                Timecuelog.cue22=art_marker(RevCue.Cue22)+1;   
                
                a13=find(ori_marker==13);
                Cue.Cue13=(a13+1)*0.5;
                RevCue.Cue13=intersect(Cue.Cue13,trial_keep);
                M.Cue13=M_data(RevCue.Cue13);
                Timecue.cue13=art_marker(RevCue.Cue13);
                Timecuelog.cue13=art_marker(RevCue.Cue13)+1;   
                
                a23=find(ori_marker==23);
                Cue.Cue23=(a23+1)*0.5;
                RevCue.Cue23=intersect(Cue.Cue23,trial_keep);
                M.Cue23=M_data(RevCue.Cue23);
                Timecue.cue23=art_marker(RevCue.Cue23);
                Timecuelog.cue23=art_marker(RevCue.Cue23)+1;   
                Data_allblock.M_data=M;
i=[];component=[];ERP_marker=[];
 for cueIn=1:6
                
                name=['Timecue.cue',num2str(Data_allblock.type(cueIn))];
                eval(['nTrials=length(',name,');']);
                eval(['Data_allblock.regressor1(1:nTrials,cueIn)=',name,';'])%The time points corresponding to all trials of each type  
                namelog=['Timecuelog.cue',num2str(Data_allblock.type(cueIn))];              
                eval(['Data_allblock.regressor2(1:nTrials,cueIn)=',namelog,';'])%The time points corresponding to all trials of each type  
                Data_allblock.nTrials(cueIn) = nTrials;%the trial of everystimtype         
 end
%%
    index=[1:1:points_keep(end)];
    difindex=setdiff(index,points_keep);
    Data_allblock.oriHbO_after=Data_allblock.HbO;
    Data_allblock.oriHbO_after(difindex,:)=0;
    Data_allblock.stdch=datastand;
    %save([fpathout,subnum,'_remmotiondata.mat'],'Data_allblock');
end