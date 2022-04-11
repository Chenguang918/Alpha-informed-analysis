function t_keep = s_artifact(data,marker,nStd,tLength)
% % % data = mean(d_blocks(1,1).HbO(:,1:22),2);
% % % marker = d_blocks(1,1).markerTime_raw - d_blocks(1,1).strp +1;
% % % nStd = 3;
% % % tLength = 40;
%%%%%this function is to detect movement and clear marker about these
%%%%%movement
%%%%data, a vector,signal,
%%%%marker,a vector, marker onset
%%%%artifact standard, delet data beyond standard
%%%%tLength, trial length for artifact
%%%%t_keep, marker indices that kept

mean_d = mean(data);%%%%mean value of data
std_d = std(data);%%%%% std value of data
upper_d = mean_d + nStd * std_d;%%%%%the max to keep
down_d = mean_d - nStd * std_d;%%%%the min to keep


right_out = find(data >= upper_d);%%%beyond the max indices
left_out = find(data <= down_d);%%%beyond the min indices

%%%%delet marker that near left_ and right_ out

if isempty(right_out) && isempty(left_out)
    t_marker = marker;
    t_raw = [1:length(marker)];%%%%%trial sequence,1,2,3,....
    t_keep = t_raw;
else
    t_diff1 = [];
    t_diff2 = [];
    t_marker = [];
    t_keep = [];
    t1 = 0;
    t2 = 0;
    for i = 1 : length(marker)
        if ~isempty(right_out)
            t_diff1 = abs(marker(i) - right_out);
            t1 = length(find(t_diff1<=tLength./2));
        end
        if ~isempty(left_out)
            t_diff2 = abs(marker(i) - left_out);
            t2 = length(find(t_diff2<=tLength./2));
        end
        
        %%%%if t1 =0 and t2 = 0, then this marker is far from movement
        if t1 == 0 && t2 == 0
            t_marker = [t_marker,marker(i)];
            t_keep = [t_keep,i];
        end
    end      
end
% % % figure
% % % plot(data)
% % % hold on
% % % plot(marker,data(marker),'r+')
% % % plot(t_marker,data(t_marker),'b+')
% % % hold off



