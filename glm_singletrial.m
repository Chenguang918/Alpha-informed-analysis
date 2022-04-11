function [beta,a1,a2] = glm_singletrial(y, onset1, onset2, duration1, duration2, modulation1, modulation2, plotindication, stat, samplingFrequency)
% function [beta, T, pvalue, covb] = glm(y, onset, duration, modulation, plotindication, doPassFilter, samplingFrequency)
% This is to do glm estimation on NIRS data
% 
% y: NIRS data, a matrix of size Nxc. Where N is the sample size and c is # of channels
% onset: onset timing (in seconds) of each type of events. It is a cell array, each component is for one type event
% duration: (optional) duration of event in seconds. It should have exactly the same
% size of onset. Duration 0 means punctuated event. Default duration is 0
% modulation: (optional) modeulation of event. It should have exactly the same size of onset. Default is 1
% plotindication (optional) 1 by default. if plot result
% stats (optional) 0:model have no constant  1:model have constant 
% samplingFrequency (optional), 10 Hz by default. The sampling frequency of
% your NIRS signal. If you are using ETG4000, it is probably 10 Hz.
%
% The returned variables are matrix of dimension cxn, where c is # of
% channels and n is the number of conditions you specified in onset.
%
% beta: beta coefficient
% T: T value of beta
% pvalue: pvalue of beta
% covb: covariance matrix of beta. This can be used to calculate the T value of
% conrast. The variance of contrast = c' * covb * c
%
% Standard SPM hrf is used to convole the event
%
% Example:
% [beta, T, pvalue, covb] = glm(hbdata{1}{1}, {onset}, {duration}, {[]});
% 
% ZhaoChenGuang
% 2017/02/16




N = size(y, 1); % sample size

% passfilter hb data
if ~exist('doPassFilter')
    doPassFilter = 1;
end

% now convolve
M = length(onset1); % # of conditions
x1 = zeros(N, M);   % standard
x2 = zeros(N, M);   %
% x1
for ii=1:M
    x1(:,ii) = eventTime2Delta(onset1{ii}, N, modulation1{ii}, duration1{ii}, 1, 1/samplingFrequency);
end
% figure
% plot(x1)
a1=x1;
for ii=1:M
    X2(:,ii) = eventTime2Delta(onset2{ii}, N, modulation2{ii}, duration2{ii}, 1, 1/samplingFrequency);
end
% figure
% plot(X2)
a2=X2;
%
switch stat
    case 0
var=[x1 X2];%
    case 1
var=[ones(size(x1)) x1 X2];%
end
for ii=1:size(y,2) % loop for each column
      HbO=y(:,ii);
      [b,bint,r,rint,stats]=regress(HbO,var);
      beta(:,ii)=b;
      %Bint(:,ii)=bint;
      %R(:,ii)=r;
      %Rint(:,ii)=rint;
      %R2(ii,:)=stats(1,1);
      %P(ii,:)=stats(1,3);
end
%
if ~exist('plotindication')
    plotindication = 0;%%%%modified by HJ default is 1
end
if(plotindication)    
    if(size(T,1) == 48)
        figure;
        plotTopoMap(T(1:24,1), '4x4', [min(T(:)) max(T(:))]);
        figure;
        plotTopoMap(T(25:48,1), '4x4', [min(T(:)) max(T(:))]);
    end
    if(size(T,1) == 44)
        figure;
        plotTopoMap(T(1:22,1), '3x5', [min(T(:)) max(T(:))]);
        figure;
        plotTopoMap(T(23:44,1), '3x5', [min(T(:)) max(T(:))]);
    end
    if(size(T,1) == 52)
        figure;
        plotTopoMap(T(1:end,1), '3x11', [min(T(:)) max(T(:))]);     
    end
    %plotTopoMap(beta(1:24), '4x4', [min(beta(:)) max(beta(:))])
    %plotTopoMap(beta(25:48), '4x4', [min(beta(:)) max(beta(:))])
end
end
