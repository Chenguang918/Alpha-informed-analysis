function delta = eventTime2Delta(eventTime, totalLength, modulation, duration, toConv, tr)
% function delta = eventTime2Delta(eventTime, totalLength, modulation, duration, toConv, tr)
%
% eventTime: a vector which contains the time of events (in second)
% totalLength: (optional) how many data points (in scan)? for example, if
% tr=2, then totalLength of 50 is indeed 50x2=100 seconds.
% modulation: (optional) modulator of hrf height. should be same length with eventTime. 
% duration: (optional) duration of events (in second)
% toConv: (optional) conv with hrf? 
%     1: (default), convolve with standard hrf
%     'h':    same with 1
%     'g':    convolve with gamma function
%     0:  don't convolve
%     '': don't convolve
% tr: (optional) 2s by default
%
% delta: the result sequence
%
% Example:
%   eventTime = [40, 80, 123];
%   d = eventTime2Delta(eventTime, 150, [1, 1.5, 0.8]);
%   plot(d);
%   figure;plot(eventTime2Delta([0 50], [100], [1 2], [30 0], 'g', 2))
%
% Xu Cui
% 2005/05/27 created
% 2005/12/28 modified (normalize)
% 2006/6/2 modified (normalize again)

if ~exist('tr')
    tr = 2; % second
end

if ~exist('toConv')
    toConv = 1;
end

if ~exist('totalLength')
    totalLength = round(max(eventTime)/tr);
end

if ~exist('duration')
    durationorheight = eventTime*0;
    duration = eventTime*0;
end
    
if exist('modulation')
    if isempty(modulation)
        modulation = ones(size(eventTime));
    elseif length(modulation) ~= length(eventTime)
        error('modulator length should be equal to eventTime length');
        return;
    end
else
    modulation = ones(size(eventTime));
end

if isempty(eventTime)
    delta = zeros(1, totalLength);
    return
end

delta = zeros(1,max(round(totalLength*tr/0.01), round(max(eventTime)/0.01+1)));
for ii=1:length(eventTime)
    if duration(ii)==0
        amplitude = modulation(ii)*100;
    else
        amplitude = modulation(ii);
    end
	delta(round(eventTime(ii)/0.01)+1:round((eventTime(ii)+duration(ii))/0.01)+1) = amplitude;
end

if ~isempty(toConv)
    if toConv == 'h' | toConv == 1
        delta = conv(delta, spm_hrf(0.01));
        %delta = delta/1.9254e-3;
    elseif toConv == 'g'
        delta = conv(delta, gammahrf([0:0.01:15]));
        %delta = delta/100;
    end
end

for ii=1:totalLength
    delta2(ii) = delta(round((ii-1)*tr/0.01)+1);
end

delta = delta2;