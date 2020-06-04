
function [bins counts] = interspike_histogram(spkTr1, spkTr2, maxInt, varargin)

% This function calsulates the inter-spike interval histogram
%
%       [hist] = spk_crosscorr(spktr1, spktr2, maxInterval);
%
% Max interval is in ms.
%
% Or with options ...
%
%       [rtn] = spk_crosscorr(spktr1, spktr2, maxInterval, 'paramName', paramValue);
%
% 'trialDur', trialLength   - Trial length, in Sec. Y Axis will be spikes/s,
%                             rather than counts.
% 'divisions', N            - N divisions in histogram. Default = 50. (For
%                             one side of histogram only).
% 'hAxis', axisHandle       - Pass an axis handle to draw into.
% 'yLabel', 1 or 0          - Draw or not Y axis labels (default: draw)
% 'plot', 1 or  0           - Draw a figure;


% Input Options %
opt.hAxis = [];
opt.trialDur = [];
opt.divisions = 50; %default 50
opt.yLabel = 1;
opt.drawPlot = 1;
for ii=1:2:length(varargin)
    opt.(varargin{ii}) = varargin{ii+1};
end

% Work in ms %
spkTr1 = spkTr1.*1000;
spkTr2 = spkTr2.*1000;

% Find all Intervals (N-(N+1), N-(N+2), N-(N+3) etc .. ) %
nSpk = length(spkTr1);
int = repmat(0,[nSpk nSpk-1]);
for ii=1:nSpk-1                         % Missing out last iteration prevents 0 intervals (N-N).
    ind = [nSpk-ii+1:nSpk 1:nSpk-ii];   % Index to shift spikes to N+nSpk position (wrapping around first-last spikes)
    spkTr2Shift = spkTr2(ind);
    intTemp = spkTr1 - spkTr2Shift;
    int(:,ii) = intTemp;
end

% Bin %
binwidth = maxInt/opt.divisions;
bins = -maxInt:binwidth:maxInt;
if isempty(spkTr1) || isempty(spkTr2)
    % This is necessary for graceful failure when nSpk=0 %
    counts = repmat(NaN, size(bins));
else
    % This is normal %
    counts = histc(reshape(int,1,[]),bins);
end
if ~isempty(opt.trialDur)
    counts = (counts ./ opt.trialDur);
end

