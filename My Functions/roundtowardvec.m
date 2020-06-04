function X=roundtowardvec(X,roundvec,type)
%function newnums=roundtowardvec(X,[roundvec],[type])
%
% This function rounds number(s) toward given values. If more than one
% number is given to round, it will return the matrix with each rounded
% value, otherwise it will return the single rounded value. It will ignore
% NaNs and return them back with NaNs.
%
% Inputs: X: the number(s) that you want rounded
%
%         roundvec:(opt) the values to round X to. If none given, it will
%           default to -inf:1:inf (and use the built in functions).
%
%         type:(opt) specifies which kind of rounding you want
%           the function to use.
%
%           Choices are: 'round' - round to nearest value
%                        'floor' - round toward -Inf
%                        'ceil'  - round toward Inf
%                        'fix'   - round toward 0
%                        'away'  - round away from 0 (ceil if positive and floor if negative)
%                     (see help files for more clarity)
%
%           If no type is given, the function will default to rounding to
%           the nearest value.
%
% Outputs: newnums: rounded values, in same shape as X input matrix
%          indices: indices of rounded values in roundvec

if nargin==0
	help roundtowardvec; %if nothing given, tell what to give
	return
elseif isempty(X)
	%if given empty, return empty without going through whole script
	return
end

if nargout>1
	error('Too many output variables are given');
end
if ~exist('type','var') || isempty(type)
	type='round';  %%round to nearest value if not specified
end
if ~exist('roundvec','var') || isempty(roundvec) || all(isnan(roundvec))
	if strcmpi(type,'round')
		%to nearest integer
		X=round(X);
	elseif strcmpi(type,'away')
		%nearest integer away from 0
		X=ceil(abs(X)).*sign(X);
	elseif strcmpi(type,'fix')
		%nearest integer toward 0
		X=fix(X);
	elseif strcmpi(type,'floor')
		%nearest integer toward -inf
		X=floor(X);
	elseif strcmpi(type,'ceil')
		%nearest integer toward inf
		X=ceil(X);
	else
		error('%sRound type not recognized. Options are:\n''round'' - round to nearest value\n''floor'' - round toward -Inf\n''ceil''  - round toward Inf\n''fix''   - round toward 0\n''away''  - round away from 0','')
	end
else
	%Ignore nan in roundvec
	roundvec(isnan(roundvec))=[];
	
	%Record which values are nan to ignore
	Xnan=isnan(X);
	
	%Hold onto size for returning value
	sz=size(X);
	
	%Calculate differences
	X=X(:);
	roundvec=roundvec(:)';
	diffs=bsxfun(@minus,X,roundvec);
	
	if strcmpi(type,'round') %to nearest value
		[~,inds]=min(abs(diffs),[],2);
		X=roundvec(inds);
	elseif strcmpi(type,'fix') %to nearest value toward 0
		
		iless=X<0;
		X(iless)=roundtowardvec(X(iless),roundvec,'ceil');
		X(~iless)=roundtowardvec(X(~iless),roundvec,'floor');
	elseif strcmpi(type,'ceil') %nearest value toward inf
		diffs(diffs>0)=nan;
		[~,inds]=min(abs(diffs),[],2);
		
		i_inf=X>max(roundvec);
		X=roundvec(inds);
		X(i_inf)=inf;
	elseif strcmpi(type,'floor') %nearest value toward -inf
		diffs(diffs<0)=nan;
		[~,inds]=min(abs(diffs),[],2);
		
		i_inf=X<min(roundvec);
		X=roundvec(inds);
		X(i_inf)=-inf;
	elseif strcmpi(type,'away') %nearest value away from 0
		
		iless=X<0;
		X(~iless)=roundtowardvec(X(~iless),roundvec,'ceil');
		X(iless)=roundtowardvec(X(iless),roundvec,'floor');
	else
		error('%sRound type not recognized. Options are:\n''round'' - round to nearest value\n''floor'' - round toward -Inf\n''ceil''  - round toward Inf\n''fix''   - round toward 0\n''away''  - round away from 0','')
	end
	
	%Return to output side
	X=reshape(X(:),sz);
	
	%Ignore nan in input dataset
	X(Xnan)=nan;
end


end