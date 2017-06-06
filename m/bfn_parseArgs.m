function [ArgStruct ArgStructNonField] = bfn_parseArgs(args, ArgStruct, varargin)
% BFN_PARSEARG2  Parsing input arguments for a function
%
% Syntax:
%   ArgStruct = bfn_parseArgs(args, ArgStruct)
%   ArgStruct = bfn_parseArgs(args, ArgStruct, 'unrecognized')
%
% Description:
%   It parses input arguments for a function. It is an
%   advanced version of the function parseArgs in the Matlab package 
%   for crosswavelet and wavelet coherence analysis which is provided 
%   by Aslak Grinsted in
%   http://www.pol.ac.uk/home/research/waveletcoherence/.
%
% Input Arguments:
%   ArgStruct - the structure full of named arguments with default values.
%   Flagtype - parameters that don't require a value. 
%     (the value will be set to 1 if it is present)
%      Aliases can be used to map one argument-name to several 
%      argstruct fields.
%
% Options:
%   unrecognized - whether to recognize any argument which does not belong
%                   to ArgStruct.
%
% Output Arguments:
%   ArgStruct - the updated structure of arguments
%
% Examples:
%   Define the acceptable named arguments and assign default values
%   Args=struct('Holdaxis',0, ...
%        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
%        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
%        'rows',[],'cols',[]); 
%
%   The capital letters define abrreviations.  
%   Eg. parseargtest('spacingvertical',0) is equivalent to  parseargtest('sv',0) 
%
%   Args=bfn_parseArgs(varargin,Args, ... % fill the arg-struct with values entered by the user
%           {'Holdaxis'}, ... %this argument has no value (flag-type)
%           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
%
%   disp(Args)
%
% References:
%   Orfanidis, S.J., Optimum Signal Processing. An Introduction. 2nd
%   Edition, Prentice-Hall, Englewood Cliffs, NJ, 1996.
%
%__________________________________________________________________________
% Original from Aslak Grinsted, Modified by Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_parseArgs.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

%% main
ArgStructNonField = {};
if ~isempty(args)
    if isstruct(args{1})
        ArgStruct = parseArgs(args, ArgStruct, varargin);    
    else
        [ArgStructNonField, args] = bfn_trim_nonchar(args);
        if ~isempty(args)
            if iscell(args{1})
                ArgStruct = parseArgs(args{1}, ArgStruct, varargin);
            elseif iscell(args)
                ArgStruct = parseArgs(args, ArgStruct, varargin);
            end
        end
    end
end

end

%% the function of parsing arguments
function ArgStruct = parseArgs(args, ArgStruct, varargin)

persistent matlabver

if isempty(matlabver)
    matlabver=ver('MATLAB');
    matlabver=str2double(matlabver.Version);
end

Aliases={};
FlagTypeParams='';

if (length(varargin)>0) 
    FlagTypeParams=lower(strvcat(varargin{1}));  %#ok
    if length(varargin)>1
        Aliases=varargin{2};
    end
end


%---------------Get "numeric" arguments
NumArgCount=1;
while (NumArgCount<=size(args,2))&&(~ischar(args{NumArgCount}))
    NumArgCount=NumArgCount+1;
end
NumArgCount=NumArgCount-1;
if (NumArgCount>0)
    ArgStruct.NumericArguments={args{1:NumArgCount}};
else
    ArgStruct.NumericArguments={};
end 


%--------------Make an accepted fieldname matrix (case insensitive)
Fnames=fieldnames(ArgStruct);
for i=1:length(Fnames)
    name=lower(Fnames{i,1});
    Fnames{i,2}=name; %col2=lower
    Fnames{i,3}=[name(Fnames{i,1}~=name) ' ']; %col3=abreviation letters (those that are uppercase in the ArgStruct) e.g. SpacingHoriz->sh
    %the space prevents strvcat from removing empty lines
    Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams,'exact')); %Does this parameter have a value?
end
FnamesFull=strvcat(Fnames{:,2}); %#ok
FnamesAbbr=strvcat(Fnames{:,3}); %#ok

if length(Aliases)>0  
    for i=1:length(Aliases)
        name=lower(Aliases{i,1});
        FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(name,FnamesFull,'exact'); %&??????? exact or not? 
        end
        Aliases{i,2}=FieldIdx;
        Aliases{i,3}=[name(Aliases{i,1}~=name) ' ']; %the space prevents strvcat from removing empty lines
        Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
    end
    %Append aliases to the end of FnamesFull and FnamesAbbr
    FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); %#ok
    FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3})); %#ok
end

%--------------get parameters--------------------
if ~isempty(args)
    if isstruct(args{1})
        a            = args{1};
        a_fieldnames = fieldnames(a);
        numFnames = size(FnamesFull,1) - 1;
        for i = 1:numFnames
            Fname       = strtrim(FnamesFull(i,:));
            FieldIdx    = strmatch(Fname,lower(a_fieldnames),'exact');

            if ~isempty(FieldIdx) 
                if matlabver>=6
                    val = a.(a_fieldnames{FieldIdx,1});
                    ArgStruct.(a_fieldnames{FieldIdx,1}) = val; %try the line below if you get an error here
                else
                    val = getfield(a,a_fieldnames{FieldIdx,1});
                    ArgStruct = setfield(ArgStruct,a_fieldnames{FieldIdx,1},val); %#ok <-works in old matlab versions
                end
                a = rmfield(a, a_fieldnames{FieldIdx,1});
            end            
        end
        
        % recognize unspecified fields
        ar = a;
        ar_fieldnames = fieldnames(ar);
        numRemainFnames = length(ar_fieldnames);        
        for FieldIdx = 1:numRemainFnames            
            if strcmpi(FlagTypeParams,'unrecognized')
                if matlabver>=6
                    val = ar.(ar_fieldnames{FieldIdx,1});
                    ArgStruct.(ar_fieldnames{FieldIdx,1}) = val; %try the line below if you get an error here
                else
                    val = getfield(ar,ar_fieldnames{FieldIdx,1});
                    ArgStruct = setfield(ArgStruct,ar_fieldnames{FieldIdx,1},val); %#ok <-works in old matlab versions
                end
            else
                fprintf('Warning: Unrecognized parameter: %s\n', ar_fieldnames{FieldIdx,1});
            end
        end

    else

        l=NumArgCount+1; 
        while (l<=length(args))
            a=args{l};
            if ischar(a)
                paramHasValue=1; % assume that the parameter has is of type 'param',value
                a=lower(a);
                FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
                if isempty(FieldIdx) 
                    FieldIdx=strmatch(a,FnamesFull,'exact'); 
                end
                if (length(FieldIdx)>1) %shortest fieldname should win 
                    [mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));%#ok
                    FieldIdx=FieldIdx(mxi);
                end
                if FieldIdx>length(Fnames) %then it's an alias type.
                    FieldIdx=Aliases{FieldIdx-length(Fnames),2}; 
                end

                if ~isempty(FieldIdx) 
                    for curField=FieldIdx' %if it is an alias it could be more than one.
                        if (Fnames{curField,4})
                            if (l+1>length(args))
                                error(['Expected a value for parameter: ' Fnames{curField,1}])
                            end
                            val=args{l+1};
                        else %FLAG PARAMETER
                            if (l<length(args)) %there might be a explicitly specified value for the flag
                                val=args{l+1};
                                if isnumeric(val)
                                    if (numel(val)==1)
                                        val=logical(val);
                                    else
                                        error(['Invalid value for flag-parameter: ' Fnames{curField,1}])
                                    end
                                else
                                    val=true;
                                    paramHasValue=0; 
                                end
                            else
                                val=true;
                                paramHasValue=0; 
                            end
                        end
                        if matlabver>=6
                            ArgStruct.(Fnames{curField,1})=val; %try the line below if you get an error here
                        else
                            ArgStruct=setfield(ArgStruct,Fnames{curField,1},val); %#ok <-works in old matlab versions
                        end
                    end
                else
                    if strcmpi(FlagTypeParams,'unrecognized')
                        if (l+1>length(args))
                            error(['Expected a value for parameter: ' a])
                        end
                        val=args{l+1};
                        if matlabver>=6
                            ArgStruct.(args{l})=val; %try the line below if you get an error here
                        else
                            ArgStruct=setfield(ArgStruct,args{l},val); %#ok <-works in old matlab versions
                        end
                    else                        
                        fprintf('Warning: Unrecognized parameter: %s\n', args{l});
                    end
                end
                l=l+1+paramHasValue; %if a wildcard matches more than one
            else
                error(['Expected a named parameter: ' num2str(a)])
            end 
        end

    end
end

ArgStruct = rmfield(ArgStruct, 'NumericArguments');

end
