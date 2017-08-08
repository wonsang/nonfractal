function [A, B] = bfn_trim_nonchar(C)
% BFN_TRIM_NONCHAR  Trimming the beginning non-char elements from a cell.
%
% Syntax:
%   [A, B] = bfn_trim_nonchar(C)
%
% Description:
%   It separates a cell into both a cell containing successive non-char
%   elements and a cell whose first element has the char format.
%
% Input Arguments:
%   C - a cell.
%
% Output Arguments:
%   A - a cell which contains successive non-char elements
%   A - a cell whose first element is formatted as char.
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_trim_nonchar.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

A = {}; B = {};

N = length(C);
for i = 1:N
    if ischar(C{i})
        B = C(i:end);
        return
    else
        A{i} = C{i};
    end    
end
