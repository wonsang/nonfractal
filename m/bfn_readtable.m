function Y = bfn_readtable(filename, varargin)
% BFN_READTABLE  Read a matrix from a file
%
% Syntax:
%   Y = bfn_readtable(filename,'Property Name 1','Property value 1',...)
%
% Description:
%   It reads a matrix from a file with predefined options. All elements
%   are recognized as real numbers.
%
% Input Arguments:
%   filename : the file name with or without path
%
% Options:
%   delimiter : the separator of elements in the matrix.
%       Defaults is ';'.
%   headerlines : the number of header lines. Default is 1.
%   skipcols : the number of skipped left columns. Default is 1.
%
% Output Arguments:
%   Y : an rows x cols matrix
%
% References:
%   Orfanidis, S.J., Optimum Signal Processing. An Introduction. 2nd
%    Edition, Prentice-Hall, Englewood Cliffs, NJ, 1996.
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% $ Id: bfn_readtable.m 0028 2013-03-12 00:13:18 brainfnet $
% Copyright (c) 2011 Brainfnet.

params = struct('delimiter'     ,';',...
                'headerlines'   ,1,...
                'skipcols'      ,1);         
params = bfn_parseArgs(varargin,params);

% count the number of rows
fid = fopen(filename,'r');
fseek(fid, 0, 'eof');
chunksize = ftell(fid);
fseek(fid, 0, 'bof');
ch = fread(fid, chunksize, '*uchar');
nrows = sum(ch == sprintf('\n')); % number of lines 
fclose(fid);

% read data
X       = textread(filename,'%s','delimiter', params.delimiter, ...
                'headerlines',params.headerlines);
rows    = nrows - params.headerlines;
cols    = size(X,1)/rows;
Y       = reshape(textread(filename,'%s','delimiter', params.delimiter, ...
                'headerlines',params.headerlines), ...
            cols, rows)';
Y       = Y(:,params.skipcols+1:end);
Y       = str2double(Y);
