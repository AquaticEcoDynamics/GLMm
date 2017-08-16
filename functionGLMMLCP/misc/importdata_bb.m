function [flddata] = importdata_bb(filename)
%Written by Brendan Busch 2 February 2015
% Imitates import data used when using a single command for reading in
% files with alternative delimiters.
% Customised for Lake Analyzer files used in the GLM MLCP

fid= fopen(filename,'rt');

% Headers
header1 = fgetl(fid);

headers = strsplit(header1,'temp');
%For MATLAB 2012 string split followed by string
%headers = strsplit('temp',header1);

headers = regexprep(headers(2:end),'\s','');
headers = regexprep(headers,',','');
headers = regexprep(headers,',','');

%Find if delimiter is space, tab or comma
sss = strfind(header1,',');

frewind(fid)
% read single line: number of x-values
if ~isempty(sss) %Comma deliminated
    textformat = ['%s',repmat('%f',1,length(headers))];
    datacell = textscan(fid,textformat,'Headerlines',1,'Delimiter',',');
    mdates = datenum(datacell{1});
    for i = 1:length(headers)
	   data(:,i) = datacell{i+1};
    end
else
    textformat = ['%s%s',repmat('%f',1,length(headers))];
    datacell = textscan(fid,textformat,'Headerlines',1);
    spaces(1:length(datacell{1}),1) = ' ';
    dates = [char(datacell{1}) spaces char(datacell{2})];
    mdates = datenum(dates);
    for i = 1:length(headers)
	   data(:,i) = datacell{i+2};
    end
end
fclose(fid);
flddata.headers = str2num(char(headers));
flddata.mdates = mdates;
flddata.data = data;