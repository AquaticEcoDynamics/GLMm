function [nml] = getGLMnml(GLMnmlFile)

% function [nml_text] = getGLMnml(GLMnmlFile)
%
% Read a GLM namelist file and pull our relevant lake metrics
%
% Inputs:
%		GLMnmlFile     : filename of GLM namelist file to read from
% Outputs
%		nml : a matlab structure that contains configuration data read from
%		the GLM name list file
%
% Uses:
%
% Written by L. Bruce 18 December 2013
% 


%Open GLM name list file
fid = fopen(GLMnmlFile,'r');

%-------------------------------------------------------------------------%
%Read in all lines -------------------------------------------------------%
%-------------------------------------------------------------------------%

%Loop through all lines in the file discarding comment lines and blanks
ii=1;
while 1
    tline = fgetl(fid);
    if tline < 0 , break, end %End of file reached
    %remove spaces
    tline = regexprep(tline,'\s','');
    %remove tabs
    tline = regexprep(tline,'\t','');
    if isempty(tline) %Blank line
        continue
    end
    if strcmp(tline(1),'!') %Comment line
        continue
    end
    %Remove anything after ! ie comments
    cmnt_i = strfind(tline,'!');
    if ~isempty(cmnt_i)
        tline = tline(1:cmnt_i(1)-1);
    end
    nml.str{ii} = tline;
    ii = ii+1;
end

fclose(fid);
nml.char = char(nml.str);

%Get first character cells
%First 2 characters
c2 = cellstr(nml.char(:,1:2));
c3 = cellstr(nml.char(:,1:3));
c5 = cellstr(nml.char(:,1:5));
c8 = cellstr(nml.char(:,1:8));
c10 = cellstr(nml.char(:,1:10));

%-------------------------------------------------------------------------%
%Read in morphometry -----------------------------------------------------%
%-------------------------------------------------------------------------%

%Number of basin values

B_indx = strcmp(c8,'bsn_vals');
tline = nml.str{B_indx};
bsn_vals = str2num(char(tline(10:end)));

H_indx = find(strcmp(c2,'H='));
A_indx = find(strcmp(c2,'A='));
%Get heights
ii = H_indx;
Hstr = [];
while 1
    tline = nml.str{ii};
    tline = regexprep(tline,'H=','');
    Hstr = [Hstr regexp(tline,',','split')];
    nml.H = str2num(char(Hstr));
    if length(nml.H) >= bsn_vals, break, end %Last area value read in
    ii = ii+1;
end

%Get areas
ii = A_indx;
Astr = [];
while 1
    tline = nml.str{ii};
    tline = regexprep(tline,'A=','');
    Astr = [Astr regexp(tline,',','split')];
    nml.A = str2num(char(Astr));
    if length(nml.A) >= bsn_vals, break, end %Last area value read in
    ii = ii+1;
end

%Basin length and width
bl_indx = strcmp(c8,'bsn_len=');
tline = nml.str{bl_indx};
nml.bsn_len = str2num(char(tline(9:end)));

bw_indx = strcmp(c8,'bsn_wid=');
tline = nml.str{bw_indx};
nml.bsn_wid = str2num(char(tline(9:end)));

%Base elev
be_indx = strcmp(c8,'base_ele');
tline = nml.str{be_indx};
nml.base_elev = str2num(char(tline(11:end)));
%Crest elev
ce_indx = strcmp(c8,'crest_el');
tline = nml.str{ce_indx};
nml.crest_elev = str2num(char(tline(12:end)));

%Get maximum depth for contour plots
nml.max_depth_plots = nml.crest_elev - nml.base_elev;


%-------------------------------------------------------------------------%
%Read in general model set up --------------------------------------------%
%-------------------------------------------------------------------------%

Kw_indx = strcmp(c2,'Kw');
tline = nml.str{Kw_indx};
nml.Kw = str2num(char(tline(4:end)));

Lat_indx = strcmp(c8,'latitude');
tline = nml.str{Lat_indx};
nml.Lat = str2num(char(tline(10:end)));

%-------------------------------------------------------------------------%
%Read in run time parameters----------------------------------------------%
%-------------------------------------------------------------------------%

%Format for time stamp
tf_indx = strcmp(c8,'timefmt=');
tline = nml.str{tf_indx};
tf = str2num(char(tline(9:end)));
%Start date
sd_indx = strcmp(c5,'start');
tline = nml.str{sd_indx};
nml.start_date = datenum(char(tline(8:17)),'yyyy-mm-dd');
%End date
if tf == 2
    fd_indx = strcmp(c5,'stop=');
    tline = nml.str{fd_indx};
    nml.end_date = datenum(char(tline(7:16)),'yyyy-mm-dd');
else
    fd_indx = strcmp(c8,'num_days');
    tline = nml.str{fd_indx};
    nml.end_date = nml.start_date + str2num(tline(10:end));
end

%-------------------------------------------------------------------------%
%Meteorological parameters------------------------------------------------%
%-------------------------------------------------------------------------%

%Met daily flag
md_indx = strcmp(c8,'subdaily');
tline = nml.str{md_indx};
nml.subdaily = tline(10:end);

%Longwave flag
lw_indx = strcmp(c5,'lw_ty');
tline = nml.str{lw_indx};
nml.lw_type = tline(9:end);

%Met file name
mf_indx = strcmp(c8,'meteo_fl');
tline = nml.str{mf_indx};
nml.metfile = tline(11:end-1);

%Wind factor
wf_indx = strcmp(c8,'wind_fac');
tline = nml.str{wf_indx};
nml.wind_factor = str2num(tline(13:end));

%-------------------------------------------------------------------------%
%Output parameters--------------------------------------------------------%
%-------------------------------------------------------------------------%

%Output frequency
ns_indx = strcmp(c5,'nsave');
tline = nml.str{ns_indx};
nml.out_freq = str2double(tline(7:end))/24;

%Lake output file name
of_indx = strcmp(c8,'csv_lake');
tline = nml.str{of_indx};
nml.lakefile = tline(17:end-1);

%-------------------------------------------------------------------------%
%Read inflow parameters---------------------------------------------------%
%-------------------------------------------------------------------------%

%Number of inflows
ni_indx = strcmp(c8,'num_infl');
tline = nml.str{ni_indx};
nml.num_inflows = str2num(tline(13:end));

%Inflow file names
if_indx = find(strcmp(c10,'inflow_fl=')==1);
%Loop through if more than one line
%nml.inflow_files = [];
num_lines = 0;
while 1
    num_lines = num_lines + 1;
    tline = nml.str{if_indx};
    if num_lines == 1
        num_splits = length(findstr(tline(11:end),'.csv'));
        inf_files = regexp(tline(11:end),',','split');
        nml.inflow_files(num_lines,1:num_splits) = inf_files(1:num_splits);
    else
        num_splits = length(findstr(tline,'.csv'));
        inf_files = regexp(tline,',','split');
        nml.inflow_files(num_lines,1:num_splits) = inf_files(1:num_splits);
    end
    nml.inf_num_splits(num_lines) = num_splits;
    if ~strcmp(tline(end),','), break, end  %end of areas list
    if_indx = if_indx + 1;
end
nml.inf_num_lines = num_lines;


%-------------------------------------------------------------------------%
%Read outflow parameters---------------------------------------------------%
%-------------------------------------------------------------------------%

%Number of outflows
no_indx = strcmp(c8,'num_outl');
tline = nml.str{no_indx};
nml.num_outflows = str2num(tline(12:end));

%Outflow file names
of_indx = find(strcmp(c10,'outflow_fl')==1);
%Loop through if more than one line
num_lines = 0;
while 1
    num_lines = num_lines + 1;
    tline = nml.str{of_indx};
    if num_lines == 1
        num_splits = length(findstr(tline(12:end),'.csv'));
        outf_files = regexp(tline(12:end),',','split');
        nml.outflow_files(num_lines,1:num_splits) = outf_files(1:num_splits);
    else
        num_splits = length(findstr(tline,'.csv'));
        outf_files = regexp(tline,',','split');
        nml.outflow_files(num_lines,1:num_splits) = outf_files(1:num_splits);
    end
    nml.outf_num_splits(num_lines) = num_splits;
    if ~strcmp(tline(end),','), break, end  %end of areas list
    of_indx = of_indx + 1;
end
nml.outf_num_lines = num_lines;

