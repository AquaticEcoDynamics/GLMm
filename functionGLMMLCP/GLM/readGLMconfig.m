function data = readGLMconfig(filename)
% function data = readGLMconfig(filename)
% Inputs:
%      filename: Name of GLM configuration name list file
%
% Outputs:
%      data: MATLAB structure containing configuration setting for GLM MCMC
%      parameter estimation run
%
%
% Written by b. Busch 1 November 2013
%
% Reads configuration information from name list including list of
% parameters used in MCMC run.

fid = fopen(filename,'rt');

data = [];

EOF=0;
while ~EOF
    tline = fgetl(fid);
    if (tline == -1), break, end %end of nml file
    if isempty(tline)
        continue
    end
    
    if strcmp(tline(1),'&')
        
        isinternal = 1;
        varname = regexprep(tline,'&','');
        while isinternal
            
            tline = fgetl(fid);
            
            if (tline == -1), break, end %end of nml file
            if isempty(tline)
                continue
            else
                if ~strcmp(tline(1),'&')
                    if ~strcmp(tline(1),'!')  %Is not a comment line
                        %tline_cmp = regexprep(tline,'\s','');
                        eval(['data.(varname).',tline]);
                    end
                else
                    varname = regexprep(tline,'&','');
                    
                end
            end
        end
    end
end
fclose(fid);