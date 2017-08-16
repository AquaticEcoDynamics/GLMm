function newGLMnml(theta,conf,params)
%
%This function takes the initial glm2.nml from nml and creates a new
%glm2.nml with parameters taken from param updated to theta determined by
%the MCMC chain
%

%-------------------------------------------------------------------------%
%First take initial glm2.nml file from nml and insert parameters----%
%-------------------------------------------------------------------------%
%
%Open GLM name list file
fid = fopen([conf.paths.working_dir,'nml/glm2_init.nml'],'r');


%Loop through line by line until reach parameter then insert parameter
%value

ii=1;
while 1
    isparam = 0;
    tline = fgetl(fid);
    if (tline == -1), break, end %end of nml file
    if isempty(tline)
        continue
    end
    if ~strcmp(tline(1),'!')  %Is not a comment line
         tline_cmp = regexprep(tline,'\s','');
         equal_i = strfind(tline_cmp,'=');
         if ~isempty(equal_i) %Is a parameter line
            %Check to see if parameter is included in list of parameters
             for par_i = 1:length(params)
                 if strcmp(params{par_i}{1},tline_cmp(1:equal_i-1));
                     nml.txt{ii} = ['   ',tline_cmp(1:equal_i-1),' = ',num2str(theta(par_i))];
                     isparam = 1;
                 end
             end
         end
    end
    %if parameter line not changed then stay as same
    if ~isparam
        nml.txt{ii} = tline;
    end
    ii = ii+1;
end
fclose(fid);

%-------------------------------------------------------------------------%
%Write new glm2.nml GLM name list file ------------------------------------%
%-------------------------------------------------------------------------%
%Open new GLM name list file
fid = fopen([conf.paths.working_dir,'nml/glm2_new.nml'],'w');
for ii = 1:length(nml.txt)
    fprintf(fid,'%s\n',nml.txt{ii});
end
fclose(fid);
