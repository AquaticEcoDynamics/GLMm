function fitparams = runGLMModelFit(theta,fld_data,params,out_name,conf)

% function nmae = runGLMModelFit(theta,fld_data,params,run_info,out_name)
%
% Inputs:
%       fld_data   : MATLAB data structure containing observed data
%		params     : list of parameters for current run
%       out_name   : text used to describe current run
%       conf       : configuration information required to run GLM
%
%       fit parameters fld = name given to field variable, sim =
%       corresponding name given to GLM simulated data
%
% Outputs:
%       fitparms = MATLAB data strycture containing all model fit
%       parameters comparing simulated to field
%
% Uses:
%      readGLMnetcdf.m
%
% Written by L. Bruce 15 August 2013
% First takes parameter values and creates a glm2.nml
% Then runs the model off remote host to create a simulation file output.nc
% Takes GLM simulated output and compares against field data for
% variable "varname".
% Calculates the following measures of best fit for varname:
%    Sample size for observed data (N)
%    Mean of observed data (mean)
%    Variance for observed data (var)
%    Normalised mean absolute error (NMAE)
%    Root square mean error (RSME)
%    Normalised mean error (NME)
%    Percent relative error (PRE)
%    Normalised root square mean error (NRSME)
%    Nash-Sutcliffe Model Efficiency (NSE)
%    Correlation coefficient (r), offsett and slope
%  Saves results to .csv file


%Path names
paths = conf.paths;

%Name of output.nc file to save sim file to
GLMfile = ['output_',out_name,'.nc'];

%-------------------------------------------------------------------------%
%First write new glm2.nml based on current parameters ---------------------%
%-------------------------------------------------------------------------%
newGLMnml(theta,conf,params)

%-------------------------------------------------------------------------%
%                        REMOTE CONNECTION                                %
%-------------------------------------------------------------------------%

if conf.config.remote_flag==1
   %-------------------------------------------------------------------------%
   %Second copy the new glm2.nml,run model from remote host, return output----%
   %-------------------------------------------------------------------------%

   run_info = conf.run_info;
      
   %Open up connection with remote host
   ssh2_conn = ssh2_config(run_info.host_name,run_info.usr_name,run_info.password);
   
   %Copy across new glm2.nml file with new parameter set
   ssh2_conn = scp_put(ssh2_conn, 'glm2_new.nml',run_info.remote_dir,'nml','glm2.nml');

   %Run the model
   ssh2_conn = ssh2_command(ssh2_conn, run_info.run_glm);
   %Copy output file to working directory
   ssh2_conn = scp_get(ssh2_conn, run_info.output_file);

   %Move output file to the Output directory then analyse model fit for new
   %parameters
   movefile('output.nc' ,['Output/',GLMfile]);
   
   %close connection when done
   ssh2_conn = ssh2_close(ssh2_conn);


%-------------------------------------------------------------------------%
%                        WINDOWS RUN                                      %
%-------------------------------------------------------------------------%
elseif conf.config.remote_flag == 0
   %-------------------------------------------------------------------------%
   %Second copy the new glm2.nml to sim folder, run model, return output------%
   %-------------------------------------------------------------------------%

   %Copy across new glm2.nml file with new parameter set
   movefile('nml/glm2_new.nml','sim/glm2.nml');

   %Run the model
   oldpth=pwd; 
   cd('sim')
   run_glm = ['!',conf.paths.run_glm];
   eval(run_glm') 
   cd(oldpth) 
   %Move output file to the Output directory then analyse model fit for new
   %parameters
   movefile('sim/output.nc' ,['Output/',GLMfile]);
   
%-------------------------------------------------------------------------%
%                        MAC RUN                                      %
%-------------------------------------------------------------------------%
elseif conf.config.remote_flag == 2
   %-------------------------------------------------------------------------%
   %Second copy the new glm2.nml to sim folder, run model, return output------%
   %-------------------------------------------------------------------------%

   %Copy across new glm2.nml file with new parameter set
   movefile('nml/glm2_new.nml','sim/glm2.nml');

   %Run the model
   oldpth=pwd; 
   cd('sim') 
   eval('!sh glm.sh') 
   cd(oldpth) 
   %Move output file to the Output directory then analyse model fit for new
   %parameters
   movefile('sim/output.nc' ,['Output/',GLMfile]);

end %Run remote or off Windows

%-------------------------------------------------------------------------%
%Third return measure of model fit to the Optimisation routine -----------%
%-------------------------------------------------------------------------%
%Determine nmae for the new output file
%At this stage just calibrating against temperature
fitparams = calcGLMModelFit(fld_data,['Output/',GLMfile],conf.config.varname,conf.config.spin_up);

%-------------------------------------------------------------------------%
%Save results ------------------------------------------------------------%
%-------------------------------------------------------------------------%

save([conf.paths.working_dir,'Results/GLM_fitparams_',out_name,'.mat'],'fitparams')
movefile('GLM_fitparams.csv',[conf.paths.working_dir,'Results/GLM_fitparams_',out_name,'.csv']);