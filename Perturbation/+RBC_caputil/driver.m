%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'RBC_caputil';
M_.dynare_version = '5.0';
oo_.dynare_version = '5.0';
options_.dynare_version = '5.0';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'epsiz'};
M_.exo_names_tex(1) = {'epsiz'};
M_.exo_names_long(1) = {'epsiz'};
M_.exo_names(2) = {'epsit'};
M_.exo_names_tex(2) = {'epsit'};
M_.exo_names_long(2) = {'epsit'};
M_.endo_names = cell(6,1);
M_.endo_names_tex = cell(6,1);
M_.endo_names_long = cell(6,1);
M_.endo_names(1) = {'c'};
M_.endo_names_tex(1) = {'c'};
M_.endo_names_long(1) = {'c'};
M_.endo_names(2) = {'n'};
M_.endo_names_tex(2) = {'n'};
M_.endo_names_long(2) = {'n'};
M_.endo_names(3) = {'u'};
M_.endo_names_tex(3) = {'u'};
M_.endo_names_long(3) = {'u'};
M_.endo_names(4) = {'k'};
M_.endo_names_tex(4) = {'k'};
M_.endo_names_long(4) = {'k'};
M_.endo_names(5) = {'z'};
M_.endo_names_tex(5) = {'z'};
M_.endo_names_long(5) = {'z'};
M_.endo_names(6) = {'theta'};
M_.endo_names_tex(6) = {'theta'};
M_.endo_names_long(6) = {'theta'};
M_.endo_partitions = struct();
M_.param_names = cell(13,1);
M_.param_names_tex = cell(13,1);
M_.param_names_long = cell(13,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'alpha'};
M_.param_names_tex(2) = {'alpha'};
M_.param_names_long(2) = {'alpha'};
M_.param_names(3) = {'dss'};
M_.param_names_tex(3) = {'dss'};
M_.param_names_long(3) = {'dss'};
M_.param_names(4) = {'B'};
M_.param_names_tex(4) = {'B'};
M_.param_names_long(4) = {'B'};
M_.param_names(5) = {'gx'};
M_.param_names_tex(5) = {'gx'};
M_.param_names_long(5) = {'gx'};
M_.param_names(6) = {'gamma'};
M_.param_names_tex(6) = {'gamma'};
M_.param_names_long(6) = {'gamma'};
M_.param_names(7) = {'rhot'};
M_.param_names_tex(7) = {'rhot'};
M_.param_names_long(7) = {'rhot'};
M_.param_names(8) = {'rhoz'};
M_.param_names_tex(8) = {'rhoz'};
M_.param_names_long(8) = {'rhoz'};
M_.param_names(9) = {'sigmat'};
M_.param_names_tex(9) = {'sigmat'};
M_.param_names_long(9) = {'sigmat'};
M_.param_names(10) = {'sigmaz'};
M_.param_names_tex(10) = {'sigmaz'};
M_.param_names_long(10) = {'sigmaz'};
M_.param_names(11) = {'mu'};
M_.param_names_tex(11) = {'mu'};
M_.param_names_long(11) = {'mu'};
M_.param_names(12) = {'phi1'};
M_.param_names_tex(12) = {'phi1'};
M_.param_names_long(12) = {'phi1'};
M_.param_names(13) = {'phi2'};
M_.param_names_tex(13) = {'phi2'};
M_.param_names_long(13) = {'phi2'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 6;
M_.param_nbr = 13;
M_.orig_endo_nbr = 6;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [1 2 3 4];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 6;
M_.eq_nbr = 6;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 4 10;
 0 5 11;
 0 6 12;
 1 7 0;
 2 8 13;
 3 9 14;]';
M_.nstatic = 0;
M_.nfwrd   = 3;
M_.npred   = 1;
M_.nboth   = 2;
M_.nsfwrd   = 5;
M_.nspred   = 3;
M_.ndynamic   = 6;
M_.dynamic_tmp_nbr = [18; 19; 4; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'Intertemporal Euler Equation' ;
  2 , 'name' , 'Intratemporal Euler Equation' ;
  3 , 'name' , 'Capital utilization' ;
  4 , 'name' , 'Budget constraint' ;
  5 , 'name' , 'Production shock' ;
  6 , 'name' , 'Utilization shock' ;
};
M_.mapping.c.eqidx = [1 2 4 ];
M_.mapping.n.eqidx = [1 2 3 4 ];
M_.mapping.u.eqidx = [1 2 3 4 ];
M_.mapping.k.eqidx = [1 2 3 4 ];
M_.mapping.z.eqidx = [1 3 4 6 ];
M_.mapping.theta.eqidx = [1 2 3 4 5 ];
M_.mapping.epsiz.eqidx = [6 ];
M_.mapping.epsit.eqidx = [5 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [4 5 6 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(6, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(13, 1);
M_.endo_trends = struct('deflator', cell(6, 1), 'log_deflator', cell(6, 1), 'growth_factor', cell(6, 1), 'log_growth_factor', cell(6, 1));
M_.NNZDerivatives = [30; 94; -1; ];
M_.static_tmp_nbr = [11; 5; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 0.988173;
beta = M_.params(1);
M_.params(2) = 0.333333;
alpha = M_.params(2);
M_.params(11) = 5;
mu = M_.params(11);
M_.params(3) = 0.01843;
dss = M_.params(3);
M_.params(5) = 1.0029;
gx = M_.params(5);
M_.params(12) = 0.0333333;
phi1 = M_.params(12);
M_.params(13) = 100;
phi2 = M_.params(13);
M_.params(6) = 1;
gamma = M_.params(6);
M_.params(7) = 0.9;
rhot = M_.params(7);
M_.params(8) = 0.9;
rhoz = M_.params(8);
M_.params(9) = 0.0096;
sigmat = M_.params(9);
M_.params(10) = 0.0072;
sigmaz = M_.params(10);
M_.params(4) = 0.334784;
B = M_.params(4);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(1) = 0.829255;
oo_.steady_state(4) = 10.5409;
oo_.steady_state(2) = 0.333333;
oo_.steady_state(3) = 1;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(10))^2;
M_.Sigma_e(2, 2) = (M_.params(9))^2;
options_.irf = 100;
options_.order = 2;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'RBC_caputil_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_caputil_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_caputil_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_caputil_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_caputil_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_caputil_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_caputil_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
