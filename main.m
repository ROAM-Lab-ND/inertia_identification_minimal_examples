% Need to run this file with its enclosing folder as the working directory
fileInfo = dir(matlab.desktop.editor.getActiveFilename);
cd(fileInfo.folder);

clear
init_path

load('CheetahSysID.mat');
model = Cheetah3LegModel();

%% Data Setup
N = length(Y);
n_links = 3; % # of conventional links
n_rotors= 3; % # of rigid-bodies modeled for rotors
n_bodies= n_links+n_rotors;
n_dofs  = length(qd{1});

% Time vector
dt = 1e-3;
t = (0:(length(q)-1))*dt;

% Option to regenerate the regressors.
% The regressors are provided in the dataset, 
% but here's how you can compute them
regenerate_regressors = 0;

if regenerate_regressors
   Y = {};
   for i =1:N
       if mod(i,100) == 0
           fprintf('%d / %d Regressors Computed\n',i,N);
       end
       [Y_i, Yrot_i] = RegressorClassical( model, q{i}, qd{i},qdd{i});
       Y{i,1} = [Y_i Yrot_i];
    end
end

% Setup Friction Regressor
B = repmat({zeros(n_dofs,n_dofs)}, N,1 );
Bc = repmat({zeros(n_dofs,n_dofs)}, N,1 );
for i = 1:N
   B{i}  = diag( qd{i} );
   Bc{i} = diag( sign( qd{i} ) );
end

% Convert cell arrays to matrices
tau_stack = cell2mat(tau);
Y_stack = cell2mat(Y);
tau_mat = cell2mat(tau')';
B_stack = cell2mat(B);
Bc_stack = cell2mat(Bc);

Y_total = [Y_stack B_stack Bc_stack];

% Optional: Select subset for training data
N_train   = N; % Full data set
Y_train   =   Y_stack(1:n_dofs*N_train,:);
B_train   =   B_stack(1:n_dofs*N_train,:);
Bc_train  =  Bc_stack(1:n_dofs*N_train,:);
tau_train = tau_stack(1:n_dofs*N_train,:);

%% Prior
% The parameter order for each link is: 
% m, h_x, h_y, h_z, I_xx, I_yy, I_zz, I_yz, I_xz, I_xy
prior_params = {};
for i = 1:n_links
    prior_params{i}         = inertiaMatToVec( model.I{i}    );
    prior_params{i+n_links} = inertiaMatToVec( model.I_rotor{i} );
end
for i = 1:n_bodies
    J_prior{i} = inertiaVecToPinertia(prior_params{i});
end

%% Conventional System ID with Entropic Regularization

fprintf('=========== Convex Approach (SDP with LMIs) =========\n');

cvx_setup;

% Clear CVX
cvx_begin
cvx_end

% Set solver
cvx_solver mosek % <- Can be changed if you don't have it.

weight_regularization = 1e-2;
use_const_pullback_approx = 0;

cvx_begin
    variable params(10,n_bodies)  % inertial params of links / rotors (units vary)
    variable b(n_dofs)            % viscous friction coefficient (Nm / (rad/s) 
    variable bc(n_dofs)           % coulomb friction coefficient (Nm)
    
    variable J(4,4,n_bodies) semidefinite % pseudo-inertia for each body
    
    expression bregman(n_bodies) % Entropic / bregman divergence for regularization
    expression e_train(3*N_train,1); % Training error (Nm)
    
    % Training error
    e_train = Y_train*params(:) + B_train*b + Bc_train*bc - tau_train;
    
    % Entropic (i.e., bregman) divergence between psuedo-inertias
    for i=1:n_bodies
        if use_const_pullback_approx == 0
            % d( J || J_prior )
            bregman(i) = -log_det( J(:,:,i) ) + log(det(J_prior{i})) ...
                            + trace(J_prior{i} \ J(:,:,i) ) - 4;
        else
            % constant pullback approximation
            M{i} = pullbackMetric(prior_params{i});
            bregman(i) = 1/2*(params(:,i) - prior_params{i})'*M{i}*(params(:,i) - prior_params{i});
        end
    end
    
    % Objective = Squared 2-norm of residual + Regularization
    minimize( 1/2*sum_square_abs(e_train)/length(e_train) ...
                    + weight_regularization*sum(bregman) )
    
    % Only constraints are to set the pseudo-inertias based on params
    subject to:
        for i=1:n_bodies
            J(:,:,i) == inertiaVecToPinertia( params(:,i) );
        end
cvx_end

tau_predict_entropic = reshape(Y_stack*params(:) + B_stack*b + Bc_stack*bc,n_dofs,N)';
plotTorquePredictions(1,'Convex, Entropic Regularized',t,tau_mat, tau_predict_entropic);


%% Unconstrained System ID with Entropic Regularization

fprintf('=========== Log Cholesky Approach (unconstrained) =========\n');

options = optimoptions('fminunc');
options.Algorithm = 'trust-region';
options.HessianFcn = 'objective';
options.SpecifyObjectiveGradient = true;
options.StepTolerance = 1.0000e-6;
options.MaxFunctionEvaluations = 132000;
options.MaxIterations = 2000;
options.Display = 'iter';

data.Y   = [Y_train B_train Bc_train];
data.tau = tau_train;
data.prior = cell2mat(prior_params);
data.gamma = weight_regularization;

x_init = zeros(66,1);    

% Main optimization step
cholParams = fminunc( @(x) objective(x, data), x_init,options );

% System Param Map converts from cholesky Parameters to conventional params
piParams = systemParamMap( cholParams, data.prior );


tau_predic_logchol = reshape(Y_total*piParams - tau_stack, n_dofs,N)';

plotTorquePredictions(2,'Log Cholesky, Entropic Regularized',t,tau_mat, tau_predict_entropic);


%% Kinematics Plotting
q_mat = cell2mat(q')';
figure(3)
clf
subplot(311);
plot(t,q_mat(:,1))
ylabel('Ab/Ad angle (rad)');
title('Leg Kinematics');
subplot(312);
plot(t,q_mat(:,2))
ylabel('Hip angle (rad)');
subplot(313);
plot(t,q_mat(:,3))
ylabel('Knee angle (rad)');
xlabel('Time (s)');

%% Helpers
function RMS= rms(varargin)
    %
    % Written by Phillip M. Feldman  March 31, 2006
    %
    % rms computes the root-mean-square (RMS) of values supplied as a
    % vector, matrix, or list of discrete values (scalars).  If the input is
    % a matrix, rms returns a row vector containing the RMS of each column.
    % David Feldman proposed the following simpler function definition:
    %
    %    RMS = sqrt(mean([varargin{:}].^2))
    %
    % With this definition, the function accepts ([1,2],[3,4]) as input,
    % producing 2.7386 (this is the same result that one would get with
    % input of (1,2,3,4).  I'm not sure how the function should behave for
    % input of ([1,2],[3,4]).  Probably it should produce the vector
    % [rms(1,3) rms(2,4)].  For the moment, however, my code simply produces
    % an error message when the input is a list that contains one or more
    % non-scalars.
    if (nargin == 0)
       error('Missing input.');
    end
    % Section 1: Restructure input to create x vector.
    if (nargin == 1)
       x= varargin{1};
    else
       for i= 1 : size(varargin,2)
          if (prod(size(varargin{i})) ~= 1)
             error(['When input is provided as a list, ' ...
                    'list elements must be scalar.']);
          end
          x(i)= varargin{i};
       end
    end
    % Section 2: Compute RMS value of x.
    RMS= sqrt (mean (x .^2) );
end

function plotTorquePredictions(figNum, name, t, tau_actual , tau_predict)
    tau_rms = rms(tau_actual - tau_predict);
    figure(figNum)
    clf
    subplot(311)
    plot(t,tau_actual(:,1)); hold on;
    plot(t,tau_predict(:,1),'r','LineWidth',1.5 )
    ylabel('Ab/Ad Torque (Nm)');
    rms_tag = sprintf('(%s) [RMS=%.2f, %.2f, %.2f (Nm)]',name, tau_rms(1), tau_rms(2), tau_rms(3));
    title(['Torque Predictions ' rms_tag]);
    legend('Measured','Predicted')

    subplot(312)
    plot(t,tau_actual(:,2)); hold on;
    plot(t,tau_predict(:,2),'r','LineWidth',1.5)
    ylabel('Hip Torque (Nm)');

    subplot(313)
    plot(t,tau_actual(:,3)); hold on;
    plot(t,tau_predict(:,3),'r','LineWidth',1.5)
    ylabel('Knee Torque (Nm)');
    xlabel('Time (s)');
end

% All of this from here on out is basically just to get analytical
% jacobians / hessians of the objective.

function [f, J, H] = objective(x, data)
    if nargout == 1
        f = regularizedLeastSquares( systemParamMap(x,data.prior), data);
    elseif nargout == 1
        % We have a map g: Log Cholesky -> standard parameters
        % We also have a objective function f : standard parameters -> Real
        % It is easy to get the jacobian and hessian of f and g. 
        % The function used below gets the jacobian and hessian of their
        % composition.
        [f, J] = composedJacobianAndHessian( @(x) regularizedLeastSquares(x, data) ...
                                           , @(x) systemParamMap(x, data.prior) , x);
    else
        [f, J, H] = composedJacobianAndHessian(  @(x) regularizedLeastSquares(x, data) ...
                                           , @(x) systemParamMap(x, data.prior) , x);
    end
end
