% Discrete time solver
% uses griddedInterpolant and fmincon to speed up the DP in comparison to
% the algorithm used in the CEMA2019 submission
% also uses parfor to speed it up dramatically
clear variables
T = 1; % length of period
pen = 300; % penalty per unit of non compliance
S_min = 0;
S_max = pen;
req = 500; % requirement for compliance
b_max = 2*req;
h = 1*req; % constant for now, can make this deterministic in time without anything changing
mu_f = 0; % drift of F_t (assumed ABM here, but more realistically a jump process)
sigma_f = 10; % volatility of F_t (assumed ABM here, but more realistically a jump process)
zeta = 0.6;
gamma = 0.6;
time_steps = 50;
dt = T/time_steps;
dS = sqrt(3 * dt) * sigma_f;
grid_points_s = ceil((S_max - S_min)/dS);
grid_points_b = 401; % have to choose this carefully so that we get enough grid points for the interpolation to not be problematic
T_grid = linspace(0, 1, time_steps+1);
psi = 0.01;
eta = 0.02;
nsim = 2000;

nu = 10; % standard deviation of stochastic generation noise


options = optimset('Display', 'off', 'MaxFunEvals', 100); % options for minimization
s = rng; % set seed
% initialization of parameters
S_grid = linspace(S_min, S_max, grid_points_s);
b_grid = linspace(0, 2*req, grid_points_b);
[bb, ss] = meshgrid(b_grid, S_grid);
[X, Y] = ndgrid(S_grid, b_grid);

pars = [time_steps pen req h mu_f sigma_f zeta gamma T psi eta];

V = zeros(time_steps, grid_points_s, grid_points_b);
gen_opt = zeros(time_steps, grid_points_s, grid_points_b);
trade_opt = zeros(time_steps, grid_points_s, grid_points_b);
% iterating through timesteps and minimizing the cost at each step
rng(s);
tic
%parpool(25)
% constraints for fmincon
A = [-1, 0];
bb = 0;
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];

for t = time_steps:-1:1
    s_noise = normrnd(0, sqrt(dt), [1, nsim]);
    e_noise = normrnd(0, nu*sqrt(dt), [1 nsim]);
    if t < time_steps
        F = griddedInterpolant(X, Y, squeeze(V(t+1, : ,:)));
    else
        F = 0;
    end
    parfor s = 1:grid_points_s
        for b = 1:grid_points_b
            if b_grid(b) + h * dt * (time_steps-t+1) >= req % guaranteed to comply with baseline generation
                b0 = h;
                t0 = -S_grid(s)/gamma; 
            elseif b_grid(b) + (h + pen / zeta) * dt * (time_steps-t+1) <= req
                % not going to comply even with extreme production
                b0 = h + pen / zeta;
                t0 = (pen - S_grid(s))/gamma;
                if t0 == 0 % helps with getting the right minimum
                    t0 = -1;
                end
            else % should roughly correspond to the 'interesting' regime
                b0 = req;
                t0 = -1;
            end
            x0 = [b0 t0];
            f = @(ctrl)runningCostAlt(ctrl, t, S_grid(s), b_grid(b), pars, F, s_noise, e_noise);
            [x, fval] = fmincon(f, x0, A, bb, Aeq, beq, lb, ub, nonlcon, options);
            V(t, s, b) = fval;
            gen_opt(t,s,b) = x(1);
            trade_opt(t,s,b) = x(2);
        end
        fprintf([ num2str(s) ' ']);
        if mod(s,10) ==0
            fprintf('\n');
        end 
    end
    fprintf([ num2str(t) '\n']);
end
toc
save('sg_pi_etamax_psimax.mat')

% runningCostAlt uses the faster method of interpolation - results in
% significant time savings

% x0 = [500 0]
% f = @(ctrl)runningCostAlt(ctrl, 49, 150, 500, pars, F, s_noise, e_noise);
% [x, fval] = fmincon(f, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
function y = runningCostAlt(ctrl, t, S, b, pars, F, s_noise, e_noise)
% actions, timestep, stock price, banking level, model parameters, old
% values, grid_points
time_steps = pars(1);
pen = pars(2);
req = pars(3);
h = pars(4);
mu_f = pars(5);
sigma_f = pars(6);
zeta = pars(7);
gamma = pars(8);
T = pars(9);
psi = pars(10);
eta = pars(11);

dt = T/time_steps;
gen = ctrl(1);
trade= ctrl(2);
    if t == time_steps
        z = 1 / 2 * zeta * max(gen - h, 0)^2 *dt + trade.*S*dt +...
            1 / 2 * gamma * trade^2 *dt + pen*max(req - max(0, gen*dt + e_noise) - trade*dt - b, 0);
        y = mean(z);
    else
        new_b = min(max(0, (b + max(0, gen*dt + e_noise) + trade*dt)), 2*req);
        new_S = max(0, min(pen, S + mu_f * dt - psi * max(0, gen*dt + e_noise) + eta * trade * dt + sigma_f * s_noise));
        running_value = 1 / 2 * zeta * max(gen - h, 0)^2 *dt + trade*S*dt +...
            1 / 2 * gamma * trade^2 *dt;
        Vq = F(new_S', new_b');
        future_value = mean(Vq); % need to interpolate this from vals(t+1,:,:)
        y = running_value + future_value;        
    end
end

