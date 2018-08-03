close all;
% gating function parameters taken from Foundations of Comp Neuroscience pg
% 32: same values Hodgkin and Huxley got experimentally from squid
% constant values on pg 32 of this textbook

% define constants

% max conductances
g_na = 120;
g_k = 36;
g_l = 0.3;

% Equilibrium Potentials
E_k = -77;
E_na = 50;
E_l = -54.4;

c_m = 1; % cannot find value for this
phi = 1; % scales for temperature
I_ex = @(t) 10 * (t >= 50) .* (t <= 100); % input current; constant in this case, but we can change that if we wish

% x = [V n m h]
diff_eqns = @(t,x) [
1/c_m .*(I_ex(t) - g_na * x(3)^3 * x(4) * (x(1) - E_na) - g_k*x(2)^4 * (x(1) - E_k) - g_l * (x(1) - E_l));
phi * (alpha_n(x(1)) * (1-x(2)) - beta_n(x(1)) * x(2));
phi * (alpha_m(x(1)) * (1-x(3)) - beta_m(x(1)) * x(3));
phi * (alpha_h(x(1)) * (1-x(4)) - beta_h(x(1)) * x(4))];

x0 = [-70 0.3 0.1 0.6]; % initial conditions; approximated off graph pg 33 FCN
tspan = [0 200];

[t_vals, y_vals] = ode45(diff_eqns, tspan, x0);

subplot(311)
plot(t_vals, y_vals(:, 1), 'r-'), hold on;
plot(t_vals, ones(1, length(t_vals)) * -55, 'k--');
ylabel("Membrane Potential (mV)");
xlabel("time (ms)");
legend({"V(t)", "threshold"});
title("Action Potential");
xlim(tspan)

subplot(312)
I_vals = I_ex(t_vals);
plot(t_vals, I_vals, 'm');
ylabel("I_{ext}(t)")
title("Applied Current I(t)");
ylim([0 1.25*max(I_vals)]);

subplot(313);
plot(t_vals, y_vals(:, 2)), hold on;
plot(t_vals, y_vals(:, 3)), hold on;
plot(t_vals, y_vals(:, 4));
title("Gating Variables")
legend({"n", "m", "h"});
xlim(tspan);

gcf