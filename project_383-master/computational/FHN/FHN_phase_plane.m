close all; clc; clear all;

format short
% y = [V, w]
FHN_eqns = @(t,y, phi, a, b, I_app) [
    y(1) - (y(1).^3)./ 3 - y(2) + I_app; 
    phi*(y(1) +a - b*y(2))
    ];

% standard conditions for oscillatory behavior
a_std = 0.7;
b_std = 0.8;
phi_std = 0.08;
I_app = 0.5;

% each pair is a different initial condition
y0 = [0 0; 0.2 0.2; 0.5 0.5; -0.1 0.5; -0.5 -0.5; -0.5 0.5; 0.5 -0.5]';
figure(1);
[fixed_points, stability] = plot_phase_plane(FHN_eqns, phi_std, a_std, b_std, I_app, y0)



% stable fixed point; stable spiral
I_app = 0;

figure(2)
[fixed_points, stability] = plot_phase_plane(FHN_eqns, phi_std, a_std, b_std, I_app, y0)


% Bistability: No oscillations (3 real fixed points)

b_std = 2;
I_app = 0.5;
figure(3)
[fixed_points, stability] = plot_phase_plane(FHN_eqns, phi_std, a_std, b_std, I_app, y0)


%% analyze bifurcation 
close all;
FHN_eqns = @(t,y, phi, a, b, I_app) [
    y(1) - (y(1).^3)./ 3 - y(2) + I_app; 
    phi*(y(1) +a - b*y(2))
    ];

a_std = 0.7;
b_std = 0.8;
phi_std = 0.08;
I_app = 0.5;
y0 = [0 0; 0.2 0.2; 0.5 0.5; -0.1 0.5; -0.5 -0.5; -0.5 0.5; 0.5 -0.5]';
figure(4)
[fixed_points, stability] = plot_phase_plane(FHN_eqns, phi_std, a_std, b_std, I_app, y0)


b_std = 1.3;
figure(5)
[fixed_points, stability] = plot_phase_plane(FHN_eqns, phi_std, a_std, b_std, I_app, y0)

b_std = 1.5;
figure(6)
[fixed_points, stability] = plot_phase_plane(FHN_eqns, phi_std, a_std, b_std, I_app, y0)
%% plotting phase plane

% standard parameters for oscillatory behavior
% to source [5]
function [fixed_points, stability] = plot_phase_plane(FHN_eqns, phi, a, b, I_app, y0)

FHN_std = @(t,y) FHN_eqns(t,y, phi, a, b, I_app);

% returns W_vals associated with given V values 
V_nullc = @(V, I_app) -V.^3/3 + V + I_app;

%returns V_values associated with given W values
W_nullc = @(W) -a + b .* W;

W_vals = -20:0.1:20;
V_vals = -20:0.1:20;

tspan = [0,100];



V_dot = @(X, Y) X - (X.^3)./ 3 - Y + I_app;
W_dot = @(X, Y) phi*(X + a - b*Y);
[X, Y] = meshgrid(-20:5:20, -20:5:20);
U = V_dot(X,Y);
V = W_dot(X,Y);

%plots nullclines and solution

subplot(121);
%quiver(X, Y, V, U, 'k'); %doesnt work for some reason
%hold on;
plot(W_nullc(W_vals), W_vals, 'm'); % W nullcline
hold on;
plot(V_vals, V_nullc(V_vals, I_app), 'b'); % V nullcline
hold on;

dims = size(y0);
for i = 1:dims(2)
    [t_vals, y_vals] = ode45(FHN_std, tspan, y0(:,i));
    plot(y_vals(:,1), y_vals(:,2), 'k'); % solution in V-W space
    hold on;
end

plot(zeros(1, 21), -10:10, 'k'); % central axis
hold on;
plot(-10:10, zeros(1,21), 'k'); % central axis
hold on;
text(-2,-0.8, strcat("I_{app} = ", num2str(I_app), ", a = ", num2str(a), ", b = ", num2str(b), ", \phi = ", num2str(phi)), 'fontsize', 20);




xlabel("V (volts)");
ylabel("W (mAmps)");
xlim([-2.5, 2.5]); %-2.5,2.5
ylim([-1,2]) %-1,2
lg = legend({"W nullcline $(\dot{W} = 0)$", "V nullcline $( \dot{V} = 0)$", "trajectory of FHN model"});
lg.Interpreter = 'latex';
lg.FontSize = 20;
lg.Location = 'best';
set(gca, 'fontsize', 20);
title("Phase Plane");

subplot(122)
plot(t_vals, y_vals(:,1), 'r');
hold on;
plot(-10:100, 0*(-10:100) - 0.8, 'k'); % central axis
hold on;
ylabel("Membrane Potential (V)");
xlabel("time (ms)");
xlim(tspan)
ylim([-5, 5]);
set(gca, 'fontsize', 20);



% Finds Fixed points under given parameters
syms V
eqn = b/3 * V.^3 + (1-b)*V + (a - I_app) == 0;
V_star = vpa(solve(eqn, V));
W_star = (V_star + a) / b;
fixed_points = [V_star, W_star];

% Determines stability of fixed points under given parameters
stability = [];
for i = 1:length(fixed_points)
    fp = fixed_points(i, :);
    if isreal(fp)
        J_f = [1-fp(1).^2, -1; phi, -b*phi];
        format short
        eigenvals = eig(J_f);
        Jf_eigenvalues = round(double(eigenvals), 5)
        if all(real(eigenvals) < 0) %stable only if all eigenvalues real componenet < 0;
            %stable
            if ~all(isreal(eigenvals))
                stability = [stability; "stable spiral"];
            else 
               stability = [stability; "stable node"]; 
            end
        elseif all(real(eigenvals) > 0)
            if ~all(isreal(eigenvals))
                stability = [stability; "unstable spiral"];
            else
                stability = [stability; "unstable node"];
            end
        else
            stability = [stability; "saddle"];
        end
    else
       stability = [stability; "N/A"]; % not stable/cannot access imaginary roots 
    end
end

fixed_points = double(fixed_points);

end