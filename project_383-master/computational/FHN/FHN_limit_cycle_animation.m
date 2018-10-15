% vector field visualization

% creates animation for vector field of FHN under standard parameters at
% different values of the input current. Shows how limit cycle
% emerges/disappears as the input current changes.

close all; clc;
FHN_eqns = @(t,y, phi, a, b, I_app) [
    y(1) - (y(1).^3)./ 3 - y(2) + I_app; 
    phi*(y(1) + a - b*y(2))
    ];

% standard parameters
a = 0.7;
b = 0.8;
phi = 0.08;

% initial conditions for trajectories in V-W space
y0 = [0 0; 0.2 0.2; 0.5 0.5; -0.1 0.5; -0.5 -0.5; -0.5 0.5; 0.5 -0.5]';
tspan = [0,100];
W_vals = -20:0.1:20;
V_vals = -20:0.1:20;

% returns W_vals associated with given V values 
V_nullc = @(V, I_app) -V.^3/3 + V + I_app;

%returns V_values associated with given W values
W_nullc = @(W) -a + b .* W;

figure(1)

for I_app = 0:0.01:2
    V_dot = @(X, Y) X - (X.^3)./ 3 - Y + I_app;
    W_dot = @(X, Y) phi*(X + a - b*Y);
    [X, Y] = meshgrid(linspace(-2.5,2.5, 20),linspace(-1,2, 20) );
    U = V_dot(X,Y);
    V = W_dot(X,Y);
    
    dims = size(y0);
    for i = 1:dims(2)
        [t_vals, y_vals] = ode45(@(t,y) FHN_eqns(t,y, phi, a,b, I_app), tspan, y0(:,i));
        plot(y_vals(:,1), y_vals(:,2), 'k'); % solution in V-W space
        hold on;
    end
    
    quiver(X,Y,U,V,'r', 'autoscalefactor', 2);
    hold on;
    
    plot(W_nullc(W_vals), W_vals, 'm'); % W nullcline
    hold on;
    plot(V_vals, V_nullc(V_vals, I_app), 'b'); % V nullcline
    hold on;
    xlim([-2.5, 2.5]);
    ylim([-1, 2]);
    title(strcat('I_{app} = ', num2str(I_app)));
    pause(0.1);
    clf;
end
    
    