clear; clc;

AU = 1.495978707e11;

% Anfangsbedingungen Komet
r0 = [50*AU; 0];
v0 = [0; 1e3];
y0 = [r0; v0];

% Zeit
tmax = 200 * 365.25 * 24 * 3600;    %200 Jahre
tspan = [0 tmax];

opts = odeset('RelTol',1e-9,'AbsTol',1e-9);

% Simulation simpleG
[tS, yS] = ode45(@kometODE_simple, tspan, y0, opts);

% Simulation fullG
[tF, yF] = ode45(@kometODE_full, tspan, y0, opts);

% Erdbahn
tE = linspace(0, tmax, 2000);
rE = earthPos(tE);

% Plot
figure;
hold on; grid on; axis equal;

plot(yS(:,1)/AU, yS(:,2)/AU, 'b', 'LineWidth', 1.5);
plot(yF(:,1)/AU, yF(:,2)/AU, 'r', 'LineWidth', 1.5);
plot(rE(1,:)/AU, rE(2,:)/AU, 'k--', 'LineWidth', 1.2);

plot(0,0,'yo','MarkerFaceColor','y'); % Sonne

xlabel('x [AE]');
ylabel('y [AE]');
legend('Komet (simpleG)', 'Komet (fullG)', 'Erde', 'Sonne');
title('Bahnkurven von Erde und Komet');