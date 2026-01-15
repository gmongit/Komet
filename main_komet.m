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
tE = linspace(0, tmax, 20000);
rE = earthPos(tE);

% Minimaler Abstand Erde–Komet

% Erdpositionen zu den jeweiligen Zeitpunkten
rE_simple = earthPos(tS');   
rE_full   = earthPos(tF');  

% Kometenpositionen
rK_simple = yS(:,1:2)';      
rK_full   = yF(:,1:2)';

% Abstände
d_simple = vecnorm(rK_simple - rE_simple, 2, 1);
d_full   = vecnorm(rK_full   - rE_full,   2, 1);

% Minima
[dmin_simple, idxS] = min(d_simple);
[dmin_full,   idxF] = min(d_full);

% Zeitpunkte
tmin_simple = tS(idxS);
tmin_full   = tF(idxF);

% Differenz
delta_d = abs(dmin_simple - dmin_full);





% Ausgabe
fprintf('Minimaler Abstand (simpleG): %.3e m (%.3f AE)\n', ...
        dmin_simple, dmin_simple/AU);

fprintf('Minimaler Abstand (fullG)  : %.3e m (%.3f AE)\n', ...
        dmin_full, dmin_full/AU);

fprintf('Differenz: %.3e m (%.3f AE)\n', ...
        delta_d, delta_d/AU);





% Plot
figure;
hold on; grid on; axis equal;

plot(yS(:,1)/AU, yS(:,2)/AU, 'b', 'LineWidth', 1.5);
plot(yF(:,1)/AU, yF(:,2)/AU, 'r', 'LineWidth', 1.5);
plot(rE(1,:)/AU, rE(2,:)/AU, 'k--', 'LineWidth', 1.2);

plot(0,0,'yo','MarkerFaceColor','y'); % Sonne

plot(rK_simple(1,idxS)/AU, rK_simple(2,idxS)/AU, 'bo', 'MarkerFaceColor','b');
plot(rK_full(1,idxF)/AU,   rK_full(2,idxF)/AU,   'ro', 'MarkerFaceColor','r');

xlabel('x [AE]');
ylabel('y [AE]');
legend('Komet (simpleG)', 'Komet (fullG)', 'Erde', 'Sonne');
title('Bahnkurven von Erde und Komet');