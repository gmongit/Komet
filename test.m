function fit_karthesisch_full_realistic()
% Kartesische Beobachtungen: Best-fit Anfangsbedingungen + minimaler Abstand zur Erde
% FULL-Modell: Sonne + Erde
% Zusatz: "realistische Bahn" per Penalty => Perihel zur Sonne darf nicht unter qMinAU fallen
%
% Erwartete CSV: komet_karthesisch.csv mit Spalten: t, x, y, (z optional)
%   - t in Sekunden
%   - x,y,(z) in Metern
% Kommentarzeilen mit % sind erlaubt.
%
% Alles in EINER Datei (inkl. earthPos/fullG/ODE), damit nix durcheinanderkommt.

clear; clc;

AU   = 1.495978707e11;       % m
DAY  = 24*3600;
YEAR = 365.25*DAY;

%% =======================
%  0) Einstellungen
%% =======================
csvFile = "komet_karthesisch.csv";

fit_r0 = false;      % true: (x0,y0) auch fitten; false: r0 = erster Messpunkt (stabiler)

qMinAU  = 4.0;        % "realistische" minimale Sonnendistanz (Perihel) in AU
lambdaQ = 1e7;        % Stärke der Penalty fürs Perihel (größer = härter)
lambdaHyp = 1e8;      % harte Penalty falls hyperbolisch (E>=0)

nFine   = 50000;      % für Plot + Min-Abstand

%% =======================
%  1) Daten robust laden
%% =======================
D = readmatrix(csvFile, "CommentStyle","%");
D = D(all(isfinite(D),2), :);

if size(D,2) < 3
    error("CSV hat zu wenige Spalten. Erwartet: t,x,y,(z).");
end

t_obs = D(:,1);
r_obs = D(:,2:3);     % 2D: x,y

% sortieren
[t_obs, idx] = sort(t_obs);
r_obs = r_obs(idx,:);

% Kurzer Import-Check (damit du sofort siehst, ob Units passen)
r_obs_norm_AU = vecnorm(r_obs,2,2)/AU;
fprintf("Import-Check:\n");
fprintf("  t-span: %.2f Jahre\n", (t_obs(end)-t_obs(1))/YEAR);
fprintf("  |r_obs|: min=%.6f AU, max=%.6f AU\n\n", min(r_obs_norm_AU), max(r_obs_norm_AU));

%% =======================
%  2) Startwerte
%% =======================
r0_guess = r_obs(1,:)';

dt = t_obs(2) - t_obs(1);
v0_guess = (r_obs(2,:)' - r_obs(1,:)') / dt;

% Skaliert fitten (stabiler)
% Wenn fit_r0=false: q = [vx0/1000; vy0/1000]
% Wenn fit_r0=true : q = [x0/AU; y0/AU; vx0/1000; vy0/1000]
if fit_r0
    q0 = [r0_guess(1)/AU; r0_guess(2)/AU; v0_guess(1)/1000; v0_guess(2)/1000];
else
    q0 = [v0_guess(1)/1000; v0_guess(2)/1000];
end

%% =======================
%  3) Optimierung (Best-fit)
%% =======================
obj = @(q) costFun(q, t_obs, r_obs, r0_guess, fit_r0, AU, qMinAU, lambdaQ, lambdaHyp);

opts = optimset('Display','iter', 'MaxIter',6000, 'MaxFunEvals',25000);
q_best = fminsearch(obj, q0, opts);

% zurückrechnen
if fit_r0
    r0 = [q_best(1)*AU; q_best(2)*AU];
    v0 = [q_best(3)*1000; q_best(4)*1000];
else
    r0 = r0_guess;
    v0 = [q_best(1)*1000; q_best(2)*1000];
end

fprintf("\n=== Best-fit Anfangsbedingungen (FULL + Penalty) ===\n");
fprintf("r0 = [%.6e; %.6e] m (= [%.6f; %.6f] AU)\n", r0(1), r0(2), r0(1)/AU, r0(2)/AU);
fprintf("v0 = [%.6e; %.6e] m/s\n", v0(1), v0(2));

% Perihelabschätzung (2-Body Sonne)
[qAU, aAU, e, E] = perihelApproxAU(r0, v0);
fprintf("Perihel approx: q=%.6f AU, a=%.6f AU, e=%.6f, E=%+.3e\n", qAU, aAU, e, E);

%% =======================
%  4) Beste Bahn integrieren
%% =======================
y0 = [r0; v0];

odeOpts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',5e6); % ~58 Tage
sol = ode113(@kometODE_full_local, [t_obs(1) t_obs(end)], y0, odeOpts);

% RMSE an Beobachtungszeiten
Yobs = deval(sol, t_obs).';      % Nx4
r_fit_obs = Yobs(:,1:2);
rmse = sqrt(mean(sum((r_fit_obs - r_obs).^2, 2)));
fprintf("\nFit-RMSE: %.6e m (= %.6f AU)\n", rmse, rmse/AU);

%% =======================
%  5) Minimaler Abstand Erde–Komet
%% =======================
tFine = linspace(t_obs(1), t_obs(end), nFine);
Yfine = deval(sol, tFine);
rC = Yfine(1:2,:);               % 2xN (m)
rE = earthPos_local(tFine);      % 2xN (m)

dist = vecnorm(rC - rE, 2, 1);
[minDist, kMin] = min(dist);

fprintf("\n=== Minimaler Abstand Komet–Erde (im Beobachtungszeitraum) ===\n");
fprintf("d_min = %.6e m (= %.6f AU)\n", minDist, minDist/AU);
fprintf("t_min = %.6e s (= %.2f Jahre nach Start)\n", tFine(kMin), tFine(kMin)/YEAR);

% Sanity: min Abstand zur Sonne (numerisch aus feiner Kurve)
rSunMinAU = min(vecnorm(rC,2,1))/AU;
fprintf("Sanity: min(|rC|) numerisch = %.6f AU\n", rSunMinAU);

%% =======================
%  6) Plot: Gesamt + Zoom + d_min markieren
%% =======================
rC_AU = rC / AU;
rE_AU = rE / AU;

figure; hold on; grid on; axis equal;
plot(r_obs(:,1)/AU, r_obs(:,2)/AU, 'ko', 'DisplayName','Beobachtungen');
plot(rC_AU(1,:), rC_AU(2,:), 'r-', 'LineWidth',1.5, 'DisplayName','Fit-Bahn (glatt)');
plot(rE_AU(1,:), rE_AU(2,:), 'b--', 'DisplayName','Erde');
plot(0,0,'yo','MarkerFaceColor','y','DisplayName','Sonne');

% d_min markieren
plot(rC_AU(1,kMin), rC_AU(2,kMin), 'rp', 'MarkerFaceColor','r', 'DisplayName','Komet bei d_{min}');
plot(rE_AU(1,kMin), rE_AU(2,kMin), 'bp', 'MarkerFaceColor','b', 'DisplayName','Erde bei d_{min}');
plot([rC_AU(1,kMin) rE_AU(1,kMin)], [rC_AU(2,kMin) rE_AU(2,kMin)], 'k-', 'LineWidth',1.2, 'DisplayName','d_{min}');

xlabel('x [AU]'); ylabel('y [AU]');
title('Kometenbahn-Fit aus Beobachtungsdaten (FULL + realistische Penalty)');
legend('Location','best');
xlim([-2 55]); ylim([-15 15]);

% Zoom um die Sonne (damit man Perihel wirklich sieht)
figure; hold on; grid on; axis equal;
plot(rC_AU(1,:), rC_AU(2,:), 'r-', 'LineWidth',1.5);
plot(rE_AU(1,:), rE_AU(2,:), 'b--');
plot(0,0,'yo','MarkerFaceColor','y');
xlabel('x [AU]'); ylabel('y [AU]');
title('Zoom um die Sonne');
xlim([-1.0 6.0]); ylim([-2.0 2.0]);

end

%% =====================================================================
%  Cost-Funktion: SSE (Beobachtung vs Modell) + Penalty für Perihel/Hyperbel
%% =====================================================================
function J = costFun(q, t_obs, r_obs, r0_fixed, fit_r0, AU, qMinAU, lambdaQ, lambdaHyp)

% Parameter zurückrechnen
if fit_r0
    r0 = [q(1)*AU; q(2)*AU];
    v0 = [q(3)*1000; q(4)*1000];
else
    r0 = r0_fixed;
    v0 = [q(1)*1000; q(2)*1000];
end

y0 = [r0; v0];

% Integration direkt auf Beobachtungszeiten (schnell & simpel)
odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',5e6);
[~, Y] = ode113(@kometODE_full_local, t_obs, y0, odeOpts);

rModel = Y(:,1:2);     % Nx2
err = (rModel - r_obs) / AU;

% SSE normalisiert
J = sum(err(:,1).^2 + err(:,2).^2);

% --- "Realistische Bahn" Penalty: Perihel zur Sonne nicht unter qMinAU
[qAU, ~, ~, E] = perihelApproxAU(r0, v0);

if qAU < qMinAU
    J = J + lambdaQ * (qMinAU - qAU)^2;
end

% --- Harte Penalty, wenn hyperbolisch (E >= 0)
if E >= 0
    J = J + lambdaHyp * (1 + E)^2;
end

end

%% =====================================================================
%  Perihel-Approximation (2-Body Sonne, 2D)
%% =====================================================================
function [qAU, aAU, e, E] = perihelApproxAU(r0, v0)
AU  = 1.495978707e11;
muS = 1.32712440018e20;

r = norm(r0);
v = norm(v0);

% 2D Drehimpuls (Skalar)
h = r0(1)*v0(2) - r0(2)*v0(1);

% spezifische Energie
E = 0.5*v^2 - muS/r;

% große Halbachse (nur sinnvoll wenn E<0)
if E < 0
    a = -muS/(2*E);
else
    a = inf;
end

% Exzentrizität
e = sqrt( max(0, 1 + (2*E*h^2)/(muS^2)) );

% Perihel
if isfinite(a)
    q = a*(1 - e);
else
    q = 0; % bei hyperbolisch: "invalid", wird sowieso hart bestraft
end

qAU = q/AU;
aAU = a/AU;
end

%% =====================================================================
%  FULL-ODE (2D): Sonne + Erde
%% =====================================================================
function dydt = kometODE_full_local(t, y)
r = y(1:2);
v = y(3:4);
a = fullG_local(r, t);
dydt = [v; a];
end

function a = fullG_local(r, t)
% Sonne
muS = 1.32712440018e20;     % m^3/s^2
rn  = norm(r);
aS  = -muS * r / rn^3;

% Erde
muE = 3.986004418e14;       % m^3/s^2
rE  = earthPos_local(t);
d   = r - rE;
dn  = norm(d);
aE  = -muE * d / dn^3;

a = aS + aE;
end

function rE = earthPos_local(t)
% Kreisbahn Erde um Sonne, 2D, heliiozentrisch
AU = 1.495978707e11;
T  = 365.25*24*3600;
omega = 2*pi/T;

x = -AU * sin(omega*t);
y =  AU * cos(omega*t);

rE = [x; y];
end
