function fit_karthesisch_full_dv_final()
% Kartesische Beobachtungen (FULL): Best-fit Bahn + minimaler Abstand zur Erde
% Robust: fitte nur Delta-v um v0_guess (statt v0 direkt)
% Realismus: Perihel-Penalty (Sonne 2-body approx)
%
% Benötigt im selben Ordner:
%   komet_karthesisch.csv
%   kometODE_full.m  (-> fullG.m, earthPos.m)

AU   = 1.495978707e11;
DAY  = 24*3600;
YEAR = 365.25*DAY;

% ================= Einstellungen =================
qMinAU  = 0.1;        % "realistisch": Perihel >= 4 AU
lambdaQ = 0;        % Penalty-Stärke
dvMax   = 200;        % |Delta-v| Grenze je Komponente (m/s) -> Fit bleibt nahe am Guess

nFine   = 50000;      % feines Sampling für d_min

% Fit-Solver (schnell & robust)
fitRelTol = 1e-6;
fitAbsTol = [1e7 1e7 1e-1 1e-1];     % [x y vx vy]
fitMaxStep = 5e6;

% Finale Integration (genauer)
finRelTol = 1e-9;
finAbsTol = [1e3 1e3 1e-4 1e-4];
finMaxStep = 5e6;

% Optimierer
opt = optimset('Display','iter','MaxIter',2000,'MaxFunEvals',8000);

% ================= 1) Daten laden =================
D = readmatrix("komet_karthesisch.csv","CommentStyle","%");
D = D(all(isfinite(D),2), :);

if size(D,2) < 3
    error("CSV hat zu wenige Spalten. Erwartet: t,x,y,(z).");
end

t_obs = D(:,1);
r_obs = D(:,2:3);  % nur x,y

[t_obs, idx] = sort(t_obs);
r_obs = r_obs(idx,:);

fprintf("Import: t-span = %.2f Jahre, min(|r_obs|)=%.6f AU\n", ...
    (t_obs(end)-t_obs(1))/YEAR, min(vecnorm(r_obs,2,2))/AU);

% ================= 2) Startwerte =================
r0 = r_obs(1,:)';   % r0 fix (stabil)

dt = t_obs(2) - t_obs(1);
v0_guess = (r_obs(2,:)' - r_obs(1,:)')/dt;

fprintf("v0_guess = [%.2f, %.2f] m/s  |v0|=%.2f m/s\n", ...
    v0_guess(1), v0_guess(2), norm(v0_guess));

% ================= 3) Fit: Delta-v =================
dv0 = [0;0];  % Start bei "keine Korrektur"

obj = @(dv) costFun_dv(dv, v0_guess, t_obs, r_obs, r0, AU, qMinAU, lambdaQ, dvMax, fitRelTol, fitAbsTol, fitMaxStep);

dv_best = fminsearch(obj, dv0, opt);

v0 = v0_guess + dv_best;

% ===== POSTCHECK: minr aus derselben Dynamik (sollte >= qMinAU sein) =====
odeOptsCheck = odeset('RelTol',fitRelTol,'AbsTol',fitAbsTol,'MaxStep',fitMaxStep);
solCheck = ode113(@kometODE_full, [t_obs(1) t_obs(end)], [r0; v0], odeOptsCheck);

tC = linspace(t_obs(1), t_obs(end), 50000);
YC = deval(solCheck, tC).';
minr_check = min(sqrt(YC(:,1).^2 + YC(:,2).^2)/AU);

fprintf("POSTCHECK minr = %.6f AU (qMinAU=%.3f)\n", minr_check, qMinAU);

fprintf("\n=== Best-fit Anfangsbedingungen (FULL + Realismus, Delta-v Fit) ===\n");
fprintf("r0 = [%.6e; %.6e] m (= [%.6f; %.6f] AU)\n", r0(1), r0(2), r0(1)/AU, r0(2)/AU);
fprintf("dv = [%.2f; %.2f] m/s\n", dv_best(1), dv_best(2));
fprintf("v0 = [%.2f; %.2f] m/s  |v0|=%.2f m/s\n", v0(1), v0(2), norm(v0));

qAU = perihelApproxAU_2D(r0, v0);
fprintf("Perihel approx (Sonne 2-body): q=%.6f AU\n", qAU);

% ================= 4) Final integrieren =================
y0 = [r0; v0];
odeOptsFin = odeset('RelTol',finRelTol,'AbsTol',finAbsTol,'MaxStep',finMaxStep);
sol = ode113(@kometODE_full, [t_obs(1) t_obs(end)], y0, odeOptsFin);

% RMSE an Beobachtungszeiten
Yobs = deval(sol, t_obs).';
r_fit = Yobs(:,1:2);
rmse = sqrt(mean(sum((r_fit - r_obs).^2,2)));

fprintf("\nFit-RMSE: %.6e m (= %.6f AU)\n", rmse, rmse/AU);

% ================= 5) Minimaler Abstand Erde–Komet =================
tFine = linspace(t_obs(1), t_obs(end), nFine);
Yfine = deval(sol, tFine);
rC = Yfine(1:2,:);
rE = earthPos(tFine);

dist = vecnorm(rC - rE, 2, 1);
[minDist, kMin] = min(dist);

fprintf("\n=== Minimaler Abstand Komet–Erde (im Beobachtungszeitraum) ===\n");
fprintf("d_min = %.6e m (= %.6f AU)\n", minDist, minDist/AU);
fprintf("t_min = %.6e s (= %.2f Jahre nach Start)\n", tFine(kMin), tFine(kMin)/YEAR);

fprintf("Sanity: min(|rC|) numerisch = %.6f AU\n", min(vecnorm(rC,2,1))/AU);
fprintf("max |y_obs| = %.3f AU\n", max(abs(r_obs(:,2)))/AU);
fprintf("max |y_fit| = %.3f AU\n", max(abs(rC(2,:)))/AU);

% ================= 6) Plot =================
figure; hold on; grid on; axis equal;
plot(r_obs(:,1)/AU, r_obs(:,2)/AU, 'ko', 'DisplayName','Beobachtungen');
plot(rC(1,:)/AU, rC(2,:)/AU, 'r-', 'LineWidth',1.5, 'DisplayName','Fit-Bahn');
plot(rE(1,:)/AU, rE(2,:)/AU, 'b--', 'DisplayName','Erde');
plot(0,0,'yo','MarkerFaceColor','y','DisplayName','Sonne');

% d_min markieren
plot(rC(1,kMin)/AU, rC(2,kMin)/AU, 'rp', 'MarkerFaceColor','r', 'DisplayName','Komet bei d_{min}');
plot(rE(1,kMin)/AU, rE(2,kMin)/AU, 'bp', 'MarkerFaceColor','b', 'DisplayName','Erde bei d_{min}');
plot([rC(1,kMin) rE(1,kMin)]/AU, [rC(2,kMin) rE(2,kMin)]/AU, 'k-', 'LineWidth',1.2, 'DisplayName','d_{min}');

xlabel('x [AU]'); ylabel('y [AU]');
title('Kometenbahn-Fit (FULL) aus kartesischen Beobachtungen (Delta-v Fit)');
legend('Location','best');
xlim([-2 55]); ylim([-15 15]);

end

% ============================================================
% Cost-Funktion: fit Delta-v (SSE + Perihel-Penalty + Robustheit)
% ============================================================
function J = costFun_dv(dv, v0_guess, t_obs, r_obs, r0, AU, qMinAU, lambdaQ, dvMax, relTol, absTol, maxStep)

dv = dv(:);

% (1) Delta-v Begrenzung: Fit soll nahe am Guess bleiben
if any(abs(dv) > dvMax)
    J = 1e9 + sum((abs(dv)-dvMax).^2);
    return;
end

% (2) v0 setzen
v0 = v0_guess + dv;
y0 = [r0; v0];


% (3) ODE Optionen
odeOpts = odeset('RelTol',relTol,'AbsTol',absTol,'MaxStep',maxStep);

warnState = warning; 
warning('off','all');

try
    % =========================================================
    % EINMAL integrieren über das ganze Intervall
    % =========================================================
    sol = ode113(@kometODE_full, [t_obs(1) t_obs(end)], y0, odeOpts);

    % Wenn Solver nicht bis zum Ende kommt -> Kandidat verwerfen
    if sol.x(end) < t_obs(end)
        warning(warnState);
        J = 1e10;
        return;
    end

    % =========================================================
    % (4) SSE an Beobachtungszeiten (deval)
    % =========================================================
    Yobs = deval(sol, t_obs).';          % Nx4
    if size(Yobs,1) ~= size(r_obs,1) || any(~isfinite(Yobs(:)))
        warning(warnState);
        J = 1e10;
        return;
    end

    rModel = Yobs(:,1:2);
    err = (rModel - r_obs)/AU;
    SSE = sum(err(:,1).^2 + err(:,2).^2);

    % =========================================================
    % (5) Realismus numerisch: min(|r|) im ganzen Zeitraum
    % =========================================================
    tC = linspace(t_obs(1), t_obs(end), 800); % 400-1000 ok
    YC = deval(sol, tC).';                    % 800x4
    rSun_AU = sqrt(YC(:,1).^2 + YC(:,2).^2)/AU;
    minr = min(rSun_AU);
    fprintf("COST minr = %.3f AU (qMinAU=%.3f)\n", minr, qMinAU);


    if minr < qMinAU
    fprintf("PENALTY HIT: minr=%.3f AU < %.3f AU\n", minr, qMinAU);
end


    if minr < qMinAU
       J = 1e10;      % Kandidat sofort verwerfen
       return;
    
end


J = SSE;           % nur noch Fit-Fehler, weil Realismus erfüllt ist

    warning(warnState);

catch
    warning(warnState);
    J = 1e10;
end
end



% ============================================================
% Perihel-Approximation (Sonne 2-body, 2D) für Realismus-Penalty
% ============================================================
function qAU = perihelApproxAU_2D(r0, v0)
AU  = 1.495978707e11;
muS = 1.32712440018e20;

r = norm(r0);
v = norm(v0);

h = r0(1)*v0(2) - r0(2)*v0(1);
E = 0.5*v^2 - muS/r;

if E >= 0
    qAU = 0;      % dann greift die Penalty sowieso stark
    return;
end

a = -muS/(2*E);
e = sqrt(max(0, 1 + (2*E*h^2)/(muS^2)));
q = a*(1 - e);

qAU = q/AU;
end
