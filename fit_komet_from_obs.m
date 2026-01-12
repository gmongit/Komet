function fit_komet_from_obs()
clear; clc;

AU = 1.495978707e11;

%% Beobachtungsdaten laden (CSV hat Header wie "% t")
T = readtable("komet_karthesisch.csv", "VariableNamingRule","preserve");

t_obs = T{:,1};          % Sekunden
r_obs = T{:,2:3};        % Meter (x,y)

% Start-Schätzer aus den ersten Punkten
r0_guess = r_obs(1,:)';
v0_guess = (r_obs(2,:)' - r_obs(1,:)') / (t_obs(2) - t_obs(1));

% Skaliert optimieren (stabiler für fminsearch):
% q = [x0/AU, y0/AU, vx0/1000, vy0/1000]
q0 = [r0_guess(1)/AU; r0_guess(2)/AU; v0_guess(1)/1000; v0_guess(2)/1000];

obj = @(q) costFun(q, t_obs, r_obs, AU);

opts = optimset('Display','iter','MaxIter',4000,'MaxFunEvals',12000);
q_best = fminsearch(obj, q0, opts);

% zurückrechnen auf SI
r0 = [q_best(1)*AU; q_best(2)*AU];
v0 = [q_best(3)*1000; q_best(4)*1000];

fprintf("\nBest-fit Anfangsbedingungen:\n");
fprintf("r0 = [%.6e; %.6e] m\n", r0(1), r0(2));
fprintf("v0 = [%.6e; %.6e] m/s\n", v0(1), v0(2));

%% Beste Bahn integrieren
y0 = [r0; v0];
tspan = [t_obs(1) t_obs(end)];
odeOpts = odeset('RelTol',1e-9,'AbsTol',1e-9);

sol = ode45(@kometODE_full, tspan, y0, odeOpts);

% Modell an Beobachtungszeiten auswerten
y_fit = deval(sol, t_obs);
r_fit = y_fit(1:2,:)';

rmse = sqrt(mean(sum((r_fit - r_obs).^2,2)));
fprintf("\nFit-RMSE: %.6e m (%.6f AU)\n", rmse, rmse/AU);

%% Minimaler Abstand Komet–Erde
tFine = linspace(t_obs(1), t_obs(end), 50000);
yFine = deval(sol, tFine);
rC = yFine(1:2,:);          % Komet
rE = earthPos(tFine);       % Erde

dist = sqrt(sum((rC - rE).^2,1));
[minDist, idx] = min(dist);

fprintf("\nMinimaler Abstand Komet–Erde:\n");
fprintf("d_min = %.6e m = %.6f AU\n", minDist, minDist/AU);
fprintf("Zeitpunkt t = %.6e s (%.2f Jahre nach Start)\n", tFine(idx), tFine(idx)/(365.25*24*3600));

%% Plot
figure; hold on; grid on; axis equal;
plot(r_obs(:,1)/AU, r_obs(:,2)/AU, 'ko', 'DisplayName','Beobachtungen');
plot(r_fit(:,1)/AU, r_fit(:,2)/AU, 'r-', 'LineWidth',1.5, 'DisplayName','Fit-Bahn');

rEplot = earthPos(tFine);
plot(rEplot(1,:)/AU, rEplot(2,:)/AU, 'b--', 'DisplayName','Erde');

plot(0,0,'yo','MarkerFaceColor','y','DisplayName','Sonne');

xlabel('x [AU]'); ylabel('y [AU]');
legend('Location','best');
title('Kometenbahn-Fit aus Beobachtungsdaten');
end

function sse = costFun(q, t_obs, r_obs, AU)
% q = [x0/AU, y0/AU, vx0/1000, vy0/1000]

r0 = [q(1)*AU; q(2)*AU];
v0 = [q(3)*1000; q(4)*1000];
y0 = [r0; v0];

tspan = [t_obs(1) t_obs(end)];
odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-8); % im Fit etwas lockerer = schneller

sol = ode45(@kometODE_full, tspan, y0, odeOpts);

yModel = deval(sol, t_obs);
rModel = yModel(1:2,:)';

err = rModel - r_obs;

% Sum of squared errors (normiert, damit Zahlen net explodieren)
sse = sum(sum((err./AU).^2));
end
