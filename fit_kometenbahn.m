function fit_kometenbahn()
% Bestimmung einer Kometenbahn aus kartesischen Beobachtungsdaten
% Ziel:
% - bestpassende Bahn durch Variation der Anfangsgeschwindigkeit
% - minimalen Abstand zwischen Komet und Erde bestimmen

% Konstanten 
% Astronomische Einheit und Zeitumrechnungen
AU   = 1.495978707e11;     % Meter
DAY  = 24*3600;            % Sekunden
YEAR = 365.25*DAY;         % Sekunden

% Einstellungen
% Minimaler akzeptierter Sonnenabstand (Plausibilitätscheck)
qMinAU = 2;

% Maximale Korrektur der Anfangsgeschwindigkeit (Stabilität)
dvMax = 2000;

% Feines Zeitraster für Abstand Erde–Komet
nFine = 50000;

% ODE-Einstellungen für den Fit (schneller, robuster)
fitRelTol  = 1e-6;
fitAbsTol  = [1e7 1e7 1e-1 1e-1];
fitMaxStep = 5e6;

% ODE-Einstellungen für die finale Bahn (genauer)
finRelTol  = 1e-7;
finAbsTol  = [1e5 1e5 1e-2 1e-2];
finMaxStep = 2e7;   % größer erlauben

% Optimierer-Einstellungen
opt = optimset('Display','iter','MaxIter',2000,'MaxFunEvals',8000);

% CSV einlesen (Zeit, x, y)
D = readmatrix("komet_karthesisch.csv","CommentStyle","%");
D = D(all(isfinite(D),2), :);

% Zeit und Position extrahieren
t_raw = D(:,1);
r_obs = D(:,2:3);

t_obs = t_raw;   % Sekunden

% Daten nach Zeit sortieren
[t_obs, idx] = sort(t_obs);
r_obs = r_obs(idx,:);

fprintf("Span Jahre: %.2f\n", (t_obs(end)-t_obs(1))/YEAR);


% Startposition: erster Beobachtungspunkt
r0 = r_obs(1,:)';

% Startgeschwindigkeit aus Differenz der ersten zwei Punkte
dt = t_obs(2) - t_obs(1);
v0_guess = (r_obs(2,:)' - r_obs(1,:)') / dt;

% Bahn-Fit
% Es wird nur eine Korrektur dv zur Anfangsgeschwindigkeit gefittet
dv0 = [0;0];

% Zielfunktion für den Optimierer
obj = @(dv) costFun_dv(dv, v0_guess, t_obs, r_obs, r0, ...
                       AU, qMinAU, dvMax, ...
                       fitRelTol, fitAbsTol, fitMaxStep);

% Optimierung starten
dv_best = fminsearch(obj, dv0, opt);

% Beste Anfangsgeschwindigkeit
v0 = v0_guess + dv_best;

% Finale Bahn
% Anfangszustand
y0 = [r0; v0];

% Finale, genaue Integration
odeOptsFin = odeset('RelTol',finRelTol,'AbsTol',finAbsTol,'MaxStep',finMaxStep);
sol = ode15s(@kometODE_full, [t_obs(1) t_obs(end)], y0, odeOptsFin);

if sol.x(end) < t_obs(end)
    error("ODE abgebrochen bei t=%.3e, benötigt bis %.3e", sol.x(end), t_obs(end));
end


% Fit-Qualität
% Modellpositionen zu Beobachtungszeiten
Yobs = deval(sol, t_obs).';
r_fit = Yobs(:,1:2);

% RMSE als Maß für die Anpassungsgüte
rmse = sqrt(mean(sum((r_fit - r_obs).^2,2)));

fprintf("Fit-RMSE: %.3e m (%.3f AU)\n", rmse, rmse/AU);

% Minimaler Abstand Erde–Komet
% Feines Zeitraster
tFine = linspace(t_obs(1), t_obs(end), nFine);

% Kometenpositionen
Yfine = deval(sol, tFine);
rC = Yfine(1:2,:);

% Erdpositionen
rE = earthPos(tFine);

% Abstände berechnen
dist = vecnorm(rC - rE, 2, 1);
[minDist, kMin] = min(dist);

fprintf("Minimaler Abstand Erde–Komet: %.3f AU\n", minDist/AU);

% Plausibilitätscheck
% Minimaler Sonnenabstand entlang der Bahn
minSunDist = min(vecnorm(rC,2,1))/AU;
fprintf("Minimaler Sonnenabstand: %.3f AU\n", minSunDist);

% Plot
figure; hold on; grid on; axis equal;

% Beobachtungsdaten
plot(r_obs(:,1)/AU, r_obs(:,2)/AU, 'ko', 'DisplayName','Beobachtungen');

% Gefittete Bahn
plot(rC(1,:)/AU, rC(2,:)/AU, 'r-', 'LineWidth',1.5, 'DisplayName','Kometenbahn');

% Erdbahn
plot(rE(1,:)/AU, rE(2,:)/AU, 'b--', 'DisplayName','Erde');

% Sonne
plot(0,0,'yo','MarkerFaceColor','y','DisplayName','Sonne');

xlabel('x [AU]');
ylabel('y [AU]');
legend('Location','best');
title('Kometenbahn aus kartesischen Beobachtungsdaten');
dt = diff(t_obs);


function J = costFun_dv(dv, v0_guess, t_obs, r_obs, r0, ...
                        AU, qMinAU, dvMax, ...
                        relTol, absTol, maxStep)

% Zielfunktion für den Bahn-Fit
% Minimiert den Positionsfehler zu den Beobachtungen
% Verwirft unplausible Bahnen (zu sonnennah)

% Sicherstellen, dass dv ein Spaltenvektor ist
dv = dv(:);

% Begrenzung von Delta-v
% Verhindert unrealistisch große Anfangsgeschwindigkeiten
if any(abs(dv) > dvMax)
    J = 1e9;
    return;
end

% Anfangszustand
% Anfangsgeschwindigkeit als Korrektur zum Schätzwert
v0 = v0_guess + dv;

% Zustandsvektor: Position und Geschwindigkeit
y0 = [r0; v0];

% ODE-Einstellungen für den Fit
odeOpts = odeset('RelTol',relTol,'AbsTol',absTol,'MaxStep',maxStep);

% Warnungen temporär ausschalten (Solver kann warnen)
warnState = warning;
warning('off','all');

try
    % Numerische Integration
    % Bewegungsgleichung des Kometen über gesamten Zeitraum
    sol = ode113(@kometODE_full, [t_obs(1) t_obs(end)], y0, odeOpts);

    % Falls der Solver vorzeitig abbricht, Lösung verwerfen
    if sol.x(end) < t_obs(end)
        warning(warnState);
        J = 1e10;
        return;
    end

    % Anpassungsfehler
    % Modellpositionen zu den Beobachtungszeiten
    Yobs = deval(sol, t_obs).';
    rModel = Yobs(:,1:2);

    % Quadratischer Fehler in astronomischen Einheiten
    err = (rModel - r_obs) / AU;
    SSE = sum(err(:,1).^2 + err(:,2).^2);

    % Plausibilitätscheck-
    % Minimaler Sonnenabstand über gesamten Zeitraum
    tC = linspace(t_obs(1), t_obs(end), 10000);
    YC = deval(sol, tC).';
    rSun_AU = sqrt(YC(:,1).^2 + YC(:,2).^2) / AU;
    minr = min(rSun_AU);

    % Bahn verwerfen, wenn sie zu sonnennah wird
    if minr < qMinAU
        warning(warnState);
        J = 1e10;
        return;
    end

    % - Zielfunktionswert
    % Nur der Fit-Fehler zählt, wenn Bahn plausibel ist
    J = SSE;

    % Warnungen wiederherstellen
    warning(warnState);

catch
    % Falls numerisch etwas schiefgeht, Lösung verwerfen
    warning(warnState);
    J = 1e10;
end



