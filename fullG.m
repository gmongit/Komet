function g = fullG(r, t)
%FULLG Gravitationsfeld (Beschleunigung) von Sonne + Erde.
%   g = fullG(r,t)
%   r ... Position (2x1 oder 1x2) im Sonnensystem [m]
%   t ... Zeitpunkt [s]
%   g ... Beschleunigung am Punkt r durch Sonne + Erde [m/s^2]
%
%   Voraussetzung: Ihr habt eine Funktion earthPos(t),
%   die die Erdposition relativ zur Sonne als 2x1 Vektor [m] liefert.

    % Konstanten (SI)
    G       = 6.67430e-11;      % [m^3/(kg*s^2)]
    M_earth = 5.9722e24;        % [kg]

    r = r(:);

    % 1) Beitrag Sonne (Sonne sitzt im Ursprung)
    g = simpleG(r);

    % 2) Beitrag Erde (Erde sitzt bei rE(t))
    rE = earthPos(t);     % <-- muss von euch kommen, in [m]
    rE = rE(:);

    dvec = r - rE;        % Vektor von Erde zum Punkt
    d    = norm(dvec);

    if d == 0
        error('fullG: Punkt liegt exakt auf der Erde (Abstand = 0).');
    end

    % Beschleunigung durch Erde addieren
    g = g - G * M_earth / d^3 * dvec;
end
