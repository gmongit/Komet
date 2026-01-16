function g = fullG(r, t)
%FULLG Beschleunigung durch Sonne + Erde (2D, SI)



    % Konstanten
    muS = 1.32712440018e20;   % GM Sonne [m^3/s^2] G*m
    muE = 3.986004418e14;     % GM Erde [m^3/s^2]  G*m

    r = r(:);

    % Sonne (im Ursprung)
    rn = norm(r);
    gS = -muS * r / rn^3;

    % Erde
    rE = earthPos(t);   
    rE = rE(:);

    dvec = r - rE;
    d    = norm(dvec);

    % (kleine Sicherheitsklammer, damit nix explodiert falls d extrem klein)
    d = max(d, 1e7);      % 10.000 km

    gE = -muE * dvec / d^3;

    g = gS + gE;
end
