function g = simpleG(r)
%SIMPLEG Gravitationsfeld (Beschleunigung) nur von der Sonne.
%   g = simpleG(r)
%   r ... Position (2x1 oder 1x2) relativ zur Sonne [m]
%   g ... Beschleunigung am Punkt r durch Sonne [m/s^2]
%
%   Richtung: immer ZUM Massezentrum => Minuszeichen!

    % Konstanten (SI)
    G      = 6.67430e-11;       % [m^3/(kg*s^2)]
    M_sun  = 1.98847e30;        % [kg]

    % Sicherstellen: Spaltenvektor
    r = r(:);

    % Abstand
    d = norm(r);

    % Schutz gegen Division durch 0
    if d == 0
        error('simpleG: r darf nicht [0;0] sein (Abstand = 0).');
    end

    % Gravitationsbeschleunigung: g = -G*M/|r|^3 * r
    g = -G * M_sun / d^3 * r;
end
