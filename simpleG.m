function g = simpleG(r)

    G      = 6.67430e-11;   %Gravitationskonstante in m³/(kg·s²)
    M_sun  = 1.98847e30;    %Masse der Sonne in kg    

   
    r = r(:); %Spaltenvektor

    % Abstand
    d = norm(r); %Abstand vom Punkt r zur Sonne(a^2 + b^2)

    if d == 0
        error('simpleG: r darf nicht [0;0] sein (Abstand = 0).');
    end

    g = -G * M_sun / d^3 * r; %je größer d desto geringer g
end
