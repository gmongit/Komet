function rE = earthPos(t)

AU = 1.495978707e11; %mittlerer Abstand Erde zu Sonne

T  = 365.25*24*3600; %Zeit Jahr in Sekunden

omega = 2*pi/T; %winkelgeschwindigkeit


x = -AU * sin(omega*t);    %x, y Position der berechnen auf Kreisbahn   
y =  AU * cos(omega*t);

rE = [x; y];


end

