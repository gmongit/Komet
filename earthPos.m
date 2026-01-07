function rE = earthPos(t)
% earthPos  Position der Erde auf Kreisbahn um die Sonne (SI-Einheiten)
% Input:
%   t  ... Zeit in Sekunden (Skalar oder Vektor)
% Output:
%   rE ... 2xN Matrix: [x; y] in Metern

AU = 1.495978707e11;          % m
T  = 365.25*24*3600;          % s
omega = 2*pi/T;

x = -AU * sin(omega*t);       % erfüllt x(0)=0
y =  AU * cos(omega*t);

rE = [x; y];
end

