function dydt = kometODE_simple(t, y)
% Zustand y = [x; y; vx; vy]

r = y(1:2);
v = y(3:4);

a = simpleG(r);

dydt = [v; a];
end