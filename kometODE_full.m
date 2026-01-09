function dydt = kometODE_full(t, y)
% Zustand y = [x; y; vx; vy]

r = y(1:2);
v = y(3:4);

a = fullG(r, t);

dydt = [v; a];
end