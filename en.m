function [et, ek, ep] = en(vel, pos, mass)
G = 6.6743e-11;
N = size(vel, 1);
steps = size(vel, 2);
ek = zeros(steps, 1);
ep = zeros(steps, 1);
et = zeros(steps, 1);
for step = 1:steps-1
    for i = 1:N
        ekin = 0.5 * mass(i) * norm(squeeze(vel(i, step, :)))^2;
        epot = 0;
        for j = i:N
            if i ~= j
                dpos = pos(j, step, :) - pos(i, step, :);
				r = sqrt(sum(dpos.^2));
                epot = epot - G*mass(j)*mass(i)/r;
            end
        end
        ek(step) = ek(step) + ekin;
        ep(step) = ep(step) + epot;
    end
    et(step) = ek(step) + ep(step);
end
hold on;
plot(et);
plot(ek);
plot(ep);

