function [et, ek, ep] = en(vel, pos, mass)
% takes the output of a simulation as input 
% vel: array of the velocities of the form array(body, timestep, component)
% pos: array of the positions of the form array(body, timestep, component)
% mass: array of the masses of the form array(body)
% calculates the kinetic, potential and total energy of the system for
% every timestep

% constants
G = 6.6743e-11;
N = size(vel, 1);
steps = size(vel, 2);

% initialization
ek = zeros(steps, 1);
ep = zeros(steps, 1);
et = zeros(steps, 1);

for step = 1:steps-1
    % calculate the kinetic energy for each body
    for i = 1:N
        ekin = 0.5 * mass(i) * norm(squeeze(vel(i, step, :)))^2;
        ek(step) = ek(step) + ekin;
    end

    % calculate the potential energy for each pair of bodies
    for i = 1:N-1
        epot = 0;
        for j = i+1:N
            dpos = pos(j, step, :) - pos(i, step, :);
			r = sqrt(sum(dpos.^2));
            epot = epot - G*mass(j)*mass(i)/r;
        end
        ep(step) = ep(step) + epot;
    end
    
    et(step) = ek(step) + ep(step);
end
hold on;
plot(et);
plot(ek);
plot(ep);

