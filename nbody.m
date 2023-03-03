function [pos, vel] = nbody(duration, dt, N, mass, pos0, vel0) 
% n body simulation
% units: 
%   s, m, kg
% input:
%   initial values and constants
% side effects:
%   creates moving figure
% output:
%   arrays containing the state of the system for each timestep
%   array containing the mass of each body

% constants
G = 6.6743e-11;
steps = ceil(duration/dt);

% arrays for positions and velocities
pos = zeros(N, steps, 3);
vel = zeros(N, steps, 3);

% initial values
vel(:, 1, :) = vel0;
pos(:, 1, :) = pos0;

% frame of reference
ref = sum(vel(:, 1, :).*mass', 1) / sum(mass);
vel(:, 1, :) = vel(:, 1, :) - ref;

% run sim
for step = 1:steps-1
    if mod(step, 100) == 0
        step
    end

    % update ith body
	for i = 1:N
		vel(i, step+1, :) = vel(i, step, :);
		for j = 1:N
			if i ~= j
                % calculate distance r = |posj - posi|
                dpos = pos(j, step, :) - pos(i, step, :);
				r = sqrt(sum(dpos.^2));
				vel(i, step + 1, :) = vel(i, step + 1, :) + dt*G*mass(j)*dpos / (r^3);
			end
        end
		pos(i, step + 1, :) = pos(i, step, :) + dt*vel(i, step + 1, :);
    end
end
