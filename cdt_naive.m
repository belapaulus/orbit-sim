function [pos, vel, dvel, timestamp] = cdt_inaccurate(duration, k, N, mass, pos0, vel0) 
% n body simulation
% input:
%   initial values and constants
% output:
%   arrays containing the state of the system for each timestep
%   array containing the mass of each body


% constants
G = 6.6743e-11;

% arrays for positions and velocities and timestamps
pos = zeros(N, 1, 3);
vel = zeros(N, 1, 3);
dvel = zeros(N, 1, 3);
timestamp = zeros(1);

% initial values
vel(:, 1, :) = vel0;
pos(:, 1, :) = pos0;

% set the frame of reference
ref = sum(vel(:, 1, :).*mass', 1) / sum(mass);
vel(:, 1, :) = vel(:, 1, :) - ref;

% run sim
step = 1;
time = 0;
while time < duration 
    if mod(step, 100) == 0
        step
    end
    
    timestamp(step) = time;
    dvel(:, step, :) = 0;

    % calculate the change in velocity for each body
    for i = 1:N-1
        for j = i+1:N
            % calculate the vector from i to j and its length
            dpos = squeeze(pos(j, step, :) - pos(i, step, :));
			r = norm(dpos);
            
            % calculate the change in velocity for i and j
			dvel(i, step, :) = squeeze(dvel(i, step, :)) + G*mass(j)*dpos/(r^3);
            dvel(j, step, :) = squeeze(dvel(j, step, :)) - G*mass(i)*dpos/(r^3);
        end
    end

    % calculate next time step 
    vnorm = zeros(N, 1);
    for i = 1:N
        vnorm(i) = norm(squeeze(vel(1, step, :))); 
    end
    dt = k/max(vnorm);

    % update the state
    for i = 1:N
        vel(i, step + 1, :) = vel(i, step, :) + dt*dvel(i, step, :);
		pos(i, step + 1, :) = pos(i, step, :) + dt*vel(i, step + 1, :);
    end

    time = time + dt;
    step = step + 1;
end