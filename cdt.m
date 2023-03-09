function [pos, vel, vel2, dvel, timestamp] = cdt(duration, k, N, mass, pos0, vel0) 
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
vel2 = zeros(N, 1, 3); % array for saving the velocities at each t + dt/2
dvel = zeros(N, 1, 3);
timestamp = zeros(1);

% initial values
vel(:, 1, :) = vel0;
pos(:, 1, :) = pos0;

% set the frame of reference
vref = sum(vel(:, 1, :).*mass', 1) / sum(mass);
vel(:, 1, :) = vel(:, 1, :) - vref;
pref = sum(pos(:, 1, :).*mass', 1) / sum(mass);
pos(:, 1, :) = pos(:, 1, :) - pref;

% calculate the change in velocity for each body
% we are doing this here once for the initial change in velocity
for i = 1:N-1
    for j = i+1:N
        % calculate the vector from i to j and its length
        dpos = squeeze(pos(j, 1, :) - pos(i, 1, :));
		r = norm(dpos);
        
        % calculate the change in velocity for i and j
		dvel(i, 1, :) = squeeze(dvel(i, 1, :)) + G*mass(j)*dpos/(r^3);
        dvel(j, 1, :) = squeeze(dvel(j, 1, :)) - G*mass(i)*dpos/(r^3);
    end
end

% run sim
step = 1;
time = 0;
while time < duration 
    if mod(step, 100) == 0
        step
    end
    
    timestamp(step) = time;

    % calculate next time step 
    vnorm = zeros(N, 1);
    for i = 1:N
        vnorm(i) = norm(squeeze(vel(1, step, :))); 
    end
    dt = k/max(vnorm);
    
    % calculate the velocity at t + dt/2
    % v(t+dt/2) = v(t) + dt/2 * F(x(t))/m
    vel2(:, step, :) = vel(:, step, :) + dt/2*dvel(:, step, :);
    % calculate the position at t + dt
    % x(t+dt) = x(t) + dt*v(t+dt/2)
    pos(:, step+1, :) = pos(:, step, :) + dt*vel2(:, step, :);

    % calculate the change in velocity for each body at the next half step
    % this will be used here but also in the next iteration    
    dvel(:, step+1, :) = 0;
    for i = 1:N-1
        for j = i+1:N
            % calculate the vector from i to j and its length
            dpos = squeeze(pos(j, step+1, :) - pos(i, step+1, :));
			r = norm(dpos);
            
            % calculate the change in velocity for i and j
			dvel(i, step+1, :) = squeeze(dvel(i, step+1, :)) + G*mass(j)*dpos/(r^3);
            dvel(j, step+1, :) = squeeze(dvel(j, step+1, :)) - G*mass(i)*dpos/(r^3);
        end
    end

    % calculate the velocity at t + dt
    % v(t+dt) = v(t+dt/2) + dt/2 * F(x(t+dt))/m
    vel(:, step+1, :) = vel2(:, step, :) + dt/2 * dvel(:, step+1, :);

    time = time + dt;
    step = step + 1;
end