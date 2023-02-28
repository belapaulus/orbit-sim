function [pos, vel, mass] = n_body_video(duration, dt, N, mass, pos0, vel0) 
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
%

% constants
G = 6.6743e-11;
steps = ceil(duration/dt);
camera = max(pos0, [], "all"); % for creating the camera 

% arrays for positions and velocities
pos = zeros(N, steps, 3);
vel = zeros(N, steps, 3);

% initial values
vel(:, 1, :) = vel0;
pos(:, 1, :) = pos0;

% frame of reference
ref = sum(vel(:, 1, :).*mass', 1) / sum(mass);
vel(:, 1, :) = vel(:, 1, :) - ref;

figure

% run sim
for step = 1:steps-1
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

    % clear the current frame
    clf

    % plot the current positions and trails
    txt = ['Step: ' num2str(step)];
    subtitle(txt);
    hold on
    for i = 1:N
	    plot3( ...
            squeeze(pos(i, 1:step, 1)), ...
            squeeze(pos(i, 1:step, 2)), ...
            squeeze(pos(i, 1:step, 3)) ...
        )
        plot3( ...
            squeeze(pos(i, step, 1)), ...
            squeeze(pos(i, step, 2)), ...
            squeeze(pos(i, step, 3)), ...
            'r*' ...
        )
    end    
    hold off

    % set axis limits
    axis equal
    xlim([-1.5*camera,1.5*camera])
    ylim([-1.5*camera,1.5*camera])
    zlim([-1.5*camera,1.5*camera])
    view([ ...
        0.5*camera, ...
        0.5*camera, ...
        0.2*camera ...
        ])
    shg
end
