% n body simulation
% units: s, m, kg

% timestep 1 day = 24 * 60 * 60s
duration = 5e7;
dt = 24 * 60 * 60;
steps = ceil(duration/dt);

% constants
G = 6.6743e-11;
N = 3;
mass = [2e30, 6e29, 4e24];

% arrays for positions and velocities
pos = zeros(steps, 3, N);
vel = zeros(steps, 3, N);

% initial values
vel(1, :, 1) = [0, -6923.076923076924, 0];
vel(1, :, 2) = [0, 30e3 - 6923.076923076924, 5000];
vel(1, :, 3) = [5000, 0, 0];
pos(1, :, 2) = [150e9, 0, 0];
pos(1, :, 3) = [0, 150e9, 0];

% frame of reference
ref = zeros(1,3);
for i = 1:N
	ref = ref + (squeeze(vel(1, :, i)).*mass(i));
end
ref = ref./sum(mass);
for i = 1:N
	vel(1, :, i) = vel(1, :, i) - ref;
end

% create video object
writerObj = VideoWriter("testvideo.mp4", 'MPEG-4'); % see VideoWriter docs for other formats
writerObj.FrameRate = 30; % set frame rate for video
open(writerObj)
figure

% run sim
for step = 1:steps-1
    % update ith body
	for i = 1:N
		vel(step + 1, :, i) = vel(step, :, i);
		for j = 1:N
			if i ~= j
                % calculate distance r = |posj - posi|
                dpos = pos(step, :, j) - pos(step, :, i);
				r = sqrt(sum(dpos.^2));
				vel(step + 1, :, i) = vel(step + 1, :, i) + dt*G*mass(j)*dpos / (r^3);

				pos(step + 1, :, i) = pos(step, :, i) + dt*vel(step + 1, :, i);
			end
		end
    end

    % clear the current frame
    clf

    % plot the current positions and trails
    hold on
    for i = 1:N
	    plot3( ...
            squeeze(pos(1:step, 1, i)), ...
            squeeze(pos(1:step, 2, i)), ...
            squeeze(pos(1:step, 3, i)) ...
        )
        plot3( ...
            squeeze(pos(step, 1, i)), ...
            squeeze(pos(step, 2, i)), ...
            squeeze(pos(step, 3, i)), ...
            'r*' ...
        )
    end    
    hold off

    % set axis limits so they don't change through movie frames
    axis equal
    xlim([-1.5*x(1),1.5*x(1)])
    ylim([-1.5*x(1),1.5*x(1)])
    zlim([-1.5*x(1),1.5*x(1)])
    view([ ...
        0.5*x(1), ...
        0.5*x(1), ...
        0.2*x(1) ...
        ])
    
    % write this figure to the video as a new frame
    currentframe = getframe();
    writeVideo(writerObj, currentframe)
end

% close video (saves video)
close(writerObj)