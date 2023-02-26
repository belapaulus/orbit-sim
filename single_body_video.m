% single body simulation
% units: s, m, kg

% timestep 1 day = 24 * 60 * 60s
% simulation length 1 year = 356 * 1 day
dt = 24 * 60 * 60;
nT = 365;

% constants
G = 6.6743e-11;
M = 1.989e30;

% arrays for position and speed
x = zeros(nT + 1);
y = zeros(nT + 1);
u = zeros(nT + 1);
v = zeros(nT + 1);

% initial position and speed
x(1) = 1.4809e11;
y(1) = 0;
u(1) = 0;
v(1) = 30e3; % sqrt(4*G*M/450e9);

% create video object
writerObj = VideoWriter("testvideo.mp4", 'MPEG-4'); % see VideoWriter docs for other formats
writerObj.FrameRate = 10; % set frame rate for video
open(writerObj)

figure
for i = 1:nT
    % update speed
    r = sqrt(x(i)^2 + y(i)^2);
    u(i + 1) = u(i) - dt * G * M * (x(i) / r^3);
    v(i + 1) = v(i) - dt * G * M * (y(i) / r^3);

    % update position
    x(i + 1) = x(i) + dt*u(i+1);
    y(i + 1) = y(i) + dt*v(i+1);

    % clear the current frame
    clf

    % plot the current position
    hold on
    plot(0,0,'r*')
    plot(x(1:i), y(1:i))
    plot(x(i), y(i), 'g*')
    hold off

    % set axis limits so they don't change through movie frames
    axis equal
    xlim([-1.5*x(1),1.5*x(1)])
    ylim([-1.5*x(1),1.5*x(1)])
    
    % write this figure to the video as a new frame
    currentframe = getframe();
    writeVideo(writerObj, currentframe)
end

% close video (saves video)
close(writerObj)