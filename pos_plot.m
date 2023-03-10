function pos_plot(pos)
% takes an array of positions of the from array(body, timestep, component)
% as input
% produces a plot of the positions
N = size(pos, 1);
hold on
for i = 1:N
	    plot3( ...
            squeeze(pos(i, :, 1)), ...
            squeeze(pos(i, :, 2)), ...
            squeeze(pos(i, :, 3)), ...
            '.'...
        )
        axis equal
end
end

