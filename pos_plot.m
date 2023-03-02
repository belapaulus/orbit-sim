function pos_plot(pos)
N = size(pos, 1);
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
