%% Connection Test
% Rows = sensors (Accel / Mag / Gyro), Cols = axes (X / Y / Z).
% Sliding time window (last WIN_SEC seconds) with animatedline ring buffers.
% Stop by closing the figure window.
%
% HyperIMU datagram order (your setup):  ACCEL , MAG , GYRO
%   ax ay az mx my mz gx gy gz
% No filtering / no calibration - this is raw-signal inspection only.

%% ----------------------------- USER CONFIG -----------------------------
PORT = 5555;                 % must match HyperIMU's UDP target port

COL_AX = 1; COL_AY = 2; COL_AZ = 3;     % accelerometer [m/s^2]
COL_MX = 4; COL_MY = 5; COL_MZ = 6;     % magnetometer  [uT]
COL_GX = 7; COL_GY = 8; COL_GZ = 9;     % gyroscope     [rad/s]

WIN_SEC   = 10;              % width of the scrolling time window [s]
DT_MAX    = 0.10;            % gap guard (unused for plotting, kept for parity)

% --------------------------------- UDP ---------------------------------
clear u
u = udpport("datagram", "LocalPort", PORT);
flush(u);

% ------------------------------ FIGURE ---------------------------------
fig = figure('Name','HyperIMU - all sensors (raw)','NumberTitle','off', ...
             'Position',[60 60 1200 760]);
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
title(tl,'Live raw sensor channels (last 10 s)');

sensorNames = {'Accel [m/s^2]','Mag [\muT]','Gyro [rad/s]'};
axisNames   = {'X','Y','Z'};
cols        = {COL_AX COL_AY COL_AZ; COL_MX COL_MY COL_MZ; COL_GX COL_GY COL_GZ};
lineClr     = [0.85 0.20 0.20; 0.20 0.65 0.25; 0.20 0.40 0.85];  % X/Y/Z colours

axArr = gobjects(3,3);
lnArr = gobjects(3,3);
colIdx = zeros(3,3);
for s = 1:3                         % sensor row
    for a = 1:3                     % axis column
        ax = nexttile(tl);  hold(ax,'on');  grid(ax,'on');
        lnArr(s,a) = animatedline(ax,'Color',lineClr(a,:),'LineWidth',1.0);
        title(ax,sprintf('%s  %s', sensorNames{s}, axisNames{a}));
        if s == 3, xlabel(ax,'t [s]'); end
        xlim(ax,[0 WIN_SEC]);
        axArr(s,a) = ax;
        colIdx(s,a) = cols{s,a};
    end
end

% ------------------------------ LIVE LOOP ------------------------------
tStart   = tic;
needCols = [COL_AX COL_AY COL_AZ COL_MX COL_MY COL_MZ COL_GX COL_GY COL_GZ];

while ishandle(fig)
    got = false;

    % drain everything queued; append every valid packet to the ring lines
    while u.NumDatagramsAvailable > 0
        tNow = toc(tStart);
        dg   = read(u, 1, "string");
        vals = str2double(split(string(dg.Data), ","));
        if numel(vals) < max(needCols) || any(isnan(vals(needCols)))
            continue                                  % malformed / short line
        end
        for s = 1:3
            for a = 1:3
                addpoints(lnArr(s,a), tNow, vals(colIdx(s,a)));
            end
        end
        got = true;
    end

    % scroll all tiles to the last WIN_SEC seconds, then redraw once
    if got
        tNow = toc(tStart);
        x0 = max(0, tNow - WIN_SEC);
        for s = 1:3
            for a = 1:3
                xlim(axArr(s,a), [x0, x0 + WIN_SEC]);
            end
        end
        drawnow limitrate
    else
        pause(0.005);
    end
end

clear u
disp("Stopped.");




