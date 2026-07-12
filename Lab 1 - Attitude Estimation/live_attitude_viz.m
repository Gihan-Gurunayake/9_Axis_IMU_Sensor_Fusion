%% Real-time roll/pitch estimation from HyperIMU (UDP) with live 3D view
% Stop the program by closing the figure window.
%
% --------------------------------------------------------------------
%  Pipeline per packet:  parse -> measure dt -> predict -> update
%  Animation is redrawn (throttled) once per outer loop pass.
% --------------------------------------------------------------------

%% ----------------------------- USER CONFIG -----------------------------
PORT = 5555;   % must match HyperIMU's UDP target port

% --- CSV column layout of each datagram ---
% EDIT to match the sensor order you enabled in HyperIMU. Default assumes
% the line is:  ax, ay, az, gx, gy, gz, [anything after is ignored]
COL_AX = 1; COL_AY = 2; COL_AZ = 3;     % accelerometer columns
COL_GX = 4; COL_GY = 5; COL_GZ = 6;     % gyroscope columns

% --- Gyro units as they arrive from HyperIMU ---
% Android's "Gyroscope" sensor outputs rad/s natively -> leave this FALSE.
% Set TRUE only if your particular stream is in deg/s.
GYRO_IN_DEG = false;

% --- Optional axis sign flips (set to -1 to invert after the static test) ---
SGN_GX = +1; SGN_GY = +1; SGN_GZ = +1;

% --- Kalman tuning (same spirit as the offline script) ---
Q_DIAG = 0.01;               % process noise on [phi bphi theta btheta]
R_DIAG = 10;                 % measurement noise on [phi_acc theta_acc]

% --- dt sanity limits (seconds) ---
DT_MIN = 1e-4;               % floor (avoids divide-by-tiny in transients)
DT_MAX = 0.10;               % packets spaced wider than this = treated as a gap

%% ------------------------------- KF INIT -------------------------------
P = eye(4);
Q = eye(4) * Q_DIAG;
R = eye(2) * R_DIAG;
x = [0; 0; 0; 0];            % state: [phi; bias_phi; theta; bias_theta] (rad)
C = [1 0 0 0; 0 0 1 0];

haveInit = false;            % seed state from the first accel sample
tPrev    = [];
phi      = 0;
theta    = 0;
sensor_dataset = [];

%% --------------------------------- UDP ---------------------------------
clear u
u = udpport("datagram", "LocalPort", PORT);
flush(u);

%% ------------------------------ 3D FIGURE ------------------------------
fig = figure('Name','Live Attitude (roll/pitch)','NumberTitle','off');
ax3 = axes('Parent',fig); hold(ax3,'on');
view(ax3, 135, 25);  axis(ax3,'equal');  grid(ax3,'on');  rotate3d(ax3,'on');
xlim(ax3,[-1.2 1.2]); ylim(ax3,[-1.2 1.2]); zlim(ax3,[-1.2 1.2]);
xlabel(ax3,'X'); ylabel(ax3,'Y'); zlabel(ax3,'Z');

% phone-like slab (half-extents): narrow in X, long in Y, thin in Z
hx = 0.35; hy = 0.70; hz = 0.05;
V = [-hx -hy -hz;  hx -hy -hz;  hx  hy -hz; -hx  hy -hz; ...
     -hx -hy  hz;  hx -hy  hz;  hx  hy  hz; -hx  hy  hz];
Fc = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

tf = hgtransform('Parent',ax3);
patch('Parent',tf,'Vertices',V,'Faces',Fc, ...
      'FaceColor',[0.20 0.50 0.90],'FaceAlpha',0.92,'EdgeColor','k');
% highlight the +Z ("screen") face so orientation is readable at a glance
patch('Parent',tf,'Vertices',V,'Faces',[5 6 7 8], ...
      'FaceColor',[0.95 0.85 0.20],'EdgeColor','k');

% fixed world reference axes
quiver3(ax3,0,0,0, 1,0,0,'r','LineWidth',1.5,'MaxHeadSize',0.5);
quiver3(ax3,0,0,0, 0,1,0,'g','LineWidth',1.5,'MaxHeadSize',0.5);
quiver3(ax3,0,0,0, 0,0,1,'b','LineWidth',1.5,'MaxHeadSize',0.5);
ttl = title(ax3,'roll = ---   pitch = ---   (waiting for packets)');

%% ------------------------------ LIVE LOOP ------------------------------
tStart = tic;
needCols = [COL_AX COL_AY COL_AZ COL_GX COL_GY COL_GZ];

while ishandle(fig)
    got = false;

    % ---- drain everything queued; filter EVERY packet (keeps gyro integ.) ----
    while u.NumDatagramsAvailable > 0
        tNow = toc(tStart);
        dg   = read(u, 1, "string");
        vals = str2double(split(string(dg.Data), ","));
        %sensor_dataset = [sensor_dataset vals];

        if numel(vals) < max(needCols) || any(isnan(vals(needCols)))
            continue                          % malformed / short line
        end

        axm = vals(COL_AX); aym = vals(COL_AY); azm = vals(COL_AZ);
        p = SGN_GX * vals(COL_GX);
        q = SGN_GY * vals(COL_GY);
        r = SGN_GZ * vals(COL_GZ);
        if GYRO_IN_DEG
            p = p*pi/180;  q = q*pi/180;  r = r*pi/180;
        end

        % accelerometer-derived angles (ratio-based -> accel units cancel)
        phi_acc   = atan2( aym, sqrt(axm^2 + azm^2));
        theta_acc = atan2(-axm, sqrt(aym^2 + azm^2));

        % --- seed filter from the first valid sample, then start integrating ---
        if ~haveInit
            x(1) = phi_acc;  x(3) = theta_acc;
            phi  = x(1);     theta = x(3);
            tPrev = tNow;  haveInit = true;  got = true;
            continue
        end

        dt = tNow - tPrev;  tPrev = tNow;
        dt = min(max(dt, DT_MIN), DT_MAX);    % clamp jitter / gaps

        % --- rebuild dt-dependent matrices each step (variable sampling) ---
        A = [1 -dt 0 0; 0 1 0 0; 0 0 1 -dt; 0 0 0 1];
        B = [dt 0; 0 0; 0 dt; 0 0];

        phi_hat = x(1);  theta_hat = x(3);
        phi_dot   = p + sin(phi_hat)*tan(theta_hat)*q + cos(phi_hat)*tan(theta_hat)*r;
        theta_dot = cos(phi_hat)*q - sin(phi_hat)*r;

        % Predict
        x = A*x + B*[phi_dot; theta_dot];
        P = A*P*A' + Q;

        % Update
        z = [phi_acc; theta_acc];
        yk = z - C*x;
        S  = R + C*P*C';
        K  = P*C'/S;     % mrdivide: stabler than inv(S)
        x  = x + K*yk;
        P  = (eye(4) - K*C)*P;

        phi = x(1);  theta = x(3);
        got = true;
    end

    % ---- throttled redraw ----
    if got
        % world-from-body with zero yaw (ZYX): R = Ry(theta) * Rx(phi)
        M = makehgtform('yrotate', theta) * makehgtform('xrotate', phi);
        set(tf, 'Matrix', M);
        set(ttl, 'String', sprintf('roll = %+6.1f\\circ   pitch = %+6.1f\\circ', ...
            phi*180/pi, theta*180/pi));
        drawnow limitrate
    else
        pause(0.005);
    end
end

clear u
disp("Stopped.");