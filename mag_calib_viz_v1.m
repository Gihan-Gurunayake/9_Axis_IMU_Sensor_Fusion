
%% Real-time ROLL/PITCH/YAW estimation from HyperIMU (UDP) with live 3D view
% Extension of the roll/pitch-only lab script:
%   SECTION A : 30 s magnetometer calibration (live ellipsoid -> A, b -> sphere)
%   SECTION B : 6-state Kalman filter [phi bphi theta btheta psi bpsi]
%               yaw measured from the CALIBRATED magnetometer, tilt-compensated
%               with the filter's own roll/pitch states. Yaw innovation is
%               wrapped to the shortest arc (required at the +/-180 deg seam).
%
% HyperIMU datagram order (per user setup):  ACCEL , MAGNETOMETER , GYRO
%   -> ax ay az mx my mz gx gy gz
% Stop the program by closing the live-attitude figure window.

%% ----------------------------- USER CONFIG -----------------------------
clear all;
PORT = 5555;                 % must match HyperIMU's UDP target port

% --- CSV column layout: accel(1-3), MAG(4-6), GYRO(7-9). NOTE: gyro moved
%     from 4-6 (old script) to 7-9 because the magnetometer sits between.
COL_AX = 1; COL_AY = 2; COL_AZ = 3;
COL_MX = 4; COL_MY = 5; COL_MZ = 6;     % magnetometer (uT on Android)
COL_GX = 7; COL_GY = 8; COL_GZ = 9;     % gyroscope    (rad/s on Android)

GYRO_IN_DEG = false;         % Android gyro is rad/s natively
SGN_GX = +1; SGN_GY = +1; SGN_GZ = +1;

% --- AXIS / SIGN ALIGNMENT (this is what fixes "orientation looks altered") ---
% All three sensors must express the SAME right-handed body frame that the
% roll/pitch/yaw math assumes: X = phone right, Y = phone top, Z = out of screen.
% HyperIMU streams raw per-sensor values; the magnetometer in particular is
% frequently permuted or sign-flipped relative to accel/gyro on a given device.
%
% MAG_MAP picks which raw mag column becomes body [X Y Z]; MAG_SGN flips signs.
% Identity = [1 2 3] / [+1 +1 +1]. See the STATIC-TEST guide at the bottom to
% determine the correct values for YOUR phone, then edit these two lines only.
MAG_MAP = [1 2 3];           % e.g. [2 1 3] swaps mag X<->Y
MAG_SGN = [-1 -1 +1];        % e.g. [+1 +1 -1] flips mag Z

% Accelerometer sign flips (leave +1 unless the static test says otherwise).
SGN_AX = +1; SGN_AY = +1; SGN_AZ = +1;

% Overall yaw sense: set -1 if yaw turns the WRONG way after axes are correct.
YAW_SENSE = -1;

% --- Magnetometer calibration ---
CAL_DURATION = 120;           % seconds of data collection
CAL_MIN_SAMPLES = 300;       % refuse to fit on less than this
CAL_REFIT_EVERY = 0.5;       % seconds between live ellipsoid re-fits

% --- Kalman tuning ---
Q_DIAG = 0.01;               % process noise, all 6 states (retune psi pair if needed)
R_TILT = 10;                 % accel roll/pitch measurement variance (as before)
R_PSI  = 10;                 % mag yaw measurement variance (TUNE: raise if yaw is twitchy)

% --- Heading reference ---
D_DECL = 0;                  % magnetic declination [rad]. Colombo ~ -1.91 deg
                             % (= -0.0333 rad) if psi must reference TRUE north.
                             % Leave 0 to reference magnetic north (recommended).

% --- Yaw measurement gate: skip mag update when field is disturbed ---
MAG_NORM_TOL = 0.25;         % accept ||m_cal|| within [1-tol, 1+tol] (H = 1)

% --- dt sanity limits (seconds) ---
DT_MIN = 1e-4;  DT_MAX = 0.10;

%% =================== SECTION A: MAGNETOMETER CALIBRATION ===============
% Rotate the phone SLOWLY through ALL orientations (figure-eights + full
% rolls). Coverage matters more than time: a lazy wave gives a degenerate
% ellipsoid and the fit will be rejected. Keep the USB cable slack away
% from the phone body while doing this.

clear u
u = udpport("datagram", "LocalPort", PORT);
flush(u);

figCal = figure('Name','Magnetometer calibration','NumberTitle','off', ...
                'Position',[80 80 1150 520]);
axRaw = subplot(1,2,1); hold(axRaw,'on'); grid(axRaw,'on'); axis(axRaw,'equal');
view(axRaw,135,25); rotate3d(axRaw,'on');
xlabel(axRaw,'m_x [\muT]'); ylabel(axRaw,'m_y [\muT]'); zlabel(axRaw,'m_z [\muT]');
hPts  = plot3(axRaw,nan,nan,nan,'.','MarkerSize',4,'Color',[0.25 0.45 0.85]);
hEll  = surf(axRaw,nan(2),nan(2),nan(2),'FaceColor',[0.9 0.4 0.2], ...
             'FaceAlpha',0.25,'EdgeColor','none');
tRaw  = title(axRaw,'collecting...');

axCal = subplot(1,2,2); hold(axCal,'on'); grid(axCal,'on'); axis(axCal,'equal');
view(axCal,135,25); rotate3d(axCal,'on');
xlabel(axCal,'x'); ylabel(axCal,'y'); zlabel(axCal,'z');
title(axCal,'calibrated (after fit)');
[sx,sy,sz] = sphere(24);
surf(axCal,sx,sy,sz,'FaceColor','none','EdgeColor',[0.7 0.7 0.7]); % unit sphere

magBuf   = zeros(0,3);
tCal     = tic;  tLastFit = -inf;
bHat = []; Mt = [];                       % fit results (raw uT coordinates)

while toc(tCal) < CAL_DURATION && ishandle(figCal)
    while u.NumDatagramsAvailable > 0
        dg   = read(u, 1, "string");
        vals = str2double(split(string(dg.Data), ","));
        if numel(vals) < COL_MZ || any(isnan(vals(COL_MX:COL_MZ))), continue, end
        mraw = vals([COL_MX COL_MY COL_MZ]);
        mbody = MAG_SGN(:) .* mraw(MAG_MAP);      % same body frame as accel/gyro
        magBuf(end+1,:) = mbody.'; %#ok<SAGROW>   % calibrate in the CORRECT frame
    end

    % ---- live ellipsoid re-fit (throttled) ----
    if size(magBuf,1) >= 100 && toc(tCal) - tLastFit > CAL_REFIT_EVERY
        tLastFit = toc(tCal);
        [bTry, MtTry, ok] = ellipsoidFitLocal(magBuf);
        if ok
            bHat = bTry;  Mt = MtTry;
            [Xe,Ye,Ze] = ellipsoidSurfLocal(bHat, Mt, sx, sy, sz);
            set(hEll,'XData',Xe,'YData',Ye,'ZData',Ze);
        end
    end

    set(hPts,'XData',magBuf(:,1),'YData',magBuf(:,2),'ZData',magBuf(:,3));
    set(tRaw,'String',sprintf('raw + fitted ellipsoid   |   N = %d   |   %4.1f s left', ...
        size(magBuf,1), max(0, CAL_DURATION - toc(tCal))));
    drawnow limitrate
end

% ---- final fit + hard validity gate ----
assert(size(magBuf,1) >= CAL_MIN_SAMPLES, ...
    'Calibration aborted: only %d samples (< %d). Check stream / columns.', ...
    size(magBuf,1), CAL_MIN_SAMPLES);
[bHat, Mt, ok] = ellipsoidFitLocal(magBuf);
assert(ok, ['Calibration aborted: data do not define an ellipsoid ', ...
    '(insufficient orientation coverage). Redo with slow full rotations.']);

% ---- recover A (H = 1 convention: calibrated field has unit norm) ----
[Vm, Lm] = eig((Mt+Mt.')/2);              % Mt = M/H^2, symmetric PD
lam      = diag(Lm);
A_cal    = Vm * diag(1./sqrt(lam)) * Vm.'; % A  = V L^(-1/2) V'
A_inv    = Vm * diag(   sqrt(lam)) * Vm.'; % A^-1 (exact closed form, no inv())

fprintf('\n================ MAGNETOMETER CALIBRATION RESULT ================\n');
fprintf('Hard-iron bias  b  [uT] :  [% .3f  % .3f  % .3f]\n', bHat);
fprintf('Soft-iron matrix A (H = 1 normalisation):\n'); disp(A_cal);
fprintf('A is the symmetric square-root representative (A recoverable only\n');
fprintf('up to a right orthogonal factor; any rotational soft-iron component\n');
fprintf('appears as a constant yaw offset).\n');

% ---- calibrated cloud on the unit sphere (adjacent plot) ----
mc = (A_inv * (magBuf - bHat.').').';
plot3(axCal, mc(:,1), mc(:,2), mc(:,3), '.', 'MarkerSize',4, ...
      'Color',[0.20 0.60 0.30]);
nrm = vecnorm(mc,2,2);
title(axCal, sprintf('calibrated: ||m_{cal}|| = %.3f \\pm %.3f', mean(nrm), std(nrm)));
set(tRaw,'String',sprintf('DONE  |  N = %d  |  b = [%.2f %.2f %.2f] \\muT', ...
    size(magBuf,1), bHat));
drawnow

%% ------------------------------- KF INIT -------------------------------
% State: [phi; bias_phi; theta; bias_theta; psi; bias_psi]   (rad)
P = eye(6);
Q = eye(6) * Q_DIAG;
x = zeros(6,1);
C3 = [1 0 0 0 0 0;                       % accel -> phi
      0 0 1 0 0 0;                       % accel -> theta
      0 0 0 0 1 0];                      % mag   -> psi ONLY (decoupling row)
R3 = diag([R_TILT R_TILT R_PSI]);
C2 = C3(1:2,:);  R2 = R3(1:2,1:2);       % tilt-only update (mag gated out)

haveInit = false;
tPrev    = [];
phi = 0; theta = 0; psi = 0;

flush(u);                                % reuse the same UDP port

% ------------------------------ 3D FIGURE ------------------------------
fig = figure('Name','Live Attitude (roll/pitch/yaw)','NumberTitle','off');
ax3 = axes('Parent',fig); hold(ax3,'on');
view(ax3, 135, 25);  axis(ax3,'equal');  grid(ax3,'on');  rotate3d(ax3,'on');
xlim(ax3,[-1.2 1.2]); ylim(ax3,[-1.2 1.2]); zlim(ax3,[-1.2 1.2]);
xlabel(ax3,'X'); ylabel(ax3,'Y'); zlabel(ax3,'Z');

hx = 0.35; hy = 0.70; hz = 0.05;
V = [-hx -hy -hz;  hx -hy -hz;  hx  hy -hz; -hx  hy -hz; ...
     -hx -hy  hz;  hx -hy  hz;  hx  hy  hz; -hx  hy  hz];
Fc = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

tf = hgtransform('Parent',ax3);
patch('Parent',tf,'Vertices',V,'Faces',Fc, ...
      'FaceColor',[0.20 0.50 0.90],'FaceAlpha',0.92,'EdgeColor','k');
patch('Parent',tf,'Vertices',V,'Faces',[5 6 7 8], ...
      'FaceColor',[0.95 0.85 0.20],'EdgeColor','k');

quiver3(ax3,0,0,0, 1,0,0,'r','LineWidth',1.5,'MaxHeadSize',0.5);
quiver3(ax3,0,0,0, 0,1,0,'g','LineWidth',1.5,'MaxHeadSize',0.5);
quiver3(ax3,0,0,0, 0,0,1,'b','LineWidth',1.5,'MaxHeadSize',0.5);
ttl = title(ax3,'roll = ---  pitch = ---  yaw = ---  (waiting)');

% ------------------------------ LIVE LOOP ------------------------------
tStart   = tic;
needCols = [COL_AX COL_AY COL_AZ COL_MX COL_MY COL_MZ COL_GX COL_GY COL_GZ];

while ishandle(fig)
    got = false;

    while u.NumDatagramsAvailable > 0
        tNow = toc(tStart);
        dg   = read(u, 1, "string");
        vals = str2double(split(string(dg.Data), ","));
        if numel(vals) < max(needCols) || any(isnan(vals(needCols))), continue, end

        axm = SGN_AX * vals(COL_AX); aym = SGN_AY * vals(COL_AY); azm = SGN_AZ * vals(COL_AZ);
        p = SGN_GX * vals(COL_GX);
        q = SGN_GY * vals(COL_GY);
        r = SGN_GZ * vals(COL_GZ);
        if GYRO_IN_DEG, p = p*pi/180; q = q*pi/180; r = r*pi/180; end

        % --- magnetometer: SAME remap as calibration, THEN calibrate ---
        mraw  = vals([COL_MX COL_MY COL_MZ]);
        mbody = MAG_SGN(:) .* mraw(MAG_MAP);
        m_cal = A_inv * (mbody - bHat);                        % ~unit norm

        % --- accelerometer tilt angles (unchanged) ---
        phi_acc   = atan2( aym, sqrt(axm^2 + azm^2));
        theta_acc = atan2(-axm, sqrt(aym^2 + azm^2));

        % --- seed all three angles from the first valid sample ---
        if ~haveInit
            x(1) = phi_acc;  x(3) = theta_acc;
            x(5) = magYawLocal(m_cal, x(1), x(3)) + D_DECL;
            phi = x(1); theta = x(3); psi = x(5);
            tPrev = tNow;  haveInit = true;  got = true;
            continue
        end

        dt = tNow - tPrev;  tPrev = tNow;
        dt = min(max(dt, DT_MIN), DT_MAX);

        % --- dt-dependent matrices: old 4x4/4x2 blocks + appended yaw block ---
        A6 = [1 -dt 0   0  0   0;
              0   1 0   0  0   0;
              0   0 1 -dt  0   0;
              0   0 0   1  0   0;
              0   0 0   0  1 -dt;
              0   0 0   0  0   1];
        B6 = [dt 0 0; 0 0 0; 0 dt 0; 0 0 0; 0 0 dt; 0 0 0];

        % --- Euler rates from CURRENT estimates (yaw row now restored) ---
        phi_hat = x(1);  theta_hat = x(3);
        cth = cos(theta_hat);
        cth = sign(cth + (cth==0)) * max(abs(cth), 0.05);   % gimbal guard
        phi_dot   = p + sin(phi_hat)*tan(theta_hat)*q + cos(phi_hat)*tan(theta_hat)*r;
        theta_dot = cos(phi_hat)*q - sin(phi_hat)*r;
        psi_dot   = (sin(phi_hat)*q + cos(phi_hat)*r) / cth;

        % Predict
        x = A6*x + B6*[phi_dot; theta_dot; psi_dot];
        P = A6*P*A6' + Q;

        % --- yaw measurement: tilt-compensate m_cal with PREDICTED phi/theta ---
        psi_mag = magYawLocal(m_cal, x(1), x(3)) + D_DECL;
        magOK   = abs(norm(m_cal) - 1) < MAG_NORM_TOL;   % disturbance gate

        % Update
        if magOK
            z  = [phi_acc; theta_acc; psi_mag];
            yk = z - C3*x;
            yk(3) = atan2(sin(yk(3)), cos(yk(3)));       % REQUIRED wrap (+/-180 seam)
            S = R3 + C3*P*C3';
            K = P*C3'/S;
            x = x + K*yk;
            P = (eye(6) - K*C3)*P;
        else                                             % field disturbed: tilt only
            z  = [phi_acc; theta_acc];
            yk = z - C2*x;
            S = R2 + C2*P*C2';
            K = P*C2'/S;
            x = x + K*yk;
            P = (eye(6) - K*C2)*P;
        end
        x(5) = atan2(sin(x(5)), cos(x(5)));              % keep psi in (-pi, pi]

        phi = x(1);  theta = x(3);  psi = x(5);
        got = true;
    end

    if got
        % world-from-body ZYX: prepend the yaw rotation to the old transform
        M = makehgtform('zrotate', YAW_SENSE*psi) * makehgtform('yrotate', theta) ...
          * makehgtform('xrotate', phi);
        set(tf, 'Matrix', M);
        set(ttl, 'String', sprintf( ...
            'roll = %+6.1f\\circ  pitch = %+6.1f\\circ  yaw = %+6.1f\\circ', ...
            phi*180/pi, theta*180/pi, psi*180/pi));
        drawnow limitrate
    else
        pause(0.005);
    end

    fprintf('m=[% .1f % .1f % .1f]  ax=[% .1f % .1f % .1f]\n', mbody, axm, aym, azm);
end

clear u
disp("Stopped.");

%% ---------------------------- LOCAL FUNCTIONS ---------------------------
function [b, Mt, ok] = ellipsoidFitLocal(Mdat)
% Constrained LS ellipsoid fit (SVD null-vector), per the derived algorithm:
% design row d(m) = [x^2 y^2 z^2 2xy 2xz 2yz 2x 2y 2z 1];  theta = V(:,end).
% Pre-centred / pre-scaled for conditioning (uT^2 vs 1 spans ~4 decades).
    mu = mean(Mdat,1);
    sc = mean(vecnorm(Mdat - mu, 2, 2));
    n  = (Mdat - mu) / sc;
    xs = n(:,1); ys = n(:,2); zs = n(:,3);
    D  = [xs.^2 ys.^2 zs.^2 2*xs.*ys 2*xs.*zs 2*ys.*zs 2*xs 2*ys 2*zs ones(size(xs))];
    [~,~,Vd] = svd(D, 'econ');
    th = Vd(:,end);
    Qm = [th(1) th(4) th(5); th(4) th(2) th(6); th(5) th(6) th(3)];
    uv = th(7:9);  s = th(10);
    ev = eig((Qm+Qm.')/2);
    if all(ev < 0), Qm = -Qm; uv = -uv; s = -s; ev = -ev; end
    ok = all(ev > 0);                       % mixed signs => not an ellipsoid
    if ~ok, b = zeros(3,1); Mt = eye(3); return, end
    bn  = -(Qm \ uv);
    Mn  = Qm / (uv.'*(Qm\uv) - s);          % = M/H^2 in normalised coords
    ok  = all(eig((Mn+Mn.')/2) > 0);
    b   = mu.' + sc*bn;                     % undo pre-conditioning
    Mt  = Mn / sc^2;                        % raw-coordinate M/H^2
end

function [Xe,Ye,Ze] = ellipsoidSurfLocal(b, Mt, sx, sy, sz)
% Map the unit sphere through the ellipsoid (m-b)' Mt (m-b) = 1 for display.
    [Vm,Lm] = eig((Mt+Mt.')/2);
    T = Vm * diag(1./sqrt(diag(Lm)));       % semi-axes = 1/sqrt(eig)
    Pts = T * [sx(:) sy(:) sz(:)].' + b;
    Xe = reshape(Pts(1,:), size(sx));
    Ye = reshape(Pts(2,:), size(sy));
    Ze = reshape(Pts(3,:), size(sz));
end

function psi_m = magYawLocal(m, phi, theta)
% Tilt-compensated heading, EXACT for the lab convention
% R_f^b = Rx(phi) Ry(theta) Rz(psi)  (verified numerically to machine precision).
    mxh = m(1)*cos(theta) + m(2)*sin(phi)*sin(theta) - m(3)*cos(phi)*sin(theta);
    myh = m(2)*cos(phi)   + m(3)*sin(phi);
    psi_m = atan2(myh, mxh);
end

% ------------------------- STATIC-TEST GUIDE (read) ---------------------
% If orientation looks altered, DO NOT random-flip signs. Identify the axis
% frame first with these one-time tests, then set the config at the top.
%
% Add a temporary line inside the live loop to print the corrected triples:
%   fprintf('a=[% .2f % .2f % .2f]  m=[% .2f % .2f % .2f]\n', axm,aym,azm, mbody);
%
% ACCEL (defines roll/pitch). Hold the phone and watch the a=[...] print:
%   flat, screen up   -> az strongly POSITIVE (~+9.8), ax~0, ay~0
%   tilt top edge up  -> ax should go NEGATIVE (that is +pitch here)
%   tilt right edge up-> ay should go POSITIVE  (that is +roll here)
%   If a sign is opposite, flip the matching SGN_AX/AY/AZ.
%
% GYRO (defines how fast angles move). Rotate about each phone axis:
%   the on-screen angle must move the SAME way as the phone. If one axis
%   moves backward, flip SGN_GX/GY/GZ for that axis.
%
% MAGNETOMETER (defines yaw). With axes above already correct:
%   Lay the phone flat, rotate it slowly clockwise (viewed from above).
%   Watch m=[...]: as you spin, mx and my should trace a circle. Point the
%   phone's +Y (top) toward magnetic north -> the horizontal mag vector should
%   lie mostly along +body-X after tilt-comp. If yaw is rotated 90 deg or
%   mirrored, the mag axes are permuted/flipped: set MAG_MAP / MAG_SGN so the
%   mag frame matches the accel frame (same X, Y, Z directions).
%   Quick check: the calibrated cloud must be a SPHERE (it will be regardless),
%   but yaw only reads correctly when MAG axes == ACCEL axes.
%
% YAW_SENSE: after all axes are correct, if the slab still yaws the wrong way
% on screen (pure mirror), set YAW_SENSE = -1. This is display-only; it does
% not touch the filter. Never use it to mask an axis mismatch.
%
% ISOLATION RULE: if the calibration ellipsoid looked clean and the sphere is
% tight (||m_cal|| ~ 1.00), the problem is NOT calibration -- it is axis/sign
% alignment above. Fix it there, not by retuning Q/R.