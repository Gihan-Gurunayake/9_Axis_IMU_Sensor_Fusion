%% Real-time quaternion MEKF attitude estimation from HyperIMU (UDP), live 3D view
% Fully coupled MARG fusion per the derivation PDF:
%   SECTION A : 30 s magnetometer calibration (live ellipsoid -> A, b -> sphere)
%   SECTION B : STILL initialisation (~1.5 s): gyro bias b0, roll/pitch from
%               accel, and the magnetic reference m_ref MEASURED from data
%               (absorbs declination + device sign convention; yaw therefore
%               reads RELATIVE to the starting pose, psi0 = 0).
%   SECTION C : MEKF loop. Nominal state = unit quaternion (body-to-nav,
%               Hamilton scalar-first) + gyro bias. 6x6 error covariance.
%               Accel and mag enter as FULL 3-axis vector measurements with
%               H = [ [h_hat]x | 0 ]. Multiplicative injection, Joseph form.
%
% HyperIMU datagram order:  ACCEL , MAG , GYRO  ->  ax ay az mx my mz gx gy gz
% Stop by closing the live-attitude figure.
%
% OPERATIONAL REQUIREMENTS (new vs the Euler version):
%   * Phone must be HELD STILL during the INIT window, or b0 and m_ref are
%     silently wrong and the filter starts poisoned.
%   * MAG_MAP / MAG_SGN must be the TRUE device axis alignment. The coupled
%     filter consumes the mag as a vector: a wrong map now corrupts roll and
%     pitch near disturbances, not just yaw.

%% ----------------------------- USER CONFIG -----------------------------
clear all;
PORT = 5555;

COL_AX = 1; COL_AY = 2; COL_AZ = 3;      % accel [m/s^2]
COL_MX = 4; COL_MY = 5; COL_MZ = 6;      % mag   [uT]
COL_GX = 7; COL_GY = 8; COL_GZ = 9;      % gyro  [rad/s]

GYRO_IN_DEG = false;
SGN_AX = +1; SGN_AY = +1; SGN_AZ = +1;
SGN_GX = +1; SGN_GY = +1; SGN_GZ = +1;
MAG_MAP = [1 2 3];
MAG_SGN = [-1 -1 +1];        % your device's verified axis alignment

% --- Magnetometer calibration ---
CAL_DURATION = 90;  CAL_MIN_SAMPLES = 300;  CAL_REFIT_EVERY = 0.5;

% --- Still-init window ---
INIT_TIME = 1.5;             % seconds; HOLD THE PHONE STILL

% --- MEKF noise parameters (TUNING KNOBS, not derived values) ---
GYRO_NOISE_STD = 0.01;       % rad/s          -> Q attitude block
BIAS_WALK_STD  = 5e-4;       % rad/s/sqrt(s)  -> Q bias block
ACC_NOISE_STD  = 0.5;        % m/s^2 (inflated: handheld motion is model error)
MAG_NOISE_STD  = 0.05;       % unit-field units

% --- Initial uncertainty ---
P0_ATT  = deg2rad(10)^2;     % rad^2 per attitude axis
P0_BIAS = 0.02^2;            % (rad/s)^2 per bias axis

% --- Measurement gates ---
ACC_GATE     = 1.0;          % accept accel update if | ||a|| - g | < this [m/s^2]
MAG_NORM_TOL = 0.25;         % accept mag  update if | ||m_cal|| - 1 | < this

GRAV = 9.81;
DT_MIN = 1e-4;  DT_MAX = 0.10;

%% =================== SECTION A: MAGNETOMETER CALIBRATION ===============
% Rotate the phone SLOWLY through ALL orientations. Coverage beats time.
clear u
u = udpport("datagram", "LocalPort", PORT);
flush(u);

figCal = figure('Name','Magnetometer calibration','NumberTitle','off', ...
                'Position',[80 80 1150 520]);
axRaw = subplot(1,2,1); hold(axRaw,'on'); grid(axRaw,'on'); axis(axRaw,'equal');
view(axRaw,135,25); rotate3d(axRaw,'on');
xlabel(axRaw,'m_x [\muT]'); ylabel(axRaw,'m_y [\muT]'); zlabel(axRaw,'m_z [\muT]');
hPts = plot3(axRaw,nan,nan,nan,'.','MarkerSize',4,'Color',[0.25 0.45 0.85]);
hEll = surf(axRaw,nan(2),nan(2),nan(2),'FaceColor',[0.9 0.4 0.2], ...
            'FaceAlpha',0.25,'EdgeColor','none');
tRaw = title(axRaw,'collecting...');

axCal = subplot(1,2,2); hold(axCal,'on'); grid(axCal,'on'); axis(axCal,'equal');
view(axCal,135,25); rotate3d(axCal,'on');
title(axCal,'calibrated (after fit)');
[sx,sy,sz] = sphere(24);
surf(axCal,sx,sy,sz,'FaceColor','none','EdgeColor',[0.7 0.7 0.7]);

magBuf = zeros(0,3);  tCal = tic;  tLastFit = -inf;
bHat = []; Mt = [];
while toc(tCal) < CAL_DURATION && ishandle(figCal)
    while u.NumDatagramsAvailable > 0
        dg = read(u,1,"string");
        vals = str2double(split(string(dg.Data), ","));
        if numel(vals) < COL_MZ || any(isnan(vals(COL_MX:COL_MZ))), continue, end
        mraw  = vals([COL_MX COL_MY COL_MZ]);
        mbody = MAG_SGN(:) .* mraw(MAG_MAP);
        magBuf(end+1,:) = mbody.'; %#ok<SAGROW>
    end
    if size(magBuf,1) >= 100 && toc(tCal)-tLastFit > CAL_REFIT_EVERY
        tLastFit = toc(tCal);
        [bTry,MtTry,ok] = ellipsoidFitLocal(magBuf);
        if ok
            bHat = bTry;  Mt = MtTry;
            [Xe,Ye,Ze] = ellipsoidSurfLocal(bHat,Mt,sx,sy,sz);
            set(hEll,'XData',Xe,'YData',Ye,'ZData',Ze);
        end
    end
    set(hPts,'XData',magBuf(:,1),'YData',magBuf(:,2),'ZData',magBuf(:,3));
    set(tRaw,'String',sprintf('raw + ellipsoid | N = %d | %4.1f s left', ...
        size(magBuf,1), max(0,CAL_DURATION-toc(tCal))));
    drawnow limitrate
end
assert(size(magBuf,1) >= CAL_MIN_SAMPLES, 'Too few samples (%d).', size(magBuf,1));
[bHat,Mt,ok] = ellipsoidFitLocal(magBuf);
assert(ok, 'Data do not define an ellipsoid: redo with full slow rotations.');

[Vm,Lm] = eig((Mt+Mt.')/2);  lam = diag(Lm);
A_cal = Vm*diag(1./sqrt(lam))*Vm.';
A_inv = Vm*diag(   sqrt(lam))*Vm.';
fprintf('\n========== MAG CALIBRATION ==========\n');
fprintf('b [uT] : [% .3f % .3f % .3f]\n', bHat);
fprintf('A (H=1):\n'); disp(A_cal);

mc = (A_inv*(magBuf - bHat.').').';
plot3(axCal,mc(:,1),mc(:,2),mc(:,3),'.','MarkerSize',4,'Color',[0.20 0.60 0.30]);
nrm = vecnorm(mc,2,2);
title(axCal,sprintf('calibrated: ||m|| = %.3f \\pm %.3f',mean(nrm),std(nrm)));
drawnow

%% ============== SECTION B: STILL INITIALISATION (~1.5 s) ================
fprintf('\nHOLD THE PHONE STILL for %.1f s (measuring gyro bias + m_ref)...\n', INIT_TIME);
flush(u);
accA = zeros(0,3); magA = zeros(0,3); gyrA = zeros(0,3);
tInit = tic;
needCols = [COL_AX COL_AY COL_AZ COL_MX COL_MY COL_MZ COL_GX COL_GY COL_GZ];
while toc(tInit) < INIT_TIME
    while u.NumDatagramsAvailable > 0
        dg = read(u,1,"string");
        vals = str2double(split(string(dg.Data), ","));
        if numel(vals) < max(needCols) || any(isnan(vals(needCols))), continue, end
        a = [SGN_AX*vals(COL_AX); SGN_AY*vals(COL_AY); SGN_AZ*vals(COL_AZ)];
        w = [SGN_GX*vals(COL_GX); SGN_GY*vals(COL_GY); SGN_GZ*vals(COL_GZ)];
        if GYRO_IN_DEG, w = w*pi/180; end
        mraw  = vals([COL_MX COL_MY COL_MZ]);
        mbody = MAG_SGN(:) .* mraw(MAG_MAP);
        accA(end+1,:) = a.';  gyrA(end+1,:) = w.'; %#ok<SAGROW>
        magA(end+1,:) = (A_inv*(mbody - bHat)).'; %#ok<SAGROW>
    end
    pause(0.005);
end
assert(size(accA,1) >= 20, 'Init window received too few packets.');
% stillness sanity check: gyro should be near zero
assert(norm(mean(gyrA,1)) < 0.15, ...
    'Init FAILED: phone was moving (mean |gyro| = %.3f rad/s). Rerun and hold still.', ...
    norm(mean(gyrA,1)));

aAvg = mean(accA,1).';   mAvg = mean(magA,1).';   bhat = mean(gyrA,1).';
phi0   = atan2( aAvg(2), sqrt(aAvg(1)^2 + aAvg(3)^2));
theta0 = atan2(-aAvg(1), sqrt(aAvg(2)^2 + aAvg(3)^2));
qhat   = eul2quatLocal(phi0, theta0, 0);              % psi0 = 0: heading is RELATIVE
m_ref  = Rquat(qhat) * (mAvg / norm(mAvg));           % measured reference (unit)
% For ABSOLUTE heading instead: set psi0 from your tilt-compensated formula
% here (one-time), THEN compute m_ref the same way.

P = blkdiag(P0_ATT*eye(3), P0_BIAS*eye(3));           % 6x6 error covariance
Ra = ACC_NOISE_STD^2 * eye(3);
Rm = MAG_NOISE_STD^2 * eye(3);
fprintf('init: roll %.1f deg, pitch %.1f deg, |b0| = %.4f rad/s, m_ref = [% .3f % .3f % .3f]\n', ...
    rad2deg(phi0), rad2deg(theta0), norm(bhat), m_ref);

%% ------------------------------ 3D FIGURE ------------------------------
fig = figure('Name','Live Attitude (MEKF, quaternion)','NumberTitle','off');
ax3 = axes('Parent',fig); hold(ax3,'on');
view(ax3,135,25); axis(ax3,'equal'); grid(ax3,'on'); rotate3d(ax3,'on');
xlim(ax3,[-1.2 1.2]); ylim(ax3,[-1.2 1.2]); zlim(ax3,[-1.2 1.2]);
xlabel(ax3,'X'); ylabel(ax3,'Y'); zlabel(ax3,'Z');

hx = 0.35; hy = 0.70; hz = 0.05;
V = [-hx -hy -hz; hx -hy -hz; hx hy -hz; -hx hy -hz; ...
     -hx -hy  hz; hx -hy  hz; hx hy  hz; -hx hy  hz];
Fc = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
tf = hgtransform('Parent',ax3);
patch('Parent',tf,'Vertices',V,'Faces',Fc, ...
      'FaceColor',[0.20 0.50 0.90],'FaceAlpha',0.92,'EdgeColor','k');
patch('Parent',tf,'Vertices',V,'Faces',[5 6 7 8], ...
      'FaceColor',[0.95 0.85 0.20],'EdgeColor','k');
quiver3(ax3,0,0,0,1,0,0,'r','LineWidth',1.5,'MaxHeadSize',0.5);
quiver3(ax3,0,0,0,0,1,0,'g','LineWidth',1.5,'MaxHeadSize',0.5);
quiver3(ax3,0,0,0,0,0,1,'b','LineWidth',1.5,'MaxHeadSize',0.5);
ttl = title(ax3,'waiting for packets');

% ============== SECTION C: MEKF LIVE LOOP ==============================
tStart = tic;  tPrev = [];
flush(u);
while ishandle(fig)
    got = false;
    while u.NumDatagramsAvailable > 0
        tNow = toc(tStart);
        dg = read(u,1,"string");
        vals = str2double(split(string(dg.Data), ","));
        if numel(vals) < max(needCols) || any(isnan(vals(needCols))), continue, end

        a = [SGN_AX*vals(COL_AX); SGN_AY*vals(COL_AY); SGN_AZ*vals(COL_AZ)];
        w = [SGN_GX*vals(COL_GX); SGN_GY*vals(COL_GY); SGN_GZ*vals(COL_GZ)];
        if GYRO_IN_DEG, w = w*pi/180; end
        mraw  = vals([COL_MX COL_MY COL_MZ]);
        m_cal = A_inv*(MAG_SGN(:).*mraw(MAG_MAP) - bHat);

        if isempty(tPrev), tPrev = tNow; continue, end
        dt = min(max(tNow - tPrev, DT_MIN), DT_MAX);  tPrev = tNow;

        % ---------------- PREDICT ----------------
        w_hat = w - bhat;
        qhat  = qmulLocal(qhat, qincLocal(w_hat, dt));   % exact increment (norm-exact)
        qhat  = qhat / norm(qhat);                       % vs floating-point drift
        Fk = [eye(3)-skewLocal(w_hat)*dt, -eye(3)*dt;
              zeros(3),                    eye(3)];
        Qk = blkdiag(GYRO_NOISE_STD^2*dt*eye(3), BIAS_WALK_STD^2*dt*eye(3));
        P  = Fk*P*Fk.' + Qk;

        % ---------------- ACCEL UPDATE (gated) ----------------
        if abs(norm(a) - GRAV) < ACC_GATE
            [qhat,bhat,P] = mekfUpdateLocal(qhat,bhat,P, a, [0;0;GRAV], Ra);
        end
        % ---------------- MAG UPDATE (gated) ----------------
        if abs(norm(m_cal) - 1) < MAG_NORM_TOL
            [qhat,bhat,P] = mekfUpdateLocal(qhat,bhat,P, m_cal/norm(m_cal), m_ref, Rm);
        end
        got = true;
    end

    if got
        Rbn = Rquat(qhat);                    % body-to-nav = world-from-body
        M = eye(4);  M(1:3,1:3) = Rbn;
        set(tf,'Matrix',M);
        [phi,theta,psi] = quat2eulLocal(qhat);
        set(ttl,'String',sprintf( ...
            'roll = %+6.1f\\circ  pitch = %+6.1f\\circ  yaw(rel) = %+6.1f\\circ  |b| = %.3f', ...
            rad2deg(phi), rad2deg(theta), rad2deg(psi), norm(bhat)));
        drawnow limitrate
    else
        pause(0.005);
    end
end
clear u
disp("Stopped.");

%% ---------------------------- LOCAL FUNCTIONS ---------------------------
function [q,b,P] = mekfUpdateLocal(q,b,P, z, ref, R)
% One full-vector MEKF update: h = R(q)' ref, H = [ [h]x | 0 ], Joseph form,
% multiplicative injection, implicit reset (error mean -> 0).
    h = Rquat(q).' * ref;
    H = [skewLocal(h), zeros(3)];
    S = H*P*H.' + R;
    K = (P*H.')/S;
    dx = K*(z - h);
    IKH = eye(6) - K*H;
    P = IKH*P*IKH.' + K*R*K.';               % Joseph form
    q = qmulLocal(q, qincLocal(dx(1:3), 1)); % exact exp-map injection
    q = q / norm(q);
    b = b + dx(4:6);
end

function q = qmulLocal(a,b)
% Hamilton product, scalar-first.
    q = [a(1)*b(1) - a(2:4).'*b(2:4);
         a(1)*b(2:4) + b(1)*a(2:4) + cross(a(2:4), b(2:4))];
end

function dq = qincLocal(v, dt)
% Unit quaternion for rotation v*dt (exp map). Safe at ||v||->0.
    ang = norm(v)*dt;
    if ang < 1e-12
        dq = [1;0;0;0];
    else
        dq = [cos(ang/2); sin(ang/2)*(v/norm(v))];
    end
end

function R = Rquat(q)
% Body-to-nav rotation matrix, Hamilton scalar-first (v^n = R q v^b).
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    R = [1-2*(q2^2+q3^2),  2*(q1*q2-q0*q3),  2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3),  1-2*(q1^2+q3^2),  2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2),  2*(q2*q3+q0*q1),  1-2*(q1^2+q2^2)];
end

function S = skewLocal(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function q = eul2quatLocal(phi, theta, psi)
% q = qz(psi) (x) qy(theta) (x) qx(phi), matching R(q)' = Rx Ry Rz (lab).
    qx = [cos(phi/2);   sin(phi/2); 0; 0];
    qy = [cos(theta/2); 0; sin(theta/2); 0];
    qz = [cos(psi/2);   0; 0; sin(psi/2)];
    q = qmulLocal(qz, qmulLocal(qy, qx));
    q = q / norm(q);
end

function [phi,theta,psi] = quat2eulLocal(q)
% Extraction matching the lab Euler convention (verified to 1e-15).
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    phi   = atan2(2*(q0*q1+q2*q3), 1-2*(q1^2+q2^2));
    theta = asin(max(-1,min(1, 2*(q0*q2-q1*q3))));   % clip vs fp overshoot
    psi   = atan2(2*(q1*q2+q0*q3), 1-2*(q2^2+q3^2));
end

function [b, Mt, ok] = ellipsoidFitLocal(Mdat)
% SVD ellipsoid fit per the calibration derivation; pre-conditioned.
    mu = mean(Mdat,1);
    sc = mean(vecnorm(Mdat - mu, 2, 2));
    n  = (Mdat - mu)/sc;
    xs=n(:,1); ys=n(:,2); zs=n(:,3);
    D = [xs.^2 ys.^2 zs.^2 2*xs.*ys 2*xs.*zs 2*ys.*zs 2*xs 2*ys 2*zs ones(size(xs))];
    [~,~,Vd] = svd(D,'econ');  th = Vd(:,end);
    Qm = [th(1) th(4) th(5); th(4) th(2) th(6); th(5) th(6) th(3)];
    uv = th(7:9);  s = th(10);
    ev = eig((Qm+Qm.')/2);
    if all(ev<0), Qm=-Qm; uv=-uv; s=-s; ev=-ev; end
    ok = all(ev>0);
    if ~ok, b=zeros(3,1); Mt=eye(3); return, end
    bn = -(Qm\uv);
    Mn = Qm/(uv.'*(Qm\uv) - s);
    ok = all(eig((Mn+Mn.')/2) > 0);
    b  = mu.' + sc*bn;
    Mt = Mn/sc^2;
end

function [Xe,Ye,Ze] = ellipsoidSurfLocal(b, Mt, sx, sy, sz)
    [Vm,Lm] = eig((Mt+Mt.')/2);
    T = Vm*diag(1./sqrt(diag(Lm)));
    Pts = T*[sx(:) sy(:) sz(:)].' + b;
    Xe = reshape(Pts(1,:),size(sx));
    Ye = reshape(Pts(2,:),size(sy));
    Ze = reshape(Pts(3,:),size(sz));
end