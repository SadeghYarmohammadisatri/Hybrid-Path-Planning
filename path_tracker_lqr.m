clc;
clear;
close all;

%% ============================================================
%  A* MOTION PLANNING + LQR MOTION CONTROL
%  Base MATLAB implementation (no toolboxes required)
% =============================================================

%% 1. Create environment

Nx = 100;
Ny = 100;

map = zeros(Ny,Nx);

%% 2. Define obstacles

map(30:50,35:55) = 1;
map(60:75,65:80) = 1;

%% 3. Define start and goal

start = [10 10];
goal  = [90 80];

%% 4. Initialize A* algorithm

gCost = inf(Ny,Nx);
fCost = inf(Ny,Nx);

visited = false(Ny,Nx);
openSet = false(Ny,Nx);

parentX = zeros(Ny,Nx);
parentY = zeros(Ny,Nx);

gCost(start(2),start(1)) = 0;

fCost(start(2),start(1)) = norm(goal-start);

openSet(start(2),start(1)) = true;

%% 5. A* search

directions = [
     1  0;
    -1  0;
     0  1;
     0 -1;
     1  1;
     1 -1;
    -1  1;
    -1 -1
];

found = false;

while any(openSet(:))

    candidates = find(openSet);

    [~,idx] = min(fCost(candidates));

    current = candidates(idx);

    [cy,cx] = ind2sub(size(map),current);

    if cx == goal(1) && cy == goal(2)
        found = true;
        break;
    end

    openSet(cy,cx) = false;
    visited(cy,cx) = true;

    for k = 1:size(directions,1)

        nx = cx + directions(k,1);
        ny = cy + directions(k,2);

        if nx < 1 || nx > Nx || ...
           ny < 1 || ny > Ny
            continue;
        end

        if map(ny,nx) == 1
            continue;
        end

        if visited(ny,nx)
            continue;
        end

        if directions(k,1) ~= 0 && directions(k,2) ~= 0

            if map(cy,nx) == 1 || map(ny,cx) == 1
                continue;
            end

        end

        stepCost = norm(directions(k,:));

        tentativeG = gCost(cy,cx) + stepCost;

        if tentativeG < gCost(ny,nx)

            parentX(ny,nx) = cx;
            parentY(ny,nx) = cy;

            gCost(ny,nx) = tentativeG;

            h = norm(goal-[nx ny]);

            fCost(ny,nx) = tentativeG+h;

            openSet(ny,nx) = true;

        end

    end

end

%% 6. Reconstruct path

if ~found
    error('No path found');
end

path = goal;

cx = goal(1);
cy = goal(2);

while ~(cx == start(1) && cy == start(2))

    px = parentX(cy,cx);
    py = parentY(cy,cx);

    path = [px py; path];

    cx = px;
    cy = py;

end

%% ============================================================
%  PATH PREPARATION
% =============================================================

%% 7. Smooth the A* path

% Remove unnecessary collinear waypoints and smooth corners
% using iterative averaging.

pathSmooth = path;

for iter = 1:100

    candidate = pathSmooth;

    for i = 2:size(pathSmooth,1)-1

        candidate(i,:) = ...
            0.5*pathSmooth(i,:) + ...
            0.25*pathSmooth(i-1,:) + ...
            0.25*pathSmooth(i+1,:);

    end

    % Reject a smoothing iteration if it creates a collision.
    % Sampling is used to check each path segment.

    valid = true;

    for i = 1:size(candidate,1)-1

        p1 = candidate(i,:);
        p2 = candidate(i+1,:);

        segmentLength = norm(p2-p1);
        nCheck = max(2,ceil(segmentLength/0.2)+1);

        for j = 0:nCheck-1

            a = j/(nCheck-1);

            p = (1-a)*p1 + a*p2;

            ix = round(p(1));
            iy = round(p(2));

            if ix < 1 || ix > Nx || ...
               iy < 1 || iy > Ny || ...
               map(iy,ix) == 1

                valid = false;
                break;

            end

        end

        if ~valid
            break;
        end

    end

    if valid
        pathSmooth = candidate;
    else
        break;
    end

end

%% 8. Resample path by arc length

dsOriginal = hypot(diff(pathSmooth(:,1)), ...
                   diff(pathSmooth(:,2)));

sOriginal = [0; cumsum(dsOriginal)];

% Remove any duplicate arc-length points
[sUnique,uniqueIdx] = unique(sOriginal,'stable');

pathUnique = pathSmooth(uniqueIdx,:);

sRef = (0:0.25:sUnique(end))';

if sRef(end) < sUnique(end)
    sRef = [sRef; sUnique(end)];
end

Xref = interp1(sUnique,pathUnique(:,1),sRef,'linear');
Yref = interp1(sUnique,pathUnique(:,2),sRef,'linear');

%% 9. Calculate heading and curvature

dx = gradient(Xref,sRef);
dy = gradient(Yref,sRef);

Psiref = unwrap(atan2(dy,dx));

dpsi = gradient(Psiref,sRef);

kappaRef = dpsi;

referencePath = [Xref Yref Psiref];

save('referencePath.mat','referencePath');

%% ============================================================
%  VEHICLE MODEL
% =============================================================

%% 10. Vehicle parameters

L = 2.8;               % Wheelbase (m)
v = 3.0;               % Constant vehicle speed (m/s)

deltaMax = deg2rad(30); % Maximum steering angle

%% 11. Simulation settings

dt = 0.05;
Tmax = 80;

N = round(Tmax/dt);

X = zeros(N+1,1);
Y = zeros(N+1,1);
psi = zeros(N+1,1);

deltaHistory = zeros(N,1);
eyHistory = zeros(N,1);
epsiHistory = zeros(N,1);
time = (0:N)'*dt;

%% 12. Initial vehicle state

X(1) = start(1);
Y(1) = start(2);

psi(1) = Psiref(1);

%% ============================================================
%  LQR CONTROLLER DESIGN
% =============================================================

%% 13. Linearized error dynamics

Ac = [
    0 v;
    0 0
];

Bc = [
    0;
    v/L
];

% Discretization using exact ZOH formulas
Ad = [
    1 v*dt;
    0 1
];

Bd = [
    v^2*dt^2/(2*L);
    v*dt/L
];

%% 14. LQR weighting matrices

Q = diag([15 8]);

R = 2;

%% 15. Solve discrete Riccati equation iteratively

P = Q;

for iter = 1:10000

    Ktemp = (R + Bd'*P*Bd) \ (Bd'*P*Ad);

    Pnew = Q + Ad'*P*Ad ...
         - Ad'*P*Bd*Ktemp;

    if norm(Pnew-P,'fro') < 1e-10
        P = Pnew;
        break;
    end

    P = Pnew;

end

K = (R + Bd'*P*Bd) \ (Bd'*P*Ad);

disp('LQR Gain Matrix:');
disp(K);

%% ============================================================
%  CLOSED-LOOP SIMULATION
% =============================================================

%% 16. Motion control

lastStep = N;

for k = 1:N

    % Find closest reference point
    distance = hypot(Xref-X(k),Yref-Y(k));

    [~,idx] = min(distance);

    % Reference state
    Xd = Xref(idx);
    Yd = Yref(idx);
    psid = Psiref(idx);

    % Lateral tracking error
    ey = -(X(k)-Xd)*sin(psid) ...
         +(Y(k)-Yd)*cos(psid);

    % Heading error
    epsi = atan2( ...
        sin(psi(k)-psid), ...
        cos(psi(k)-psid));

    % Error state
    errorState = [ey; epsi];

    % Feedforward steering
    deltaFF = atan(L*kappaRef(idx));

    % LQR feedback
    deltaFB = -K*errorState;

    % Total steering command
    deltaCmd = deltaFF + deltaFB;

    % Steering saturation
    delta = max(-deltaMax,min(deltaMax,deltaCmd));

    % Save controller outputs
    deltaHistory(k) = delta;
    eyHistory(k) = ey;
    epsiHistory(k) = epsi;

    % Kinematic bicycle model
    X(k+1) = X(k) + v*cos(psi(k))*dt;

    Y(k+1) = Y(k) + v*sin(psi(k))*dt;

    psi(k+1) = psi(k) ...
             + v/L*tan(delta)*dt;

    % Stop when vehicle reaches the goal
    goalDistance = hypot( ...
        X(k+1)-goal(1), ...
        Y(k+1)-goal(2));

    if goalDistance < 1.0

        lastStep = k;
        break;

    end

end

%% 17. Trim simulation arrays

X = X(1:lastStep+1);
Y = Y(1:lastStep+1);
psi = psi(1:lastStep+1);

deltaHistory = deltaHistory(1:lastStep);
eyHistory = eyHistory(1:lastStep);
epsiHistory = epsiHistory(1:lastStep);

time = time(1:lastStep+1);

%% ============================================================
%  RESULTS
% =============================================================

%% 18. Plot motion planning and vehicle trajectory

figure;

imagesc([1 Nx],[1 Ny],map);

set(gca,'YDir','normal');

colormap(gray);

hold on;

plot(path(:,1),path(:,2), ...
    'r--','LineWidth',1.5);

plot(Xref,Yref, ...
    'b','LineWidth',2);

plot(X,Y, ...
    'g','LineWidth',2);

plot(start(1),start(2), ...
    'ko','MarkerFaceColor','g');

plot(goal(1),goal(2), ...
    'ko','MarkerFaceColor','r');

xlabel('X Position (m)');
ylabel('Y Position (m)');

title('A* Motion Planning + LQR Motion Control');

legend('A* Path','Smoothed Reference', ...
       'Vehicle Trajectory','Start','Goal');

axis equal;
grid on;

%% 19. Lateral tracking error

figure;

plot(time(1:end-1),eyHistory, ...
    'LineWidth',2);

xlabel('Time (s)');
ylabel('Lateral Error (m)');

title('LQR Lateral Tracking Error');

grid on;

%% 20. Heading error

figure;

plot(time(1:end-1),rad2deg(epsiHistory), ...
    'LineWidth',2);

xlabel('Time (s)');
ylabel('Heading Error (deg)');

title('LQR Heading Tracking Error');

grid on;

%% 21. Steering angle

figure;

plot(time(1:end-1),rad2deg(deltaHistory), ...
    'LineWidth',2);

xlabel('Time (s)');
ylabel('Steering Angle (deg)');

title('LQR Steering Control');

grid on;

%% 22. Tracking performance

RMSE_y = sqrt(mean(eyHistory.^2));

MaxError = max(abs(eyHistory));

fprintf('\nTracking Performance:\n');

fprintf('Lateral RMSE: %.3f m\n',RMSE_y);

fprintf('Maximum lateral error: %.3f m\n',MaxError);

fprintf('Final distance to goal: %.3f m\n', ...
    hypot(X(end)-goal(1),Y(end)-goal(2)));