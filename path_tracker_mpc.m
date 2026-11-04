clc;
clear;
close all;

%% ============================================================
%  A* MOTION PLANNING + MPC PATH TRACKING
%
%  Motion Planner : A* Algorithm
%  Controller     : Linear Model Predictive Control (MPC)
%  Vehicle Model  : Kinematic Bicycle Model
%  Optimization   : Projected Gradient Descent
%
%  No additional MATLAB toolboxes required.
% =============================================================

%% ============================================================
%  PART 1: ENVIRONMENT AND MOTION PLANNING
% =============================================================

%% 1. Create environment

Nx = 100;
Ny = 100;

map = zeros(Ny,Nx);

%% 2. Define obstacles

map(30:50,35:55) = 1;
map(60:75,65:80) = 1;

%% 3. Define start and goal

start = [10 10];      % [X Y] (m)
goal  = [90 80];      % [X Y] (m)

%% 4. Inflate obstacles for vehicle clearance

clearance = 3;        % Grid cells; 1 cell = 1 m

inflatedMap = conv2( ...
    map,ones(2*clearance+1),'same') > 0;

if inflatedMap(start(2),start(1)) || ...
   inflatedMap(goal(2),goal(1))

    error('Start or goal is inside an inflated obstacle.');

end

%% 5. Initialize A* algorithm

gCost = inf(Ny,Nx);
fCost = inf(Ny,Nx);

visited = false(Ny,Nx);
openSet = false(Ny,Nx);

parentX = zeros(Ny,Nx);
parentY = zeros(Ny,Nx);

gCost(start(2),start(1)) = 0;

fCost(start(2),start(1)) = norm(goal-start);

openSet(start(2),start(1)) = true;

%% 6. Define movement directions

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

%% 7. A* search

found = false;

while any(openSet(:))

    % Find node with minimum fCost
    candidates = find(openSet);

    [~,idx] = min(fCost(candidates));

    current = candidates(idx);

    [cy,cx] = ind2sub(size(map),current);

    % Check goal
    if cx == goal(1) && cy == goal(2)

        found = true;
        break;

    end

    openSet(cy,cx) = false;
    visited(cy,cx) = true;

    % Explore neighbors
    for j = 1:size(directions,1)

        nx = cx + directions(j,1);
        ny = cy + directions(j,2);

        % Boundary check
        if nx < 1 || nx > Nx || ...
           ny < 1 || ny > Ny

            continue;

        end

        % Obstacle check
        if inflatedMap(ny,nx)

            continue;

        end

        % Visited check
        if visited(ny,nx)

            continue;

        end

        % Prevent diagonal corner cutting
        if directions(j,1) ~= 0 && ...
           directions(j,2) ~= 0

            if inflatedMap(cy,nx) || ...
               inflatedMap(ny,cx)

                continue;

            end

        end

        % Movement cost
        stepCost = norm(directions(j,:));

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

%% 8. Reconstruct A* path

if ~found

    error('No collision-free A* path found.');

end

path = goal;

cx = goal(1);
cy = goal(2);

while ~(cx == start(1) && cy == start(2))

    px = parentX(cy,cx);
    py = parentY(cy,cx);

    path = [px py;path];

    cx = px;
    cy = py;

end

disp('A* motion planning completed successfully.');

%% ============================================================
%  PART 2: PATH SMOOTHING
% =============================================================

%% 9. Smooth A* path

pathSmooth = path;

for iter = 1:150

    candidate = pathSmooth;

    candidate(2:end-1,:) = ...
        0.5*pathSmooth(2:end-1,:) + ...
        0.25*pathSmooth(1:end-2,:) + ...
        0.25*pathSmooth(3:end,:);

    % Accept only collision-free smoothing updates
    if isPathFree(candidate,inflatedMap)

        pathSmooth = candidate;

    else

        break;

    end

end

%% 10. Calculate path arc length

segmentLength = hypot( ...
    diff(pathSmooth(:,1)), ...
    diff(pathSmooth(:,2)));

sOriginal = [0;cumsum(segmentLength)];

[sOriginal,uniqueIdx] = ...
    unique(sOriginal,'stable');

pathSmooth = pathSmooth(uniqueIdx,:);

%% 11. Resample reference path

ds = 0.25;           % Reference spacing (m)

sRef = (0:ds:sOriginal(end))';

if sRef(end) < sOriginal(end)

    sRef(end+1) = sOriginal(end);

end

Xref = interp1( ...
    sOriginal,pathSmooth(:,1),sRef,'pchip');

Yref = interp1( ...
    sOriginal,pathSmooth(:,2),sRef,'pchip');

% If interpolation crosses an obstacle,
% use linear interpolation instead.

if ~isPathFree([Xref Yref],inflatedMap)

    Xref = interp1( ...
        sOriginal,pathSmooth(:,1),sRef,'linear');

    Yref = interp1( ...
        sOriginal,pathSmooth(:,2),sRef,'linear');

end

if ~isPathFree([Xref Yref],inflatedMap)

    error('Reference path is not collision-free.');

end

%% 12. Calculate reference heading

dX = gradient(Xref,sRef);
dY = gradient(Yref,sRef);

PsiRef = unwrap(atan2(dY,dX));

%% 13. Calculate reference curvature

KappaRef = gradient(PsiRef,sRef);

%% 14. Save reference path

referencePath = [Xref Yref PsiRef];

save('referencePath.mat','referencePath');

%% ============================================================
%  PART 3: VEHICLE PARAMETERS
% =============================================================

%% 15. Vehicle parameters

L = 2.8;               % Wheelbase (m)

Vx = 3.0;              % Longitudinal velocity (m/s)

deltaMax = deg2rad(30); % Maximum steering angle (rad)

deltaRateMax = deg2rad(45); % Maximum steering rate (rad/s)

%% ============================================================
%  PART 4: MPC CONTROLLER DESIGN
% =============================================================

%% 16. MPC parameters

Ts = 0.1;              % Sampling time (s)

Np = 15;               % Prediction horizon

Nc = Np;               % Control horizon

% Prediction time = Np*Ts = 1.5 seconds

%% 17. MPC weighting matrices

Q = diag([18 10]);

R = 0.7;

Rd = 8;

% Q  : Tracking error penalty
% R  : Steering correction penalty
% Rd : Steering correction variation penalty

%% 18. Linearized vehicle error model

% State vector:
%
% x = [ey;
%      epsi]
%
% ey   = lateral tracking error (m)
% epsi = heading tracking error (rad)
%
% Control input:
%
% u = steering correction (rad)

Ac = [
    0 Vx;
    0 0
];

Bc = [
    0;
    Vx/L
];

%% 19. Discrete-time model

Ad = [
    1 Vx*Ts;
    0 1
];

Bd = [
    Vx^2*Ts^2/(2*L);
    Vx*Ts/L
];

%% 20. Build prediction matrices

F = zeros(2*Np,2);

G = zeros(2*Np,Nc);

for i = 1:Np

    F(2*i-1:2*i,:) = Ad^i;

    for j = 1:min(i,Nc)

        G(2*i-1:2*i,j) = ...
            Ad^(i-j)*Bd;

    end

end

%% 21. Build MPC cost matrices

Qbar = kron(eye(Np),Q);

Rbar = R*eye(Nc);

% Difference matrix
D = eye(Nc);

for i = 2:Nc

    D(i,i-1) = -1;

end

H = G'*Qbar*G + Rbar + Rd*(D'*D);

% Projected-gradient step size
stepSize = 1/(2*max(eig(H)));

maxIterations = 100;

optimizationTolerance = 1e-5;

%% ============================================================
%  PART 5: SIMULATION INITIALIZATION
% =============================================================

%% 22. Simulation settings

Tmax = 90;

N = ceil(Tmax/Ts);

time = (0:N)'*Ts;

%% 23. Initialize vehicle states

X = zeros(N+1,1);
Y = zeros(N+1,1);
Psi = zeros(N+1,1);

X(1) = start(1);
Y(1) = start(2);

Psi(1) = PsiRef(1);

%% 24. Initialize tracking variables

deltaHistory = zeros(N,1);

eyHistory = zeros(N,1);

epsiHistory = zeros(N,1);

distanceHistory = zeros(N,1);

%% 25. Initialize MPC variables

U = zeros(Nc,1);

previousSteering = 0;

idxLast = 1;

lastStep = N;

reachedGoal = false;
collisionDetected = false;

%% ============================================================
%  PART 6: CLOSED-LOOP MPC PATH TRACKING
% =============================================================

%% 26. Start simulation

for k = 1:N

    %% 26.1 Find closest reference point

    searchIndices = idxLast:min( ...
        idxLast+120,length(sRef));

    distance = hypot( ...
        Xref(searchIndices)-X(k), ...
        Yref(searchIndices)-Y(k));

    [~,localIdx] = min(distance);

    idx = searchIndices(localIdx);

    idxLast = idx;

    %% 26.2 Reference vehicle states

    Xdesired = Xref(idx);

    Ydesired = Yref(idx);

    PsiDesired = PsiRef(idx);

    %% 26.3 Calculate tracking errors

    dx = X(k)-Xdesired;

    dy = Y(k)-Ydesired;

    % Lateral tracking error
    ey = -dx*sin(PsiDesired) + ...
          dy*cos(PsiDesired);

    % Heading tracking error
    epsi = atan2( ...
        sin(Psi(k)-PsiDesired), ...
        cos(Psi(k)-PsiDesired));

    %% 26.4 MPC state vector

    x0 = [
        ey;
        epsi
    ];

    %% 26.5 Future reference positions

    futureDistance = ...
        sRef(idx)+Vx*Ts*(0:Np-1)';

    futureDistance = min( ...
        sRef(end),futureDistance);

    %% 26.6 Future path curvature

    futureCurvature = interp1( ...
        sRef,KappaRef,futureDistance,'linear');

    %% 26.7 Feedforward steering

    deltaFF = atan(L*futureCurvature);

    deltaFF = max( ...
        -deltaMax,min(deltaMax,deltaFF));

    %% 26.8 Build optimization cost

    % Penalize changes in steering correction
    % relative to the previous applied correction.

    d0 = zeros(Nc,1);

    d0(1) = previousSteering-deltaFF(1);

    f = G'*Qbar*(F*x0) - Rd*(D'*d0);

    %% 26.9 Shift previous MPC solution

    U = [U(2:end);U(end)];

    %% 26.10 MPC optimization

    for iteration = 1:maxIterations

        % Cost-function gradient
        gradientJ = 2*(H*U+f);

        % Gradient descent
        Unew = U-stepSize*gradientJ;

        %% Apply steering constraints

        previous = previousSteering;

        for j = 1:Nc

            % Total steering command
            desiredSteering = deltaFF(j)+Unew(j);

            % Steering-rate constraint
            lowerBound = max( ...
                -deltaMax, ...
                previous-deltaRateMax*Ts);

            upperBound = min( ...
                deltaMax, ...
                previous+deltaRateMax*Ts);

            % Apply constraints
            desiredSteering = min( ...
                upperBound, ...
                max(lowerBound,desiredSteering));

            % Convert back to correction
            Unew(j) = desiredSteering-deltaFF(j);

            previous = desiredSteering;

        end

        %% Check optimization convergence

        if norm(Unew-U) < optimizationTolerance

            U = Unew;
            break;

        end

        U = Unew;

    end

    %% 26.11 Apply first optimal steering command

    delta = deltaFF(1)+U(1);

    previousSteering = delta;

    %% 26.12 Save controller outputs

    deltaHistory(k) = delta;

    eyHistory(k) = ey;

    epsiHistory(k) = epsi;

    distanceHistory(k) = hypot(dx,dy);

    %% ========================================================
    %  PART 7: VEHICLE DYNAMICS
    % =========================================================

    %% 27. Kinematic bicycle model

    yawRate = Vx/L*tan(delta);

    % Midpoint integration
    PsiMid = Psi(k)+0.5*Ts*yawRate;

    X(k+1) = X(k)+Vx*Ts*cos(PsiMid);

    Y(k+1) = Y(k)+Vx*Ts*sin(PsiMid);

    Psi(k+1) = Psi(k)+Ts*yawRate;

    %% 28. Collision detection

    vehicleSegment = [
        X(k)   Y(k);
        X(k+1) Y(k+1)
    ];

    if ~isPathFree(vehicleSegment,inflatedMap)

        collisionDetected = true;

        warning('Vehicle entered an inflated obstacle.');

        lastStep = k;

        break;

    end

    %% 29. Check goal

    goalDistance = hypot( ...
        X(k+1)-goal(1), ...
        Y(k+1)-goal(2));

    if goalDistance < 1.0

        reachedGoal = true;

        lastStep = k;

        break;

    end

end

%% ============================================================
%  PART 8: PROCESS RESULTS
% =============================================================

%% 30. Trim simulation results

X = X(1:lastStep+1);

Y = Y(1:lastStep+1);

Psi = Psi(1:lastStep+1);

deltaHistory = deltaHistory(1:lastStep);

eyHistory = eyHistory(1:lastStep);

epsiHistory = epsiHistory(1:lastStep);

distanceHistory = distanceHistory(1:lastStep);

time = time(1:lastStep+1);

%% 31. Calculate performance metrics

RMSE_y = sqrt(mean(eyHistory.^2));

MAE_y = mean(abs(eyHistory));

MaxError = max(abs(eyHistory));

RMSE_psi = sqrt(mean(epsiHistory.^2));

MaxSteering = max(abs(rad2deg(deltaHistory)));

FinalDistance = hypot( ...
    X(end)-goal(1), ...
    Y(end)-goal(2));

%% 32. Display results

fprintf('\n========================================\n');

fprintf('A* + MPC PATH TRACKING RESULTS\n');

fprintf('========================================\n');

fprintf('Goal reached: %s\n',mat2str(reachedGoal));

fprintf('Collision detected: %s\n', ...
    mat2str(collisionDetected));

fprintf('Simulation time: %.2f s\n',time(end));

fprintf('Lateral RMSE: %.4f m\n',RMSE_y);

fprintf('Lateral MAE: %.4f m\n',MAE_y);

fprintf('Maximum lateral error: %.4f m\n',MaxError);

fprintf('Heading RMSE: %.4f deg\n', ...
    rad2deg(RMSE_psi));

fprintf('Maximum steering angle: %.2f deg\n', ...
    MaxSteering);

fprintf('Final distance to goal: %.3f m\n', ...
    FinalDistance);

fprintf('========================================\n');

%% ============================================================
%  PART 9: VISUALIZATION
% =============================================================

%% 33. Plot planned path and vehicle trajectory

figure('Name','A* + MPC Motion Planning and Control');

imagesc([1 Nx],[1 Ny],inflatedMap);

set(gca,'YDir','normal');

colormap(gray);

hold on;

h1 = plot( ...
    path(:,1),path(:,2), ...
    'r--','LineWidth',1.5);

h2 = plot( ...
    Xref,Yref, ...
    'b','LineWidth',2);

h3 = plot( ...
    X,Y, ...
    'g','LineWidth',2);

h4 = plot( ...
    start(1),start(2), ...
    'go','MarkerFaceColor','g');

h5 = plot( ...
    goal(1),goal(2), ...
    'ro','MarkerFaceColor','r');

xlabel('X Position (m)');

ylabel('Y Position (m)');

title('A* Motion Planning + MPC Path Tracking');

legend( ...
    [h1 h2 h3 h4 h5], ...
    {'A* Path', ...
     'Smoothed Reference Path', ...
     'MPC Vehicle Trajectory', ...
     'Start', ...
     'Goal'}, ...
    'Location','best');

axis equal;
grid on;

%% 34. Plot lateral tracking error

figure('Name','Lateral Tracking Error');

plot( ...
    time(1:end-1),eyHistory, ...
    'b','LineWidth',2);

xlabel('Time (s)');

ylabel('Lateral Error (m)');

title('MPC Lateral Tracking Error');

grid on;

%% 35. Plot heading tracking error

figure('Name','Heading Tracking Error');

plot( ...
    time(1:end-1),rad2deg(epsiHistory), ...
    'r','LineWidth',2);

xlabel('Time (s)');

ylabel('Heading Error (deg)');

title('MPC Heading Tracking Error');

grid on;

%% 36. Plot steering angle

figure('Name','MPC Steering Control');

stairs( ...
    time(1:end-1),rad2deg(deltaHistory), ...
    'k','LineWidth',2);

hold on;

yline(rad2deg(deltaMax),'r--');

yline(-rad2deg(deltaMax),'r--');

xlabel('Time (s)');

ylabel('Steering Angle (deg)');

title('MPC Steering Angle');

grid on;

%% 37. Plot vehicle heading

figure('Name','Vehicle Heading');

plot( ...
    time,rad2deg(unwrap(Psi)), ...
    'LineWidth',2);

xlabel('Time (s)');

ylabel('Vehicle Heading (deg)');

title('Vehicle Yaw Angle');

grid on;

%% ============================================================
%  PART 10: SAVE RESULTS
% =============================================================

%% 38. Save simulation data

results.time = time;

results.X = X;
results.Y = Y;
results.Psi = Psi;

results.delta = deltaHistory;

results.ey = eyHistory;
results.epsi = epsiHistory;

results.referencePath = referencePath;

results.RMSE_y = RMSE_y;
results.MAE_y = MAE_y;
results.MaxError = MaxError;

results.reachedGoal = reachedGoal;
results.collisionDetected = collisionDetected;

save('AStar_MPC_Results.mat','results');

disp('MPC simulation completed.');
disp('Results saved to AStar_MPC_Results.mat');

%% ============================================================
%  LOCAL FUNCTION: PATH COLLISION CHECK
% =============================================================

function free = isPathFree(points,occupiedMap)

    [nY,nX] = size(occupiedMap);

    free = true;

    for i = 1:size(points,1)-1

        p0 = points(i,:);

        p1 = points(i+1,:);

        segmentLength = norm(p1-p0);

        nSamples = max( ...
            2,ceil(segmentLength/0.15)+1);

        for j = 0:nSamples-1

            alpha = j/(nSamples-1);

            p = (1-alpha)*p0 + alpha*p1;

            ix = round(p(1));
            iy = round(p(2));

            if ix < 1 || ix > nX || ...
               iy < 1 || iy > nY || ...
               occupiedMap(iy,ix)

                free = false;

                return;

            end

        end

    end

end