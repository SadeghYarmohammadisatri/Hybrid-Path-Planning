clc;
clear;
close all;

%% 1. Create occupancy map
map = binaryOccupancyMap(100,100,1);

%% 2. Add obstacles
[x1,y1] = meshgrid(35:55,30:50);
setOccupancy(map,[x1(:),y1(:)],1);

[x2,y2] = meshgrid(65:80,60:75);
setOccupancy(map,[x2(:),y2(:)],1);

%% 3. Account for vehicle clearance
inflate(map,2);

%% 4. Create vehicle configuration space
ss = stateSpaceSE2;

ss.StateBounds = [
    map.XWorldLimits;
    map.YWorldLimits;
    -pi pi
];

%% 5. Collision validator
validator = validatorOccupancyMap(ss);
validator.Map = map;
validator.ValidationDistance = 0.1;

%% 6. Create Hybrid A* planner
planner = plannerHybridAStar(validator, ...
    MinTurningRadius=6, ...
    MotionPrimitiveLength=3);

%% 7. Define start and goal
startPose = [10 10 0];
goalPose  = [90 80 0];

%% 8. Generate path
pathObj = plan(planner,startPose,goalPose);

path = pathObj.States;

if isempty(path)
    error('No path found.');
end

%% 9. Extract reference path
Xref = path(:,1);
Yref = path(:,2);
Psiref = path(:,3);

%% 10. Visualize
figure;
show(map);
hold on;

plot(Xref,Yref,'b','LineWidth',2);

plot(startPose(1),startPose(2),'go', ...
    'MarkerFaceColor','g');

plot(goalPose(1),goalPose(2),'ro', ...
    'MarkerFaceColor','r');

title('Hybrid A* Motion Planning');
xlabel('X (m)');
ylabel('Y (m)');
grid on;
axis equal;