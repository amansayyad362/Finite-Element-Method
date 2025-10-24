% Create PDE model
model = createpde('thermal','steadystate');  % Steady-state thermal problem

% Define geometry: rectangle [0,2] x [0,3]
rect = [3 4 0 2 2 0 0 0 3 3]';  % Rectangle in decsg format
gd = rect;
ns = char('R1');
ns = ns';
sf = 'R1';
g = decsg(gd,sf,ns);
geometryFromEdges(model,g);

% Generate mesh
meshHmax = 0.05;   % approximate element size
generateMesh(model,'Hmax',meshHmax);


% Set thermal properties (analogous to conductivity h)
thermalProperties(model,'ThermalConductivity',50);

% Set internal heat source (analogous to q0)
internalHeatSource(model,500);

% Apply Dirichlet boundary conditions (temperature = 100) on all edges
thermalBC(model,'Edge',1:model.Geometry.NumEdges,'Temperature',100);

% Solve
result = solve(model);

% Extract temperature
T = result.Temperature;

disp(['Max T: ', num2str(max(T))])
disp(['Min T: ', num2str(min(T))])