%% Number of polymers and monomers per polymer
num_polymers = 80;
M = 10; % Monomers per polymer
total_monomers = num_polymers * M; % Total number of monomers

% Physical parameters
sigma = 1;
ks = 100; % Spring constant
kb = 0; % Bending constant
kr = 100; % Repulsion constant
gamma = 1; % Damping coefficient
v0 = 1; % Active velocity
D = 1; % Diffusion coefficient
DR = 1; % Rotational diffusion coefficient
sigmasq = sigma * sigma;
box_size = sqrt(num_polymers) * M * sigma;%box size

% Define a cutoff distance for neighbour interactions
cutoff_distance = 2 * sigma; % 可以自己调节

% Initialize neighbour lists
neighbor_lists = cell(total_monomers, 1);

% Initialize positions with no overlap
positions = initialize_no_overlap(num_polymers, M, sigma);

% Perform energy minimization to remove any potential overlaps
positions = minimize_energy_fire(positions, neighbor_lists, num_polymers, M, ks, kb, kr, sigma, cutoff_distance);

% Initialize orientations for all monomers
theta = 2 * pi * rand(total_monomers, 1); % Random initial orientations

% Setup for visualization and recording
figure;
axis equal;
hold on;
video = VideoWriter('polymer_motion.avi');
open(video);

% Simulation parameters
dt = 0.001;
numSteps = 10000;
all_positions = zeros(total_monomers, 2, numSteps); % Store all positions at each step

%======================新加的部分！！！！！！！！！！



% Main simulation loop
for step = 1:numSteps

    % Update neighbour lists (not every step, depending on your system's dynamics)
    if mod(step, 10) == 0 % Example: update every 10 steps
        for i = 1:total_monomers
           neighbor_lists{i} = find_neighbours(i, positions, cutoff_distance, total_monomers);
        end
    end
    
    
    if step < numSteps
        cla; % Clear the plot to draw new frame
    end
    
    all_positions(:, :, step) = positions; % Store current positions
    

     % Calculate forces for all pairs of monomers, not just adjacent ones
    dUs = zeros(total_monomers, 2); % Spring forces
    dUb = zeros(total_monomers, 2); % Bending forces
    dUr = zeros(total_monomers, 2); % Repulsion forces

    % Spring force calculations
    for p = 1:num_polymers
        for i = 1:(M-1)
            monomer_idx = (p-1)*M + i;
            next_monomer_idx = monomer_idx + 1;
            direction = (positions(next_monomer_idx, :) - positions(monomer_idx, :)) / norm(positions(next_monomer_idx, :) - positions(monomer_idx, :));
            dUs(monomer_idx, :) = dUs(monomer_idx, :) + ks * (norm(positions(next_monomer_idx, :) - positions(monomer_idx, :)) - sigma) * direction;
            dUs(next_monomer_idx, :) = dUs(next_monomer_idx, :) - ks * (norm(positions(next_monomer_idx, :) - positions(monomer_idx, :)) - sigma) * direction;
        end
    end

    
% Bending energy
epsilon = 1e-6;  % Threshold to avoid division by zero

for i = 2:M-1
    vec1 = positions(i, :) - positions(i-1, :);
    vec2 = positions(i+1, :) - positions(i, :);
    
    % Calculate the cosine of the bending angle
    cos_theta = dot(vec1, vec2) / (norm(vec1) * norm(vec2));
    
    % Check if the angle is too small
    if abs(cos_theta) > 1 - epsilon
        % The angle is too small, set the gradient to zero
        continue;
    end
    
    % Calculate the actual bending angle
    theta_i = acos(cos_theta);
    
    % Calculate the derivative of the bending energy with respect to theta_i
    dUb_dtheta_i = -kb * (theta_i);
    
    % Calculate the gradients of cos_theta with respect to positions
    nabla_ri_cos_theta = (vec1 / norm(vec1) - cos_theta * vec2 / norm(vec2)^2) / sqrt(1 - cos_theta^2);
    nabla_riplus1_cos_theta = (vec2 / norm(vec2) - cos_theta * vec1 / norm(vec1)^2) / sqrt(1 - cos_theta^2);
    nabla_riminus1_cos_theta = -nabla_ri_cos_theta - nabla_riplus1_cos_theta;
    
    % Update the bending energy gradient
    dUb(i, :) = dUb(i, :) + dUb_dtheta_i * nabla_ri_cos_theta;
    dUb(i+1, :) = dUb(i+1, :) + dUb_dtheta_i * nabla_riplus1_cos_theta;
    dUb(i-1, :) = dUb(i-1, :) + dUb_dtheta_i * nabla_riminus1_cos_theta;
end



    
    % Repulsion force calculations using neighbour lists
    for i = 1:M
        for j = i+1:M
            xij = positions(i,1) - positions(j,1);
            if(abs(xij)<sigma)
                yij = positions(i,2) - positions(j,2);
                if(abs(yij)<sigma)
                    rijsq = xij^2 + yij^2;
                    if rijsq < sigmasq
                        rij = sqrt(rijsq);
                        dUr(i, :) = dUr(i, :) - kr * (rij - sigma) * [xij, yij] / rij;
                        dUr(j, :) = dUr(j, :) + kr * (rij - sigma) * [xij, yij] / rij;
                    end
                end
            end
        end
    end
    
    % Total force on each monomer
    F = dUs + dUb + dUr;

    % Update positions using the total force and orientation vectors
    n = [cos(theta), sin(theta)]; % Orientation vectors
    positions = positions + dt * (F / gamma + v0 * n) + sqrt(2 * D) * randn(total_monomers, 2) * sqrt(dt);

    % Apply periodic boundary conditions
    positions(i, :) = mod(positions(i, :), box_size);

    % Update orientations with random noise
    theta = theta + sqrt(2 * DR) * randn(total_monomers, 1) * sqrt(dt);


    % Visualization and video writing code
    if(~mod(step,10))
        for p = 1:num_polymers
            polymer_positions = positions((p-1)*M+1:p*M, :);
            plot(polymer_positions(:,1), polymer_positions(:,2), 'o-'); % Plot polymer as a line with circles at monomers
        end
        title(['Step: ', num2str(step)]);
        drawnow;
        frame = getframe(gcf);
        writeVideo(video, frame);
    end
end

close(video); % Close the video writer after the loop




% Function to initialize positions of multiple polymers with no overlap
% and at different angles
function positions = initialize_no_overlap1(num_polymers, M, sigma)
    positions = zeros(num_polymers * M, 2); % Initialize all positions
    % Assuming a 2D square box for simplicity of placement
    box_size = sqrt(num_polymers) * M * sigma; % Estimate a box size

    for p = 1:num_polymers
        % Randomly place the first monomer of each polymer
        placed = false;
        while ~placed
            pos = rand(1, 2) * (box_size - sigma) + sigma / 2;
            placed = true;
            % Check for overlap with other polymers' first monomers
            for other_p = 1:(p-1)
                other_pos = positions((other_p-1) * M + 1, :);
                if norm(pos - other_pos) < sigma
                    placed = false;
                    break;
                end
            end
        end
        % Assign the position of the first monomer of the current polymer
        positions((p-1) * M + 1, :) = pos;

        % Randomly choose an angle for the polymer
        angle = rand() * 2 * pi; 
        direction = [cos(angle), sin(angle)]; % Direction vector

        % Place subsequent monomers in a chain at the chosen angle
        for m = 2:M
            angle = angle + (rand() - 0.5) * 0.2; % Small random deviation
            direction = [cos(angle), sin(angle)];
            positions((p-1) * M + m, :) = positions((p-1) * M + m - 1, :) + sigma * direction;
        end
    end
end




%% Energy minimization function to remove overlaps and minimize potential energy
function positions = minimize_energy_fire(positions, neighbor_lists, num_polymers, M, ks, kb, kr, sigma, cutoff_distance)
    % Directly defined FIRE parameters
    dt = 0.1; % Time step for FIRE
    alpha = 0.1; % Starting alpha for FIRE
    f_dec = 0.5; % Factor for decreasing time step
    f_inc = 1.1; % Factor for increasing time step
    f_alpha = 0.99; % Factor for decreasing alpha
    N_min = 5; % Minimum number of steps before updating dt and alpha
    F_max = 1e-3; % Maximum force for termination
    max_iterations = 10000; % Example maximum number of iterations for the loop

    total_monomers = num_polymers * M;
    velocities = zeros(total_monomers, 2); % Initialize velocities

    % Main optimization loop
    for iteration = 1:max_iterations
        F = calculate_total_forces(positions, neighbor_lists, num_polymers, M, ks, kb, kr, sigma, cutoff_distance);
        power = sum(sum(F .* velocities));

        if power > 0
            velocities = (1 - alpha) * velocities + alpha * (F / norm(F(:))) * norm(velocities(:));
            if iteration > N_min
                dt = min(dt * f_inc, dt); % Assuming dt_max is the same as initial dt
                alpha = alpha * f_alpha;
            end
        else
            velocities = zeros(size(velocities));
            dt = dt * f_dec;
            alpha = 0.1;
        end

        positions = positions + velocities * dt + F * dt^2 / 2;
        velocities = velocities + F * dt / 2;

        % Convergence check (can be based on max force, total energy, etc.)
        max_force = max(sqrt(sum(F.^2, 2)));
        if max_force < F_max
            break;
        end
    end
end


% Function to find neighbours of a given monomer
function neighbours = find_neighbours(monomer_idx, positions, cutoff_distance, total_monomers)
    distances = sqrt(sum((positions - positions(monomer_idx,:)).^2, 2));
    neighbours = find(distances < cutoff_distance & distances > 0 & (1:total_monomers)' ~= monomer_idx);
end

function F = calculate_total_forces(positions, neighbor_lists, num_polymers, M, ks, kb, kr, sigma, cutoff_distance)
    total_monomers = size(positions, 1);
    dUs = zeros(total_monomers, 2); % Spring forces
    dUb = zeros(total_monomers, 2); % Bending forces
    dUr = zeros(total_monomers, 2); % Repulsion forces

    % Spring force calculations
    for p = 1:num_polymers
        for i = 1:(M-1)
            monomer_idx = (p-1)*M + i;
            next_monomer_idx = monomer_idx + 1;
            direction = (positions(next_monomer_idx, :) - positions(monomer_idx, :)) / norm(positions(next_monomer_idx, :) - positions(monomer_idx, :));
            dUs(monomer_idx, :) = dUs(monomer_idx, :) + ks * (norm(positions(next_monomer_idx, :) - positions(monomer_idx, :)) - sigma) * direction;
            dUs(next_monomer_idx, :) = dUs(next_monomer_idx, :) - ks * (norm(positions(next_monomer_idx, :) - positions(monomer_idx, :)) - sigma) * direction;
        end
    end

    % Bending energy
    epsilon = 1e-6;
 
        for i = 2:M-1
            vec1 = positions(i, :) - positions(i-1, :);
            vec2 = positions(i+1, :) - positions(i, :);
            
            % Calculate the cosine of the bending angle
            cos_theta = dot(vec1, vec2) / (norm(vec1) * norm(vec2));
            
            % Check if the angle is too small
            if abs(cos_theta) > 1 - epsilon
                % The angle is too small, set the gradient to zero
                continue;
            end
            
            % Calculate the actual bending angle
            theta_i = acos(cos_theta);
            
            % Calculate the derivative of the bending energy with respect to theta_i
            dUb_dtheta_i = -kb * (theta_i);
            
            % Calculate the gradients of cos_theta with respect to positions
            nabla_ri_cos_theta = (vec1 / norm(vec1) - cos_theta * vec2 / norm(vec2)^2) / sqrt(1 - cos_theta^2);
            nabla_riplus1_cos_theta = (vec2 / norm(vec2) - cos_theta * vec1 / norm(vec1)^2) / sqrt(1 - cos_theta^2);
            nabla_riminus1_cos_theta = -nabla_ri_cos_theta - nabla_riplus1_cos_theta;
            
            % Update the bending energy gradient
            dUb(i, :) = dUb(i, :) + dUb_dtheta_i * nabla_ri_cos_theta;
            dUb(i+1, :) = dUb(i+1, :) + dUb_dtheta_i * nabla_riplus1_cos_theta;
            dUb(i-1, :) = dUb(i-1, :) + dUb_dtheta_i * nabla_riminus1_cos_theta;
        end

    % Repulsion force calculations using neighbor lists
    for i = 1:total_monomers
        for j = neighbor_lists{i}'
            rij = positions(j,:) - positions(i,:);
            r = norm(rij);
            if r < cutoff_distance && r > 0 % Prevent division by zero
                dUr(i, :) = dUr(i, :) - kr * (rij / r) * (cutoff_distance - r);
                dUr(j, :) = dUr(j, :) + kr * (rij / r) * (cutoff_distance - r);
            end
        end
    end

    % Total force on each monomer
    F = dUs + dUb + dUr;
end


% Polymers are not all wingly straight line, 弯曲的 不是统一直线的
function positions = initialize_no_overlap2(num_polymers, M, sigma)
    positions = zeros(num_polymers * M, 2); % Initialize all positions
    % Assuming a 2D square box for simplicity of placement
    box_size = sqrt(num_polymers) * M * sigma; % Estimate a box size

    for p = 1:num_polymers
        for m = 1:M
            placed = false;
            while ~placed
                if m == 1
                    % Randomly place the first monomer of the current polymer
                    pos = rand(1, 2) * (box_size - sigma) + sigma / 2;
                else
                    % Place subsequent monomers in a chain at the chosen angle
                    angle = rand() * 2 * pi + (rand() - 0.5) * 0.2; % Small random deviation
                    direction = [cos(angle), sin(angle)];
                    pos = positions((p-1) * M + m - 1, :) + sigma * direction;
                end

                % Check for overlap with all previously placed monomers
                placed = true;
                for other_monomer_idx = 1:(p-1)*M + m - 1
                    if norm(pos - positions(other_monomer_idx, :)) < sigma
                        placed = false;
                        break;
                    end
                end
            end
            % Assign the position of the current monomer of the current polymer
            positions((p-1) * M + m, :) = pos;
        end
    end
end




function cellIndex = getCellIndex(position, cellSize, numCells)
    % Convert a position to a cell index in the cell-linked list
    cellCoords = ceil(position / cellSize);
    cellIndex = min(max(cellCoords, 1), numCells); % Ensure within bounds
end

function neighboringCells = getNeighboringCells(cellIndex, numCells)
    % Get indices of neighboring cells (including the cell itself)
    neighboringCells = [];
    for dx = -1:1
        for dy = -1:1
            neighborCell = cellIndex + [dx, dy];
            if all(neighborCell > 0 & neighborCell <= numCells)
                neighboringCells = [neighboringCells; neighborCell];
            end
        end
    end
end

function positions = initialize_no_overlap(num_polymers, M, sigma)
    positions = zeros(num_polymers * M, 2); % Initialize all positions
    box_size = sqrt(num_polymers) * M * sigma; % Estimate a box size
    cellSize = 2 * sigma; % Cell size based on interaction range
    numCells = ceil(box_size / cellSize);
    cells = cell(numCells, numCells); % Initialize cell-linked list

    for p = 1:num_polymers
        for m = 1:M
            placed = false;
            while ~placed
                if m == 1
                    pos = rand(1, 2) * (box_size - sigma) + sigma / 2;
                else
                    angle = rand() * 2 * pi;
                    direction = [cos(angle), sin(angle)];
                    pos = positions((p-1) * M + m - 1, :) + sigma * direction;
                end
                
                % Check for overlap using cell-linked list
                cellIdx = getCellIndex(pos, cellSize, numCells);
                neighboringCellsIdx = getNeighboringCells(cellIdx, numCells);
                placed = true;
                for idx = 1:size(neighboringCellsIdx, 1)
                    cellMonomers = cells{neighboringCellsIdx(idx, 1), neighboringCellsIdx(idx, 2)};
                    for otherIdx = cellMonomers'
                        if norm(pos - positions(otherIdx, :)) < sigma
                            placed = false;
                            break;
                        end
                    end
                    if ~placed
                        break;
                    end
                end

                if placed
                    positions((p-1) * M + m, :) = pos;
                    cells{cellIdx(1), cellIdx(2)} = [cells{cellIdx(1), cellIdx(2)}; (p-1)*M + m];
                end
            end
        end
    end
end

