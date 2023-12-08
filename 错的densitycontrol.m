%有neighbur,cell-linking , periodic, random orientation,density control
%=======================================================================
%% Number of polymers and monomers per polymer
num_polymers = 40;
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
%box_size = sqrt(num_polymers) * M * sigma;
box_size = 100;%box size

% Define a cutoff distance for neighbour interactions
cutoff_distance = 2 * sigma; % 可以自己调节

% Initialize neighbour lists
neighbor_lists = cell(total_monomers, 1);

% Initialize positions with no overlap
positions = initialize_no_overlap(num_polymers, M, sigma);

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
numSteps = 1000;
all_positions = zeros(total_monomers, 2, numSteps); % Store all positions at each step

%======================新加的部分！！！！！！！！！！



% Main simulation loop

% Initialize parameters
density_check_interval = 10; % Interval to check and adjust density

% Initial density parameters
target_density = 0.4; % Starting target density, e.g., 40%
final_density = 0.9; % Final target density, e.g., 90%
density_increment = 0.1; % Increment in density for each step
shrink_factor = 0.9; % Factor by which the box is shrunk each time

% Main simulation loop
for step = 1:numSteps
    % Update neighbour lists (not every step, depending on your system's dynamics)
    if mod(step, 50) == 0 % Example: update every 10 steps
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

    

    %Bending energy
    epsilon = 1e-6; % Threshold to avoid division by zero
    
    for p = 1:num_polymers
        for i = 2:(M-1)
            monomer_idx = (p-1)*M + i;
            prev_monomer_idx = monomer_idx - 1;
            next_monomer_idx = monomer_idx + 1;
    
            vec1 = positions(monomer_idx, :) - positions(prev_monomer_idx, :);
            vec2 = positions(next_monomer_idx, :) - positions(monomer_idx, :);
    
            cos_theta = dot(vec1, vec2) / (norm(vec1) * norm(vec2));
    
            if abs(cos_theta) > 1 - epsilon
                continue;
            end
                theta_i = acos(cos_theta);
                dUb_dtheta_i = -kb * (theta_i);
    
                nabla_ri_cos_theta = (vec1 / norm(vec1) - cos_theta * vec2 / norm(vec2)^2) / sqrt(1 - cos_theta^2);
                nabla_riplus1_cos_theta = (vec2 / norm(vec2) - cos_theta * vec1 / norm(vec1)^2) / sqrt(1 - cos_theta^2);
                nabla_riminus1_cos_theta = -nabla_ri_cos_theta - nabla_riplus1_cos_theta;
    
                dUb(monomer_idx, :) = dUb(monomer_idx, :) + dUb_dtheta_i * nabla_ri_cos_theta;
                dUb(next_monomer_idx, :) = dUb(next_monomer_idx, :) + dUb_dtheta_i * nabla_riplus1_cos_theta;
                dUb(prev_monomer_idx, :) = dUb(prev_monomer_idx, :) + dUb_dtheta_i * nabla_riminus1_cos_theta;
        end
    end



    % Repulsion force calculations using neighbour lists
    for i = 1:total_monomers
        for j_idx = 1:length(neighbor_lists{i})
            j = neighbor_lists{i}(j_idx); % Index of the neighboring monomer
            rij = positions(j,:) - positions(i,:);
            r = norm(rij);
            if r < cutoff_distance && r > 0 % Avoid division by zero
                dUr(i, :) = dUr(i, :) - kr * (rij / r) * (cutoff_distance - r);
                dUr(j, :) = dUr(j, :) + kr * (rij / r) * (cutoff_distance - r);
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


     overlap_threshold = sigma; % Define an overlap threshold
    % Check and adjust density at specific intervals
    if mod(step, density_check_interval) == 0
        current_density = calculate_density(num_polymers, M, sigma, box_size)

        % If the current density is less than the target density
        if current_density < target_density
            % Shrink the box to increase the density
            [positions, box_size] = shrink_box(positions, box_size, shrink_factor);

            % Apply the FIRE algorithm to minimize potential energy and remove overlaps
            positions = minimize_energy_fire(positions, neighbor_lists, num_polymers, M, ks, kb, kr, sigma, cutoff_distance,box_size);

        % Check for and handle overlaps
        for i = 1:total_monomers
            for j = (i+1):total_monomers
                distance = norm(positions(i,:) - positions(j,:));
                if distance < overlap_threshold
                    % Calculate a displacement vector
                    displacement = (overlap_threshold - distance) * (positions(j,:) - positions(i,:)) / distance;
                    
                    % Reposition the monomers
                    positions(i,:) = positions(i,:) - displacement / 2;
                    positions(j,:) = positions(j,:) + displacement / 2;
                end
            end
        end


            % Recalculate density after shrinking the box
            current_density = calculate_density(num_polymers, M, sigma, box_size);

            % Update the target density for the next iteration, if needed
            if current_density >= target_density
                target_density = min(target_density + density_increment, final_density);
            end
        end
    end


 
    % Visualization and video writing code
    if(~mod(step,10))
        cla; % Clear the plot to draw new frame
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


%% Energy minimization function to remove overlaps and minimize potential energy

function positions = minimize_energy_fire(positions, neighbor_lists, num_polymers, M, ks, kb, kr, sigma, cutoff_distance,box_size)
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
        F = calculate_total_forces(positions, neighbor_lists, num_polymers, M, ks, kb, kr, sigma, cutoff_distance,box_size);
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

      
        % Update positions with periodic boundary conditions
        positions = positions + velocities * dt + F * dt^2 / 2;
        for i = 1:size(positions, 1)
            positions(i, :) = mod(positions(i, :), box_size); % Apply PBC in both dimensions
        end


        % Update velocities
        velocities = velocities + F * dt / 2;

        % Convergence check (can be based on max force, total energy, etc.)
        max_force = max(sqrt(sum(F.^2, 2)));
        if max_force < F_max
            break;
        end
    end
end

% Function to find neighbours of a given monomer
function neighbours = find_neighbours(monomer_idx, positions, cutoff_distance, total_monomers,box_size)
    distances = sqrt(sum((positions - positions(monomer_idx,:)).^2, 2));
    neighbours = find(distances < cutoff_distance & distances > 0 & (1:total_monomers)' ~= monomer_idx);
end

function F = calculate_total_forces(positions, neighbor_lists, num_polymers, M, ks, kb, kr, sigma, cutoff_distance,box_size)
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
    for p = 1:num_polymers
        for i = 2:M-1
            monomer_idx = (p-1)*M + i;
            vec1 = positions(monomer_idx, :) - positions(monomer_idx-1, :);
            vec2 = positions(monomer_idx+1, :) - positions(monomer_idx, :);
            cos_theta = dot(vec1, vec2) / (norm(vec1) * norm(vec2));

            if abs(cos_theta) > 1 - epsilon
                continue;
            end

            theta_i = acos(cos_theta);

            dUb_dtheta_i = -kb * (theta_i);

            nabla_ri_cos_theta = (vec1 / norm(vec1) - cos_theta * vec2 / norm(vec2)^2) / sqrt(1 - cos_theta^2);
            nabla_riplus1_cos_theta = (vec2 / norm(vec2) - cos_theta * vec1 / norm(vec1)^2) / sqrt(1 - cos_theta^2);
            nabla_riminus1_cos_theta = -nabla_ri_cos_theta - nabla_riplus1_cos_theta;

            dUb(monomer_idx, :) = dUb(monomer_idx, :) + dUb_dtheta_i * nabla_ri_cos_theta;
            dUb(monomer_idx+1, :) = dUb(monomer_idx+1, :) + dUb_dtheta_i * nabla_riplus1_cos_theta;
            dUb(monomer_idx-1, :) = dUb(monomer_idx-1, :) + dUb_dtheta_i * nabla_riminus1_cos_theta;
        end
    end

    % Repulsion force calculations using neighbor lists
    for i = 1:total_monomers
        for j = neighbor_lists{i}'
            rij = positions(j,:) - positions(i,:);

            % Apply periodic boundary conditions
            for dim = 1:2
                if rij(dim) > box_size / 2
                    rij(dim) = rij(dim) - box_size;
                elseif rij(dim) < -box_size / 2
                    rij(dim) = rij(dim) + box_size;
                end
            end

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

% this is the inital of the polymer is not a wiggle straight line
% function positions = initialize_no_overlap(num_polymers, M, sigma)
%     positions = zeros(num_polymers * M, 2); % Initialize all positions
%     box_size = 100; % Estimate a box size
%     cellSize = 2 * sigma; % Cell size based on interaction range
%     numCells = ceil(box_size / cellSize);
%     cells = cell(numCells, numCells); % Initialize cell-linked list
% 
%     for p = 1:num_polymers
%         for m = 1:M
%             placed = false;
%             while ~placed
%                 if m == 1
%                     pos = rand(1, 2) * (box_size - sigma) + sigma / 2;
%                 else
%                     angle = rand() * 2 * pi;
%                     direction = [cos(angle), sin(angle)];
%                     pos = positions((p-1) * M + m - 1, :) + sigma * direction;
%                 end
% 
%                 % Check for overlap using cell-linked list
%                 cellIdx = getCellIndex(pos, cellSize, numCells);
%                 neighboringCellsIdx = getNeighboringCells(cellIdx, numCells);
%                 placed = true;
%                 for idx = 1:size(neighboringCellsIdx, 1)
%                     cellMonomers = cells{neighboringCellsIdx(idx, 1), neighboringCellsIdx(idx, 2)};
%                     for otherIdx = cellMonomers'
%                         if norm(pos - positions(otherIdx, :)) < sigma
%                             placed = false;
%                             break;
%                         end
%                     end
%                     if ~placed
%                         break;
%                     end
%                 end
% 
%                 if placed
%                     positions((p-1) * M + m, :) = pos;
%                     cells{cellIdx(1), cellIdx(2)} = [cells{cellIdx(1), cellIdx(2)}; (p-1)*M + m];
%                 end
%             end
%         end
%     end
% end



% The function to simulate a wormlike chain of active beads (monomers) for each 
% polymer, while also generating polymers at different angles, we need to 
% make a few adjustments to your existing code. The goal is to ensure each 
% polymer is initialized as a slightly wiggly line (to simulate the wormlike
% chain) and each polymer starts at a different random angle.


function positions = initialize_no_overlap(num_polymers, M, sigma)
    positions = zeros(num_polymers * M, 2); % Initialize all positions
    box_size = 100; % Estimate a box size
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



%% Density Check
function density = calculate_density(num_polymers, M, sigma, box_size)
    % Volume occupied by each monomer (assuming circular area)
    monomer_volume = pi * (sigma/2)^2;
    
    % Total volume occupied by all monomers
    total_monomer_volume = num_polymers * M * monomer_volume;
    
    % Volume of the box
    box_volume = box_size^2;
    
    % Density calculation
    density = total_monomer_volume / box_volume;
end

function [new_positions, new_box_size] = shrink_box(positions, box_size, shrink_factor)
    % Shrink the box
    new_box_size = box_size * shrink_factor;
    
    % Scale the positions of the monomers
    new_positions = positions * shrink_factor;
end

