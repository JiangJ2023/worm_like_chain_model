%% Number of polymers and monomers per polymer
num_polymers = 100;
M = 10; % Monomers per polymer
total_monomers = num_polymers * M; % Total number of monomers

% Physical parameters
sigma = 1;
ks = 100; % Spring constant
kb = 100; % Bending constant
kr = 100; % Repulsion constant
gamma = 1; % Damping coefficient
v0 = 1; % Active velocity
D = 1; % Diffusion coefficient
DR = 1; % Rotational diffusion coefficient

% Initialize positions with no overlap
positions = initialize_no_overlap(num_polymers, M, sigma);

% Perform energy minimization to remove any potential overlaps
positions = minimize_energy(positions, num_polymers, M, sigma, ks);

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
% Define a cutoff distance for neighbour interactions
cutoff_distance = 2 * sigma; % 可以自己调节

% Initialize neighbour lists
neighbour_lists = cell(total_monomers, 1);


% Main simulation loop
for step = 1:numSteps

    % Update neighbour lists (not every step, depending on your system's dynamics)
    if mod(step, 10) == 0 % Example: update every 10 steps
        for i = 1:total_monomers
            neighbour_lists{i} = find_neighbours(i, positions, cutoff_distance, total_monomers);
        end
    end
    
    
    if step < numSteps
        cla; % Clear the plot to draw new frame
    end
    
    all_positions(:, :, step) = positions; % Store current positions
    
   % Calculate pairwise differences using broadcasting
   % Pairwise Distance Calculation: Uses broadcasting (repmat) to compute distances between all pairs of monomers.
    dx = repmat(positions(:,1), 1, total_monomers) - repmat(positions(:,1).', total_monomers, 1);
    dy = repmat(positions(:,2), 1, total_monomers) - repmat(positions(:,2).', total_monomers, 1);
    dr = sqrt(dx.^2 + dy.^2); % Pairwise distances


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
       epsilon = 1e-6;
    for i = 2:M-1
        vec1 = positions(i, :) - positions(i-1, :);
        vec2 = positions(i+1, :) - positions(i, :);
        cos_theta = dot(vec1, vec2) / (norm(vec1) * norm(vec2));

        if abs(cos_theta) > 1 - epsilon
            continue;
        end

        theta_i = acos(cos_theta);

        dUb_dtheta_i = -kb * (theta_i);

        nabla_ri_cos_theta = (vec1 / norm(vec1) - cos_theta * vec2 / norm(vec2)^2) / sqrt(1 - cos_theta^2);
        nabla_riplus1_cos_theta = (vec2 / norm(vec2) - cos_theta * vec1 / norm(vec1)^2) / sqrt(1 - cos_theta^2);
        nabla_riminus1_cos_theta = -nabla_ri_cos_theta - nabla_riplus1_cos_theta;

        dUb(i, :) = dUb(i, :) + dUb_dtheta_i * nabla_ri_cos_theta;
        dUb(i+1, :) = dUb(i+1, :) + dUb_dtheta_i * nabla_riplus1_cos_theta;
        dUb(i-1, :) = dUb(i-1, :) + dUb_dtheta_i * nabla_riminus1_cos_theta;
    end
    
    % Repulsion force calculations using neighbour lists
    for i = 1:total_monomers
        for j = neighbour_lists{i}'
            rij = positions(j,:) - positions(i,:);
            r = norm(rij);
            % Calculate forces only if within cutoff distance
            if r < cutoff_distance && r > 0 % Prevent division by zero
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
function positions = initialize_no_overlap(num_polymers, M, sigma)
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
        % Place subsequent monomers in a chain
        for m = 2:M
            positions((p-1) * M + m, :) = positions((p-1) * M + m - 1, :) + [sigma, 0];
        end
    end
end

%% Energy minimization function to remove overlaps and minimize potential energy
function positions = minimize_energy(positions, num_polymers, M, sigma, ks)
    total_monomers = num_polymers * M;
    for iteration = 1:10000
        F = zeros(total_monomers, 2);
        % Calculate repulsive forces between all pairs of monomers
        for i = 1:total_monomers
            for j = i+1:total_monomers
                rij = positions(j,:) - positions(i,:);
                r = norm(rij);
                % Apply forces if monomers are within a critical distance
                if r < 2 * sigma % Using 2 * sigma to consider the diameter of monomers
                    Fij = ks * (r - 2 * sigma) * (rij / r);
                    F(i,:) = F(i,:) - Fij;
                    F(j,:) = F(j,:) + Fij;
                end
            end
        end
        
        % Stop the minimization if forces are small
        if max(max(abs(F))) < 1e-3
            break;
        end
        
        % Update positions with a simple steepest descent approach
        positions = positions - 0.01 * F; 
    end
end

% Function to find neighbours of a given monomer
function neighbours = find_neighbours(monomer_idx, positions, cutoff_distance, total_monomers)
    distances = sqrt(sum((positions - positions(monomer_idx,:)).^2, 2));
    neighbours = find(distances < cutoff_distance & distances > 0 & (1:total_monomers)' ~= monomer_idx);
end
