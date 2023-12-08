%% Worm like chain model
% We consider a wormlike chain, consisting of M active beads of 
% radius \Sigma at positions ri,  1; 2; ...;M. The potential energy 
% Upol= Us+ Ub+ Ur of the chain has three contributions.

% Equations of Motion:
% Overdamped Langevin equation for the position of each bead.
% Stochastic differential equation for the orientation of each bead.

%% Model Parameters
M = 20; % Number of active beads
sigma = 1; % Diameter of each monomer
ks = 100; % Spring constant for stretching
kb = 100; % Bending constant
kr = 100; % Repulsion constant
gamma = 1; % Damping coefficient
v0 = 1; % Active velocity
D = 1; % Diffusion coefficient
DR = 1; % Rotational diffusion coefficient



%% Initialize positions, velocities, and orientations

% To ensure that the monomers are closer together initially, we can modify
% the initialization of the positions such that they are closer to 
% each other. One way to do this is to initialize them along a line or a 
% slight curve, rather than randomly in the 2D space.

%In this initialization, the x positions of the monomers are 
% spaced by 0.1 units (which is the value of sigma), ensuring they
% are close together. The y positions have a slight random deviation to 
% give some randomness in the vertical direction.

% positions: This initializes the 2D positions of the M beads. Each bead has an (x, y) coordinate.
% theta: This initializes the orientation angles for the M beads. Each bead has an orientation in the range [0, 2*pi].

% Initialize positions along a line with slight randomness in y-direction.

% Function explainations：
% The cumsum function in MATLAB computes the cumulative sum of elements 
% along a given dimension of an array. The cumulative sum is the running 
% total of the sum of elements.

% Why we define positions?
% It initializes the positions matrix such that the x-coordinates of the 
% monomers are spaced approximately 0.1 units apart, and the y-coordinates
% have a slight random deviation. This ensures that the monomers are 
% initialized close to each other in a sort of slightly wiggly line.
positions = cumsum([sigma*ones(M,1), zeros(M,1)], 1);

%这个是缠在一起的，The orientation, represented by the variable theta, 
% describes the direction in which each monomer (or bead) of the polymer 
% is pointing. In a 2D plane, the orientation can be any angle in the range
% [0, 2*pi] radians, which corresponds to a full circle. 
% positions = rand(M, 2); % Random initial positions
theta = 2 * pi * rand(M, 1); % Random initial orientations


% Visualization setup
figure;
%axis([-M*sigma M*sigma -M*sigma M*sigma]);
%axis([0 2 -1 2]); %fixed rectangle boundary
axis equal
hold on;


% Video setup
video = VideoWriter('polymer_motion.avi');
%video.FrameRate = 15; 
open(video);


%% Time evolution
dt = 0.001; % Time step
numSteps = 1000;

% 加入MSD：Initialize a matrix to store positions at each time step
all_positions = zeros(M, 2, numSteps);

%% Calculated parameters
sigmasq = sigma * sigma;

for step = 1:numSteps
    % Clear for each frame except the last one
    if step < numSteps
        cla;
    end

     % Store positions at each time step
    all_positions(:, :, step) = positions;
%% Compute pairwise distances and differences

% pdist: Computes the pairwise distance between rows of a matrix.
% squareform: Converts a vector of distances into a square, symmetric matrix.
% reshape: Reshapes a matrix into a specified size.

% Mathematical Explanation:

%rij: This matrix provides the scalar distances between any pairs of 
%     monomers. It doesn't matter whether they are adjacent or not. 
%     It's a symmetric matrix where rij(i,j) gives the distance between 
%     monomer i and monomer j.

% dr: This matrix provides the direction and magnitude of the vector 
%     pointing from one monomer to another for any pair of monomers. 
%     It's not limited to adjacent monomers. For instance, dr(:,i,j) 
%     gives the vector pointing from monomer j to monomer i.



    dr = positions' - reshape(positions', [2, 1, M]);
    % imagesc(rij);

    %% Stretching force for each bead
    % The stretching force between two consecutive monomers is proportional to the difference between their current distance and the equilibrium distance sigma. The direction of the force is along the line connecting the two monomers.
    % The force on monomer i due to monomer i+1 is opposite in direction to the force on monomer i+1 due to monomer i.
    % F = -grad(Us);
    dUs = zeros(M, 2);
    for i = 1:M-1
        direction = (positions(i+1, :) - positions(i, :)) / norm(positions(i+1, :) - positions(i, :));
        dUs(i,   :) = dUs(i,   :) + ks * (norm(positions(i+1, :) - positions(i, :)) - sigma) * direction;
        dUs(i+1, :) = dUs(i+1, :) - ks * (norm(positions(i+1, :) - positions(i, :)) - sigma) * direction;
    end
  

%% Bending energy gradient with stability check
epsilon = 1e-6;  % Threshold to avoid division by zero
dUb = zeros(M, 2);  % Initialize bending energy gradient

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


    
    %% Repulsion energy gradient
    % The squeeze function in MATLAB is used to remove singleton dimensions 
    % from an array. A singleton dimension is a dimension for which the 
    % size is 1.

    % dr is a 3D matrix of shape [2, M, M], where M is the number of 
    % monomers. dr(:, i, j) is a 2x1 vector representing the difference 
    % in positions between monomer j and monomer i.

    % When you use squeeze(dr(:, i, j)), it removes the singleton 
    % dimensions and returns a 2-element vector. Essentially, 
    % it converts the 2x1 matrix into a simple 2-element row vector.

    %the entire expression Kr * ... calculates the repulsion force
    % contribution from monomer j to monomer i based on their distance 
    % and direction. This force is then added to the existing dUr(i, :),
    % which accumulates the repulsion forces from all other monomers
    % on monomer i.
    dUr = zeros(M, 2);
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
    
    % Total force
    F = dUs + dUb + dUr;
   

    %% Update positions using overdamped Langevin equation
    n = [cos(theta), sin(theta)];
    positions = positions + dt * (F/gamma + v0 * n) + sqrt(2*D) * randn(M, 2) * sqrt(dt);
   
    %% Update orientations
    xi = randn(M, 1);
    theta = theta + sqrt(2*DR) * xi * sqrt(dt);
    

    %% Visualization
    if(~mod(step,10))
        h = plot(positions(:,1), positions(:,2), '-');
        for iM=1:M
            h(iM+1) = rectangle('position',[positions(iM,1)-sigma/2 positions(iM,2)-sigma/2 sigma sigma],'curvature',[1 1]);
        end
        title(['Step: ', num2str(step)]);
        drawnow;
        
        % Save frame
        frame = getframe(gcf);
        writeVideo(video, frame);

        % pause();
    
        % Clear for next frame
        %cla;
    end
end

close(video);
