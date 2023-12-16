
function main
    % Parameters
    noOfPointsInSolution = 6; % Define number of points in the solution path
    popSize = 10; % Set population size
    maxGenerations = 10;  % Set maximum number of generations (iterations)
    startPoint = [1, 1];   % Define the starting point of the path
    finishPoint = [500, 500]; % Define the end point of the path 
    tournamentSize = 2;  % Set the size for tournament selection
    mutationScale = 0.25; % Define scale for Gaussian mutation
    geneRange=[1,500];  % Define the range of gene values in the chromosome
    mutationRate = 0.1;  % Set the mutation rate
   
    % Read and Generate Map
    map = imbinarize(imread('random_map.bmp'));

    % Start timer
    tic;

    % User Input for Selection, Crossover, and Mutation Methods
    selectionMethod = input('Enter selection method (0 for RWS, 1 for Tournament, 2 for Rank-Based): ');
    crossoverMethod = input('Enter crossover method (0 for Method 1, 1 for Method 2): ');
    mutationMethod = input('Enter mutation method (0 for Method 1, 1 for Method 2): ');

    % Initialize Population
    population = initializePopulation(popSize, noOfPointsInSolution, startPoint, finishPoint);

    % Main GA Loop
    for gen = 1:maxGenerations
        % Fitness Evaluation
        fitness = evaluateFitness(population, map, startPoint, finishPoint);

        % Selection
        selected = selectPopulation(population, fitness, selectionMethod, tournamentSize);

        % Crossover
        offspring = crossover(selected, crossoverMethod);
       
        % Mutation
        mutatedOffspring = mutate(offspring, mutationMethod, mutationRate, geneRange, mutationScale);

        % Create New Generation
        population = [selected; mutatedOffspring];

        % Find the best path and its fitness in the current generation
        [bestPath, bestFitness] = findBestPath(population, fitness);
    end

    % Display Final Best Path
    displayFinalPath(map, bestPath, startPoint, finishPoint);

    % Stop timer and calculate execution time
    executionTime = toc;
    fprintf('Execution Time: %f seconds\n', executionTime); % Print the execution time

    % Calculate Total Euclidean Distance of the Optimal Path
    totalDistance = calculateTotalDistance(bestPath);
    fprintf('Total Distance: %f\n', totalDistance); % Print the total distance  
end

function population = initializePopulation(popSize, noOfPointsInSolution, startPoint, finishPoint)
    % Initialize an empty population matrix
    population = zeros(popSize, 2 * (noOfPointsInSolution + 2)); % +2 for start and finish points
   
    for i = 1:popSize
        % Generate random points within the map bounds for intermediate waypoints
        intermediatePoints = randi([1, 500], 1, 2 * noOfPointsInSolution);

        % Each row in population represents a path including start and finish points
        population(i, :) = [startPoint, intermediatePoints, finishPoint];
    end
end

% Evaluate the fitness of path to find the shortest path and avoid obstacles
function fitness = evaluateFitness(population, map, startPoint, finishPoint)
    numIndividuals = size(population, 1); % Get the number of individuals in the population
    fitness = zeros(numIndividuals, 1); % Initialize fitness array to zeros
    mapSize = size(map); % Retrieve the size of the map for bounds checking

    for i = 1:numIndividuals % Loop over each individual in the population
        % Create the path from the individual's genes
        waypoints = [startPoint; reshape(population(i, 1:end-2), 2, [])'; finishPoint];
        pathLength = 0; % Initialize path length to 0
        obstacleHits = 0; % Initialize obstacle hit count to 0
        totalAngleChange = 0; % Initialize total angle change to 0
        outOfBoundsPenalty = 0; % Initialize out-of-bounds penalty to 0

        for j = 1:(size(waypoints, 1) - 1) % Loop over each segment in the path
            % Calculate the Euclidean distance between waypoints
            segmentLength = norm(waypoints(j + 1, :) - waypoints(j, :));
            pathLength = pathLength + segmentLength; % Add segment length to total path length

            % Check for obstacle hits using Bresenham's line algorithm
            obstacleHits = obstacleHits + bresenham(waypoints(j, :), waypoints(j + 1, :), map);

            % Penalize path segments that go out of bounds
            if any(waypoints(j, :) < 1 | waypoints(j, :) > mapSize)
                outOfBoundsPenalty = outOfBoundsPenalty + 10000; % Add penalty for out-of-bounds
            end

            % Calculate angle change for path smoothness
            if j > 1 % Skip the first waypoint which doesn't have a previous segment
                prevVector = waypoints(j, :) - waypoints(j - 1, :); % Vector of previous segment
                currVector = waypoints(j + 1, :) - waypoints(j, :); % Vector of current segment
                % Calculate angle between vectors
                angle = acosd(dot(prevVector, currVector) / (norm(prevVector) * norm(currVector)));
                if ~isnan(angle) % Check for NaN which occurs if vectors are identical
                    totalAngleChange = totalAngleChange + angle; % Add angle to total angle change
                end
            end
        end

        % Calculate fitness value for the individual
        fitness(i) = 1 / (pathLength + obstacleHits * 1000000000 + totalAngleChange + outOfBoundsPenalty);
    end
end

% a function to help draw straightlines 
function hitCount = bresenham(p1, p2, map)
    hitCount = 0; % Initialize hit counter to 0
    x1 = round(p1(1)); % Round the first x-coordinate to the nearest integer
    y1 = round(p1(2)); % Round the first y-coordinate to the nearest integer
    x2 = round(p2(1)); % Round the second x-coordinate to the nearest integer
    y2 = round(p2(2)); % Round the second y-coordinate to the nearest integer

    % Clamp coordinates to ensure they are within the map boundaries
    x1 = max(1, min(x1, size(map, 2)));
    y1 = max(1, min(y1, size(map, 1)));
    x2 = max(1, min(x2, size(map, 2)));
    y2 = max(1, min(y2, size(map, 1)));

    dx = abs(x2 - x1); % Calculate the absolute difference in x
    dy = -abs(y2 - y1); % Calculate the negative absolute difference in y
    sx = sign(x2 - x1); % Determine the sign of the difference in x
    sy = sign(y2 - y1); % Determine the sign of the difference in y
    err = dx + dy; % Initialize the error term for the Bresenham algorithm

    % Bresenham's line algorithm
    while true
        % Check if the current coordinate is within the map bounds
        if x1 >= 1 && x1 <= size(map, 2) && y1 >= 1 && y1 <= size(map, 1)
            if map(y1, x1) == 0 % If the current coordinate hits an obstacle
                hitCount = hitCount + 1; % Increment hit counter
            end
        end

        % Check if the endpoint is reached
        if x1 == x2 && y1 == y2
            break; % Break the loop if the endpoint is reached
        end

        e2 = 2 * err; % Calculate doubled error term
        if e2 >= dy % Adjust error term and coordinate if necessary
            err = err + dy;
            x1 = x1 + sx;
        end
        if e2 <= dx % Adjust error term and coordinate if necessary
            err = err + dx;
            y1 = y1 + sy;
        end
    end
end

% user has an option between three different selection methods
function selected = selectPopulation(population, fitness, selectionMethod, tournamentSize)
    numIndividuals = size(population, 1);
    selected = zeros(size(population));
   
    if selectionMethod == 0 % Roulette Wheel Selection
      selected = rouletteWheelSelection(population, fitness);
    elseif selectionMethod == 1
        % Tournament Selection (to be implemented)
      selected = tournamentSelection(population, fitness, tournamentSize);

    elseif selectionMethod ==2   % Rank-Based Selection (to be implemented)
        selected = rankBasedSelection(population, fitness);

    end
end

% implementation of the roulette wheel selection method
function selected = rouletteWheelSelection(population, fitness)
    numIndividuals = size(population, 1); % Get the number of individuals in the population
    selected = zeros(size(population)); % Initialize a matrix to store the selected individuals

    cumulativeFitness = cumsum(fitness); % Compute the cumulative sum of fitness values
    totalFitness = cumulativeFitness(end); % Get the total fitness for normalization

    normalizedFitness = cumulativeFitness / totalFitness; % Normalize fitness values to sum to 1

    for i = 1:numIndividuals % Iterate over each individual
        r = rand; % Generate a random number between 0 and 1
        idx = binarySearch(normalizedFitness, r); % Find the index of the individual to select
        selected(i, :) = population(idx, :); % Add the selected individual to the selected matrix
    end
    fprintf('roulette chosen'); % Print a confirmation message
end

% we can use a binary search to produce optimal solution
function idx = binarySearch(cumulative, value)
    left = 1; % Initialize the left boundary of the search
    right = length(cumulative); % Initialize the right boundary of the search
    
    while left < right % Continue the loop until the search boundaries meet
        mid = left + floor((right - left) / 2); % Calculate the midpoint of the current search boundary
        if cumulative(mid) < value % If the value at midpoint is less than the target
            left = mid + 1; % Move the left boundary to the right
        else % Otherwise, the value at midpoint is greater than or equal to the target
            right = mid; % Move the right boundary to the left
        end
    end
    idx = left; % The target value is not less than the value at the left boundary
end

% implementation of the rankbased selection method
function selected = rankBasedSelection(population, fitness)
    numIndividuals = size(population, 1); % Get the number of individuals in the population
    selected = zeros(size(population)); % Initialize a matrix to store the selected individuals

    [~, sortedIndices] = sort(fitness, 'descend'); % Sort fitness in descending order and get indices
    ranks = 1:numIndividuals; % Assign ranks to each individual

    transformedRanks = 1 ./ ranks; % Transform ranks by taking their reciprocal
    selectionProbabilities = transformedRanks / sum(transformedRanks); % Normalize transformed ranks to sum to 1

    for i = 1:numIndividuals % Iterate over each individual
        idx = find(rand <= cumsum(selectionProbabilities), 1, 'first'); % Find the index to select based on rank probabilities
        selected(i, :) = population(sortedIndices(idx), :); % Add the selected individual to the selected matrix
    end
    fprintf('rank based selection works'); % Print a confirmation message
end

% implementation of the tournament selection method
function selected = tournamentSelection(population, fitness, tournamentSize)
    numIndividuals = size(population, 1); % Get the number of individuals in the population
    selected = zeros(size(population)); % Initialize a matrix to store the selected individuals

    for i = 1:numIndividuals % Iterate over each individual for selection
        indices = randperm(numIndividuals, tournamentSize); % Randomly select 'tournamentSize' individuals for the tournament
        [~, bestIdx] = max(fitness(indices)); % Find the index of the individual with the highest fitness in the tournament
        selected(i, :) = population(indices(bestIdx), :); % Add the fittest individual to the selected matrix
    end
    fprintf('tournament selection works\n'); % Print a confirmation message that the function has completed
end

% user has 2 options for crossover method twoPoint and uniform crossover
function offspring = crossover(selected, crossoverMethod)
    numPairs = size(selected, 1) / 2; % Calculate the number of pairs for crossover
    offspring = zeros(size(selected)); % Initialize offspring matrix
    
    for i = 1:numPairs % Iterate over each pair
        idx1 = 2*i-1; % Index of the first parent
        idx2 = 2*i; % Index of the second parent
        parent1 = selected(idx1, :); % Extract first parent
        parent2 = selected(idx2, :); % Extract second parent
        
        % Perform crossover based on the specified method
        switch crossoverMethod
            case 0
                offspringPair = twoPointCrossover(parent1, parent2); % Two-point crossover
            case 1
                offspringPair = uniformCrossover(parent1, parent2); % Uniform crossover
        end
        
        % Ensure the end point of each offspring is [500, 500]
        offspringPair(1, end-1:end) = [500, 500];
        offspringPair(2, end-1:end) = [500, 500];
        
        offspring(idx1:idx2, :) = offspringPair; % Place the offspring in the matrix
    end
end

% A method that ensures diversity in population
function offspring = uniformCrossover(parent1, parent2)
    chromosomeLength = length(parent1); % Get the length of the chromosome
    mask = rand(1, chromosomeLength) > 0.5; % Generate a mask for crossover

    % Create offspring by mixing genes based on the mask
    offspring1 = parent1;
    offspring2 = parent2;
    offspring1(mask) = parent2(mask);
    offspring2(mask) = parent1(mask);

    offspring = [offspring1; offspring2]; % Combine offspring into a single matrix
end

function offspring = twoPointCrossover(parent1, parent2)
    chromosomeLength = length(parent1); % Get the length of the chromosome
    points = sort(randperm(chromosomeLength-1, 2)); % Randomly choose two points for crossover
    crossoverPoint1 = points(1); % First crossover point
    crossoverPoint2 = points(2); % Second crossover point
    
    % Create offspring by combining segments from both parents
    offspring = [parent1(1:crossoverPoint1), parent2(crossoverPoint1+1:crossoverPoint2), parent1(crossoverPoint2+1:end);
                 parent2(1:crossoverPoint1), parent1(crossoverPoint1+1:crossoverPoint2), parent2(crossoverPoint2+1:end)];
end

% user has option between different mutation methods
function mutatedPopulation = mutate(population, mutationMethod, mutationRate, geneRange, mutationScale)
    numIndividuals = size(population, 1); % Determine the number of individuals in the population

    for i = 1:numIndividuals % Loop through each individual in the population
        % Check the mutation method to apply
        if mutationMethod == 0
            % Apply simple random mutation to the individual
            individual = simpleRandomMutation(population(i, :), mutationRate, geneRange);
        elseif mutationMethod == 1
            % Apply Gaussian mutation to the individual
            individual = gaussianMutation(population(i, :), mutationRate, mutationScale);
        end
        
        % Ensure the end point of each individual's path is [500,500]
        individual(end-1:end) = [500, 500];
        
        % Update the individual in the population with the mutated individual
        population(i, :) = individual;
    end

    % Return the mutated population
    mutatedPopulation = population;
end

% helps to maintain genetic diversity in a population
function mutated = simpleRandomMutation(individual, mutationRate, geneRange)
    mutationMask = rand(size(individual)) <= mutationRate;  % Create a mask for mutation
    randomValues = randi(geneRange, size(individual)); % Generate random values for mutation
    mutatedIndividual = individual; % Copy of the individual

    individual(mutationMask) = randomValues(mutationMask); % Apply mutation
     % Ensure mutatedIndividual stays within bounds
    mutatedIndividual = max(1, min(mutatedIndividual, 500));
    mutated = individual;
end

% mutation derived from Gaussian (normal) distribution
function mutated = gaussianMutation(individual, mutationRate, mutationScale)
    mutationMask = rand(size(individual)) <= mutationRate; % Create a mask for mutation
    gaussianPerturbations = mutationScale * randn(size(individual)); % Generate Gaussian perturbations
    mutatedIndividual = individual; % Copy of the individual
    individual(mutationMask) = individual(mutationMask) + gaussianPerturbations(mutationMask); % Apply mutation
    % Ensure mutatedIndividual stays within bounds
    mutatedIndividual = max(1, min(mutatedIndividual, 500));
    mutated = individual;
end

% we find the best path of the population based on fitness
function [bestPath, bestFitness] = findBestPath(population, fitness)
    [bestFitness, bestIdx] = max(fitness); % Find the best fitness and its index
    bestPath = population(bestIdx, :); % Extract the best path
end

% calculate euclidian distance
function totalDistance = calculateTotalDistance(path)
    totalDistance = 0; % Initialize total distance to zero

    % Check if path is a vector and reshape it if necessary
    if isvector(path)
        path = reshape(path, [], 2);
    end

    % Calculate the total distance of the path
    for i = 1:(size(path, 1) - 1)
        pointA = path(i, :); % Current point
        pointB = path(i + 1, :); % Next point
        totalDistance = totalDistance + sqrt(sum((pointB - pointA) .^ 2)); % Euclidean distance
    end
end

% display the final path to the screen
function displayFinalPath(map, solution, startPoint, finishPoint)
    % Convert solution to path
    path = [startPoint; reshape(solution, 2, [])'; finishPoint];

    % Clear current figure
    clf;

    % Show map
    imshow(map);
    hold on; % Hold on to plot the path on top of the map

    % Draw a rectangle around the map
    rectangle('position', [1, 1, size(map, 2) - 1, size(map, 1) - 1], 'edgecolor', 'k');

    % Plot the path
    line(path(:, 2), path(:, 1), 'Color', 'red', 'LineWidth', 2);

    hold off; % Release the hold
end
