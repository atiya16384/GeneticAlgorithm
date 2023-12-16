% MATLAB script to generate a binary image with circles, squares, rectangles, and triangles

% Set the size of the image
imageSize = 500;

% Create a blank white image
binaryImage = ones(imageSize);

% Number of random shapes
numShapes = 30;

% Generate random shapes
for i = 1:numShapes
    % Random shape type: 1 for rectangle, 2 for circle, 3 for triangle
    shapeType = randi([1, 3]);

    % Random position
    posX = randi([1, imageSize]);
    posY = randi([1, imageSize]);

    % Random size
    sizeX = randi([20, 50]); % Adjust the range as needed
    sizeY = randi([20, 50]); % Adjust the range as needed

    % Add the shape to the image
    if shapeType == 1
        % Rectangle
        binaryImage(posY:posY+sizeY, posX:posX+sizeX) = 0;
    elseif shapeType == 2
        % Circle
        [X, Y] = meshgrid(1:imageSize);
        mask = (X - posX).^2 + (Y - posY).^2 <= (sizeX / 2)^2;
        binaryImage(mask) = 0;
    else
        % Triangle
        triangle = poly2mask([posX, posX+sizeX, posX+sizeX/2], ...
                             [posY, posY, posY+sizeY], imageSize, imageSize);
        binaryImage(triangle) = 0;
    end
end

binaryImage = binaryImage(1:imageSize,1:imageSize);

% Ensure there is a clear path from (1,1) to (500,500)
binaryImage(imageSize-10:imageSize, 1:imageSize) = 1;
binaryImage(1:imageSize, 1:10) = 1;

% Display the generated image
imshow(binaryImage);
title('Binary Image with Random Shapes and a Clear Path');

% Save the image to a file if needed
imwrite(binaryImage, 'random_map.bmp');
