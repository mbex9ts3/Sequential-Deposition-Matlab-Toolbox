function [centroids, radii, index] = drop_roll_zones(maskIn, distMod, param1, param2, minDg, maxDg, vizFlag)
% Modification of the classic 'drop and roll' algorithm for layered 
% granular media generation with conditional grain size 

% Takes container bounds and particle size distribution as input and
% returns the centroids and radii of the generated particles

% Collision detection solutions are close form for speed
% Circle circle intersections solved with a modified (vectorized) version
% of the native Matlab circcirc.m function for speed

% INPUTS: 
% y_bounds - container y boundaries [ymin ymax]
% x_bounds - container x boundaries [xmin xmax]
% layers - nx1 vector containing max layer heights (termination criteria)
%        - values should grow sequentially and not exceed y_bounds(2)
% distMod - nx1 cell of strings defining distributional forms associated 
%         - each layer. Supported distributions:
%         - 'Normal' (param1 = mean / param2 = standard deviation)
%         - 'Uniform' (param1 = minimum / param2 = maximum)
%         - 'Exponential' (param1 = mean / param2 = [])
%         - 'Lognormal' (param1 = mean / param2 = standard deviation)
%         - 'Weibull' (param1 = scale parameter / param2 = shape parameter)
%         - 'Unimodal' (param1 = value / param2 = nan)
% param1 - nx1 distribution parameter 1 (see above for definitions)
% param2 - nx1 distribution parameter 2 (see above for definitions)
% param1_sc - nx1 scaling factor for param1 (vertical gradients)
%           - set as 1 for uniform parameters in each layer
%           - applied sequentially / set to 1 if 
% param2_sc - nx1 scaling factor for param2 (vertical gradients)
%           - set as 1 for uniform parameters in each layer
% minDg - nx1 minimuum grainsize for each layer
% maxDg - nx1 maximum grainsize for each layer

% OUTPUTS: 
% centroids - generated circle centroids (nx2)
% radii - generated circle radii (nx1)

% Author: Thomas Seers: thomas.seers@qatar.tamu.edu (2020)

% setup container (image limts)
x_bounds = [1 size(maskIn, 2)];
y_bounds = [1 size(maskIn, 1)];

% for checks on zonation
if min(min(maskIn)) == false
    maskTest = maskIn + 1;
else 
    maskTest = maskIn;
end

% initialize
nZones = size(unique(maskIn),1);
idx = randperm(nZones, 1);
distroIn = distMod{idx};
param1In = param1(idx);
param2In = param2(idx);
param1_sc = 1;
param2_sc = 1;
maxDgIn = maxDg(idx);
minDgIn = minDg(idx);
if strcmp(distroIn, 'Exponential') == true
    expoFlag = true;
else
    expoFlag = false;
end

if strcmp(distroIn, 'Unimodal') == true
    uniFlag = true;
else
    uniFlag = false;
end

% expo flag (one input)
if expoFlag == true
    circIn = random(distroIn, param1In);
else if uniFlag == true
        circIn = param1In;
    else
        circIn = random(distroIn, param1In, param2In);
    end
end

if uniFlag == false
    while circIn < minDgIn || circIn > maxDgIn
        if expoFlag == true
            circIn = random(distroIn, param1In);
        else
            circIn = random(distroIn, param1In, param2In);
        end
    end
end
circR = circIn/2;

% x cent
x_cent = rand*x_bounds(2);
while x_cent < x_bounds(1) + circR || x_cent > x_bounds(2) - circR
    x_cent = rand*x_bounds(2);
end

circC_T = [x_cent y_bounds(1) + circR circR]; % 'drop'
circB_T = [x_cent - circR x_cent + circR]; % bounds in x

% test if circ is in the correct region
currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
index = [];
if isequal(currentVal, idx)
    
    % append to list
    circC = circC_T;
    circB = circB_T;
    index = vertcat(index, currentVal);
else
    while ~isequal(currentVal, idx)
        % expo flag (one input)
        if expoFlag == true
            circIn = random(distroIn, param1In);
        else if uniFlag == true
                circIn = param1In;
            else
                circIn = random(distroIn, param1In, param2In);
            end
        end
        
        if uniFlag == false
            while circIn < minDgIn || circIn > maxDgIn
                if expoFlag == true
                    circIn = random(distroIn, param1In);
                else
                    circIn = random(distroIn, param1In, param2In);
                end
            end
        end
        circR = circIn/2;
        
        % x cent
        x_cent = rand*x_bounds(2);
        while x_cent < x_bounds(1) + circR || x_cent > x_bounds(2) - circR
            x_cent = rand*x_bounds(2);
        end
        
        circC_T = [x_cent y_bounds(1) + circR circR]; % 'drop'
        circB_T = [x_cent - circR x_cent + circR]; % bounds in x
        
        % test if circ is in the correct region
        currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
        
        if isequal(currentVal, idx)
            
            % append to list
            circC = circC_T;
            circB = circB_T;
            index = vertcat(index, currentVal);
        end
    end
end

particleCount = 0;
count = 0;
nZones = size(unique(maskIn),1);
% set termination criteria max(circC(:,2)) > y_bounds(2)
% for i = 1:50000
rejectN = 0; % count consecutive rejected cases to identify degenerate geometry /stuck cases
while max(circC(:,2)) < y_bounds(2)
% for i = 1:9000
    count = count+1;
    particleCount = particleCount+1;
    disp(['fitting particle ' num2str(size(circC,1))]);
    disp(['current total particle count =  ' num2str(particleCount)]);
    yFlag = false; % to test if centroid falls outside the mask in the loop
    
    % initialize
    idx = randperm(nZones, 1);
    distroIn = distMod{idx};
    param1In = param1(idx);
    param2In = param2(idx);
    param1_scIn = 1;
    param2_scIn = 1;
    maxDgIn = maxDg(idx);
    minDgIn = minDg(idx);
    
    if strcmp(distroIn, 'Exponential') == true
        expoFlag = true;
    else
        expoFlag = false;
    end
    
    if strcmp(distroIn, 'Unimodal') == true
        uniFlag = true;
    else
        uniFlag = false;
    end
    
    % apply scaling
    if param1In > minDgIn && param1In < maxDgIn
        param1In = param1In*param1_scIn;
        param2In = param2In*param2_scIn;
    end
    
    % expo flag (one input)
    if expoFlag == true
        circIn = random(distroIn, param1In);
    else if uniFlag == true
            circIn = param1In;
        else
            circIn = random(distroIn, param1In, param2In);
        end
    end
    
    if uniFlag == false
        while circIn < minDgIn || circIn > maxDgIn
            if expoFlag == true
                circIn = random(distroIn, param1In);
            else
                circIn = random(distroIn, param1In, param2In);
            end
        end
    end
    circR = circIn/2;
    
    % x cent
    x_cent = rand*x_bounds(2);
    while x_cent < x_bounds(1) + circR || x_cent > x_bounds(2) - circR
        x_cent = rand*x_bounds(2);
    end
    
    % test if any circs are in the path of x_cent +/- circR
    x_path = [x_cent - circR x_cent + circR];
    logi = true(size(circC, 1),1);
    logi(circB(:,2) < x_path(:,1)) = false; % lower
    logi(circB(:,1) > x_path(:,2)) = false; % lower
    
    if sum(logi) == 0 % no objects in the path
        circC_T = [x_cent y_bounds(1) + circR circR]; % 'drop'
        circB_T = [x_cent - circR x_cent + circR]; % bounds in x
        
        % test if circ is in the correct region
        if rejectN < 10
            currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
        else % reduce acceptence criteria to any orthogonal edge vertices
            
            % build mask
            ang=0:0.1:2*pi;
            xp=circR*sin(ang);
            yp=circR*cos(ang);
            poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
            BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
            maskT = maskIn;
            maskT(BW == false) = false;
            if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                currentVal = idx;
            end
        end
        
        if isequal(currentVal, idx)
            
            % append to list
            circC = vertcat(circC, circC_T);
            circB = vertcat(circB, circB_T);
            index = vertcat(index, currentVal);
            rejectN = 0;
           
            % test if
            stopFlag = true; % stop criteria met
        else
            stopFlag = nan; % resting particle rejected
            disp('particle outside the ROI');
            rejectN = rejectN+1;
        end
        
    else
        % initialize stopFlag
        stopFlag = false;

        % find first intersecting circle along V = [0 -1]
        circsTest = circC(logi,:);
        
        % Solve quadratic equation for minimum distance to the intersected 
        % circs: d (possibly zero, one, or two solutions)
        % px^2 + (py - d)^2 = (r1 + r2)^2
        % (py - d)^2 = (r1 + r2)^2 - px^2
        % d = py +/- Sqrt((r1 + r2)^2 - px^2)
        r1 = repmat(circR, size(circsTest,1),1);
        r2 = circC(logi ,3);
        p1x = repmat(x_cent, size(circsTest,1),1);
        p1y = repmat(y_bounds(2), size(circsTest,1),1);
        px = p1x - circC(logi,1);
        py = p1y - circC(logi,2);
        d1 = py + sqrt((r1 + r2).^2 - px.^2);
        d2 = py - sqrt((r1 + r2).^2 - px.^2);
        
        % collect  results
        minD = min([d1 d2], [], 2);
        [mV I] = min(minD);
        
        % select intersecting circle
        circHit = circsTest(I,:);
        
        % translate             
        circC_T = [x_cent y_bounds(2) - minD(I) circR]; % 'drop'
        
        % initialize rotation around circHit ('roll')
        % nb 3 post rotation termination conditions exits:
        % (1) circC_T is stranded, whereby its centroid is between the
        %     centroids of two adjacent circles
        % (2) circC_T comes to rest against the side wall of the container
        % (3) circC_T comes to rest against the base of the container
        
        % establish sense of rotation - if p1 > p2; R is clockwise; else
        % R is counter clockwise
        if circC_T(1) < circHit(1) % anticlockwise
            Rflag = false;
        else 
            Rflag = true; % clockwise
        end
        
        % check for container wall intersections (test)
        r3 = circC_T(3)+circHit(3);       % check wall (dependant upon Rflag)
        if Rflag == false % check left
            Mtan = x_bounds(1) + circC_T(3);
            [xoutW,youtW] = linecirc(inf,Mtan,circHit(1),circHit(2),r3);
        else
            Mtan = x_bounds(2) - circC_T(3); % check right
            [xoutW,youtW] = linecirc(inf,Mtan,circHit(1),circHit(2),r3);
        end
        if size(xoutW,1) > 1 % if tangent
            xoutW = transpose(xoutW);
            youtW = transpose(youtW);
        end
        
        % remove lower contact
        [Wmax, WI] = max(youtW);
        xoutW = xoutW(WI);
        youtW = youtW(WI);
        
        % check for container basal intersection
        Mtan = y_bounds(1) + circC_T(3);
        [xoutB,youtB] = linecirc(0,Mtan,circHit(1),circHit(2),r3);
        if size(xoutB,1) > 1 % if tangent
            xoutB = transpose(xoutB);
            youtB = transpose(youtB);
        end
        
        % remove lower contact
        [Bmax, maxI] = max(xoutB);
        [Bmax, minI] = min(xoutB);

        if Rflag == false % check left side: remove right intersection
            xoutB = xoutB(minI);
            youtB = youtB(minI);
        else % check right side: remove left intersction
            xoutB = xoutB(maxI);
            youtB = youtB(maxI);
        end
                    
        % locate all collisions on the arc path around circHit
        % find potential hits
        x1 = repmat(circHit(1),size(circC,1),1);
        x2 = circC(:,1);
        y1 = repmat(circHit(2),size(circC,1),1);
        y2 = circC(:,2);
        r1 = repmat(circHit(3)+circC_T(3),size(circC,1),1);
        r2 = circC(:,3)+repmat(circC_T(3),size(circC,1),1);
        
        % vectorized implementation of circcirc.m
        [xout, yout] = circcircvect(x1,y1,r1,x2,y2,r2); % problem
        
        % remove hits on along the opposite arc-path
        if Rflag == true % clockwise
            yout(xout < circC_T(1)) = nan;
            xout(xout < circC_T(1)) = nan;
        else % anticlockwise
            yout(xout > circC_T(1)) = nan;
            xout(xout > circC_T(1)) = nan;
        end
        
        %test all hits to find the highest centroid
        int1 = [xout(:,1) yout(:,1)];
        int2 = [xout(:,2) yout(:,2)];

        % determine highest centroid contact in the container (if exist)
        contain = vertcat([xoutW youtW], [xoutB,youtB]);
        if sum(sum(~isnan(contain)))
            % collision with the container detected
            [mCont, iCont] = max(contain(:,2)); % highest centroid
        else
            iCont = nan;
        end
                
        % circle test
        if sum(sum(isnan(int1))) < size(int1,1)*size(int1,2)  || sum(sum(isnan(int2))) < size(int2,1)*size(int2,2) % circle intersections located
           [mCirc1, iCirc1] = max(int1(:,2)); % highest centroid
           [mCirc2, iCirc2] = max(int2(:,2)); % highest centroid
        else
            iCirc1 = nan;
            iCirc2 = nan;
        end
        
        % determine system state
        if isnan(iCont) % container not detected
            if isnan(iCirc1) && isnan(iCirc2) % degenerate case
                disp('degenerate case');
                continue
            else if mCirc1 > mCirc2 || isnan(mCirc2) % first centroid is higher
                    circC_T = [int1(iCirc1,1) int1(iCirc1,2) circR]; % 'roll'
                    circB_T = [int1(iCirc1,1) - circR int1(iCirc1,1) + circR]; % bounds in x
                    iCirc = iCirc1;
                else
                    circC_T = [int2(iCirc2,1) int2(iCirc2,2) circR]; % 'roll'
                    circB_T = [int2(iCirc2,1) - circR int2(iCirc2,1) + circR]; % bounds in x
                    iCirc = iCirc2;
                end
            end
            
            % determine if stop condition (local depression) is met
            if circHit(1) < circC(iCirc,1)
                xMinTest = circHit(1);
                xMaxTest = circC(iCirc,1);
            else
                xMinTest = circC(iCirc,1);
                xMaxTest = circHit(1);
            end
            if circC_T(1) > xMinTest && circC_T(1) < xMaxTest % new position lies between centroids
                
                % check the particle is not out of bounds
                if circC_T(2) > y_bounds(2)
                    break
                end
                % test if circ is in the correct region
                if rejectN < 10
                    currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
                else % reduce acceptence criteria to any orthogonal edge vertices
                    
                    % build mask
                    ang=0:0.1:2*pi;
                    xp=circR*sin(ang);
                    yp=circR*cos(ang);
                    poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
                    BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
                    maskT = maskIn;
                    maskT(BW == false) = false;
                    if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                        currentVal = idx;
                    end
                end

                
                if isequal(currentVal, idx)
                    
                    % append to list
                    circC = vertcat(circC, circC_T);
                    circB = vertcat(circB, circB_T);
                    index = vertcat(index, currentVal);
                    rejectN = 0;

                    % test if
                    stopFlag = true; % stop criteria met
                else
                    stopFlag = nan; % resting particle rejected
                    disp('particle outside the ROI');
                    rejectN = rejectN+1;
                end
            end
            
        else if isnan(iCirc1) && isnan(iCirc2) % no circles in the arc-path
                circC_T = [contain(iCont,:) circR]; % 'roll'
                circB_T = [circC_T(1) - circR circC_T(1) + circR]; % bounds in x
                
                % check the particle is not out of bounds
                if circC_T(2) > y_bounds(2)
                    break
                end
                
                % test if circ is in the correct region
                if rejectN < 10
                    currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
                else % reduce acceptence criteria to any orthogonal edge vertices
                    
                    % build mask
                    ang=0:0.1:2*pi;
                    xp=circR*sin(ang);
                    yp=circR*cos(ang);
                    poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
                    BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
                    maskT = maskIn;
                    maskT(BW == false) = false;
                    if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                        currentVal = idx;
                    end
                end
                
                if isequal(currentVal, idx)
                    
                    % append to list
                    circC = vertcat(circC, circC_T);
                    circB = vertcat(circB, circB_T);
                    index = vertcat(index, currentVal);
                    rejectN = 0;
                    
                    % test if
                    stopFlag = true; % stop criteria met
                else
                    stopFlag = nan; % resting particle rejected
                    disp('particle outside the ROI');
                    rejectN = rejectN+1;

                end
                
            else % both container and circles lie in the arc-path
                % determine whether circ or container gives highest
                % centroid
                [mAll, iAll] = max([mCont, mCirc1, mCirc2]);
                if iAll == 1 % container is first collision
                    circC_T = [contain(iCont,:) circR]; % 'roll'
                    circB_T = [circC_T(1) - circR circC_T(1) + circR]; % bounds in x
                    
                    % check the particle is not out of bounds
                    if circC_T(2) > y_bounds(2)
                        break
                    end
                    
                    % test if circ is in the correct region
                    if rejectN < 10
                        currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
                    else % reduce acceptence criteria to any orthogonal edge vertices
                        
                        % build mask
                        ang=0:0.1:2*pi;
                        xp=circR*sin(ang);
                        yp=circR*cos(ang);
                        poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
                        BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
                        maskT = maskIn;
                        maskT(BW == false) = false;
                        if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                            currentVal = idx;
                        end
                    end

                    
                    if isequal(currentVal, idx)
                        
                        % append to list
                        circC = vertcat(circC, circC_T);
                        circB = vertcat(circB, circB_T);
                        index = vertcat(index, currentVal);
                        rejectN = 0;

                        % test if
                        stopFlag = true; % stop criteria met
                    else
                        stopFlag = nan; % resting particle rejected
                        disp('particle outside the ROI');
                        rejectN = rejectN+1;
                    end
                    
                else % circle is the first collision
                    if iAll == 2
                        circC_T = [int1(iCirc1,1) int1(iCirc1,2) circR]; % 'roll'
                        circB_T = [int1(iCirc1,1) - circR int1(iCirc1,1) + circR]; % bounds in x
                        iCirc = iCirc1;
                    else
                        circC_T = [int2(iCirc2,1) int2(iCirc2,2) circR]; % 'roll'
                        circB_T = [int2(iCirc2,1) - circR int2(iCirc2,1) + circR]; % bounds in x
                        iCirc = iCirc2;
                    end
                    
                    % determine if stop condition (local depression) is met
                    if circHit(1) < circC(iCirc,1)
                        xMinTest = circHit(1);
                        xMaxTest = circC(iCirc,1);
                    else
                        xMinTest = circC(iCirc,1);
                        xMaxTest = circHit(1);
                    end
                    if circC_T(1) > xMinTest && circC_T(1) < xMaxTest % new position lies between centroids
                        
                        % check the particle is not out of bounds
                        if circC_T(2) > y_bounds(2)
                            break
                        end
                        
                        % test if circ is in the correct region
                        if rejectN < 10
                            currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
                        else % reduce acceptence criteria to any orthogonal edge vertices
                            
                            % build mask
                            ang=0:0.1:2*pi;
                            xp=circR*sin(ang);
                            yp=circR*cos(ang);
                            poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
                            BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
                            maskT = maskIn;
                            maskT(BW == false) = false;
                            if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                                currentVal = idx;
                            end
                        end
                        
                        if isequal(currentVal, idx)
                            
                            % append to list
                            circC = vertcat(circC, circC_T);
                            circB = vertcat(circB, circB_T);
                            index = vertcat(index, currentVal);
                            rejectN = 0;

                            % test if
                            stopFlag = true; % stop criteria met
                        else
                            stopFlag = nan; % resting particle rejected
                            disp('particle outside the ROI');
                            rejectN = rejectN+1;

                        end
                    end
                end
            end
        end
    end

    % test if stop criteria is met: break or conditional loop
    if stopFlag == true || isnan(stopFlag)
        continue
    else
        circHit = circC(iCirc,:); % new rotational centroid and radius
        
        % iterate until stable location is found
        while stopFlag == false
            
            % initialize rotation around circHit ('roll')
            % nb 3 post rotation termination conditions exits:
            % (1) circC_T is stranded, whereby its centroid is between the
            %     centroids of two adjacent circles
            % (2) circC_T comes to rest against the side wall of the container
            % (3) circC_T comes to rest against the base of the container
            
            
            % establish sense of rotation - if p1 > p2; R is clockwise; else
            % R is counter clockwise
            if circC_T(1) < circHit(1) % anticlockwise
                Rflag = false;
            else
                Rflag = true; % clockwise
            end
        
            % check for container wall intersections (test)
            r3 = circC_T(3)+circHit(3);       % check wall (dependant upon Rflag)
            if Rflag == false % check left
                Mtan = x_bounds(1) + circC_T(3);
                [xoutW,youtW] = linecirc(inf,Mtan,circHit(1),circHit(2),r3);
            else
                Mtan = x_bounds(2) - circC_T(3); % check right
                [xoutW,youtW] = linecirc(inf,Mtan,circHit(1),circHit(2),r3);
            end
            if size(xoutW,1) > 1 % if tangent
                xoutW = transpose(xoutW);
                youtW = transpose(youtW);
            end
            
            % remove lower contact
            [Wmax, WI] = max(youtW);
            xoutW = xoutW(WI);
            youtW = youtW(WI);
            
            % check for container basal intersection
            Mtan = y_bounds(1) + circC_T(3);
            [xoutB,youtB] = linecirc(0,Mtan,circHit(1),circHit(2),r3);
            if size(xoutB,1) > 1 % if tangent
                xoutB = transpose(xoutB);
                youtB = transpose(youtB);
            end
            
            % remove lower contact
            [Bmax, maxI] = max(xoutB);
            [Bmax, minI] = min(xoutB);
            
            if Rflag == false % check left side: remove right intersection
                xoutB = xoutB(minI);
                youtB = youtB(minI);
            else % check right side: remove left intersction
                xoutB = xoutB(maxI);
                youtB = youtB(maxI);
            end
            
            % locate all collisions on the arc path around circHit
            % find potential hits
            x1 = repmat(circHit(1),size(circC,1),1);
            x2 = circC(:,1);
            y1 = repmat(circHit(2),size(circC,1),1);
            y2 = circC(:,2);
            r1 = repmat(circHit(3)+circC_T(3),size(circC,1),1);
            r2 = circC(:,3)+repmat(circC_T(3),size(circC,1),1);
            
            % vectorized implementation of circcirc.m
            [xout, yout] = circcircvect(x1,y1,r1,x2,y2,r2);
            
            % remove hits on along the opposite arc-path
            if Rflag == true % clockwise
                yout(xout < circC_T(1)) = nan;
                xout(xout < circC_T(1)) = nan;
            else % anticlockwise
                yout(xout > circC_T(1)) = nan;
                xout(xout > circC_T(1)) = nan;
            end

            %test all hits to find the highest centroid
            int1 = [xout(:,1) yout(:,1)];
            int2 = [xout(:,2) yout(:,2)];
                        
            % check for circC_T centroid and cut
            tol = 0.0000000001; % numerical error set to 1x10^-10
            Lia1 = ismembertol(int1,circC_T(1:2), tol, 'ByRows',true);
            Lia2 = ismembertol(int2,circC_T(1:2), tol, 'ByRows',true);
            int1(Lia1,:) = nan;
            int2(Lia2,:) = nan;
            
            % determine highest centroid contact in the container (if exist)
            contain = vertcat([xoutW youtW], [xoutB,youtB]);
            if sum(sum(~isnan(contain)))
                % collision with the container detected
                [mCont, iCont] = max(contain(:,2)); % highest centroid
            else
                iCont = nan;
            end
            
            % circle test
            if sum(sum(isnan(int1))) < size(int1,1)*size(int1,2)  || sum(sum(isnan(int2))) < size(int2,1)*size(int2,2) % circle intersections located
                [mCirc1, iCirc1] = max(int1(:,2)); % highest centroid
                [mCirc2, iCirc2] = max(int2(:,2)); % highest centroid
            else
                iCirc1 = nan;
                iCirc2 = nan;
            end
            
            % determine system state
            if isnan(iCont) % container not detected
                if isnan(iCirc1) && isnan(iCirc2)
                    break
                else if mCirc1 > mCirc2 || isnan(mCirc2) % first centroid is higher
                        circC_T = [int1(iCirc1,1) int1(iCirc1,2) circR]; % 'roll'
                        circB_T = [int1(iCirc1,1) - circR int1(iCirc1,1) + circR]; % bounds in x
                        iCirc = iCirc1;
                    else
                        circC_T = [int2(iCirc2,1) int2(iCirc2,2) circR]; % 'roll'
                        circB_T = [int2(iCirc2,1) - circR int2(iCirc2,1) + circR]; % bounds in x
                        iCirc = iCirc2;
                    end
                end
                
                % determine if stop condition (local depression) is met
                if circHit(1) < circC(iCirc,1)
                    xMinTest = circHit(1);
                    xMaxTest = circC(iCirc,1);
                else
                    xMinTest = circC(iCirc,1);
                    xMaxTest = circHit(1);
                end
                
                if circC_T(1) > xMinTest && circC_T(1) < xMaxTest % new position lies between centroids
                      
                    % check the particle is not out of bounds
                    if circC_T(2) > y_bounds(2)
                        yFlag = true;
                        break
                    end
                    
                    % test if circ is in the correct region
                    if rejectN < 10
                        currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
                    else % reduce acceptence criteria to any orthogonal edge vertices
                        
                        % build mask
                        ang=0:0.1:2*pi;
                        xp=circR*sin(ang);
                        yp=circR*cos(ang);
                        poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
                        BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
                        maskT = maskIn;
                        maskT(BW == false) = false;
                        if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                            currentVal = idx;
                        end
                    end
                    
                    
                    if isequal(currentVal, idx)
                        
                        % append to list
                        circC = vertcat(circC, circC_T);
                        circB = vertcat(circB, circB_T);
                        index = vertcat(index, currentVal);
                        rejectN = 0;

                        % test if
                        stopFlag = true; % stop criteria met
                    else
                        stopFlag = nan; % resting particle rejected
                        disp('particle outside the ROI');
                        rejectN = rejectN+1;
                    end
                end
                
            else if isnan(iCirc1) && isnan(iCirc2) % no circles in the arc-path
                    circC_T = [contain(iCont,:) circR]; % 'roll'
                    circB_T = [circC_T(1) - circR circC_T(1) + circR]; % bounds in x
                    
                    % check the particle is not out of bounds
                    if circC_T(2) > y_bounds(2)
                        yFlag = true;
                        break
                    end
                    
                    % test if circ is in the correct region
                    if rejectN < 10
                        currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
                    else % reduce acceptence criteria to any orthogonal edge vertices
                        
                        % build mask
                        ang=0:0.1:2*pi;
                        xp=circR*sin(ang);
                        yp=circR*cos(ang);
                        poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
                        BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
                        maskT = maskIn;
                        maskT(BW == false) = false;
                        if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                            currentVal = idx;
                        end
                    end
                        
                    
                    if isequal(currentVal, idx)
                        
                        % append to list
                        circC = vertcat(circC, circC_T);
                        circB = vertcat(circB, circB_T);
                        index = vertcat(index, currentVal);
                        rejectN = 0;

                        % test if
                        stopFlag = true; % stop criteria met
                    else
                        stopFlag = nan; % resting particle rejected
                        disp('particle outside the ROI');
                        rejectN = rejectN+1;
                    end
                    
                else % both container and circles lie in the arc-path
                    % determine whether circ or container gives highest
                    % centroid
                    [mAll, iAll] = max([mCont, mCirc1, mCirc2]);
                    if iAll == 1 % container is first collision
                        circC_T = [contain(iCont,:) circR]; % 'roll'
                        circB_T = [circC_T(1) - circR circC_T(1) + circR]; % bounds in x
                        
                        % check the particle is not out of bounds
                        if circC_T(2) > y_bounds(2)
                            yFlag = true;
                            break
                        end
                        
                        % test if circ is in the correct region
                        if rejectN < 10
                            currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
                        else % reduce acceptence criteria to any orthogonal edge vertices
                            
                            % build mask
                            ang=0:0.1:2*pi;
                            xp=circR*sin(ang);
                            yp=circR*cos(ang);
                            poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
                            BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
                            maskT = maskIn;
                            maskT(BW == false) = false;
                            if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                                currentVal = idx;
                            end
                        end
                        
                        if isequal(currentVal, idx)
                            
                            % append to list
                            circC = vertcat(circC, circC_T);
                            circB = vertcat(circB, circB_T);
                            index = vertcat(index, currentVal);
                            rejectN = 0;

                            % test if
                            stopFlag = true; % stop criteria met
                        else
                            stopFlag = nan; % resting particle rejected
                            disp('particle outside the ROI');
                            rejectN = rejectN+1;
                        end
                        
                    else % circle is the first collision
                        if iAll == 2
                            circC_T = [int1(iCirc1,1) int1(iCirc1,2) circR]; % 'roll'
                            circB_T = [int1(iCirc1,1) - circR int1(iCirc1,1) + circR]; % bounds in x
                            iCirc = iCirc1;
                        else
                            circC_T = [int2(iCirc2,1) int2(iCirc2,2) circR]; % 'roll'
                            circB_T = [int2(iCirc2,1) - circR int2(iCirc2,1) + circR]; % bounds in x
                            iCirc = iCirc2;
                        end
                        
                        % determine if stop condition (local depression) is met
                        if circHit(1) < circC(iCirc,1)
                            xMinTest = circHit(1);
                            xMaxTest = circC(iCirc,1);
                        else
                            xMinTest = circC(iCirc,1);
                            xMaxTest = circHit(1);
                        end
                        
                        if circC_T(1) > xMinTest && circC_T(1) < xMaxTest % new position lies between centroids
                            
                            % check the particle is not out of bounds
                            if circC_T(2) > y_bounds(2)
                                yFlag = true;
                                break
                            end
                            
                            % test if circ is in the correct region
                            if rejectN < 10
                                currentVal = maskTest(round(circC_T(2)), round(circC_T(1)));
                            else % reduce acceptence criteria to any orthogonal edge vertices
                                
                                % build mask
                                ang=0:0.1:2*pi;
                                xp=circR*sin(ang);
                                yp=circR*cos(ang);
                                poly = [transpose(circC_T(2)+xp) transpose(circC_T(1)+yp)];
                                BW = poly2mask(poly(:,1),poly(:,2),size(maskIn, 1),size(maskIn, 2));
                                maskT = maskIn;
                                maskT(BW == false) = false;
                                if size(find(maskT == idx), 1) > sum(sum(BW))*0.3333 % accept if 1/3 overlap
                                    currentVal = idx;
                                end
                            end
                        
                            if isequal(currentVal, idx)
                                
                                % append to list
                                circC = vertcat(circC, circC_T);
                                circB = vertcat(circB, circB_T);
                                index = vertcat(index, currentVal);
                                rejectN = 0;

                                % test if
                                stopFlag = true; % stop criteria met
                            else
                                stopFlag = nan; % resting particle rejected
                                disp('particle outside the ROI');
                                rejectN = rejectN+1;
                            end
                        end
                    end
                end
            end
            
            % test if stop criteria is met: break or conditional loop
            if stopFlag == true || isnan(stopFlag)
                continue
            else
                circHit = circC(iCirc,:); % new rotational centroid and radius
            end
            
        end
        if yFlag == true
            break
        end
    end
end

% output
centroids = circC(:, 1:2);
radii = circC(:, 3);

% visualization routine
if vizFlag == true
    % visualization loop
    Grains = cell(size(radii, 1),1);
    for i = 1:size(radii, 1)
        ang=0:0.01:2*pi;
        xp=radii(i)*cos(ang);
        yp=radii(i)*sin(ang);
        Grains{i} = [transpose(centroids(i,1)+xp) transpose(centroids(i,2)+yp)];
        Grains{i} = Grains{i}(1:5:end,:);
        plot(centroids(i,1)+xp, centroids(i,2)+yp);
        hold on
    end
    
    % colormap (prop = grain radii)
    elAsc = radii - min(radii);
    sc = 1/max(elAsc);
    elAsc = elAsc*sc;
    elAsc = elAsc*254;
    elAsc = elAsc+1;
    elAsc = round(elAsc);
    cm = jet(255);
    for i = 1:size(radii,1)
        Cols = cm(elAsc(i),:);
        fill(Grains{i}(:,1),Grains{i}(:,2),Cols); hold on;
    end
    axis equal;
end


