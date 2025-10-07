function [centroids, radii] = drop_roll(x_bounds, y_bounds, distMod, param1, param2, minDg, maxDg, vizFlag)
% test data
% x_bounds = [0 1]
% y_bounds = [0 3]
% distMod = 'Normal'
% param1 = 0.05
% param2 = 0.05
% minDg = 0.01
% maxDg = 0.1
% vizFlag = true
% Implementation of the classic 'drop and roll' algorithm for granular media
% generation

% Takes container bounds and particle size distribution as input and
% returns the centroids and radii of the generated particles

% Collision detection solutions are close form for speed
% Circle circle intersections solved with a modified (vectorized) version
% of the native Matlab circcirc.m function for speed

% INPUTS: 
% y_bounds - container y boundaries [ymin ymax]
% x_bounds - container x boundaries [xmin xmax]
% distMod - string defining distributional form. Supported distributions:
%         - 'Normal' (param1 = mean / param2 = standard deviation)
%         - 'Uniform' (param1 = minimum / param2 = maximum)
%         - 'Exponential' (param1 = mean / param2 = [])
%         - 'Lognormal' (param1 = mean / param2 = standard deviation)
%         - 'Weibull' (param1 = scale parameter / param2 = shape parameter)
%         - 'Mono' - (param1 = dG / param2 = dG)
% param1 - distribution parameter 1 (see above for definitions)
% param2 - distribution parameter 2 (see above for definitions)
% minDg - minimuum grainsize (scalar)
% maxDg - maximum grainsize (scalar)

% OUTPUTS: 
% centroids - generated circle centroids (nx2)
% radii - generated circle radii (nx1)

% Author: Thomas Seers: thomas.seers@qatar.tamu.edu (2025)

% test distribution input
testString = [{'Exponential'}; {'Normal'}; {'Uniform'};  {'Lognormal'}; {'Weibull'}; {'Mono'}];
logiSt = false(6,1);
for i = 1:6
    logiSt(i) = strcmp(testString{i}, distMod);
end

if ~isequal(sum(logiSt), 1)
    error('Input distribution not recognized');
end

% expo flag (one input)
if logiSt(1) == true
    distroFlag = 1;
else if logiSt(6) == true
        distroFlag = -1;
    else
        distroFlag = 0;
    end
end

% initialize
% expo flag (one input)
if distroFlag == 1
    circIn = random(distMod, param1);
else if distroFlag == -1
        circIn = param1;
    else
        circIn = random(distMod, param1, param2);
    end
end

if ~isequal(distroFlag, -1)
    while circIn < minDg || circIn > maxDg
        if distroFlag == 1
            circIn = random(distMod, param1);
        else
            circIn = random(distMod, param1, param2);
        end
    end
end
circR = circIn/2;

% x cent
x_cent = rand*x_bounds(2);
while x_cent < x_bounds(1) + circR || x_cent > x_bounds(2) - circR
    x_cent = rand*x_bounds(2);
end

circC = [x_cent y_bounds(1) + circR circR];
circB = [x_cent - circR x_cent + circR]; % bounds in x


% set termination criteria max(circC(:,2)) > y_bounds(2)
while max(circC(:,2)) < y_bounds(2)
    
    % initialize
    % expo flag (one input)
    if distroFlag == 1
        circIn = random(distMod, param1);
    else if distroFlag == -1
            circIn = param1;
        else
            circIn = random(distMod, param1, param2);
        end
    end
    
    if ~isequal(distroFlag, -1)
        while circIn < minDg || circIn > maxDg
            if distroFlag == 1
                circIn = random(distMod, param1);
            else
                circIn = random(distMod, param1, param2);
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
        
        % append to list
        circC = vertcat(circC, circC_T);
        circB = vertcat(circB, circB_T);
        stopFlag = true; % stop criteria met
        
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
        
        % d = (p1y - p2y) + sqrt((r1 + r2)^2 - (p1x - p2x)^2)
         
       
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
            if isnan(iCirc1) && isnan(iCirc2)
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
                % append to list
                circC = vertcat(circC, circC_T);
                circB = vertcat(circB, circB_T);
                stopFlag = true;
            end
            
        else if isnan(iCirc1) && isnan(iCirc2) % no circles in the arc-path
                circC_T = [contain(iCont,:) circR]; % 'roll'
                circB_T = [circC_T(1) - circR circC_T(1) + circR]; % bounds in x
                
                % stop condition met: append to list
                circC = vertcat(circC, circC_T);
                circB = vertcat(circB, circB_T);
                stopFlag = true;
                
            else % both container and circles lie in the arc-path
                % determine whether circ or container gives highest
                % centroid
                [mAll, iAll] = max([mCont, mCirc1, mCirc2]);
                if iAll == 1 % container is first collision
                    circC_T = [contain(iCont,:) circR]; % 'roll'
                    circB_T = [circC_T(1) - circR circC_T(1) + circR]; % bounds in x
                    
                    % stop condition met: append to list
                    circC = vertcat(circC, circC_T);
                    circB = vertcat(circB, circB_T);
                    stopFlag = true;
                    
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
                        % append to list
                        circC = vertcat(circC, circC_T);
                        circB = vertcat(circB, circB_T);
                        stopFlag = true;
                    end
                end
            end
        end
    end

    % test if stop criteria is met: break or conditional loop
    if stopFlag == true
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
                    % append to list
                    circC = vertcat(circC, circC_T);
                    circB = vertcat(circB, circB_T);
                    stopFlag = true;
                end
                
            else if isnan(iCirc1) && isnan(iCirc2) % no circles in the arc-path
                    circC_T = [contain(iCont,:) circR]; % 'roll'
                    circB_T = [circC_T(1) - circR circC_T(1) + circR]; % bounds in x
                    
                    % stop condition met: append to list
                    circC = vertcat(circC, circC_T);
                    circB = vertcat(circB, circB_T);
                    stopFlag = true;
                    
                else % both container and circles lie in the arc-path
                    % determine whether circ or container gives highest
                    % centroid
                    [mAll, iAll] = max([mCont, mCirc1, mCirc2]);
                    if iAll == 1 % container is first collision
                        circC_T = [contain(iCont,:) circR]; % 'roll'
                        circB_T = [circC_T(1) - circR circC_T(1) + circR]; % bounds in x
                        
                        % stop condition met: append to list
                        circC = vertcat(circC, circC_T);
                        circB = vertcat(circB, circB_T);
                        stopFlag = true;
                        
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
                            % append to list
                            circC = vertcat(circC, circC_T);
                            circB = vertcat(circB, circB_T);
                            stopFlag = true;
                        end
                    end
                end
            end
            
            % test if stop criteria is met: break or conditional loop
            if stopFlag == true
                continue
            else
                circHit = circC(iCirc,:); % new rotational centroid and radius
            end
            
        end
    end
end

% output
centroids = circC(:, 1:2);
radii = circC(:, 3);

% visualization routine
if vizFlag == true
    if distroFlag == -1
        for i = 1:size(radii, 1)
            ang=0:0.01:2*pi;
            xp=radii(i)*cos(ang);
            yp=radii(i)*sin(ang);
            Grains{i} = [transpose(centroids(i,1)+xp) transpose(centroids(i,2)+yp)];
            Grains{i} = Grains{i}(1:5:end,:);
            plot(centroids(i,1)+xp, centroids(i,2)+yp, 'k');
            hold on
        end
        
        Cols = [0.8 0.7 0.2];
        for i = 1:size(radii,1)
            fill(Grains{i}(:,1),Grains{i}(:,2),Cols); hold on;
        end
        axis equal;
        
    else
        % visualization loop
        Grains = cell(size(radii, 1),1);
        for i = 1:size(radii, 1)
            ang=0:0.01:2*pi;
            xp=radii(i)*cos(ang);
            yp=radii(i)*sin(ang);
            Grains{i} = [transpose(centroids(i,1)+xp) transpose(centroids(i,2)+yp)];
            Grains{i} = Grains{i}(1:5:end,:);
            plot(centroids(i,1)+xp, centroids(i,2)+yp, 'k');
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
            if strcmp(distMod, 'Mono') == true
                Cols = 'r';
            else
                Cols = cm(elAsc(i),:);
            end
            fill(Grains{i}(:,1),Grains{i}(:,2),Cols); hold on;
        end
        axis equal;
    end
end







