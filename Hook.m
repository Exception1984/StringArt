classdef Hook < handle
    properties
        width
        pos2d
        rotZAngleRadians
        
        p1
        p2
        p3
        p4
    end
    
    methods
        function obj = Hook(width, pos2d, rotZAngleRadians)
            obj.width = width;
            obj.pos2d = pos2d;
            obj.rotZAngleRadians = rotZAngleRadians;
            
            % Canonical position values
            % p1 ----- p4
            % |        |
            % |        |
            % |        |
            % p2-------p3
            
            hw = 0.5 * width;
            canP1 = [-hw hw]';
            canP2 = [-hw -hw]';
            canP3 = [hw -hw]';
            canP4 = [hw hw]';
            
            % Compute Corner Points
            rotMatrix = rotationMatrixZ(rotZAngleRadians);
            
            obj.p1 = pos2d + rotMatrix * canP1;
            obj.p2 = pos2d + rotMatrix * canP2;
            obj.p3 = pos2d + rotMatrix * canP3;
            obj.p4 = pos2d + rotMatrix * canP4;
        end
        
        function [AB, BA, DiagAB, DiagBA] = computeStrings(obj, hookB, stringWidth)
            aPoints = [obj.p1 obj.p2 obj.p3 obj.p4];
            bPoints = [hookB.p1 hookB.p2 hookB.p3 hookB.p4];
            points = [aPoints bPoints];
            
            bestLengthAB = [];
            bestLengthBA = [];
            bestLengthDiagAB = [];
            bestLengthDiagBA = [];
            
            r = 0.5 * stringWidth;
            rSquared = r * r;
            numRem = 3;
            
            thresh = 1.0e-8;
            
            for i = 1 : 4
                bPoint = bPoints(:, i);
                K = convhull([aPoints(1, :) bPoint(1)], [aPoints(2, :) bPoint(2)]);
                
                isA = K < 5;
                jumpIndices = find(isA ~= [isA(2:end)' isA(1)]');
                
                % Find direction of strings
                if K(jumpIndices(1)) < 5
                    % First jump is from A to B
                    
                    aToBIndex = K(jumpIndices(1));
                    bToAIndex = K(jumpIndices(2) + 1);
                    
                    ab = [points(:, aToBIndex) bPoint]';
                    ba = [bPoint points(:, bToAIndex)]';
                else
                    % First jump is from B to A
                    
                    aToBIndex = K(jumpIndices(2));
                    bToAIndex = K(jumpIndices(1) + 1);
                    
                    ab = [points(aToBIndex) bPoint]';
                    ba = [bPoint points(bToAIndex)]';
                end
                
                % Hesse Normal Form: n * x - d = 0
                
                remB = bPoints;
                remB(:, i) = [];
                
                % Analyze AB
                dirAB = ab(2, :) - ab(1, :);
                lengthAB = norm(dirAB);
                dirAB = dirAB / lengthAB;
                n = [dirAB(2) -dirAB(1)];
                d = n * bPoints(:, i);
                
                indicatorB = n * remB - d;
                indicatorB(abs(indicatorB) < thresh) = 0;
                
                if sum(indicatorB >= 0) == numRem || sum(indicatorB <= 0) == numRem
                    % ab is tangent to Hook B, test if it is AB or DiagAB
                    remA = aPoints;
                    remA(:, aToBIndex) = [];
                
                    indicatorA = n * remA - d;
                    indicatorA(abs(indicatorA) < thresh) = 0;
                    
                    if sum(indicatorA >= 0) == sum(indicatorB >= 0)
                        % AB
                        if isempty(bestLengthAB) || bestLengthAB > lengthAB
                            bestLengthAB = lengthAB;
                            AB = ab + 0.5 * stringWidth * repmat(n, 2, 1);
                        end
                    else
                        % DiagAB
                        if isempty(bestLengthDiagAB) || bestLengthDiagAB > lengthAB
                            bestLengthDiagAB = lengthAB;
                            d = norm(ab(1, :) - ab(2, :));
                            l = sqrt(0.25 * d^2 - rSquared);
                            intPoint = 0.5 * sum(ab);

                            alpha = asin(2 * r / d);
                            f = tan(alpha);

                            diagABDir = dirAB - f * n;
                            diagABDir = diagABDir / norm(diagABDir);

                            a = intPoint - l * diagABDir;
                            b = intPoint + l * diagABDir;

                            DiagAB = [a; b];
                        end
                    end
                end
                
                % Analyze BA
                dirBA = ba(2, :) - ba(1, :);
                lengthBA = norm(dirBA);
                dirBA = dirBA / lengthBA;
                n = [dirBA(2) -dirBA(1)];
                d = n * bPoints(:, i);
                
                indicatorB = n * remB - d;
                indicatorB(abs(indicatorB) < thresh) = 0;
                
                if sum(indicatorB >= 0) == numRem || sum(indicatorB <= 0) == numRem
                    % ab is tangent to Hook B, test if it is AB or DiagAB
                    remA = aPoints;
                    remA(:, bToAIndex) = [];
                
                    indicatorA = n * remA - d;
                    indicatorA(abs(indicatorA) < thresh) = 0;
                    
                    if sum(indicatorA >= 0) == sum(indicatorB >= 0)
                        % BA
                        if isempty(bestLengthBA) || bestLengthBA > lengthBA
                            bestLengthBA = lengthBA;
                            BA = ba + 0.5 * stringWidth * repmat(n, 2, 1);
                        end
                    else
                        % DiagAB
                        if isempty(bestLengthDiagBA) || bestLengthDiagBA > lengthBA
                            bestLengthDiagBA = lengthBA;
                            d = norm(ba(1, :) - ba(2, :));
                            l = sqrt(0.25 * d^2 - rSquared);
                            intPoint = 0.5 * sum(ba);

                            alpha = asin(2 * r / d);
                            f = tan(alpha);

                            diagBADir = dirBA + f * n;
                            diagBADir = diagBADir / norm(diagBADir);

                            a = intPoint + l * diagBADir;
                            b = intPoint - l * diagBADir;

                            DiagBA = [b; a];
                        end
                    end
                end
            end
        end
        
        function intersection = intersectsString(obj, p1, p2)
            thresh = 1.0e-8;
            stringDir = p2 - p1;
            stringDir = stringDir / norm(stringDir);
            
            n = [-stringDir(2) stringDir(1)];
            
            % x * n - d = 0
            d = n * p1;
            
            indicator = n * [obj.p1 obj.p2 obj.p3 obj.p4] - d;
            indicator(abs(indicator) < thresh) = 0;
                
            if sum(indicator >= 0) == 4 || sum(indicator <= 0) == 4
                intersection = false;
            else
                intersection = true;
            end
        end
        
        function angleRad = angleToXAxis(obj, vec)
            normVec = vec / norm(vec);
            xVec = [1 0]';
            
            angleRad = -acos(normVec' * xVec);
            
            % Evaluate quadrant and correct angle
            if normVec(2) < 0
                angleRad = -angleRad;
            end                    
        end
    end
end