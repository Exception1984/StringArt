function rotMatrixZ = rotationMatrixZ(angleRadians)
    rotMatrixZ = [cos(angleRadians) -sin(angleRadians);
                  sin(angleRadians)  cos(angleRadians)];
end