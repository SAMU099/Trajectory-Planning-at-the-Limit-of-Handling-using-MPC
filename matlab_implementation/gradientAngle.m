function grad = gradientAngle(angles)
%GRADIENTANGLE(angles) calculates the gradient of an array of angles.
% This is different than the standard gradient of an array because of the
% periodicity of 2*pi. The closest angle difference in a period is chosen.

n = length(angles);
grad(1) = angdiff(angles(1),angles(2));
for i=2:n-1
    grad(i) = 0.5*angdiff(angles(i-1),angles(i+1));
end
grad(n) = angdiff(angles(n-1),angles(n));
end

