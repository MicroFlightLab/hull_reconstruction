function beta=signedAngleTwo3DVectors(startVec,endVec,rotateVec,oneVecFlag)
% Coded by Seth B. Wagenman, Riverside Research Institute, July 2020
% Calculate the angle (in radians) between two three-D vectors defined as
% the angle from startVec to endVec in the plane defined by its normal
% vector being either the cross product or negation of startVec x endVec.
%
% Pass 1 for 'oneVecFlag' if the first three inputs are single vectors.
%
% Usage details for the 'rotateVec' argument are as follows:
% rotateVec is a user-specified (normal) vector pointing from the side of
% the plane that startVec and endVec lie in, which points in the direction
% about which the angle is measured.  If you do not know what this normal
% vector is, then pass an empty vector and the function will substitute
% their cross product, but in this case the angles will always be positive,
% so the user must employ some other test external to the function to
% determine when the angle should be positive or negative based on the
% directions of input vectors startVec and endVec, and adjust accordingly.


	function unitLengthVec3D=unit3DNonZeroVectorByRow(vectorMatrix,oneVecFlag)
	% Normalize three-dimensional vector(s) to have Euclidean length of one.
	% Fails to process vector(s) of zero length; specify if only passing one.
	% Coded by Seth B. Wagenman, Riverside Research Institute, August 2020

	% Transform inputs if necessary and/or ask for user clarification of inputs
	[firstDim, secondDim] = size(vectorMatrix);
	rotationNeededPriorToReturn = 0;
	if ~oneVecFlag && (firstDim == secondDim) && (firstDim == 3)
		vectorDimensionUnclearMessage = ...
			'Is the vector dimension along the first axis?  Enter y or n:  ';
		userSpecifiedVectorDimension = 0;
		while ~userSpecifiedVectorDimension
			userInputVectorDimension = input(vectorDimensionUnclearMessage);
			if userInputVectorDimension == 'y'
				vectorMatrix = vectorMatrix.';
				lengthDim = secondDim;
				vecDim = firstDim;
				rotationNeededPriorToReturn = 1;
				userSpecifiedVectorDimension = 1;
			elseif userInputVectorDimension == 'n'
				lengthDim = firstDim;
				vecDim = secondDim;
				userSpecifiedVectorDimension = 1;
			else
				disp('Try entering eitehr lowercase y OR lowercase n ONLY...')
				disp('HINT: CAPS lock might be active on your keyboard...')
			end
		end
	elseif oneVecFlag && (firstDim == 3)
		vectorMatrix = vectorMatrix.';
		lengthDim = secondDim;
		vecDim = firstDim;
		rotationNeededPriorToReturn = 1;
	elseif (firstDim ~= 3) && (secondDim ~= 3)
		error('At least one dimension of input vector(s) must be three...')
	else
		lengthDim = firstDim;
		vecDim = secondDim;
	end

	% unitLengthVec3D is what the function returns
	EuclideanNorm3D = sqrt(sum(vectorMatrix.^2, 2));
	if any(EuclideanNorm3D <= eps) % Ensure vector does not have zero length
		error('Found a vector of length zero; cannot continue processing.')
	else
		unitLengthVec3D = vectorMatrix./ EuclideanNorm3D;
	end

	% Return a unit vector with the same orientation as the input vectorMatrix
	if rotationNeededPriorToReturn
		unitLengthVec3D = unitLengthVec3D.'
	end

	end % End of function


% Check for different sized inputs and other usage errors
if size(startVec) ~= size(endVec)
    error('Input vectors must have same dimensions, in same order...')
elseif ~isempty(rotateVec)
	if ~oneVecFlag
		if size(startVec) ~= size(rotateVec)
			error('Input vectors must have same dimensions, in same order...')
		end
		if any(abs(dot(rotateVec, cross(startVec, endVec), 2)) < 2.5 * eps)
			error('rotateVec is nearly coplanar with startVec and endVec.')
		end
	else
		if abs(dot(rotateVec, cross(startVec, endVec))) < 2.5 * eps
			error('rotateVec is nearly coplanar with startVec and endVec.')
		end
	end
end

firstDim = size(startVec, 1);
secondDim = size(startVec, 2);

% Check if user passed only one set of vectors, and/or rotated them
singleInputsRotated = oneVecFlag && (firstDim == 1) && (secondDim == 3);
multipleInputsRotated = ~oneVecFlag && (firstDim == 3) && (secondDim > 1);
if (singleInputsRotated || multipleInputsRotated)
    startVec = startVec.';
    endVec = endVec.';
end
longDimension = size(startVec, 1);

% If the user passes an empty rotateVec, use startVec x endVec
temp = cross(startVec, endVec); % temp is normal to both startVec & endVec
if isempty(rotateVec) % take unit normal to the plane in default direction
    rotateVec = unit3DNonZeroVectorByRow(temp, oneVecFlag);
else % Use cross product (or negation if it points in specified direction)
	if ~oneVecFlag % Flip the cross products/normal vectors if necessary
		negationVector = zeros(longDimension, 1);
		for row = 1:longDimension
			negationVector(row) = ...
				sign(dot(rotateVec(row, :), temp(row, :)));
			if negationVector(row) < 0
				temp(row, :) = -temp(row, :);
			end
		end
	else % Flip the cross product/normal vector if necessary
		if sign(dot(rotateVec, temp)) < 0
			temp = -temp
		end
	end
	rotateVec = unit3DNonZeroVectorByRow(temp, oneVecFlag);
end

% For derivation of the following formulas, see StackOverflow answer at:
% https://stackoverflow.com/a/33920320/12763497 or search for the term
% "signed angle between two 3D vectors with the same origin within the same
% plane" in the search bar and look for the 2015 answer by Adrian Leonhard.
% That proof is the basis of the corresponding post in math.stackexchange:
% https://math.stackexchange.com/a/3759166/747468, which explains how to
% calculate the four-quadrant arctangent based on Adrian Leonhard's proof.

% Numerator of tan(beta) = sin(beta)
if ~oneVecFlag
    sineBeta = dot(rotateVec, cross(startVec, endVec), 2);
else
    sineBeta = dot(rotateVec, cross(startVec, endVec));
end

% Denominator of tan(beta) = cos(beta)
if ~oneVecFlag
    cosineBeta = dot(startVec, endVec, 2);
else
    cosineBeta = dot(startVec, endVec);
end

% Two-argument (four-quadrant) arctangent takes sine & cosine as arguments
beta = atan2(sineBeta, cosineBeta);

acos(dot(startVec', endVec'))
2*pi - acos(dot(startVec', endVec'))

v_rot = Functions.rodrigues_rot(startVec,endVec,acos(dot(startVec', endVec'))')


end % End of function
