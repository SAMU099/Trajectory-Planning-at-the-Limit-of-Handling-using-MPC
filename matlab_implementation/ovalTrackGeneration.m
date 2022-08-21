function Track = ovalTrackGeneration(TrackStraightLength, TrackCornerRadius, TrackWidth, TrackResApprox)
%OVALTRACKGENERATION(TrackStraightLength, TrackCornerRadius, TrackWidth, TrackResApprox) generates a
% parametrized track with two straights and two corners. The inputs are the length of the straights,
% the radius of the corner, the width of the track, and the approximate spacial resolution
% (approximate because the length of the track has to be perfectly
% divisible by it).

%% Track Resolution
TrackLength = 2*pi*TrackCornerRadius + 2*TrackStraightLength;
TrackNOPoints = round(TrackLength/TrackResApprox/2)*2;
TrackRes = TrackLength/TrackNOPoints;
CornerRes = TrackRes/TrackCornerRadius;
TrackSegmentsRemainder(1) = mod(TrackStraightLength/2,TrackRes);
TrackSegmentsRemainder(2) = mod(pi*TrackCornerRadius+TrackSegmentsRemainder(1),TrackRes);
TrackSegmentsRemainder(3) = mod(TrackStraightLength+TrackSegmentsRemainder(2),TrackRes);
TrackSegmentsRemainder(4) = mod(pi*TrackCornerRadius+TrackSegmentsRemainder(3),TrackRes);
TrackSegmentsRemainder(5) = mod(TrackStraightLength/2+TrackSegmentsRemainder(4),TrackRes);

%% Track Profile
% Calculate track coordinates
HalfTrackStraight1Coord = [zeros([floor(TrackStraightLength/2/TrackRes)+1,1]), (0:TrackRes:TrackStraightLength/2)'] + ...
    [TrackCornerRadius, 0];
HalfTrackCornerCoord = [cos(CornerRes-TrackSegmentsRemainder(1)/TrackCornerRadius:CornerRes:pi)', sin(CornerRes-TrackSegmentsRemainder(1)/TrackCornerRadius:CornerRes:pi)'] * TrackCornerRadius + ...
    [0, TrackStraightLength/2];
% Note: because of the calculation of the resolution (see TrackNOPoint)
% there is a point exactly at half track, so we can flip the straight
% already defined
HalfTrackStraight2Coord = flipud(HalfTrackStraight1Coord(2:end,:)) - [2*TrackCornerRadius, 0];
HalfTrackCoord = [HalfTrackStraight1Coord; HalfTrackCornerCoord; HalfTrackStraight2Coord];
Track.CentrelineCoord = [HalfTrackCoord; -1*HalfTrackCoord];
% Calculate track orientation/direction
TrackGradient = gradient(Track.CentrelineCoord')';
Track.Heading = atan2(TrackGradient(:,2),TrackGradient(:,1));
% Calculate track curvature
Track.Curvature = gradientAngle(Track.Heading')'./TrackRes;
% Calculate right edge
Track.Width = TrackWidth * ones(length(Track.Heading),1);
Track.LeftEdgeCoord = Track.CentrelineCoord + Track.Width .* [cos(Track.Heading+pi/2), sin(Track.Heading+pi/2)];
Track.RightEdgeCoord = Track.CentrelineCoord + Track.Width .* [cos(Track.Heading-pi/2), sin(Track.Heading-pi/2)];
% Add useful data to the structure
Track.Res = TrackRes;
Track.NOPoints = TrackNOPoints;

end

