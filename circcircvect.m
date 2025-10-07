function [xout, yout] = circcircvect(x1,y1,r1,x2,y2,r2)
%CIRCCIRC  Intersections of circles in Cartesian plane (vector input)
%
%  [xout,yout] = CIRCCIRC(x1,y1,r1,x2,y2,r2) finds the points
%  of intersection (if any), given two lists of circles, each defined by
%  centers and radii in x-y coordinates.  In general, two points are
%  returned.  When the circles do not intersect or are identical,
%  NaNs are returned.  When the two circles are tangent, two identical
%  points are returned.  All inputs must be scalars.
%
%  See also LINECIRC.

% Copyright 1996-2007 The MathWorks, Inc.
% Written by:  E. Brown, E. Byrns / modified by T. D. Seers 2020

r3 = sqrt((x2-x1).^2 + (y2-y1).^2);
anought=atan2((y2-y1),(x2-x1));
 
%Law of cosines
aone=acos(-((r2.^2-r1.^2-r3.^2)./(2*r1.*r3)));

alpha1 = anought+aone;
alpha2 = anought-aone;

xout=[x1 x1]+[r1 r1].*cos([alpha1 alpha2]);
yout=[y1 y1]+[r1 r1].*sin([alpha1 alpha2]);

indx1=find(r3>r1+r2);  % too far apart to intersect
indx2=find(r2>r3+r1);  % circle one completely inside circle two
indx3=find(r1>r3+r2);  % circle two completely inside circle one
indx4=find((r3<10*eps)&(abs(r1-r2)<10*eps)); % circles identical
indx=[indx1(:);indx2(:);indx3(:);indx4(:)];

% Replace complex results (no intersection or identical)
% with NaNs.
if ~isempty(indx)
    xout(indx,:) = NaN;    yout(indx,:) = NaN;
end
