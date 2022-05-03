function distance=point_to_line_distance(pt, v1, v2)
%Calculate distance between a point and a line in 2D or 3D.
% syntax:
% distance = point_to_line(pt, v1, v2)
% pt is a nx3 matrix with xyz coordinates for n points
% v1 and v2 are vertices on the line (each 1x3)
% d is a nx1 vector with the orthogonal distances
%
% 2D inputs are extended to 3D by setting all z-values to 0, all inputs should be either nx2 or nx3
% (1x2 or 1x3 for v1 and v2).
%
% The actual calculation is a slightly edit version of this line (which only works for 3D points):
% distance=norm(cross(v1-v2,pt-v2))/norm(v1-v2)
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.3.2
% Date:    2020-06-28
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})
%check inputs
if nargin~=3
    error('HJW:point_to_line_distance:nargin',...
        'Incorrect number of inputs, expected 3.');
end
if ~isnumeric(pt) || ~any(size(pt,2)==[2 3]) || any(size(pt)==0)
    error('HJW:point_to_line_distance:pt_type_size',...
        'First input (pt) is not numeric or has an incorrect shape.')
end
if ~isnumeric(v1) || numel(v1)~=size(pt,2)
    error('HJW:point_to_line_distance:v_type_size',...
        ['Second input (v1) is not numeric or has an incorrect size.',char(10),...
        'Expected 1x3 or 1x2, which should match the first input.']) %#ok<CHARTEN>
end
if ~isnumeric(v2) || numel(v2)~=size(pt,2)
    error('HJW:point_to_line_distance:v_type_size',...
        ['Third input (v2) is not numeric or has an incorrect size.',char(10),...
        'Expected 1x3 or 1x2, which should match the first input.']) %#ok<CHARTEN>
end
%prepare inputs
v1=v1(:)';%force 1x3 or 1x2
v2=v2(:)';%force 1x3 or 1x2
if length(v1)==2,v1(3)=0;  end%extend 1x2 to 1x3 if needed
if length(v2)==2,v2(3)=0;  end%extend 1x2 to 1x3 if needed
if size(pt,2)==2,pt(1,3)=0;end%extend nx2 to nx3 if needed
v1_ = repmat(v1,size(pt,1),1);
v2_ = repmat(v2,size(pt,1),1);
%actual calculation
a = v1_ - v2_;
b = pt - v2_;
distance = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
%this is equivalent to the following line for a single point
%distance=norm(cross(v1-v2,pt-v2))/norm(v1-v2)

% Note, we'll want an alpha parameter for the location on the line too

end