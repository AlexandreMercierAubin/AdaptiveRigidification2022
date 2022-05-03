function distance=point_to_line_segment_distance(pt, v1, v2)
%assume nx2,1x2,1x2 inputs
distance=zeros(size(pt,1),1);
%type 1 are all points closer to the end-point than to the rest of the line
%type 2 is closer to the start-point
%type 3 is the rest
%for type 3 we can use the point_to_line_distance FEX submission
%https://www.mathworks.com/matlabcentral/fileexchange/64396-point-to-line-distance
%dsp is the distance between the start and points
%dep is the distance between the end and points
%dse is the distance between start and end
%dpl is the distance between P and the projection on the entire line
dsp=sqrt( (pt(:,1)-v1(1)).^2 + (pt(:,2)-v1(2)).^2 );% + (pt(:,3)-v1(3)).^2 );
dep=sqrt( (pt(:,1)-v2(1)).^2 + (pt(:,2)-v2(2)).^2 );% + (pt(:,3)-v2(3)).^2 );
dse=sqrt( sum( (v1-v2).^2 ) );

dpl = point_to_line_distance(pt, v1, v2);

type1= sqrt(dse^2+dep.^2)<=dsp ;
type2= sqrt(dse^2+dsp.^2)<=dep ;
type3= ~type1 & ~type2;

distance(type1)=dep(type1);
distance(type2)=dsp(type2);
distance(type3)=dpl(type3);

end