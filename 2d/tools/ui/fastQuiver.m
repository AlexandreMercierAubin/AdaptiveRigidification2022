function h = fastQuiver( position, size, varargin )
    
    if ( numel(varargin) > 0 )
        color = varargin{1};
    else
        color = [1,0,0];
    end
    if( numel(varargin) > 1 )
        scale = varargin{2};
    else
        scale = 0.1;
    end

    dest = position + size;
    ortho = [-size(2), size(1)];
    dir1 = scale * (ortho - size);
    dir2 = scale * (-ortho - size);

    x1 = [ position(1), dest(1) + dir1(1), dest(1) + dir2(1) ];      
    x2 = [ dest(1), dest(1), dest(1) ];
    y1 = [ position(2), dest(2) + dir1(2), dest(2) + dir2(2) ];      
    y2 = [ dest(2), dest(2), dest(2) ];
    
    h = plot( [x1;x2], [y1;y2], 'color', color, 'LineWidth', 2 );
    
%     plot([position(1); dest(1)], [position(2); dest(2)], 'color', color);
%     plot([dest(1) + dir1(1) * 0.1; dest(1)], [dest(2) + dir1(2) * 0.1; dest(2)], 'color', color);
%     plot([dest(1) + dir2(1) * 0.1; dest(1)], [dest(2) + dir2(2) * 0.1; dest(2)], 'color', color);
end

