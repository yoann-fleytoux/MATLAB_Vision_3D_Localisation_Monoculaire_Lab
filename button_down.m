function button_down(src,evnt)

persistent indx;
global App3d_ini;
colork={'cyan' 'magenta' 'yellow' 'red' 'green' 'blue'};

if isempty(indx)
   indx = 1;
end

% The x y z points of the clicked line
xl = get(src, 'Xdata');
yl = get(src, 'Ydata');
zl = get(src, 'Zdata');

text((xl(1)+xl(2))/2+5, (yl(1)+yl(2))/2+5, (zl(1)+zl(2))/2, num2str(indx),'color',colork{indx});
set(src, 'color', colork{indx});
indx = indx + 1;

App3d_ini = [ App3d_ini; xl yl zl ];

end
