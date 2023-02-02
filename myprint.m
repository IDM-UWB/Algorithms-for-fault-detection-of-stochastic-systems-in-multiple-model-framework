function myprint(width,h,ratio)
%MYPRINT   Print a figure into pdf and eps files
%   MYPRINT(WIDTH,H,RATIO) Prints a figure with the handle H, into a files,
%   where WIDTH is width in cm and ratio is ratio of sides
%
%   Copyright (C) 2022  Ivo Puncochar
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

inch2cm = 2.54; % one inch is 2.54 cm
if nargin<3
    ratio = 6/8;
    if nargin<2
        h = gcf;
        if nargin<1
            width = 8*inch2cm;
        end
    end
end
width = width/inch2cm;
figure(h)
%set(gca,'Position', [0.13 0.11 0.775 0.815])
set(h,'PaperUnits', 'inches');
set(h,'PaperSize', width*[1 ratio]);
set(h,'PaperPosition', width*[0 0 1 ratio]);
filename = get(h,'Name');

% First save figure
saveas(h,filename,'fig')

% then print it to pdf
print(h,'-dpdf',filename)
%saveas(h,filename,'pdf')

% print it to eps
print(h,'-depsc2',filename)

% print it to wmf
print(h,'-dmeta',filename)

% print it to 
print(h,'-dpng',filename)