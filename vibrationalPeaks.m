clc; clear; close all;

% Plot the band heads of the various nitrogen systems (i.e. Find the
%   position of the rotationless vibrational trasistion). Plot as a function
%   of the up vibrational level vs wavelength position. 
%
%       PROCEDURE:
%           1. Define the upper and lower vibrational levels of interest.
%           2. Create a MeshGrid matrix to evaluate all possible
%           transtion between vibrational levels of the two electronic
%           systems.
%           3. Pass the meshGrid variables to roVib and obtain all the
%           wavelengths.
%               - NOTE:     roVib can allow for NxN vibrational or
%               rotational number inputs. However, the other will need to
%               be an integer number. 
%           4. Reshape the NxN matricies into a Mx1 matrix for scatter
%           plotting.
%           5. Create annotation lables for each transition.
%           6. Create a Scatter plot with labels.
%           7. Repeat process for any other system.

v1  = 0:3;
v2  = 0:6;
J1  = 0;
J2  = 0;

figure
hold on;

% FIRST NEGATIVE SYSTEM
[V1,V2] = meshgrid(v1,v2);
wl      = roVib('1-',V1, V2, J1, J2, 0);

V1      = reshape(V1,[],1);
V2      = reshape(V2,[],1);
wl      = reshape(wl,[],1);
txt     = cell(size(V1));

for i   = 1:length(V1)
    txt(i)  = {sprintf('(%d, %d)',V1(i),V2(i))};
end

scatter(wl,V1,'filled','DisplayName', 'First Negative System');
text(wl,V1+0.1,txt,'FontWeight','bold','FontSize',8,...
        'HorizontalAlignment','center');

% SECOND POSITIVE SYSTEM
[V1,V2] = meshgrid(v1,v2);
wl      = roVib('2+',V1, V2, J1, J2, 0);

V1      = reshape(V1,[],1);
wl      = reshape(wl,[],1);

scatter(wl,V1+0.25,'filled','DisplayName', 'Second Positive System');
text(wl,V1+0.35,txt,'FontWeight','bold','FontSize',8,...
        'HorizontalAlignment','center');

% FIRST POSITIVE SYSTEM
[V1,V2] = meshgrid(v1,v2);
wl      = roVib('1+',V1, V2, J1, J2, 0);

V1      = reshape(V1,[],1);
V2      = reshape(V2,[],1);
wl      = reshape(wl,[],1);

scatter(wl,V1,'filled','DisplayName', 'First Positive System');
text(wl,V1+0.25,txt,'FontWeight','bold','FontSize',8,...
        'HorizontalAlignment','center');

axis([200 1100 -1 max(v1)+2])
legend('show');
