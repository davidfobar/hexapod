% createLink: Creates a structure to represent a robot's link based on its DH parameters
%
% link = createLink(a, d, alpha, theta, offset, centOfMass, mass, inertia, isRotary)
% This function creates a structure for a robot's link based on its Denavit-Hartenberg (DH)
% parameters and other physical properties of the link. This structure contains important
% information such as the joint angle (theta), offset, center of mass, mass, inertia, and
% whether the joint is rotary or prismatic.
%
% link = the created structure representing the link.
%
% a = DH parameter a (meters)
% d = DH parameter d (meters)
% alpha = DH parameter α (radians)
% theta = joint angle provided by encoder (radians)
% offset = offset between encoder orientation and DH zero-angle (radians or meters for prismatic)
% centOfMass = position of the link's center of mass
% mass = link mass (kg)
% inertia = link mass moment of inertia (kg m^2)
% isRotary = 1 if the joint is rotational, 0 if prismatic, -1 if static
%
% Raul
% 10871796
% MEGN 544
% 10/18/25

function link = createLink(a, d, alpha, theta, offset, centOfMass, mass, inertia, isRotary)
       % Create the link structure
    link.a = a;                     % DH parameter a (meters)
    link.d = d;                     % DH parameter d (meters)
    link.alpha = alpha;             % DH parameter alpha (radians)
    link.theta = theta;             % Joint angle provided by encoder (radians)
    link.offset = offset;           % Offset between encoder orientation and DH zero-angle
    link.centOfMass = centOfMass;   % Center of mass of the link
    link.mass = mass;               % Link mass 
    % (kg)
    link.inertia = inertia;         % Link mass moment of inertia (kg m^2)
    link.isRotary = isRotary;       % 1 if joint is rotational, 0 if prismatic, -1 if static

        % Handle empty theta or d based on joint type (rotary/prismatic/static)
    if isRotary == 1  % Rotary joint
        if isempty(theta)
            link.theta = 0;  % Default value for theta if empty
        end
    elseif isRotary == 0  % Prismatic joint
        if isempty(d)
            link.d = 0;  % Default value for d if empty
        end
    elseif isRotary == -1  % Static joint
        link.theta = theta;  % Static joints have no variable motion, changed from "NaN"
        link.d = d;         %changed from "NaN"
    end


end

