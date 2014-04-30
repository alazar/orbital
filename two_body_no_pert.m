function state_vector = two_body_no_pert(m_planet, m_satellite, r_0, v_0, sim_time)

%TWO_BODY_NO_PERT solves the two body problem without including any
%perturbative inputs. The assumptions here is a satellite moving around a
%perfectly round planet. position and velocity components are w.r.t. an
%inertial frame (planet centered inertial [PCI])
%
% inputs:
%   m_planet - mass of planet (kg)
%   m_satellite - mass of satellite (kg)
%   r_0 - initial position of satellite [PCI] (km/s)
%   v_0 - initial velocity of satellite [PCI] (km/s)
%   sim_time - simulation time (s)
%
% outputs:
%   state_vector - column vector of [x y z vx vy vz]' at end of sim time
%   
    G = 6.6742e-20; 
    mu = G*(m_planet + m_satellite); % Gravitational parameter
    state_vector = rkf_integrate(@state_transition, [0 sim_time], [r_0 v_0]);
end

function dXdt = state_transition(X, t)
% our state transitions function: dX/dt = f(X, t)
% X is column vector of [x y z dx/dt dy/dt dz/dt]

    dXdt = zeros(size(X));
    dXdt(1) = X(4);
    dXdt(2) = X(5); % dX/dt is velocity
    dXdt(3) = X(6);
    
    r = sqrt(sum(X(1:3).^2));
    
    dXdt(4) = -mu*x/r^3; % acceleration components 
    dXdt(5) = -mu*y/r^3;
    dXdt(6) = -mu*z/r^3;

end