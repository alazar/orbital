function state_matrix = RK_integrate(state_transition_fn, start_state, integration_time, tableau_id)

%RK_INTEGRATE does a numerical integration of the ODE dX/dt = f(X, t)
%where X is a vector, given the start state and interval of integrations:
%
% inputs:
%   state_transition: handle to equation for dX/dt = f(t, x)
%   start_state: vector of start state
%   integration_time: time interval over which to do the integration
%   tableau_id(opt): index specifying which Butcher Tableau to use for
%   integration.

    RK4_matrix = [0 0 0 0; ...
                  1/2 0 0 0; ...
                  0 1/2 0 0;...
                  0 0 1 0];
    RK4_weights = [0 1/2 1/2 0]; 
    RK4_nodes = [1/6 1/3 1/3 1/6];
    
    
    a = RK4_matrix;
    b = RK4_weights;
    c = RK4_nodes; 
    s = 4;
    
    b = repmat(b, length(start_state + 1), []);
    b(end, :) = 0;
    
    t = 0;
    h = 1;
    
    state_matrix = nan(length(start_state), ceil(integration_time/h));
    state_matrix(:, 1) = start_state;
    index = 1;
    
    while t < integration_time
        state_matrix(:, index + 1) = next_state_vec(state_transition_fn, a, b, c, h, s, t, state_matrix(:, index));
        t = t + h;
        index = index + 1;
    end
    
end

function state_vec = next_state_vec(state_transition_fn, a, b, c, h, s, t, current_state)

    k = zeros(length(current_state), s+1);
    
    k(:,1) = h*feval(state_transition_fn, current_state, t)';
    for k_index = 2:s
       t_ = t + c(k_index)*h;
       temp = a(k_index, 1:k_index-1);
       if length(temp) > 1
           weighed_ks = sum(transpose(repmat(temp, length(current_state), []).*k(:,1:k_index-1)))';
       else
           weighed_ks = temp.*k(:,1:k_index-1);
       end
       k(:, k_index) = h*feval(state_transition_fn, current_state + weighed_ks, t_);
    end
    
    state_vec = current_state + sum(transpose(b.*k(:,1:end-1)))';
        
end