% Please note that this file is written specificallyf or cantilever beam
% problem

% bc_node_x = this vector stores the nodes where x dof is fixed
% bc_node_y = this vector stores the nodes where y dof is fixed


bc_node_x = [] ;
bc_node_y = [] ;

number_of_nodes_with_x_bc =  0 ;
number_of_nodes_with_y_bc =  0 ;

for i = 1:nno
    if abs(Xn(i,1)) <= eps % if x coordinate is zero
        bc_node_x = [bc_node_x ; i];
        number_of_nodes_with_x_bc = number_of_nodes_with_x_bc + 1 ;
    end
end

for i = 1:nno
    if abs(Xn(i,2)) <= eps  % if y coordinate is zero
        bc_node_y = [bc_node_y ; i];
        number_of_nodes_with_y_bc = number_of_nodes_with_y_bc + 1 ;
    end
end


% Form the boundary condition array. ir contains the nodes where the
% displacements are unknown. To get this array we first have to get the 
% dofs in the x and y directions. These are stored in arrays ir_x and ir_y.

% ir = zeros(number_of_nodes_with_x_bc+number_of_nodes_with_y_bc, 1) ;

% Form the boundary condition array in x direction

ir_x =  zeros(number_of_nodes_with_x_bc,1) ;

for i = 1: number_of_nodes_with_x_bc
    ir_x(i,1) = 2 * bc_node_x(i,1) - 1 ;
end

% Form the boundary condition array in y direction
ir_y =  zeros(number_of_nodes_with_y_bc,1) ;
for i = 1: number_of_nodes_with_y_bc
    ir_y(i,1) = 2 * bc_node_y(i,1) ;
end
 

 ir_disp = [ir_y ; ir_x]; 


% Using setdiff function of MATLAB to get the list of dofs where the
% displacement is unknown
ir = setdiff(1:ndof,ir_disp)' ;