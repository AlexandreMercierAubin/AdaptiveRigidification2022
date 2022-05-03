classdef contactInfo < handle
    %CONTACTINFO Contains information for each contact found
    %   The class can contain point position along with the indices of DOFs
    %   in the master list when all meshe dofs are concatenated.  The
    %   concatenated dof lists allow for contacts to be trivially warm
    %   started across time steps (i.e., contact between the same point and
    %   edge).  
    
    properties
        point       % position of the contact
        contactID = int64(0); % bit contatenation of nodes
        % list of nodes, up to 3, in the full list of nodes of all meshes,
        % each index must be less than 2^20 to fit in this contcatenation
        
        % lambda was stored here previously, but it is slow to access, so
        % the vector of all lambdas is instead kept in the quicksolvecache.
        
        normal
        tangent
        tangent2
        velocity    % if the constraint has time dependence, then this 
                    % should go on the LRHS of the schur complement contact
                    % solve (i.e., desired velocity in contact frame)
                    
        pointAlpha
        normalAlpha
        tangentAlpha
        pointAlphaInv
        normalAlphaInv
        tangentAlphaInv
        color = NaN
                   
        contactFinderId
        frictionCoefficient
        meshID
    end
    
    methods
        function obj = contactInfo( point, normal, tangent, mu, nodeIndices, contactFinderId ,tangent2, meshID)
            if nargin < 7
                tangent2 = [0,0,0];
            end
            if nargin < 8
                meshID = 1; 
            end
            
            %CONTACTINFO Construct an instance of this class
            obj.point = point;
            obj.normal = normal;
            obj.tangent = tangent;
            obj.tangent2 = tangent2;
            obj.frictionCoefficient = mu;
            obj.meshID = meshID;
            assert( numel(nodeIndices) <= 3 );
            
            tmp = zeros(1,3,'int64'); obj.contactFinderId = contactFinderId;
            tmp( 1:numel(nodeIndices) ) = int64(nodeIndices);
            
            obj.contactID = tmp(1) + bitshift(tmp(2),20) + bitshift(tmp(3),40);
            obj.contactFinderId = contactFinderId;
            obj.velocity = [0;0];
            % note lambda will be set after the solve, while the contact
            % info is created when the contact is found.
        end        
    end
end