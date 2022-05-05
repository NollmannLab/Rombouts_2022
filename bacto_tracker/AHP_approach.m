
function [Rank] = AHP_approach(Cand,Cell)

%% Function that applies Analytical Hierarchy Processing method (from decision making theory)
%% to identify the optimal cell based on a given set of criteria
%% ------------------------------------------------------------------------------------------

% Input of this function is:
%   - Parameters of cell to be tracked
%   - Cell IDs of cell to be tracked and of candidates
%   - Parameters in matrix of all candidates

% AHP structure:
% GOAL: Find optimal candidate to be next cell in the track
% Criteria: Cell length - Area under mask - Overlap of masks in t=0 and t=1
% Candidates: Given as input to this function

%% ------------------------------------------------------------------------------------------

if isempty(Cand)
    Rank = [];
elseif size(Cand,1)==1
    Rank = horzcat(Cand(:,1),1);
else
    % STEP 1: Define the criteria and their respective priorities
    % -----------------------------------------------------------
    % We start with equal priorities
    % Sum of priorities for criteria needs to be equal to 1
    
    Prio_length = 1/3;
    Prio_area = 1/3;
    Prio_overlap = 1/3;
    
    
    
    % STEP 2: Set-up the matrix of scores
    % -----------------------------------
    % We define the scores by the following procedure:
    %   - For each parameter for a given criterion for each candidate we
    %     calculte 'how far the candidate is from the actual cell' by
    %     exp(-(abs((Cell-Cand)/Cell)*3))
    %   - Set-up of the pairwise matrix
    %   - Calculate the scores for in the matrix (by dividing the values of y)
    
    Total_y= [];
    Total_y_length = [];
    Total_y_area = [];
    
    for i = 1 : size(Cand,1)
        for j = 2 : size(Cand,2)-1
            y = exp(-(abs((Cell(:,j)-Cand(i,j))/Cell(:,j))*3));
            
            if j == 2
                Total_y_length = vertcat(Total_y_length,y);
            elseif j == 3
                Total_y_area = vertcat(Total_y_area,y);
            end
        end
    end
    
    % for cell length
    PWM_length = zeros(size(Cand,1));
    for n = 1 : size(PWM_length)
        for m = 1 : size(PWM_length)
            PWM_length(n,m) = Total_y_length(n)/Total_y_length(m);
        end
    end
    
    % for cell area
    PWM_area = zeros(size(Cand,1));
    for n = 1 : size(PWM_area)
        for m = 1 : size(PWM_area)
            PWM_area(n,m) = Total_y_area(n)/Total_y_area(m);
        end
    end
    
    % STEP 3: Calculate normalisation factor and local priorities
    % -----------------------------------------------------------
 
    Length = [];
    for col = 1 : size(PWM_length,2)
        
        Norm_length = sum(PWM_length(:,col));
        Local_length = [];
        
        for l = 1 : size(PWM_length,1);
            L = PWM_length(l,col)/Norm_length;
            Local_length = vertcat(Local_length,L);
        end
        
        Length = horzcat(Length,Local_length);
    end
    
    Avg_prio_length = [];
    for l = 1 : size(Length,1)
        Avg_length = mean(Length(l,:));
        Avg_prio_length = vertcat(Avg_prio_length, Avg_length);
    end
    
    Area = [];
    for col = 1 : size(PWM_area,2)
        
        Norm_area = sum(PWM_area(:,col));
        Local_area = [];
        
        for l = 1 : size(PWM_area,1);
            L = PWM_area(l,col)/Norm_area;
            Local_area = vertcat(Local_area,L);
        end
        
        Area = horzcat(Area,Local_area);
    end
    
    Avg_prio_area = [];
    for l = 1 : size(Area,1)
        Avg_area = mean(Area(l,:));
        Avg_prio_area = vertcat(Avg_prio_area, Avg_area);
    end
    
    Avg_prio_overlap = [];    
    Norm_overlap = sum(Cand(:,4));
    for l = 1 : size(Cand,1)
        L = Cand(l,4)/Norm_overlap;
        if isnan(L)
            L = 1;
        end
       
        Avg_prio_overlap = vertcat(Avg_prio_overlap, L);
    end
    
    
    % STEP 4: Calculate the global priorities and rank candidates
    % -----------------------------------------------------------
    Rank = [];
    R = [];
    for c = 1 : size(Cand,1)
        C_length = Prio_length*Avg_prio_length(c);
        C_area = Prio_area*Avg_prio_area(c);
        C_overlap = Prio_overlap*Avg_prio_overlap(c);
        C = C_length+C_area+C_overlap;
        R = vertcat(R,C);
    end
    
    Rank = horzcat(Cand(:,1),R);
    
end

