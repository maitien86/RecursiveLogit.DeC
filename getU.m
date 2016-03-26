%   Get Utility
%%
function Ufull = getU(x, isLS)

    global incidenceFull;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;  
    global GLS;
    global isFixedUturn;
    global memShip;    
    mu = 1;
    u1 = x(1) * EstimatedTime;
    u2 = x(2) * TurnAngles;
    u3 = x(3) * LeftTurn;
    if isFixedUturn == false
        u4 = x(4) * Uturn;
        if isLS == true
            u5 = x(5) * GLS;
        else
            u5 = 0;
        end
    else 
        u4 = -20 * Uturn;
        if isLS == true
            u5 = x(4) * GLS;
        else
            u5 = 0;
        end
    end
    u = (u1 + u2 + u3 + u4 + u5);
    Ufull = incidenceFull .* u;
end


