function [value,isterminal,direction] = zillEvent( t , x , rc )
% detects a contact event in a Zilli system simulation

value = sqrt(x(1)*x(1)+x(2)*x(2))-rc;
 direction=0; % for a finite impact simulation this should be 0 to detect 
                %transitions both in and out of contact
                % for a rigid impact system this should be 1 , to only
                % detect the intial impactt (otherwise odd nonphysical
                % effects occur)
isterminal=true ;
end

