%% cyclical signal detector. Looks forward to see if signal returns to current state with a max departure threshold. 


function [DATA,n_departures] = depart_detector( frames_present, sig, params )

%Length of signal
L = length(sig);
%Threshold of significant change. 
thresh = params.delS_threshold;

%Data structure for resets. 
DATA = struct('peak',[],'frames',[],'start',[],'stop',[],'net_delS',[]);

%Counter
c=0;

start_idx = 1;

while start_idx <= length(sig) - params.min_window_size
    
    %Initialize first step. 
    delS = sig(start_idx + 1 ) - sig(start_idx );
    dir = sign(delS);
    
    %Now see if things continue in same direction. 
    for stop_idx = start_idx + 2 : length(sig)
        
        %Next increment change in signal. 
        delS = sig(stop_idx) - sig(stop_idx-1);
        
        %Does the signal change direction? 
        condition_1 = -dir*delS > params.buff;
        
        %Does the signal stagnate? Can only test if stop_idx +1 >
        %params.stagnation_window
        if( (stop_idx-start_idx) >= params.stagnation_window)
            test_range = stop_idx - params.stagnation_window : stop_idx + params.stagnation_window;
            %Clip if it goes beyong signal. 
            test_range = intersect(test_range, 1:length(sig));
            test_sig   = sig( test_range )';
            %Fit slope through range. 
            P = polyfit(test_range,test_sig,1);
            %Is slope greater than threshold. 
            condition_2 = abs(P(1)) < params.stagnation_threshold;
        else
            condition_2 = 0;
        end
        
        
        if( condition_1 || condition_2 )
            %Check to see if magnitude / duration
            %are satisfactory. 
            net_delS = sig(stop_idx) - sig(start_idx );
            net_delT = stop_idx -1 - start_idx;
            if( abs(net_delS)>=thresh && net_delT >= params.min_window_size)
                c=c+1;
                %Figure out whether to stop at stop_idx or one before
                %depending on direction from stop-1 to stop? 
                stop_dir = sign(  sig(stop_idx) - sig(stop_idx-1) );
                if(stop_dir==dir)
                    stop_frame = stop_idx;
                else
                    stop_frame = stop_idx - 1;
                end
                DATA(c).frames = start_idx:stop_frame;
                DATA(c).net_delS = net_delS;
                
                %If there was a valid departure, then shift start to the
                %stop point.
                start_idx = stop_frame;
            else
                %There was no valid departure, so only increment start by 1
                %so we don't miss any dynamics. 
                start_idx = start_idx + 1;
            end

            break
            
        %Otherwise, the departure continues. 
        end
           
        
    end
    
    %If the departure went until end of signal, check if magnitude/duration
    %are ok, then break. 
    if( stop_idx == length(sig))
        
        net_delS = sig(stop_idx) - sig(start_idx );
        net_delT = stop_idx - start_idx;
        if( abs(net_delS)>=thresh && net_delT >= params.min_window_size)
            c=c+1;
            DATA(c).frames = start_idx:stop_idx;
            DATA(c).net_delS = net_delS;
        end
        
        break
        
    end
    
end

                
n_departures = c;            

        
     

end