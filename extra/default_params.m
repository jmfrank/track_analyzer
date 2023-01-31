function out_params = default_params(seg_type, t)



switch seg_type
    
    
    case {'cells','foci'}
        
        % Set gui fields. 
        out_params.bg=0;

        out_params.tile_size=25;
        out_params.ClipLimit=0.05;

        out_params.alpha=0;
        out_params.beta=0.1;
        out_params.gamma=8;
        out_params.delta=1;
        out_params.epsilon=0;
        out_params.sigmagradient=[5,5,2];

        out_params.median_filter = [5,5,1];

        out_params.simple_threshold = 100;
        out_params.percentile= 78*ones(1,t);
        out_params.I_sm_sigma=5;
        out_params.Hdepth=1000;
        % Now tack on a filtering step. 
        out_params.AbsMinVol=600;
        out_params.AbsMaxVol=3000;
        out_params.imclose_r=5;
        out_params.MeanFilterSensitivity=0.4;
        out_params.MeanFilterNeighborhood=[5,5,3];
        out_params.OutlierThreshold=60000;
        out_params.h_min_depth=1.5;
        out_params.merge_dist=20;
        out_params.WaterShedMaxVol = 3000;
        
        
        
    case 'spots'
        out_params.smooth_sigma=[3,3,1];
        out_params.spot_sigma = [3,3,1];
        out_params.local_thresh_percentile(1) = 90;
        out_params.local_thresh_percentile(2) = 90;
        
end

end