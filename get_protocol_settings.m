function model_settings = get_protocol_settings(nslices)
    % Create a struct with model and protocol settings

    
    model_settings = struct;

    model_settings.theta = 80; % flip angle (degrees)
    model_settings.TR = 2; % repitition time (s)
    model_settings.TE = 30/1000; % echo time (s)
    model_settings.T1 = 4.7; % 4.2734; %4.3; % (s) Normal Values of Magnetic Relaxation Parameters of Spine 
    model_settings.T2 = 1.6; % 1.5776; %2.00; % (s) Components with the Synthetic MRI Sequence
    
    model_settings.ds = 10; % Nr of dummy scans
    model_settings.n_slices = nslices; % Nr of excited slices of consideration
    
    model_settings.L = 3.4; % slice thickness (mm)
    model_settings.sp = 0.5; % spacing (mm)
    model_settings.stp = 20; % number of points per mm for mesh
    
    model_settings.xmin = -(model_settings.n_slices-1)*(model_settings.L+model_settings.sp) + model_settings.sp; % min position
    model_settings.xmax = model_settings.n_slices*(model_settings.L+model_settings.sp) - model_settings.sp; % max position
    
    
    model_settings.x0 = 0:(model_settings.L+model_settings.sp):(model_settings.n_slices-1)*(model_settings.L + model_settings.sp); % ground position slices (mm)
    model_settings.dt = 2/37; % Time per slice (s)
    model_settings.Nt = 170; %60; %170; % Nr of excitations per slice
    model_settings.t0 = 30; % time before excitations
    
    
    t_ex = zeros(model_settings.Nt + model_settings.ds, model_settings.n_slices); % time at excitation (s)
    for i = 1:model_settings.n_slices
        j = floor(i/2);
        if (mod(i, 2) == 0)
            t_ex(:,i) = model_settings.t0 + (18 + j)*model_settings.dt + model_settings.TR*(0:(model_settings.ds + model_settings.Nt-1));
        else
            t_ex(:,i) = model_settings.t0 + j*model_settings.dt + model_settings.TR*(0:(model_settings.ds + model_settings.Nt-1));
        end
    end

    model_settings.t_ex = t_ex;


end