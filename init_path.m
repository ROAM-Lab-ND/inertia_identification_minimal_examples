function init_path()

    for dir = split(path(),':')
        if matches(dir,'spatial_v2')
            rmpath(dir);
        end
    end

    spatial_v2_path= [pwd '/spatial_v2_extended'];
    err_msg = 'You need to run `git submodule update --init --recursive` to load the submodule first';

    assert(exist(spatial_v2_path,'dir') > 0,err_msg);

    addpath(genpath(pwd));

    try
        cvx;
    catch
        assert(1==1,'CVX not detected. Please install it first.')
    end
end

