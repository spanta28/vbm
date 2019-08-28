# web_vbm
Runs vbm pipeline using singularity image 

Sample_input_data is provided in this repo. 

Optional parameters info:
The optional params info is detailed in vbm optional params description.xlsx

json input file is required. 

Ex json input to run_script.py: 
{"data":[["/home/vagrant/input1","sub1","sess1"],["/home/vagrant/input2","sub2","sess2"]], "options":{"smoothing_x_mm":10.0,"smoothing_y_mm":10.0,"smoothing_z_mm":10.0,"smoothing_implicit_masking":false,"reorient_params_x_mm":0, "reorient_params_y_mm":0, "reorient_params_z_mm":0,"reorient_params_pitch":0, "reorient_params_roll":0, "reorient_params_yaw":0,"BIAS_REGULARISATION":0.0001,"FWHM_GAUSSIAN_SMOOTH_BIAS":60,"affine_regularization":"mni","warping_regularization":[0, 1e-3, 0.5, 0.05, 0.2], "sampling_distance":3,"mrf_weighting":1.0,"cleanup":1}}
Note: 
trendscenter/vbm_web_image hosted on dockerhub is private. So contact me with your docker ID so that I can give you access to it. 

You can run the code in foll. ways:
1) Build dockerfile file from here 
docker build -t vbm_image_name Dir_containing_Dockerfile


Sometimes you may create a docker container on your local machine and when you move it elsewhere you may not have permissions to run somethings inside that container. So go inside the container and run "chmod -R 755 /" or equivalent to make sure your required software is executable.

To change the cache dir location , incase you dont have permissions to write to the default
in bash:
export SINGULARITY_CACHEDIR="/user/container"

same with writing to default /tmp folder
export SINGULARITY_TMPDIR="/user/container/tmp"

then run singularity build vbm.img trendscenter/vbm_web_image
if you want to create the vbm.img as writable image for editing , then :
singularity build --sandbox vbm.img docker://trendscenter/vbm_web_image
or if you have permissions to use --writable then:
singularity build --writable vbm.img docker://trendscenter/vbm_web_image

This is what I used and it works fine.

Sometimes due to changes in open source code of packages being used in this container, things may not work as expected. 
If you want to use the image that I already built then proceed to step 2) 

2) use trendscenter/vbm_web_image located on dockerhub

To change the cache dir location , incase you dont have permissions to write to the default 
in bash:
export SINGULARITY_CACHEDIR="/export/user/container"

same with writing to default /tmp folder
export SINGULARITY_TMPDIR="/user/container/tmp"

singularity build vbm.img docker://trendscenter/vbm_web_image

if you want to create the vbm.img as writable image for editing , then :
singularity build --sandbox vbm.img docker://trendscenter/vbm_web_image
or if you have permissions to use --writable then:
singularity build --writable vbm.img docker://trendscenter/vbm_web_image

This is what I used and it works fine.

3) Once downloaded to local machine:

3.a) You can use run_script.py and supply required args to have the script call singularity image and run the pipeline or do it manually as in step 3.b)

 python run_script.py -h
usage: run_script.py [-h] -json JSON -o OUTPUT -tmp TMP -image_path IMAGE_PATH
                     [-g SLURM] [-type {nifti,dicoms}]

Run vbm pipeline on mprage T1w dicoms

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -json JSON, --json JSON
                        Enter the full path of location of json input file . 
  -o OUTPUT, --output OUTPUT
                        Existing output directory to store vbm outputs
  -tmp TMP              Nipype tmp directory, it can be same as output
                        directory
  -image_path IMAGE_PATH
                        Enter full path to singularity vbm image name

optional arguments:
  -g SLURM, --slurm SLURM
                        Use slurm, True or False, default=False
  -type {nifti,dicoms}, --data_type {nifti,dicoms}
                        input data type. Ex:'nifti','dicoms'


python run_script.py -image_path /container/gsu_vbm.img/ -i /container/input -o /container/output -json /container/json_data/input.json -tmp /output -sub sub1 -sess sess1

vbm pipeline is running...
vbm output directory: /out/sub1/sess1/anat
VBM preprocessing completed.

3.b) Singularity command
singularity exec --contain --workdir NIPYPE_TMP_DIR -B INPUT_DIR:/data,OUTPUT_DIR:/out vbm.img python3 /computation/run_vbm.py -sub SUBJECT_ID -sess SESSION_ID -json '/json_data/input.json'

    -c/--contain        This option disables the sharing of filesystems on
                        your host (e.g. /dev, $HOME and /tmp).
                        
    --workdir stores tmp files that singulairty image will write , in this case nipype(package that implements vbm pipeline)         writes to this. Nipype will try to write to host /tmp or var_tmp and will have permission issue. Hence usage of this           option.
    
    -B Binds path on host to container
    
    Variables:
    NIPYPE_TMP_DIR : Directory to store nipype tmp files, can be anywhere that the image has write access to. Could be stored     in sessions output dir in the logs as well.
    INPUT_DIR : Full path to directory where input data is located : Data can be a nifti file or dicoms
    The script assumes there is only one nifti file
    OUTPUT_DIR : Directory where outputs will be written to
    
    Required arguments:
    Each input scan is from sepecific subject and sessionthat The foll. arguemnts are used to create output directories using     these args
    SUBJECT_ID : subject id 
    SESSION_ID : session id
    json : full path to accessible json file for the singularity image
    
The optional params info is detailed in vbm optional params description.xlsx
    
    
    
Sample run:
Input dir
ls /input/
Re_mprage_5e_RMS_0003_0001.nii

Command executed
singularity exec --contain --workdir /test/tmp -B /input:/data,/output:/out gsu_vbm.img python3 /computation/run_vbm.py -sub sub-02 -sess -json '/json_data/input.json'

vbm pipeline is running...
/out/sub-02/session-01/anat
VBM preprocessing completed.

Output data
ls /output/
outputs_description.txt  quality_control_readme.txt  sub-02

ls /output/sub-02/session-01/anat/
Re_mprage_5e_RMS_0003_0001.nii	vbm_spm12

ls /output/sub-02/session-01/anat/vbm_spm12/
c1Re.nii  c4Re.nii  mwc1Re.nii	mwc4Re.nii  Re.nii	 smwc2Re.nii  smwc5Re.nii  swc2Re.nii  swc5Re.nii	   wc1Re.nii  wc3Re.nii  wc6Re.nii
c2Re.nii  c5Re.nii  mwc2Re.nii	mwc5Re.nii  Re_seg8.mat  smwc3Re.nii  smwc6Re.nii  swc3Re.nii  swc6Re.nii	   wc1Re.png  wc4Re.nii
c3Re.nii  c6Re.nii  mwc3Re.nii	mwc6Re.nii  smwc1Re.nii  smwc4Re.nii  swc1Re.nii   swc4Re.nii  vbm_corr_value.txt  wc2Re.nii  wc5Re.nii
