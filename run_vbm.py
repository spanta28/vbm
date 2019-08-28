#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer includes the interface adapter(IA) for parsing input args to read in pre-processing structural nifti T1w scans
This layer sends the output to vbm_use_cases_layer with the appropriate inputs to run the pipeine using nipype interface

Sample run for input data of nifti 
python3 /computation/run_vbm.py -sub sub1 -sess sess1 -smooth 6
"""
import contextlib,traceback


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr
    e.g.:
    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')
    """

    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, 'w')
        os.dup2(dest_file.fileno(), stdchannel.fileno())

        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()

import json, argparse, getopt, re
import warnings, os, glob, sys
import nibabel as nib
from  math import radians as rad


## Load Nipype spm interface ##
from nipype.interfaces import spm

import vbm_use_cases_layer

# Create a dict to store all paths to softwares,templates & store parameters, names of output files

template_dict = {
    'spm_version':
    '12.7507',
    'matlab_cmd':
    '/opt/spm12/run_spm12.sh /opt/mcr/v95 script',
    'spm_path':
    '/opt/spm12/fsroot',
    'tpm_path':
    '/opt/spm12/fsroot/spm/spm12/tpm/TPM.nii',
    'transf_mat_path':
    os.path.join('/computation', 'transform.mat'),
    'dicom_dir':
    '',
    'dicom_output_dir':
    '',
    'scan_type':
    'T1w',
    'subject':
    '',
    'session':
    '',
    'reorient_params_x_mm': 0,
    'reorient_params_y_mm': 0,
    'reorient_params_z_mm': 0,
    'reorient_params_pitch': 0,
    'reorient_params_roll': 0,
    'reorient_params_yaw': 0,
    'reorient_params_x_scaling': 1,
    'reorient_params_y_scaling': 1,
    'reorient_params_z_scaling': 1,
    'reorient_params_x_affine': 0,
    'reorient_params_y_affine': 0,
    'reorient_params_z_affine': 0,
    'bounding_box':
    '',
    'FWHM_SMOOTH': [10, 10, 10],
    'BIAS_REGULARISATION':
    0.0001,
    'FWHM_GAUSSIAN_SMOOTH_BIAS':
    60,
    'affine_regularization': 'mni',
    'warping_regularization': [0, 1e-3, 0.5, 0.05, 0.2],
    'sampling_distance': 3.0,
    'mrf_weighting':
    1.0,
    'cleanup':
    1,
    'implicit_masking':
    False,
    'correlation_value':
    0.90,
    'vbm_output_dirname':
    'vbm_spm12',
    'output_zip_dir':
    'vbm_outputs',
    'qa_flagged_filename':
    'QA_flagged_subjects.txt',
    'display_image_name':
    'wc1Re.png',
    'display_pngimage_name':
    'Gray matter (Normalized)',
    'cut_coords': (0, 0, 0),
    'display_nifti':
    'wc1Re.nii',
    'qc_nifti':
    'swc1*nii',
    'vbm_qc_filename':
    'vbm_corr_value.txt',
    'outputs_manual_name':
    'outputs_description.txt',
    'flag_warning':
    ' QC warning: Only half of the data is pre-processed or passed the QA, please check the data!',
    'bids_outputs_manual_content':
    "Prefixes descriptions for segmented images:c1-Gray matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
    "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
    "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'nifti_outputs_manual_content':
    "Prefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
    "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
    "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'dicoms_outputs_manual_content':
    "Prefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
    "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
    "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'qc_readme_name':
    'quality_control_readme.txt',
    'qc_readme_content':
    "In each subject's anat/vbm_spm12 directory,vbm_corr_value.txt gives the correlation value of the swc1*nii file with spm12 adult brain tissue probability maps (TPM.nii) file from SPM12 toolbox"
    "\nIf your subjects are kids/adolescent, take this into considersation as correlation value may not pass threshold. For an adult study, its suggested that scans with correlation value <=0.91 "
    "\nshould be manually looked into for possible reorientation to the ac-pc line. Subjects that do not pass this QA metric will be saved in vbm_outputs/QA_flagged_subjects.txt"
}
"""
More info. on keys in template_dict

spm_path is path to spm software inside docker . 
SPM is Statistical Parametric Mapping toolbox for matlab 
Info. from http://www.fil.ion.ucl.ac.uk/spm/
"Statistical Parametric Mapping refers to the construction and assessment of spatially extended statistical processes used to test hypotheses about functional imaging data. 
These ideas have been instantiated in software that is called SPM.
The SPM software package has been designed for the analysis of brain imaging data sequences. 
The sequences can be a series of images from different cohorts, or time-series from the same subject. 
The current release is designed for the analysis of fMRI, PET, SPECT, EEG and MEG."

tpm_path is the path where the SPM structural template nifti file is stored
This file is used to :
1) Perform segmentation in the VBM pipeline
2) Compute correlation value to smoothed, warped grey matter from output of pipeline, which is stored in the vbm_qc_filename

transf_mat_path is the path to the transformation matrix used in running the reorient step of the pipeline
scan_type is the type of structural scans on which is accepted by this pipeline
FWHM_SMOOTH is the full width half maximum smoothing kernel value in mm in x,y,z directions
vbm_output_dirname is the name of the output directory to which the outputs from this pipeline are written to
vbm_qc_filename is the name of the VBM quality control text file , which is placed in vbm_output_dirname

For nifti files , it is assumed that they are T1w (T1 weighted) type of scans
FWHM_SMOOTH is an optional parameter that can be passed to this script
"""
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")


def software_check():
    """This function returns the spm standalone version installed inside the docker
    """
    spm.SPMCommand.set_mlab_paths(
        matlab_cmd=template_dict['matlab_cmd'], use_mcr=True)
    return (spm.SPMCommand().version)


def args_parser():
    parser = argparse.ArgumentParser(
        description='Run vbm pipeline on mprage T1w dicoms')
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-sub', required=True, nargs=1, type=str, help='Enter subject_id')
    required.add_argument(
        '-sess', required=True, nargs=1, type=str, help='Enter session')
    required.add_argument(
        '-json',
        '--json',
        nargs=1,
        type=str,
        required=True,
        help=
        'Enter the path to json file. Ex: -json "/data/input.json"'
    )

    return parser.parse_args()

def data_parser(args):
    if 'smoothing_x_mm' in args['options']:
         template_dict['FWHM_SMOOTH'][0]= float(args['options']['smoothing_x_mm'])
    if 'smoothing_y_mm' in args['options']:
         template_dict['FWHM_SMOOTH'][1]= float(args['options']['smoothing_y_mm'])
    if 'smoothing_z_mm' in args['options']:
        template_dict['FWHM_SMOOTH'][2] = float(args['options']['smoothing_z_mm'])

    if 'smoothing_implicit_masking' in args['options']:
        template_dict['implicit_masking']=args['options']['smoothing_implicit_masking']

    if 'reorient_params_x_mm' in args['options']:
        template_dict['reorient_params_x_mm'] = float(args['options']['reorient_params_x_mm'])
    if 'reorient_params_y_mm' in args['options']:
        template_dict['reorient_params_y_mm'] = float(args['options']['reorient_params_y_mm'])
    if 'reorient_params_z_mm' in args['options']:
        template_dict['reorient_params_z_mm'] = float(args['options']['reorient_params_z_mm'])
    if 'reorient_params_pitch' in args['options']:
        template_dict['reorient_params_pitch'] = float((args['options']['reorient_params_pitch']))
    if 'reorient_params_roll' in args['options']:
        template_dict['reorient_params_roll'] = float((args['options']['reorient_params_roll']))
    if 'reorient_params_yaw' in args['options']:
        template_dict['reorient_params_yaw'] = float((args['options']['reorient_params_yaw']))
    if 'reorient_params_x_scaling' in args['options']:
        template_dict['reorient_params_x_scaling'] = float(args['options']['reorient_params_x_scaling'])
    if 'reorient_params_y_scaling' in args['options']:
        template_dict['reorient_params_y_scaling'] = float(args['options']['reorient_params_y_scaling'])
    if 'reorient_params_z_scaling' in args['options']:
        template_dict['reorient_params_z_scaling'] = float(args['options']['reorient_params_z_scaling'])
    if 'reorient_params_x_affine' in args['options']:
        template_dict['reorient_params_x_affine'] = float(args['options']['reorient_params_x_affine'])
    if 'reorient_params_y_affine' in args['options']:
        template_dict['reorient_params_y_affine'] = float(args['options']['reorient_params_y_affine'])
    if 'reorient_params_z_affine' in args['options']:
        template_dict['reorient_params_z_affine'] = float(args['options']['reorient_params_z_affine'])

    if 'BIAS_REGULARISATION' in args['options']:
        template_dict['BIAS_REGULARISATION']=float(args['options']['BIAS_REGULARISATION'])

    if 'FWHM_GAUSSIAN_SMOOTH_BIAS' in args['options']:
        template_dict['FWHM_GAUSSIAN_SMOOTH_BIAS']=args['options']['FWHM_GAUSSIAN_SMOOTH_BIAS']

    if 'affine_regularization' in args['options']:
        template_dict['affine_regularization']=args['options']['affine_regularization']

    if 'warping_regularization' in args['options']:
        if len(args['options']['warping_regularization'])==5:
            template_dict['warping_regularization']=args['options']['warping_regularization']

    if 'sampling_distance' in args['options']:
        template_dict['sampling_distance']=float(args['options']['sampling_distance'])

    if 'mrf_weighting' in args['options']:
        template_dict['mrf_weighting']=float(args['options']['mrf_weighting'])

    if 'cleanup' in args['options']:
        template_dict['cleanup']=int(args['options']['cleanup'])

    if 'registration_template' in args['options']:
        if os.path.isfile(args['options']['registration_template']) and (str(
                ((nib.load(template_dict['tpm_path'])).shape)) == str(
            ((nib.load(args['options']['registration_template'])).shape))):
            template_dict['tpm_path'] = args['options']['registration_template']


def convert_reorientparams_save_to_mat_script():
    from pathlib2 import Path
    import shutil
    shutil.copy('/computation/convert_to_mat_file_template.m', '/computation/convert_to_mat_file.m')
    path = Path('/computation/convert_to_mat_file.m')
    text = path.read_text()
    text = text.replace('x_mm', str(template_dict['reorient_params_x_mm']))
    text = text.replace('y_mm', str(template_dict['reorient_params_y_mm']))
    text = text.replace('z_mm', str(template_dict['reorient_params_z_mm']))
    text = text.replace('pitch', str(template_dict['reorient_params_pitch']))
    text = text.replace('roll', str(template_dict['reorient_params_roll']))
    text = text.replace('yaw', str(template_dict['reorient_params_yaw']))
    text = text.replace('x_scaling', str(template_dict['reorient_params_x_scaling']))
    text = text.replace('y_scaling', str(template_dict['reorient_params_y_scaling']))
    text = text.replace('z_scaling', str(template_dict['reorient_params_z_scaling']))
    text = text.replace('x_affine', str(template_dict['reorient_params_x_affine']))
    text = text.replace('y_affine', str(template_dict['reorient_params_y_affine']))
    text = text.replace('z_affine', str(template_dict['reorient_params_z_affine']))
    path.write_text(text)
    # Run reorient.m script using spm12 standalone and Matlab MCR
    with stdchannel_redirected(sys.stderr, os.devnull):
        spm.SPMCommand.set_mlab_paths(matlab_cmd='/opt/spm12/run_spm12.sh /opt/mcr/v95 script /computation/convert_to_mat_file.m',
                                  use_mcr=True)

if __name__ == '__main__':


    try:
        # Code to check input args
        args = args_parser()
        with open(args.json[0], 'r') as f:
            json_string = json.load(f)
        data_parser(json_string)
    except Exception as e:
        sys.stderr.write('Unable to read input args/data. Error_log:'+str(e)+str(traceback.format_exc()))

    # Code to convert input reorient params options to transform.mat for reorientation, try,except,pass statement is just to supress the error its printing out.
    # Have not been able to come up with a more elegant solution yet. Already tried to supress warnings etc.
    try:
        convert_reorientparams_save_to_mat_script()
    except:
        pass

    try:
        # Check if spm is running
        with stdchannel_redirected(sys.stderr, os.devnull):
            spm_check = software_check()
        if spm_check != template_dict['spm_version']:
            raise EnvironmentError("spm unable to start in vbm docker")

        read_data = '/data'
        output_dir = '/out'

        if args.sub[0] and args.sess[0]:
            template_dict['subject'] = args.sub[0]
            template_dict['session'] = args.sess[0]
            # Check if data has nifti files
            if os.path.isfile(glob.glob(os.path.join(
                    read_data, '*.nii*'))[0]) and os.access(output_dir, os.W_OK):
                nifti_file = glob.glob(os.path.join(read_data, '*.nii*'))[0]
                computation_output = vbm_use_cases_layer.setup_pipeline(
                    data=nifti_file,
                    write_dir=output_dir,
                    data_type='nifti',
                    **template_dict)
                sys.stdout.write(computation_output)
            # Check if inputs are dicoms
            elif os.path.isdir(read_data) and os.listdir(read_data) and os.access(
                    output_dir, os.W_OK):
                dicom_header_info = os.popen('strings' + ' ' + glob.glob(
                    os.path.join(read_data, '*'))[0] + '|grep DICM').read()
                if 'DICM' in dicom_header_info:
                    computation_output = vbm_use_cases_layer.setup_pipeline(
                        data=read_data,
                        write_dir=output_dir,
                        data_type='dicoms',
                        **template_dict)
                    sys.stdout.write(computation_output)
            else:
                sys.stdout.write("No data found")
    except Exception as e:
        sys.stderr.write('Unable to pre-process data. Error_log:'+str(e)+str(traceback.format_exc()))
