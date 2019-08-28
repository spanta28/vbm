#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer defines the nodes of pre-processing pipeline
"""

import nipype.pipeline.engine as pe
import nipype.interfaces.spm as spm
spm.terminal_output = 'file'
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
from nipype.interfaces.spm.utils import DicomImport


## 1 Dicom convert node & settings ##
class Dicom_convert:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=DicomImport(), name='converter')
        self.node.inputs.in_files = template_dict['dicom_dir']
        self.node.inputs.output_dir = template_dict['dicom_output_dir']
        self.node.run()


## 2 Reorientation node & settings ##
class Reorient:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.ApplyTransform(), name='reorient')
        self.node.inputs.mat = template_dict['transf_mat_path']
        self.node.inputs.paths = template_dict['spm_path']


## 3 Segementation Node and settings ##
class Segment:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.NewSegment(), name='segmentation')
        self.node.inputs.paths = template_dict['spm_path']
        self.node.inputs.channel_info = (
            template_dict['BIAS_REGULARISATION'],
            template_dict['FWHM_GAUSSIAN_SMOOTH_BIAS'], (False, False))
        self.node.inputs.affine_regularization=template_dict['affine_regularization']
        self.node.inputs.warping_regularization=template_dict['warping_regularization']
        self.node.inputs.sampling_distance=template_dict['sampling_distance']
        self.node.inputs.mrf_weighting=template_dict['mrf_weighting']
        self.node.inputs.cleanup=template_dict['cleanup']

def transform_list(normalized_class_images):
    return [each[0] for each in normalized_class_images]


class List_Normalized_Images:
    def __init__(self):
        self.node = pe.Node(
            interface=Function(
                input_names='normalized_class_images',
                output_names='list_norm_images',
                function=transform_list),
            name='List_normalized_images')


## 4 Smoothing Node & Settings ##
class Smooth:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.Smooth(), name='smoothing')
        self.node.inputs.paths = template_dict['spm_path']
        self.node.inputs.fwhm = template_dict['FWHM_SMOOTH']
        self.node.inputs.implicit_masking=template_dict['implicit_masking']

## 5 Datsink Node that collects segmented, smoothed files and writes to temp_write_dir ##
class Datasink:
    def __init__(self):
        self.node = pe.Node(interface=DataSink(), name='sinker')