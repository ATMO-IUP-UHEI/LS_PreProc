#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:08:14 2022

@author: guanter
"""

import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import georasters as gr   
import spectral.io.envi as envi

def read_ang (tree, lab_ang):
    
    for fact in tree.iter(tag = lab_ang):
        try:
            ang1 = np.float(fact.find('upper_left').text)
            ang2 = np.float(fact.find('upper_right').text)
            ang3 = np.float(fact.find('lower_left').text)
            ang4 = np.float(fact.find('lower_right').text)
        except:
            pass

    return np.mean([ang1, ang2, ang3, ang4])


def read_xml (file_xml, swir_flg):
    
    tree = ET.parse(file_xml) #read in the XML

    
    if swir_flg:
        lab_str_m = 'swir'    
    else:    
        lab_str_m = 'vnir'    
    
    wl_center = []
    wl_fwhm = []
    gain_arr = []
    offs_arr = []
    for fact in tree.iter(tag = 'bandID'):
        try:
            wvl = np.float(fact.find('wavelengthCenterOfBand').text)
            fwhm = np.float(fact.find('FWHMOfBand').text)
            gain = np.float(fact.find('GainOfBand').text)
            offset = np.float(fact.find('OffsetOfBand').text)

            wl_center = np.append(wl_center, wvl)
            wl_fwhm = np.append(wl_fwhm, fwhm)
            gain_arr = np.append(gain_arr, gain)
            offs_arr = np.append(offs_arr, offset)
        except:
            pass

    #plt.figure()
    #plt.plot(gain_arr)
    #plt.figure()
    #plt.plot(offs_arr)
    
    
    for fact in tree.iter(tag = lab_str_m + 'ProductQuality'):  
        try:
            num_bd_sp = (np.array(fact.find('numChannelsExpected').text)).astype(int)
        except:
            pass

    num_bd = len(wl_center)
    
    if lab_str_m == 'swir':
        wl_center = wl_center[num_bd-num_bd_sp:]
        wl_fwhm = wl_fwhm[num_bd-num_bd_sp:]
        gain_arr = gain_arr[num_bd-num_bd_sp:]
        offs_arr = offs_arr[num_bd-num_bd_sp:]
    else:   
        wl_center = wl_center[0:num_bd_sp]
        wl_fwhm = wl_fwhm[0:num_bd_sp]
        gain_arr = gain_arr[0:num_bd_sp]
        offs_arr = offs_arr[0:num_bd_sp]

    
    sza = 90. - read_ang (tree, 'sunElevationAngle')
    saa = read_ang (tree, 'sunAzimuthAngle')
    vaa = read_ang (tree, 'sceneAzimuthAngle')
    vza = read_ang (tree, 'acrossOffNadirAngle')

    for fact in tree.iter(tag = 'specific'):
        try:
            hsf = np.float(fact.find('meanGroundElevation').text)
        except:
            pass

    return wl_center, wl_fwhm, hsf, sza, saa, vaa, vza, gain_arr, offs_arr 

def import_enmap (file_img_env, file_img_xml, swir_flg, tiff_flg):
    
    #sc_coef = 1.e-2
    sc_coef = 1.e+3
    
    wl_center, wl_fwhm, hsf, sza, saa, vaa, vza, gain_arr, offs_arr = read_xml (file_img_xml, swir_flg)
    
    if tiff_flg:
        im_tmp = np.transpose(gr.load_tiff(file_img_env), (2, 1, 0))
    else:    
        img = envi.open(file_img_env)
        im_tmp = np.transpose(img.open_memmap(writeable = True), (1, 0, 2))

    ncols, nrows, num_bd = im_tmp.shape

    im_ini = np.zeros([ncols, nrows, num_bd])        
    for bd in range(0, len(wl_center)):        
        im_ini[:, :, bd] = (im_tmp[:, :, bd] * gain_arr[bd] + offs_arr[bd])* sc_coef
        
    im_tmp = 0
    
    return im_ini, wl_center, wl_fwhm, hsf, sza, saa, vza, vaa

swir_flg = True

path_dat = '/mnt/LARS_NAS/Hyperspectral/EnMAP/turkmenistan/dims_op_oc_oc-en_700289847_1/ENMAP.HSI.L1B/ENMAP-HSI-L1BDT0000004147_03-2022-10-02T07:48:37.859_2022-10-07T12:41:04/ENMAP01-____L1B-DT0000004147_20221002T074837Z_003_V010106_20221007T093455Z/'
img_str = 'ENMAP01-____L1B-DT0000004147_20221002T074837Z_003_V010106_20221007T093455Z'          
file_img_env = path_dat + img_str + '-SPECTRAL_IMAGE_SWIR' + '.TIF'
file_xml = path_dat + img_str + '-METADATA.XML'

ltoa_img, wl_center, wl_fwhm, hsf, sza, saa, vza, vaa = import_enmap (file_img_env, file_xml, swir_flg, tiff_flg)
