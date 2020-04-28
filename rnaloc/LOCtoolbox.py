#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:43:47 2018

@author: fmueller
"""

# Imports
#import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pyplot as plt

import sys
import os
from rnaloc import annotationImporter
from rnaloc import maskGenerator
from rnaloc import FQtoolbox
from rnaloc import utils
import re
from scipy import ndimage
import json
import time
import base64
from read_roi import read_roi_file
import numpy as np


def calc_nuclear_enrichment(FQ_file, binsHist, channels={'nuclei':''}, show_plots=False, Zrange=None, dZ=2, plot_callback=None, progress_callback=None, log_callback=None):
    """
    Enrichment at the nuclear MEMBRANE
    Function uses annotations generated in FIJI and creates mask based
    on the specified parameters.
    """

    # Get input args. Has to be FIRST call!
    input_args = locals()

    # Get folder
    drive, path_and_file = os.path.splitdrive(FQ_file)
    path_results, file_results = os.path.split(path_and_file)
    file_base, ext = os.path.splitext(file_results)

    path_save = os.path.join(drive,path_results, file_base, 'NucEnvelopDist_{}'.format(time.strftime("%y%m%d-%H%M", time.localtime())))
    if not os.path.isdir(path_save):
        os.makedirs(path_save)

    # *************************************************************************
    # Load nuclear outline
    utils.log_message(f'Reading segmentation masks of nuclei', callback_fun = log_callback)

    # Generate binary masks for a selected data-set
    binaryGen = maskGenerator.BinaryMaskGenerator(erose_size=5,
                                                  obj_size_rem=500,
                                                  save_indiv=True,
                                                  progress_callback=progress_callback)

    # Import annotations for nuclei
    path_annot = os.path.join(drive,path_results,'zstack_segmentation')
    folderImporter = annotationImporter.FolderImporter(channels=channels,
                                                       data_category={'roi':''},
                                                       annot_ext='__RoiSet.zip',
                                                       progress_callback=progress_callback)
    annotDict = folderImporter.load(path_annot)

    # The generate function uses as an input the sub-dictionary for one data-category and one channel
    annotatFiles = annotDict['roi']['nuclei']
    mask_dict_nuclei = binaryGen.generate(annotatFiles)

    keys_delete = []
    for k, v in annotatFiles.items():

        # Make sure that key in present in mask, otherwise delete
        if k in mask_dict_nuclei:
            v.update(mask_dict_nuclei[k])

        else:
            keys_delete.append(k)

    # *************************************************************************
    #  Load embryo outline
    file_embryo = os.path.join(drive,path_results,'embryo_contour.roi')

    # Convert dictionary & get size of ROIS
    roi_import = read_roi_file(file_embryo)
    roi_dict = {}
    roi_dict['embryo'] = {}
    roi_dict['embryo']['type'] = roi_import['embryo_contour']['type']
    roi_dict['embryo']['pos'] = np.column_stack((roi_import['embryo_contour']['y'], roi_import['embryo_contour']['x']))

    # Assemble structure to call mask generator
    image_size = annotatFiles.values().__iter__().__next__()['image'].shape
    image_fake = np.zeros(image_size, dtype=np.uint8)

    annotat_files_embryo = {}
    annotat_files_embryo['embryo_contour'] = {}
    annotat_files_embryo['embryo_contour']['roi'] = roi_dict
    annotat_files_embryo['embryo_contour']['image'] = image_fake

    mask_dict_embryo = binaryGen.generate(annotat_files_embryo)
    mask_embryo = mask_dict_embryo['embryo_contour']['mask_fill']

    # *************************************************************************
    #  Load and analyze FQ results

    fq_dict = FQtoolbox.read_FQ_matlab(FQ_file)
    spots_all = FQtoolbox.get_rna(fq_dict)

    # Z position in pixel
    Zrna = np.divide(spots_all[:, [2]], fq_dict['settings']['microscope']['pix_z']).astype(int)

    # Other parameters for calculation
    dist_membr_RNA = np.array([])
    dist_membr_pix = np.array([])
    idx = 0

    # Loop over all z-slices
    utils.log_message(f'Loop over z-slices', callback_fun = log_callback)

    N_annot = len(annotatFiles)
    for idx_file, (k_annot, v_annot) in enumerate(annotatFiles.items()):

        # Indicate progress via callback
        if progress_callback:
            perc = int(100*(idx_file+1)/(N_annot))
            progress_callback({"task":"analyze_slices","text":f"{perc}%, {k_annot}","progress":perc})
        else:
            print(f'Slice: {k_annot}')

        # Get Z coordinate
        m = re.search('.*_Z([0-9]*)\.tif',k_annot)
        Zmask = int(m.group(1))

        # Check if outside of specified z range
        if Zrange is not None:
            if not (Zrange[0]==0 and Zrange[1]==0):
                if (Zmask < Zrange[0]) or (Zmask > Zrange[1]):
                    print(f'Z-slice outside of specified range: {Zmask}')
                    continue

        # Get z-range for loop
        if dZ==0:
            Zloop = np.logical_and(Zrna <= Zmask + dZ,Zrna >= Zmask - dZ).flatten()
        else:
            Zloop = (Zrna == Zmask).flatten()
            
        spots_loop = spots_all[Zloop,:]
        spots_loop_XY = np.divide(spots_loop[:,[0,1]],fq_dict['settings']['microscope']['pix_xy']).astype(int)

        # Distance transform
        dist_nuc_outside = ndimage.distance_transform_edt(~v_annot['mask_fill'])
        dist_nuc_inside = ndimage.distance_transform_edt(v_annot['mask_fill'])  # Negate mask
        dist_nuc = dist_nuc_outside - dist_nuc_inside

        # Indices have to be inverted to access array
        dist_nuc_RNA_loop = dist_nuc[spots_loop_XY[:, 0], spots_loop_XY[:, 1]]

        # Get distance from membrane for all pixel in the cell
        dist_membr_pix_loop = dist_nuc[mask_embryo.astype(bool)]

        # Save values
        if idx == 0:
            dist_membr_RNA = np.copy(dist_nuc_RNA_loop)
            dist_membr_pix = np.copy(dist_membr_pix_loop)
        else:
            dist_membr_RNA = np.append(dist_membr_RNA, dist_nuc_RNA_loop, axis=0)
            dist_membr_pix = np.append(dist_membr_pix, dist_membr_pix_loop, axis=0)
        idx += 1

    # *************************************************************************
    #  Load and analyze FQ results
    utils.log_message(f'Analyzing smFISH data', callback_fun = log_callback)

    x_ticks = binsHist[:-1]
    width = 0.9*np.diff(binsHist)
    width_full = np.diff(binsHist)
    center = (binsHist[:-1] + binsHist[1:]) / 2

    #histRNA_all, bins = np.histogram(dist_membr_RNA,binsHist ,density=False)
    #histPIX_all, bins = np.histogram(dist_membr_pix,binsHist ,density=False)
    #histRNA_norm = np.divide(histRNA_all,histPIX_all)
    #counts_total = np.nansum(np.multiply(histRNA_norm,width_full))
    #histRNA_norm = np.divide(histRNA_norm,counts_total)

    # Distance of all pixel from distance map
    hist_pix_all, bins = np.histogram(dist_membr_pix, binsHist, density=False)
    hist_pix_all_norm = hist_pix_all/hist_pix_all.sum()

    # Distance of all RNAs
    hist_rna_all, bins = np.histogram(dist_membr_RNA, binsHist, density=False)
    hist_rna_all_norm = hist_rna_all/hist_rna_all.sum()

    # Re-normalize RNA distances
    hist_rna_all_norm_pix = np.divide(hist_rna_all_norm, hist_pix_all_norm)
    hist_rna_all_norm_pix = np.nan_to_num(hist_rna_all_norm_pix)

    # Save file with histogram
    histo_dist = {'center': center,
                  'width': width_full,
                  'bins': binsHist,
                  'x_ticks': x_ticks,
                  'hist_pix_all': hist_pix_all,
                  'hist_pix_all_norm': hist_pix_all_norm,
                  'hist_rna_all': hist_rna_all,
                  'hist_rna_all_norm': hist_rna_all_norm,
                  'hist_rna_all_norm_pix': hist_rna_all_norm_pix}

    # Save entire analysis results as json
    input_args.pop('show_plot', None)
    input_args.pop('plot_callback', None)
    input_args.pop('progress_callback', None)
    input_args.pop('log_callback', None)

    analysis_results = {'args': input_args,
                        'histogram': histo_dist,}

    name_json = os.path.join(path_save, 'DataAnalysis.json')

    with open(name_json, 'w') as fp:
        json.dump(analysis_results, fp, sort_keys=True, indent=4, cls=utils.NumpyEncoder)

    # Save histogram of pooled data as csv
    name_csv = os.path.join(path_save, '_HistogramDistances.csv')
    histo_dist.pop('bins', None)
    csv_header = ';'.join(histo_dist.keys())
    hist_values = np.array( list(histo_dist.values())).transpose()
    np.savetxt(name_csv, hist_values, delimiter=";", fmt='%f', header=csv_header, comments='')

    # *************************************************************************
    # Plot results
    # Don't show plots when they are saved
    if not show_plots:
        plt.ioff()

    fig1, ax = plt.subplots(3,1)
    ax[0].bar(center, hist_rna_all, align='center', width=width)
    ax[0].set_xlabel('Dist to nuc')
    ax[0].set_ylabel('# RNAs')
    ax[0].set_xticks(x_ticks)
    ax[0].set_xticklabels(x_ticks.astype(int))

    ax[1].bar(center, hist_pix_all, align='center', width=width)
    ax[1].set_xlabel('Dist to nuc')
    ax[1].set_ylabel('# pixels')
    ax[1].set_xticks(x_ticks)
    ax[1].set_xticklabels(x_ticks.astype(int))

    ax[2].bar(center, hist_rna_all_norm_pix, align='center', width=width)
    ax[2].set_xlabel('Distance to nuc')
    ax[2].set_ylabel('RNA counts [a.u.]')
    ax[2].set_xticks(x_ticks)
    ax[2].set_xticklabels(x_ticks.astype(int))

    ax[0].title.set_text('RNAs')
    ax[1].title.set_text('All pixel in images')
    ax[2].title.set_text('RNAs renormalized with pixels')

    plt.tight_layout()

    # Save and close if display not required
    name_save = os.path.join(path_save, '_DistanceEnrichmentSummary.png')
    plt.savefig(name_save, dpi=300)

    if not show_plots:
        plt.close()

    if plot_callback:
        plot_callback(name_save)

    if progress_callback:
        with open(name_save, 'rb') as f:
            data = f.read()
            result = base64.b64encode(data).decode('ascii')
            img_url = 'data:image/png;base64,' + result
            progress_callback({"task": "show_results",
                               "src": img_url})


def process_file(FQ_file, bin_prop = (0,90,20), channels={'cells':'C3-'},data_category={'roi':''},annotation_extension ='__RoiSet.zip',img_extension='.tif',show_plots = False,Zrange=None,dZ=2,plot_callback=None,progress_callback=None,log_callback=None):
    '''
    Enrichment along the CELL MEMBRANE
    Function uses annotations generated in FIJI and creates mask based
    on the specified parameters.

    Args:

        plot_callback ... callback function to plot results (e.g. in ImJoy interface)

        Zrange [tuple, 2 elements]. [Optional] Tuple specifying minimum and maximum z-value that is considered
        in analysis.

        bin_prop [Tuple, 3 elements]. Specifies the bins for the histograms (min, max,delta).

    '''
    
    # Get input args. Has to be FIRST call!
    input_args = locals()
    input_args['plot_callback'] = str(input_args['plot_callback'])
    input_args['progress_callback'] = str(input_args['progress_callback'])
    input_args['log_callback'] = str(input_args['log_callback'])
        
        
    # Show function name and input arguments
    function_name = sys._getframe(  ).f_code.co_name
    utils.log_message(f"Function ({function_name}) called with:\n {str(input_args)} ", callback_fun=log_callback)

    # Make sure input args are correct - assignments with 0 can come from ImJoy
    if Zrange[0] ==0 or Zrange[1] ==0:
        Zrange = None

    if bin_prop[1] == 0 or bin_prop[1] == 0:
        bin_prop = (0,90,20)

    ## Prepare folder to save results
    drive, path_and_file = os.path.splitdrive(FQ_file)
    path_results, file_results = os.path.split(path_and_file)
    file_base, ext = os.path.splitext(file_results)

    path_save = os.path.join(drive, path_results, file_base, 'MembDist_{}'.format(time.strftime("%y%m%d-%H%M", time.localtime())))
    utils.log_message(f"Save results in folder: {path_save}", callback_fun=log_callback)
    if not os.path.isdir(path_save):
        os.makedirs(path_save)

    ## Open FQ results
    fq_dict = FQtoolbox.read_FQ_matlab(FQ_file)
    spots_all = FQtoolbox.get_rna(fq_dict)
    Zrna = spots_all[:,[18]]

    # Open annotations
    utils.log_message(f'\n== Open annotations', callback_fun = log_callback)
    
    if 'RoiSet.zip' in annotation_extension:

        path_annot = os.path.join(drive, path_results,'zstack_segmentation')
        utils.log_message(f"Opening annotation folder: {path_annot}", callback_fun=log_callback)
        
        folderImporter = annotationImporter.FolderImporter(channels=channels,
                                                           data_category=data_category,
                                                           annot_ext=annotation_extension,
                                                           progress_callback=progress_callback)
        annotDict = folderImporter.load(path_annot)
        str_size = str(annotDict['roi_size'])
        utils.log_message(f'Average roi size: {str_size}', callback_fun=log_callback)
        
        #utils.log_message(f'{annotDict['roi_size']}', callback_fun=log_callback)
        # Generate binary masks for a selected data-set
        binaryGen = maskGenerator.BinaryMaskGenerator(erose_size=5,
                                                      obj_size_rem=500,
                                                      save_indiv=True,
                                                      progress_callback=progress_callback)

        # The generate function uses as an input the sub-dictionary for one data-category and one channel
        annotatFiles = annotDict['roi']['cells']
        maskDict = binaryGen.generate(annotatFiles)

        # Use a loop and the update function to add the mask dictionary to the loaded annotation dictionary
        #utils.log_message(str(maskDict.keys()), callback_fun=log_callback)
        utils.log_message(f"Keys of mask dictionary: {str(maskDict.keys())} ", callback_fun=log_callback)
                                     
        keys_delete = []
        for k, v in annotatFiles.items():

            # Make sure that key in present in mask, otherwise delete
            if k in maskDict:
                v.update(maskDict[k])

            else:
                keys_delete.append(k)

        utils.log_message(f"Deleted keys of mask dictionary: {str(keys_delete)} ", callback_fun=log_callback)

    # Bins of histogram
    binsHist = np.arange(bin_prop[0],bin_prop[1],bin_prop[2])
    width = 0.8 * (binsHist[1] - binsHist[0])
    center = (binsHist[:-1] + binsHist[1:]) / 2

    # Other parameters for calculation
    dist_membr_RNA = np.array([])
    dist_membr_pix = np.array([])
    idx = 0

    # Loop over all z-slices
    utils.log_message(f'\n==  Loop over slices', callback_fun = log_callback)
    hist_slice ={}
    N_annot = len(annotatFiles)

    for idx_file, (k_annot, v_annot) in enumerate(annotatFiles.items()):

        # Indicate progress via callback
        if progress_callback:
            perc = int(100*(idx_file+1)/(N_annot))
            progress_callback({"task":"analyze_slices","text":f"{perc}%, {k_annot}","progress":perc})
        else:
            #print(f'Slice {idx_file+1}/{N_annot}: {k_annot}')
            utils.log_message(f'Slice {idx_file+1}/{N_annot}: {k_annot}', callback_fun = log_callback)
        # Get Z coordinate
        m = re.search('.*_Z([0-9]*)\.tif',k_annot)
        Zmask = int(m.group(1))

        # Check if outside of specified z range
        if Zrange is not None:
            if (Zmask <= Zrange[0]) or (Zmask >= Zrange[1]):
                print('Slice outside of range')
                continue

        # Get z-range for loop
        Zloop = np.logical_and(Zrna <= Zmask + dZ,Zrna >= Zmask - dZ).flatten()
        spots_loop = spots_all[Zloop,:]
        spots_loop_XY = spots_loop[:,[16, 17]].astype(int)

         # Distance transform
        dist_membr = ndimage.distance_transform_edt(~v_annot['mask_edge'])  # Negate mask

        # Indices have to be inversed to access array
        dist_membr_RNA_loop = dist_membr[spots_loop_XY[:,0],spots_loop_XY[:,1]]

        # Get distance from membrane for all pixel in the cell
        mask_all = v_annot['mask_fill'] + v_annot['mask_edge']
        dist_membr_pix_loop = dist_membr[mask_all]

        # Save values
        if idx == 0:
            dist_membr_RNA = np.copy(dist_membr_RNA_loop)
            dist_membr_pix = np.copy(dist_membr_pix_loop)
        else:
            dist_membr_RNA = np.append(dist_membr_RNA,dist_membr_RNA_loop,axis=0)
            dist_membr_pix = np.append(dist_membr_pix,dist_membr_pix_loop,axis=0)
        idx+=1

        # Calculate histograms
        histRNA, bins = np.histogram(dist_membr_RNA_loop,binsHist ,density=False)
        histpix, bins = np.histogram(dist_membr_pix_loop,binsHist ,density=False)

        histRNAnorm = histRNA/histRNA.sum()
        histpixnorm = histpix/histpix.sum()

        histRNAnormPix = np.divide(histRNAnorm,histpixnorm)
        histRNAnormPix = np.nan_to_num(histRNAnormPix)

        hist_plot = {'width':width,'center':center,'bins':bins,
                     'histRNA':histRNA,'histpix':histpix,
                     'histRNAnormPix':histRNAnormPix}
        hist_slice[f'Z{Zmask}'] = hist_plot

        # Plot results
        name_save = os.path.join(path_save,f'Z-{Zmask}.png')
        plot_results_slice(Zmask,v_annot,mask_all,spots_loop_XY,dist_membr,hist_plot,name_save,show_plots)



    # Analyze all slices
    histRNA_all, bins = np.histogram(dist_membr_RNA,binsHist ,density=False)
    histRNA_all_norm = histRNA_all/histRNA_all.sum()

    # Renormalize with pixel counts
    histpix_all, bins = np.histogram(dist_membr_pix,binsHist ,density=False)
    histpix_all_norm = histpix_all/histpix_all.sum()

    hist_RNA_all_normPix = np.divide(histRNA_all_norm,histpix_all_norm)
    hist_RNA_all_normPix = np.nan_to_num(hist_RNA_all_normPix)

    hist_plot_all = {'width':width,'center':center, 'bins':bins,
                     'histRNA_all':histRNA_all,
                     'histRNA_all_norm':histRNA_all_norm,
                     'histpix_all_norm':histpix_all_norm,
                     'hist_RNA_all_normPix':hist_RNA_all_normPix}
    name_save = os.path.join(path_save,'_DistanceEnrichmentSummary.png')
    plot_results_all(hist_plot_all,name_save,show_plots)

    if plot_callback:
        plot_callback(name_save)

    if progress_callback:
        with open(name_save, 'rb') as f:
            data = f.read()
            result = base64.b64encode(data).decode('ascii')
            imgurl = 'data:image/png;base64,' + result
            progress_callback({"task":"show_results","src":imgurl})

    # Save entire analysis results as json (remove callback functions, since those cause errors)
    input_args.pop('show_plot', None)
    input_args.pop('plot_callback', None)
    input_args.pop('progress_callback', None)
    input_args.pop('log_callback', None)
    
    analysis_results = {'args': input_args,
                        'hist_all': hist_plot_all,
                        'hist_slice': hist_slice}

    name_json = os.path.join(path_save, 'DataAll.json')

    with open(name_json, 'w') as fp:
        json.dump(analysis_results, fp,sort_keys=True, indent=4, cls=utils.NumpyEncoder)

    # Save histogram of pooled data as csv
    name_csv = os.path.join(path_save, '_HistogramPooled.csv')
    hist_plot_all.pop('bins', None)
    hist_plot_all.pop('width', None)
    csv_header  = ';'.join(hist_plot_all.keys())
    hist_values = np.array( list(hist_plot_all.values())).transpose()
    np.savetxt(name_csv, hist_values, delimiter=";",fmt='%f',header=csv_header,comments='')

    return analysis_results


def plot_results_all(hist_plot,name_save = None, show_plot=False):

    # Don't show plots when they are saved
    if not show_plot:
        plt.ioff()

    # Get parameters to plot histogram
    center = hist_plot['center']
    width = hist_plot['width']
    histRNA_all = hist_plot['histRNA_all']
    histRNA_all_norm = hist_plot['histRNA_all_norm']
    histpix_all_norm = hist_plot['histpix_all_norm']
    hist_RNA_all_normPix = hist_plot['hist_RNA_all_normPix']

    # Plot results
    fig1, ax = plt.subplots(2,2,num='Distance comparison')
    fig1.set_size_inches((8,8))

    ax[0][0].bar(center, histRNA_all, align='center', width=width)
    ax[0][0].set_xlabel('Distance cell cortex [pix]')
    ax[0][0].set_ylabel('# RNAs')
    ax[0][0].set_xticks(center)
    ax[0][0].set_xticklabels(center.astype(int))

    ax[0][1].bar(center, histRNA_all_norm, align='center', width=width/2)
    ax[0][1].set_xlabel('Distance cell cortex [pix]')
    ax[0][1].set_ylabel('Frequency')
    ax[0][1].set_xticks(center)
    ax[0][1].set_xticklabels(center.astype(int))

    ax[1][0].bar(center, histpix_all_norm, align='center', width=width)
    ax[1][0].set_xlabel('Distance cell cortex [pix]')
    ax[1][0].set_ylabel('Renormalized counts')
    ax[1][0].set_xticks(center)
    ax[1][0].set_xticklabels(center.astype(int))

    ax[1][1].bar(center, hist_RNA_all_normPix, align='center', width=width)
    ax[1][1].set_xlabel('Distance cell cortex [pix]')
    ax[1][1].set_ylabel('Renormalized counts')
    ax[1][1].set_xticks(center)
    ax[1][1].set_xticklabels(center.astype(int))


    ax[0][0].title.set_text('RNA-absolute counts')
    ax[0][1].title.set_text('RNA-normalized counts')
    ax[1][0].title.set_text('All pixel-renormalized')
    ax[1][1].title.set_text('RNA renormalized with pixels')

    plt.tight_layout()

    if name_save:
        plt.savefig(name_save,dpi=200)

    if not show_plot:
        plt.close()



def plot_results_slice(Zmask,mask,mask_all,spots_loop_XY,dist_membr,hist_plot,name_save=None,show_plot=False):

    if not show_plot:
         plt.ioff()

    # Find min and max values for plotting
    pad = 10
    indMaskAx0 = np.argwhere(mask_all.sum(axis=0))
    minAx0 = indMaskAx0[0]-pad
    maxAx0 = indMaskAx0[-1]+pad

    indMaskAx1 = np.argwhere(mask_all.sum(axis=1))
    minAx1 = indMaskAx1[0]-pad
    maxAx1 = indMaskAx1[-1]+pad

    # Set distance outside of cell to 0 for better plotting
    dist_membr_plot = np.copy(dist_membr)
    dist_membr_plot[np.logical_not(mask_all)] = 0

    # Get parameters to plot histogram
    center = hist_plot['center']
    width = hist_plot['width']
    histRNA = hist_plot['histRNA']
    histpix = hist_plot['histpix']
    histRNAnormPix = hist_plot['histRNAnormPix']

    # Generate plot
    fig1, ax = plt.subplots(2,3,num='Distance to cell membrane analysis. Z={}'.format(Zmask))
    fig1.set_size_inches((13,6))

    ax[0][0].imshow(mask['image'],cmap="hot")
    ax[0][0].get_xaxis().set_visible(False)
    ax[0][0].get_yaxis().set_visible(False)
    ax[0][0].set_xlim(minAx0, maxAx0)
    ax[0][0].set_ylim(minAx1, maxAx1)

    ax[0][1].imshow(mask_all,cmap="hot")
    ax[0][1].get_xaxis().set_visible(False)
    ax[0][1].get_yaxis().set_visible(False)
    ax[0][1].set_xlim(minAx0, maxAx0)
    ax[0][1].set_ylim(minAx1, maxAx1)

    imgdum = ax[0][2].imshow(dist_membr_plot,cmap="hot")
    ax[0][2].set_xlim(minAx0, maxAx0)
    ax[0][2].set_ylim(minAx1, maxAx1)

    ax[0][2].get_xaxis().set_visible(False)
    ax[0][2].get_yaxis().set_visible(False)
    FQtoolbox.colorbar(imgdum)

    for kROI, vROI in mask['roi'].items():
        roi_pos = vROI['pos']
        ax[0][2].plot(roi_pos[:,1],roi_pos[:,0],'b-')

    ax[0][2].scatter(spots_loop_XY[:,1],spots_loop_XY[:,0],color='g',s=4)


    ax[1][0].bar(center, histRNA, align='center', width=width)
    ax[1][0].set_xticks(center)
    ax[1][0].set_xticklabels(center.astype(int))
    ax[1][0].set_xlabel('Distance from cell cortex [pixel]')
    ax[1][0].set_ylabel('# RNAs')

    ax[1][1].bar(center, histpix, align='center', width=width)
    ax[1][1].set_xticks(center)
    ax[1][1].set_xticklabels(center.astype(int))
    ax[1][1].set_xlabel('Distance from cell cortex [pixel]')
    ax[1][1].set_ylabel('# pixels')

    ax[1][2].bar(center, histRNAnormPix, align='center', width=width)
    ax[1][2].set_xticks(center)
    ax[1][2].set_xticklabels(center.astype(int))
    ax[1][2].set_xlabel('Distance from cell cortex [pixel]')
    ax[1][2].set_ylabel('Renormalized frequency')

    # Set titles
    ax[0][0].title.set_text('Cell cortex')
    ax[0][1].title.set_text('Cell mask')
    ax[0][2].title.set_text('Distance transform')

    ax[1][0].title.set_text('RNAs')
    ax[1][1].title.set_text('All pixel')
    ax[1][2].title.set_text('Renormalized RNA distance')

    plt.tight_layout()

    if name_save:
        plt.savefig(name_save,dpi=200)

    if not show_plot:
        plt.close()
