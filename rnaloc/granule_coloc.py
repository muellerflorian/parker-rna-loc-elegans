# Imports
import os
import numpy as np
from skimage import io
import matplotlib.pyplot as plt
from skimage.feature import blob_dog, blob_log
from math import sqrt
import pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import ConvexHull, Delaunay
from skimage.morphology import dilation
from matplotlib import colors
from matplotlib_venn import venn2
import seaborn as sns

# Own libraries
from rnaloc.utils import create_folder, create_mask_sphere, cmap_create_random
from rnaloc.FQtoolbox import read_FQ_matlab, get_rna


def analyze_coloc(clusters_label, granules_label, df_clusters, df_granules, img_granules, img_clusters, show_plot=True, path_save=''):
    """
    Analyze coloc of DB-scan cluster, and detected blob granules.

    Args:
        clusters_label (3d np array): 3D image with cluster labels.
        granules_label (3d np array): 3D image with granule labels.
        df_clusters (Pandas DF): summarizes cluster properties.
        df_granules (Pandas DF): summarizes granules properties.
        img_granules (2D np array): 2D image of granules.
        img_clusters (2D np array): 2D image of clusters.
        show_plot (Bool): Show plot.
        path_save (string): path to save results. No results will be saved if empty.
    Returns:
        df_clusters (Pandas DF): summarizes cluster properties.
        df_granules (Pandas DF): summarizes granules properties.
    """

    # **** Analyze overlap

    # Overlap mask
    mask_overlap = clusters_label.astype('bool') & granules_label.astype('bool')

    # Overlap of RNA clusters
    df_clusters['coloc'] = False
    RNA_cluster_overlap = clusters_label[mask_overlap]
    RNA_cluster_idx_overlap = np.unique(RNA_cluster_overlap)
    idx_coloc = np.isin(df_clusters['index'], RNA_cluster_idx_overlap)
    df_clusters.loc[idx_coloc, 'coloc'] = True
    print(f'RNA clusters with granules: {RNA_cluster_idx_overlap}')

    # Overlap of P-granules
    df_granules['coloc'] = False
    granules_overlap = granules_label[mask_overlap]
    granules_idx_overlap = np.unique(granules_overlap)
    idx_coloc = np.isin(df_granules['index'], granules_idx_overlap)
    df_granules.loc[idx_coloc, 'coloc'] = True
    print(f'Granules with RNA clusters: {granules_idx_overlap}')

    # Save data-frames
    if path_save:
        df_clusters.to_csv(os.path.join(path_save, f'RNA_clusters.csv'), sep=';')
        df_granules.to_csv(os.path.join(path_save, f'P_granules.csv'), sep=';')

    # *** Plot results of granule and RNA cluster detection: color-coded masks

    # Create mask with granules being co-localized (pixel value 2), or not (pixel value 1)
    img_granules_label_MIP = granules_label.max(axis=0)
    img_granules_label_MIP_coloc = np.copy(img_granules_label_MIP)
    img_granules_label_MIP_coloc[img_granules_label_MIP > 0] = 1

    for row in df_granules.itertuples(index=True, name='Pandas'):
        index = getattr(row, 'index')
        coloc = getattr(row, 'coloc')
        if coloc:
            img_granules_label_MIP_coloc[img_granules_label_MIP == index] = 2

    # Create mask with RNA clusters being co-localized (pixel value 2), or not (pixel value 1)
    img_RNA_clusters_label_MIP = clusters_label.max(axis=0)
    img_RNA_clusters_label_MIP_coloc = np.copy(img_RNA_clusters_label_MIP)
    img_RNA_clusters_label_MIP_coloc[img_RNA_clusters_label_MIP > 0] = 1

    for row in df_clusters.itertuples(index=True, name='Pandas'):
        index = getattr(row, 'index')
        coloc = getattr(row, 'coloc')
        if coloc:
            img_RNA_clusters_label_MIP_coloc[img_RNA_clusters_label_MIP == index] = 2

    # Show figure
    cmap_random = cmap_create_random()  # Random look-up table
    cmap_3colors = colors.ListedColormap(('black', 'magenta', 'green'))

    # Show plot
    if not show_plot:
        plt.ioff()

    fig, ax = plt.subplots(2, 3)
    fig.set_size_inches((10, 10))

    ax[0][0].imshow(np.log10(img_granules), cmap="Greys_r")
    ax[0][0].get_xaxis().set_visible(False)
    ax[0][0].get_yaxis().set_visible(False)
    ax[0][0].set_title('Granules [log 10 scale]')

    ax[0][1].imshow(np.log10(img_granules), cmap="Greys_r")
    ax[0][1].imshow(img_granules_label_MIP, cmap=cmap_random)
    ax[0][1].get_xaxis().set_visible(False)
    ax[0][1].get_yaxis().set_visible(False)
    ax[0][1].set_title('Detected granules')

    ax[0][2].imshow(img_granules_label_MIP_coloc, cmap=cmap_3colors)
    ax[0][2].get_xaxis().set_visible(False)
    ax[0][2].get_yaxis().set_visible(False)
    ax[0][2].set_title('Granules: green-coloc; magenta-not coloc')

    ax[1][0].imshow(img_clusters, cmap="Greys_r")
    ax[1][0].get_xaxis().set_visible(False)
    ax[1][0].get_yaxis().set_visible(False)
    ax[1][0].set_title('FISH')

    ax[1][1].imshow(img_clusters, cmap="Greys_r")
    ax[1][1].imshow(img_RNA_clusters_label_MIP, cmap=cmap_random)
    ax[1][1].get_xaxis().set_visible(False)
    ax[1][1].get_yaxis().set_visible(False)
    ax[1][1].set_title('Detected clusters')

    ax[1][2].imshow(img_RNA_clusters_label_MIP_coloc, cmap=cmap_3colors)
    ax[1][2].get_xaxis().set_visible(False)
    ax[1][2].get_yaxis().set_visible(False)
    ax[1][2].set_title('Clusters: green-coloc; magenta-not coloc')

    plt.tight_layout()
    if path_save:
        plt.savefig(os.path.join(path_save, f'COMPARE_granules_clusters.png'), dpi=300)

    if show_plot:
        plt.show()
    else:
        plt.close()

    # *** Summarize as a Venn diagram
    fig, ax = plt.subplots(1, 1, num='Venn diagram')
    fig.set_size_inches((4, 4))
    plt.clf()

    venn2(subsets=(df_clusters.shape[0], df_granules.shape[0], df_clusters['coloc'].sum()),
          set_labels=('RNA clusters', 'P-granules'))
    plt.tight_layout()
    if path_save:
        plt.savefig(os.path.join(path_save, f'VENN_diagram_overlap.png'), dpi=300)

    if show_plot:
        plt.show()
    else:
        plt.close()

    # *** Show RNA cluster/P-granule size as a function of being co-localized or not
    fig, ax = plt.subplots(1, 2)
    fig.set_size_inches((4, 4))

    # Channel 1 as box plot
    sns.boxplot(x="coloc", y="n_rna", data=df_clusters, ax=ax[0])
    ax[0].set_title('RNA clusters')
    ax[0].set_xlabel('Colocalized?')
    ax[0].set_ylabel('# RNAs per clusters')
    ax[0].set_xticklabels(['NO', 'YES'])

    # Channel 2 as box plot
    sns.boxplot(x="coloc", y="int_mean", data=df_granules, ax=ax[1])
    ax[1].set_title('P granules')
    ax[1].set_xlabel('Colocalized?')
    ax[1].set_ylabel('Average intensity')
    ax[1].set_xticklabels(['NO', 'YES'])
    plt.tight_layout()
    if path_save:
        plt.savefig(os.path.join(path_save, f'Intensity_vs_clustering.png'), dpi=300)

    if show_plot:
        plt.show()
    else:
        plt.close()


    return df_clusters, df_granules


def dbscan_analyze(clusters, img, pos, show_plot=True, path_save=''):
    """
    Analyze results of DB scan.

    Args:
        clusters (1D array: For each position, index of corresponding cluster. -1 if not assigned.
        img (3d np array): 3D where clusters where called.
        pos (3D np array): Detected 3D positions in pixel.
        show_plot (Bool): Show plot.
        path_save (string): path to save results. No results will be saved if empty.
    Returns:
        img_clusters_label (3D np): image with labels of detected clusters.
        df_clusters (Pandas DF): summarizes cluster properties.

    """
    #%% Loop over RNA cluster results and create mask image

    # Get index of all clusters
    clusters_idx, cluster_size = np.unique(clusters, return_counts=True)
    clusters_idx = np.delete(clusters_idx, np.where(clusters_idx == -1))   # Delete the ones without cluster assignment

    # Prepare variables
    img_clusters_label = np.zeros(img.shape).astype('int16')  # Image with masks of RNA clusters
    cluster_properties = np.zeros((clusters_idx.shape[0], 5))           # Numpy array to store cluster properties

    for loop_idx, cluster_idx in np.ndenumerate(clusters_idx):

        # Get RNA positions in cluster
        rna_in_cluster = (clusters == cluster_idx)
        rna_loop = pos[rna_in_cluster, :]
        print(f'Analyzing cluster {cluster_idx} with {rna_loop.shape[0]} RNAs')

        # Hull with triangulation
        try:
            hull = ConvexHull(rna_loop)
            tri = Delaunay(rna_loop[hull.vertices, :],qhull_options = 'QL')

            # Create grid and check which one are in hull
            ymax, xmax, zmax = rna_loop.max(axis=0)
            ymin, xmin, zmin = rna_loop.min(axis=0)

            yy,xx,zz = np.meshgrid(np.arange(ymin, ymax+1, 1),
                                   np.arange(xmin, xmax+1, 1),
                                   np.arange(zmin, zmax+1, 1))

            points_test = np.empty((yy.flatten().shape[0],3))
            points_test[:, 0] = yy.flatten()
            points_test[:, 1] = xx.flatten()
            points_test[:, 2] = zz.flatten()

            in_cluster = tri.find_simplex(points_test)
            index_in = np.where(in_cluster > -1)

            # Order in Python is ZYX
            img_clusters_label[
                            points_test[index_in, 2].astype('int'),
                            points_test[index_in, 0].astype('int'),
                            points_test[index_in, 1].astype('int')] = cluster_idx+1  # Cluster idx starts at 0
        except:
            print('No triangulated surface could be created. Will use RNA positions directly.')
            # Order in Python is ZYX
            img_clusters_label[
                                   rna_loop[:, 2].astype('int'),
                                   rna_loop[:, 0].astype('int'),
                                   rna_loop[:, 1].astype('int')] = cluster_idx+1  # Cluster idx starts at 0

        cluster_properties[loop_idx, [0, 1, 2]] = rna_loop.mean(axis=0).astype('int')
        cluster_properties[loop_idx, 3] = rna_loop.shape[0]
        cluster_properties[loop_idx, 4] = loop_idx[0]+1

    # Create dataframe with results
    df_clusters = pd.DataFrame({
                    'index': cluster_properties[:, 4].astype('int'),
                    'z': cluster_properties[:, 0].astype('int'),
                    'y': cluster_properties[:, 1].astype('int'),
                    'x': cluster_properties[:, 2].astype('int'),
                    'n_rna':  cluster_properties[:, 3].astype('int'),
    })

    # Opening with small ball to smooth clusters
    img_clusters_label = dilation(img_clusters_label, selem=create_mask_sphere(1))

    if path_save:
        io.imsave(os.path.join(path_save, 'img_clusters_label.tif'), img_clusters_label)

        io.imsave(os.path.join(path_save, 'img_clusters_mask.tif'), img_clusters_label.astype('bool').astype('int16'))

    return img_clusters_label, df_clusters


def dbscan_apply(pos, img_2d, eps, min_samples, show_plot=True, name_save=''):
    """
    Perform DB scan.

    Args:
        pos (3D np array): Positions
        img_2d (2d np array): 2D image to show overlay
        eps (float): Maximum distance.
        min_samples (integer): Minimum number of elements per cluster.
        show_plot (Bool): Show plot.
        name_save (string): name to save results image. No results will be saved if empty.
    Returns:
        clusters (1D array): For each position, index of corresponding cluster. -1 if not assigned.
    """

    # Perform DB scan
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    clusters = dbscan.fit_predict(pos)
    cmap_random = cmap_create_random()  # Random look-up table

    # ** Plot the cluster assignment

    # Show plot
    if not show_plot:
        plt.ioff()

    fig, ax = plt.subplots(1, 3, figsize=(6, 8))
    ax[0].imshow(img_2d, cmap="Greys_r")
    ax[0].get_xaxis().set_visible(False)
    ax[0].get_yaxis().set_visible(False)
    ax[0].set_title('FISH')

    ax[1].scatter(pos[:, 1], pos[:, 0], s=0.5)
    plt.xlabel("X")
    plt.ylabel("Y")
    ax[1].get_xaxis().set_visible(False)
    ax[1].get_yaxis().set_visible(False)
    ax[1].set_aspect('equal', 'box')
    ax[1].invert_yaxis()
    ax[1].set_title('All detection')

    ax[2].scatter(pos[:, 1], pos[:, 0], c=clusters, cmap=cmap_random, s=1)
    plt.xlabel("X")
    plt.ylabel("Y")
    ax[2].set_aspect('equal', 'box')
    ax[2].get_xaxis().set_visible(False)
    ax[2].get_yaxis().set_visible(False)
    ax[2].invert_yaxis()
    ax[2].set_title('DB-scan results')

    plt.tight_layout()
    if name_save:
        plt.savefig(name_save, bbox_inches='tight', pad_inches=0, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()

    return clusters


def detect_blobs_3d(img, blobs_threshold, blobs_method,show_plot=False, path_save=''):
    """
    Detect blobs in 3D images.
    https://scikit-image.org/docs/dev/auto_examples/features_detection/plot_blob.html

    Args:
        img (2D np array): 2D image
        blobs_threshold (float): detection threshold
        blobs_method (string): blob detection methods. Supported are
                log: Laplacian of Gaussian
                dog: Difference of Gaussian.
        show_plot (Bool): Show plot.
        path_save (string): path to save results. No results will be saved if empty.
    Returns:
        img_granules_label (3D np array): label image of detected granules.
        df_granules (Pandas df): summarizing detected granules.

    """

    # ** Actual blob detection
    if blobs_method == 'log':
        blobs_3D = blob_log(img, min_sigma=2, max_sigma=10, num_sigma=10, threshold=blobs_threshold)
    if blobs_method == 'dog':
        blobs_3D = blob_dog(img, max_sigma=10, threshold=blobs_threshold)

    # Calculate radius
    blobs_3D[:, 3] = blobs_3D[:, 3] * sqrt(2)

    # ** Summarize and get 3D images of labels and mask
    granules_intensity = np.zeros((blobs_3D.shape[0], 1))
    granules_idx = np.zeros((blobs_3D.shape[0], 1))

    n_z, n_y, n_x = img.shape
    img_granules_label = np.zeros((n_z, n_y, n_x)).astype('uint16')

    #print('Image size)')
    #print(img_granules_label.shape)
    for idx, blob in enumerate(blobs_3D):
        z, y, x, r = blob[[0, 1, 2, 3]].astype('int')
        print(f'R,Z,Y,X: {r}, {z},{y},{x}')

        # Get default range in image and mask
        z_min = z-r-1
        z_min_mask = 0
        z_max = z+r
        z_max_mask = 2 * r + 1

        y_min = y-r-1
        y_min_mask = 0
        y_max = y+r
        y_max_mask = 2 * r + 1

        x_min = x-r-1
        x_min_mask = 0
        x_max = x + r
        x_max_mask = 2 * r + 1

        #print(f'Z before: range image: {z_max - z_min}, range mask: {z_max_mask - z_min_mask}')
        #print(f'Y before: range image: {y_max - y_min}, range mask: {y_max_mask - y_min_mask}')
        #print(f'X before: range image: {x_max - x_min}, range mask: {x_max_mask - x_min_mask}')

        # Correct max values
        if z_max > n_z-1:
            z_max_mask = z_max_mask - (z_max-n_z)
            z_max = n_z

        if y_max > n_y-1:
            y_max_mask = y_max_mask - (y_max-n_y)
            y_max = n_y

        if x_max > n_x-1:
            x_max_mask = x_max_mask - (x_max-n_x)
            x_max = n_x

        # Correct min values
        if z_min < 0:
            z_min_mask = 0 + (-z_min)
            z_min = 0

        if y_min < 0:
            y_min_mask = 0 + (-y_min)
            y_min = 0

        if x_min < 0:
            x_min_mask = 0 + (-x_min)
            x_min = 0

        #print(f'Z: range image: {z_max - z_min}, range mask: {z_max_mask - z_min_mask}')
        #print(f'Y: range image: {y_max - y_min}, range mask: {y_max_mask - y_min_mask}')
        #print(f'X: range image: {x_max - x_min}, range mask: {x_max_mask - x_min_mask}')


        # Spherical mask
        sphere_loop = create_mask_sphere(r)
        sphere_loop_add = sphere_loop[z_min_mask:z_max_mask, y_min_mask:y_max_mask, x_min_mask:x_max_mask]

        # Label image
        img_granules_label[z_min:z_max, y_min:y_max, x_min:x_max] = (idx + 1) * sphere_loop_add

        # Mask and average intensity of region
        img_mask_loop = np.zeros(img.shape).astype('uint16')
        img_mask_loop[z_min:z_max, y_min:y_max, x_min:x_max] = sphere_loop_add
        granules_intensity[idx, 0] = img[img_mask_loop.astype('bool')].mean()
        granules_idx[idx, 0] = idx + 1

    img_granules_mask = img_granules_label.astype('bool').astype('int16')

    if path_save:
        io.imsave(os.path.join(path_save, 'img_granules_label.tif'), img_granules_label)
        io.imsave(os.path.join(path_save, 'img_granules_mask.tif'), img_granules_mask)

    # Data-frame with results
    df_granules = pd.DataFrame({
        'index': granules_idx[:, 0].astype('int'),
        'z': blobs_3D[:, 0].astype('int'),
        'y': blobs_3D[:, 1].astype('int'),
        'x': blobs_3D[:, 2].astype('int'),
        'radius': blobs_3D[:, 3],
        'int_mean': granules_intensity[:, 0].astype('int'),

    })

    # ** Plot results
    img_MIP = img.max(axis=0)

    # Show plot
    if not show_plot:
        plt.ioff()

    fig, ax = plt.subplots(1, 2, figsize=(6, 6))
    ax[0].imshow(np.log10(img_MIP), cmap="Greys_r")
    ax[0].set_axis_off()
    ax[0].set_title('P-granules [log10]')

    ax[1].imshow(np.log10(img_MIP), cmap="Greys_r")
    for blob in blobs_3D:
        z, y, x, r = blob
        c = plt.Circle((x, y), r, color='red', linewidth=2, fill=False)
        ax[1].add_patch(c)
    ax[1].set_axis_off()
    ax[1].set_title(f'Detection in 3D [{blobs_method}-{blobs_threshold}]')

    plt.tight_layout()

    if path_save:
        plt.savefig(os.path.join(path_save, f'3D_blob_detection__{blobs_method}-{blobs_threshold}.png'), dpi=300)

    if show_plot:
        plt.show()
    else:
        plt.close()

    return img_granules_label, df_granules


def detect_blobs_2d(img, blobs_threshold, blobs_method,show_plot=False, name_save=''):
    """
    Detect blobs in 2D images.
    https://scikit-image.org/docs/dev/auto_examples/features_detection/plot_blob.html
    
    Args:
        img (2D np array): 2D image
        blobs_threshold (float): detection threshold
        blobs_method (string): blob detection methods. Supported are
                log: Laplacian of Gaussian
                dog: Difference of Gaussian.
        show_plot (Bool): Show plot. 
        name_save (string): name to save plot.
    Returns:
        path_project (dictionary: updated dictionoary different project folders

    """

    if blobs_method == 'log':
        blobs_2D = blob_log(img, min_sigma=2, max_sigma=10, num_sigma=10, threshold=blobs_threshold)
    if blobs_method == 'dog':
        blobs_2D = blob_dog(img, max_sigma=10, threshold=blobs_threshold)

    # Calculate radius
    blobs_2D[:, 2] = blobs_2D[:, 2] * sqrt(2)

    # Show plot
    if not show_plot:
        plt.ioff()

    fig, ax = plt.subplots(1, 2, figsize=(6, 8))
    ax[0].imshow(np.log10(img), cmap="Greys_r")
    ax[0].set_axis_off()
    ax[0].set_title('P-granules [log10]')

    ax[1].imshow(np.log10(img), cmap="Greys_r")
    for blob in blobs_2D:
        y, x, r = blob
        c = plt.Circle((x, y), r, color='red', linewidth=2, fill=False)
        ax[1].add_patch(c)
    ax[1].set_axis_off()
    ax[1].set_title(f'Detection in 2D [{blobs_method}-{blobs_threshold}]')

    plt.tight_layout()

    if name_save:
        plt.savefig(name_save, bbox_inches='tight', pad_inches=0, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def get_folders(path_project, channel_identifier):
    """
    Determine folders for projects.
    * path_FQ = '/Volumes/PILON_HD2/fmueller/Documents/Data/smFISH_collaborations/Dylan/190729_PGranules_Detection/data/L4440_wDMP04_glh-1-21/C2/results_GMM'
    * path_img_RNA      >> By default one level up from FQ results.
    * path_img_granules >> By default two levels up from FQ results.
    * path_save         >> By default two levels up from FQ results in a folder called 'coloc_granules_CHANNELIDENT'

    Args:
        path_project (dictionary): contains different project folders
        channel_identifier (tuple): first element is identifier of FISh image, second of granules channels

    Returns:
        path_project (dictionary: updated dictionoary different project folders

    """

    # Get different folders
    path_FQ = path_project['FQ']
    path_save = path_project['save']
    path_img_RNA = path_project['img_granules']
    path_img_granules = path_project['img_RNA']

    # Get path to save images
    if ('..' in path_save) | (not path_save):
        path_save = os.path.abspath(os.path.join(path_FQ, '..', '..', 'coloc_granules', f'Channel_{channel_identifier[0]}'))
        create_folder(path_save)

    # Get path with images
    if ('..' in path_img_RNA) | (not path_img_RNA):
        path_img_RNA = os.path.abspath(os.path.join(path_FQ, '..'))

    if ('..' in path_img_granules) | (not path_img_granules):
        path_img_granules = os.path.abspath(os.path.join(path_FQ, '..', '..'))

    # Check all directories
    if not os.path.exists(path_img_granules):
        raise NotADirectoryError(f'Directory with granule images does not exist: {path_img_granules}',path_img_granules)

    if not os.path.exists(path_img_RNA):
        raise NotADirectoryError(f'Directory with FISH image does not exist: {path_img_RNA}',path_img_RNA)

    path_project = {
        'FQ': path_FQ,
        'save': path_save,
        'img_granules': path_img_granules,
        'img_RNA': path_img_RNA,
    }

    return path_project


def read_data(file_FQ, channel_identifier, path_project):
    """
    Read all data.

    Args:
        file_FQ (string): filename of FQ results.
        channel_identifier (tuple): first element is identifier of FISh image, second of granules channels
        path_project (dictionary): contains different project folders
    Returns:
        path_project (dictionary: updated dictionoary different project folders

    """

    # Read FQ results
    file_load = os.path.join(path_project['FQ'], file_FQ)
    print(f'Loading file: {file_load}')
    fq_dict = read_FQ_matlab(file_load)

    # Position in pixel
    spots_all = get_rna(fq_dict)
    spots_all_pix = np.copy(spots_all)
    spots_all_pix[:, [0, 1]] = np.divide(spots_all[:, [0, 1]], fq_dict['settings']['microscope']['pix_xy']).astype(int)
    spots_all_pix[:, 2] = np.divide(spots_all[:, 2], fq_dict['settings']['microscope']['pix_z']).astype(int)
    spots_all_pix = spots_all_pix.astype('int')

    # Load FISH image
    file_img_FISH = fq_dict['file_names']['smFISH']
    try:
        img_FISH = io.imread(os.path.join(path_project['img_RNA'], file_img_FISH))
    except IOError as e:
        print(f'Unable to open image of granules {file_img_FISH}')

    # Granule image
    file_img_granules = file_img_FISH.replace(channel_identifier[0], channel_identifier[1])
    if file_img_granules == file_img_FISH:
        raise UserWarning(f'Filename for granules identical with filename for FISH: {file_img_FISH}', file_img_FISH)

    try:
        img_granules = io.imread(os.path.join(path_project['img_granules'], file_img_granules))
    except IOError as e:
        print(f'Unable to open image of granules {file_img_granules}')


    return spots_all, spots_all_pix, img_FISH, img_granules