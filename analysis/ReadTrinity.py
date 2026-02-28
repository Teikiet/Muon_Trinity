#Read the eventio file
from eventio import IACTFile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
def read_cherenkov_hits(input_file, plot=True):
    with IACTFile(input_file) as f:
            events = iter(f)
            event = next(events)
            #print(event.header)
            azimuth_rad = (event.header['azimuth']) #event azimuth
            zenith_rad = (event.header['zenith']) #event zenith
            X = event.photon_bunches[0]['x'] #in cm
            Y = event.photon_bunches[0]['y'] #in cm
            T = event.photon_bunches[0]['time'] #in ns
            cos_X = event.photon_bunches[0]['cx'] #cosine of the x direction
            cos_Y = event.photon_bunches[0]['cy'] #cosine of the y direction
            wavelength = event.photon_bunches[0]['wavelength'] #in nm
            zem = event.photon_bunches[0]['zem']  # in cm from Earth's center to the photon emission point
            zem = zem * 1e-2   # Convert from cm to m
            zem =-6371e3 + zem  # Convert to height above Earth's surface
            ##print("zenith", np.degrees(zenith_rad))
            ##print("azimuth", np.degrees(azimuth_rad))
            ##print("max zem", np.max(zem))
    Slant = (zem-2944) / np.cos(zenith_rad)  # Calculate slant depth
    Weight = np.ones(len(X))*5
    cos_Z = np.sqrt(1 - cos_X**2 - cos_Y**2)  # Calculate cosine of the z direction
    Theta_z = np.degrees(np.arccos(cos_Z))
    Phi = np.degrees(np.arctan2(cos_Y, cos_X))
    ##print(np.nanmean(Theta_z))
    ##print(np.nanmean(Phi))
    if plot:
        plt.figure(figsize=(6, 5))

        h = plt.hist2d(
        X/1e2,
        Y/1e2,
        bins=100,          # adjust binning as needed
        cmap="viridis",
        #range=[[-15, 15], [-15, 15]],
        #range=[[-0.65, 0.65], [-0.65, 0.65]],
        weights=Weight,
        norm=LogNorm()  # use logarithmic color scale
        )
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.title("Cherenkov Hits: 2D Histogram (x vs y)")
        plt.colorbar(h[3], label="Counts")
        plt.tight_layout()
        plt.show()
    return X, Y, T, cos_X, cos_Y, wavelength, zem, Theta_z, Phi, Slant, Weight


#GROPT Reader
import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
# Open ROOT file
#file = uproot.open("/uufs/chpc.utah.edu/common/home/u1520754/Detector_Sims/Test_Dir/GROPT/Tilt_91.560000/diffuse_BSM_1.000_E_8.000000_Ze_45.000000_Az_0.000000_0.root")  # Local file example
# List all contents
def plot_GROPT(in_file, plot=True):
    file = uproot.open(in_file)
    ##print(file.keys())

    # Access a tree
    tree = file['T1;1']  # Replace with your tree name

    # See available branches
    #print(tree.keys())

    # Read all data into arrays
    data = tree.arrays()

    # Convert to pandas DataFrame (optional)
    
    df = tree.arrays(library="pd")

    # Or to numpy
    arrays = tree.arrays(library="np")
    cosx = arrays['photonDcosX'][0]
    cosy = arrays['photonDcosY'][0]
    photonX_np = arrays['photonX'][0]
    photonY_np = arrays['photonY'][0]
    T = arrays['time'][0]
    if plot:
        plt.hist2d(photonX_np, photonY_np, bins=100, range=[[-50,50],[-50,50]], cmap='viridis')#,norm=LogNorm())
        plt.colorbar(label='Counts')
        plt.xlabel('Photon X Position (mm)')
        plt.ylabel('Photon Y Position (mm)')
        plt.title('2D Histogram of Photon Positions')
        plt.show()
    return cosx, cosy, photonX_np, photonY_np, T

#CARE reader function

import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Open ROOT file
#file = uproot.open("/uufs/chpc.utah.edu/common/home/u1520754/Detector_Sims/Test_Dir/CARE/Tilt_91.560000/diffuse_BSM_1.000_E_8.000000_Ze_45.000000_Az_0.000000_0.root")

def plot_CARE(in_file, plot=True):
    file = uproot.open(in_file)
    # List all contents
    #print(file.keys())

    # Access a tree
    tree = file['Events/T0;1']  # Replace with your tree name

    # See available branches
    #print(tree.keys())

    # Load the data
    arrays = tree.arrays(library="np")

    # Assuming you read `vPEInPixel`, and pixel positions are stored or known
    vPEInPixel = arrays["vPEInPixel"][0]  # Photoelectrons for each pixel
    #vPEInPixel = np.flip(vPEInPixel, axis=0)  # Assuming the data is flipped
    # Pixel configuration from your dataset
    # Pixel configuration from your dataset
    pixel_config = [ #I need to regenerate it again in the future
        # Pixel ID, x, y (In mm, extracted from your input data)
        (0, -46.875, -46.875), (1, -46.875, -40.625), 
        (2, -46.875, -34.375), (3, -46.875, -28.125),
        (4, -40.625, -46.875), (5, -40.625, -40.625), 
        (6, -40.625, -34.375), (7, -40.625, -28.125),
        (8, -34.375, -46.875), (9, -34.375, -40.625),
        (10, -34.375, -34.375), (11, -34.375, -28.125),
        (12, -28.125, -46.875), (13, -28.125, -40.625),
        (14, -28.125, -34.375), (15, -28.125, -28.125),
        (16, -21.875, -46.875), (17, -21.875, -40.625),
        (18, -21.875, -34.375), (19, -21.875, -28.125),
        (20, -15.625, -46.875), (21, -15.625, -40.625),
        (22, -15.625, -34.375), (23, -15.625, -28.125),
        (24, -9.375, -46.875), (25, -9.375, -40.625),
        (26, -9.375, -34.375), (27, -9.375, -28.125),
        (28, -3.125, -46.875), (29, -3.125, -40.625),
        (30, -3.125, -34.375), (31, -3.125, -28.125),
        (32, 3.125, -46.875), (33, 3.125, -40.625),
        (34, 3.125, -34.375), (35, 3.125, -28.125),
        (36, 9.375, -46.875), (37, 9.375, -40.625),
        (38, 9.375, -34.375), (39, 9.375, -28.125),
        (40, 15.625, -46.875), (41, 15.625, -40.625),
        (42, 15.625, -34.375), (43, 15.625, -28.125),
        (44, 21.875, -46.875), (45, 21.875, -40.625),
        (46, 21.875, -34.375), (47, 21.875, -28.125),
        (48, 28.125, -46.875), (49, 28.125, -40.625),
        (50, 28.125, -34.375), (51, 28.125, -28.125),
        (52, 34.375, -46.875), (53, 34.375, -40.625),
        (54, 34.375, -34.375), (55, 34.375, -28.125),
        (56, 40.625, -46.875), (57, 40.625, -40.625),
        (58, 40.625, -34.375), (59, 40.625, -28.125),
        (60, 46.875, -46.875), (61, 46.875, -40.625),
        (62, 46.875, -34.375), (63, 46.875, -28.125),
        (64, -46.875, -21.875), (65, -46.875, -15.625),
        (66, -46.875, -9.375), (67, -46.875, -3.125),
        (68, -40.625, -21.875), (69, -40.625, -15.625),
        (70, -40.625, -9.375), (71, -40.625, -3.125),
        (72, -34.375, -21.875), (73, -34.375, -15.625),
        (74, -34.375, -9.375), (75, -34.375, -3.125),
        (76, -28.125, -21.875), (77, -28.125, -15.625),
        (78, -28.125, -9.375), (79, -28.125, -3.125),
        (80, -21.875, -21.875), (81, -21.875, -15.625),
        (82, -21.875, -9.375), (83, -21.875, -3.125),
        (84, -15.625, -21.875), (85, -15.625, -15.625),
        (86, -15.625, -9.375), (87, -15.625, -3.125),
        (88, -9.375, -21.875), (89, -9.375, -15.625),
        (90, -9.375, -9.375), (91, -9.375, -3.125),
        (92, -3.125, -21.875), (93, -3.125, -15.625),
        (94, -3.125, -9.375), (95, -3.125, -3.125),
        (96, 3.125, -21.875), (97, 3.125, -15.625),
        (98, 3.125, -9.375), (99, 3.125, -3.125),
        (100, 9.375, -21.875), (101, 9.375, -15.625),
        (102, 9.375, -9.375), (103, 9.375, -3.125),
        (104, 15.625, -21.875), (105, 15.625, -15.625),
        (106, 15.625, -9.375), (107, 15.625, -3.125),
        (108, 21.875, -21.875), (109, 21.875, -15.625),
        (110, 21.875, -9.375), (111, 21.875, -3.125),
        (112, 28.125, -21.875), (113, 28.125, -15.625),
        (114, 28.125, -9.375), (115, 28.125, -3.125),
        (116, 34.375, -21.875), (117, 34.375, -15.625),
        (118, 34.375, -9.375), (119, 34.375, -3.125),
        (120, 40.625, -21.875), (121, 40.625, -15.625),
        (122, 40.625, -9.375), (123, 40.625, -3.125),
        (124, 46.875, -21.875), (125, 46.875, -15.625),
        (126, 46.875, -9.375), (127, 46.875, -3.125),
        (128, -46.875, 3.125), (129, -46.875, 9.375),
        (130, -46.875, 15.625), (131, -46.875, 21.875),
        (132, -40.625, 3.125), (133, -40.625, 9.375),
        (134, -40.625, 15.625), (135, -40.625, 21.875),
        (136, -34.375, 3.125), (137, -34.375, 9.375),
        (138, -34.375, 15.625), (139, -34.375, 21.875),
        (140, -28.125, 3.125), (141, -28.125, 9.375),
        (142, -28.125, 15.625), (143, -28.125, 21.875),
        (144, -21.875, 3.125), (145, -21.875, 9.375),
        (146, -21.875, 15.625), (147, -21.875, 21.875),
        (148, -15.625, 3.125), (149, -15.625, 9.375),
        (150, -15.625, 15.625), (151, -15.625, 21.875),
        (152, -9.375, 3.125), (153, -9.375, 9.375),
        (154, -9.375, 15.625), (155, -9.375, 21.875),
        (156, -3.125, 3.125), (157, -3.125, 9.375),
        (158, -3.125, 15.625), (159, -3.125, 21.875),
        (160, 3.125, 3.125), (161, 3.125, 9.375),
        (162, 3.125, 15.625), (163, 3.125, 21.875),
        (164, 9.375, 3.125), (165, 9.375, 9.375),
        (166, 9.375, 15.625), (167, 9.375, 21.875),
        (168, 15.625, 3.125), (169, 15.625, 9.375),
        (170, 15.625, 15.625), (171, 15.625, 21.875),
        (172, 21.875, 3.125), (173, 21.875, 9.375),
        (174, 21.875, 15.625), (175, 21.875, 21.875),
        (176, 28.125, 3.125), (177, 28.125, 9.375),
        (178, 28.125, 15.625), (179, 28.125, 21.875),
        (180, 34.375, 3.125), (181, 34.375, 9.375),
        (182, 34.375, 15.625), (183, 34.375, 21.875),
        (184, 40.625, 3.125), (185, 40.625, 9.375),
        (186, 40.625, 15.625), (187, 40.625, 21.875),
        (188, 46.875, 3.125), (189, 46.875, 9.375),
        (190, 46.875, 15.625), (191, 46.875, 21.875),
        (192, -46.875, 28.125), (193, -46.875, 34.375),
        (194, -46.875, 40.625), (195, -46.875, 46.875),
        (196, -40.625, 28.125), (197, -40.625, 34.375),
        (198, -40.625, 40.625), (199, -40.625, 46.875),
        (200, -34.375, 28.125), (201, -34.375, 34.375),
        (202, -34.375, 40.625), (203, -34.375, 46.875),
        (204, -28.125, 28.125), (205, -28.125, 34.375),
        (206, -28.125, 40.625), (207, -28.125, 46.875),
        (208, -21.875, 28.125), (209, -21.875, 34.375),
        (210, -21.875, 40.625), (211, -21.875, 46.875),
        (212, -15.625, 28.125), (213, -15.625, 34.375),
        (214, -15.625, 40.625), (215, -15.625, 46.875),
        (216, -9.375, 28.125), (217, -9.375, 34.375),
        (218, -9.375, 40.625), (219, -9.375, 46.875),
        (220, -3.125, 28.125), (221, -3.125, 34.375),
        (222, -3.125, 40.625), (223, -3.125, 46.875),
        (224, 3.125, 28.125), (225, 3.125, 34.375),
        (226, 3.125, 40.625), (227, 3.125, 46.875),
        (228, 9.375, 28.125), (229, 9.375, 34.375),
        (230, 9.375, 40.625), (231, 9.375, 46.875),
        (232, 15.625, 28.125), (233, 15.625, 34.375),
        (234, 15.625, 40.625), (235, 15.625, 46.875),
        (236, 21.875, 28.125), (237, 21.875, 34.375),
        (238, 21.875, 40.625), (239, 21.875, 46.875),
        (240, 28.125, 28.125), (241, 28.125, 34.375),
        (242, 28.125, 40.625), (243, 28.125, 46.875),
        (244, 34.375, 28.125), (245, 34.375, 34.375),
        (246, 34.375, 40.625), (247, 34.375, 46.875),
        (248, 40.625, 28.125), (249, 40.625, 34.375),
        (250, 40.625, 40.625), (251, 40.625, 46.875),
        (252, 46.875, 28.125), (253, 46.875, 34.375),
        (254, 46.875, 40.625), (255, 46.875, 46.875)
    ]


    # Extract the x, y coordinates and ids for all pixels
    pixel_ids = [p[0] for p in pixel_config]
    x_coords = [p[1] for p in pixel_config]
    y_coords = [p[2] for p in pixel_config]

    # Normalize coordinates to start from (0,0)
    min_x, min_y = min(x_coords), min(y_coords)
    max_x, max_y = max(x_coords), max(y_coords)
    pixel_pitch = 6.25  # Spacing between pixels

    # Compute grid dimensions
    num_cols = int((max_x - min_x) / pixel_pitch) + 1
    num_rows = int((max_y - min_y) / pixel_pitch) + 1

    # Map pixel IDs to their respective (row, col) indices in the grid
    pixel_map = {}
    for pix_id, x, y in pixel_config:
        col = int((x - min_x) / pixel_pitch)
        row = int((y - min_y) / pixel_pitch)
        pixel_map[pix_id] = (row, col)

    # Create an empty 2D matrix for the pixel signals
    camera_matrix = np.zeros((num_rows, num_cols))

    # Populate the matrix using `vPEInPixel` values
    for pix_id, pe_count in enumerate(vPEInPixel):
        row, col = pixel_map[pix_id]
        camera_matrix[row, col] = pe_count

    # Plot the camera signal
    if plot:
        plt.imshow(np.flip(camera_matrix, axis=(0)), cmap='viridis', extent=(min_x, max_x, min_y, max_y))
        plt.colorbar(label='Photoelectron Counts')
        plt.title('Camera Pixel Signals')
        plt.xlabel('Pixel X Position (mm)')
        plt.ylabel('Pixel Y Position (mm)')
        plt.show()
    # all 256 FADC traces
    all_vFADCTrace = []
    for i in range(256):
        try:
            # Get FADC trace for pixel i
            trace_key = f"vFADCTraces{i}"
            if trace_key in arrays:
                vFADCTrace = arrays[trace_key][0]
                all_vFADCTrace.append(vFADCTrace)
        except:
            continue
    all_vFADCTrace = np.array(all_vFADCTrace)
    time_ns = np.arange(len(all_vFADCTrace[0])) * 10 - 170 # 10 ns per sample and trigger at 170 ns (sample 17)
    dc_2_pe = 24.1
    pedestal_dc = 500
    total_trace = np.sum(all_vFADCTrace - pedestal_dc, axis=0)
    max_trace = np.max(all_vFADCTrace - pedestal_dc, axis=0)
    total_pe = np.clip(total_trace, 0, None) / dc_2_pe
    max_pe = np.clip(max_trace, 0, None) / dc_2_pe
    return camera_matrix, time_ns, total_pe, max_pe, arrays, pixel_config, min_x, max_x, min_y, max_y, pixel_pitch
