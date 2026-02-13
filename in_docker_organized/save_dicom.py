# from __future__ import print_function
# import sys, time, os
# import numpy as np
# import pydicom
# import nibabel as nib
# from pydicom.dataset import Dataset, FileDataset
# from pydicom.uid import CTImageStorage
# from tqdm import tqdm

# pixel_dtypes = {
#     "float64": np.float64,
#     "float32": np.float32,
#     "uint16": np.uint16
# }

# def generate_series_metadata():
#     """Generate consistent UIDs and metadata for the series"""
#     series_uid = pydicom.uid.generate_uid()
#     study_uid = pydicom.uid.generate_uid()
#     frame_uid = pydicom.uid.generate_uid()
    
#     return {
#         'SeriesInstanceUID': series_uid,
#         'StudyInstanceUID': study_uid,
#         'FrameOfReferenceUID': frame_uid,
#         'SeriesDate': time.strftime("%Y%m%d"),
#         'SeriesTime': time.strftime("%H%M%S"),
#         'SeriesDescription': 'Converted from NIFTI',
#         'PatientName': 'SinglePatient',
#         'PatientID': 'SYNTHETIC123',
#         'StudyID': 'Study1',
#         'PatientBirthDate': '20000101',
#         'PatientSex': 'O',
#         'AccessionNumber': 'ACC0001'
#     }

# def analyze_data_range(arr):
#     """Analyze array data to determine optimal scaling parameters"""
#     data_min = float(np.min(arr))
#     data_max = float(np.max(arr))
#     data_mean = float(np.mean(arr))
#     data_std = float(np.std(arr))
    
#     # Print detailed statistics
#     print(f"Array statistics - min: {data_min:.6f}, max: {data_max:.6f}, mean: {data_mean:.6f}, std: {data_std:.6f}")
    
#     return data_min, data_max, data_mean, data_std

# def convert_slice(arr, output_path, index=0, series_metadata=None, spacing=None):
#     """Convert single slice to DICOM with enhanced float64 handling"""
#     # Analyze the data range
#     data_min, data_max, data_mean, data_std = analyze_data_range(arr)
    
#     # Create file meta
#     file_meta = pydicom.Dataset()
#     file_meta.MediaStorageSOPClassUID = CTImageStorage
#     file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()
#     file_meta.ImplementationClassUID = pydicom.uid.PYDICOM_IMPLEMENTATION_UID
#     file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian

#     # Create dataset
#     ds = FileDataset(output_path, {}, 
#                     file_meta=file_meta,
#                     preamble=b"\0" * 128)
    
#     # Enhanced float64/float32 handling for optimal display
#     if arr.dtype in [np.float64, np.float32]:
#         print(f"Converting {arr.dtype} data to DICOM format")
        
#         # Calculate optimal scaling parameters based on data range
#         data_range = data_max - data_min
        
#         if data_range > 0:
#             if data_min >= -1024 and data_max <= 3071:  # Standard CT range
#                 # Use direct conversion for CT-like data
#                 print("Data appears to be in CT Hounsfield units range")
#                 arr_scaled = arr.astype(np.int16)
#                 slope = 1.0
#                 intercept = 0.0
#             elif data_max <= 1.0 and data_min >= 0.0:
#                 # Normalized data (0-1 range) - scale to full 16-bit range
#                 print("Data appears to be normalized (0-1 range)")
#                 arr_scaled = (arr * 65535).astype(np.uint16)
#                 slope = 1.0 / 65535
#                 intercept = 0.0
#                 ds.PixelRepresentation = 0  # unsigned
#             else:
#                 # Custom scaling to maximize precision
#                 # Calculate optimal scaling to use most of int16 range
#                 # Avoid using max range to prevent overflow issues
#                 print("Using custom scaling for float data")
#                 target_range = 65000.0  # slightly less than full int16 range (65536)
                
#                 if abs(data_min) > abs(data_max):
#                     # Data has larger negative values
#                     slope = abs(data_min) / 32000.0
#                 else:
#                     # Data has larger positive values or is balanced
#                     slope = data_range / target_range
                
#                 intercept = data_min
                
#                 # Apply the scaling transformation with rounding for better precision
#                 arr_scaled = np.round((arr - intercept) / slope).astype(np.int16)
                
#                 # Double-check conversion
#                 max_scaled = np.max(arr_scaled)
#                 min_scaled = np.min(arr_scaled)
#                 print(f"Scaled data range: {min_scaled} to {max_scaled}")
                
#                 # Verify that original values can be recovered
#                 sample_orig = arr.flat[0]
#                 sample_scaled = arr_scaled.flat[0]
#                 sample_recovered = sample_scaled * slope + intercept
#                 print(f"Sample value - original: {sample_orig:.6f}, recovered: {sample_recovered:.6f}, diff: {abs(sample_orig-sample_recovered):.6f}")
#         else:
#             # Handle flat data (all values the same)
#             print("Warning: Data has no range (all values are the same)")
#             arr_scaled = np.zeros(arr.shape, dtype=np.int16)
#             slope = 1.0
#             intercept = data_min
#     else:
#         # Handle integer types with appropriate casting
#         if arr.dtype == np.uint16:
#             # Already in DICOM-compatible format
#             arr_scaled = arr.copy()
#             slope = 1.0
#             intercept = 0.0
#         else:
#             # Scale other integer types to int16 range if needed
#             data_min = float(np.min(arr))
#             data_max = float(np.max(arr))
            
#             if data_max > 32767 or data_min < -32768:
#                 # Need rescaling
#                 data_range = data_max - data_min
#                 slope = data_range / 65000.0
#                 intercept = data_min
#                 arr_scaled = np.round((arr - intercept) / slope).astype(np.int16)
#             else:
#                 # Can be directly converted to int16
#                 arr_scaled = arr.astype(np.int16)
#                 slope = 1.0
#                 intercept = 0.0
    
#     # Set essential DICOM attributes
#     ds.SOPClassUID = CTImageStorage
#     ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
#     ds.Modality = "CT"
    
#     # Set image details
#     ds.SamplesPerPixel = 1
#     ds.PhotometricInterpretation = "MONOCHROME2"
#     ds.Rows = arr_scaled.shape[0]
#     ds.Columns = arr_scaled.shape[1]
#     ds.BitsAllocated = 16
#     ds.BitsStored = 16
#     ds.HighBit = 15
#     ds.PixelRepresentation = 1  # 1 for signed
#     ds.PixelData = arr_scaled.tobytes()
    
#     # Set DICOM rescale parameters - crucial for proper display
#     ds.RescaleSlope = str(slope)
#     ds.RescaleIntercept = str(intercept)
    
#     # Set appropriate window values for better default display
#     if data_range > 0:
#         # For custom data: use 4 standard deviations for window width
#         if data_std > 0:
#             # Center on mean, width covers 4 standard deviations
#             center = data_mean
#             width = data_std * 4.0
#         else:
#             # Fallback: use min/max
#             center = (data_max + data_min) / 2
#             width = data_max - data_min
            
#         # For CT-like data: use standard CT windowing
#         if data_min >= -1024 and data_max <= 3071:
#             # Default to soft tissue window
#             center = 40
#             width = 400
        
#         print(f"Setting window center: {center:.2f}, window width: {width:.2f}")
#         ds.WindowCenter = str(center)
#         ds.WindowWidth = str(width)
    
#     # Set consistent series information
#     if series_metadata:
#         for key, value in series_metadata.items():
#             setattr(ds, key, value)
    
#     # Set instance-specific attributes
#     ds.InstanceNumber = str(index)
    
#     # Set spatial position information
#     if spacing:
#         ds.PixelSpacing = [str(spacing[0]), str(spacing[1])]
#         ds.SliceThickness = str(spacing[2])
#         ds.SpacingBetweenSlices = str(spacing[2])
        
#         # Set slice location (essential for correct 3D positioning)
#         ds.SliceLocation = str(index * spacing[2])
#         ds.ImagePositionPatient = [str(0), str(0), str(index * spacing[2])]
#     else:
#         ds.PixelSpacing = ["1", "1"]
#         ds.SliceThickness = "1"
#         ds.ImagePositionPatient = [str(0), str(0), str(index)]
    
#     # Patient orientation (standard axial orientation)
#     ds.ImageOrientationPatient = ['1', '0', '0', '0', '1', '0']
#     ds.PatientPosition = 'HFS'  # Head First-Supine
    
#     ds.save_as(output_path)

# def nifti2dicom_1file(nifti_path, out_dir, is_2d=False):
#     """Convert single nifti file to DICOM series"""
#     print(f"Converting {nifti_path} to DICOM...")
    
#     # Load nifti file
#     nifti_img = nib.load(nifti_path)
#     nifti_array = nifti_img.get_fdata()
    
#     # Print information about the array
#     print(f"Nifti array shape: {nifti_array.shape}, dtype: {nifti_array.dtype}")
#     print(f"Min value: {np.min(nifti_array)}, Max value: {np.max(nifti_array)}")
    
#     # Get spacing information from nifti header
#     spacing = list(nifti_img.header.get_zooms())
#     if len(spacing) < 3:
#         spacing.append(1.0)  # Add slice thickness for 2D
    
#     print(f"Pixel spacing: {spacing}")
    
#     # Generate consistent series metadata
#     series_metadata = generate_series_metadata()
    
#     if is_2d or len(nifti_array.shape) < 3:
#         # Handle 2D image case
#         if len(nifti_array.shape) > 2:
#             print("Warning: 3D array found but is_2d=True, using first slice only")
#             nifti_array = nifti_array[:,:,0]
        
#         convert_slice(nifti_array, os.path.join(out_dir, 'slice0000.dcm'), 0, 
#                      series_metadata, spacing)
#         print(f"Created 1 DICOM file")
#     else:
#         # Handle 3D image case
#         print(f"Creating {nifti_array.shape[2]} DICOM slices...")
        
#         for slice_idx in tqdm(range(nifti_array.shape[2])):
#             output_path = os.path.join(out_dir, f'slice{slice_idx:04d}.dcm')
#             convert_slice(nifti_array[:,:,slice_idx], output_path, slice_idx,
#                         series_metadata, spacing)
        
#         print(f"Created {nifti_array.shape[2]} DICOM files")

# # Main script
# if len(sys.argv) < 4:
#     print(f"Usage: python {__file__} <input_nifti> <output_directory> <is_2D> [pixel_type]")
#     sys.exit(1)

# input_nifti = sys.argv[1]
# output_dir = sys.argv[2]
# is_2D = bool(int(sys.argv[3]))

# os.makedirs(output_dir, exist_ok=True)
# nifti2dicom_1file(input_nifti, output_dir, is_2D)

