# Khufai_Formation_Oman_EBSD_PIC_SEM_data_Shuram_excursion

This repository contains analytical data supporting the study of dolomite fabrics in ooids formed during the onset of the Shuram excursion from the uppermost Khufai Formation, Huqf region, Oman. Two additional dolomite samples were compared from the Buah and Birba Formations. Samples were prepared as thin sections (0.25 μm thickness) and epoxy-mounted blocks. EBSD, SEM-EDS, and EPMA samples were polished with colloidal silica and carbon-coated prior to analysis, PIC samples were polished to 50 nm with Al₂O₃. Acid etching was completed on polished thin sections with 10% HCl for 30 seconds prior to SEM analysis. XRD analysis was completed with small powder quantities with a zero-background silicon chip.

Repository Structure

README.md: This file, providing an overview of the repository and dataset descriptions.

Data Folders:

data/EBSD/

Electron Backscatter Diffraction (EBSD) data: Contains raw and processed EBSD datasets from three samples: Khufai Formation ooids (MD6 258.6m), spherulitic dolomite from conophyton stromatolite (SB1_CON), and dolomite cements from pisolite complex (SB1_ONC). Data collected using ThermoFisher Helios Hydra 5 SEM with Oxford Instruments Symmetry S3 EBSD detector at 25 kV, 51 nA beam current. Three maps were generated with a with 250 nm step size and one map was generated with a 500 nm step size. Maps provide full 3-D crystal orientations at micron-scale resolution. Data processed with MTEX toolbox for MATLAB, providing quantitative crystallographic fabric analysis.

data/EBSD/code/

MATLAB scripts for EBSD data analysis and figure generation: Contains MATLAB scripts utilizing the open-source MTEX toolbox for crystallographic analysis. Scripts generate inverse pole figure maps, pole figures, crystal orientation statistics, and c-axis to elongation angle calculations. Analysis includes texture strength measurements and comparative crystallographic fabric assessment.
MTEX Software: The analysis uses the MTEX toolbox for MATLAB, available at: https://mtex-toolbox.github.io/ (Bachmann et al., 2010)

data/EPMA/

Electron Probe Microanalysis data: Quantitative elemental spot analyses and maps from JEOL JXA-8200 electron microprobe (15 kV, 20 nA, 1 μm beam). Includes Ca, Mg, Fe, Si, Sr, and Mn elemental maps and spot analysis concentrations with detection limits of 177-300 ppm, calibrated using calcite, dolomite, and rhodochrosite standards. Elemental maps show concentric banding patterns in ooid cortices, silica cements, nuclei types, etc.

data/PIC/

Polarization-dependent Imaging Contrast (PIC) mapping data: High-resolution (20 nm resolution, 56 nm pixels) crystal orientation maps using synchrotron-based soft X-ray imaging with variable linear polarization. Maps crystal c-axis orientations in 3D for dolomite samples, providing nanoscale details of crystal morphologies and orientations. Data were collected on beamline 11.0.1.1 at the Advanced Light Source, Lawrence Berkeley National Laboratory, using X-ray PhotoEmission Electron spectroMicroscopy (X-PEEM). Synchrotron-based technique using variable polarization soft X-rays to map crystal c-axis orientations at 20 nm resolution. Complements EBSD with nanoscale details of crystal arrangements.

data/SEM_acid_etched/

Scanning Electron Microscopy images of acid-etched surfaces: Secondary electron images of samples etched with 10% HCl for 30 seconds using Hitachi TM4000II microscope (10 kV, 8.5 mm working distance). Images reveal crystal morphologies, boundaries, and microtextures at sub-micron resolution across ooids from sample MD6 258.6m and comparative samples SB1_CON and SB1_ONC.

data/SEM_BSE/

Backscattered Electron (BSE) imaging data: BSE images showing compositional contrast and crystal relationships. Collected simultaneously with EDS mapping at 15 kV accelerating voltage. Images highlight dolomite crystal arrangements, cement relationships, porosity, and inclusion distributions.

data/photomicrographs/

Optical microscopy images: Plane-polarized light (PPL) and cross-polarized light (XPL) photomicrographs taken with an AxioImager M2 showing ooid textures, crystal fabrics, nuclei types, and cement relationships. Images document radial crystal arrangements, cortical banding, and optical properties of dolomite crystals. In XPL, it is possible to distinguish silica cements from the glass slide.

data/stitched_photomicrograph/

Large-scale composite microscopy images: Complete thin-section panoramic images created from multiple 5× magnification photomicrographs using an AxioImager M2. Provides a comprehensive view of sample textures, ooid distribution, and spatial relationships between analytical regions for sample MD6 258.6m.

data/UV_Fluorescence/

Ultraviolet fluorescence microscopy images: UV-induced fluorescence images taken on the AxioImager M2 revealing organic matter distribution, cement generations, and diagenetic features. Helps distinguish primary depositional fabrics from secondary alteration.

data/Field_photos/

Field context images: Outcrop photographs showing stratigraphic relationships, depositional structures, and sample collection locations from Mukhaibah Dome and Khufai Dome localities in the Huqf region, Oman.

data/XPL_movies/

Cross-polarized light rotation movies: Time-lapse microscopy videos showing optical behavior during stage rotation on a rotating stage installed on an AxioImager M2. Demonstrates extinction patterns, crystal orientations, and optical continuity in dolomite ooids and comparative samples. The stage is always rotated clockwise.  

data/XRD/

X-ray Diffraction patterns: Powder diffraction data from PANalytical X'Pert Pro (Cu radiation, 45 kV, 40 mA). Scanned 5-70° 2θ with 0.008° step size. Includes dolomite ordering calculations from (015)/(110) peak intensity ratios and phase identification.

Citation
Wilcots, J., Bergmann, K.B., Gilbert, P.U.P.A., Cross, A., 2025. Khufai Formation Oman EBSD PIC SEM Data - Shuram Excursion. Zenodo. https://doi.org/10.5281/zenodo.16415784

References
Bachmann, F., Hielscher, R., & Schaeben, H. (2010). Texture analysis with MTEX--free and open source software toolbox. Solid State Phenomena, 160, 63-68.

Data Availability
This repository is public and contains all supporting analytical data for the study of Ediacaran dolomite ooid crystallography and its implications for interpreting the Shuram excursion onset.
