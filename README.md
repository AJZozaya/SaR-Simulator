# 'SARrawSim' is a Matlab-based program for generating raw data from a general pulsed linear frequency-modulated synthetic aperture radar (SAR) operating in stripmap mode.

# Author: Alfonso J. Zozaya S.

# Code Summary Description

The code consists of six (06) sections. The first section, Probing Scene Definition and Drawing, defines the probing area. In section 2, Radar Parameters Definition, radar parameters, including the transmitted pulse and time span for A-scope visualization, are established. Section 3, named Random Targets Generation, randomly generates a set of NofT point targets. In the following section, Radar Travel Animation and Raw Data Generation, a loop animates the radar travel while checking for targets inside the half-power beam-width antenna pattern (equivalent to detecting targets inside the antenna footprint). Subsequently, indexes of targets within the antenna footprint are extracted (equivalent to target ranging). The echo is then created, Gaussian noise is added, and the result is visualized in A-scope type. The next section, Raw Image, renders a raw data image using the pcolor Matlab function. The final section, Completion of Metadata Data, consolidates metadata and raw data in a mat file named rawdata.mat for further processing.

# Usage

The code can be run as it is in its present form, but the user, according to their own requirements, can edit the sections: Probing Scene Definition and Drawing, Radar Parameters Definition, and Random Targets Generation. In any case, the radar parameters have to be defined based on slant range dr and cross-range du resolutions. The following two reccomendations are made: (1) From the range resolution dr, secondary values are calculated, including the bandwidth B = c/2dr, where c represents the speed of light. Additionally, the pulse duration is determined as tau=100/B, the chirp rate as k= B/tau, the fast-time sampling frequency as Fs=osfB, and its corresponding period as Ts= 1/Fs, where osf denotes the oversampling factor. The pulse to be transmitted is already defined based on parameters tau and B. (2) From the cross-resolution du=lambda0/2Dtheta, the horizontal half-power antenna beamwidth Dtheta, and carrier frequency f0 are determined. Then, from the Doppler frequency bandwidth BD=2v Dtheta/lambda0, the platform speed v and the pulse repetition interval tR=1/osfBD are specified. All values established at this stage constitute the metadata.

# 'Focusing' is a Matlab-based program for focusing raw data from a general pulsed linear frequency-modulated synthetic aperture radar (SAR) operating in stripmap mode using the Range-Doppler Algorithm (RDA).

# Author: Alfonso J. Zozaya S.

# Code Summary Description

This code loads the file 'data.mat' generated from 'SARrawSim' and applies the Range-Doppler Algorithm (RDA) to the raw data. The RDA is executed based on metadata contained in the file, resulting in a single complex radar image. The script successively calls the functions: 'eR=RaNGeC(P, e, NofS)', '[ERA, ERCA, fD]=RCMC(eR, NofR, dr, tR, r, f0, v)', and 'eRA=azimuthC(ERCA, v, fD, r, lambda0)'. These functions, included in the code itself, perform the core processes of the RDA, i.e., focusing the data in range, correcting range-cell migration, and focusing the data in azimuth, respectively. The outcome is the conversion of the raw data into a single look complex (SLC) radar image.

# Usage

Run the 'Focusing' script after running 'SARrawSim'.

# Contact

a.zozayas@utem.cl
