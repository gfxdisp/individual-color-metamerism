The experiment data is organized as follows:

1). (Obs responses) 
Individual participant responses are collected in .csv files named Corrected_???.mat, where ??? is the three letter observer ID assigned during the experiment. 
The files contain a matrix that is 6x22. Each row contains data for a different display pair. These pairs are, by row numbers:

1 - (Display) 1 and (Display) 2
2 - 1 and 3
3 - 1 and 4
4 - 2 and 3
5 - 2 and 4
6 - 3 and 4

To create a matrix containing these display pairs use Matlab command nchoosek(1:4,2).

Display numbers are:
1: LG C2
2: Sony X310
3: Samsung Laser Projector
4: ASUS VG246

Each column contains data for a specific color. The colors are presented in the same order as in the calibration files. Colors 2 and 5 are not evaluated on display
4 because they are outside of the gamut.
The first 11 columns contain the data from the first half of the experiment, where observers evaluated all colors on all display pairs once. Columns 12 to 22 contain
data from the second half of the experiment, where each evaluation was repeated. The evaluations match positions in both submatrices, i.e. position (1,1) contains the
answer to the first evaluation of color 1 on displays 1 and 2, position (1,12) contains the second evaluation of the same color on the same displays.

2). (Spectra)
Folder contains .csv and .mat files containing the same data: the spectra of individual colors used in the experiment. Each matrix has 11 rows, one for each color used in the experiment. Each column contains spectral radiance measured on the specified display at each wavelength, from 380nm to 780 nm at 1nm intervals. The wavelengths are specified in the wavelengths.(csv|mat) file. Data is in absolute radiance values (W * sr^-1 *m^-2).

3). (AuxData)
Additional data needed to calculate the observer metamerism indices. 

ciexyz31_1.csv - 1931 CIE Standard 2 deg observer tabulated at 1nm

IndividualObs_2deg - Contains 1000 simulated observers from Yuta Asano's PhD dissertation, tabulated from 390nm to 780nm at 5nm intervals.