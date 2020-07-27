# Making music with Matlab

This is a Music generator program written in Matlab and built using Matlab R2020a, tested successfully on macOS Mojave.
The application generates a 10 seconds audio clip for the Mortal Kombat theme song.  

## Usage Instructions

1. To Compile the .m file : open the .m file on matlab and hit run
2. To change the name of the output files generated simply replace in the code. 
3. The output songs are stored in the default folder of your Matlab. Eg. /HamadElKahza/downloads in my case

## Features


*	Plots the time and frequency domain components of the music created; 
*	Plots the Spectogram and other components of the work (Attached is a project description for your reference)
*	Generates 3 different audio files :
    1. The original Mortal Kombat theme song named MortalKombat.wav
    2. The song created adding the White Gaussian noise, the file is named MortalKombatWithNoise.wav
    3. The filtered song using a low pass filter, file named MortalKombatFiltered.wav
