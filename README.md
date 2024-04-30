# FrameExtractor

Those small python macros help to extract info from a gwf file and store it into a ROOT file (for the moment)

## STEP 1: Data recovery

First of all you need to recover some data recorded from interferomaters. Below are the data location at **CCIN2P3**:

- **Standard strain**: /sps/virgo/USERS/mbta/O4/preprocessing/chunk01/mbtaD/prod
- **Standard strain + injections**: /sps/virgo/USERS/mbta/O4/preprocessing/chunk01/mbtaI/prod

In this case the data are stored in **gwf** files, grouped by paquets of 100000s per file. The file with injections normally contains all the necessary info.

## STEP 2: Get the info out of gwf file

One will process the data with the macro **gwf_extractor.py**. Suppose that you want to extract 100 seconds at a frequency of 1024 Hz of the file **input.gwf**, the command will be:

**python gwf_extractor.py —input input.gwf —output out.root -d 100 -fe 1024**

This will produce a **ROOT** file containing 2 trees: 

- **strain** will contain the non-gated data from Hanford and Livingston, both with and without soft injections.
- **injections** will contain the global MC truth of the injections.

## STEP 3: Whitening and PSDs production

A second macro will start from this file to produce the PSDs and the whitened data:

**python whiten.py -i out.root -o out_whiten.root**

The final file contains the following trees:

- **strain** will contain the non-gated data from Hanford and Livingston, both with and without soft injections, whitened.
- **injections** will contain the global MC truth of the injections.
- **PSDs**: the PSDs computed over the data, one every 50s by default
