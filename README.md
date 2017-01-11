# rtx2mnc
Converts a RT-structure from MIRADA to a MNC label file with the contours set to 1.
The program takes as input the volume file used to create the RTx struct, e.g. a PET image, in MNC format as well. The output label file will match the input volumes image dimensions.

## Usage:
<pre><code>
rtx2mnc < VOLUME.mnc > < RTx > < out_label.mnc >
      	
      	< VOLUME.mnc > is the file which the RTx was defined on.
      	             The out_label will have the dimensions and voxel size matching this file.
      	< RTx > is the RT struct in DICOM format.
      	< out_label.mnc > is the resulting MINC file with the contours in the RT file set to 1.
</code></pre>

## Requires/Dependencies:
 - MINC tools
 - DCMTK