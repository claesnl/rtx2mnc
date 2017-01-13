#include <minc2.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <vector>
#include <math.h>
#include <map>
#include <ParseArgv.h>
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmdata/dcpxitem.h"

using namespace std;

/* Set variables for handling the input volumes */
mihandle_t	container_volume, label_volume;
midimhandle_t dimensions[3], *dimensions_new;
char 		*pname;
char 		**infiles, **outfiles;

/* Varies for the RT struct and MINC volume */
int 		i,j,slice;
misize_t start[3], count[3];
unsigned int  sizes[3];
double        *slab;
double        *slab_new;
double        world_location[3];
double        dvoxel_location[3];
unsigned long voxel_location[3];
map< int,vector< vector<double> > > contours;

int main(int argc, char **argv)
{
	
	/* Input handling */
	pname = argv[0];
	if(argc < 2)
	{
      	(void) fprintf(stderr, "\nUsage: %s <VOLUME.mnc> <RTx> <out_label.mnc>\n", pname);
      	(void) fprintf(stderr, "       <VOLUME.mnc>\tis the file which the RTx was defined on.\n");
      	(void) fprintf(stderr, "                   \tThe out_label will have the dimensions and voxel size matching this file.\n");
      	(void) fprintf(stderr, "       <RTx>\t\tis the RT struct in DICOM format.\n");
      	(void) fprintf(stderr, "       <out_label.mnc>\tis the resulting MINC file with the contours in the RT file set to 1.\n\n");
		exit(EXIT_FAILURE);
	}
	outfiles = &argv[argc-1];
	infiles = &argv[1];

	/* READ DICOM RT STRUCT 
	 * The file is formatted as (one) ROIContourSequence
	 *						 -> (Abitrary number of) ContourSequence
	 *						 -> -> ContourNumber, NumberOfContourPoints and ContourData
	 * Contour data is a double-string consisting of x\y\z\[..]\x\y\z coordinates.
	 *
	*/

	fprintf(stderr, "Loading RTx file: %s\n",infiles[1]);

	DcmFileFormat fileformat;
  	if (fileformat.loadFile(infiles[1]).good()){
  		// Load the DICOM-RT dataset
  		DcmDataset *dataset = fileformat.getDataset();
  		// Find the ROIContourSequence tag
  		DcmSequenceOfItems *pROIContourSequnce = NULL;
  		if(dataset->findAndGetSequence(DCM_ROIContourSequence, pROIContourSequnce, true).good()){
  			// Assuming only 1 - get it.
  			DcmDataset *ROIContourSequence = (DcmDataset *) pROIContourSequnce->getItem(0);
  			// Find the ContourSequence tags
  			DcmSequenceOfItems *pContourSequnce = NULL;
  			if(ROIContourSequence->findAndGetSequence(DCM_ContourSequence, pContourSequnce, true).good()){
  				// Get the number of contour sequences in the file
  				Uint16 numimages = pContourSequnce->card();

  				fprintf(stderr, "\tFound %d contour data entries in the RTx file.\n",numimages);

  				for( int iter=0; iter<numimages; iter++){

  					// Get each ContourSequence at a time
  					DcmItem *ContourSequence = pContourSequnce->getItem(iter);
  					
  					// For future reference, get the contour Number and number of points in sequence
  					const char *number = NULL;
  					ContourSequence->findAndGetString(DCM_ContourNumber, number, 4);
  					string s_number(number);
  					const char *npoints = NULL;
  					ContourSequence->findAndGetString(DCM_NumberOfContourPoints, npoints, 4);
  					string s_npoints(npoints);

  					// Get the contour data points
  					vector< vector<double> > xyz_coordinates;
					const char *cdata = NULL;
  					ContourSequence->findAndGetString(DCM_ContourData, cdata, 4);
  					string s_cdata(cdata);

  					// Extract each combination of x,y,z coordinates and insert them in map
  					size_t pos = 0;
  					string token;
  					string delimiter = "\\";
  					int xyz_counter = 0;
  					int coordinate_counter = 0;
  					while((pos = s_cdata.find(delimiter)) != string::npos){
  						token = s_cdata.substr(0,pos);

  						// If 0, we are at x again, add new coordinate point to vector
  						if(xyz_counter == 0){
  							vector<double> new_point(1,3);
  							xyz_coordinates.push_back(new_point);
  						}

  						// Add coordinate point (x, y or z) to vector
  						xyz_coordinates[coordinate_counter].push_back(atof(token.c_str()));

  						// If 3, we have added z, reset to x for next iteration
  						if(xyz_counter+1 == 3){
  							xyz_counter = 0;
  							coordinate_counter += 1;
  						} else {
  							xyz_counter += 1;
  						}

  						s_cdata.erase(0,pos+delimiter.length());
  					}
  					// Add final z that is left in the sequence
  					xyz_coordinates[coordinate_counter].push_back(atof(s_cdata.c_str()));
  					// Add all points from contour to map of contours
  					contours.insert( pair< int,vector< vector<double> > >(iter,xyz_coordinates));

  				}

  			}

  		}

  	}
	
	/* Load the MINC container file from which we will obtain 
	 * information such as voxel size and image dimensions.
	 *
	 * We can only work with minc 2.0 files */

  	fprintf(stderr, "Loading PET file: %s\n",infiles[0]);

	system("mkdir /tmp/minc_plugins/");
	std::string in1 = "mincconvert -clobber ";
	in1 += infiles[0];
	in1 += " -2 /tmp/minc_plugins/in.mnc";
	system(in1.c_str());
	

	if (miopen_volume("/tmp/minc_plugins/in.mnc", MI2_OPEN_RDWR, &container_volume) != MI_NOERROR)
	{
		fprintf(stderr, "Error opening input file: %s.\n",infiles[0]);
		exit(EXIT_FAILURE);
	}
	
	// allocate new dimensions 
	dimensions_new = (midimhandle_t *) malloc(sizeof(midimhandle_t) * 3);

	// get the dimension sizes
	miget_volume_dimensions(container_volume, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL, MI_DIMORDER_FILE, 3, dimensions);
	misize_t sizes_tmp[3];
	miget_dimension_sizes(dimensions, 3, sizes_tmp);
	for (i=0; i < 3; i++) {
		sizes[i] = (unsigned int) sizes_tmp[i];
		count[i] = (unsigned long) sizes[i];
		micopy_dimension(dimensions[i], &dimensions_new[i]);
	}

	fprintf(stderr, "\tDimensions of PET file: %dx%dx%d\n",sizes[0],sizes[1],sizes[2]);

	// allocate memory for the hyperslab - make it size of entire volume
	start[0] = start[1] = start[2] = 0;
	slab_new = (double *) malloc(sizeof(double) * sizes[0] * sizes[1] * sizes[2]);

	// Initialize new hyperslab for the label file with zeros everywhere
  	for (i=0; i < sizes[0] * sizes[1] * sizes[2]; i++) {
		slab_new[i] = 0.0;
	}
	
	/* Loop through the RT struct map, set new label to 1 for every pixel on the contours */
	map< int , vector < vector<double> > >::iterator it = contours.begin();
  	for(it=contours.begin(); it!= contours.end(); ++it){
  		vector< vector<double> > contour = it->second;
		for(int cord_i = 0; cord_i < contour.size(); ++cord_i){

			// Convert to voxel coordinate
			world_location[0] = -contour[cord_i][1]; // manual -1* (why?)
			world_location[1] = -contour[cord_i][2]; // manual -1* (why?)
			world_location[2] = contour[cord_i][3];
			miconvert_world_to_voxel(container_volume, world_location, dvoxel_location);
			for (int vox_i=0; vox_i<3; vox_i++)
				voxel_location[vox_i] = (unsigned long) dvoxel_location[vox_i];
			// Calculate the position in the hyperslab 
			int contour_voxel = voxel_location[0]*sizes[1]*sizes[2] + voxel_location[1]*sizes[2] + voxel_location[2];
			// Set new label value to 1
			slab_new[contour_voxel] = 1.0;
		}
	}

	fprintf(stderr, "Writing volume and data to file: %s\n",outfiles[0]);

	// create the new volume 
	if (micreate_volume("/tmp/minc_plugins/label_volume.mnc", 3, dimensions_new, MI_TYPE_UBYTE,
	MI_CLASS_REAL, NULL, &label_volume) != MI_NOERROR) {
		fprintf(stderr, "Error creating new volume\n");
		return(1);
	}
	
	// create the data for the new volume 
	if (micreate_volume_image(label_volume) != MI_NOERROR) {
		fprintf(stderr, "Error creating volume data\n");
		return(1);
	}
		
	// set valid and real range 
	miset_volume_valid_range(label_volume, 255, 0);
	miset_volume_range(label_volume, 32, 0);

	// write the modified hyperslab to the file 
	if (miset_real_value_hyperslab(label_volume, MI_TYPE_DOUBLE,
		start, count, slab_new) != MI_NOERROR) {
		fprintf(stderr, "Error setting hyperslab\n");
		return(1);
	}

	// closes the volume and makes sure all data is written to file 
	miclose_volume(container_volume);
	miclose_volume(label_volume);

	// free memory 
	free(dimensions_new);
	free(slab);
	free(slab_new);
	
	// Convert back to minc 1 volume and save to wanted location
	std::string s_outfiles = "mincconvert -clobber /tmp/minc_plugins/label_volume.mnc ";
	s_outfiles += outfiles[0];
	system(s_outfiles.c_str());
	
	system("rm -rf /tmp/minc_plugins/");

	fprintf(stderr, "Finished rtx2mnc\n");

	return(0);
}
