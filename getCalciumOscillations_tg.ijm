setBatchMode(true);
run("Close All");

inputFolder = "";
outputFolder = "";

if (!(File.exists(outputFolder))) {
	File.makeDirectory(outputFolder);
}

// this is an initial radius - the contrast will be calculated between 0.5r and 1.5r.
r = 10;
// this is for Gaussian blurring of the image
rGaussian = 1;

processFolder(inputFolder, outputFolder);

function processFolder(inputFolder, outputFolder) {
	list = getFileList(inputFolder);
	for (i=0; i<list.length; i++) {
		processImage(inputFolder, list[i], outputFolder);
	}
}

function processImage(inputFolder, imageFilename, outputFolder) {
	
	if (endsWith(imageFilename, "tif")){
		
		temp = split(imageFilename,".");
		filenamePrefix = temp[0];
		print(filenamePrefix);
		
		listOutputFolder = outputFolder + filenamePrefix +"/";
		if (!(File.exists(listOutputFolder))) {
			File.makeDirectory(listOutputFolder);
		}
		
		open( inputFolder+imageFilename );
		win0 = getTitle();

		// get the dimensions of the stack
		Stack.getDimensions(width, height, channels, slices, frames);
		// make a new image with the same dimensions to record the segmented dots
		dotImage = filenamePrefix + "-dots.tif";
		newImage(dotImage, "8-bit black", width, height, 1);

		totalSliceNumber = nSlices;
		print(totalSliceNumber);
		selectWindow(win0);
		run("Z Project...", "projection=[Average Intensity]");
		rename("Average intensity projection");
		winAverage = getTitle();
		selectWindow(win0);
		run("Z Project...", "projection=[Max Intensity]");
		rename("Max intensity projection");
		winMIP = getTitle();
		imageCalculator("Subtract create", winMIP, winAverage);
		rename("Diff Max-Average");
		winDiff = getTitle();
		run("Gaussian Blur...", "sigma="+rGaussian);
		getStatistics(area, mean, min, max, std);

		selectWindow(winDiff);
		run("Duplicate...", " ");
		rRollingBall = floor(3*r) + 1;
		run("Subtract Background...", "rolling="+rRollingBall);
		rename("BG subtracted");
		winBGsub = getTitle();

		// may need to reduce the noise level.
		noiseLevel = 2.0*std;
		run("Clear Results");
		selectWindow(winDiff);
		run("Find Maxima...", "noise="+noiseLevel+" output=List exclude");

		maximaListFilename = outputFolder + filenamePrefix + "-maximaList.txt";
		saveAs("Results", maximaListFilename);
		selectWindow("Results");
		run("Close");//this is how to close non-image windows
		
		filestring=File.openAsString(maximaListFilename); 
		rows=split(filestring, "\n");
		x=newArray(rows.length);
		y=newArray(rows.length);
		
		for(i=1; i<rows.length; i++){
			
			columns=split(rows[i],"\t");
			x[i-1]=parseInt(columns[1]);
			y[i-1]=parseInt(columns[2]);
			
			selectWindow(winBGsub);
			results = estimateRadius(x[i-1], y[i-1], r);
			rMaxC = results[0];
			maxC = results[1];
			print(maxC);
			
			// contrast definition is (meanIn - meanOut) / (meanIn + meanOut)
			// 0.75 is equivalent to a 7-fold difference in to out
			// 0.3333 is equivalent to a 2-fold difference in to out
			contrastThreshold = 0.75;
			//contrastThreshold = 0.0;
			if ( maxC < contrastThreshold ){
				continue;
			}
			else {
				selectWindow(dotImage);
				setColor(255);
				fillOval(x[i-1]-rMaxC+1, y[i-1]-rMaxC+1, rMaxC*2-1, rMaxC*2-1);

				outFile = File.open(listOutputFolder + filenamePrefix + "-totalIntensitySeries-" + i + ".txt");
				print(outFile, "frame number" + "\t" + "total intensity");
				selectWindow(win0);
				for(j=1; j<slices+1; j++){
					setSlice(j);
					totalIntensity = getTotalIntensity(x[i-1], y[i-1], rMaxC);
					//print(totalIntensity);
					print(outFile, j + "\t" + totalIntensity);
				}
				File.close(outFile);
			}
		}
		selectWindow(dotImage);
		save(outputFolder + dotImage);
		run("Close All");
	}
	else {
		print(imageFilename + "is not a tif file!");
	}
	//run("Close All");
}

function getMeanIntensity(x, y, r) {
	// Calculates the mean intensity of an oval centered at (x,y) with radius r
	//makeOval(x-r, y-r, r*2+1, r*2+1);// makes the center to be exactly the specified pixel.
	makeOval(x-r+1, y-r+1, r*2-1, r*2-1);// makes the center to be exactly the specified pixel.
	getStatistics(areaIn, meanIntensityIn);
	run("Select None");
	return meanIntensityIn;
}

function getTotalIntensity(x, y, r) {
	// Calculates the total intensity of an oval centered at (x,y) with radius r,
	// where the local background is estimated in a 2-pixel width ring.
	//makeOval(x-r, y-r, r*2+1, r*2+1);// makes the center to be exactly the specified pixel.
	makeOval(x-r+1, y-r+1, r*2-1, r*2-1);// makes the center to be exactly the specified pixel.
	getStatistics(areaIn, meanIntensityIn);
	run("Enlarge...", "enlarge=2");// a 2-pixel width ring around the circle is considered as local background
	getStatistics(areaOut, meanIntensityOut);
	run("Select None");
	meanIntensityBg = (meanIntensityOut * areaOut - meanIntensityIn * areaIn) / (areaOut - areaIn);
	totalIntensity = (meanIntensityIn - meanIntensityBg) * areaIn;
	return totalIntensity;
}

function getContrastV0(x, y, r) {
	// The contrast here is defined as (I1 - I2) / (I1 + I2).
	// Therefore, a contrast of 0.33 corresponds to 2-fold difference,
	// whereas a contrast of 0.5 correspoinds to 3-fold difference.
	makeOval(x-r, y-r, r*2+1, r*2+1);// makes the center to be exactly the specified pixel.
	//makeOval(x-r+1, y-r+1, r*2-1, r*2-1);// makes the center to be exactly the specified pixel.
	getStatistics(areaIn, meanIntensityIn);
	run("Enlarge...", "enlarge=2");// a 2-pixel width ring around the circle is considered as local background
	getStatistics(areaOut, meanIntensityOut);
	run("Select None");
	meanIntensityBg = (meanIntensityOut * areaOut - meanIntensityIn * areaIn) / (areaOut - areaIn);
	contrast = (meanIntensityIn - meanIntensityBg) / (meanIntensityIn + meanIntensityBg);
	return contrast;
}

function getContrast(x, y, r) {
	// The contrast here is defined as (I1 - I2) / (I1 + I2).
	// I1 is the mean intensity of a 1-pixel width inner ring at r
	// I2 is the mean intensity of a 1-pixel width outer ring at r
	// Therefore, a contrast of 0.33 corresponds to 2-fold difference,
	// whereas a contrast of 0.5 correspoinds to 3-fold difference.
	makeOval(x-r, y-r, r*2+1, r*2+1);// makes the center to be exactly the specified pixel.
	getStatistics(area0, mean0);
	//makeOval(x-r+1, y-r+1, r*2-1, r*2-1);// makes the center to be exactly the specified pixel.
	run("Enlarge...", "enlarge=-1");// a 1-pixel width inner ring
	getStatistics(areaIn, meanIn);
	run("Enlarge...", "enlarge=2");// a 1-pixel width outer ring
	getStatistics(areaOut, meanOut);
	run("Select None");
	I1 = ( mean0 * area0 - meanIn * areaIn ) / (area0 - areaIn);
	I2 = ( meanOut * areaOut - mean0 * area0 ) / (areaOut - area0);
	contrast = ( I1 - I2 ) / ( I1 + I2 );
	return contrast;
}

function estimateRadius(x, y, r) {
	// Loop through every radius ranging from 0.5*r to 2*r, get the contrast.
	// The estimated radius is defined as the radius maximizing the contrast.
	rMin = floor(r/2);
	//rMin = 2;
	rMax = floor(1.5*r);
	N = rMax - 1;// the total number of radius points to search for
	contrastList = newArray(N);
	for(r=rMin; r < rMax+1; r++){
		contrast = getContrast(x, y, r);
		contrastList[r-rMin] = contrast;
	}
	Array.getStatistics(contrastList, minC, maxC, meanC, stdDevC);
	ranks = Array.rankPositions(contrastList);
	rMaxC = rMin + ranks[ranks.length-1];
	results = newArray(rMaxC, maxC);
	return results;
}