#include "correlation.h"

/************************ TODO 2 **************************/
/*
 *	INPUT:
 *		origImg:		the original image,
 *		imgWidth:		the width of the image
 *		imgHeight:		the height of the image
 *						the image is arranged such that
 *						origImg[3*(row*imgWidth+column)+0],
 *						origImg[3*(row*imgWidth+column)+1],
 *						origImg[3*(row*imgWidth+column)+2]
 *						are R, G, B values for pixel at (column, row).
 *
 *      kernel:			the 2D filter kernel,
 *		knlWidth:		the width of the kernel
 *		knlHeight:		the height of the kernel
 *
 *		scale, offset:  after correlating the kernel with the origImg,
 *						each pixel should be divided by scale and then added by offset
 *
 *		selection:      a byte array of the same size as the image,
 *						indicating where in the original image should be filtered, e.g.,
 *						selection[k] == 1 ==> pixel k should be filtered
 *                      selection[k] == 0 ==> pixel k should NOT be filtered
 *                      a special case is selection is a NULL pointer, which means all the pixels should be filtered.
 *
 *  OUTPUT:
 *		rsltImg:		the filtered image of the same size as original image.
 *						it is a valid pointer ( allocated already ).
 */

void image_filter(double* rsltImg, const unsigned char* origImg, const unsigned char* selection,
                  int imgWidth, int imgHeight,
                  const double* kernel, int knlWidth, int knlHeight,
                  double scale, double offset)
{

printf("Filtering image\n");

for (int i = 0; i < imgWidth; i++) {
	for (int j = 0; j < imgHeight; j++) {
		double pixel[3] = {0.0, 0.0, 0.0};
		pixel_filter(pixel, i, j, origImg, imgWidth, imgHeight, kernel, knlWidth, knlHeight, scale, offset);
		for (int c = 0; c < 3; c++) {
			rsltImg[3*(j*imgWidth+i)+c] = pixel[c];
		}
	}
}

}

/************************ END OF TODO 2 **************************/


/************************ TODO 3 **************************/
/*
 *	INPUT:
 *      x:				a column index,
 *      y:				a row index,
 *		origImg:		the original image,
 *		imgWidth:		the width of the image
 *		imgHeight:		the height of the image
 *						the image is arranged such that
 *						origImg[3*(row*imgWidth+column)+0],
 *						origImg[3*(row*imgWidth+column)+1],
 *						origImg[3*(row*imgWidth+column)+2]
 *						are R, G, B values for pixel at (column, row).
 *
 *      kernel:			the 2D filter kernel,
 *		knlWidth:		the width of the kernel
 *		knlHeight:		the height of the kernel
 *
 *		scale, offset:  after correlating the kernel with the origImg,
 *						the result pixel should be divided by scale and then added by offset
 *
 *  OUTPUT:
 *		rsltPixel[0], rsltPixel[1], rsltPixel[2]:
 *						the filtered pixel R, G, B values at row y , column x;
 */

void pixel_filter(double rsltPixel[3], int x, int y, const unsigned char* origImg, int imgWidth, int imgHeight,
                  const double* kernel, int knlWidth, int knlHeight,
                  double scale, double offset)
{

//printf("Filtering pixel %u, %u\n", x, y);

for (int c = 0; c < 3; c++) {
	for (int i = 0; i < knlWidth; i++) {
		for (int j = 0; j < knlHeight; j++) {
			int col = i + x - ((knlHeight - 1)/2);
			int row = j + y - ((knlWidth - 1)/2);
			if (col >= imgWidth) {
				col = x;
			}
			if (row >= imgHeight) {
				row = y;
			}
			//printf("Fetching pixel %u, %u\n", col, row);
			rsltPixel[c] += kernel[j*knlWidth+i] * origImg[3*(row*imgWidth+col)+c];
		}
	}
	rsltPixel[c] /= scale;
	rsltPixel[c] += offset;
}

}

/************************ END OF TODO 3 **************************/

