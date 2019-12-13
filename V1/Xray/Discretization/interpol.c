/*****************************************************************************
 *	Date: January 5, 2009
 *----------------------------------------------------------------------------
 *	This C program is based on the following paper:
 *		P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
 *		IEEE Transactions on Medical Imaging,
 *		vol. 19, no. 7, pp. 739-758, July 2000.
 *----------------------------------------------------------------------------
 *	Philippe Thevenaz
 *	EPFL/STI/IMT/LIB/BM.4.137
 *	Station 17
 *	CH-1015 Lausanne VD
 *----------------------------------------------------------------------------
 *	phone (CET):	+41(21)693.51.61
 *	fax:			+41(21)693.37.01
 *	RFC-822:		philippe.thevenaz@epfl.ch
 *	X-400:			/C=ch/A=400net/P=switch/O=epfl/S=thevenaz/G=philippe/
 *	URL:			http://bigwww.epfl.ch/
 *----------------------------------------------------------------------------
 *	This file is best viewed with 4-space tabs (the bars below should be aligned)
 *	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|
 *  |...|...|...|...|...|...|...|...|...|...|...|...|...|...|...|...|...|...|
 ****************************************************************************/

/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#include "mex.h"

/*****************************************************************************
 *	Interpolation function
 ****************************************************************************/

extern double	InterpolatedValue
				(
					double	*Bcoeff,	/* input B-spline array of coefficients */
					mwSize	Width,		/* width of the image */
					mwSize	Height,		/* height of the image */
					double	x,			/* x coordinate where to interpolate */
					double	y,			/* y coordinate where to interpolate */
					mwSize	SplineDegree/* degree of the spline model */
				)
{ /* begin InterpolatedValue */

	double	*p;
	double	xWeight[10], yWeight[10];
	double	interpolated;
	double	w, w2, w4, t, t0, t1;
	mwSize	xIndex[10], yIndex[10];
	mwSize	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L;
	mwSize	i, j, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (mwSize)floor(x) - SplineDegree / 2L;
		j = (mwSize)floor(y) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	else {
		i = (mwSize)floor(x + 0.5) - SplineDegree / 2L;
		j = (mwSize)floor(y + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 0L:
			/* x */
			xWeight[0] = 1.0;
			/* y */
			yWeight[0] = 1.0;
			break;
		case 1L:
			/* x */
			w = x - (double)xIndex[0];
			xWeight[1] = w;
			xWeight[0] = 1.0 - xWeight[1];
			/* y */
			w = y - (double)yIndex[0];
			yWeight[1] = w;
			yWeight[0] = 1.0 - yWeight[1] ;
			break;
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			break;
		case 6L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= xWeight[0] / 720.0;
			xWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			xWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w
				* (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			xWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			xWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w
				* (1.0 + w)))))) / 48.0;
			xWeight[6] = 1.0 / 2.0 + w;
			xWeight[6] *= xWeight[6] * xWeight[6];
			xWeight[6] *= xWeight[6] / 720.0;
			xWeight[5] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[6];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= yWeight[0] / 720.0;
			yWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			yWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w
				* (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			yWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			yWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w
				* (1.0 + w)))))) / 48.0;
			yWeight[6] = 1.0 / 2.0 + w;
			yWeight[6] *= yWeight[6] * yWeight[6];
			yWeight[6] *= yWeight[6] / 720.0;
			yWeight[5] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] - yWeight[3]
				- yWeight[4] - yWeight[6];
			break;
		case 7L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			xWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				* (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
			xWeight[2] = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			xWeight[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2
				* (-4.0 + w)))) / 144.0;
			xWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
			xWeight[5] = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
			xWeight[7] = w2;
			xWeight[7] *= xWeight[7] * xWeight[7];
			xWeight[7] *= w / 5040.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[5] - xWeight[7];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			yWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				* (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
			yWeight[2] = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			yWeight[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2
				* (-4.0 + w)))) / 144.0;
			yWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
			yWeight[5] = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
			yWeight[7] = w2;
			yWeight[7] *= yWeight[7] * yWeight[7];
			yWeight[7] *= w / 5040.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] - yWeight[3]
				- yWeight[4] - yWeight[5] - yWeight[7];
			break;
		case 8L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] / 40320.0;
			w2 = w * w;
			xWeight[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))
				* (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w
				* (-3.0 + w)))) / 5040.0;
			xWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
			xWeight[3] = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
			xWeight[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			xWeight[5] = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w
				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
			xWeight[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))
				* (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w
				* (3.0 + w)))) / 5040.0;
			xWeight[8] = 1.0 / 2.0 + w;
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8] / 40320.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[5] - xWeight[7] - xWeight[8];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] / 40320.0;
			w2 = w * w;
			yWeight[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))
				* (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w
				* (-3.0 + w)))) / 5040.0;
			yWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
			yWeight[3] = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
			yWeight[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			yWeight[5] = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w
				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
			yWeight[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))
				* (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w
				* (3.0 + w)))) / 5040.0;
			yWeight[8] = 1.0 / 2.0 + w;
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8] / 40320.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] - yWeight[3]
				- yWeight[4] - yWeight[5] - yWeight[7] - yWeight[8];
			break;
		case 9L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * (1.0 - w) / 362880.0;
			xWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			xWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			xWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			xWeight[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			xWeight[5] = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			xWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			xWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			xWeight[9] = w2 * w2;
			xWeight[9] *= xWeight[9] * w / 362880.0;
			xWeight[8] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[5] - xWeight[6] - xWeight[7] - xWeight[9];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * (1.0 - w) / 362880.0;
			yWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			yWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			yWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			yWeight[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			yWeight[5] = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			yWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			yWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			yWeight[9] = w2 * w2;
			yWeight[9] *= yWeight[9] * w / 362880.0;
			yWeight[8] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] - yWeight[3]
				- yWeight[4] - yWeight[5] - yWeight[6] - yWeight[7] - yWeight[9];
			break;
		default:
			mexErrMsgTxt("Invalid spline degree");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	for (k = 0L; k <= SplineDegree; k++) {
		yIndex[k] = (Width == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Width2 * ((-yIndex[k]) / Width2))
			: (yIndex[k] - Width2 * (yIndex[k] / Width2)));
		if (Width <= yIndex[k]) {
			yIndex[k] = Width2 - yIndex[k];
		}
			
		xIndex[k] = (Height == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Height2 * ((-xIndex[k]) / Height2))
			: (xIndex[k] - Height2 * (xIndex[k] / Height2)));
		if (Height <= xIndex[k]) {
			xIndex[k] = Height2 - xIndex[k];
		}
	}

	/* perform interpolation */
	interpolated = 0.0;
	for (j = 0L; j <= SplineDegree; j++) {
		// p = Bcoeff + (ptrdiff_t)(yIndex[j] * Width); // row major
		p = Bcoeff + (ptrdiff_t)(yIndex[j] * Height); // column major
		w = 0.0;
		for (i = 0L; i <= SplineDegree; i++) {
			w += xWeight[i] * p[xIndex[i]];
		}
		interpolated += yWeight[j] * w;
	}

	return(interpolated);
} /* end InterpolatedValue */


/*****************************************************************************
 *	MEX entry point function
 ****************************************************************************/



void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	if ( !( nrhs == 3 || nrhs == 4) ) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
				"3 or 4 inputs required.");
	}
	
	double* Bcoeff = mxGetPr( prhs[0] ); // column major (from MATLAB)
	
	mwSize Width, Height;
	mwSize SplineDegree;
	
	if (nrhs == 3) { // 1D
		
		Height = fmax( mxGetN(prhs[0]), mxGetM(prhs[0]) );
		Width = 1;
		
		mwSize M = mxGetM(prhs[1]);
		mwSize N = mxGetN(prhs[1]);
		
		double *x = mxGetPr( prhs[1] );
		double y = 0;
		
		SplineDegree = mxGetScalar( prhs[2] );
		
		plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
		double* result = mxGetPr( plhs[0] );
		
		for (mwSize i = 0; i < M*N; ++i) {
			result[i] = InterpolatedValue(Bcoeff, Width, Height, x[i], y, SplineDegree);
		}
		
		
	} else if (nrhs == 4) { // 2D

		Width = mxGetN(prhs[0]);
		Height = mxGetM(prhs[0]);
		
		if (mxGetM(prhs[1]) != mxGetM(prhs[2]) || mxGetN(prhs[1]) != mxGetN(prhs[2])) {
			mexErrMsgTxt("x and y dims must match");
		}
		
		mwSize M = mxGetM(prhs[1]);
		mwSize N = mxGetN(prhs[1]);
		double *x = mxGetPr( prhs[1] );
		double *y = mxGetPr( prhs[2] );
		
		SplineDegree = mxGetScalar( prhs[3] ); 
		
		plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
		double* result = mxGetPr( plhs[0] );
		
		for (mwSize i = 0; i < M*N; ++i) {
			result[i] = InterpolatedValue(Bcoeff, Width, Height, x[i], y[i], SplineDegree);
		}
	}
	
	
	/*mexPrintf("x %g\n", x);
	mexPrintf("y %g\n", y);
	mexPrintf("Width %d\n", Width);
	mexPrintf("Height %d\n", Height);
	mexPrintf("SplineDegree %d\n", SplineDegree);*/

	


}