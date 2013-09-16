/* iScissor.cpp */
/* Main file for implementing project 1.  See TODO statments below
 * (see also correlation.cpp and iScissor.h for additional TODOs) */

#include <assert.h>
#include <cmath>

#include "correlation.h"
#include "iScissor.h"

const double linkLengths[8] = { 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2 };

// two inlined routines that may help;

inline Node& NODE(Node* n, int i, int j, int width)
{
    return *(n + j * width + i);
}

inline unsigned char PIXEL(const unsigned char* p, int i, int j, int c, int width)
{
    return *(p + 3 * (j * width + i) + c);
}

/************************ TODO 1 ***************************
 *InitNodeBuf
 *	INPUT:
 *		img:	a RGB image of size imgWidth by imgHeight;
 *		nodes:	a allocated buffer of Nodes of the same size, one node corresponds to a pixel in img;
 *  OUTPUT:
 *      initializes the column, row, and linkCost fields of each node in the node buffer.
 */

void InitNodeBuf(Node* nodes, const unsigned char* img, int imgWidth, int imgHeight)
{
	double *newImgDoub = new double[imgWidth * imgHeight * 3];
	unsigned char *newImgChar = new unsigned char[imgWidth * imgHeight * 3];

	double blurKernel[] = {
		0.1111, 0.1111, 0.1111 ,
		0.1111, 0.1112, 0.1111 ,
		0.1111, 0.1111, 0.1111
	};

	image_filter(newImgDoub, img, NULL, imgWidth, imgHeight, blurKernel, 3, 3, 1, 0);
	
	for (int a = 0; a < imgWidth * imgHeight * 3; a++) {
		newImgChar[a] = ceil(newImgDoub[a] - 0.5);
	}

	double max = 0.0;

	for (int i = 0; i < imgWidth; i++) {
		for (int j = 0; j < imgHeight; j++) {
			NODE(nodes, i, j, imgWidth).column = i;
			NODE(nodes, i, j, imgWidth).row = j;
			for (int d = 0; d < 8; d++) {
				double pixel[] = {0.0, 0.0, 0.0};
				pixel_filter(pixel, i, j, newImgChar, imgWidth, imgHeight, kernels[d], 3, 3, 1, 0);
				double avg = (abs(pixel[0]) + abs(pixel[1]) + abs(pixel[2]))/3.0;
				if (d % 2 == 1) {
					avg = SQRT2 * avg;
				}
				if (avg > max) {
					max = avg;
				}
				NODE(nodes, i, j, imgWidth).linkCost[d] = avg;
			}
		}
	}

	for (int i = 0; i < imgWidth; i++) {
		for (int j = 0; j < imgHeight; j++) {
			for (int d = 0; d < 8; d++) {
				NODE(nodes, i, j, imgWidth).linkCost[d] = max - NODE(nodes, i, j, imgWidth).linkCost[d];
			
			}
		}
	}

}
/************************ END OF TODO 1 ***************************/

static int offsetToLinkIndex(int dx, int dy)
{
    int indices[9] = { 3, 2, 1, 4, -1, 0, 5, 6, 7 };
    int tmp_idx = (dy + 1) * 3 + (dx + 1);
    assert(tmp_idx >= 0 && tmp_idx < 9 && tmp_idx != 4);
    return indices[tmp_idx];
}

/************************ TODO 4 ***************************
 *LiveWireDP:
 *	INPUT:
 *		seedX, seedY:	seed position in nodes
 *		nodes:			node buffer of size width by height;
 *      width, height:  dimensions of the node buffer;
 *		selection:		if selection != NULL, search path only in the subset of nodes[j*width+i] if selection[j*width+i] = 1;
 *						otherwise, search in the whole set of nodes.
 *		numExpanded:		compute only the first numExpanded number of nodes; (for debugging)
 *	OUTPUT:
 *		computes the minimum path tree from the seed node, by assigning
 *		the prevNode field of each node to its predecessor along the minimum
 *		cost path from the seed to that node.
 */
void LiveWireDP(int seedX, int seedY, Node* nodes, int width, int height, const unsigned char* selection, int numExpanded)
{
    seedXX = seedX;
    seedYY = seedY;
    CTypedPtrHeap<Node> *pq = new CTypedPtrHeap<Node>();
    
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            NODE(nodes, i, j, width).state = INITIAL;
            NODE(nodes, i, j, width).prevNode = NULL;
            NODE(nodes, i, j, width).distance = 0;
        }
    }

    pq->Insert(&NODE(nodes, seedX, seedY, width));

    int curNumExpanded = 0;

    while (!pq->IsEmpty() && curNumExpanded < numExpanded) {
        Node& q = *pq->ExtractMin();
        q.state = EXPANDED;
        for (int i = 0; i < 8; ++i)
        {
            int offsetX, offsetY;
            q.nbrNodeOffset(offsetX, offsetY, i);

            int neighborCol = offsetX + q.column;
            int neighborRow = offsetY + q.row;

            if (neighborRow < 0 || neighborRow >= height || neighborCol < 0 || neighborCol >= width || (selection != NULL && selection[neighborRow*width+neighborCol] == 0)) continue;
            Node& r = NODE(nodes, neighborCol, neighborRow, width);

            if (r.state == INITIAL)
            {
                r.totalCost = q.totalCost + q.linkCost[i];
                r.state = ACTIVE;
                r.prevNode = &q;
                r.distance = q.distance + linkLengths[i];
                pq->Insert(&r);
            }
            else if (r.state == ACTIVE)
            {
                if ((q.totalCost + q.linkCost[i] < r.totalCost) || (q.totalCost + q.linkCost[i] == r.totalCost && r.distance > q.distance + linkLengths[i]))
                {
                    r.totalCost = q.totalCost + q.linkCost[i];
                    r.prevNode = &q;
                    r.distance = q.distance + linkLengths[i];
                }
            }
        }
        curNumExpanded ++;
    }
}
/************************ END OF TODO 4 ***************************/

/************************ TODO 5 ***************************
 *MinimumPath:
 *	INPUT:
 *		nodes:				a node buffer of size width by height;
 *		width, height:		dimensions of the node buffer;
 *		freePtX, freePtY:	an input node position;
 *	OUTPUT:
 *		insert a list of nodes along the minimum cost path from the seed node to the input node.
 *		Notice that the seed node in the buffer has a NULL predecessor.
 *		And you want to insert a *pointer* to the Node into path, e.g.,
 *		insert nodes+j*width+i (or &(nodes[j*width+i])) if you want to insert node at (i,j), instead of nodes[nodes+j*width+i]!!!
 *		after the procedure, the seed should be the head of path and the input code should be the tail
 */

void MinimumPath(CTypedPtrDblList <Node>* path, int freePtX, int freePtY, Node* nodes, int width, int height)
{
    Node* seed = &NODE(nodes, seedXX, seedYY, width);
    Node* curNode = &NODE(nodes, freePtX, freePtY, width);

    path->AddHead(curNode);
    while (curNode->prevNode) {
        path->AddHead(curNode->prevNode);
        curNode = curNode->prevNode;
    }
    path->AddHead(seed);

}
/************************ END OF TODO 5 ***************************/

/************************ An Extra Credit Item ***************************
 *SeedSnap:
 *	INPUT:
 *		img:				a RGB image buffer of size width by height;
 *		width, height:		dimensions of the image buffer;
 *		x,y:				an input seed position;
 *	OUTPUT:
 *		update the value of x,y to the closest edge based on local image information.
 */

void SeedSnap(int& x, int& y, unsigned char* img, int width, int height)
{
	double *dx = new double[width * height * 3];
	double *dy = new double[width * height * 3];
	
	unsigned char *selection = new unsigned char[width * height];

	int max = 0;
	int maxCol;
	int maxRow;

	double dxKernel[] = {
		-0.8889, 0.1111, 1.1111,
		-1.8889, 0.1112, 2.1111,
		-0.8889, 0.1111, 1.1111
	};

	double dyKernel[] = {
		1.1111, 2.1111, 1.1111,
		0.1111, 0.1112, 0.1111,
		-0.8889, -1.8889, -0.8889
	};

	for (int i = -10; i <= 10; i++) {
		for (int j = -10; j <= 10; j++) {
			selection[j*width+i] = 1;
		}
	}

	image_filter(dx, img, selection, width, height, dxKernel, 3, 3, 9, 0);
	image_filter(dy, img, selection, width, height, dyKernel, 3, 3, 9, 0);

	for (int i = -10; i <= 10; i++) {
		for (int j = -10; j <= 10; j++) {
			int col = i + x;
			int row = j + y;
			
			if (col >= 0 && col < width && row >= 0 && row < height) {
				int dxVal = 0;
				int dyVal = 0;
				int deriv = 0;
				double multiplier = 1 - (sqrt((i*i) + (j*j))/30);
				for (int c = 0; c < 3; c++) {
					dxVal += abs(dx[3*(row*width+col)+c]);
					dyVal += abs(dy[3*(row*width+col)+c]);
				}
				if (dxVal > dyVal) {
					deriv = dxVal;
				} else {
					deriv = dyVal;
				}
				deriv = deriv * multiplier;
				if (deriv > max) {
					max = deriv;
					maxCol = col;
					maxRow = row;
				}
			}
		}
	}

	printf("You clicked on %u, %u, seed snapped to max value of %u found at %u, %u\n", x, y, max, maxCol, maxRow);

	x = maxCol;
	y = maxRow;
}

//generate a cost graph from original image and node buffer with all the link costs;
void MakeCostGraph(unsigned char* costGraph, const Node* nodes, const unsigned char* img, int imgWidth, int imgHeight)
{
    int graphWidth = imgWidth * 3;
    int graphHeight = imgHeight * 3;
    int dgX = 3;
    int dgY = 3 * graphWidth;

    int i, j;
    for (j = 0; j < imgHeight; j++) {
        for (i = 0; i < imgWidth; i++) {
            int nodeIndex = j * imgWidth + i;
            int imgIndex = 3 * nodeIndex;
            int costIndex = 3 * ((3 * j + 1) * graphWidth + (3 * i + 1));

            const Node* node = nodes + nodeIndex;
            const unsigned char* pxl = img + imgIndex;
            unsigned char* cst = costGraph + costIndex;

            cst[0] = pxl[0];
            cst[1] = pxl[1];
            cst[2] = pxl[2];

            //r,g,b channels are grad info in seperate channels;
            DigitizeCost(cst	   + dgX, node->linkCost[0]);
            DigitizeCost(cst - dgY + dgX, node->linkCost[1]);
            DigitizeCost(cst - dgY      , node->linkCost[2]);
            DigitizeCost(cst - dgY - dgX, node->linkCost[3]);
            DigitizeCost(cst	   - dgX, node->linkCost[4]);
            DigitizeCost(cst + dgY - dgX, node->linkCost[5]);
            DigitizeCost(cst + dgY	   ,  node->linkCost[6]);
            DigitizeCost(cst + dgY + dgX, node->linkCost[7]);
        }
    }
}

