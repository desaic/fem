
#include "KdTree.h"

#define SWAP_POINTS(a,b) \
			KdTreePoint tmp = points[a];\
		    points[a] = points[b];\
		    points[b] = tmp;


// ******************
// global definitions
// ******************
double*						g_queryOffsets = new double[3];
double					g_queryPosition[3];

KdTree::KdTree(const double *positions, int idim, const unsigned int nOfPositions, const unsigned int maxBucketSize) {
  m_dim = idim;
	m_bucketSize			= maxBucketSize;
	m_positions				= positions;
	m_nOfPositions			= nOfPositions;
	m_points				= new KdTreePoint[nOfPositions];
	m_nOfFoundNeighbours	= 0;
	m_nOfNeighbours			= 0;
	m_queryPriorityQueue	= new PQueue();
	for (unsigned int i=0; i<nOfPositions; i++) {
    for (int icoord=0; icoord<m_dim; icoord++)
    {
		  m_points[i].pos[icoord] = positions[m_dim*i+icoord];
    }
		m_points[i].index = i;
	}
	computeEnclosingBoundingBox(m_boundingBoxLowCorner, m_boundingBoxHighCorner);
	m_root = new KdNode();
	double maximum[3], minimum[3];
	getSpread(m_points, nOfPositions, maximum, minimum);
	createTree(*m_root, 0, nOfPositions, maximum, minimum);
	setNOfNeighbours(1);
}


KdTree::~KdTree() {
	delete m_root;
	delete[] m_points;
	delete m_queryPriorityQueue;
}

void KdTree::computeEnclosingBoundingBox(double *lowCorner, double *hiCorner) {
	const double *tmp;
	copyPosition(hiCorner, m_points[0].pos);
  copyPosition(lowCorner, m_points[0].pos);

	for (unsigned int i=1; i<m_nOfPositions; i++) {
		tmp = m_positions + m_dim*i;
		if (hiCorner[0] < tmp[0]) {
			hiCorner[0] = tmp[0];
		}
		else if (lowCorner[0] > tmp[0]) {
			lowCorner[0] = tmp[0];
		}
		if (hiCorner[1] < tmp[1]) {
			hiCorner[1] = tmp[1];
		}
		else if (lowCorner[1] > tmp[1]) {
			lowCorner[1] = tmp[1];
		}
    if (m_dim>2)
    {
      if (hiCorner[2] < tmp[2]) {
        hiCorner[2] = tmp[2];
      }
      else if (lowCorner[2] > tmp[2]) {
        lowCorner[2] = tmp[2];
      }
    }
  }		 
}

double KdTree::computeBoxDistance(const double *q, const double *lo, const double *hi) {
	register double dist = 0.0;
	register double t;

	if (q[0] < lo[0]) {
		t = lo[0] - q[0];
		dist = t*t;
	}
	else if (q[0] > hi[0]) {
		t = q[0] - hi[0];
		dist = t*t;
	}
	if (q[1] < lo[1]) {
		t = lo[1] - q[1];
		dist += t*t;
	}
	else if (q[1] > hi[1]) {
		t = q[1] - hi[1];
		dist += t*t;
	}
  if (m_dim>2)
  {
    if (q[2] < lo[2]) {
      t = lo[2] - q[2];
      dist += t*t;
    }
    else if (q[2] > hi[2]) {
      t = q[2] - hi[2];
      dist += t*t;
    }
  }

  return dist;
}

void KdTree::queryPosition(const double *position) {
	if (m_neighbours.size() == 0) {
		return;
	}
	g_queryOffsets[0] = g_queryOffsets[1] = g_queryOffsets[2] = 0.0;
	m_queryPriorityQueue->init();
	m_queryPriorityQueue->insert(-1, DBL_MAX);
	copyPosition(g_queryPosition, position);
	double dist = computeBoxDistance(position, m_boundingBoxLowCorner, m_boundingBoxHighCorner);
	m_root->queryNode(dist, m_queryPriorityQueue);
	
	if (m_queryPriorityQueue->getMax().index == -1) {
		m_queryPriorityQueue->removeMax();
	}

	m_nOfFoundNeighbours = m_queryPriorityQueue->getNofElements();

	for(int i=m_nOfFoundNeighbours-1; i>=0; i--) {
		m_neighbours[i] = m_queryPriorityQueue->getMax();
		m_queryPriorityQueue->removeMax();
	}
}

void KdTree::queryRange(const double *position, const double maxSqrDistance) {
	if (m_neighbours.size() == 0) {
		return;
	}
	g_queryOffsets[0] = g_queryOffsets[1] = g_queryOffsets[2] = 0.0;
	m_queryPriorityQueue->init();
	m_queryPriorityQueue->insert(-1, maxSqrDistance);
	copyPosition(g_queryPosition, position);

	double dist = computeBoxDistance(position, m_boundingBoxLowCorner, m_boundingBoxHighCorner);	
	m_root->queryNode(dist, m_queryPriorityQueue);

	if (m_queryPriorityQueue->getMax().index == -1) {
		m_queryPriorityQueue->removeMax();
	}

	m_nOfFoundNeighbours = m_queryPriorityQueue->getNofElements();

	for(int i=m_nOfFoundNeighbours-1; i>=0; i--) {
		m_neighbours[i] = m_queryPriorityQueue->getMax();
		m_queryPriorityQueue->removeMax();
	}
}

void KdTree::setNOfNeighbours (const unsigned int newNOfNeighbours) {
	if (newNOfNeighbours != m_nOfNeighbours) {
		m_nOfNeighbours = newNOfNeighbours;
		m_queryPriorityQueue->setSize(m_nOfNeighbours);
		m_nOfNeighbours = newNOfNeighbours;
		m_neighbours.resize(m_nOfNeighbours);
		m_nOfFoundNeighbours = 0;
	}
}


void KdTree::createTree(KdNode &node, int start, int end, double *maximum, double*minimum) {
	int	mid;

	int n = end-start;
	double diff[3];
  for (int icoord=0; icoord<m_dim; icoord++)
  {
    diff[icoord] = maximum[icoord] - minimum[icoord];
  }
	short dim;
	// get longest axe
	if (diff[0] > diff[1]) {
		if (m_dim < 3 || diff[0] > diff[2]) {
			dim = 0;	//x-axe is longest axe
		}
		else {
			dim = 2;	// z-axe is longest axe
		}
	}
	else {
		if (m_dim < 3 || diff[1] > diff[2]) {
			dim = 1;	// y-axe is longest axe
		}
		else {
			dim = 2;	// z-axe is longest axe
		}
	}
	
	node.m_dim = (unsigned char)dim;
	double bestCut = (maximum[dim]+minimum[dim])/2.0f;
	double min, max;
	getMinMax(m_points+start, n, dim, min, max);	// find min/max coordinates
	if (bestCut < min)		// slide to min or max as needed
		node.m_cutVal = min;
    else if (bestCut > max)
		node.m_cutVal = max;
    else
		node.m_cutVal = bestCut;

    int br1, br2;
	splitAtMid(m_points+start, n, dim, node.m_cutVal, br1, br2);	// permute points accordingly

	if (bestCut < min) mid = start+1;
    else if (bestCut > max) mid = end-1;
    else if (br1 > n/2.0) mid = start+br1;
    else if (br2 < n/2.0) mid = start+br2;
    else mid = start + (n>>1);

	TreeNode** childNodes = new TreeNode*[2];
	node.m_children = childNodes;
	if (mid-start <= m_bucketSize) {
		// new leaf
		KdLeaf* leaf = new KdLeaf();
		node.m_children[0] = leaf;
		leaf->m_points = (m_points+start);
		leaf->m_nOfElements = mid-start;
    leaf->m_dim = m_dim;
	}
	else {
		// new node
		KdNode* childNode = new KdNode();
		node.m_children[0] = childNode;
		double oldMax = maximum[dim];
		maximum[dim] = node.m_cutVal;
		createTree(*childNode, start, mid, maximum, minimum);
		maximum[dim] = oldMax;
	}
	
	if (end-mid <= m_bucketSize) {
		// new leaf
		KdLeaf* leaf = new KdLeaf();
		node.m_children[1] = leaf;
		leaf->m_points = (m_points+mid);
		leaf->m_nOfElements = end-mid;
    leaf->m_dim = m_dim;
	}
	else {
		// new node
		minimum[dim] = node.m_cutVal;
		KdNode* childNode = new KdNode();
		node.m_children[1] = childNode;
		createTree(*childNode, mid, end, maximum, minimum);
	}
}

void KdTree::getSpread(KdTreePoint* points, int nOfPoints, double *maximum, double *minimum) {
	double *pos;
  copyPosition(maximum, points->pos);
  copyPosition(minimum, points->pos);
	points++;
	for (int i = 1; i < nOfPoints; i++) {
		pos = points->pos;
		if (pos[0] < minimum[0]) {
			minimum[0] = pos[0];
		}
		if (pos[0] > maximum[0]) {
			maximum[0] = pos[0];
		}
		if (pos[1] < minimum[1]) {
			minimum[1] = pos[1];
		}
		if (pos[1] > maximum[1]) {
			maximum[1] = pos[1];
		}
    if (m_dim>2)
    {
      if (pos[2] < minimum[2]) {
        minimum[2] = pos[2];
      }
      if (pos[2] > maximum[2]) {
        maximum[2] = pos[2];
      }
    }
    points++;
  }
}

void KdTree::getMinMax(KdTreePoint *points, int nOfPoints, int dim, double &min, double &max) {
	min = points->pos[dim];
	max = points->pos[dim];
	points++;
	for (int i=1; i<nOfPoints; i++) {
		if (points->pos[dim] < min) {
			min = points->pos[dim];
		}
		else if (points->pos[dim] > max) {
			max = points->pos[dim];
		}
		points++;
	}
}


void KdTree::splitAtMid(KdTreePoint *points, int nOfPoints, int dim, double cutVal, int &br1, int &br2) {
    int l = 0;
    int r = nOfPoints-1;
    for(;;) {				// partition points[0..n-1] about the cut value
		while (l < nOfPoints && points[l].pos[dim] < cutVal) {
			l++;
		}
		while (r >= 0 && points[r].pos[dim] >= cutVal) {
			r--;
		}
		if (l > r) 
			break;
		SWAP_POINTS(l,r);
		l++; 
		r--;
    }
    br1 = l;			// now: points[0..br1-1] < cutVal <= points[br1..n-1]
    r = nOfPoints-1;
    for(;;) {				// partition points[br1..n-1] about cutVal
		while (l < nOfPoints && points[l].pos[dim] <= cutVal) { 
			l++;
		}
		while (r >= br1 && points[r].pos[dim] > cutVal) {
			r--;
		}
		if (l > r) 
			break;
		SWAP_POINTS(l,r);
		l++; 
		r--;
    }
    br2 = l;			// now: points[br1..br2-1] == cutVal < points[br2..n-1]
}

void KdNode::queryNode(double rd, PQueue* queryPriorityQueue) {
	register double old_off = g_queryOffsets[m_dim];
	register double new_off = g_queryPosition[m_dim] - m_cutVal;
	if (new_off < 0) {
		m_children[0]->queryNode(rd, queryPriorityQueue);
		rd = rd - SQR(old_off) + SQR(new_off);
		if (rd < queryPriorityQueue->getMaxWeight()) {
			g_queryOffsets[m_dim] = new_off;
			m_children[1]->queryNode(rd, queryPriorityQueue);
			g_queryOffsets[m_dim] = old_off;
		}
	}
	else {
		m_children[1]->queryNode(rd, queryPriorityQueue);
		rd = rd - SQR(old_off) + SQR(new_off);
		if (rd < queryPriorityQueue->getMaxWeight()) {
			g_queryOffsets[m_dim] = new_off;
			m_children[0]->queryNode(rd, queryPriorityQueue);
			g_queryOffsets[m_dim] = old_off;
		}
	}

}

void KdLeaf::queryNode(double rd, PQueue* queryPriorityQueue) {
	double sqrDist;
	//use pointer arithmetic to speed up the linear traversing
	KdTreePoint* point = m_points;
	for (register unsigned int i=0; i<m_nOfElements; i++) {
		sqrDist = 0;
    for (int icoord=0; icoord<m_dim; icoord++)
    {
      double diff = point->pos[icoord] - g_queryPosition[icoord];
      sqrDist += diff*diff;
    }
		if (sqrDist < queryPriorityQueue->getMaxWeight()) {
			queryPriorityQueue->insert(point->index, sqrDist);
		}
		point++;
	}		
}

