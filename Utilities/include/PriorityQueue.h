
#ifndef __PRIORITYQUEUE_H_
#define __PRIORITYQUEUE_H_

#include <vector>
#include <assert.h>

/**
 * Fix sized maximum priority queue 
 *
 * @author Richard Keiser
 * @version 2.0
 */
template<class Index, class Weight>
class MaxPriorityQueue {

	
public:
	typedef struct el{
		Index	index;
		Weight weight;
	} Element;

	/** 
	 * set maximum size of queue
	 *
	 * @param size
	 *		  size of the queue
	 */
	void setSize(int size) {
		m_queue.resize(size+1);
		m_nOfElements = size;
		init();
	}

	/**
	 * initialize queue with zero elements
	 */
	inline void init() {
		m_current = 0;
	}

	/** 
	 * returns true if queue is empty
	 */
	inline bool isEmpty() {
		return m_current == 0;
	}

	/** 
	 * returns true if queue is full
	 */
	inline bool isFull() {
		return m_current==(m_nOfElements);
	}

	/** 
	 * returns number of elements in queue
	 */
	inline int getNofElements() {
		return m_current;
	}

	/**
	 * inserts a new element in O(log n)
	 * An element consists of an index and a weight
	 */
	inline void insert(Index index, Weight weight) {
		if (isFull()) {
			// if full replace max
			m_queue[1].index = index;
			m_queue[1].weight = weight;
			restore(1, m_nOfElements);
		} 
		else {
			m_current++;
			m_queue[m_current].index = index;
			m_queue[m_current].weight = weight;
			register int i=m_current;
			while(i>1 && (m_queue[i].weight > m_queue[i>>1].weight)) {
				swapElements(i, i>>1);
				i >>= 1;
			}
		}
	}

	/**
	 * gets the element with maximal weight in O(1)
	 */
	inline Element getMax() {
		assert(!isEmpty());
		return m_queue[1];
	}

	/**
	 * gets the index of the element with maximal weight in O(1)
	 */
	inline Index getMaxIndex() {
		assert(!isEmpty());
		return m_queue[1].index;
	}

	/** 
	 * gets the maximal weight in O(1)
	 */
	inline Weight getMaxWeight() {
		assert(!isEmpty());
		return m_queue[1].weight;
	}

	/**
	 * deletes the element with maximal weight in O(log n)
	 */
	inline void removeMax() {
		assert(!isEmpty());
		m_queue[1] = m_queue[m_current];
		restore(1, m_current);
		m_current--;
	}

	/**
	 * gets all elements in the queue
	 */
	inline Element* getQueue() {
		return &(*(m_queue.begin() + 1));
	}


protected:
	inline void restore(register int L, register int R) {
		int i, j;
		i = L;
		while (i <= (R>>1)) {
			if( 2*i < R && m_queue[2*i+1].weight > m_queue[2*i].weight) {
				j = 2*i+1;
			}
			else {
				j = 2*i;
			}
			if (m_queue[j].weight > m_queue[i].weight) {
				swapElements(j, i);
				swap(&i, &j);

			}
			else {
				i = R;
			}
		}
	}

private:
    std::vector<Element> m_queue;
	int m_nOfElements;
	int m_current;

	inline void swapElements(const int i, const int j) {
		Element tmp;
		tmp = m_queue[i];
		m_queue[i] = m_queue[j];
		m_queue[j] = tmp;
	}

	inline void swap(int *i, int *j) {
		int tmp = *i;
		*i = *j;
		*j = tmp;
	}
			
};

/**
 * Fix sized minimum priority queue 
 *
 * @author Richard Keiser
 * @version 2.0
 */
template<class Index, class Weight>
class MinPriorityQueue {

	
public:
	typedef struct el{
		Index	index;
		Weight	weight;
	} Element;

	/** 
	 * set maximum size of queue
	 *
	 * @param size
	 *		  size of the queue
	 */
	void setSize(int size) {
		m_queue.resize(size+1);
		m_nOfElements = size;
		init();
	}

	/**
	 * initialize queue with zero elements
	 */
	inline void init() {
		m_current = 0;
	}

	/** 
	 * returns true if queue is empty
	 */
	inline bool isEmpty() {
		return m_current == 0;
	}

	/** 
	 * returns true if queue is full
	 */
	inline bool isFull() {
		return m_current==(m_nOfElements);
	}

	/** 
	 * returns number of elements in queue
	 */
	inline int getNofElements() {
		return m_current;
	}

	/**
	 * inserts a new element in O(log n)
	 * An element consists of an index and a weight
	 */
	inline void insert(Index index, Weight weight) {
		if (isFull()) {
			// if full replace max
			m_queue[1].index = index;
			m_queue[1].weight = weight;
			restore(1, m_nOfElements);
		} 
		else {
			m_current++;
			m_queue[m_current].index = index;
			m_queue[m_current].weight = weight;
			register int i=m_current;
			while(i>1 && (m_queue[i].weight < m_queue[i>>1].weight)) {
				swapElements(i, i>>1);
				i >>= 1;
			}
		}
	}

	/**
	 * gets the element with minimal weight in O(1)
	 */
	inline Element getMin() {
		assert(!isEmpty());
		return m_queue[1];
	}

	/**
	 * gets the index of the element with minimum weight in O(1)
	 */
	inline Index getMinIndex() {
		assert(!isEmpty());
		return m_queue[1].index;
	}

	/** 
	 * gets the minimal weight in O(1)
	 */
	inline Weight getMinWeight() {
		assert(!isEmpty());
		return m_queue[1].weight;
	}

	/**
	 * deletes the element with minimal weight in O(log n)
	 */
	 inline void removeMin() {
		assert(!isEmpty());
		m_queue[1] = m_queue[m_current];
		restore(1, m_current);
		m_current--;
	}

	 /**
	  * gets the element in the queue
	  */
	inline Element* getQueue() {
		return &(*(m_queue.begin() + 1));
	}


protected:
	inline void restore(register int L, register int R) {
		int i, j;
		i = L;
		while (i <= (R>>1)) {
			if( 2*i < R && m_queue[2*i+1].weight < m_queue[2*i].weight) {
				j = 2*i+1;
			}
			else {
				j = 2*i;
			}
			if (m_queue[j].weight < m_queue[i].weight) {
				swapElements(j, i);
				swap(&i, &j);

			}
			else {
				i = R;
			}
		}
	}

private:
    std::vector<Element> m_queue;
	int m_nOfElements;
	int m_current;

	inline void swapElements(const int i, const int j) {
		Element tmp;
		tmp = m_queue[i];
		m_queue[i] = m_queue[j];
		m_queue[j] = tmp;
	}

	inline void swap(int *i, int *j) {
		int tmp = *i;
		*i = *j;
		*j = tmp;
	}
			
};

#endif

