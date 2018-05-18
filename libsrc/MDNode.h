#ifndef __vagabond__MDNode__
#define __vagabond__MDNode__

/**
 * \class MDNode
 * \brief Accumulates information in a drill-down tree to create a lookup table
 * or output to a CSV file.
 **/ 

#include <cstdarg>

class CSV;
class Plucker;

class MDNode
{
public:
	MDNode(int dims);
	~MDNode();
	
	/** Prior to splitting the node, set the minimum and maximum of each
	* dimension */
	void setDimension(int dim, double min, double max);

	/** Make sure all dimensions have been set prior to this */
	void splitNode(int divisions, int remaining);
	
	/** Generate a Plucker object and store it locally for this node.
	* @param subtract value to subtract from all weights to calculate the
	* Plucker chances. */
	void makePlucker(int dim, double value, double subtract);
	
	/** Get the plucker object from MDNode.
	* \param give_up if 1, will set plucker pointer to zero without release
	* \return Plucker pointer. */ 
	Plucker *getPlucker(int give_up)
	{
		Plucker *pluck = _plucker;
		
		if (give_up)
		{
			_plucker = NULL;		
		}

		return pluck;
	}

	/** Add to CSV object - to be called from CSV */
	void addToCSV(CSV *csv);

	void addChildrenToCSV(CSV *csv);
	
	double getValue()
	{
		return _value;
	}

	/** Add a value to increment the corresponding child node.
	* \param value The number to increment the node by.
	* \param count Should be equal to the number of dimensions of the array.
	* \param ... Values of dimensions to specify appropriate node.
	*/
	void addToNode(double value, int count, ...)
	{
		va_list arguments;
		va_start(arguments, count);
		double *dimvals = new double[count];

		for (int i = 0; i < count; i++)
		{
			double dimval = (double)(va_arg(arguments, double));
			dimvals[i] = dimval;
		}

		va_end(arguments);
		
		addToNode(dimvals, value);

		delete [] dimvals;
	}

	size_t nodeCount()
	{
		if (_nodes == NULL) return 0;
		
		return pow(2, _dims);
	}
	
	MDNode *node(int i)
	{
		return _nodes[i];
	}
	
	double aveDimension(int dim)
	{
		return (_mins[dim] + _maxes[dim]) / 2;
	}
private:
	MDNode *findNode(double *dimvals);
	Plucker *_plucker;
	void addNodesToPlucker(Plucker *pluck, int dim,
	                       double value, double subtract);
	void addToNode(double *dimvals, double value);
	int _dims;
	
	double *_mins;
	double *_maxes;
	
	MDNode **_nodes;
	
	double _value;
	double _weight;
};


#endif
