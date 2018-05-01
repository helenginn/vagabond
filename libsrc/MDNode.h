#ifndef __vagabond__MDNode__
#define __vagabond__MDNode__

/**
 * \class MDNode
 * \brief Accumulates information in a drill-down tree.
 **/ 

#include <cstdarg>

class CSV;

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
	
	/** Add to CSV object */
	void addToCSV(CSV *csv);
	
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
private:
	void addToNode(double *dimvals, double value);
	void addChildrenToCSV(CSV *csv);
	int _dims;
	
	double *_mins;
	double *_maxes;
	
	MDNode **_nodes;
	
	double _value;
	double _weight;
};


#endif
