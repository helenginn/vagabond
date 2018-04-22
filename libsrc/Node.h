#ifndef __vagabond__Node__
#define __vagabond__Node__

typedef struct Node
{
	double xMin;
	double xMax;
	double yMin;
	double yMax;
	Node *nextNodes[4];
	double value;
	int count;
} Node;

Node *malloc_node();
void free_node(Node *node);

void prepare_node(Node *node, int divisions, double xMin, double yMin,
                  double xMax, double yMax);
void add_to_node(Node *node, double x, double y, double z);

#endif
