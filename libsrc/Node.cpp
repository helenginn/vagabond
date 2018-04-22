#include "Node.h"
#include <climits>
#include <cstring>
#include <memory>

Node *malloc_node()
{
	Node *node = (Node *)calloc(sizeof(Node), 1);
	return node;	
}

void free_node(Node *node)
{
	if (node == NULL)
	{
		return;
	}
	
	for (int i = 0; i < 4; i++)
	{
		if (node->nextNodes[i] != NULL)
		{
			free_node(node->nextNodes[i]);
			node->nextNodes[i] = NULL;
		}
	}
	
	free(node);
	node = NULL;
}

void split_node(Node *node, int divisions)
{
	if (divisions <= 0)
	{
		return;
	}

	double xMid = (node->xMax + node->xMin) / 2;
	double yMid = (node->yMax + node->yMin) / 2;
	
	Node *topLeft = malloc_node();
	topLeft->xMin = node->xMin;
	topLeft->xMax = xMid;
	topLeft->yMin = node->yMin;
	topLeft->yMax = yMid;
	node->nextNodes[0] = topLeft;
	
	Node *topRight = malloc_node();
	topRight->xMin = xMid;
	topRight->xMax = node->xMax;
	topRight->yMin = node->yMin;
	topRight->yMax = yMid;
	node->nextNodes[1] = topRight;
	
	Node *bottomLeft = malloc_node();
	bottomLeft->xMin = node->xMin;
	bottomLeft->xMax = xMid;
	bottomLeft->yMin = yMid;
	bottomLeft->yMax = node->yMax;
	node->nextNodes[2] = bottomLeft;
	
	Node *bottomRight = malloc_node();
	bottomRight->xMin = xMid;
	bottomRight->xMax = node->xMax;
	bottomRight->yMin = yMid;
	bottomRight->yMax = node->yMax;
	node->nextNodes[3] = bottomRight;
	
	for (int i = 0; i < 4; i++)
	{
		split_node(node->nextNodes[i], divisions - 1);
	}
}

void prepare_node(Node *node, int divisions, double xMin, double yMin,
                  double xMax, double yMax)
{
	node->xMin = xMin;
	node->yMin = yMin;
	node->xMax = xMax;
	node->yMax = yMax;
	
	split_node(node, divisions);
}

void add_to_node(Node *node, double x, double y, double z)
{
	if (node->nextNodes[0] == NULL)
	{
		node->value += z;
		node->count++;
		return;
	}
	
	double xMid = (node->xMax + node->xMin) / 2;
	double yMid = (node->yMax + node->yMin) / 2;

	int quarter = -1;

	if (x < xMid && y < yMid)
	{
		quarter = 0;
	}
	if (x > xMid && y < yMid)
	{
		quarter = 1;
	}
	else if (x < xMid && y > yMid)
	{
		quarter = 2;
	}
	else if (x > xMid && y > yMid)
	{
		quarter = 3;
	}
	
	if (quarter < 0)
	{
		return;
	}

	add_to_node(node->nextNodes[quarter], x, y, z);
}


