//
//  Density2GL.h
//  VagabondViewer
//
//  Created by Helen Ginn on 22/7/2018.
//  Copyright © 2018 Ginn. All rights reserved.
//

#include <mutex>
#include "GLObject.h"
#include <vector>
#include "../libsrc/shared_ptrs.h"

struct ivec
{
	int x;
	int y;
	int z;
};

typedef enum
{
	DensityWeighted,
	DensityDifference,
	DensityOriginal,
	DensityDiffWithOriginal
} DensityType;

typedef std::map<int, std::map<int, int> > IntMap;

class GLKeeper;

class Density2GL : public GLObject
{
public:
	Density2GL();
	
	virtual void render();
	
	void makeNewDensity(CrystalPtr crystal = CrystalPtr());
	void nudgeDensity(int dir);
	
	void setOrigDensity()
	{
		_dType = DensityOriginal;
	}
	
	void setDiffDensity(bool val = true)
	{
		_dType = DensityDifference;
		_threshold = 3;
	}
	
	void setVisible(bool vis)
	{
		_visible = vis;
	}
	
	void recalculate()
	{
		_recalculate = true;
	}
	
	void setRenderType(GLuint type)
	{
		_renderType = type;
	}
	
	void setKeeper(GLKeeper *keeper)
	{
		_keeper = keeper;
	}
protected:
	virtual void bindTextures();
private:
	DensityType _dType;
	
	bool isDifferenceDensity()
	{
		return (_dType == DensityDiffWithOriginal || 
		        _dType == DensityDifference);
	}
	
	FFTPtr getFFT();
	CrystalPtr _crystal;
	GLKeeper *_keeper;
	void calculateContouring(CrystalPtr crystal);
	void makeUniformGrid();
	void setupIndexTable();
	void getSigma(FFTPtr fft);
	int _recalculate;
	std::mutex _renderLock;
	IntMap _flips;
	IntMap _allBits;
	
	vec3 _offset;
	vec3 getCentreOffset();
	double _resolution;
	double _threshold;
	double _sigma;
	double _mean;
	ivec _dims;
	bool _visible;
	std::vector<std::vector<GLuint> > _cubeIndices;
};
