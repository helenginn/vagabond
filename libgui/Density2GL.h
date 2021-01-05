//
//  Density2GL.h
//  VagabondViewer
//
//  Created by Helen Ginn on 22/7/2018.
//  Copyright © 2018 Ginn. All rights reserved.
//

#include <mutex>
#include <h3dsrc/SlipObject.h>
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

class VagabondGLWidget;

class Density2GL : public SlipObject
{
public:
	Density2GL();
	
	virtual void render(SlipGL *sender);
	
	void makeNewDensity(CrystalPtr crystal = CrystalPtr());
	void nudgeDensity(int dir);
	
	void setOrigDensity()
	{
		_dType = DensityOriginal;
	}
	
	void setDiffWithOrigDensity()
	{
		_dType = DensityDiffWithOriginal;
		_imag = 1;
		_threshold = 3;
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
	
	void setKeeper(VagabondGLWidget *keeper)
	{
		_keeper = keeper;
	}
protected:
	virtual void bindTextures();
	virtual void extraUniforms();
private:
	DensityType _dType;
	
	bool isDifferenceDensity();
	
	VagFFTPtr getFFT();
	CrystalPtr _crystal;
	VagabondGLWidget *_keeper;
	void calculateContouring(CrystalPtr crystal);
	void makeUniformGrid();
	void setupIndexTable();
	void getSigma(VagFFTPtr fft);
	int _recalculate;
	int _imag;
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
