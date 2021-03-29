// vagabond
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#ifndef __vagabond__BlobMesh__
#define __vagabond__BlobMesh__

#include <h3dsrc/Mesh.h>
#include "../libsrc/shared_ptrs.h"

class Blob;

class BlobMesh : public Mesh
{
public:
	BlobMesh(SlipObject *p, int tri);
	BlobPtr toBlob();
	
	void setBlob(Blob *b);
	
	void multiplyScale(double mult);
	
	void deleteBlob(CrystalPtr cryst);
	void blobToCrystal(CrystalPtr cryst);
	
	static BlobMesh *meshFromBlob(SlipObject *parent, Blob *b);
	
	bool hasBlob()
	{
		return (_blob != NULL);
	}

private:
	BlobPtr _blob;

};

#endif
