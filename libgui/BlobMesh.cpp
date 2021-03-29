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

#include "BlobMesh.h"
#include "../libsrc/Blob.h"
#include "../libsrc/Crystal.h"

BlobMesh::BlobMesh(SlipObject *p, int tri) : Mesh(p, tri)
{
	_wrapCycles = 20;
	_smoothCycles = 2;
}

BlobPtr BlobMesh::toBlob()
{
	if (!_blob)
	{
		_blob = BlobPtr(new Blob());
	}
	
	_blob->clear();
	
	changeToTriangles();

	for (size_t i = 0; i < vertexCount(); i++)
	{
		vec3 v = vec_from_pos(_vertices[i].pos);
		_blob->addVertex(v);
	}

	for (size_t i = 0; i < indexCount(); i++)
	{
		unsigned long idx = index(i);
		_blob->addIndex(idx);
	}
	
	changeToLines();
	
	return _blob;
}

void BlobMesh::blobToCrystal(CrystalPtr cryst)
{
	BlobPtr blob = toBlob();
	cryst->addBlob(blob);
}

void BlobMesh::deleteBlob(CrystalPtr cryst)
{
	if (_blob != NULL)
	{
		cryst->removeBlob(_blob);
	}
}

BlobMesh *BlobMesh::meshFromBlob(SlipObject *parent, Blob *b)
{
	BlobMesh *bm = new BlobMesh(parent, 1);
	parent->setCustomMesh(bm);
	
	bm->setBlob(b);
	bm->changeSelectionResize(1);
	bm->setSelectable(true);
	
	return bm;
}

void BlobMesh::setBlob(Blob *b)
{
	clearVertices();
	changeToTriangles();
	_blob = b->shared_from_this();

	for (size_t i = 0; i < b->vertexCount(); i++)
	{
		addVertex(b->vertex(i));
	}

	for (size_t i = 0; i < b->indexCount(); i++)
	{
		addIndex(b->index(i));
	}
	
	calculateNormals();

	changeToLines();

	setColour(0., 0., 0.);
}

void BlobMesh::multiplyScale(double mult)
{
	if (_blob != NULL)
	{
		_blob->multiplyScale(mult);
	}
}
