//
//  Vagabond2GL.cpp
//  VagabondViewer
//
//  Created by Helen Ginn on 03/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "Vagabond2GL.h"
#include "../libsrc/Options.h"
#include "../libsrc/Crystal.h"
#include "../libsrc/Bond.h"
#include "../libsrc/Atom.h"
#include "../libsrc/Element.h"

void Vagabond2GL::updateAtoms()
{
    for (AtomMap::iterator it = _atomMap.begin(); it != _atomMap.end(); it++)
    {
        AtomPtr atom = it->first;
        std::pair<int, int> pair = it->second;
        int conformer = pair.first;
        int vertex = pair.second;

        if (!atom->getModel()->isBond()) continue;

        std::vector<vec3> majBonds, minBonds;
        getPositions(atom, &minBonds, &majBonds);

        if (majBonds.size() >= conformer || minBonds.size() >= conformer)
        {
            continue;
        }

        vec3 majStart = majBonds[conformer];
        vec3 minStart = minBonds[conformer];

        _vertices[vertex].pos[0] = majStart.x;
        _vertices[vertex].pos[1] = majStart.y;
        _vertices[vertex].pos[2] = majStart.z;
        vertex++;

        for (int i = 0; i < 2; i++)
        {
            _vertices[vertex + i].pos[0] = (minStart.x + majStart.x) / 2;
            _vertices[vertex + i].pos[1] = (minStart.y + majStart.y) / 2;
            _vertices[vertex + i].pos[2] = (minStart.z + majStart.z) / 2;
        }

        vertex += 2;
        _vertices[vertex].pos[0] = minStart.x;
        _vertices[vertex].pos[1] = minStart.y;
        _vertices[vertex].pos[2] = minStart.z;

    }
}

void Vagabond2GL::getPositions(AtomPtr atom, std::vector<vec3> *min,
                               std::vector<vec3> *maj)
{
    ModelPtr minBond = atom->getModel();
    ModelPtr majBond = (ToBondPtr(minBond))->getMajor()->getModel();

    if (majBond->isBond())
    {
        *maj = ToBondPtr(majBond)->fishPositions();
    }

    *min = ToBondPtr(minBond)->fishPositions();
}

bool Vagabond2GL::shouldGetBonds()
{
    if (!_moleculeMap.size())
    {
        return true;
    }

    _renders = 0;

    OptionsPtr globalOptions = Options::getRuntimeOptions();

    for (int i = 0; i < globalOptions->crystalCount(); i++)
    {
        CrystalPtr crystal = globalOptions->getCrystal(i);

        for (int j = 0; j < crystal->moleculeCount(); j++)
        {
            MoleculePtr molecule = crystal->molecule(j);

            if (!_moleculeMap.count(molecule))
            {
                return true;
            }

            int expected = _moleculeMap[molecule];

            int existing = 0;

            for (int k = 0; k < molecule->atomCount(); k++)
            {
                if (molecule->atom(k)->getModel()->isBond())
                {
                    existing++;
                }
            }

            if (expected != existing)
            {
                return true;
            }
        }
    }

    return false;
}

void Vagabond2GL::setVertexColour(AtomPtr atom, Vertex *vertex)
{
    vertex->color[0] = 150. / 255.;
    vertex->color[1] = 150. / 255.;
    vertex->color[2] = 150. / 255.;
    vertex->color[3] = 1.0;

    if (atom->getElement()->getSymbol() == "O")
    {
        vertex->color[0] = 1.0;
        vertex->color[1] = 0.0;
        vertex->color[2] = 0.0;
    }

    if (atom->getElement()->getSymbol() == "N")
    {
        vertex->color[0] = 104. / 255.;
        vertex->color[1] = 139. / 255.;
        vertex->color[2] = 255. / 255.;
    }

    if (atom->getElement()->getSymbol() == "S")
    {
        vertex->color[0] = 255. / 255.;
        vertex->color[1] = 255. / 255.;
        vertex->color[2] = 0. / 255.;
    }
}

int Vagabond2GL::processMolecule(MoleculePtr molecule)
{
    GLuint count = (int)_vertices.size();
    int bonds = 0;

    for (int i = 0; i < molecule->atomCount(); i++)
    {
        AtomPtr atom = molecule->atom(i);

        if (atom->getModel() && atom->getModel()->isBond())
        {
            ToBondPtr(atom->getModel())->useMutex();
        }
    }

    for (int i = 0; i < molecule->atomCount(); i++)
    {
        AtomPtr atom = molecule->atom(i);
        AtomPtr major = ToBondPtr(atom->getModel())->getMajor();

        if (atom->getElement()->electronCount() <= 1)
        {
            continue;
        }

        if (atom->getModel()->isBond())
        {
            std::vector<vec3> majBonds, minBonds;
            getPositions(atom, &minBonds, &majBonds);

            for (int j = 0; j < majBonds.size(); j += 1)
            {
                vec3 majStart = majBonds[j];
                if (majBonds.size() <= j || minBonds.size() <= j)
                {
                    break;
                }


                vec3 minStart = minBonds[j];
                vec3 normal = vec3_subtract_vec3(minStart, majStart);
                vec3_set_length(&normal, 1);
                GLfloat glNorm[3];
                glNorm[0] = normal.x;
                glNorm[1] = normal.y;
                glNorm[2] = normal.z;

                Vertex vertex;
                vertex.pos[0] = majStart.x;
                vertex.pos[1] = majStart.y;
                vertex.pos[2] = majStart.z;
                memcpy(vertex.normal, &glNorm, 3 * sizeof(GLfloat));
                setVertexColour(major, &vertex);
                _vertices.push_back(vertex);

                /* Mid point */
                vertex.pos[0] = (minStart.x + majStart.x) / 2;
                vertex.pos[1] = (minStart.y + majStart.y) / 2;
                vertex.pos[2] = (minStart.z + majStart.z) / 2;
                memcpy(vertex.normal, &glNorm, 3 * sizeof(GLfloat));
                _vertices.push_back(vertex);
                setVertexColour(atom, &vertex);
                _vertices.push_back(vertex);

                /* Minor atom */
                vertex.pos[0] = minStart.x;
                vertex.pos[1] = minStart.y;
                vertex.pos[2] = minStart.z;
                memcpy(vertex.normal, &glNorm, 3 * sizeof(GLfloat));
                _vertices.push_back(vertex);
                _indices.push_back(count);
                _indices.push_back(count + 1);
                _indices.push_back(count + 2);
                _indices.push_back(count + 3);

                _atomMap[atom] = std::make_pair(j, count);
                count += 4;
            }

            if (minBonds.size() && majBonds.size())
            {
                bonds++;
            }
        }
    }

    return bonds;
}

void Vagabond2GL::findAtoms()
{
    OptionsPtr globalOptions = Options::getRuntimeOptions();

    _vertices.clear();
    _indices.clear();
    _atomMap.clear();
    _moleculeMap.clear();

    if (globalOptions)
    {
        for (int i = 0; i < globalOptions->crystalCount(); i++)
        {
            CrystalPtr crystal = globalOptions->getCrystal(i);

            for (int j = 0; j < crystal->moleculeCount(); j++)
            {
                MoleculePtr molecule = crystal->molecule(j);

                int expected = processMolecule(molecule);

                _moleculeMap[molecule] = expected;
            }
        }
    }
}

void Vagabond2GL::render()
{
    if (shouldGetBonds())
    {
        findAtoms();
        rebindProgram();
    }
    else
    {
        updateAtoms();
    }

    GLObject::render();
}
