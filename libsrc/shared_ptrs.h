//
//  shared_ptrs.h
//  vagabond
//
//  Created by Helen Ginn on 11/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef vagabond_shared_ptrs_h
#define vagabond_shared_ptrs_h

#include <memory>

class cFFTW3d;
typedef std::shared_ptr<cFFTW3d> FFTPtr;

class Molecule;
typedef std::shared_ptr<Molecule> MoleculePtr;

class Atom;
typedef std::shared_ptr<Atom> AtomPtr;


class Model;
class Absolute;

typedef std::shared_ptr<Absolute> AbsolutePtr;
typedef std::shared_ptr<Model> ModelPtr;

#endif
