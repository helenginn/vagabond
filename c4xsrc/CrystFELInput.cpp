// cluster4x
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

#include "CrystFELInput.h"
#include "Group.h"
#include "MtzFile.h"
#include "MtzFFT.h"
#include <FileReader.h>
#include "AveDiffraction.h"
#include "libccp4/csymlib.h"
#include <crystfel/stream.h>
#include <crystfel/image.h>
#include <crystfel/symmetry.h>
#include <crystfel/reflist-utils.h>

CrystFELInput::CrystFELInput(std::string streams, std::string geom, 
                             std::string spg, double res)
{
	_max = 0;
	_skip = 0;
	_geom = geom;
	_spg = spg;
	_res = res;
	
	_streams = split(streams, ',');
}

static RefList *apply_max_adu(RefList *list, double max_adu)
{
	RefList *nlist;
	Reflection *refl;
	RefListIterator *iter;

	nlist = reflist_new();
	if ( nlist == NULL ) return NULL;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		if ( get_peak(refl) < max_adu ) {
			signed int h, k, l;
			get_indices(refl, &h, &k, &l);
			Reflection *nrefl = add_refl(nlist, h, k, l);
			copy_data(nrefl, refl);
		}
	}
	reflist_free(list);
	return nlist;
}

std::vector<struct image> CrystFELInput::loadStream(std::string str)
{
	Stream *stream = open_stream_for_read(str.c_str());

	std::vector<struct image> some;
	struct image *next;
	
	some.resize(some.size() + 1);
	next = &some[some.size() - 1];
	next->det = _detector;

	int n_crystals = 0;
	int n_crystals_seen = 0;
	Crystal **crystals = NULL;
	double max_adu = +INFINITY;

	char *sym_str = NULL;
	SymOpList *sym;
	
	std::vector<std::string> contents;
	
	if (_preload.length())
	{
		contents = split(get_file_contents(_preload), ',');
	}

	if ( sym_str == NULL ) sym_str = strdup("1");
	pointgroup_warning(sym_str);
	sym = get_pointgroup(sym_str);

	while (true)
	{
//		ReadStreamFlags f = STREAM_READ_REFLECTIONS | STREAM_READ_UNITCELL;
		if (read_chunk(stream, next) != 0 )
		{
			break;
		}

		struct image *cur = &some[some.size() - 1];
		RefList *as;

		for (int i = 0; i<cur->n_crystals; i++)
		{
			n_crystals_seen++;
			
			if (n_crystals_seen < _skip)
			{
				continue;
			}
			
			Crystal *cr;
			Crystal **crystals_new;
			RefList *cr_refl;
			struct image *image;

			crystals_new = (Crystal **)realloc(crystals,
			                      (n_crystals+1)*sizeof(Crystal *));
			if ( crystals_new == NULL ) {
				ERROR("Failed to allocate memory for crystal "
				      "list.\n");
				return std::vector<struct image>();
			}
			crystals = crystals_new;
			crystals[n_crystals] = cur->crystals[i];
			cr = crystals[n_crystals];

			image = (struct image *)malloc(sizeof(struct image));
			if ( image == NULL ) {
				ERROR("Failed to allocatea memory for image.\n");
				return std::vector<struct image>();
			}

			crystal_set_image(cr, image);
			*image = *cur;
			image->n_crystals = 1;
			image->crystals = &crystals[n_crystals];

			/* This is the raw list of reflections */
			cr_refl = crystal_get_reflections(cr);

			cr_refl = apply_max_adu(cr_refl, max_adu);

			as = asymmetric_indices(cr_refl, sym);
			crystal_set_reflections(cr, as);
			crystal_set_user_flag(cr, 0);
			reflist_free(cr_refl);

			n_crystals++;

			if (n_crystals % 1000 == 0)
			{
				std::cout << "." << std::flush;
			}
		}

		if (n_crystals_seen < _skip)
		{
			continue;
		}
		
		if (contents.size() > 0)
		{
			std::string fn = next->filename;
			std::string bn = getBaseFilename(fn);

			/* not found */
			if (std::find(contents.begin(), contents.end(), bn)
			    == contents.end())
			{
				some.pop_back();
			}
		}

		some.resize(some.size() + 1);
		next = &some[some.size() - 1];
		next->det = _detector;
		next->div = NAN;
		next->bw = NAN;

		if (_max > 0 && n_crystals > _max)
		{
			break;
		}
	}
	
	std::cout << "Looked through " << n_crystals << " crystals from " <<
	"stream: " << str << std::endl;
	
	some.pop_back();
	
	close_stream(stream);
	
	return some;
}

std::vector<double> unitCellFor(Crystal *c)
{
	UnitCell *cell = crystal_get_cell(c);
	std::vector<double> uc = std::vector<double>(6, 0.);
	cell_get_parameters(cell, &uc[0], &uc[1], &uc[2],
	                    &uc[3], &uc[4], &uc[5]);
	uc[0] *= 1e10;
	uc[1] *= 1e10;
	uc[2] *= 1e10;
	uc[3] *= 180 / M_PI;
	uc[4] *= 180 / M_PI;
	uc[5] *= 180 / M_PI;

	return uc;
}

Group *CrystFELInput::process()
{
	AveDiffraction::setShouldScale(false);

	Group *g = new Group(NULL);
	int count = 0;
	
	std::vector<struct image> images;
	_detector = get_detector_geometry(_geom.c_str(), NULL);

	for (size_t i = 0; i < _streams.size(); i++)
	{
		std::cout << "Loading " << _streams[i] << std::endl;
		std::vector<struct image> some = loadStream(_streams[i]);
		images.reserve(some.size() + images.size());
		images.insert(images.end(), some.begin(), some.end());
	}
	
	int mh = 0; int mk = 0; int ml = 0;
	/* find maximum vals for h,k,l */
	for (size_t i = 0; i < images.size(); i++)
	{
		struct image *im = &images[i];
		std::string name = std::string(im->filename) + i_to_str(i);

		for (int j = 0; j < im->n_crystals; j++)
		{
		count++;
			Crystal *c = im->crystals[j];
			std::vector<double> uc = unitCellFor(c);
			mat3x3 m = mat3x3_from_unit_cell(&uc[0]);
			m = mat3x3_inverse(m);
			m = mat3x3_transpose(m);
			RefListIterator *it;
			RefList *refls = crystal_get_reflections(c);
			Reflection *ref = first_refl(refls, &it);

			while (true)
			{
				ref = (next_refl(ref, it));
				
				if (ref == NULL)
				{
					break;
				}

				int h, k, l;
				get_indices(ref, &h, &k, &l);
				vec3 hkl = make_vec3(h, k, l);
				mat3x3_mult_vec(m, &hkl);
				
				double r = 1 / vec3_length(hkl);
				if (r < _res)
				{
					continue;
				}
				
				h = abs(h); k = abs(k); l = abs(l);
				
				mh = std::max(h, mh);
				mk = std::max(k, mk);
				ml = std::max(l, ml);
			}
		}
	}
	
	mh *= 2; mk *= 2; ml *= 2;
	VagFFTPtr fftTemplate = VagFFTPtr(new VagFFT(mh, mk, ml));
	count = 0;

	int num = atoi(_spg.c_str());
	if (num == 0)
	{
		std::cout << "Please enter space group, --spg=<num>" << std::endl;
		exit(1);
	}

	CSym::CCP4SPG *spg = CSym::ccp4spg_load_by_ccp4_num(num);
	
	for (size_t i = 0; i < images.size(); i++)
	{
		struct image *im = &images[i];
		std::string name = std::string(im->filename) + "_" + i_to_str(i);

		for (int j = 0; j < im->n_crystals; j++)
		{
			count++;
			Crystal *c = im->crystals[j];
			std::vector<double> uc = unitCellFor(c);
			mat3x3 m = mat3x3_from_unit_cell(&uc[0]);
			m = mat3x3_inverse(m);
			m = mat3x3_transpose(m);

			std::string subname = name + "_" + i_to_str(j);
			std::string nickname = getBaseFilename(name);
			MtzFile *file = new MtzFile(name);
			file->setPanddaName(subname);
			file->setMetadata(nickname);

			MtzFFTPtr mtz = MtzFFTPtr(new MtzFFT(g, *fftTemplate));
			mtz->multiplyAll(NAN);
			mtz->setMtzFile(file);
			mtz->setUnitCell(uc);
			mtz->setSpaceGroup(spg);

			RefListIterator *it;
			RefList *refls = crystal_get_reflections(c);
			Reflection *ref = first_refl(refls, &it);

			while (true)
			{
				ref = (next_refl(ref, it));
				
				if (ref == NULL)
				{
					break;
				}

				int h, k, l, _h, _k, _l;
				get_indices(ref, &h, &k, &l);
				CSym::ccp4spg_put_in_asu(spg, h, k, l, &_h, &_k, &_l);
				vec3 hkl = make_vec3(h, k, l);
				mat3x3_mult_vec(m, &hkl);
				
				double r = 1 / vec3_length(hkl);
				if (r < _res)
				{
					continue;
				}

				double intensity = get_intensity(ref);
				
				mtz->setReal(_h, _k, _l, intensity);
				mtz->setImag(_h, _k, _l, 0);
			}

			mtz->applySymmetry(true, -1, false);
			g->addMtz(mtz);
		}
	}

	return g;
}
