// Vagabond
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

#ifndef __Vagabond__GeomVariant__
#define __Vagabond__GeomVariant__

class GeomVariant
{
public:
	static double *pointerToNonPGIVNonXProCNA()
	{
		return &_npgiv_nxpro_cna[0];
	}

	static double *pointerToNonPGIVNonXProNAB()
	{
		return &_npgiv_nxpro_nab[0];
	}

	static double *pointerToNonPGIVNonXProBAC()
	{
		return &_npgiv_nxpro_bac[0];
	}

	static double *pointerToNonPGIVNonXProNAC()
	{
		return &_npgiv_nxpro_nac[0];
	}

	static double *pointerToNonPGIVNonXProACN()
	{
		return &_npgiv_nxpro_acn[0];
	}

	static double *pointerToNonPGIVNonXProACO()
	{
		return &_npgiv_nxpro_aco[0];
	}

	static double *pointerToNonPGIVNonXProOCN()
	{
		return &_npgiv_nxpro_ocn[0];
	}

	static double *pointerToIleValNonXProCNA()
	{
		return &_ileval_nxpro_cna[0];
	}

	static double *pointerToIleValNonXProNAB()
	{
		return &_ileval_nxpro_nab[0];
	}

	static double *pointerToIleValNonXProBAC()
	{
		return &_ileval_nxpro_bac[0];
	}

	static double *pointerToIleValNonXProNAC()
	{
		return &_ileval_nxpro_nac[0];
	}

	static double *pointerToIleValNonXProACN()
	{
		return &_ileval_nxpro_acn[0];
	}

	static double *pointerToIleValNonXProACO()
	{
		return &_ileval_nxpro_aco[0];
	}

	static double *pointerToIleValNonXProOCN()
	{
		return &_ileval_nxpro_ocn[0];
	}

	/* for a given pointer, get phi and psi-dependent angle */
	static double getAngle(double *table, double phi, double psi);

	static double interpolateAngle(double *table, double phi, double psi);
private:
	/* 36 x 36 values lookup tables, psi = fast, phi/psi variation 10 deg
	 * from -180 to +170 */
	
	static double _npgiv_nxpro_cna[];
	static double _npgiv_nxpro_nab[];
	static double _npgiv_nxpro_nac[];
	static double _npgiv_nxpro_bac[];
	static double _npgiv_nxpro_acn[];
	static double _npgiv_nxpro_aco[];
	static double _npgiv_nxpro_ocn[];
	
	static double _ileval_nxpro_cna[];
	static double _ileval_nxpro_nab[];
	static double _ileval_nxpro_nac[];
	static double _ileval_nxpro_bac[];
	static double _ileval_nxpro_acn[];
	static double _ileval_nxpro_aco[];
	static double _ileval_nxpro_ocn[];

	static double _gly_nxpro_cna[];
	static double _gly_nxpro_nab[];
	static double _gly_nxpro_nac[];
	static double _gly_nxpro_bac[];
	static double _gly_nxpro_acn[];
	static double _gly_nxpro_aco[];
	static double _gly_nxpro_ocn[];

	static double _pro_nxpro_cna[];
	static double _pro_nxpro_nab[];
	static double _pro_nxpro_nac[];
	static double _pro_nxpro_bac[];
	static double _pro_nxpro_acn[];
	static double _pro_nxpro_aco[];
	static double _pro_nxpro_ocn[];

  static double _npgiv_xpro_cna[];
  static double _npgiv_xpro_nab[];
  static double _npgiv_xpro_nac[];
  static double _npgiv_xpro_bac[];
  static double _npgiv_xpro_aco[];
  static double _npgiv_xpro_acn[];
  static double _npgiv_xpro_ocn[];

  static double _ileval_xpro_cna[];
  static double _ileval_xpro_nab[];
  static double _ileval_xpro_nac[];
  static double _ileval_xpro_bac[];
  static double _ileval_xpro_aco[];
  static double _ileval_xpro_acn[];
  static double _ileval_xpro_ocn[];

  static double _gly_xpro_cna[];
  static double _gly_xpro_nab[];
  static double _gly_xpro_nac[];
  static double _gly_xpro_bac[];
  static double _gly_xpro_aco[];
  static double _gly_xpro_acn[];
  static double _gly_xpro_ocn[];

  static double _pro_xpro_cna[];
  static double _pro_xpro_nab[];
  static double _pro_xpro_nac[];
  static double _pro_xpro_bac[];
  static double _pro_xpro_aco[];
  static double _pro_xpro_acn[];
  static double _pro_xpro_ocn[];
};

#endif
