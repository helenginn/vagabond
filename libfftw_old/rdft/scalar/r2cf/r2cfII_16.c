/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Sun Oct 29 08:18:05 EDT 2017 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_r2cf.native -fma -compact -variables 4 -pipeline-latency 4 -n 16 -name r2cfII_16 -dft-II -include rdft/scalar/r2cfII.h */

/*
 * This function contains 66 FP additions, 48 FP multiplications,
 * (or, 18 additions, 0 multiplications, 48 fused multiply/add),
 * 32 stack variables, 7 constants, and 32 memory accesses
 */
#include "rdft/scalar/r2cfII.h"

static void r2cfII_16(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DK(KP198912367, +0.198912367379658006911597622644676228597850501);
     DK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DK(KP668178637, +0.668178637919298919997757686523080761552472251);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ivs, R1 = R1 + ivs, Cr = Cr + ovs, Ci = Ci + ovs, MAKE_VOLATILE_STRIDE(64, rs), MAKE_VOLATILE_STRIDE(64, csr), MAKE_VOLATILE_STRIDE(64, csi)) {
	       E T5, TZ, TB, TT, Tr, TK, Tu, TJ, Ti, TH, Tl, TG, Tc, T10, TE;
	       E TU;
	       {
		    E T1, TR, T4, TS, T2, T3;
		    T1 = R0[0];
		    TR = R0[WS(rs, 4)];
		    T2 = R0[WS(rs, 2)];
		    T3 = R0[WS(rs, 6)];
		    T4 = T2 - T3;
		    TS = T2 + T3;
		    T5 = FNMS(KP707106781, T4, T1);
		    TZ = FNMS(KP707106781, TS, TR);
		    TB = FMA(KP707106781, T4, T1);
		    TT = FMA(KP707106781, TS, TR);
	       }
	       {
		    E Tn, Ts, Tq, Tt, To, Tp;
		    Tn = R1[WS(rs, 7)];
		    Ts = R1[WS(rs, 3)];
		    To = R1[WS(rs, 1)];
		    Tp = R1[WS(rs, 5)];
		    Tq = To - Tp;
		    Tt = To + Tp;
		    Tr = FMA(KP707106781, Tq, Tn);
		    TK = FMA(KP707106781, Tt, Ts);
		    Tu = FNMS(KP707106781, Tt, Ts);
		    TJ = FMS(KP707106781, Tq, Tn);
	       }
	       {
		    E Te, Tj, Th, Tk, Tf, Tg;
		    Te = R1[0];
		    Tj = R1[WS(rs, 4)];
		    Tf = R1[WS(rs, 2)];
		    Tg = R1[WS(rs, 6)];
		    Th = Tf - Tg;
		    Tk = Tf + Tg;
		    Ti = FNMS(KP707106781, Th, Te);
		    TH = FMA(KP707106781, Tk, Tj);
		    Tl = FNMS(KP707106781, Tk, Tj);
		    TG = FMA(KP707106781, Th, Te);
	       }
	       {
		    E T8, TC, Tb, TD;
		    {
			 E T6, T7, T9, Ta;
			 T6 = R0[WS(rs, 5)];
			 T7 = R0[WS(rs, 1)];
			 T8 = FMA(KP414213562, T7, T6);
			 TC = FNMS(KP414213562, T6, T7);
			 T9 = R0[WS(rs, 3)];
			 Ta = R0[WS(rs, 7)];
			 Tb = FMA(KP414213562, Ta, T9);
			 TD = FMS(KP414213562, T9, Ta);
		    }
		    Tc = T8 - Tb;
		    T10 = TD - TC;
		    TE = TC + TD;
		    TU = T8 + Tb;
	       }
	       {
		    E Td, T13, Tw, T14, Tm, Tv;
		    Td = FMA(KP923879532, Tc, T5);
		    T13 = FNMS(KP923879532, T10, TZ);
		    Tm = FMA(KP668178637, Tl, Ti);
		    Tv = FMA(KP668178637, Tu, Tr);
		    Tw = Tm - Tv;
		    T14 = Tm + Tv;
		    Cr[WS(csr, 6)] = FNMS(KP831469612, Tw, Td);
		    Ci[WS(csi, 5)] = FNMS(KP831469612, T14, T13);
		    Cr[WS(csr, 1)] = FMA(KP831469612, Tw, Td);
		    Ci[WS(csi, 2)] = -(FMA(KP831469612, T14, T13));
	       }
	       {
		    E Tx, T11, TA, T12, Ty, Tz;
		    Tx = FNMS(KP923879532, Tc, T5);
		    T11 = FMA(KP923879532, T10, TZ);
		    Ty = FNMS(KP668178637, Tr, Tu);
		    Tz = FNMS(KP668178637, Ti, Tl);
		    TA = Ty - Tz;
		    T12 = Tz + Ty;
		    Cr[WS(csr, 5)] = FNMS(KP831469612, TA, Tx);
		    Ci[WS(csi, 1)] = FMA(KP831469612, T12, T11);
		    Cr[WS(csr, 2)] = FMA(KP831469612, TA, Tx);
		    Ci[WS(csi, 6)] = FMS(KP831469612, T12, T11);
	       }
	       {
		    E TF, TX, TM, TY, TI, TL;
		    TF = FMA(KP923879532, TE, TB);
		    TX = FNMS(KP923879532, TU, TT);
		    TI = FNMS(KP198912367, TH, TG);
		    TL = FMA(KP198912367, TK, TJ);
		    TM = TI + TL;
		    TY = TL - TI;
		    Cr[WS(csr, 7)] = FNMS(KP980785280, TM, TF);
		    Ci[WS(csi, 3)] = FMA(KP980785280, TY, TX);
		    Cr[0] = FMA(KP980785280, TM, TF);
		    Ci[WS(csi, 4)] = FMS(KP980785280, TY, TX);
	       }
	       {
		    E TN, TV, TQ, TW, TO, TP;
		    TN = FNMS(KP923879532, TE, TB);
		    TV = FMA(KP923879532, TU, TT);
		    TO = FMA(KP198912367, TG, TH);
		    TP = FNMS(KP198912367, TJ, TK);
		    TQ = TO - TP;
		    TW = TO + TP;
		    Cr[WS(csr, 4)] = FNMS(KP980785280, TQ, TN);
		    Ci[WS(csi, 7)] = FNMS(KP980785280, TW, TV);
		    Cr[WS(csr, 3)] = FMA(KP980785280, TQ, TN);
		    Ci[0] = -(FMA(KP980785280, TW, TV));
	       }
	  }
     }
}

static const kr2c_desc desc = { 16, "r2cfII_16", {18, 0, 48, 0}, &GENUS };

void X(codelet_r2cfII_16) (planner *p) {
     X(kr2c_register) (p, r2cfII_16, &desc);
}

#else

/* Generated by: ../../../genfft/gen_r2cf.native -compact -variables 4 -pipeline-latency 4 -n 16 -name r2cfII_16 -dft-II -include rdft/scalar/r2cfII.h */

/*
 * This function contains 66 FP additions, 30 FP multiplications,
 * (or, 54 additions, 18 multiplications, 12 fused multiply/add),
 * 32 stack variables, 7 constants, and 32 memory accesses
 */
#include "rdft/scalar/r2cfII.h"

static void r2cfII_16(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP555570233, +0.555570233019602224742830813948532874374937191);
     DK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DK(KP195090322, +0.195090322016128267848284868477022240927691618);
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ivs, R1 = R1 + ivs, Cr = Cr + ovs, Ci = Ci + ovs, MAKE_VOLATILE_STRIDE(64, rs), MAKE_VOLATILE_STRIDE(64, csr), MAKE_VOLATILE_STRIDE(64, csi)) {
	       E T5, T11, TB, TV, Tr, TK, Tu, TJ, Ti, TH, Tl, TG, Tc, T10, TE;
	       E TS;
	       {
		    E T1, TU, T4, TT, T2, T3;
		    T1 = R0[0];
		    TU = R0[WS(rs, 4)];
		    T2 = R0[WS(rs, 2)];
		    T3 = R0[WS(rs, 6)];
		    T4 = KP707106781 * (T2 - T3);
		    TT = KP707106781 * (T2 + T3);
		    T5 = T1 + T4;
		    T11 = TU - TT;
		    TB = T1 - T4;
		    TV = TT + TU;
	       }
	       {
		    E Tq, Tt, Tp, Ts, Tn, To;
		    Tq = R1[WS(rs, 7)];
		    Tt = R1[WS(rs, 3)];
		    Tn = R1[WS(rs, 1)];
		    To = R1[WS(rs, 5)];
		    Tp = KP707106781 * (Tn - To);
		    Ts = KP707106781 * (Tn + To);
		    Tr = Tp - Tq;
		    TK = Tt - Ts;
		    Tu = Ts + Tt;
		    TJ = Tp + Tq;
	       }
	       {
		    E Te, Tk, Th, Tj, Tf, Tg;
		    Te = R1[0];
		    Tk = R1[WS(rs, 4)];
		    Tf = R1[WS(rs, 2)];
		    Tg = R1[WS(rs, 6)];
		    Th = KP707106781 * (Tf - Tg);
		    Tj = KP707106781 * (Tf + Tg);
		    Ti = Te + Th;
		    TH = Tk - Tj;
		    Tl = Tj + Tk;
		    TG = Te - Th;
	       }
	       {
		    E T8, TC, Tb, TD;
		    {
			 E T6, T7, T9, Ta;
			 T6 = R0[WS(rs, 1)];
			 T7 = R0[WS(rs, 5)];
			 T8 = FNMS(KP382683432, T7, KP923879532 * T6);
			 TC = FMA(KP382683432, T6, KP923879532 * T7);
			 T9 = R0[WS(rs, 3)];
			 Ta = R0[WS(rs, 7)];
			 Tb = FNMS(KP923879532, Ta, KP382683432 * T9);
			 TD = FMA(KP923879532, T9, KP382683432 * Ta);
		    }
		    Tc = T8 + Tb;
		    T10 = Tb - T8;
		    TE = TC - TD;
		    TS = TC + TD;
	       }
	       {
		    E Td, TW, Tw, TR, Tm, Tv;
		    Td = T5 - Tc;
		    TW = TS + TV;
		    Tm = FMA(KP195090322, Ti, KP980785280 * Tl);
		    Tv = FNMS(KP980785280, Tu, KP195090322 * Tr);
		    Tw = Tm + Tv;
		    TR = Tv - Tm;
		    Cr[WS(csr, 4)] = Td - Tw;
		    Ci[WS(csi, 7)] = TR + TW;
		    Cr[WS(csr, 3)] = Td + Tw;
		    Ci[0] = TR - TW;
	       }
	       {
		    E Tx, TY, TA, TX, Ty, Tz;
		    Tx = T5 + Tc;
		    TY = TV - TS;
		    Ty = FNMS(KP195090322, Tl, KP980785280 * Ti);
		    Tz = FMA(KP980785280, Tr, KP195090322 * Tu);
		    TA = Ty + Tz;
		    TX = Tz - Ty;
		    Cr[WS(csr, 7)] = Tx - TA;
		    Ci[WS(csi, 3)] = TX + TY;
		    Cr[0] = Tx + TA;
		    Ci[WS(csi, 4)] = TX - TY;
	       }
	       {
		    E TF, T12, TM, TZ, TI, TL;
		    TF = TB + TE;
		    T12 = T10 - T11;
		    TI = FMA(KP831469612, TG, KP555570233 * TH);
		    TL = FMA(KP831469612, TJ, KP555570233 * TK);
		    TM = TI - TL;
		    TZ = TI + TL;
		    Cr[WS(csr, 6)] = TF - TM;
		    Ci[WS(csi, 2)] = T12 - TZ;
		    Cr[WS(csr, 1)] = TF + TM;
		    Ci[WS(csi, 5)] = -(TZ + T12);
	       }
	       {
		    E TN, T14, TQ, T13, TO, TP;
		    TN = TB - TE;
		    T14 = T10 + T11;
		    TO = FNMS(KP555570233, TJ, KP831469612 * TK);
		    TP = FNMS(KP555570233, TG, KP831469612 * TH);
		    TQ = TO - TP;
		    T13 = TP + TO;
		    Cr[WS(csr, 5)] = TN - TQ;
		    Ci[WS(csi, 1)] = T13 + T14;
		    Cr[WS(csr, 2)] = TN + TQ;
		    Ci[WS(csi, 6)] = T13 - T14;
	       }
	  }
     }
}

static const kr2c_desc desc = { 16, "r2cfII_16", {54, 18, 12, 0}, &GENUS };

void X(codelet_r2cfII_16) (planner *p) {
     X(kr2c_register) (p, r2cfII_16, &desc);
}

#endif
