// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_INSTANTIATE
#define H_GMSHFEM_INSTANTIATE

#define _GET_ARGS_1(a, ...) a
#define _GET_ARGS_2(a, b, ...) b
#define _GET_ARGS_3(a, b, c, ...) c
#define _GET_ARGS_4(a, b, c, d, ...) d
#define _GET_ARGS_5(a, b, c, d, e, ...) e
#define _GET_ARGS_6(a, b, c, d, e, f, ...) f
#define _GET_ARGS_7(a, b, c, d, e, f, g, ...) g
#define _GET_ARGS_8(a, b, c, d, e, f, g, h, ...) h
#define _GET_ARGS_9(a, b, c, d, e, f, g, h, i, ...) i
#define _GET_ARGS_10(a, b, c, d, e, f, g, h, i, j, ...) j

#define GET_ARGS(I, ...) _GET_ARGS_##I(__VA_ARGS__)

#define TEMPLATE_ARGS(...) __VA_ARGS__

// ************************
// Class
// ************************

// 1 arg
#define _INSTANTIATE_CLASS_N(N, name, argsA) template class name< GET_ARGS(N, argsA) >;

#define _INSTANTIATE_CLASS_1(name, argsA) _INSTANTIATE_CLASS_N(1, name, TEMPLATE_ARGS(argsA))
#define _INSTANTIATE_CLASS_2(name, argsA)             \
  _INSTANTIATE_CLASS_N(2, name, TEMPLATE_ARGS(argsA)) \
  _INSTANTIATE_CLASS_1(name, TEMPLATE_ARGS(argsA));
#define _INSTANTIATE_CLASS_3(name, argsA)             \
  _INSTANTIATE_CLASS_N(3, name, TEMPLATE_ARGS(argsA)) \
  _INSTANTIATE_CLASS_2(name, TEMPLATE_ARGS(argsA));
#define _INSTANTIATE_CLASS_4(name, argsA)             \
  _INSTANTIATE_CLASS_N(4, name, TEMPLATE_ARGS(argsA)) \
  _INSTANTIATE_CLASS_3(name, TEMPLATE_ARGS(argsA));
#define _INSTANTIATE_CLASS_5(name, argsA)             \
  _INSTANTIATE_CLASS_N(5, name, TEMPLATE_ARGS(argsA)) \
  _INSTANTIATE_CLASS_4(name, TEMPLATE_ARGS(argsA));

#define INSTANTIATE_CLASS(name, sizeA, argsA) _INSTANTIATE_CLASS_##sizeA(name, TEMPLATE_ARGS(argsA))

// 2 args
#define _INSTANTIATE_CLASS_M_N(M, N, name, argsA, argsB) template class name< GET_ARGS(M, argsA), GET_ARGS(N, argsB) >;

#define _INSTANTIATE_CLASS_2_N(N, name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_M_N(1, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_M_N(2, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_N(N, name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_2_N(N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))    \
  _INSTANTIATE_CLASS_M_N(3, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_M_N(4, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_5_N(N, name, argsA, argsB)                         \
  _INSTANTIATE_CLASS_4_N(N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_M_N(5, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))

#define _INSTANTIATE_CLASS_2_1(name, argsA, argsB) _INSTANTIATE_CLASS_2_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_2_2(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_2_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_2_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_2_3(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_2_N(3, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_2_2(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_2_4(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_2_N(4, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_2_3(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))

#define _INSTANTIATE_CLASS_4_1(name, argsA, argsB) _INSTANTIATE_CLASS_4_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_2(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_3(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(3, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_2(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_4(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(4, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_3(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))

#define _INSTANTIATE_CLASS_4_1(name, argsA, argsB) _INSTANTIATE_CLASS_4_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_2(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_3(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(3, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_2(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_4(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(4, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_3(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_5(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(5, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_4(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_6(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(6, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_5(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_7(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(7, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_6(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_8(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(8, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_7(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_4_9(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_4_N(9, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_4_8(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))

#define _INSTANTIATE_CLASS_5_1(name, argsA, argsB) _INSTANTIATE_CLASS_5_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_5_2(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_5_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_5_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_5_3(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_5_N(3, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_5_2(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_5_4(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_5_N(4, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_5_3(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define _INSTANTIATE_CLASS_5_5(name, argsA, argsB)                            \
  _INSTANTIATE_CLASS_5_N(5, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB)) \
  _INSTANTIATE_CLASS_5_4(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))

#define INSTANTIATE_CLASS_2(name, sizeA, sizeB, argsA, argsB) _INSTANTIATE_CLASS_##sizeA##_##sizeB(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))

// 3 args
#define _INSTANTIATE_CLASS_O_M_N(O, M, N, name, argsA, argsB, argsC) template class name< GET_ARGS(O, argsA), GET_ARGS(M, argsB), GET_ARGS(N, argsC) >;

#define _INSTANTIATE_CLASS_4_M_N(M, N, name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_O_M_N(1, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_O_M_N(2, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_O_M_N(3, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_O_M_N(4, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))

#define _INSTANTIATE_CLASS_4_1_N(N, name, argsA, argsB, argsC) \
  _INSTANTIATE_CLASS_4_M_N(1, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_2_N(N, name, argsA, argsB, argsC)                                        \
  _INSTANTIATE_CLASS_4_1_N(N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_M_N(2, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_3_N(N, name, argsA, argsB, argsC)                                        \
  _INSTANTIATE_CLASS_4_2_N(N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_M_N(3, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_4_N(N, name, argsA, argsB, argsC)                                        \
  _INSTANTIATE_CLASS_4_3_N(N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_M_N(4, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))

#define _INSTANTIATE_CLASS_4_1_1(name, argsA, argsB, argsC) _INSTANTIATE_CLASS_4_1_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_1_2(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_1_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_1_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))

#define _INSTANTIATE_CLASS_4_2_1(name, argsA, argsB, argsC) _INSTANTIATE_CLASS_4_2_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_2_2(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_2_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_2_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_2_3(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_2_N(3, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_2_2(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_2_4(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_2_N(4, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_2_3(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))

#define _INSTANTIATE_CLASS_4_3_1(name, argsA, argsB, argsC) _INSTANTIATE_CLASS_4_3_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_3_2(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_3_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_3_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_3_3(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_3_N(3, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_3_2(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_3_4(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_3_N(4, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_3_3(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))

#define _INSTANTIATE_CLASS_4_4_1(name, argsA, argsB, argsC) _INSTANTIATE_CLASS_4_4_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_4_2(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_4_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_4_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_4_3(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_4_N(3, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_4_2(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define _INSTANTIATE_CLASS_4_4_4(name, argsA, argsB, argsC)                                           \
  _INSTANTIATE_CLASS_4_4_N(4, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC)) \
  _INSTANTIATE_CLASS_4_4_3(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))

#define INSTANTIATE_CLASS_3(name, sizeA, sizeB, sizeC, argsA, argsB, argsC) _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))

// 4 args
#define _INSTANTIATE_CLASS_P_O_M_N(P, O, M, N, name, argsA, argsB, argsC, argsD) template class name< GET_ARGS(P, argsA), GET_ARGS(O, argsB), GET_ARGS(M, argsC), GET_ARGS(N, argsD) >;

#define _INSTANTIATE_CLASS_4_O_M_N(O, M, N, name, argsA, argsB, argsC, argsD)                                                          \
  _INSTANTIATE_CLASS_P_O_M_N(1, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_P_O_M_N(2, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_P_O_M_N(3, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_P_O_M_N(4, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))

#define _INSTANTIATE_CLASS_4_1_M_N(M, N, name, argsA, argsB, argsC, argsD) \
  _INSTANTIATE_CLASS_4_O_M_N(1, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define _INSTANTIATE_CLASS_4_3_M_N(M, N, name, argsA, argsB, argsC, argsD)                                                          \
  _INSTANTIATE_CLASS_4_1_M_N(M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))    \
  _INSTANTIATE_CLASS_4_O_M_N(2, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_4_O_M_N(3, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))

#define _INSTANTIATE_CLASS_4_1_2_N(N, name, argsA, argsB, argsC, argsD)                                                          \
  _INSTANTIATE_CLASS_4_1_M_N(1, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_4_1_M_N(2, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define _INSTANTIATE_CLASS_4_3_4_N(N, name, argsA, argsB, argsC, argsD)                                                          \
  _INSTANTIATE_CLASS_4_3_M_N(1, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_4_3_M_N(2, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_4_3_M_N(3, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_4_3_M_N(4, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))

#define _INSTANTIATE_CLASS_4_1_2_1(name, argsA, argsB, argsC, argsD) \
  _INSTANTIATE_CLASS_4_1_2_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define _INSTANTIATE_CLASS_4_1_2_2(name, argsA, argsB, argsC, argsD)                                                       \
  _INSTANTIATE_CLASS_4_1_2_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_4_1_2_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define _INSTANTIATE_CLASS_4_3_4_1(name, argsA, argsB, argsC, argsD) \
  _INSTANTIATE_CLASS_4_3_4_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define _INSTANTIATE_CLASS_4_3_4_2(name, argsA, argsB, argsC, argsD)                                                       \
  _INSTANTIATE_CLASS_4_3_4_1(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD)) \
  _INSTANTIATE_CLASS_4_3_4_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))

#define INSTANTIATE_CLASS_4(name, sizeA, sizeB, sizeC, sizeD, argsA, argsB, argsC, argsD) _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))

// 5 args
#define _INSTANTIATE_CLASS_Q_P_O_M_N(Q, P, O, M, N, name, argsA, argsB, argsC, argsD, argsE) template class name< GET_ARGS(Q, argsA), GET_ARGS(P, argsB), GET_ARGS(O, argsC), GET_ARGS(M, argsD), GET_ARGS(N, argsE) >;

#define _INSTANTIATE_CLASS_4_P_O_M_N(P, O, M, N, name, argsA, argsB, argsC, argsD, argsE)                                                                         \
  _INSTANTIATE_CLASS_Q_P_O_M_N(1, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_Q_P_O_M_N(2, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_Q_P_O_M_N(3, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_Q_P_O_M_N(4, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))

#define _INSTANTIATE_CLASS_4_1_O_M_N(O, M, N, name, argsA, argsB, argsC, argsD, argsE) \
  _INSTANTIATE_CLASS_4_P_O_M_N(1, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))
#define _INSTANTIATE_CLASS_4_3_O_M_N(O, M, N, name, argsA, argsB, argsC, argsD, argsE)                                                                         \
  _INSTANTIATE_CLASS_4_1_O_M_N(O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))    \
  _INSTANTIATE_CLASS_4_P_O_M_N(2, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_4_P_O_M_N(3, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))

#define _INSTANTIATE_CLASS_4_1_2_M_N(M, N, name, argsA, argsB, argsC, argsD, argsE)                                                                         \
  _INSTANTIATE_CLASS_4_1_O_M_N(1, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_4_1_O_M_N(2, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))
#define _INSTANTIATE_CLASS_4_3_4_M_N(M, N, name, argsA, argsB, argsC, argsD, argsE)                                                                         \
  _INSTANTIATE_CLASS_4_3_O_M_N(1, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_4_3_O_M_N(2, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_4_3_O_M_N(3, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_4_3_O_M_N(4, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))

#define _INSTANTIATE_CLASS_4_1_2_1_N(N, name, argsA, argsB, argsC, argsD, argsE) \
  _INSTANTIATE_CLASS_4_1_2_M_N(1, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))
#define _INSTANTIATE_CLASS_4_3_4_1_N(N, name, argsA, argsB, argsC, argsD, argsE) \
  _INSTANTIATE_CLASS_4_3_4_M_N(1, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))

#define _INSTANTIATE_CLASS_4_1_2_1_2(name, argsA, argsB, argsC, argsD, argsE)                                                                         \
  _INSTANTIATE_CLASS_4_1_2_1_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_4_1_2_1_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))
#define _INSTANTIATE_CLASS_4_3_4_1_2(name, argsA, argsB, argsC, argsD, argsE)                                                                         \
  _INSTANTIATE_CLASS_4_3_4_1_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE)) \
  _INSTANTIATE_CLASS_4_3_4_1_N(2, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))

#define INSTANTIATE_CLASS_5(name, sizeA, sizeB, sizeC, sizeD, sizeE, argsA, argsB, argsC, argsD, argsE) _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD##_##sizeE(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE))

// 7 args
#define _INSTANTIATE_CLASS_S_R_Q_P_O_M_N(S, R, Q, P, O, M, N, name, argsA, argsB, argsC, argsD, argsE, argsF, argsG) \
  template class name< GET_ARGS(S, argsA), GET_ARGS(R, argsB), GET_ARGS(Q, argsC), GET_ARGS(P, argsD), GET_ARGS(O, argsE), GET_ARGS(M, argsF), GET_ARGS(N, argsG) >;

#define _INSTANTIATE_CLASS_4_R_Q_P_O_M_N(R, Q, P, O, M, N, name, argsA, argsB, argsC, argsD, argsE, argsF, argsG)                                                                                                       \
  _INSTANTIATE_CLASS_S_R_Q_P_O_M_N(1, R, Q, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG)) \
  _INSTANTIATE_CLASS_S_R_Q_P_O_M_N(2, R, Q, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG)) \
  _INSTANTIATE_CLASS_S_R_Q_P_O_M_N(3, R, Q, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG)) \
  _INSTANTIATE_CLASS_S_R_Q_P_O_M_N(4, R, Q, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG))

#define _INSTANTIATE_CLASS_4_1_Q_P_O_M_N(Q, P, O, M, N, name, argsA, argsB, argsC, argsD, argsE, argsF, argsG) \
  _INSTANTIATE_CLASS_4_R_Q_P_O_M_N(1, Q, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG))

#define _INSTANTIATE_CLASS_4_1_1_P_O_M_N(P, O, M, N, name, argsA, argsB, argsC, argsD, argsE, argsF, argsG) \
  _INSTANTIATE_CLASS_4_1_Q_P_O_M_N(1, P, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG))

#define _INSTANTIATE_CLASS_4_1_1_2_O_M_N(O, M, N, name, argsA, argsB, argsC, argsD, argsE, argsF, argsG)                                                                                                       \
  _INSTANTIATE_CLASS_4_1_1_P_O_M_N(1, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG)) \
  _INSTANTIATE_CLASS_4_1_1_P_O_M_N(2, O, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG))

#define _INSTANTIATE_CLASS_4_1_1_2_1_M_N(M, N, name, argsA, argsB, argsC, argsD, argsE, argsF, argsG) \
  _INSTANTIATE_CLASS_4_1_1_2_O_M_N(1, M, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG))

#define _INSTANTIATE_CLASS_4_1_1_2_1_1_N(N, name, argsA, argsB, argsC, argsD, argsE, argsF, argsG) \
  _INSTANTIATE_CLASS_4_1_1_2_1_M_N(1, N, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG))

#define _INSTANTIATE_CLASS_4_1_1_2_1_1_1(name, argsA, argsB, argsC, argsD, argsE, argsF, argsG) \
  _INSTANTIATE_CLASS_4_1_1_2_1_1_N(1, name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG))

#define INSTANTIATE_CLASS_7(name, sizeA, sizeB, sizeC, sizeD, sizeE, sizeF, sizeG, argsA, argsB, argsC, argsD, argsE, argsF, argsG) _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD##_##sizeE##_##sizeF##_##sizeG(name, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD), TEMPLATE_ARGS(argsE), TEMPLATE_ARGS(argsF), TEMPLATE_ARGS(argsG))

// ************************
// Function
// ************************

#define TEMPLATE_PARAM_1 T_TypeA
#define TEMPLATE_PARAM_2 T_TypeB
#define TEMPLATE_PARAM_3 T_TypeC
#define TEMPLATE_PARAM_4 T_TypeD
#define TEMPLATE_CLASS_PARAM_1 T_TypeClassA

#define TEMPLATE_PARAMS(...) __VA_ARGS__
#define TEMPLATE_RETURN(...) __VA_ARGS__

#define PACK(...) __VA_ARGS__

// 1 arg
#define _FAKE_CLASS(ret, namespace, name, typeA, params, tag) \
  template< typeA T_TypeA >                                   \
  class FAKE_CLASS_##name##_##tag                             \
  {                                                           \
    auto fake()                                               \
    {                                                         \
      ret (*pf)(params) = namespace##name< T_TypeA >;         \
      return pf;                                              \
    }                                                         \
  };

#define _FAKE_CLASS_OPLL(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                  \
  class FAKE_CLASS_operatorLL##_##tag                        \
  {                                                          \
    auto fake()                                              \
    {                                                        \
      ret (*pf)(params) = namespace##operator<< < T_TypeA >; \
      return pf;                                             \
    }                                                        \
  };

#define _FAKE_CLASS_OPP(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                 \
  class FAKE_CLASS_operatorP##_##tag                        \
  {                                                         \
    auto fake()                                             \
    {                                                       \
      ret (*pf)(params) = namespace##operator+< T_TypeA >;  \
      return pf;                                            \
    }                                                       \
  };

#define _FAKE_CLASS_OPM(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                 \
  class FAKE_CLASS_operatorM##_##tag                        \
  {                                                         \
    auto fake()                                             \
    {                                                       \
      ret (*pf)(params) = namespace##operator-< T_TypeA >;  \
      return pf;                                            \
    }                                                       \
  };

#define _FAKE_CLASS_OPT(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                 \
  class FAKE_CLASS_operatorT##_##tag                        \
  {                                                         \
    auto fake()                                             \
    {                                                       \
      ret (*pf)(params) = namespace##operator*< T_TypeA >;  \
      return pf;                                            \
    }                                                       \
  };

#define _FAKE_CLASS_OPD(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                 \
  class FAKE_CLASS_operatorD##_##tag                        \
  {                                                         \
    auto fake()                                             \
    {                                                       \
      ret (*pf)(params) = namespace##operator/< T_TypeA >;  \
      return pf;                                            \
    }                                                       \
  };

#define _FAKE_CLASS_OPV(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                 \
  class FAKE_CLASS_operatorV##_##tag                        \
  {                                                         \
    auto fake()                                             \
    {                                                       \
      ret (*pf)(params) = namespace##operator%< T_TypeA >;  \
      return pf;                                            \
    }                                                       \
  };

#define _FAKE_CLASS_OPE(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                 \
  class FAKE_CLASS_operatorE##_##tag                        \
  {                                                         \
    auto fake()                                             \
    {                                                       \
      ret (*pf)(params) = namespace##operator=< T_TypeA >;  \
      return pf;                                            \
    }                                                       \
  };

#define _FAKE_CLASS_OPPP(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                  \
  class FAKE_CLASS_operatorPP##_##tag                        \
  {                                                          \
    auto fake()                                              \
    {                                                        \
      ret (*pf)(params) = namespace##operator+< T_TypeA >;   \
      return pf;                                             \
    }                                                        \
  };

#define _FAKE_CLASS_OPMM(ret, namespace, typeA, params, tag) \
  template< typeA T_TypeA >                                  \
  class FAKE_CLASS_operatorMM##_##tag                        \
  {                                                          \
    auto fake()                                              \
    {                                                        \
      ret (*pf)(params) = namespace##operator-< T_TypeA >;   \
      return pf;                                             \
    }                                                        \
  };

#define INSTANTIATE_FCT(ret, namespace, name, tag, sizeA, typeA, argsA, params)           \
  _FAKE_CLASS(TEMPLATE_RETURN(ret), namespace, name, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_##name##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPLL(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPLL(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorLL##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPP(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPP(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorP##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPM(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPM(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorM##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPT(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPT(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorT##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPD(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPD(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorD##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPV(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPV(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorV##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPE(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPE(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorE##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPPP(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPPP(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorPP##_##tag, TEMPLATE_ARGS(argsA))
#define INSTANTIATE_OPMM(ret, namespace, tag, sizeA, typeA, argsA, params)               \
  _FAKE_CLASS_OPMM(TEMPLATE_RETURN(ret), namespace, typeA, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA(FAKE_CLASS_operatorMM##_##tag, TEMPLATE_ARGS(argsA))


// 2 args
#define _FAKE_CLASS_1_1(ret, namespace, className, name, typeClassA, typeA, params, tag) \
  template< typeClassA T_TypeClassA, typeA T_TypeA >                                     \
  class FAKE_CLASS_##className##_##name##_##tag                                          \
  {                                                                                      \
    auto fake()                                                                          \
    {                                                                                    \
      typedef ret (className< T_TypeClassA >::*Fct)(params);                             \
      Fct pf = &className< T_TypeClassA >::template name< T_TypeA >;                     \
      return pf;                                                                         \
    }                                                                                    \
  };

#define _FAKE_CLASS_1_1_CONST(ret, namespace, className, name, typeClassA, typeA, params, tag) \
  template< typeClassA T_TypeClassA, typeA T_TypeA >                                           \
  class FAKE_CLASS_##className##_##name##_##tag                                                \
  {                                                                                            \
    auto fake()                                                                                \
    {                                                                                          \
      typedef ret (className< T_TypeClassA >::*Fct)(params) const;                             \
      Fct pf = &className< T_TypeClassA >::template name< T_TypeA >;                           \
      return pf;                                                                               \
    }                                                                                          \
  };

#define _FAKE_CLASS_OPLL_1_1(ret, namespace, className, typeClassA, typeA, params, tag) \
  template< typeClassA T_TypeClassA, typeA T_TypeA >                                    \
  class FAKE_CLASS_##className##_##operatorLL##_##tag                                   \
  {                                                                                     \
    auto fake()                                                                         \
    {                                                                                   \
      typedef ret (className< T_TypeClassA >::*Fct)(params);                            \
      Fct pf = &className< T_TypeClassA >::template operator<< < T_TypeA >;             \
      return pf;                                                                        \
    }                                                                                   \
  };

#define _FAKE_CLASS_2(ret, namespace, name, typeA, typeB, params, tag) \
  template< typeA T_TypeA, typeB T_TypeB >                             \
  class FAKE_CLASS_##name##_##tag                                      \
  {                                                                    \
    auto fake()                                                        \
    {                                                                  \
      typedef ret (*Fct)(params);                                      \
      Fct pf = namespace##name< T_TypeA, T_TypeB >;                    \
      return pf;                                                       \
    }                                                                  \
  };

#define _FAKE_CLASS_OPP_2(ret, namespace, typeA, typeB, params, tag) \
  template< typeA T_TypeA, typeB T_TypeB >                           \
  class FAKE_CLASS_operatorP##_##tag                                 \
  {                                                                  \
    auto fake()                                                      \
    {                                                                \
      ret (*pf)(params) = namespace##operator+< T_TypeA, T_TypeB >;  \
      return pf;                                                     \
    }                                                                \
  };

#define _FAKE_CLASS_OPM_2(ret, namespace, typeA, typeB, params, tag) \
  template< typeA T_TypeA, typeB T_TypeB >                           \
  class FAKE_CLASS_operatorM##_##tag                                 \
  {                                                                  \
    auto fake()                                                      \
    {                                                                \
      ret (*pf)(params) = namespace##operator-< T_TypeA, T_TypeB >;  \
      return pf;                                                     \
    }                                                                \
  };

#define _FAKE_CLASS_OPT_2(ret, namespace, typeA, typeB, params, tag) \
  template< typeA T_TypeA, typeB T_TypeB >                           \
  class FAKE_CLASS_operatorT##_##tag                                 \
  {                                                                  \
    auto fake()                                                      \
    {                                                                \
      ret (*pf)(params) = namespace##operator*< T_TypeA, T_TypeB >;  \
      return pf;                                                     \
    }                                                                \
  };

#define _FAKE_CLASS_OPD_2(ret, namespace, typeA, typeB, params, tag) \
  template< typeA T_TypeA, typeB T_TypeB >                           \
  class FAKE_CLASS_operatorD##_##tag                                 \
  {                                                                  \
    auto fake()                                                      \
    {                                                                \
      ret (*pf)(params) = namespace##operator/< T_TypeA, T_TypeB >;  \
      return pf;                                                     \
    }                                                                \
  };

#define _FAKE_CLASS_OPV_2(ret, namespace, typeA, typeB, params, tag) \
  template< typeA T_TypeA, typeB T_TypeB >                           \
  class FAKE_CLASS_operatorV##_##tag                                 \
  {                                                                  \
    auto fake()                                                      \
    {                                                                \
      ret (*pf)(params) = namespace##operator%< T_TypeA, T_TypeB >;  \
      return pf;                                                     \
    }                                                                \
  };

#define INSTANTIATE_CLASS_FCT(ret, namespace, className, name, tag, sizeClassA, typeClassA, argsClassA, sizeA, typeA, argsA, params) \
  _FAKE_CLASS_1_1(TEMPLATE_RETURN(ret), namespace, className, name, typeClassA, typeA, TEMPLATE_PARAMS(params), tag)                 \
  _INSTANTIATE_CLASS_##sizeClassA##_##sizeA(FAKE_CLASS_##className##_##name##_##tag, TEMPLATE_ARGS(argsClassA), TEMPLATE_ARGS(argsA))
#define INSTANTIATE_CLASS_FCT_CONST(ret, namespace, className, name, tag, sizeClassA, typeClassA, argsClassA, sizeA, typeA, argsA, params) \
  _FAKE_CLASS_1_1_CONST(TEMPLATE_RETURN(ret), namespace, className, name, typeClassA, typeA, TEMPLATE_PARAMS(params), tag)                 \
  _INSTANTIATE_CLASS_##sizeClassA##_##sizeA(FAKE_CLASS_##className##_##name##_##tag, TEMPLATE_ARGS(argsClassA), TEMPLATE_ARGS(argsA))
#define INSTANTIATE_CLASS_OPLL(ret, namespace, className, tag, sizeClassA, typeClassA, argsClassA, sizeA, typeA, argsA, params) \
  _FAKE_CLASS_OPLL_1_1(TEMPLATE_RETURN(ret), namespace, className, typeClassA, typeA, TEMPLATE_PARAMS(params), tag)             \
  _INSTANTIATE_CLASS_##sizeClassA##_##sizeA(FAKE_CLASS_##className##_##operatorLL##_##tag, TEMPLATE_ARGS(argsClassA), TEMPLATE_ARGS(argsA))

#define INSTANTIATE_FCT_2(ret, namespace, name, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, params) \
  _FAKE_CLASS_2(TEMPLATE_RETURN(ret), namespace, name, typeA, typeB, TEMPLATE_PARAMS(params), tag)     \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB(FAKE_CLASS_##name##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define INSTANTIATE_OPP_2(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, params) \
  _FAKE_CLASS_OPP_2(TEMPLATE_RETURN(ret), namespace, typeA, typeB, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB(FAKE_CLASS_operatorP##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define INSTANTIATE_OPM_2(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, params) \
  _FAKE_CLASS_OPM_2(TEMPLATE_RETURN(ret), namespace, typeA, typeB, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB(FAKE_CLASS_operatorM##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define INSTANTIATE_OPT_2(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, params) \
  _FAKE_CLASS_OPT_2(TEMPLATE_RETURN(ret), namespace, typeA, typeB, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB(FAKE_CLASS_operatorT##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define INSTANTIATE_OPD_2(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, params) \
  _FAKE_CLASS_OPD_2(TEMPLATE_RETURN(ret), namespace, typeA, typeB, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB(FAKE_CLASS_operatorD##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))
#define INSTANTIATE_OPV_2(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, params) \
  _FAKE_CLASS_OPV_2(TEMPLATE_RETURN(ret), namespace, typeA, typeB, TEMPLATE_PARAMS(params), tag) \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB(FAKE_CLASS_operatorV##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB))

// 3 args
#define _FAKE_CLASS_3(ret, namespace, name, typeA, typeB, typeC, params, tag) \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC >                     \
  class FAKE_CLASS_##name##_##tag                                             \
  {                                                                           \
    auto fake()                                                               \
    {                                                                         \
      typedef ret (*Fct)(params);                                             \
      Fct pf = namespace##name< T_TypeA, T_TypeB, T_TypeC >;                  \
      return pf;                                                              \
    }                                                                         \
  };

#define _FAKE_CLASS_OPP_3(ret, namespace, typeA, typeB, typeC, params, tag)  \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC >                    \
  class FAKE_CLASS_operatorP##_##tag                                         \
  {                                                                          \
    auto fake()                                                              \
    {                                                                        \
      ret (*pf)(params) = namespace##operator+< T_TypeA, T_TypeB, T_TypeC >; \
      return pf;                                                             \
    }                                                                        \
  };

#define _FAKE_CLASS_OPM_3(ret, namespace, typeA, typeB, typeC, params, tag)  \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC >                    \
  class FAKE_CLASS_operatorM##_##tag                                         \
  {                                                                          \
    auto fake()                                                              \
    {                                                                        \
      ret (*pf)(params) = namespace##operator-< T_TypeA, T_TypeB, T_TypeC >; \
      return pf;                                                             \
    }                                                                        \
  };

#define _FAKE_CLASS_OPT_3(ret, namespace, typeA, typeB, typeC, params, tag)  \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC >                    \
  class FAKE_CLASS_operatorT##_##tag                                         \
  {                                                                          \
    auto fake()                                                              \
    {                                                                        \
      ret (*pf)(params) = namespace##operator*< T_TypeA, T_TypeB, T_TypeC >; \
      return pf;                                                             \
    }                                                                        \
  };

#define _FAKE_CLASS_OPD_3(ret, namespace, typeA, typeB, typeC, params, tag)  \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC >                    \
  class FAKE_CLASS_operatorD##_##tag                                         \
  {                                                                          \
    auto fake()                                                              \
    {                                                                        \
      ret (*pf)(params) = namespace##operator/< T_TypeA, T_TypeB, T_TypeC >; \
      return pf;                                                             \
    }                                                                        \
  };

#define _FAKE_CLASS_OPV_3(ret, namespace, typeA, typeB, typeC, params, tag)  \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC >                    \
  class FAKE_CLASS_operatorV##_##tag                                         \
  {                                                                          \
    auto fake()                                                              \
    {                                                                        \
      ret (*pf)(params) = namespace##operator%< T_TypeA, T_TypeB, T_TypeC >; \
      return pf;                                                             \
    }                                                                        \
  };

#define INSTANTIATE_FCT_3(ret, namespace, name, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, params) \
  _FAKE_CLASS_3(TEMPLATE_RETURN(ret), namespace, name, typeA, typeB, typeC, TEMPLATE_PARAMS(params), tag)                   \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC(FAKE_CLASS_##name##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define INSTANTIATE_OPP_3(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, params) \
  _FAKE_CLASS_OPP_3(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, TEMPLATE_PARAMS(params), tag)               \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC(FAKE_CLASS_operatorP##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define INSTANTIATE_OPM_3(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, params) \
  _FAKE_CLASS_OPM_3(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, TEMPLATE_PARAMS(params), tag)               \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC(FAKE_CLASS_operatorM##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define INSTANTIATE_OPT_3(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, params) \
  _FAKE_CLASS_OPT_3(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, TEMPLATE_PARAMS(params), tag)               \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC(FAKE_CLASS_operatorT##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define INSTANTIATE_OPD_3(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, params) \
  _FAKE_CLASS_OPD_3(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, TEMPLATE_PARAMS(params), tag)               \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC(FAKE_CLASS_operatorD##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))
#define INSTANTIATE_OPV_3(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, params) \
  _FAKE_CLASS_OPV_3(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, TEMPLATE_PARAMS(params), tag)               \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC(FAKE_CLASS_operatorV##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC))

// 4 args
#define _FAKE_CLASS_4(ret, namespace, name, typeA, typeB, typeC, typeD, params, tag) \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC, typeD T_TypeD >             \
  class FAKE_CLASS_##name##_##tag                                                    \
  {                                                                                  \
    auto fake()                                                                      \
    {                                                                                \
      typedef ret (*Fct)(params);                                                    \
      Fct pf = namespace##name< T_TypeA, T_TypeB, T_TypeC, T_TypeD >;                \
      return pf;                                                                     \
    }                                                                                \
  };

#define _FAKE_CLASS_OPP_4(ret, namespace, typeA, typeB, typeC, typeD, params, tag)    \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC, typeD T_TypeD >              \
  class FAKE_CLASS_operatorP##_##tag                                                  \
  {                                                                                   \
    auto fake()                                                                       \
    {                                                                                 \
      ret (*pf)(params) = namespace##operator+< T_TypeA, T_TypeB, T_TypeC, T_TypeD >; \
      return pf;                                                                      \
    }                                                                                 \
  };

#define _FAKE_CLASS_OPM_4(ret, namespace, typeA, typeB, typeC, typeD, params, tag)    \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC, typeD T_TypeD >              \
  class FAKE_CLASS_operatorM##_##tag                                                  \
  {                                                                                   \
    auto fake()                                                                       \
    {                                                                                 \
      ret (*pf)(params) = namespace##operator-< T_TypeA, T_TypeB, T_TypeC, T_TypeD >; \
      return pf;                                                                      \
    }                                                                                 \
  };

#define _FAKE_CLASS_OPT_4(ret, namespace, typeA, typeB, typeC, typeD, params, tag)    \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC, typeD T_TypeD >              \
  class FAKE_CLASS_operatorT##_##tag                                                  \
  {                                                                                   \
    auto fake()                                                                       \
    {                                                                                 \
      ret (*pf)(params) = namespace##operator*< T_TypeA, T_TypeB, T_TypeC, T_TypeD >; \
      return pf;                                                                      \
    }                                                                                 \
  };

#define _FAKE_CLASS_OPD_4(ret, namespace, typeA, typeB, typeC, typeD, params, tag)    \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC, typeD T_TypeD >              \
  class FAKE_CLASS_operatorD##_##tag                                                  \
  {                                                                                   \
    auto fake()                                                                       \
    {                                                                                 \
      ret (*pf)(params) = namespace##operator/< T_TypeA, T_TypeB, T_TypeC, T_TypeD >; \
      return pf;                                                                      \
    }                                                                                 \
  };

#define _FAKE_CLASS_OPV_4(ret, namespace, typeA, typeB, typeC, typeD, params, tag)    \
  template< typeA T_TypeA, typeB T_TypeB, typeC T_TypeC, typeD T_TypeD >              \
  class FAKE_CLASS_operatorV##_##tag                                                  \
  {                                                                                   \
    auto fake()                                                                       \
    {                                                                                 \
      ret (*pf)(params) = namespace##operator%< T_TypeA, T_TypeB, T_TypeC, T_TypeD >; \
      return pf;                                                                      \
    }                                                                                 \
  };

#define INSTANTIATE_FCT_4(ret, namespace, name, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, sizeD, typeD, argsD, params) \
  _FAKE_CLASS_4(TEMPLATE_RETURN(ret), namespace, name, typeA, typeB, typeC, typeD, TEMPLATE_PARAMS(params), tag)                                 \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD(FAKE_CLASS_##name##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define INSTANTIATE_OPP_4(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, sizeD, typeD, argsD, params) \
  _FAKE_CLASS_OPP_4(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, typeD, TEMPLATE_PARAMS(params), tag)                             \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD(FAKE_CLASS_operatorP##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define INSTANTIATE_OPM_4(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, sizeD, typeD, argsD, params) \
  _FAKE_CLASS_OPM_4(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, typeD, TEMPLATE_PARAMS(params), tag)                             \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD(FAKE_CLASS_operatorM##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define INSTANTIATE_OPT_4(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, sizeD, typeD, argsD, params) \
  _FAKE_CLASS_OPT_4(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, typeD, TEMPLATE_PARAMS(params), tag)                             \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD(FAKE_CLASS_operatorT##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define INSTANTIATE_OPD_4(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, sizeD, typeD, argsD, params) \
  _FAKE_CLASS_OPD_4(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, typeD, TEMPLATE_PARAMS(params), tag)                             \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD(FAKE_CLASS_operatorD##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))
#define INSTANTIATE_OPV_4(ret, namespace, tag, sizeA, typeA, argsA, sizeB, typeB, argsB, sizeC, typeC, argsC, sizeD, typeD, argsD, params) \
  _FAKE_CLASS_OPV_4(TEMPLATE_RETURN(ret), namespace, typeA, typeB, typeC, typeD, TEMPLATE_PARAMS(params), tag)                             \
  _INSTANTIATE_CLASS_##sizeA##_##sizeB##_##sizeC##_##sizeD(FAKE_CLASS_operatorV##_##tag, TEMPLATE_ARGS(argsA), TEMPLATE_ARGS(argsB), TEMPLATE_ARGS(argsC), TEMPLATE_ARGS(argsD))

#endif // H_GMSHFEM_INSTANTIATE
