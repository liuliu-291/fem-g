# GetDP - Copyright (C) 1997-2020 P. Dular and C. Geuzaine, University of Liege
#
# See the LICENSE.txt file for license information. Please report all
# issues on https://gitlab.onelab.info/getdp/getdp/issues.

set(SRC
  Bessel.cpp
)

if(ENABLE_BESSEL)
  list(APPEND SRC BesselLib.f)
endif(ENABLE_BESSEL)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
append_gmshfem_contrib(bessel "${SRC};${HDR}")
