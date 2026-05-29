#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include "cctk.h"

#ifndef PPT_HFILE
#define PPT_HFILE

namespace ParameterTools
{
	CCTK_INT findArrayIndex(const char* givenvars, const char* findval);
	std::vector<CCTK_INT> findArrayIndices(const char* givenvars, const char* findvals);
	CCTK_INT findArrayIndexFromString(std::vector<std::string> givenvec, const char* findval);
	std::vector<std::string> splitOnSpace(const char* cstr);
	const char* stitchGroupStringWithIndex(const char* groupname, CCTK_INT index);
	CCTK_INT setRealParameter(const char* paramname, const char* thorn, const CCTK_REAL value);
	//CCTK_INT appendToStringParameter(const char* paramname, const char* thorn, const char* appends);
}

#endif
