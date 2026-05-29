#include "ParameterTools.hh"

using namespace ParameterTools;


std::vector<std::string> ParameterTools::splitOnSpace(const char* cstr)
{

	std::stringstream stream;
	stream << cstr;
	std::string token;
	std::vector<std::string> vec;
	while(std::getline(stream, token, ' '))
	{
		if(token.length() != 0)
		{
			CCTK_VInfo(CCTK_THORNSTRING, "splitting %s", token.c_str());
			vec.push_back(token);
		}
	}

	return vec;
}

CCTK_INT ParameterTools::findArrayIndex(const char* givenvars, const char* findval)
{
	std::vector<std::string> givenvec = splitOnSpace(givenvars);
	return findArrayIndexFromString(givenvec, findval);
}
CCTK_INT ParameterTools::findArrayIndexFromString(std::vector<std::string> givenvec, const char* findval)
{
	std::string target_lower(findval);
	std::transform(target_lower.begin(), target_lower.end(), target_lower.begin(), ::tolower);
	for(std::vector<std::string>::iterator it = givenvec.begin(); it != givenvec.end(); ++it)
	{
		CCTK_INT dist = std::distance(givenvec.begin(), it);
		std::string this_lower(*it);
		std::transform(this_lower.begin(), this_lower.end(), this_lower.begin(), ::tolower);
		CCTK_VInfo(CCTK_THORNSTRING, "comparing %s to %s", it->c_str(), target_lower.c_str());
		if(target_lower.compare(this_lower) == 0)
		{
			CCTK_VInfo(CCTK_THORNSTRING, "match at %d", dist);
			return dist;
		}
		else
		{
			CCTK_VInfo(CCTK_THORNSTRING, "no match between %s and %s", this_lower.c_str(), target_lower.c_str());
		}
	}

	//if no matches
	CCTK_VInfo(CCTK_THORNSTRING, "could not match %s", target_lower.c_str());
	return -1;
}

std::vector<CCTK_INT> ParameterTools::findArrayIndices(const char* givenvars, const char* findvals)
{
	std::vector<std::string> findvec = splitOnSpace(findvals);
	std::vector<std::string> givenvec = splitOnSpace(givenvars);
	std::vector<CCTK_INT> retval;
	for(std::vector<std::string>::iterator it = findvec.begin(); it != findvec.end(); ++it)
	{
		retval.push_back(findArrayIndexFromString(givenvec, it->c_str()));
	}
	return retval;

}

const char* ParameterTools::stitchGroupStringWithIndex(const char* groupname, CCTK_INT index)
{
	std::stringstream stream;
	stream << groupname << "[" << index << "]";
	return stream.str().c_str();
}
	
/*

CCTK_INT ParameterTools::appendToStringParameter(const char* paramname, const char* thorn, const char* appends)
{
	std::stringstream stream;
	char* valstring = CCTK_ParameterValString(paramname, thorn);
	//no leading space is necessary
	if(CCTK_Equals(valstring, ""))
	{
		stream << appends;
	}
	else
	{
		stream << valstring << " " << appends;
	}

	free(valstring);
	CCTK_VInfo(CCTK_THORNSTRING, "Appending %s to parameter %s::%s", appends, thorn, paramname);
	CCTK_INT ierr = CCTK_ParameterSet(paramname, thorn, stream.to_str().c_str());
	return ierr;
}
*/

CCTK_INT ParameterTools::setRealParameter(const char* paramname, const char* thorn, const CCTK_REAL value)
{
	std::stringstream stream;
	stream << value;
	return CCTK_ParameterSet(paramname, thorn, stream.str().c_str());
}
