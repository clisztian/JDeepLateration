#include "StdAfx.h"
#include <math.h>
#include "CIRHelper.h"
#include "LineHelper.h"


CCIRHelper::CCIRHelper(void)
{
}


CCIRHelper::~CCIRHelper(void)
{
}

int CCIRHelper::ParseCIRString(char* cStr, double* cirArray)
{
	// extracts the CIR information out of the string representation, returns the size of the detected CIR
	int		icc;
	// prepare loop
	icc		= 0;
	cStr		= strtok(cStr, "|");
	// loop
	while (cStr) {
		cirArray[icc++]	= atof(cStr);
		cStr	= strtok(NULL, "|");
	}
	// all done
	return icc;
}

int CCIRHelper::GetCIRRange(CFile* pFile, double* pMinTime, double* pMaxTime, double* pResol)
{
	CLineHelper		lineHelper;
	int				iRef;
	CString			strLine;
	CString			strWork;
	CStringArray	itemArray;
	DWORD				dwTimestamp;
	int				nDaySeconds;
	int				nItem;
	double			dBestResol;
	double			dTestResol;
	double			dBestError;
	double			dError;
	double			dCIRStartTime;
	double			dCIREndTime;
	double			dCIRResolution;

	// read until "RANGE" has been found
	while (1) {
		strLine	= lineHelper.ReadLine(pFile);
		if (strLine.IsEmpty())
			return -1;
		nItem	= lineHelper.ExtractParameters(strLine, &dwTimestamp, &nDaySeconds, &itemArray);
		if (nItem) {
			if (itemArray[3] == "RANGE")
				break;
		}
	}
	// get the params from the range
	strWork				= itemArray[4];		// start of range
	iRef					= strWork.Find("=");
	dCIRStartTime		= atof(strWork.Mid(iRef+1));
	strWork				= itemArray[5];		// end of range
	iRef					= strWork.Find("=");
	dCIREndTime			= atof(strWork.Mid(iRef+1));
	strWork				= itemArray[6];		// resolution
	iRef					= strWork.Find("=");
	dCIRResolution		= atof(strWork.Mid(iRef+1));
	// make sure the CIR resolution fits one of the sampling rates
	dBestError	= 1.0E99;
	dTestResol	= 1.0;
	while (dTestResol > 0.06) {
		dError	= fabs(dTestResol - dCIRResolution);
		if (dError < dBestError) {
			dBestError	= dError;
			dBestResol	= dTestResol;
		}
		dTestResol	*= 0.5;
	}
	dCIRResolution	= dBestResol;
	// keep params
	*pMinTime	= dCIRStartTime;
	*pMaxTime	= dCIREndTime;
	*pResol		= dCIRResolution;
	// all done
	return 0;
}