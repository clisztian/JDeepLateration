#include "StdAfx.h"
#include <math.h>
#include "NavigationHelper.h"
#include "LineHelper.h"

CNavigationHelper::CNavigationHelper(void)
{
	m_dwSize								= 0;
	m_dwAllocation						= 0;
	m_dwDistanceArrayAllocation	= 0;
	m_dwDistanceArraySize			= 0;
}


CNavigationHelper::~CNavigationHelper(void)
{
	if (m_dwAllocation)
		free (m_NavigationInfo);
	if (m_dwDistanceArrayAllocation)
		free (m_DistanceArray);
}

void CNavigationHelper::Reset()
{
	// reset all pending information
	m_dwSize						= 0;
	m_dwDistanceArraySize	= 0;
	m_nFirstDaySeconds		= 0;
	m_nLastDaySeconds			= 0;
	m_dwFirstTimestamp		= 0;
	m_dwLastTimestamp			= 0;
	m_dSlope						= 0.0;
	m_dOffset					= 0.0;
	m_dXSum						= 0.0;
	m_dYSum						= 0.0;
	m_dXYSum						= 0.0;
	m_dX2Sum						= 0.0;
}

DWORD CNavigationHelper::CollectNavigation(CString strFile)
{
	BOOL		ok;
	DWORD		dwSize;
	CFile		theFile;

	// try to open the file
	ok	= theFile.Open(strFile, CFile::modeRead, NULL);
	if (ok) {
		dwSize	= CollectNavigation(&theFile);
		theFile.Close();
	}
	else
		dwSize	= 0;
	// all done
	return dwSize;
}

void CNavigationHelper::AddNavigationEntry(CStringArray* pItemArray)
{
	// adds a single entry (string) from the log file
	CLineHelper		lineHelper;
	int				nbItem;
	int				iRef;
	CString			strWork;
	double			dXVal;
	double			dYVal;

	NavigationInfo	navInfo;
	// check the count
	nbItem	= pItemArray->GetSize();
	// analyse items
	if (nbItem == 7) {
		if (pItemArray->GetAt(3) == "NAVIGATION") {
			navInfo.dwTimestamp	= atoi(pItemArray->GetAt(0));
			navInfo.nDaySeconds	= lineHelper.ExtractNumber(pItemArray->GetAt(4));
			navInfo.dAltitude		= 500.0;
			strWork					= pItemArray->GetAt(5);
			iRef						= strWork.Find("=");
			strWork					= strWork.Mid(iRef+1);
			navInfo.dLongitude	= atof(strWork);
			strWork					= pItemArray->GetAt(6);
			iRef						= strWork.Find("=");
			strWork					= strWork.Mid(iRef+1);
			navInfo.dLatitude		= atof(strWork);
			Wgs2ChCoor(&navInfo);
			// add to container
			if (m_dwSize == m_dwAllocation) {
				if (m_dwAllocation) {
					m_dwAllocation	+= 32;
					m_NavigationInfo	= (NavigationInfo*) realloc(m_NavigationInfo, m_dwAllocation*sizeof(NavigationInfo));
				}
				else {
					m_dwAllocation		= 256;
					m_NavigationInfo	= (NavigationInfo*) calloc(m_dwAllocation, sizeof(NavigationInfo));
				}
			}
			m_NavigationInfo[m_dwSize]	= navInfo;
			m_dwSize++;
			// params for linear regression
			dXVal	= navInfo.dwTimestamp;
			dYVal	= navInfo.nDaySeconds;
			m_dXSum	+= dXVal;
			m_dYSum	+= dYVal;
			m_dXYSum	+= dXVal*dYVal;
			m_dX2Sum	+= dXVal*dXVal;
		}
	}
}

DWORD CNavigationHelper::CollectNavigation(CFile* pFile)
{
	// collects all navigation information from the file
	CLineHelper	lineHelper;
	CString		inpLine;
	CString		workString;

	CStringArray	itemArray;
	int				iref;
	int				nbItem;

	double			dXVal;
	double			dYVal;

	NavigationInfo	navInfo;
	// makes sure the file starts at the beginning
	pFile->SeekToBegin();
	// reset all pending information
	Reset();
	// collection loop
	while (1) {
		// read line
		inpLine	= lineHelper.ReadLine(pFile);
		if (inpLine.IsEmpty())
			break;
		// analyse items
		nbItem	= lineHelper.ExtractParameters(inpLine, &navInfo.dwTimestamp, &navInfo.nDaySeconds, &itemArray);
		if (nbItem == 7) {
			if (itemArray[3] == "NAVIGATION") {
				navInfo.dwTimestamp	= atoi(itemArray[0]);
				navInfo.nDaySeconds	= lineHelper.ExtractNumber(itemArray[4]);
				navInfo.dAltitude		= 500.0;
				workString				= itemArray[5];
				iref						= workString.Find("=");
				workString				= workString.Mid(iref+1);
				navInfo.dLongitude	= atof(workString);
				workString				= itemArray[6];
				iref						= workString.Find("=");
				workString				= workString.Mid(iref+1);
				navInfo.dLatitude		= atof(workString);
				Wgs2ChCoor(&navInfo);
				// add to container
				if (m_dwSize == m_dwAllocation) {
					if (m_dwAllocation) {
						m_dwAllocation	+= 32;
						m_NavigationInfo	= (NavigationInfo*) realloc(m_NavigationInfo, m_dwAllocation*sizeof(NavigationInfo));
					}
					else {
						m_dwAllocation		= 256;
						m_NavigationInfo	= (NavigationInfo*) calloc(m_dwAllocation, sizeof(NavigationInfo));
					}
				}
				m_NavigationInfo[m_dwSize]	= navInfo;
				m_dwSize++;
				// params for linear regression
				dXVal	= navInfo.dwTimestamp;
				dYVal	= navInfo.nDaySeconds;
				m_dXSum	+= dXVal;
				m_dYSum	+= dYVal;
				m_dXYSum	+= dXVal*dYVal;
				m_dX2Sum	+= dXVal*dXVal;
			}
		}
	}
	// keep limits
	if (m_dwSize) {
		m_nFirstDaySeconds	= m_NavigationInfo[0].nDaySeconds;
		m_nLastDaySeconds		= m_NavigationInfo[m_dwSize-1].nDaySeconds;
		m_dwFirstTimestamp	= m_NavigationInfo[0].dwTimestamp;
		m_dwLastTimestamp		= m_NavigationInfo[m_dwSize-1].dwTimestamp;
		// linear regression
		m_dSlope	= m_dwSize*m_dXYSum - m_dXSum*m_dYSum;
		m_dSlope	= m_dSlope / (m_dwSize*m_dX2Sum - m_dXSum*m_dXSum);
		m_dOffset	= (m_dYSum - m_dSlope*m_dXSum) / m_dwSize;
	}
	// all done
	return m_dwSize;
}

DWORD CNavigationHelper::Finalyse()
{
	// completes the collection process and prepares all params for analysis, returns the nb of valid navigation entries
	if (m_dwSize) {
		m_nFirstDaySeconds	= m_NavigationInfo[0].nDaySeconds;
		m_nLastDaySeconds		= m_NavigationInfo[m_dwSize-1].nDaySeconds;
		m_dwFirstTimestamp	= m_NavigationInfo[0].dwTimestamp;
		m_dwLastTimestamp		= m_NavigationInfo[m_dwSize-1].dwTimestamp;
		// linear regression
		m_dSlope	= m_dwSize*m_dXYSum - m_dXSum*m_dYSum;
		m_dSlope	= m_dSlope / (m_dwSize*m_dX2Sum - m_dXSum*m_dXSum);
		m_dOffset	= (m_dYSum - m_dSlope*m_dXSum) / m_dwSize;
	}
	// all done
	return m_dwSize;
}


void CNavigationHelper::Wgs2ChCoor(NavigationInfo* pInfo)
{
	// converts the WGS84 navigation info to swiss coor's

	double	lambda;
	double	phi;
	double	dVal;

	// normalise longitude and latitude
	phi		= (3600.0*pInfo->dLatitude - 169028.66) / 10000.0;
	lambda	= (3600.0*pInfo->dLongitude - 26782.5) / 10000.0;
	// calculation of EW coor
	dVal		=	600072.37 + 211455.93 * lambda;
	dVal		-=	10938.51 * lambda * phi;
	dVal		-= 0.36 * lambda * phi * phi;
	dVal		-= 44.54 * lambda * lambda * lambda;
	pInfo->dwWECoor	= DWORD(dVal + 0.5);
	// calculation of SN coor
	dVal		=	200147.07 + 308807.95 * phi;
	dVal		+=	3745.25 * lambda * lambda;
	dVal		+= 76.63 * phi * phi;
	dVal		-= 194.56 * lambda * lambda * phi;
	dVal		+= 119.79 * phi * phi * phi;
	pInfo->dwSNCoor	= DWORD(dVal + 0.5);
	// calculation of height
	dVal		= pInfo->dAltitude - 49.55;
	dVal		+= 2.73 * lambda;
	dVal		+= 6.94 * phi;
	pInfo->dAltitude	= dVal;
	// all done

}

void CNavigationHelper::ChCoor2Wgs(NavigationInfo* pInfo)
{
	// conversion from CH coor to WGS84
	double	x;
	double	y;
	double	dVal;

	x	= (pInfo->dwSNCoor - 200000.0) / 1000000.0;
	y	= (pInfo->dwWECoor - 600000.0) / 1000000.0;

	dVal	= 2.6779094 + 4.728982 * y;
	dVal	+= 0.791484 * y * x;
	dVal	+= 0.1306 * y * x * x;
	dVal	-= 0.0436 * y * y * y;
	pInfo->dLongitude	= dVal * 100.0 / 36.0;

	dVal	= 16.9023892 + 3.238272 * x;
	dVal	-= 0.270978 * y * y;
	dVal	-= 0.002528 * x * x;
	dVal	-= 0.0447 * y * y * x;
	dVal	-= 0.014 * x * x * x;
	pInfo->dLatitude	= dVal * 100.0 / 36.0;

	pInfo->dAltitude	= pInfo->dAltitude + 49.55 - 12.6*y - 22.64*x;
	// all done
}

int CNavigationHelper::GetNavigationAtTimestamp(DWORD dwTimestamp, int* pIndex, NavigationInfo* pInfo)
{
	// calculates interpolated navigation information for a given timestamp
	// returns 0 if the input value is within range
	DWORD		dwLo;
	DWORD		dwHi;
	int		nIndex;
	double	dScaler;
	// size check
	if (m_dwSize == 0) {
		memset (pInfo, 0, sizeof(NavigationInfo));
		return -1;
	}
	// low limit test
	if (dwTimestamp <= m_dwFirstTimestamp) {
		memcpy (pInfo, &m_NavigationInfo[0], sizeof(NavigationInfo));
		*pIndex	= 0;
		if (dwTimestamp == m_dwFirstTimestamp)
			return 0;
		else
			return -1;
	}
	// upper limit test
	if (dwTimestamp >= m_dwLastTimestamp) {
		memcpy (pInfo, &m_NavigationInfo[m_dwSize-1], sizeof(NavigationInfo));
		*pIndex	= m_dwSize - 1;
		if (dwTimestamp == m_dwLastTimestamp)
			return 0;
		else
			return -1;
	}
	// check for valid initial search index
	nIndex	= *pIndex;
	if ((nIndex < 0) || (DWORD(nIndex) > m_dwSize))
		nIndex	= 0;
	// search operation
	dwLo	= nIndex;
	dwHi	= __min(dwLo+1, m_dwSize-1);
	while ((dwTimestamp < m_NavigationInfo[dwLo].dwTimestamp) || (dwTimestamp > m_NavigationInfo[dwHi].dwTimestamp)) {
		if (dwTimestamp < m_NavigationInfo[dwLo].dwTimestamp) {
			dwLo--;
			dwHi--;
		}
		else {
			dwLo++;
			dwHi++;
			dwHi = __min(dwHi, m_dwSize + 1);
		}
	}
	// interpolation
	dScaler					= double(dwTimestamp - m_NavigationInfo[dwLo].dwTimestamp) / double(m_NavigationInfo[dwHi].dwTimestamp - m_NavigationInfo[dwLo].dwTimestamp);
	pInfo->dwTimestamp	= dwTimestamp;
	pInfo->dAltitude		= m_NavigationInfo[dwLo].dAltitude + dScaler*(m_NavigationInfo[dwHi].dAltitude - m_NavigationInfo[dwLo].dAltitude);
	pInfo->dLatitude		= m_NavigationInfo[dwLo].dLatitude + dScaler*(m_NavigationInfo[dwHi].dLatitude - m_NavigationInfo[dwLo].dLatitude);
	pInfo->dLongitude		= m_NavigationInfo[dwLo].dLongitude + dScaler*(m_NavigationInfo[dwHi].dLongitude - m_NavigationInfo[dwLo].dLongitude);
	pInfo->nDaySeconds	= int(m_NavigationInfo[dwLo].nDaySeconds + dScaler*(m_NavigationInfo[dwHi].nDaySeconds - m_NavigationInfo[dwLo].nDaySeconds) + 0.5);
	// get CH Coor's
	Wgs2ChCoor(pInfo);
	// keep index
	*pIndex	= dwLo;
	// all done
	return 0;
}

int CNavigationHelper::GetNavigationAtDaySeconds(int nDaySeconds, int* pIndex, NavigationInfo* pInfo)
{
	// extracts the navigation to specific day seconds
	//	the return value is 0, if the input is within the time span
	DWORD		dwLo;
	DWORD		dwHi;
	int		nIndex;
	int		nSize;
	double	scaler;
	// size check
	if (m_dwSize == 0) {
		memset (pInfo, 0, sizeof(NavigationInfo));
		*pIndex	= 0;
		return -1;
	}
	// lo limit test
	if (nDaySeconds <= m_nFirstDaySeconds) {
		memcpy (pInfo, &m_NavigationInfo[0], sizeof(NavigationInfo));
		*pIndex	= 0;
		if (nDaySeconds == m_nFirstDaySeconds)
			return 0;
		else
			return -1;
	}
	// hi limit test
	if (nDaySeconds >= m_nLastDaySeconds) {
		memcpy (pInfo, &m_NavigationInfo[m_dwSize-1], sizeof(NavigationInfo));
		*pIndex	= m_dwSize - 1;
		if (nDaySeconds == m_nLastDaySeconds)
			return 0;
		else
			return -1;
	}
	// interpolation section
	nIndex	= *pIndex;
	nSize		= m_dwSize;
	if ((nIndex < 0) || (nIndex >= nSize)) {
		dwLo	= 0;
		dwHi	= __min(1, nSize-1);
	}
	else {
		dwLo	= __min(nIndex, nSize - 2);
		dwHi	= dwLo + 1;
	}
	while ((nDaySeconds < m_NavigationInfo[dwLo].nDaySeconds) || (nDaySeconds > m_NavigationInfo[dwHi].nDaySeconds)) {
		if (nDaySeconds < m_NavigationInfo[dwLo].nDaySeconds) {
			dwLo--;
			dwHi--;
		}
		else {
			dwLo++;
			dwHi	= __min(dwHi+1, m_dwSize-1);
		}
	}
	scaler = double(nDaySeconds - m_NavigationInfo[dwLo].nDaySeconds) / double(m_NavigationInfo[dwHi].nDaySeconds - m_NavigationInfo[dwLo].nDaySeconds);
	pInfo->dAltitude	= m_NavigationInfo[dwLo].dAltitude + scaler * (m_NavigationInfo[dwHi].dAltitude - m_NavigationInfo[dwLo].dAltitude);
	pInfo->dLatitude	= m_NavigationInfo[dwLo].dLatitude + scaler * (m_NavigationInfo[dwHi].dLatitude - m_NavigationInfo[dwLo].dLatitude);
	pInfo->dLongitude	= m_NavigationInfo[dwLo].dLongitude + scaler * (m_NavigationInfo[dwHi].dLongitude - m_NavigationInfo[dwLo].dLongitude);
	pInfo->nDaySeconds	= nDaySeconds;
	// get CH corr's
	Wgs2ChCoor(pInfo);
	// keep index
	*pIndex	= dwLo;
	// all done
	return 0;
}

int CNavigationHelper::GetDaySeconds(DWORD dwTimestamp)
{
	double	dVal;

	dVal	= m_dSlope*dwTimestamp + m_dOffset;
	return int(dVal + 0.5);
}

DWORD CNavigationHelper::GetTimestamp(int nDaySeconds)
{
	double	dVal;

	dVal	= (nDaySeconds - m_dOffset) / m_dSlope;
	return DWORD(dVal + 0.5);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// distance handling

void CNavigationHelper::BuildDistanceArray()
{
	// creates an array with the distance from start for the given navigation sequence
	DWORD		dw;
	int		nXVal1;
	int		nYVal1;
	int		nXVal2;
	int		nYVal2;
	double	dDist;
	double	dSum;

	// memory handling
	if (m_dwDistanceArrayAllocation < m_dwSize) {
		if (m_dwDistanceArrayAllocation)
			free (m_DistanceArray);
		m_dwDistanceArrayAllocation	= m_dwSize;
		m_DistanceArray					= (double*) calloc(m_dwDistanceArrayAllocation, sizeof(double));
	}
	// prepare the loop
	dDist						= 0.0;
	m_DistanceArray[0]	= 0.0;
	nXVal1					= m_NavigationInfo[0].dwWECoor;
	nYVal1					= m_NavigationInfo[0].dwSNCoor;
	// collection loop
	for (dw=1; dw<m_dwSize; dw++) {
		nXVal2	= m_NavigationInfo[dw].dwWECoor;
		nYVal2	= m_NavigationInfo[dw].dwSNCoor;
		dSum		= (nXVal2 - nXVal1) * (nXVal2 - nXVal1);
		dSum		+= (nYVal2 - nYVal1) * (nYVal2 - nYVal1);

		dDist						+= sqrt(dSum);
		m_DistanceArray[dw]	= dDist;
		// keep values
		nXVal1	= nXVal2;
		nYVal1	= nYVal2;
	}
	// keep size
	m_dwDistanceArraySize	= m_dwSize;
	// all done
}

int CNavigationHelper::GetDistanceAtTimestamp(DWORD dwTimestamp, int* pIndex, double* pDistance)
{
	// calculates the distance from start for a given timestamp, returns 0 if the information is valid
	DWORD		dwLo;
	DWORD		dwHi;
	int		nIndex;
	double	dScaler;
	// size check
	if (m_dwDistanceArraySize == 0) {
		*pDistance	= 0.0;
		return -1;
	}
	// low limit test
	if (dwTimestamp <= m_dwFirstTimestamp) {
		*pDistance	= m_DistanceArray[0];
		*pIndex		= 0;
		if (dwTimestamp == m_dwFirstTimestamp)
			return 0;
		else
			return -1;
	}
	// upper limit test
	if (dwTimestamp >= m_dwLastTimestamp) {
		*pDistance	= m_DistanceArray[m_dwDistanceArraySize-1];
		*pIndex	= m_dwSize - 1;
		if (dwTimestamp == m_dwLastTimestamp)
			return 0;
		else
			return -1;
	}
	// check for valid initial search index
	nIndex	= *pIndex;
	if ((nIndex < 0) || (DWORD(nIndex) > m_dwSize))
		nIndex	= 0;
	// search operation
	dwLo	= nIndex;
	dwHi	= __min(dwLo+1, m_dwSize-1);
	while ((dwTimestamp < m_NavigationInfo[dwLo].dwTimestamp) || (dwTimestamp > m_NavigationInfo[dwHi].dwTimestamp)) {
		if (dwTimestamp < m_NavigationInfo[dwLo].dwTimestamp) {
			dwLo--;
			dwHi--;
		}
		else {
			dwLo++;
			dwHi++;
			dwHi = __min(dwHi, m_dwSize + 1);
		}
	}
	// interpolation
	dScaler			= double(dwTimestamp - m_NavigationInfo[dwLo].dwTimestamp) / double(m_NavigationInfo[dwHi].dwTimestamp - m_NavigationInfo[dwLo].dwTimestamp);
	*pDistance		= m_DistanceArray[dwLo] + dScaler*(m_DistanceArray[dwHi] - m_DistanceArray[dwLo]);
	// keep index
	*pIndex	= dwLo;
	// all done
	return 0;
}
