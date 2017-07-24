#pragma once

struct SPeakEntry {
	int		nIndex;
	DWORD		dwTimestamp;
	double	dAbsTime;
	double	dPeakTime;
	double	dLevel;
	double	dRSSI;
	double	dSNR;
	double	dFreqError;
};

struct SCIREntry {
	SPeakEntry	sPeakEntry;
	ULONGLONG	ulFilePosition;
	int			nStringLen;
	double		dCRT;
};


class CCIRHelper
{
public:
	CCIRHelper(void);
	~CCIRHelper(void);

	int	ParseCIRString(char* cStr, double* cirArray);
	int	GetCIRRange(CFile* pFile, double* pMinTime, double* pMaxTime, double* pResol);


};

