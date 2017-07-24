#pragma once

struct NavigationInfo
{
	DWORD		dwTimestamp;
	int		nDaySeconds;
	double	dLongitude;
	double	dLatitude;
	double	dAltitude;
	DWORD		dwWECoor;
	DWORD		dwSNCoor;
};


class CNavigationHelper
{
public:
	CNavigationHelper(void);
	~CNavigationHelper(void);

	void	Reset();
	void	AddNavigationEntry(CStringArray* pItemArray);
	DWORD	Finalyse();

	DWORD	CollectNavigation(CString strFile);
	DWORD	CollectNavigation(CFile* pFile);
	void	BuildDistanceArray();
	int	GetNavigationAtTimestamp(DWORD dwTimestamp, int* pIndex, NavigationInfo* pInfo);
	int	GetNavigationAtDaySeconds(int nDaySeconds, int* pIndex, NavigationInfo* pInfo);
	int	GetDistanceAtTimestamp(DWORD dwTimestamp, int* pIndex, double* pDistance);
	void	Wgs2ChCoor(NavigationInfo* pInfo);
	void	ChCoor2Wgs(NavigationInfo* pInfo);
	int	GetDaySeconds(DWORD dwTimestamp);
	DWORD	GetTimestamp(int nDaySeconds);

public:
	DWORD					m_dwSize;
	NavigationInfo*	m_NavigationInfo;
	DWORD					m_dwFirstTimestamp;
	DWORD					m_dwLastTimestamp;
	int					m_nFirstDaySeconds;
	int					m_nLastDaySeconds;

	// distance handling
	DWORD					m_dwDistanceArrayAllocation;
	DWORD					m_dwDistanceArraySize;
	double*				m_DistanceArray;

	double				m_dSlope;
	double				m_dOffset;

private:
	DWORD					m_dwAllocation;

	double				m_dXSum;
	double				m_dYSum;
	double				m_dXYSum;
	double				m_dX2Sum;
};

