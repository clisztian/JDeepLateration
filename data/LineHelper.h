#pragma once
class CLineHelper
{
public:
	CLineHelper(void);
	~CLineHelper(void);

public:
	CString	ReadLine(CFile* pFile);
	int		ExtractParameters(CString workString, DWORD* pTimestamp, int* pDaySeconds, CStringArray* pStringArray);
	int		ExtractNumber(CString workString);
	CString	GetTimeString(int daySeconds);
	CString	GetTimeString(double daySeconds);

private:
	DWORD		m_dwReadBufferAllocation;
	DWORD		m_dwReadBufferWritePosition;
	DWORD		m_dwReadBufferReadPosition;
	DWORD		m_dwReadBufferSize;
	DWORD		m_dwTempBufferSize;

	char*		m_ReadBuffer;
	char*		m_TempBuffer;
};

