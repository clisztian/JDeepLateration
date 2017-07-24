#include "StdAfx.h"
#include "LineHelper.h"


CLineHelper::CLineHelper(void)
{
	m_dwReadBufferAllocation		= 0;
	m_dwReadBufferSize				= 0;
	m_dwReadBufferWritePosition	= 0;
	m_dwReadBufferReadPosition		= 0;
	m_dwTempBufferSize				= 0;
}


CLineHelper::~CLineHelper(void)
{
	if (m_dwReadBufferAllocation)
		free (m_ReadBuffer);
	if (m_dwTempBufferSize)
		free (m_TempBuffer);
}

CString CLineHelper::ReadLine(CFile* pFile)
{
	// reads an entire line from the file

	CString	workString;
	DWORD		dwC;
	DWORD		dwS;
	DWORD		dwT;
	DWORD		dwTempSize;
	DWORD		dwStartPos;
	DWORD		dwEndPos;

	// initial read buffer handling
	if (m_dwReadBufferAllocation == 0) {
		m_dwReadBufferAllocation		= 262144;
		m_ReadBuffer						= (char*) calloc(m_dwReadBufferAllocation, 1);
		m_dwReadBufferReadPosition		= 0;
		m_dwReadBufferWritePosition	= 0;
		m_dwReadBufferSize				= 0;
		m_dwTempBufferSize				= 65536;
		m_TempBuffer						= (char*) calloc(m_dwTempBufferSize, 1);
	}
	// prepare read loop
	workString.Empty();
	dwC			= 0;
	dwStartPos	= m_dwReadBufferReadPosition;
	dwEndPos		= dwStartPos;
	dwTempSize	= m_dwReadBufferSize;
	// read loop
	while (1) {
		// check for available data
		if (dwTempSize == 0) {
			// read an additional block of information
			dwS	= __min(m_dwTempBufferSize, m_dwReadBufferAllocation - m_dwReadBufferSize);
			dwS	= pFile->Read(m_TempBuffer, dwS);
			if (dwS == 0)
				return workString;
			// copy to buffer
			if ((m_dwReadBufferWritePosition + dwS) > m_dwReadBufferAllocation) {
				dwT	= m_dwReadBufferAllocation - m_dwReadBufferWritePosition;
				memcpy (&m_ReadBuffer[m_dwReadBufferWritePosition], m_TempBuffer, dwT);
				memcpy (m_ReadBuffer, &m_TempBuffer[dwT], dwS - dwT);
			}
			else 
				memcpy (&m_ReadBuffer[m_dwReadBufferWritePosition], m_TempBuffer, dwS);
			// update write position
			m_dwReadBufferWritePosition	= (m_dwReadBufferWritePosition + dwS) % m_dwReadBufferAllocation;
			// update size in any case
			m_dwReadBufferSize	+= dwS;
			dwTempSize				+= dwS;
		}
		// search
		if (m_ReadBuffer[dwEndPos] == 10)
			break;
		// next element
		dwEndPos	= (dwEndPos + 1) % m_dwReadBufferAllocation;
		dwTempSize--;
	}
	// test version
	DWORD		dwTemp;
	if (dwEndPos > dwStartPos)
		dwTemp	= dwEndPos - dwStartPos;
	else
		dwTemp	= dwEndPos + m_dwReadBufferAllocation - dwStartPos;
	TRACE ("START: %u; STOP: %u; SIZE: %u\r\n", dwStartPos, dwEndPos, dwTemp);
	// copy operation
	memset (m_TempBuffer, 0, m_dwTempBufferSize);
	if (dwEndPos > dwStartPos) {
		dwS	= dwEndPos - dwStartPos - 1;	// skip final CR LF
		memcpy (m_TempBuffer, &m_ReadBuffer[dwStartPos], dwS);
	}
	else {
		dwS	= dwEndPos + m_dwReadBufferAllocation - dwStartPos - 1;	// skip final CR LF
		dwT	= m_dwReadBufferAllocation - dwStartPos;
		memcpy (m_TempBuffer, &m_ReadBuffer[dwStartPos], dwT);
		if (dwS > dwT)
			memcpy (&m_TempBuffer[dwT], m_ReadBuffer, dwS - dwT);
	}
	// adjust read buffer params
	m_dwReadBufferReadPosition	= (dwEndPos + 1) % m_dwReadBufferAllocation;
	m_dwReadBufferSize			-= dwS + 2;		// include CR LF
	// set the string
	workString.SetString(m_TempBuffer);
	// all done
	return workString;
}

int CLineHelper::ExtractNumber(CString workString)
{
	// extracts a single number from the input string
	int	i;
	int	lll;

	CString	numberString;
	char		ccc;

	// prepare loop
	lll	= workString.GetLength();
	i		= 0;
	numberString.Empty();

	// loop
	while (i<lll) {
		ccc	= workString[i];
		if ((ccc >= '0') && (ccc <= '9'))
			numberString	= numberString + ccc;
		else {
			if (!numberString.IsEmpty())
				i	= lll;	// to quit the loop
		}
		i++;
	}
	if (numberString.IsEmpty())
		return -1;
	else
		return atoi(numberString);
}

int CLineHelper::ExtractParameters(CString workString, DWORD* pTimestamp, int* pDaySeconds, CStringArray* pStringArray)
{
	// extracts all params separated by TAB from the input string
	//	returns the nb of paramaters found
	//	extracts vital timing information

	int		tempSize;
	CString	dayString;
	CString	strWork;
	char*		cWork;

	// reset information
	pStringArray->RemoveAll();
	// copy information since the string will suffer with strtok
	strWork	= workString;
	cWork		= strWork.GetBuffer();
	// prepare search loop
	cWork	= strtok(cWork, "\t");
	while (cWork) {
		pStringArray->Add(cWork);
		cWork	= strtok(NULL, "\t");
	}
	// size check
	tempSize = pStringArray->GetSize();
	if (tempSize < 4) {
		*pTimestamp		= 0;
		*pDaySeconds	= 0;
		return 0;	// invalid log entry
	}
	// extract the time stamp
	*pTimestamp	= atoi(pStringArray->GetAt(0));
	// extract day seconds
	dayString		= pStringArray->GetAt(2);	// skip the date
	*pDaySeconds	= 3600*atoi(dayString.Left(2)) + 60*atoi(dayString.Mid(3, 2)) + atoi(dayString.Right(2));
	// all done
	return tempSize;
}

CString CLineHelper::GetTimeString(int daySeconds)
{
	int	hh;
	int	mm;
	int	ss;

	CString	outString;

	hh	= daySeconds / 3600;
	mm	= (daySeconds - hh*3600) / 60;
	ss	= daySeconds % 60;

	outString.Format("%02d:%02d:%02d", hh, mm, ss);
	return outString;
}

CString CLineHelper::GetTimeString(double daySeconds)
{
	int	intSeconds;
	int	hh;
	int	mm;
	int	ss;

	double	rSecs;
	int		dSecs;
	CString	outString;

	intSeconds	= int(daySeconds + 1.0E-6);

	hh	= intSeconds/3600;
	mm	= (intSeconds - hh*3600) / 60;
	ss	= intSeconds % 60;
	// remove day breaks
	hh	%= 24;

	rSecs	= daySeconds - intSeconds;
	rSecs	= __max(0.0, rSecs);
	dSecs	= int(10.0*rSecs);

	outString.Format("%02d:%02d:%02d.%d", hh, mm, ss, dSecs);
	return outString;
}