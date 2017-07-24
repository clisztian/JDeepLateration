// CIRDisplayDlg.cpp: Implementierungsdatei
//

#include "stdafx.h"
#include <math.h>
#include "UniversalChannelMeas.h"
#include "CIRDisplayDlg.h"
#include "ColorMapper.h"
#include "afxdialogex.h"
#include "FileSelectionHelper.h"
#include "LineHelper.h"
#include "DopplerDisplayDlg.h"
#include "SimulationDlg.h"


// CCIRDisplayDlg-Dialogfeld

IMPLEMENT_DYNAMIC(CCIRDisplayDlg, CDialog)

CCIRDisplayDlg::CCIRDisplayDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CCIRDisplayDlg::IDD, pParent)
	, m_StrINPFILE(_T(""))
	, m_StrMINMEASTIME(_T(""))
	, m_StrMAXMEASTIME(_T(""))
	, m_StrCYCLEPERIOD(_T(""))
	, m_StrCYCLEOFFSET(_T(""))
	, m_StrSNRTHRESHOLD(_T(""))
	, m_StrMINPLACEMENT(_T(""))
	, m_StrMAXPLACEMENT(_T(""))
	, m_StrMINLEVEL(_T(""))
	, m_StrMAXLEVEL(_T(""))
	, m_CheckMAINPEAK(FALSE)
	, m_CheckSTARTSTOP(FALSE)
	, m_StrFREQINDEX(_T(""))
	, m_StrCAPTURERANGE(_T(""))
	, m_CheckMULTIFREQ(FALSE)
	, m_CheckDISTVIEW(FALSE)
{
	// prepare a peak entry array
	m_dwEntryAllocation	= 2048;
	m_dwEntrySize			= 0;
	m_PeakEntries			= (SPeakEntry*) calloc(m_dwEntryAllocation, sizeof(SPeakEntry));
	m_CIREntries			= (SCIREntry*) calloc(m_dwEntryAllocation, sizeof(SCIREntry));
	// fill the color table
	m_ColorTable[0]	= RGB(1, 1, 1);
	m_ColorTable[1]	= RGB(255, 0, 0);
	m_ColorTable[2]	= RGB(0, 255, 0);
	m_ColorTable[3]	= RGB(0, 0, 255);
	m_ColorTable[4]	= RGB(255, 255, 0);
	m_ColorTable[5]	= RGB(0, 255, 255);
	m_ColorTable[6]	= RGB(255, 128, 0);
	m_ColorTable[7]	= RGB(255, 0, 128);
}

CCIRDisplayDlg::~CCIRDisplayDlg()
{
	if (m_dwEntryAllocation) {
		free (m_PeakEntries);
		free (m_CIREntries);
	}
	if (m_pSingleCIRDlg)
		delete(m_pSingleCIRDlg);
}

void CCIRDisplayDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_INPFILE, m_StrINPFILE);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_MINTIME, m_StrMINMEASTIME);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_MAXTIME, m_StrMAXMEASTIME);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_CYCLEPERIOD, m_StrCYCLEPERIOD);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_CYCLEOFFSET, m_StrCYCLEOFFSET);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_SNRTHRESHOLD, m_StrSNRTHRESHOLD);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_MINPLACEMENT, m_StrMINPLACEMENT);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_MAXPLACEMENT, m_StrMAXPLACEMENT);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_MINLEVEL, m_StrMINLEVEL);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_MAXLEVEL, m_StrMAXLEVEL);
	DDX_Control(pDX, IDC_PEAKDISPLAY_DLG_GRAP_PEAKDISPLAY, m_GrapPEAKDISPLAY);
	DDX_Control(pDX, IDC_PEAKDISPLAY_DLG_GRAP_SCALING, m_GrapSCALING);
	DDX_Check(pDX, IDC_CIRDISPLAY_DLG_CHECK_MAINPEAK, m_CheckMAINPEAK);
	DDX_Check(pDX, IDC_CIRDISPLAY_DLG_CHECK_STARTSTOPCORR, m_CheckSTARTSTOP);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_FREQINDEX, m_StrFREQINDEX);
	DDX_Text(pDX, IDC_CIRDISPLAY_DLG_EDIT_CAPTURE, m_StrCAPTURERANGE);
	DDX_Check(pDX, IDC_CIRDISPLAY_DLG_CHECK_MULTIFREQ, m_CheckMULTIFREQ);
	DDX_Check(pDX, IDC_CIRDISPLAY_DLG_CHECK_DISTVIEW, m_CheckDISTVIEW);
}


BEGIN_MESSAGE_MAP(CCIRDisplayDlg, CDialog)
	ON_BN_CLICKED(IDC_CIRDISPLAY_DLG_BTN_INPFILE, &CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnInpfile)
	ON_BN_CLICKED(IDC_CIRDISPLAY_DLG_BTN_UPDATE, &CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnUpdate)
	ON_BN_CLICKED(IDC_CIRDISPLAY_DLG_BTN_CRTVIEW, &CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnCrtview)
	ON_BN_CLICKED(IDC_CIRDISPLAY_DLG_BTN_DOPPLERVIEW, &CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnDopplerview)
	ON_STN_DBLCLK(IDC_PEAKDISPLAY_DLG_GRAP_PEAKDISPLAY, &CCIRDisplayDlg::OnStnDblclickPeakdisplayDlgGrapPeakdisplay)
	ON_MESSAGE(WM_USER_FLOATING_DLG_CLOSE, OnFloatingDlgClose)
	ON_MESSAGE(WM_USER_FLOATING_DLG_STEP, OnFloatingDlgStep)
	ON_BN_CLICKED(IDC_CIRDISPLAY_DLG_BTN_SIMULATION, &CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnSimulation)
END_MESSAGE_MAP()


// CCIRDisplayDlg-Meldungshandler

BOOL CCIRDisplayDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  hier zusätzliche Initialisierung hinzufügen.
	m_dwEntrySize	= 0;
	m_GrapSCALING.InitGraph();
	m_pSingleCIRDlg	= NULL;
	CheckRadioButton(IDC_CIRDISPLAY_DLG_RADIO_ABSOLUTE_PEAKLEVEL, IDC_CIRDISPLAY_DLG_RADIO_RELATIVE_PEAKLEVEL, IDC_CIRDISPLAY_DLG_RADIO_ABSOLUTE_PEAKLEVEL);
	UpdateData(FALSE);

	return TRUE;  // return TRUE unless you set the focus to a control
	// AUSNAHME: OCX-Eigenschaftenseite muss FALSE zurückgeben.
}



void CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnInpfile()
{
	// TODO: Fügen Sie hier Ihren Kontrollbehandlungscode für die Benachrichtigung ein.
	BOOL						ok;
	CFileSelectionHelper	helper;
	CString					tempName;
	CString					tempDir;
	DWORD						dwC;

	// get the new file
	ok	= helper.SelectOpenFileName("", "Log-Files","log", &tempName, &tempDir);
	if (ok) {
		UpdateData(TRUE);
		m_StrINPFILE	= tempName;
		UpdateData(FALSE);
		// parse the file
		m_dwEntrySize	= ParseLogFile();
		// match corresponding CIR's
		dwC	= ParseCIRFile();
	}
}

void CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnCrtview()
{
	// TODO: Fügen Sie hier Ihren Kontrollbehandlungscode für die Benachrichtigung ein.
	CCRTDisplayDlg	dlg;

	// copy existing params to the dlg
	dlg.m_pNavHelper			= &m_NavHelper;
	dlg.m_dwDataSize			= m_dwCIRCount;
	dlg.m_CIRArray				= m_CIREntries;
	dlg.m_strInpFile			= m_StrINPFILE;
	// launch the dlg
	dlg.DoModal();
}

void CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnDopplerview()
{
	// TODO: Fügen Sie hier Ihren Kontrollbehandlungscode für die Benachrichtigung ein.
	CDopplerDisplayDlg	dlg;

	// copy existing params to the dlg
	dlg.m_pNavHelper			= &m_NavHelper;
	dlg.m_dwDataSize			= m_dwEntrySize;
	dlg.m_PeakEntryArray		= m_PeakEntries;
	// launch the dlg
	dlg.DoModal();
}

void CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnUpdate()
{
	// TODO: Fügen Sie hier Ihren Kontrollbehandlungscode für die Benachrichtigung ein.

	int				nRet;
	int				nNavIndex;
	DWORD				dw;
	DWORD				dwC;
	DWORD				displayMode;
	DWORD				graphicsMode;
	CColorMapper	colorMapper;

	int		offset[4];

	double*		xArray;
	double*		yArray;
	COLORREF*	cArray;

	BOOL		bRelMode;

	double	dMinTime;
	double	dMaxTime;
	double	dMinLevel;
	double	dMaxLevel;
	double	dTime;
	double	dSlope;
	double	dPeriod;
	double	dOffset;
	double	dMinSNR;
	double	dRecentTime;
	double	dRefValue;
	double	dCIRRange;
	double	dLimit;
	double	dDist;
	int		nIndex;

	// check data availability
	if (m_dwEntrySize < 1) {
		AfxMessageBox("No valid data available");
		return;
	}
	// get the current settings
	UpdateData(TRUE);
	dMinTime				= atof(m_StrMINMEASTIME);
	dMaxTime				= atof(m_StrMAXMEASTIME);
	dMinLevel			= atof(m_StrMINLEVEL);
	dMaxLevel			= atof(m_StrMAXLEVEL);
	dPeriod				= atof(m_StrCYCLEPERIOD)/1.0E6;		// calculation in sec, display in usec's
	dOffset				= atof(m_StrCYCLEOFFSET)/1.0E6;		// calculation in sec, dsiplay in usec's
	dMinSNR				= atof(m_StrSNRTHRESHOLD);
	if (m_StrFREQINDEX.IsEmpty())
		nIndex	= -1;
	else
		nIndex	= atoi(m_StrFREQINDEX);
	if (m_CheckSTARTSTOP)
		dSlope	= m_dStartStopSlope;
	else
		dSlope	= 0.0;
	if (m_StrCAPTURERANGE.IsEmpty())
		dCIRRange	= 1000.0;
	else
		dCIRRange	= atof(m_StrCAPTURERANGE);
	// check for valid distance information
	if (m_NavHelper.m_dwDistanceArraySize < 2)
		m_CheckDISTVIEW	= FALSE;
	// prepare arrays to hold the entire data
	xArray	= (double*) calloc(m_dwEntrySize, sizeof(double));
	yArray	= (double*) calloc(m_dwEntrySize, sizeof(double));
	cArray	= (COLORREF*) calloc(m_dwEntrySize, sizeof(COLORREF));
	// prepare color mapping
	colorMapper.SetRange(dMinLevel, dMaxLevel);
	// collection loop
	dwC					= 0;
	dRecentTime			= 0.0;
	dRefValue			= 0.0;
	nNavIndex			= -1;
	bRelMode		= GetCheckedRadioButton(IDC_CIRDISPLAY_DLG_RADIO_ABSOLUTE_PEAKLEVEL, IDC_CIRDISPLAY_DLG_RADIO_RELATIVE_PEAKLEVEL) == IDC_CIRDISPLAY_DLG_RADIO_RELATIVE_PEAKLEVEL;
	for (dw=0; dw<m_dwEntrySize; dw++) {
		dTime	= m_PeakEntries[dw].dAbsTime;
		// check frequency index (-1 indicates all frequencies)
		if (nIndex != -1) {
			if (m_PeakEntries[dw].nIndex != nIndex)
				continue;
		}
		// keep the reference value for the first item
		if (dTime != dRecentTime) {
			// track first elements
			dRecentTime	= dTime;
			dLimit		= m_PeakEntries[dw].dLevel - dCIRRange;
			// keep the distance from start, if requested
			if (m_CheckDISTVIEW)
				nRet	= m_NavHelper.GetDistanceAtTimestamp(m_PeakEntries[dw].dwTimestamp, &nNavIndex, &dDist);
			// keep power of first element as reference
			if (bRelMode)
				dRefValue	= m_PeakEntries[dw].dLevel;
		}
		else {
			if (m_CheckMAINPEAK)
				continue;
		}
		// additional checks
		if ((m_PeakEntries[dw].dSNR > dMinSNR) && (m_PeakEntries[dw].dLevel >= dLimit)) {
			if (m_CheckDISTVIEW)
				yArray[dwC]	= dDist/1000.0;		// display in km's
			else
				yArray[dwC]	= dTime;
			// calculate peak time
			dTime	= m_PeakEntries[dw].dPeakTime;
			// apply offset
			dTime	+= dOffset;
			// drift compensation
			dTime	-= dSlope * (m_PeakEntries[dw].dAbsTime - m_PeakEntries[0].dAbsTime);
			// remove period
			dTime	= fmod (dTime, dPeriod);
			// keep value
			xArray[dwC]	= 1.0E6 * dTime;		// calculation in sec, display in usec's
			// multi color handling
			if (m_CheckMULTIFREQ)
				cArray[dwC]	= m_ColorTable[m_PeakEntries[dw].nIndex];
			else
				cArray[dwC]	= colorMapper.GetColor(m_PeakEntries[dw].dLevel - dRefValue);
			dwC++;
		}
	}
	// prepare graphics
	graphicsMode	= GRAPHICS_MODE_POINT | GRAPHICS_MODE_XY;
	displayMode		= DISPLAY_MODE_LABELS;
	offset[0]	= 10;
	offset[1]	= 10;
	offset[2]	= 8;
	offset[3]	= 4;
	m_GrapPEAKDISPLAY.InitGraph(dwC, 0, graphicsMode, displayMode, offset);
	if (m_CheckDISTVIEW)
		m_GrapPEAKDISPLAY.SetLabels("Cycle placement [usec]", "Distance [km]");
	else
		m_GrapPEAKDISPLAY.SetLabels("Cycle placement [usec]", "Test time [sec]");
	m_GrapPEAKDISPLAY.xyl[0]	= atof(m_StrMINPLACEMENT);
	m_GrapPEAKDISPLAY.xyl[2]	= atof(m_StrMAXPLACEMENT);
	m_GrapPEAKDISPLAY.xyl[1]	= dMinTime;
	m_GrapPEAKDISPLAY.xyl[3]	= dMaxTime;
	m_GrapPEAKDISPLAY.SetMaxNbOfAnnotations(12, 8);
	// finally set the points
	m_GrapPEAKDISPLAY.SetData(dwC, xArray, yArray, cArray);
	m_GrapPEAKDISPLAY.OnPaint();
	// update scaling display
	m_GrapSCALING.Update(dMinLevel, dMaxLevel);
	// all done
	free (xArray);
	free (yArray);
	free (cArray);
}

void CCIRDisplayDlg::OnBnClickedCirdisplayDlgBtnSimulation()
{
	// TODO: Fügen Sie hier Ihren Kontrollbehandlungscode für die Benachrichtigung ein.
	CSimulationDlg	dlg;

	// set the params
	dlg.m_dwEntrySize	= m_dwEntrySize;
	dlg.m_PeakEntries	= m_PeakEntries;
	// show the dlg
	dlg.DoModal();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// floating dlg handling

LRESULT CCIRDisplayDlg::OnFloatingDlgClose(WPARAM wParam, LPARAM lParam)
{
	if (m_pSingleCIRDlg) {
		m_pSingleCIRDlg->Close();
		delete (m_pSingleCIRDlg);
		m_pSingleCIRDlg	= NULL;
	}
	return 0;
}

LRESULT CCIRDisplayDlg::OnFloatingDlgStep(WPARAM wParam, LPARAM lParam)
{
	// single CIR dlg notified a step
	DWORD	dwC;
	if (lParam > 0) {
		// means step up
		dwC	= m_dwSingleCIRIndex + m_nFreqs;
		if (dwC < m_dwCIRCount) {
			m_dwSingleCIRIndex	= dwC;
			ShowSingleCIRAt(m_dwSingleCIRIndex);
		}
	}
	else {
		dwC	= m_nFreqs;
		if (m_dwSingleCIRIndex >= dwC) {
			m_dwSingleCIRIndex	-= m_nFreqs;
			ShowSingleCIRAt(m_dwSingleCIRIndex);
		}
	}
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main file parsing

double CCIRDisplayDlg::GetStartStopSlope()
{
	// returns the slope out of the first and the last measurements
	double	dStartAbsTime;
	double	dStartDelayTime;
	double	dEndAbsTime;
	double	dEndDelayTime;
	double	diff;

	// get the params of the beginning
	dStartAbsTime		= m_PeakEntries[0].dAbsTime;
	dStartDelayTime	= fmod(dStartAbsTime, m_dPeriod);
	// get the last time
	dEndAbsTime		= m_PeakEntries[m_dwEntrySize-1].dAbsTime;
	dEndDelayTime	= fmod(dEndAbsTime,  m_dPeriod);
	// remove wrap around
	diff	= dEndDelayTime - dStartDelayTime;
	while (fabs(diff) > m_dPeriod/2.0) {
		if (diff > 0.0)
			dEndDelayTime	-= m_dPeriod;
		else
			dEndDelayTime	+= m_dPeriod;
		diff	= dEndDelayTime - dStartDelayTime;
	}
	// final slope calculation
	return (dEndDelayTime - dStartDelayTime) / (dEndAbsTime - dStartAbsTime);
}


double CCIRDisplayDlg::GetDValue(CString str)
{
	int	iRef;

	iRef	= str.Find("=");
	return atof(str.Mid(iRef + 1));
}

void CCIRDisplayDlg::ParsePeaks(CString str, DWORD dwTimestamp, int nIndex, double dAbsTime, double dRSSI, double dSNR, int nPeaks)
{
	// parses a single peak string
	int			i;
	int			iRef;
	double		dLevel;
	double		dDelay;
	double		dFreqError;
	SPeakEntry	tempEntry;

	// keep limits
	m_dMinTime	= __min(m_dMinTime, dAbsTime);
	m_dMaxTime	= __max(m_dMaxTime, dAbsTime);
	m_dMinLevel	= __min(dRSSI, m_dMinLevel);
	m_dMaxLevel	= __max(dRSSI, m_dMaxLevel);

	// prepare standard values for the entry
	tempEntry.dwTimestamp	= dwTimestamp;
	tempEntry.nIndex			= nIndex;
	tempEntry.dAbsTime		= dAbsTime;
	tempEntry.dRSSI			= dRSSI;
	tempEntry.dSNR				= dSNR;
	// trim the peak string
	iRef	= str.Find("=");
	str	= str.Mid(iRef + 1);
	// loop through the peaks
	for (i=0; i<nPeaks; i++) {
		// extract the level
		iRef		= str.Find("|");
		dLevel	= atof(str.Left(iRef));
		str		= str.Mid(iRef + 1);
		// extract the delay
		iRef		= str.Find("|");
		dDelay	= atof(str.Left(iRef));
		str		= str.Mid(iRef + 1);
		// extract the frequency error
		iRef		= str.Find("|");
		if (iRef == -1)
			dFreqError	= atof(str);
		else {
			dFreqError	= atof(str.Left(iRef));
			str			= str.Mid(iRef + 1);
		}
		// finalyse entry values
		tempEntry.dLevel		= dLevel;
		tempEntry.dPeakTime	= dAbsTime + dDelay / 1.0E6;
		tempEntry.dFreqError	= dFreqError;
		// memory handling
		if (m_dwEntrySize == m_dwEntryAllocation) {
			m_dwEntryAllocation	+= 256;
			m_PeakEntries			= (SPeakEntry*) realloc(m_PeakEntries, m_dwEntryAllocation*sizeof(SPeakEntry));
			m_CIREntries			= (SCIREntry*) realloc(m_CIREntries, m_dwEntryAllocation*sizeof(SCIREntry));
		}
		// keep this entry
		m_PeakEntries[m_dwEntrySize]	= tempEntry;
		m_dwEntrySize++;
		// keep reference peak
		if (i == 0)
			m_CIREntries[m_dwCIRCount++].sPeakEntry	= tempEntry;
	}
	// all done
}

DWORD CCIRDisplayDlg::ParseLogFile()
{
	// collects all reasonable information from the file, returns the nb of valid entries

	BOOL						ok;
	CFile						logFile;
	CLineHelper				helper;
	CStringArray			itemArray;
	int						nItem;
	int						nIndex;
	int						icc;
	CString					tempItem;
	CString					workString;
	double					dAbsTime;
	double					dRSSI;
	double					dSNR;
	int						nPeaks;
	int						nMaxIndex;

	DWORD		dwTimestamp;
	int		nDaySeconds;
	double	dBandwidth;
	DWORD		dwN;

	// reset all information
	m_dwEntrySize	= 0;
	m_dMinTime		= 1.0E99;
	m_dMaxTime		= -1.0;
	m_dMinLevel		= 1.0E99;
	m_dMaxLevel		= -1.0E99;
	m_dwCIRCount	= 0;
	nMaxIndex		= -1;
	m_nFreqs			= 0;
	m_NavHelper.Reset();
	// try to open the file
	ok	= logFile.Open(m_StrINPFILE, CFile::modeRead, NULL);
	if (!ok)
		return 0;
	// read the header
	workString	= helper.ReadLine(&logFile);
	if (workString.IsEmpty())
		return 0;
	// get the initial values from the file
	while (1) {
		workString	= helper.ReadLine(&logFile);
		if (workString.IsEmpty()) {
			logFile.Close();
			return 0;
		}
		nItem			= helper.ExtractParameters(workString, &dwTimestamp, &nDaySeconds, &itemArray);
		if (nItem >3)
			if (itemArray[3] == "SETUP")
				break;
	}
	// get the bandwidth
	tempItem		= itemArray[5];
	dBandwidth	= GetDValue(tempItem);
	// convert to period
	m_dPeriod	= 65536 / dBandwidth;		// period in usec's
	m_dPeriod	/= 1.0E6;						// conversion to sec's
	// CIR collection loop
	while (1) {
		workString	= helper.ReadLine(&logFile);
		if (workString.IsEmpty())
			break;
		nItem	= helper.ExtractParameters(workString, &dwTimestamp, &nDaySeconds, &itemArray);
		// check for navigation
		if (nItem) {
			if (itemArray[3] == "NAVIGATION")
				m_NavHelper.AddNavigationEntry(&itemArray);
			else {
				if (itemArray[3] == "CHANNEL") {
					if (nItem == 9) {
						nIndex	= 0;
						icc		= 4;		// old file without index
					}
					else {
						nIndex	= int(GetDValue(itemArray[4]) + 0.5);
						icc		= 5;		// new file with index
					}
					dAbsTime	= GetDValue(itemArray[icc++]);
					dRSSI		= GetDValue(itemArray[icc++]);
					dSNR		= GetDValue(itemArray[icc++]);
					nPeaks	= helper.ExtractNumber(itemArray[icc++]);
					// finally parse the peaks
					ParsePeaks(itemArray[icc], dwTimestamp, nIndex, dAbsTime, dRSSI, dSNR, nPeaks);
					// keep the highest frequency index
					if (nIndex > nMaxIndex)
						nMaxIndex	= nIndex;
				}
			}
		}
	}
	// parsing complete
	logFile.Close();
	// keep freq count
	m_nFreqs	= nMaxIndex + 1;
	// finalyse the navigation information
	dwN	= m_NavHelper.Finalyse();
	if (dwN)
		m_NavHelper.BuildDistanceArray();
	// update annotations
	m_dMinLevel	= 10 * int(0.1*m_dMinLevel) - 10;
	m_dMaxLevel	= 10 * int(0.1*m_dMaxLevel);
	m_dMinTime	= int(m_dMinTime);
	m_dMaxTime	= int(m_dMaxTime) + 1;
	m_StrMINMEASTIME.Format("%.1f", m_dMinTime);
	m_StrMAXMEASTIME.Format("%.1f", m_dMaxTime);
	m_StrMINLEVEL.Format("%.1f", m_dMinLevel);
	m_StrMAXLEVEL.Format("%.1f", m_dMaxLevel);
	m_StrCYCLEOFFSET	= "0.0";
	m_StrCYCLEPERIOD.Format("%.1f", m_dPeriod*1.0E6);	// conversion to usec's
	m_StrSNRTHRESHOLD.Empty();
	m_StrMINPLACEMENT = "0.0";
	m_StrMAXPLACEMENT	= m_StrCYCLEPERIOD;
	UpdateData(FALSE);
	// drift compensation
	m_dStartStopSlope		= GetStartStopSlope();
	// all done
	return m_dwEntrySize;
}

DWORD CCIRDisplayDlg::ParseCIRFile()
{
	// collects all primary information from the CIR file, returns the nb of CIR entries found
	BOOL			ok;
	DWORD			dwC;
	int			nBufferSize;
	int			nPrevious;
	int			nBufferPosition;
	int			nCIRStart;
	int			nCIREnd;
	int			nIndex;

	int			nOffset;
	char			tBuffer[8192];
	char			ccc;
	CFile			cirFile;
	CCIRHelper	cirHelper;

	// open the file
	nIndex	= m_StrINPFILE.ReverseFind('.');
	if (nIndex == -1)
		return 0;
	ok	= cirFile.Open(m_StrINPFILE.Left(nIndex) + ".cir", CFile::modeRead, NULL);
	if (!ok)
		return 0;
	// prepare the loop
	nOffset				= 0;
	dwC					= 0;
	nCIRStart			= -1;
	nCIREnd				= -1;
	nBufferSize			= 0;
	nBufferPosition	= 0;
	nPrevious			= 0;
	nIndex				= 0;
	// collection loop
	while (1) {
		// check whether enough information is available
		if (nBufferSize == 0) {
			nOffset	+= nPrevious;
			nBufferSize = cirFile.Read(tBuffer, 8192);
			if (nBufferSize == 0)
				break;	// end of file reached
			nPrevious			= nBufferSize;
			nBufferPosition	= 0;
		}
		// analysis loop
		ccc	= tBuffer[nBufferPosition];
		switch (ccc) {
		case 'x':
			nIndex	= 1;
			break;
		case '=':
			if (nIndex == 1) {
				// set the starting point for the CIR
				nCIRStart	= nOffset + nBufferPosition + 2;	// skip equal sign and frequency index character
			}
			break;
		case 13:
			nIndex	= -1;
			break;
		case 10:
			if ((nIndex == -1) && (nCIRStart > 0)) {
				// set the end point and finalyze the entry
				nCIREnd	= nOffset + nBufferPosition - 2;		// skip the finaly <CR><LF>
				m_CIREntries[dwC].ulFilePosition	= nCIRStart;
				m_CIREntries[dwC].nStringLen		= nCIREnd - nCIRStart + 1;
				nCIRStart								= -1;
				dwC++;
			}
			break;
		}
		// next element
		nBufferPosition++;
		nBufferSize--;
	}
	// finally get the CIR information
	cirFile.SeekToBegin();
	nIndex	= cirHelper.GetCIRRange(&cirFile, &m_dCIRStartTime, &m_dCIREndTime, &m_dCIRResolution);
	// reset the display indicator
	m_dwSingleCIRIndex	= 0;
	// clean up
	cirFile.Close();
	// all done
	if (nIndex)
		return 0;
	else
		return dwC;
}

void CCIRDisplayDlg::OnStnDblclickPeakdisplayDlgGrapPeakdisplay()
{
	// TODO: Fügen Sie hier Ihren Kontrollbehandlungscode für die Benachrichtigung ein.
	CPoint					pt;
	CRect						windowRect;
	double					dYValue;
	DWORD						dw;
	DWORD						dwRef;

	// get the current position of the cursor
	GetCursorPos(&pt);
	// get the position in the grap
	m_GrapPEAKDISPLAY.GetWindowRect(&windowRect);
	// calculate the position
	pt.x	-= windowRect.left;
	pt.y	-= windowRect.top;
	// check whether the position is within the viewport
	if ((pt.x < m_GrapPEAKDISPLAY.viewport.left) || (pt.x > m_GrapPEAKDISPLAY.viewport.right))
		return;
	if ((pt.y < m_GrapPEAKDISPLAY.viewport.top) || (pt.y > m_GrapPEAKDISPLAY.viewport.bottom))
		return;
	// calculate the position using the display xyl
	dYValue	= m_GrapPEAKDISPLAY.xyl[3] - (m_GrapPEAKDISPLAY.xyl[3] - m_GrapPEAKDISPLAY.xyl[1]) * (pt.y - m_GrapPEAKDISPLAY.viewport.top) / double(m_GrapPEAKDISPLAY.viewport.Height());
	// locate the nearest CIR information
	if (m_CheckDISTVIEW) {
		dwRef	= m_NavHelper.m_dwDistanceArraySize;
		for (dw=0; dw<m_NavHelper.m_dwDistanceArraySize; dw++) {
			if (m_NavHelper.m_DistanceArray[dw] >= dYValue) {
				dwRef	= dw;
				dw		= m_NavHelper.m_dwDistanceArraySize;		// to quit the loop
			}
		}
		if (dwRef == m_NavHelper.m_dwDistanceArraySize)
			return;
		dYValue	= m_PeakEntries[dwRef].dAbsTime;
	}
	// time based search in any case
	dwRef	= m_dwCIRCount;
	for (dw=0; dw<m_dwCIRCount; dw++) {
		if (m_CIREntries[dw].sPeakEntry.dAbsTime >= dYValue) {
			dwRef	= dw;
			dw		= m_dwCIRCount;	// to quit the loop;
		}
	}
	if (dwRef == m_dwCIRCount)
		return;
	// simple display
	ShowSingleCIRAt(dwRef);
	// keep placement
	m_dwSingleCIRIndex	= dwRef;
}

void CCIRDisplayDlg::ShowSingleCIRAt(DWORD dwIndex)
{
	// displays the CIR's starting at the given position
	BOOL						ok;
	CString					strFile;
	CString					str;
	CFile						cirFile;
	DWORD						dwWantedMask;
	DWORD						dwCurrentMask;
	DWORD						dwM;
	DWORD						dwC;
	int						nRet;
	int						nCIRSize;
	double					dPeriod;
	double					dOffset;
	double					dCenterTime;
	double*					cirArray;
	double*					xVec;
	CCIRHelper				cirHelper;
	int						offset[4]	= {8, 8, 2, 2};
		
	// get basic information
	UpdateData(TRUE);
	dPeriod	= atof(m_StrCYCLEPERIOD);
	dOffset	= atof(m_StrCYCLEOFFSET);
	// create the file name
	nRet	= m_StrINPFILE.ReverseFind('.');
	if (nRet == -1)
		return;
	ok	= cirFile.Open(m_StrINPFILE.Left(nRet) + ".cir", CFile::modeRead, NULL);
	if (!ok)
		return;
	// prepare the buffer to hold the information
	nCIRSize		= int((m_dCIREndTime - m_dCIRStartTime) / m_dCIRResolution + 1.5);
	cirArray		= (double*) calloc(nCIRSize, sizeof(double));
	xVec			= (double*) calloc(nCIRSize, sizeof(double));
	// read operation for the first CIR
	nRet	= BuildSingleCIR(&cirFile, nCIRSize, dwIndex, dPeriod, dOffset, xVec, cirArray);
	if (nRet) {
		free (cirArray);
		free (xVec);
		cirFile.Close();
		return;
	}
	// create the params to collect the CIR's at other frequencies
	dwWantedMask	= (1 << m_nFreqs);
	dwCurrentMask	= (1 << m_CIREntries[dwIndex].sPeakEntry.nIndex);
	dCenterTime		= m_CIREntries[dwIndex].sPeakEntry.dAbsTime;
	// check whether the grap has to be created
	if (m_pSingleCIRDlg == NULL) {
		m_pSingleCIRDlg = new CFloatingDlg();
		m_pSingleCIRDlg->Create(this, 0, "Single CIR display");
		m_pSingleCIRDlg->m_GrapField.InitGraph(nCIRSize, 0, GRAPHICS_MODE_XY, DISPLAY_MODE_LABELS, offset);
		m_pSingleCIRDlg->m_StrMINX		= m_StrMINPLACEMENT;
		m_pSingleCIRDlg->m_StrMAXX		= m_StrMAXPLACEMENT;
		m_pSingleCIRDlg->m_StrMINY	= "-40.0";
		m_pSingleCIRDlg->m_StrMAXY	= "0.0";
		m_pSingleCIRDlg->UpdateData(FALSE);
		m_pSingleCIRDlg->m_GrapField.SetMaxNbOfAnnotations(12, 8);
		m_pSingleCIRDlg->m_GrapField.labelX	= "Delay [usec]";
		m_pSingleCIRDlg->m_GrapField.labelY	= "Rel level [dB]";
	}
	// create the header
	str.Format("CIR at %.6f sec", m_CIREntries[dwIndex].sPeakEntry.dAbsTime);
	m_pSingleCIRDlg->SetWindowText(str);
	// prepare the window
	m_pSingleCIRDlg->SetSize(m_nFreqs, nCIRSize, m_ColorTable);
	m_pSingleCIRDlg->SetData(m_CIREntries[dwIndex].sPeakEntry.nIndex, xVec, cirArray);
	// try to collect the other freq's
	dwC	= dwIndex;
	while (dwWantedMask != dwCurrentMask) {
		if (dwC)
			dwC--;
		else
			break;
		if ((dCenterTime - m_CIREntries[dwC].sPeakEntry.dAbsTime) > 2.0)
			break;
		dwM	= (1 << m_CIREntries[dwC].sPeakEntry.nIndex);
		if ((dwM & dwCurrentMask) == 0) {
			nRet	= BuildSingleCIR(&cirFile, nCIRSize, dwC, dPeriod, dOffset, xVec, cirArray);
			if (nRet == 0) {
				dwCurrentMask	|= dwM;
				m_pSingleCIRDlg->SetData(m_CIREntries[dwC].sPeakEntry.nIndex, xVec, cirArray);
			}
		}
	}
	dwC	= dwIndex;
	while (dwWantedMask != dwCurrentMask) {
		dwC++;
		if (dwC >= m_dwCIRCount)
			break;
		if ((m_CIREntries[dwC].sPeakEntry.dAbsTime - dCenterTime) > 2.0)
			break;
		dwM	= (1 << m_CIREntries[dwC].sPeakEntry.nIndex);
		if ((dwM & dwCurrentMask) == 0) {
			nRet	= BuildSingleCIR(&cirFile, nCIRSize, dwC, dPeriod, dOffset, xVec, cirArray);
			if (nRet == 0) {
				dwCurrentMask	|= dwM;
				m_pSingleCIRDlg->SetData(m_CIREntries[dwC].sPeakEntry.nIndex, xVec, cirArray);
			}
		}
	}
	// show information in any case
	m_pSingleCIRDlg->OnBnClickedFloatingDlgBtnUpdate();
	// clean up
	free (cirArray);
	free (xVec);
	cirFile.Close();
	// all done
}



int CCIRDisplayDlg::BuildSingleCIR(CFile* pFile, int nCIRSize, DWORD dwIndex, double dPeriod, double dOffset, double* xData, double* yData)
{
	// builds a single time corrected CIR at the given index, returns 0 if the information is valid
	int			i;
	int			nStringSize;
	int			nSize;
	char*			cString;
	double		dTime;
	CCIRHelper	cirHelper;

	nStringSize	= m_CIREntries[dwIndex].nStringLen;
	cString		= (char*) calloc(nStringSize + 1, 1);
	// read operation
	pFile->Seek(m_CIREntries[dwIndex].ulFilePosition, CFile::begin);
	nSize	= pFile->Read(cString, nStringSize);
	if (nSize != nStringSize) {
		free (cString);
		return -1;
	}
	nSize	= cirHelper.ParseCIRString(cString, yData);
	// string not used anymore
	free (cString);
	if (nSize != nCIRSize)
		return -1;
	// assemble the x-vector
	dTime	= m_CIREntries[dwIndex].sPeakEntry.dAbsTime * 1.0E6;
	dTime	= fmod(dTime, dPeriod);
	dTime	+= dOffset;
	dTime	+= m_dCIRStartTime;
	for (i=0; i<nCIRSize; i++) {
		xData[i]	= dTime;
		dTime		+= m_dCIRResolution;
	}
	// all done
	return 0;
}