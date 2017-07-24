#pragma once
#include "afxwin.h"
#include "CRTDisplayDlg.h"
#include "PeakDisplayGraph.h"
#include "Scaling3D.h"
#include "FloatingDlg.h"



// CCIRDisplayDlg-Dialogfeld

class CCIRDisplayDlg : public CDialog
{
	DECLARE_DYNAMIC(CCIRDisplayDlg)

public:
	CCIRDisplayDlg(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CCIRDisplayDlg();

// Dialogfelddaten
	enum { IDD = IDD_CIRDISPLAY_DLG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	DECLARE_MESSAGE_MAP()

	virtual BOOL OnInitDialog();
	afx_msg void OnBnClickedCirdisplayDlgBtnInpfile();
	afx_msg void OnBnClickedCirdisplayDlgBtnUpdate();
	afx_msg void OnBnClickedCirdisplayDlgBtnCrtview();
	afx_msg void OnBnClickedCirdisplayDlgBtnDopplerview();
	afx_msg void OnStnDblclickPeakdisplayDlgGrapPeakdisplay();
	afx_msg void OnBnClickedCirdisplayDlgBtnSimulation();

public:
	CString m_StrINPFILE;
	CString m_StrMINMEASTIME;
	CString m_StrMAXMEASTIME;
	CString m_StrCYCLEPERIOD;
	CString m_StrCYCLEOFFSET;
	CString m_StrSNRTHRESHOLD;
	CString m_StrMINPLACEMENT;
	CString m_StrMAXPLACEMENT;
	CString m_StrMINLEVEL;
	CString m_StrMAXLEVEL;
	CString m_StrFREQINDEX;
	CString m_StrCAPTURERANGE;
	

	CPeakDisplayGraph m_GrapPEAKDISPLAY;
	CScaling3D			m_GrapSCALING;
	BOOL					m_CheckMAINPEAK;
	BOOL					m_CheckSTARTSTOP;
	BOOL					m_CheckMULTIFREQ;
	BOOL					m_CheckDISTVIEW;

private:
	DWORD		ParseLogFile();
	DWORD		ParseCIRFile();
	void		ParsePeaks(CString str, DWORD dwTimestamp, int nIndex, double dAbsTime, double dRSSI, double dSNR, int nPeaks);
	int		BuildSingleCIR(CFile* pFile,  int nCIRSize, DWORD dwIndex, double dPeriod, double dOffset, double* xData, double* yData);
	void		ShowSingleCIRAt(DWORD dwIndex);
	double	GetDValue(CString str);
	double	GetStartStopSlope();

	LRESULT	OnFloatingDlgClose(WPARAM wParam, LPARAM lParam);
	LRESULT	OnFloatingDlgStep(WPARAM wParam, LPARAM lParam);

private:
	CNavigationHelper	m_NavHelper;
	int					m_nFreqs;
	DWORD					m_dwEntryAllocation;
	DWORD					m_dwEntrySize;
	DWORD					m_dwCIRCount;

	COLORREF	m_ColorTable[8];

	SPeakEntry*	m_PeakEntries;
	SCIREntry*	m_CIREntries;

	CFloatingDlg*	m_pSingleCIRDlg;
	DWORD				m_dwSingleCIRIndex;


	double	m_dMinTime;
	double	m_dMaxTime;
	double	m_dMinLevel;
	double	m_dMaxLevel;
	double	m_dPeriod;
	double	m_dStartStopSlope;
	double	m_dOffset;

	double	m_dCIRStartTime;
	double	m_dCIREndTime;
	double	m_dCIRResolution;
public:
	
};
