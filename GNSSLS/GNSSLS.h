

#ifndef GNSSLS_H
#define GNSSLS_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include <regex>

#include "rinex.h"
#include "NRinexUtils.h"

using namespace Eigen;
using namespace std;
using namespace NGSrinex;


class GNSSLS
{
public:	
	const double PI = (3.141592653589793238463);
	const double RAD2ARC = (3600 * 180 / PI);
	const double RAD2DEG = (180 / PI);
	GNSSLS();
	//Readers
	void opensatpos(string satpos);
	void readl(NGSrinex::ObsEpoch  &currentRinexObs);
	void readEphem(NGSrinex::PRNBlock  &currentRinexEph);
	void readsatpos();
	//Setters
	void setIter(int num_iter_temp);
	void setTol(double tol_temp);
	void setSatPos(MatrixXd satPostemp);
	void setObsSats(VectorXd raw_ObsSats);
	void setA(const MatrixXd &A_temp);
	void setAPrior(double prior);
	void setL(const VectorXd & l_temp);
	void setRawL(const VectorXd & l_temp);
	void setX0(const VectorXd & x0_temp);
	void setxh(const VectorXd & xh_temp);
	void setCl(const MatrixXd & Cl_diag);
	void setP(const VectorXd diag_p);
	void setfx(const MatrixXd &fx_temp);
	void setXYZTrue(double x_t, double y_t, double z_t);
	void setElevationWeight();
	void setSingleDiffBasePos(double pBase_);
	void correctObs();
	void debug();
	void quiet();
	//Getters
	MatrixXd getsatPos();
	VectorXd getxh();
	VectorXd getx0();
	VectorXd getlh();
	VectorXd getd();
	VectorXd getw();
	double getAPost();
	MatrixXd getP();
	MatrixXd getA();
	MatrixXd getN();
	MatrixXd getNi();
	MatrixXd getCl();
	MatrixXd getClh();
	MatrixXd getCvh();
	MatrixXd getCxh();
	VectorXd getv();
	VectorXd getL();
	VectorXd getfx();
	VectorXd getobsSats();
	vector<vector<double>> getSatPos();
	vector<double> getEpochSum();
	vector<vector<double>> getEphem();
	VectorXd computeLatLong(double x, double y, double z);
	void clear();

	double dist(double xi, double yi, double zi, double xj, double yj, double zj);

	void computeA();
	void computeLS();
	void computeSP();
	void computeSatPos(double interval, int duration, bool updateepoch);
	void computeSingleDifference(string roverfilename, string basefilename);

	//initializers
	MatrixXd A, Cl, P, N, U, Ni, R;
	VectorXd l,raw_l, lh, x0, xh, xh_loc, fx, d, w, v;

private:

	int num_iter, n, u, r, obs;
	double tol, aprior, apost, DOF, obsTime, num_sats, epoch, pBase, pRov, sat_elev, x_true, y_true, z_true, p_true, l_true, h_true;
	double xDOP, yDOP, zDOP, tDOP, hDOP, vDOP, pDOP, gDOP;
	bool isquiet, isdebug, obsCorrected, isDiff, satElevWeight ;
	ifstream satposfile;
	vector<double> epochSum, ephem_row;
	vector<vector<double>> ephem, calcSatPos;

	MatrixXd Clh, Cxh, Cvh, y, satPos, Qloc;
	VectorXd obsSats, raw_obsSats, curviCoords;

	//Checkers
	void checkObstoSat();


	//computers
	void computeaPost();
	void computeP();
	void computeU();
	void computeN();
	void computeNi();
	void computeResiduals();
	void computeDelta();
	void computeClh();
	void computeCvh();
	void computeCxh();
	void computeDOP();
	MatrixXd computeInverse(MatrixXd &M);
	void computefx();
	void computeMisclosure();
	void computeEpochSum();
	void computeElevationWeight();
	void computeLocalCoords(double X, double Y, double Z);
	void computeSatPosParams(double tot, double tk, double a, double Deltan, double Mo, double e, double LilOmega, double Cuc, double Cus, double Crc, double Cis, double Cic, double Crs, double Idot, double I0, double BigOmegadot, double BigOmega, double toe);
	int findNearestEpoch(double transtime);
	MatrixXd computeRotMatrix(double p, double l, double h);

	double getMaxDelta();

	void updateObservations();
	void updateUnknowns();
	void updateCoords();


};

#endif