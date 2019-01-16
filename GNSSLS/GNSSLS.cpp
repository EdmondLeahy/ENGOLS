#include "GNSSLS.h"


GNSSLS::GNSSLS()
{
	num_iter = 0;
	n = 0;
	u = 0;
	clear();
	isdebug = false;
	isquiet = false;
	aprior = 1;
}

void GNSSLS::clear()
{
	A = MatrixXd::Zero(1, 1);
	Cl = MatrixXd::Zero(1, 1);
	Cxh = MatrixXd::Zero(1, 1);
	xh = VectorXd::Zero(1);
	w = VectorXd::Zero(1);
	fx = VectorXd::Zero(1);
	l = VectorXd::Zero(1);
	isdebug = false;
	isquiet = false;
	obsCorrected = false;
	satElevWeight = false;
}

void GNSSLS::opensatpos(string satpos)
{
	//Open the sats file
	satposfile.open(satpos);
	if (!satposfile.is_open()) { 
		fprintf(stdout, "\nUnable to open %s. Program is exiting....", satpos); 
		exit; 
	}
	
}

void GNSSLS::readl(NGSrinex::ObsEpoch  &currentRinexObs)
{
	
	//initialize number of satellites
	raw_l = VectorXd::Zero(currentRinexObs.getNumSat());
	raw_obsSats = VectorXd::Zero(currentRinexObs.getNumSat());
	//Get the observation time
	obsTime = currentRinexObs.getEpochTime().GetGPSTime().secsOfWeek;

	//For all satellites
	for (unsigned short i = 0; i < currentRinexObs.getNumSat(); ++i) {
		// get the satellite observations for the current satellite
		NGSrinex::SatObsAtEpoch satObs = currentRinexObs.getSatListElement(i);
		//We are dealing with GPS only
		if (satObs.satCode != 'G') { continue; }
		// for all observations from the current satellite (MAXOBSTYPES is defined in rinex.h)...

		//Place the SV# in the obsSat matrix
		raw_obsSats(i) = double(satObs.satNum);

		for (unsigned short j = 0; j < MAXOBSTYPES; ++j)
		{
			// Pseudorange measurements on L1 only (denoted C1; defined in rinex.h)
			if (satObs.obsList[j].obsType != C1)
				continue;

			// make sure the observation is present
			if (!satObs.obsList[j].obsPresent)
				continue;

			// Place pseudorange in observation matrix
			raw_l(i) = satObs.obsList[j].observation;
		}
	}
	//number of gps obs
	obs = l.rows();
	
}

void GNSSLS::readEphem(NGSrinex::PRNBlock &currentRinexEph)
{
	ephem_row.resize(18);

	// check for unhealthy satellite
	if (currentRinexEph.getSvHealth() != 0) {
		return;
	}

	else {
	
		//Put the data into the vector
		ephem_row[0] = currentRinexEph.getSatellitePRN();  //PRN
		ephem_row[1] = currentRinexEph.getIode();  //PRN
		ephem_row[2] = currentRinexEph.getBigOmega();      // longitude of ascending node at reference time [rad]
		ephem_row[3] = currentRinexEph.getBigOmegaDot();   // rate of longitude of ascending node [rad/s]
		ephem_row[4] = currentRinexEph.getCic();           // harmonic correction term for inclination (1 of 2) [rad]
		ephem_row[5] = currentRinexEph.getCis();           // harmonic correction term for inclination (2 of 2) [rad]
		ephem_row[6] = currentRinexEph.getCrc();           // harmonic correction term for radius (1 of 2) [m]
		ephem_row[7] = currentRinexEph.getCrs();           // harmonic correction term for radius (2 of 2) [m]
		ephem_row[8] = currentRinexEph.getCuc();           // harmonic correction term for argument of perigee (1 of 2) [rad]
		ephem_row[9] = currentRinexEph.getCus();           // harmonic correction term for argument of perigee (2 of 2) [rad]
		ephem_row[10] = currentRinexEph.getDeltan();        // correction to mean motion [rad/s]
		ephem_row[11] = currentRinexEph.getEccen();         // eccentricity [unitless]
		ephem_row[12] = currentRinexEph.getIdot();          // inclination rate [rad/s]
		ephem_row[13] = currentRinexEph.getIo();            // inclincation at reference time [rad]
		ephem_row[14] = currentRinexEph.getLilOmega();      // argument of perigee [rad]
		ephem_row[15] = currentRinexEph.getMo();            // mean anomaly at reference time [rad]
		ephem_row[16] = currentRinexEph.getSqrtA();         // square root of semi-major axis of orbit [sqrt-m]
		ephem_row[17] = currentRinexEph.getToe();           // reference time, i.e., 'time of ephemeris' [s]
		ephem.push_back(ephem_row);	
	}
	
}

void GNSSLS::readsatpos()
{
	string line;
	//Clear the buffer, in case of previous attemps
	satposfile.clear();
	while (getline(satposfile, line)) {
		stringstream linestream(line);
		linestream >> epoch >> num_sats;
		//check if it is the correct epoch
		if (abs(epoch - obsTime) < 0.01) {
			//Correct epoch
			if (isdebug) { cout << "Found correct epoch!\n"; }
			satPos = MatrixXd::Zero(num_sats, 5);
			for (int i = 0; i < num_sats; i++) {
				getline(satposfile, line);
				stringstream linestream(line);
				linestream >> satPos(i, 0) >> satPos(i, 1) >> satPos(i, 2) >> satPos(i, 3) >> satPos(i, 4);
			}

		}
		break;
	}



}

void GNSSLS::computeSingleDifference(string Roverfilename, string Basefilename)
{
	
	//ifstream Observations;
	//vector<double> satPosTemp;
	//double tempx, tempy, tempz;
	//Observations.open(Roverfilename);
	//double dataEpoch = 0;
	//double numSatsEpoch = 0;
	//string line;
	//while (getline(Observations, line)) {
	//	Observations >> dataEpoch >> numSatsEpoch;
	//	for (int i = 0; i < numSatsEpoch; i++) {
	//		Observations >> tempx >> tempy >> tempz >> 
	//	}
	//
	//
	//
	//}

}


void GNSSLS::setIter(int num_iter_temp)
{
	num_iter = num_iter_temp;
}

void GNSSLS::setTol(double tol_temp)
{
	tol = tol_temp;
}

void GNSSLS::setSatPos(MatrixXd satPostemp)
{
	satPos = satPostemp;
}

void GNSSLS::setObsSats(VectorXd rawObsSats)
{
	raw_obsSats = rawObsSats;
}


void GNSSLS::setAPrior(double prior)
{
	aprior = prior;
}

void GNSSLS::setL(const VectorXd & l_temp)
{
	l = l_temp;
}

void GNSSLS::setRawL(const VectorXd & l_temp)
{
	raw_l = l_temp;
}

void GNSSLS::setX0(const VectorXd & x0_temp)
{
	x0 = x0_temp;
}

void GNSSLS::setxh(const VectorXd & xh_temp)
{

}

void GNSSLS::setCl(const MatrixXd & Cl_diag)
{
	Cl = Cl_diag.asDiagonal();
}

void GNSSLS::setP(const VectorXd diag_p)
{
	P = diag_p.diagonal();
}

void GNSSLS::setA(const MatrixXd & A_temp)
{
	A = A_temp;
}

void GNSSLS::setfx(const MatrixXd & fx_temp)
{
	fx = fx_temp;
}

void GNSSLS::setXYZTrue(double x_t, double y_t, double z_t)
{
	x_true = x_t;
	y_true = y_t;
	z_true = z_t;



}

void GNSSLS::setElevationWeight()
{
	satElevWeight = true;
}

void GNSSLS::setSingleDiffBasePos(double pBase_)
{
	isDiff = true;
	pBase = pBase_;
}

void GNSSLS::correctObs()
{
	for (int i = 0; i < num_sats; i++) {
		l(i) = l(i) - satPos(i, 4);
	}

	obsCorrected = true;
}


void GNSSLS::debug()
{
	isdebug = true;
	isquiet = false;
}

void GNSSLS::quiet()
{
	isquiet = true;
}

MatrixXd GNSSLS::getsatPos()
{
	return satPos;
}

VectorXd GNSSLS::getxh()
{
	return xh;
}

VectorXd GNSSLS::getx0()
{
	return x0;
}

VectorXd GNSSLS::getlh()
{
	return lh;
}

VectorXd GNSSLS::getd()
{
	return d;
}

VectorXd GNSSLS::getw()
{
	return w;
}


double GNSSLS::getAPost()
{
	return apost;
}

MatrixXd GNSSLS::getP()
{
	return P;
}

MatrixXd GNSSLS::getA()
{
	return A;
}

MatrixXd GNSSLS::getN()
{
	return N;
}

MatrixXd GNSSLS::getNi()
{
	return Ni;
}

MatrixXd GNSSLS::getCl()
{
	return Cl;
}

MatrixXd GNSSLS::getClh()
{
	return Clh;
}

MatrixXd GNSSLS::getCvh()
{
	return Cvh;
}

MatrixXd GNSSLS::getCxh()
{
	return Cxh;
}

VectorXd GNSSLS::getv()
{
	return v;
}

VectorXd GNSSLS::getL()
{
	return l;
}

VectorXd GNSSLS::getfx()
{
	return fx;
}

VectorXd GNSSLS::getobsSats()
{
	return obsSats;
}

vector<vector<double>> GNSSLS::getSatPos()
{
	return calcSatPos;
}

vector<double> GNSSLS::getEpochSum()
{
	return epochSum;
}

vector<vector<double>> GNSSLS::getEphem()
{
	return ephem;
}

void GNSSLS::computeLS()
{
	//Compute non-iterative LS
	//Set initial parameters
	if (!isquiet) { cout << "\nRunning linear LeastSquares estimation......"; }
	computeA();
	n = A.rows();
	u = A.cols();
	//r = computeRank(A);
	computeP();
	computeDelta();
	updateUnknowns();
	computeResiduals();
	computeaPost();
	computeClh();
	computeCxh();
	updateObservations();

	if (!isquiet) { cout << "FINISHED\n"; }

}

void GNSSLS::computeSP()
{

	if (!isquiet) { fprintf(stdout, "\nRunning Non-Linear LeastSquares estimation for max %i iterations......\n", num_iter); }
	num_sats = raw_l.rows();
	//Need to set some initial things:	
	checkObstoSat();//Check the number of observations, make sure they match number of satelites
	if (!obsCorrected) { correctObs(); }//correct the psr observations
	n = l.rows();
	u = x0.rows();
	r = u;
	DOF = n - u;
	computeP();

	if (!isquiet) { cout << "This network has " << DOF << " degrees of freedom\n"; }

	bool Delta = false;
	int i = 0;
	while ((!Delta) && (i < num_iter)) {
		//update iteration counter
		i++;
		if (isquiet != true && isdebug) {

			cout << "-------------------------------------------\n" << endl;
			cout << "START ITERATION: " << i << endl;
			cout << "x0:\n" << x0 << endl;
			cout << "P:\n" << P << endl;
			cout << "SatPos:\n" << satPos << endl;
			cout << "L:\n" << l << endl;
			cout << "Epoch:\n" << epoch << endl;
			cout << "ObsTime:\n" << obsTime << endl;
			cout << "Sats:\n" << num_sats << endl;
			cout << "n:\n" << n << endl;
			cout << "u:\n" << u << endl;
		}
		
		//compute design matrix for LS Calculations
		computeA();

		//Compute Elevation Wise Weight
		if (satElevWeight) { computeElevationWeight(); }

		//compute the Delta
		computeDelta();

		//update the obs and unkns
		updateUnknowns();
		//Check the tolerance
		if (getMaxDelta() <= tol) {
			Delta = true;
		}

		if (isquiet != true && isdebug) {
			//TESTING
			cout << "L:\n" << l << endl;
			cout << "v:\n" << v << endl;
			cout << "Fx:\n" << fx << endl;
			cout << "d:\n" << d << endl;
			cout << "w:\n" << w << endl;
			cout << "A:\n" << A << endl;
			cout << "xh:\n" << xh << endl;
			//cout << "P:\n" << P << endl;
			//cout << "Cl:\n" << Cl << endl;
			cout << "N:\n" << N << endl;
			cout << "U:\n" << U << endl;
			cout << "n: " << n << endl;
			cout << "u: " << u << endl;
		}
	}
	if (!isquiet) { cout << "FINISHED\n"; }
	if (!isquiet) { cout << "Computation took: " << i << " iterations \n" << endl; }


	computeResiduals();//vh
	updateObservations();//lh
	computeaPost();//Apost
	computeCxh();//Cx


	computeClh();//Clh
	computeCvh();//Cvh

	computeEpochSum();



}
double GNSSLS::getMaxDelta()
{
	double max = 0.0;

	for (int i = 0; i < d.size(); i++) {
		if (abs(d(i)) > max) {
			max = abs(d(i));
		}

	}
	if (isquiet != true) { cout << "Max Delta for it: " << max << endl; }
	return max;
}

MatrixXd GNSSLS::computeInverse(MatrixXd & M)
{
	MatrixXd temp_I;
	double rows = M.rows();
	double cols = M.cols();
	if (rows == cols) {
		if (abs(M.determinant()) > 1 * pow(10, -20)) {
			//cout << "Invertable matrix. General inverse method applied" << temp_I << endl;
			temp_I = M.inverse();
		}
		else {
			cout << "Matrix is not invertable. Check observations\n" << endl;
		}
	}
	return temp_I;
}
double GNSSLS::dist(double xi, double yi, double zi, double xj, double yj, double zj)
{
	double dist = sqrt(pow(xi - xj, 2) + pow(yi - yj, 2) + pow(zi - zj, 2));
	return dist;
}
void GNSSLS::checkObstoSat()
{
	//We cannot assume that the data is in the same order
	//We use only ranges that have matching sat positions
	l = VectorXd::Zero(num_sats);
	obsSats = VectorXd::Zero(num_sats);
	for (int i = 0; i < num_sats; i++) {
		if (raw_obsSats[i] == satPos(i, 0)) {
			//the data is in same location
			l(i) = raw_l(i);
		}
		else {
			//data is not in correct location
			for (int j = 0; j < raw_obsSats.rows(); j++) {
				if (raw_obsSats[j] == satPos(i, 0)) {
					l(i) = raw_l(j);
					obsSats[i] = satPos(i, 0);
				}
				else {
					continue;
				}
			}
			if (l(i) == 0) {
				//We did not find matching range
				fprintf(stdout, "\n!No matching range for SV#: %f{1} !\n", satPos(i, 0));
				exit(1);
			}

		}

	}

}
void GNSSLS::computeaPost()
{
	VectorXd temp;
	temp = (v.transpose() * P * v) / (n - r);
	apost = pow(temp[0], 0.5);
}
void GNSSLS::computeP()
{
	P = pow(aprior, 2) * Cl.inverse();
}

void GNSSLS::computeU()
{
	computeMisclosure();
	U = A.transpose() * P * w;
}

void GNSSLS::computeN()
{
	N = A.transpose() * P * A;
	computeNi();

}

void GNSSLS::computeNi()
{
	Ni = computeInverse(N);
}

void GNSSLS::updateObservations()
{
	lh = l + v;
}

void GNSSLS::updateUnknowns()
{
	xh = x0 + d;
	x0 = xh;
	updateCoords();
}

void GNSSLS::computeResiduals()
{

	v = A * d - w;
}

void GNSSLS::computeMisclosure()
{
	computefx();
	w = fx - l;
}

void GNSSLS::computeEpochSum()
{
	//The summary of the epoch
	//Formatted the following way:
	//| #SV used in sol| X | Y | Z | dt | Lat |  Long | Height | N | E | H |GDOP | PDOP | HDOP | VDOP | apost |

	computeLocalCoords(x0(0), x0(1), x0(2));
	computeDOP();
	epochSum.resize(17);
	epochSum[0] = num_sats;
	epochSum[1] = xh[0];
	epochSum[2] = xh[1];
	epochSum[3] = xh[2];
	epochSum[4] = xh[3];
	epochSum[5] = curviCoords[0] * RAD2DEG;
	epochSum[6] = curviCoords[1] * RAD2DEG;
	epochSum[7] = curviCoords[2];
	epochSum[8] = xh_loc[0];
	epochSum[9] = xh_loc[1];
	epochSum[10] = xh_loc[2];
	epochSum[11] = xh_loc[3];
	epochSum[12] = gDOP;
	epochSum[13] = pDOP;
	epochSum[14] = hDOP;
	epochSum[15] = vDOP;
	epochSum[16] = apost;


}

void GNSSLS::computeElevationWeight()
{

	VectorXd TrueLatLong, satVectorXYZ, satVectorENU;
	double az;
	//Compute Lat Long Height
	VectorXd TruelatLong = computeLatLong(x_true, y_true, z_true);
	p_true = TruelatLong(0);
	l_true = TruelatLong(0);
	h_true = TruelatLong(0);

	//Compute Rotation matrix
	R = computeRotMatrix(p_true, l_true, h_true);
	//Calculate weight for each satellite
	satVectorXYZ = VectorXd::Zero(4);
	satVectorENU = VectorXd::Zero(4);
	P = MatrixXd::Zero(satPos.rows(), satPos.rows());
	for (int i = 0; i < satPos.rows(); i++) {
		//find the vector to the satellite
		satVectorXYZ(0) = satPos(i, 1) - x0(0);
		satVectorXYZ(1) = satPos(i, 2) - x0(1);
		satVectorXYZ(2) = satPos(i, 3) - x0(2);
		satVectorXYZ(3) = 0;

		//cout << "\nR:\n" << R << "\nVecXYZ:\n" << satVectorXYZ << endl;

		satVectorENU = R*satVectorXYZ;

		az = atan(satVectorXYZ(0) / satVectorXYZ(1));

		P(i, i) = pow(sin(az), 2);

		//cout << "\nP:\n" << P << endl;
}




}

void GNSSLS::computeLocalCoords(double X, double Y, double Z)
{
	//---------------Puts the coordinates in the local frame----------------
	MatrixXd R;

	//Need Lat and Long coords
	curviCoords = computeLatLong(X, Y, Z);

	//Rotation matrix
	R = computeRotMatrix(curviCoords[0], curviCoords[1],curviCoords[2]);
	//Compute Local Coordinates
	xh_loc = R*x0;
	Qloc = R*Ni*R.transpose();


}

void GNSSLS::computeSatPosParams(double tot, double tk, double a, double Deltan, double Mo, double e, double LilOmega, double Cuc, double Cus, double Crc, double Cis, double Cic, double Crs, double Idot, double I0, double BigOmegadot, double BigOmega, double toe)
{

	vector<double> tempSatPos(4);
	double n0, Mk, Ek, Vk, Ekd, Phik, dUk, dRk, dIk, Uk, Rk, x, y, n, Ik, BigOmegak, Xsat, Ysat, Zsat;
	//Constants
	double mu = 3.986005*pow(10, 14);	// Gravitational constant (m^3/s^2) used in mean motion computation
	const double we = 7.292115147e-5;	// Mean rotation rate of Earth (rad/s)
	n0 = sqrt((mu / pow(pow(a, 2), 3))); //Mean Motion
	n = n0 + Deltan; //Correction to mean motion
	Mk = Mo + n*tk; //Anomoly at time of transmission
	//Mk must be [0,2*PI]
	while (Mk >= 2.0*PI || Mk<0.0) {
		if (Mk>2.0*PI) {
			Mk -= 2.0*PI;
		}
		else {
			Mk += 2.0*PI;
		}

	}
	//itteratively Solve Eccentricity anomoly
	Ek = Mk;//set Ek initially
	for (int k = 0; k<6; k++) {
		Ek = Mk + e*sin(Ek); //Eccentricity
	}
	//Compute True anomoly
	Vk = atan2((sqrt(1 - pow(e, 2)) * sin(Ek)), (cos(Ek) - e));
	//Compute how much the eccentric anomoly differs
	Ekd = acos((e + cos(Vk)) / (1.0 + e*cos(Vk)));

	Phik = LilOmega + Vk;//Argument of latitude

	//compute the corrections to Keplerian Orbit
	dUk = Cuc*cos(2.0*Phik) + Cus*sin(2.0*Phik); //Along track corrections
	dRk = Crc*cos(2.0*Phik) + Crs*sin(2.0*Phik); //Across track corrections
	dIk = Cic*cos(2.0*Phik) + Cis*sin(2.0*Phik); //Across track corrections

												 //Compute corrected values
	Uk = Phik + dUk; //Corrected argument of latitude
	Rk = pow(a,2)*(1.0 - e*cos(Ek)) + dRk; // Corrected radius
	Ik = I0 + (Idot)*tk + dIk; // Corrected inclination

	//Correct longitude of ascending node
	BigOmegak = (BigOmega - we*toe) + (BigOmegadot - we)*tk; // Correct longitude of ascending node
	//BigOmegak = fmod(BigOmegak, 2.0*PI); //returns floating point remainder

	//Compute coordinates in orbital plane
	x = Rk*cos(Uk);
	y = Rk*sin(Uk);

	//compute coordinates in 
	Xsat = x*cos(BigOmegak) - y*sin(BigOmegak)*cos(Ik);
	Ysat = x*sin(BigOmegak) + y*cos(BigOmegak)*cos(Ik);
	Zsat = y*sin(Ik);
	tempSatPos[0] = tot;
	tempSatPos[1] = Xsat;
	tempSatPos[2] = Ysat;
	tempSatPos[3] = Zsat;
	calcSatPos.push_back(tempSatPos);


}

int GNSSLS::findNearestEpoch(double transtime)
{
	double ind = 0;
	double min = 1.0*pow(10,5);
	for (int i = 0; i < ephem.size(); i++) {

		//cout << abs(transtime - ephem[i][17]) << endl;
		if (abs(transtime - ephem[i][17]) < min) {
			//cout << abs(transtime - ephem[i][17]) << endl;
			min = abs(transtime - ephem[i][17]);
			ind = i;
		}
	}
	return ind;
}

void GNSSLS::computeSatPos(double interval, int duration, bool epochupdate)
{
	//Compute the position of the satellites based on the interval given (minutes)
	calcSatPos.clear();
	int counter = 0;
	int nearest_epoch = 0;
	double num_e, n, tk, a, n0, mu, Mk, Ek, Ekd, Vk, Phik, Uk, Rk, dUk, dRk, dIk, tot, IODE, BigOmega, BigOmegadot, Cic, Cis, Crc, Crs, Cuc, Cus, Deltan, e, Idot, Io, Mo, LilOmega, SqrtA, Toe;
	double x, y;
	vector<double> ephem_current; 
	vector<double> tempSatPos(2);
	//double n0, Mk, Ek, Vk, Ekd, Phik, dUk, dRk, dIk, Uk, Rk, x, y;
	//Constants
	mu = 3.986005*pow(10, 14);

	interval = interval * 60;//Into seconds
	num_e = ephem.size();
	int delimit = (duration*60*60) / interval;//How many intervals in between ephemeris published


	for (int j = 0; j < delimit; j++) {
		counter++;
		tot = (counter - 1)*interval;
		//if epochupdate, or first iteration, find nearest epoch
		if (epochupdate || j == 0) {
			nearest_epoch = findNearestEpoch(tot);
			//cout << "*********New Ephem\n";
			ephem_current = ephem[nearest_epoch];
			//Set ephem variables
			IODE = ephem_current[1];
			BigOmega = ephem_current[2];
			BigOmegadot = ephem_current[3];
			Cic = ephem_current[4];
			Cis = ephem_current[5];
			Crc = ephem_current[6];
			Crs = ephem_current[7];
			Cuc = ephem_current[8];
			Cus = ephem_current[9];
			Deltan = ephem_current[10];
			e = ephem_current[11];
			Idot = ephem_current[12];
			Io = ephem_current[13];
			LilOmega = ephem_current[14];
			Mo = ephem_current[15];
			a = ephem_current[16];
			Toe = ephem_current[17];
		}
		tk = tot - Toe;
		cout << "tk: " << tk << "Nearest Epoch: " << nearest_epoch << endl;
		computeSatPosParams(tot, tk, a, Deltan, Mo, e, LilOmega, Cuc, Cus, Crc, Cis, Cic, Crs, Idot, Io, BigOmegadot, BigOmega, Toe);

	}
}



MatrixXd GNSSLS::computeRotMatrix(double p, double l, double h)
{
	//-----Computes rotation matrix based on radians -------------
	MatrixXd R(4, 4);

	R(0, 0) = -1*sin(l);			R(0, 1) = cos(l);					R(0, 2) = 0;		R(0, 3) = 0;
	R(1, 0) = -1*sin(p)*cos(l);		R(1, 1) = -1 * sin(p)*sin(l);		R(1, 2) = cos(p) ;	R(1, 3) = 0;
	R(2, 0) = cos(p)*cos(l);		R(2, 1) = cos(p)*sin(l);			R(2, 2) = sin(p);	R(2, 3) = 0;
	R(3, 0) = 0;					R(3, 1) = 0;						R(3, 2) = 0;		R(3, 3) = 1;

	return R;

}

VectorXd GNSSLS::computeLatLong(double x, double y, double z)
{
	//----------------Compute latitude and longitude from cartesian coords
	//----------------returns in RADIANS


	VectorXd LatLong(3);
	double phi, lam, h, d, d1,d2, e, a, w, f, N,p ;
	//Reference ellipsoid coordinates
	a = 6378137.0;//m
	w = 7292115 * 10 ^ -11;// rad/s
	f = 1 / 298.257223563;
	e = sqrt(2*f-pow(f,2));
	p = sqrt(pow(x, 2) + pow(y, 2));

	//Lambda
	lam = atan2 (y , x);//radians

	d = 1;
	h = 0;
	N = 1;
	//iterate
	while (abs(d) > 0.0001) {
		d1 = N;
		phi = atan((z / p)*pow((1 - pow(e, 2)*(N / (N + h))), -1));//radians
		N = a / (sqrt(1 - pow(e, 2)*pow(sin(phi),2)));
		h = (p / cos(phi))-N;
		d2 = N;
		d = d1 - d2;
	}

	LatLong[0] = phi;
	LatLong[1] = lam;
	LatLong[2] = h;
	return LatLong;


}

void GNSSLS::computeDelta()
{
	//d = VectorXd::Zero(u);
	computeN();
	computeU();
	d = -1 * Ni * U;
}


void GNSSLS::computeClh()
{
	Clh = apost * A * Ni * A.transpose();
}

void GNSSLS::computeCvh()
{
	Cvh = apost * P.inverse() - Clh;

}

void GNSSLS::computeCxh()
{
	Cxh = apost * Ni;
}
void GNSSLS::computeDOP()
{
	xDOP = sqrt(Ni(0, 0));//qxx
	yDOP = sqrt(Ni(1, 1));//qyy
	zDOP = sqrt(Ni(2, 2));//qzz
	tDOP = sqrt(Ni(3, 3));//qdtdt
	pDOP = sqrt(pow(xDOP, 2) + pow(yDOP, 2) + pow(zDOP, 2));//position
	gDOP = sqrt(pow(xDOP, 2) + pow(yDOP, 2) + pow(zDOP, 2) + pow(tDOP, 2));//geometry
	hDOP = sqrt(Qloc(0,0)+Qloc(1,1));
	vDOP = sqrt(Qloc(2,2));

}
void GNSSLS::computefx()
{	
	fx = VectorXd::Zero(num_sats);
	for (int i = 0; i < num_sats; i++) {
		fx[i] = dist(x0[0], x0[1], x0[2], satPos(i, 1), satPos(i, 2), satPos(i, 3)) - x0[3];
	}
}

void GNSSLS::computeA()
{
	//Compute the design matrix A
	//The columns are:
	//|		xRec	|	yRec	|	zRec	|	dtRec	|
	//by number of satellites (Obesvations = equations)

	double satDist = 0;

	A = MatrixXd::Zero(num_sats, 4);
	for (int i = 0; i < num_sats; i++) {
		satDist = dist(x0(0), x0(1), x0(2), satPos(i,1), satPos(i,2), satPos(i,3));
		A(i, 0) = (x0(0) - satPos(i,1)) / satDist;//d/dxr
		A(i, 1) = (x0(1) - satPos(i,2)) / satDist;//d/dyr
		A(i, 2) = (x0(2) - satPos(i,3)) / satDist;//d/dyr
		A(i, 3) = -1;//d/dt
	}




}


void GNSSLS::updateCoords()
{

}