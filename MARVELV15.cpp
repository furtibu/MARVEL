// MARVELV15.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <Eigen/Eigen>
#include "tokenst.h"
#include <time.h>
#include <Windows.h>
#include <unordered_map>

using Eigen::SparseMatrixBase;
using Eigen::SparseView;
using namespace Eigen;
using namespace std;

typedef Eigen::SparseMatrix<long double> SpMat;
typedef Eigen::Matrix<long double, -1, 1> VectorLD;
typedef Eigen::Triplet<double> T;

int iZero = -1;
int iZeroSN1 = -1;
class Node;
class Edge;

int maxprec = 9;

#define MAX 10000

class Node
{
public:
	Node(string id, int nb)
		: id(id), mynb(nb), previous(NULL), distanceFromStart(MAX), energy(-1.0), visit(false), SNid(-1), avg(-1.00), avgunc(-1.0), trav_ind(0)
	{

	}

	void addtoTrs(Edge* trans)
	{
		myTrs.push_back(trans);
	}

	void addtoReference(string ref)
	{
		References.push_back(ref);
	}

	vector<Edge*>& getmyTR() { return myTrs; }

	int getNbofTrans() { return (int)myTrs.size(); }

	void SetEnergy(double ener) { energy = ener; }
	long double GetEnergy() { return energy; }
	bool HasEnergy() { if (energy == -1.0) { return false; } else { return true; } }

public:
	string id;
	string databaseid;
	string sym;
	int		mynb;
	int		nbinSN;
	Node* previous;
	int distanceFromStart;
	long double energy;
	long double testener;
	long double unc;
	long double uncRoland;
	long double unc2;
	long double avgunc;
	vector<Edge*> myTrs;
	vector<string> References;
	bool visit;
	long double avg;
	int SNid;
	int parity;
	int co_id;
	double distance;
	Edge* bond;
	Node* ancestor;
	int block_ind;
	int trav_ind;
	vector<Edge*> my_sp_trS;
};

bool firstNamePredicate(Node* a, Node* b)
{
	return a->GetEnergy() < b->GetEnergy();
}

bool SNnumber(Node* a, Node* b)
{
	return a->SNid < b->SNid;
}



class Edge
{
public:
	Edge(Node* node1, Node* node2, long double Freq, long double Unc, string Origunc, string Unit, string refi)
		: node1(node1), node2(node2), freq(Freq), unc(Unc), origunc(Origunc), gi(1.00 / Unc), Ref(refi), unc2(Unc*Unc),
		verybad(false), delnow(false), unit(Unit), mayDecrase(false), doDecrase(false)
	{


	}
	bool Connects(Node* node1, Node* node2)
	{
		return (
			(node1 == this->node1 &&
				node2 == this->node2) ||
				(node1 == this->node2 &&
					node2 == this->node1));
	}

	Node* getOtherNode(Node* node)
	{
		if (node == node1)
		{
			return node2;
		}
		else
		{
			return node1;
		}
	}

	bool AmIUp(Node* node)
	{
		if (node == node1)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	void SetDiff(long double difference) { diff = difference; }
	long double GetDiff() { return diff; }

public:
	Node* node1;
	Node* node2;
	int distance;
	long double freq;
	long double unc;
	long double testfreq;
	double detect;
	string origunc; 
	long double diff;
	long double gi;
	string Ref;
	string tag;
	long double unc2;
	bool verybad;
	bool inserted;
	bool explored;
	double weight;
	int block_ind;
	string ofreq;
	string ounc;
	string unit;
	bool delnow;
	int bindex;
	
	bool mayDecrase;
	bool doDecrase;
	int nbDigit;
	double optunc;

	
};

bool sortBindex(Edge* a, Edge* b)
{
	return a->bindex > b->bindex;
}

struct badMARVEL
{
	long double origunc;
	long double optimunc;
	string ref;
	string title;
};

void DFS(Node* node);
void DFS_MAP(Node* node, int SNid);
int MARVEL(vector<Node*> &energiesinSN, vector<Edge*> &transitions);
bool StringToInt(const string &s, int &i);
int mix(int a, int b, int c);

double get_wall_time() {
	LARGE_INTEGER time, freq;
	if (!QueryPerformanceFrequency(&freq)) {
		//  Handle error
		return 0;
	}
	if (!QueryPerformanceCounter(&time)) {
		//  Handle error
		return 0;
	}
	return (double)time.QuadPart / freq.QuadPart;
}


//
// Dijkstra algorithm to generate a shortest path spanning tree (SPT)
//
void dijkstra(vector<Node*> elS)
{
	struct dijkstra_aux
	{
		static void dijkstra_map(Node* start_el, int co_id, int &trav_ind)
		{
			Node *el, *neig_el;
			Edge* tr;

			priority_queue < pair<double, pair<int, Node*>>,
				vector<pair<double, pair<int, Node*>>>,
				greater < pair<double, pair<int, Node*>> >> prior_queue;

			int n;

			double sum;
			double distance;

			start_el->ancestor = NULL;
			start_el->bond = NULL;
			start_el->distance = 0.0;
			start_el->trav_ind = ++trav_ind;

			n = 0;

			prior_queue.push(make_pair(
				start_el->distance, make_pair(n, start_el)));

			while (prior_queue.size() != 0)
			{
				el = prior_queue.top().second.second;
				distance = prior_queue.top().first;
				prior_queue.pop();
				if (distance > el->distance) continue;
				el->trav_ind = ++trav_ind;
				el->co_id = co_id;
				for (int j = 0; j < (int)el->myTrs.size(); j++)
				{
					tr = el->myTrs[j];
					if (!tr->inserted) continue;
					neig_el = tr->getOtherNode(el);
					sum = el->distance + tr->weight;
					if (sum < neig_el->distance)
					{
						neig_el->distance = sum;
						neig_el->ancestor = el;
						neig_el->bond = tr;
						n++;
						prior_queue.push(make_pair(
							neig_el->distance, make_pair(n, neig_el)));
					}
				}
			}
		}
	};

	Node* el; //pointer for an EL
	Edge* tr; //pointer for a TR

	int co_id; //component identifier
	int trav_ind; //traversal index

	double max_dist; //maximum distance

	unordered_map<Node*, bool> um_viS; //unordered map to store
									 //whether the individual
									 //ELs are visited

									 //
									 // Determine the maximum distance from the root
									 //
	max_dist = -1.0;
	for (int i = 0; i < (int)elS.size(); i++)
	{
		el = elS[i];
		el->co_id = -1;
		for (int j = 0; j < (int)el->myTrs.size(); j++)
		{
			tr = el->myTrs[j];
			tr->explored = false;
			if (tr->inserted && tr->weight > max_dist)
				max_dist = tr->weight;
		}
	}
	max_dist *= (double)elS.size();

	//
	// Set the distances of the ELs
	// from their roots to maximal
	//
	for (int i = 0; i < (int)elS.size(); i++)
		elS[i]->distance = max_dist;

	//
	// Perform the Dijkstra traversal
	//
	trav_ind = co_id = -1;
	for (int i = 0; i < (int)elS.size(); i++)
	{
		el = elS[i];
		if (el->co_id == -1)
		{
			co_id++;
			dijkstra_aux::dijkstra_map(el, co_id, trav_ind);
		}
		else
		{
			el->bond->explored = true;
		}
	}
}

class Seg
{
public: 
	Seg(string iSeg, string iunit, double iesu) { tag = iSeg; unit = iunit; esu = iesu; prevesu = 0.0; calcesu = 0.0;}
	string tag;
	string unit;
	double esu;
	double calcesu;
	double prevesu;
	vector<double> diff;

	void CalcESU()
	{
		prevesu = calcesu;
		double avg = 0.0;
		double sdv = 0.0;
		for (auto db : diff)
		{
			avg += db;
		}

		avg /= diff.size();

		for (auto db : diff)
		{
			sdv += pow(db - avg, 2.0);
		}
		sdv /= (diff.size() - 1);
		sdv = sqrt(sdv);

		//cout << avg << "  " << sdv << endl;

		for (int i = 0; i < diff.size(); )
		{
			if (diff[i] >  avg + 3.0*sdv)
			{
				diff.erase(diff.begin() + i);
			}
			else
			{
				++i;
			}
		}

		avg = 0.0;
		sdv = 0.0;

		for (auto db : diff)
		{
			avg += db;
		}

		avg /= diff.size();

		for (auto db : diff)
		{
			sdv += pow(db - avg, 2.0);
		}
		sdv /= (diff.size() - 1);
		sdv = sqrt(sdv);

		calcesu = avg;

	}


};



int main()
{
	ofstream out;
	out.open("EnergyLevels.txt");
	out.setf(ios_base::fixed, ios_base::floatfield);

	ofstream outE;
	outE.open("EnergyLevelsInTransitions.txt");
	outE.setf(ios_base::fixed, ios_base::floatfield);

	ifstream data;
	data.open("transitions_HO37Cl.txt", ios::in); //transitions_H217O_W2020_Reass.txt

	if (!data) cerr << "Can't open transitions file!";

	ifstream segfile;
	segfile.open("segment_N2O_0414.txt", ios::in); //transitions_h216o_2017nov30.txt

	if (!segfile) cerr << "Can't open segments file!";

	ofstream outT;
	outT.open("CheckTransitions.txt");
	outT.setf(ios_base::fixed, ios_base::floatfield);

	ofstream outB;
	outB.open("BadTransitions.txt");
	outB.setf(ios_base::fixed, ios_base::floatfield);

	ofstream outN;
	outN.open("NewTransitions.txt");
	outN.setf(ios_base::fixed, ios_base::floatfield);

	ofstream outCl;
	outCl.open("Old_format.txt");
	outCl.setf(ios_base::fixed, ios_base::floatfield);


	vector<Node*> energies;
	vector<Node*> energiesInSN;
	vector<Edge*> transitions;

	long double freq, unc, sumwx, sumw, esu, csu, diff;
	string v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, origunc, sunc;
	string v1l, v2l, v3l, v4l, v5l, v6l, v7l, v8l, v9l, v10l, v11l, v12l, v13l;
	string Ref, assignup, assignlo, grade, unit, symup, symlow,  tempRef, vibup, viblow;
	bool haveup = false;
	bool havelow = false;
	bool mayDecrase;
	int ilow, iup, q;
	int id = 0;
	string line, sfreq, nfreq;
	vector<string> tokens;
	int linenb = 1;

	int uniq1up, uniq2up, uniqup;
	int uniq1lo, uniq2lo, uniqlo;
	typedef vector<int> IntContainer;
	typedef IntContainer::iterator IntIterator;
	IntContainer uniqlabels;
	int iv1, iv2, iv3, iv4, iv5, iv6;
	int iJup, ikaup, ikcup, iJlo, ikalo, ikclo;
	int ok = 0;
	int bindex, nDigit;

	map<string, Node*> map_ass;
	map<string, Node*>::iterator up_map_ass, low_map_ass, iter;
	map<string, Seg*> segments;
	map<string, Seg*>::iterator segmentsIter;
	
	Node* energyup;
	Node*	energylow;

	int NQN = 6;
	int minSize = 100;
	bool oldFormat = true;
	bool hasSegmentfile = false;

	if (hasSegmentfile)
	{
		while (!segfile.eof())
		{
			//segfile >> Ref >> unit >> esu >> freq >> freq >> freq >> freq >> freq >> freq >> freq >> freq;  //H2CO
			//segfile >> Ref >> unit >> freq >> freq >> freq >> freq >> freq >> freq;  // H216O
			segfile >> Ref >> unit;
			//segfile >> Ref >> unit>>esu>>csu;
			//cout << Ref << endl;

			//Ref.erase(remove(Ref.begin(), Ref.end(), '#'), Ref.end());

			Seg *seg = new Seg(Ref, unit, esu);
			segments[Ref] = seg;

			//segESU.insert(make_pair(Ref, vector<double>(1, 0.0)));
			//system("pause");
		}

		cout << "segment file reading is done!" << endl;
	}

	while (!data.eof())
	{

		getline(data, line);

		if (line.find("&") != std::string::npos) {
			continue;
		}
		
		tokens.clear();
		tokenize(tokens, line);
		int size = tokens.size();

		
		linenb++;

		if (oldFormat == false)
		{
			
			freq = atof(tokens[0].c_str());
			sfreq = tokens[0];
			nfreq = sfreq.substr(sfreq.find(".") + 1);
			nDigit = nfreq.length();
			
			origunc = tokens[1].c_str();
			sunc = tokens[2].c_str();
			//sunc.erase(remove_if(sunc.begin(), sunc.end(), not1(ptr_fun(::isdigit))), sunc.end());
			unc = atof(sunc.c_str());

			assignup = "";
			assignlo = "";
			vibup = "";
			viblow = "";

			for (int i = 3; i < (NQN + 3); i++)
				assignup += tokens[i] + " ";

			for (int i = 3; i < (NQN ); i++)
				vibup += tokens[i] + " ";

			for (int i = (NQN + 3); i < (2 * NQN + 3); i++)
				assignlo += tokens[i] + " ";

			for (int i = (NQN + 3); i < (2 * NQN ); i++)
				viblow += tokens[i] + " ";

			//cout << viblow << endl;

			if (unc == 0.00)
				unc = 0.000001;

		  //
			

			Ref = tokens[2 * NQN + 3];
			tempRef = Ref;
			size_t offset = tempRef.find(".");
			tempRef.erase(offset, 15);

			//cout << Ref << " " << endl;

			

			segmentsIter = segments.find(tempRef);

			unit = "cm-1";
			mayDecrase = true;

			if (segmentsIter != segments.end())
			{
				if (segmentsIter->second->unit == "MHz")
				{
					freq = freq*1E6 / 2.99792458E10;
					unc = unc*1E6 / 2.99792458E10;
					unit = "MHz";
					mayDecrase = false;
				}

				if (segmentsIter->second->unit == "THz")
				{
					freq = freq*1E12 / 2.99792458E10;
					unc = unc*1E12 / 2.99792458E10;
					unit = "THz";
					mayDecrase = false;
				}

				if (segmentsIter->second->unit == "GHz")
				{
					freq = freq*1E9 / 2.99792458E10;
					unc = unc*1E9 / 2.99792458E10;
					unit = "GHz";
					mayDecrase = false;
				}

				if (segmentsIter->second->unit == "kHz")
				{
					freq = freq*1E3 / 2.99792458E10;
					unc = unc*1E3 / 2.99792458E10;
					unit = "kHz";
					mayDecrase = false;
				}

				if (segmentsIter->second->unit == "Hz")
				{
					freq = freq / 2.99792458E10;
					unc = unc / 2.99792458E10;
					unit = "Hz";
					mayDecrase = false;
				}
			}
			else {
				cout << "Missing segment: " << tempRef << endl;
				system("pause");
				exit(-1);
			}

			bool convert = true;

			if (convert)
			{
				

				outCl << fixed << setprecision(8) << freq << " " << scientific << setprecision(3) << unc << " " << assignup << " " << assignlo << " " << Ref << endl;
			}

			

		}
		else {
			freq = atof(tokens[0].c_str());
			unc = atof(tokens[1].c_str());
			//sunc = tokens[2].c_str();
			//sunc.erase(remove_if(sunc.begin(), sunc.end(), not1(ptr_fun(::isdigit))), sunc.end());
			origunc = tokens[1].c_str();

			assignup = "";
			assignlo = "";

			for (int i = 2; i < (NQN + 2); i++)
				assignup += tokens[i] + " ";

			for (int i = (NQN + 2); i < (2 * NQN + 2); i++)
				assignlo += tokens[i] + " ";

			if (unc == 0.00)
				unc = 0.00000001;

			//

			//unc = 0.005;
		

			Ref = tokens[2 * NQN + 2];
			// cout << Ref << " " << endl;

			tempRef = Ref;
			size_t offset = tempRef.find(".");
			tempRef.erase(offset, 15);

			unit = "cm-1";

			if (hasSegmentfile)
			{

				segmentsIter = segments.find(tempRef);

				unit = "cm-1";

				if (segmentsIter != segments.end())
				{
					if (segmentsIter->second->unit == "MHz")
					{
						freq = freq*1E6 / 2.99792458E10;
						unc = unc*1E6 / 2.99792458E10;
						unit = "MHz";
					}

					if (segmentsIter->second->unit == "THz")
					{
						freq = freq*1E12 / 2.99792458E10;
						unc = unc*1E12 / 2.99792458E10;
						unit = "THz";
					}

					if (segmentsIter->second->unit == "GHz")
					{
						freq = freq*1E9 / 2.99792458E10;
						unc = unc*1E9 / 2.99792458E10;
						unit = "GHz";
					}

					if (segmentsIter->second->unit == "kHz")
					{
						freq = freq*1E3 / 2.99792458E10;
						unc = unc*1E3 / 2.99792458E10;
						unit = "kHz";
					}

					if (segmentsIter->second->unit == "Hz")
					{
						freq = freq / 2.99792458E10;
						unc = unc / 2.99792458E10;
						unit = "Hz";
					}
				}
				else {
					cout << "Missing segment: " << tempRef << endl;
					system("pause");
					exit(-1);
				}
			}

		}

		//out<<Ref<<endl;
		//cout<<Ref<<endl;

		if (freq < 0.0)
		{
			//outN << tokens[0] << "  " << origunc << "  " << origunc << "  " << assignup << "  " << assignlo << "  " << Ref << endl;
			outN<<tokens[0]<<"  "<<scientific<<setprecision(3)<<unc<<"  "<<assignup<<"  "<<assignlo<<"  "<<Ref<<endl;

			continue;
		}

		if (assignup == assignlo)
		{
			out << "problem in line:" << linenb - 1 << " " << assignup << " " << assignlo;
			exit(-1);
		}


		/*iv3 = atoi (tokens[4].c_str());
		ikaup = atoi (tokens[6].c_str());
		ikcup = atoi (tokens[7].c_str());
		q = iv3+ikaup+ikcup;

		if(q%2 == 0)
		{
		if(ikcup%2 == 0)
		{
		symup = "A1";
		} else {
		symup = "A2";
		}
		} else {
		if(ikcup%2 == 0)
		{
		symup = "B2";
		} else {
		symup = "B1";
		}

		}


		iv3 = atoi (tokens[10].c_str());
		ikaup = atoi (tokens[12].c_str());
		ikcup = atoi (tokens[13].c_str());
		q = iv3+ikaup+ikcup;

		if(q%2 == 0)
		{
		if(ikcup%2 == 0)
		{
		symlow = "A1";
		} else {
		symlow = "A2";
		}
		} else {
		if(ikcup%2 == 0)
		{
		symlow = "B2";
		} else {
		symlow = "B1";
		}

		}*/



		// up
		up_map_ass = map_ass.find(assignup);
		if (up_map_ass == map_ass.end())
		{
			energyup = new Node(assignup, id);
			energyup->sym = symup;
			id++;
			map_ass[assignup] = energyup;
			energies.push_back(energyup);
		}
		else
		{
			energyup = up_map_ass->second;
			energyup->sym = symup;
		}

		// low
		low_map_ass = map_ass.find(assignlo);
		if (low_map_ass == map_ass.end())
		{
			energylow = new Node(assignlo, id);
			energylow->sym = symlow;
			id++;
			map_ass[assignlo] = energylow;
			energies.push_back(energylow);
		}
		else
		{
			energylow = low_map_ass->second;
			energylow->sym = symlow;
		}

		Edge *trans = new Edge(energyup, energylow, freq, unc, origunc, unit, Ref);
		trans->ofreq = tokens[0];
		trans->ounc = tokens[1];
		trans->mayDecrase = mayDecrase;
		trans->nbDigit = nDigit;
		trans->tag = tempRef;


		trans->bindex = ceil(-1 * log10(unc));
		energyup->addtoTrs(trans);
		energylow->addtoTrs(trans);
		transitions.push_back(trans);

	}

	//for ( iter = map_ass.begin(); iter != map_ass.end(); iter++ )
	//{
	//	outE << iter->first  // string (key)
	//			 
	//			  << endl ;
	//}

	cout << "read is done!" << endl;

	bool havemarvel;
	string title;

	int SNid = 1;
	for (int i = 0; i < (int)energies.size(); i++)
	{
		if (energies[i]->visit == false)
		{
			DFS_MAP(energies[i], SNid);
			SNid++;
		}
	}

	//cout << SNid << endl;


	stable_sort(energies.begin(), energies.end(), SNnumber);

	cout << "DFS is done!" << endl;
	int northo = 0;
	for (int i = 1; i < SNid; i++)
	{
		if (energiesInSN.empty() == false)
			energiesInSN.clear();

		int idSN = 0;
		for (int j = 0; j < (int)energies.size(); j++)
		{
			if (energies[j]->SNid == i)
			{
				energies[j]->nbinSN = idSN;
				idSN++;
				energiesInSN.push_back(energies[j]);

			}
		}
		//cout<<"energiesInSN = "<<energiesInSN.size()<<endl;


		if (energiesInSN.size() >= minSize)
		{
			int iternumber = MARVEL(energiesInSN, transitions);

			sort(energiesInSN.begin(), energiesInSN.end(), firstNamePredicate);

			// final uncertainties

			Edge* tr;
			for (int i = 0; i < (int)transitions.size(); i++)
			{
				
				tr = transitions[i];
				tr->inserted = true;
				//cout << scientific << setprecision(3) << tr->unc << endl;
				tr->weight = tr->unc*tr->unc;
			}

			/*for (int j = 1; j < (int)energiesInSN.size(); j++)
			{
				if (energiesInSN[j]->energy == 0.00)
				{
					energiesInSN[j]->unc = 0.000000;
					energiesInSN[j]->uncRoland = 0.000000;
				}
			}*/

			Node* el;
			Node* act_el;
			energiesInSN[0]->unc = 0.0;
			energiesInSN[0]->uncRoland = 0.0;
			dijkstra(energiesInSN);
			for (int j = 1; j < (int)energiesInSN.size(); j++)
			{
				el = energiesInSN[j];
				el->uncRoland = sqrt(el->distance);
				//cout << el->id << " " << scientific << setprecision(3) << el->uncRoland << endl;
				act_el = el;
				while (act_el->ancestor != NULL)
				{
					
					tr = act_el->bond;
					
					el->my_sp_trS.push_back(tr);
					act_el = act_el->ancestor;
				}
			}

			

			long double uncchange;
			long double actenergy;
			double maxprecizitas = 0.0;
			int actbi;
			int bi;

			bool essential;
			int nbVeryBad = 0;
			for (int ii = 0; ii < energiesInSN.size(); ii++)
			{

				outT << ii + 1 << ")  " << energiesInSN[ii]->id << ">  " << energiesInSN[ii]->sym << "  "<<fixed << setprecision(8) << energiesInSN[ii]->energy << "  "<<scientific << setprecision(3) << energiesInSN[ii]->uncRoland << endl;
				
				nbVeryBad = 0;
				for (int j = 0; j < energiesInSN[ii]->myTrs.size(); j++)
				{
					if (energiesInSN[ii]->myTrs[j]->node1->energy != -1 && energiesInSN[ii]->myTrs[j]->node2->energy != -1)
					{
						diff = 0.0;

						maxprecizitas = sqrt(pow(energiesInSN[ii]->myTrs[j]->node1->uncRoland, 2.0) + pow(energiesInSN[ii]->myTrs[j]->node2->uncRoland, 2.0));
						
						essential = true;
						if (maxprecizitas >= energiesInSN[ii]->myTrs[j]->unc)
							essential = true;
						
						if (energiesInSN[ii]->myTrs[j]->AmIUp(energiesInSN[ii]))
						{
							actenergy = energiesInSN[ii]->myTrs[j]->freq + energiesInSN[ii]->myTrs[j]->getOtherNode(energiesInSN[ii])->energy;
						}
						else
						{
							actenergy = energiesInSN[ii]->myTrs[j]->getOtherNode(energiesInSN[ii])->energy - energiesInSN[ii]->myTrs[j]->freq;
						}

						diff = actenergy - energiesInSN[ii]->energy;

						if (essential)
						{
							outT << energiesInSN[ii]->myTrs[j]->Ref << "  " << fixed << setprecision(6)<< energiesInSN[ii]->myTrs[j]->freq << "  " << fixed << setprecision(8) << actenergy
								<< scientific << setprecision(3) << "  " << energiesInSN[ii]->myTrs[j]->unc << "  " << diff
								<< "   " << energiesInSN[ii]->myTrs[j]->getOtherNode(energiesInSN[ii])->id
								<< fixed << setprecision(8) << "  " << energiesInSN[ii]->myTrs[j]->getOtherNode(energiesInSN[ii])->energy << " "
								<< scientific << setprecision(3) << maxprecizitas;
						}
						else {
							outT << energiesInSN[ii]->myTrs[j]->Ref << "  " << fixed << setprecision(6) << energiesInSN[ii]->myTrs[j]->freq << "  " << fixed << setprecision(8) << actenergy
								<< scientific << setprecision(3) << "  " << energiesInSN[ii]->myTrs[j]->unc << "  " << diff
								<< "   " << energiesInSN[ii]->myTrs[j]->getOtherNode(energiesInSN[ii])->id
								<< fixed << setprecision(8) << "  " << energiesInSN[ii]->myTrs[j]->getOtherNode(energiesInSN[ii])->energy << " "
								<< scientific << setprecision(3) << maxprecizitas;
						}
						//energiesInSN[ii]->myTrs[j]->trialE = actenergy;

						if (diff < 0.0)
							diff *= -1;

						energiesInSN[ii]->myTrs[j]->optunc = diff;

						if (energiesInSN[ii]->myTrs.size() > 3 && diff < energiesInSN[ii]->myTrs[j]->unc  && energiesInSN[ii]->myTrs[j]->mayDecrase == true)
							energiesInSN[ii]->myTrs[j]->doDecrase = true;


						uncchange = abs(diff - energiesInSN[ii]->myTrs[j]->unc);

						if (abs(diff) > energiesInSN[ii]->myTrs[j]->unc ) //  && diff > 0.000001
						{
							//outB << energiesInSN[ii]->myTrs[j]->Ref << "  " << diff / energiesInSN[ii]->myTrs[j]->unc << endl;

							double ratio = diff / energiesInSN[ii]->myTrs[j]->unc;

							if (ratio > 10)
							{
								if (ratio > 100)
								{
									if (ratio > 1000)
									{
										if (diff > 1.0)
										{
											if (diff > 100.0)
											{
												outT << " VERY DEL_100";
												nbVeryBad++;
											}
											else {
												if (diff > 10.0)
												{
													outT << " VERY DEL_100";
													nbVeryBad++;
												}
												else {
													outT << " VERY DEL_10";
													nbVeryBad++;
												}
											}
										}
										else {
											outT << " VERY FUCK_1000";
											nbVeryBad++;
										}
									}
									else {
										outT << " VERY FUCK_100";
										nbVeryBad++;
									}
								}
								else {
									outT << " VERY BAD";
									nbVeryBad++;
								}
							}
							else {
								outT << " BAD";
							}

							actbi = abs((int)floor(log10(diff)));
							bi = abs((int)floor(log10(energiesInSN[ii]->myTrs[j]->unc)));

							if((bi - actbi) >1)
								outT << " CRITICAL ";

							if (diff > 0.05)
							{
								outT << " ERROR";
							}
							outT << " " << endl;
						}
						else
						{
							outT << endl;
						}

					}


					//
				} // myTrans

				  //outT<<energiesInSN[ii]->myTrs.size() <<" "<<nbVeryBad<<endl<<endl;
				if (energiesInSN[ii]->myTrs.size() == nbVeryBad)
				{
					outT << "ALL VERY BAD " << endl << endl;
				}

				
				outT << endl;
			} // enenrgies


			for (int j = 0; j < (int)energiesInSN.size(); j++)
			{
				out << energiesInSN[j]->id << "\t" << setprecision(maxprec) << energiesInSN[j]->energy << "\t" << setprecision(maxprec) << "\t" << energiesInSN[j]->uncRoland << "\t" << energiesInSN[j]->getNbofTrans() << endl; // <<energiesInSN[j]->unc

			}

		}

	}



	//
	// NEW TRANSITIONS
	//

	long double newunc = 0;
	long double actdiff, unc1;

	for (int i = 0; i < (int)transitions.size(); i++)
	{
		if (transitions[i]->node1->energy != -1 && transitions[i]->node2->energy != -1)
		{
			newunc = transitions[i]->unc;


			if (transitions[i]->AmIUp(transitions[i]->node1))
			{
				actdiff = abs(transitions[i]->freq - (transitions[i]->node1->energy - transitions[i]->node2->energy));
			}
			else {
				actdiff = abs(transitions[i]->freq - (transitions[i]->node2->energy - transitions[i]->node1->energy));
			}

			/*if (actdiff < 0.1*transitions[i]->unc && actdiff > 0.05)   //
			{
				//newunc = 1.1 * actdiff;
				outB << transitions[i]->ofreq << "  " << transitions[i]->origunc << "  " << scientific << setprecision(3) << actdiff << "  " << transitions[i]->node1->id << "  " << transitions[i]->node2->id << "  " << transitions[i]->Ref << endl;

				//cout << "van" << endl;
			}*/


			// Increase

			if (actdiff > transitions[i]->unc && actdiff > 0.01)   //
			{
				newunc = 1.1 * actdiff;
				outB<< transitions[i]->ofreq << "  " << transitions[i]->origunc << "  " << scientific << setprecision(3) << newunc << "  " << transitions[i]->node1->id << "  " << transitions[i]->node2->id << "  " << transitions[i]->Ref << endl;
				//std::cout << '\7';
				cout << "van" << endl;
			}

			
			unit = transitions[i]->unit;

			if (unit == "MHz")
				newunc = newunc*2.99792458E10 / 1E6;

			if (unit == "THz")
				newunc = newunc*2.99792458E10 / 1E12;

			if (unit == "GHz")
				newunc = newunc*2.99792458E10 / 1E9;

			if (unit == "kHz")
				newunc = newunc*2.99792458E10 / 1E3;

			if (unit == "Hz")
				newunc = newunc*2.99792458E10;


			int deleteTRs = 0;

			//
			// delete transitions
			//

			//outN << transitions[i]->ofreq << "  " << transitions[i]->origunc << "  " << scientific << setprecision(3) << newunc << "  " << transitions[i]->node1->id << "  " << transitions[i]->node2->id << "  " << transitions[i]->Ref << endl;
			outN<<transitions[i]->ofreq<<"  " << scientific << setprecision(3) << newunc <<"  "<<transitions[i]->node1->id<<"  "<<transitions[i]->node2->id<<"  "<<transitions[i]->Ref<<endl;
		}
		else {
			//outN << transitions[i]->ofreq << "  " << transitions[i]->origunc << "  " << transitions[i]->ounc << "  " << transitions[i]->node1->id << "  " << transitions[i]->node2->id << "  " << transitions[i]->Ref << endl;
			outN<<transitions[i]->ofreq<<"  " << scientific << setprecision(3) <<transitions[i]->origunc<<"  "<<transitions[i]->node1->id<<"  "<<transitions[i]->node2->id<<"  "<<transitions[i]->Ref<<endl;
		}

	}


	//
	// transitios for Gephi
	//

	/*ofstream outG;
	outG.open("GephiTRs.csv");
	outG.setf(ios_base::fixed,ios_base::floatfield);

	for(int i = 0; i < (int)transitions.size(); i++)
	{
	outG<<transitions[i]->node1->mynb<<";"<<transitions[i]->node2->mynb<<endl;
	}*/


}

int MARVEL(vector<Node*> &energiesInSN, vector<Edge*> &transitions)
{
	int nbEner = (int)energiesInSN.size();

	VectorLD Y(nbEner), x(nbEner); //, diagonal(nbEner);
	SpMat A(nbEner, nbEner);

	//A.reserve(VectorXi::Constant(nbEner,nbEner));

	//MatrixXd matest(nbEner,nbEner);

	vector<Node* >::iterator itEL;
	vector<Edge* >::iterator itTR;

	int row = 0, col = 0;
	long double diag = 0.00;
	long double sum = 0.00, sumtr = 0.00, uncsqr, tr;
	int nnz = 0;
	int i = 0;
	bool amiup;
	Node* actenergy;
	Edge* acttr;

	cout << "MATRIX building..." << endl;
	double wall10 = get_wall_time();

	int nonzero = 0;
	int rownb = 0;
	int maxerownb = 0;

	for (itEL = energiesInSN.begin(); itEL != energiesInSN.end(); itEL++)
	{
		actenergy = *itEL;
		rownb = 0;
		for (itTR = actenergy->myTrs.begin(); itTR != actenergy->myTrs.end(); itTR++)
		{
			rownb++;
			nonzero++;
		}
		if (maxerownb < rownb)
			maxerownb = rownb;

		nonzero++;
	}

	A.reserve(VectorXi::Constant(nbEner, maxerownb));

	//for(int i = 0; i < nbEner; i++)
	for (itEL = energiesInSN.begin(); itEL != energiesInSN.end(); itEL++)
	{
		actenergy = *itEL;
		Y[i] = 0.00;
		x[i] = 0.00;
		diag = 0.00;
		row = actenergy->nbinSN;
		sumtr = 0.00;
		//for(int j = 0; j < (int) actenergy->myTrs.size(); j++)
		for (itTR = actenergy->myTrs.begin(); itTR != actenergy->myTrs.end(); itTR++)
		{

			acttr = *itTR;
			uncsqr = acttr->unc*acttr->unc;
			Node* neigenergy = acttr->getOtherNode(actenergy);

			col = neigenergy->nbinSN;

			A.coeffRef(row, col) += -1.00 / (uncsqr);

			sum += -1.00 / uncsqr;  //*actenergy->myTrs[j]->unc

			amiup = acttr->AmIUp(actenergy);
			diag += 1.00 / uncsqr;  // *actenergy->myTrs[j]->unc
			sum += -1.00 / uncsqr;  // *actenergy->myTrs[j]->unc

			tr = acttr->freq / uncsqr;  //*actenergy->myTrs[j]->unc
			if (!amiup)
				tr *= -1.00;
			sumtr += tr;
		}
		Y[i] += sumtr;
		
		A.insert(row, row) = diag;
		i++;

	}

	//cout << "sum = " << Y[0] << endl;

	double wall11 = get_wall_time();
	cout << "Wall Time = " << wall11 - wall10 << endl;

	double wall0 = get_wall_time();

	SimplicialLDLT<Eigen::SparseMatrix<long double> > chol(A);

	if (chol.info() != Eigen::Success)
	{
		cout << " wrong decomp "<<chol.info() << endl;
		system("pause");
		exit;
	}
	x = chol.solve(Y);
	if (chol.info() != Eigen::Success)
	{
		cout << " wrong solution " << chol.info() << endl;
		system("pause");
		exit;
	}


	double wall1 = get_wall_time();
	cout << "Wall Time = " << wall1 - wall0 << endl;
	
	long double minimum = x[0];

	for (int i = 0; i < nbEner; i++)
	{
		if (x[i] < minimum)
			minimum = x[i];
	}


	long double wavg1, wavg2, actunc;

	for (int i = 0; i < nbEner; i++)
	{
		if (minimum < 0.0)
		{
			energiesInSN[i]->SetEnergy(x[i] + abs(minimum));
		}
		else
		{
			energiesInSN[i]->SetEnergy(x[i] - minimum);
		}


	}

	//stable_sort(energiesInSN.begin(), energiesInSN.end(), firstNamePredicate);


	// erzekenyseg analizis
	//vector<double> sdiff;
	//vector<Edge*> sbad;

	bool runSensi = false;


	if (runSensi)
	{

		ofstream outS;
		outS.open("Erzekenyseg.txt");
		outS.setf(ios_base::fixed, ios_base::floatfield);

		int SNid = energiesInSN[0]->SNid;
		double actdiff, mydiff;

		for (auto atr : transitions)
		{
			atr->testfreq = 0.000;
		}

		for (auto atr : transitions)
		{
			if (atr->node1->SNid != SNid)
				continue;

			atr->testfreq = 10.0*atr->unc;
			i = 0;
			for (itEL = energiesInSN.begin(); itEL != energiesInSN.end(); itEL++)
			{
				actenergy = *itEL;
				Y[i] = 0.00;
				sumtr = 0.00;
				for (itTR = actenergy->myTrs.begin(); itTR != actenergy->myTrs.end(); itTR++)
				{
					acttr = *itTR;
					uncsqr = acttr->unc*acttr->unc;
					amiup = acttr->AmIUp(actenergy);

					tr = acttr->testfreq / uncsqr;
					if (!amiup)
						tr *= -1.00;
					sumtr += tr;
				}
				Y[i] += sumtr;

				i++;

			}

			x = chol.solve(Y);
			if (chol.info() != Eigen::Success)
			{
				cout << " wrong solution " << chol.info() << endl;
				system("pause");
				exit;
			}

			long double minimum = x[0];

			for (int i = 0; i < nbEner; i++)
			{
				if (x[i] < minimum)
					minimum = x[i];
			}

			for (int i = 0; i < nbEner; i++)
			{
				if (minimum < 0.0)
				{
					energiesInSN[i]->testener = x[i] + abs(minimum);
				}
				else
				{
					energiesInSN[i]->testener = x[i] - minimum;
				}
			}



			//sdiff.clear();
			//sbad.clear();
			//for (auto btr : transitions)
			//{
			//	if (btr->node1->SNid != SNid)
			//		continue;

			//	if (btr->node1->testener != -1 && btr->node2->testener != -1)
			//	{

			//		if (btr->AmIUp(btr->node1))
			//		{
			//			actdiff = abs(btr->testfreq - (btr->node1->testener - btr->node2->testener));
			//		}
			//		else {
			//			actdiff = abs(btr->testfreq - (btr->node2->testener - btr->node1->testener));
			//		}

			//		if (btr == atr)
			//			mydiff = actdiff;

			//		if (actdiff > btr->unc)
			//		{
			//			sdiff.push_back(actdiff);
			//			sbad.push_back(btr);
			//			//outS << btr->Ref << " " << scientific << setprecision(3) << btr->unc << " " << actdiff << endl;
			//		}
			//	}
			//}

			//if (sbad.empty() == true)
			//{
			//	atr->detect = 0.0;
			//}
			//else {
			//	atr->detect = mydiff * 100 / (10.0*atr->unc);
			//}

			if (atr->AmIUp(atr->node1))
			{
				actdiff = abs(atr->testfreq - (atr->node1->testener - atr->node2->testener));
			}
			else {
				actdiff = abs(atr->testfreq - (atr->node2->testener - atr->node1->testener));
			}

			atr->detect = actdiff * 100 / (10.0*atr->unc);

			outS << atr->Ref << " " << fixed << setprecision(3) << atr->detect << endl;

			atr->testfreq = 0.00;
		}

	}



	return 0;
}

void DFS(Node* node)
{
	node->visit = true;
	for (int i = 0; i < node->getNbofTrans(); i++)
	{
		if (node->myTrs[i]->getOtherNode(node)->visit == false)
			DFS(node->myTrs[i]->getOtherNode(node));
	}
}

void DFS_MAP(Node* node, int SNid)
{
	node->visit = true;
	node->SNid = SNid;
	for (int i = 0; i < node->getNbofTrans(); i++)
	{
		if (node->myTrs[i]->getOtherNode(node)->visit == false)
			DFS_MAP(node->myTrs[i]->getOtherNode(node), SNid);
	}
}

bool StringToInt(const string &s, int &i)
{
	istringstream myStream(s);

	if (myStream >> i)
		return true;
	else
		return false;
}

int mix(int a, int b, int c)
{
	a = a - b;  a = a - c;  a = a ^ (c >> 13);
	b = b - c;  b = b - a;  b = b ^ (a << 8);
	c = c - a;  c = c - b;  c = c ^ (b >> 13);
	a = a - b;  a = a - c;  a = a ^ (c >> 12);
	b = b - c;  b = b - a;  b = b ^ (a << 16);
	c = c - a;  c = c - b;  c = c ^ (b >> 5);
	a = a - b;  a = a - c;  a = a ^ (c >> 3);
	b = b - c;  b = b - a;  b = b ^ (a << 10);
	c = c - a;  c = c - b;  c = c ^ (b >> 15);
	return c;
}


