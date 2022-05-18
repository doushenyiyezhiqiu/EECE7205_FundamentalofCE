#include<iostream>
#include<vector>
#include<limits.h>
#include"Matrix.h"
#include<ctime>

using namespace std;
using namespace Numeric_lib;

struct task {
	int num; // the number of the task
	bool ct;  // judge whether the task is a cloud task
	double pri; // 
	int FTl;  // the finish time of the task in a core
	int FTWS; // the finish time of the task in sending
	int FTC;   // the finish time of the task in cloud
	int FTWR;   // the finish time of the task in receiving
	int RTl; // the earlist  time that the task can start in local core
	int RTWS; // the earlist  time that the task can start in sending
	int RTC;  // the earlist  time that the task can start in cloud
	int RTWR; //the earlist  time that the task can start in receiving
	int ST; // the task's actual start time
	int chan; // illustrate which channel the task operate (local core = 0,1,2, cloud=3)
	bool exit;  //whether it is an exit task
	bool entry;  // whether it is an entry task
	int ready1;
	int ready2;
};

// Phase one in step one: primary assignment
void primary(vector<task>&ini, Matrix<int, 2>&ta,int t)
{
	int min;
	unsigned int i;
	unsigned int j;
	for (i = 0; i < ta.dim1(); i++)
	{
		ini[i].num = i + 1;
		min = ta(i, 0);
		for (j = 0; j < ta.dim2(); j++)
			if (ta(i, j) < min)
				min = ta(i, j);
		if (min > t)
			ini[i].ct = 1;
		else
			ini[i].ct = 0;
	}
}

// Phase two in step one: task prioritizing 
void prioritize(vector<task>&ini, Matrix<int, 2>&ta, Matrix<int, 2>&G,int t)
{
	unsigned int i;
	unsigned int j,m;
	int k ;
	double w;
	double max;
	m = ini.size() - 1;
	for (i = 0; i < ini.size(); i++)
	{
		k = 0;
		for (j = 0; j < G.dim2(); j++)
			if (G(m - i, j) == 1)
				k = k + 1;
		if (k == 0)
			ini[m - i].exit = 1;
		k = 0;
		for (j = 0; j < G.dim2(); j++)
			if (G(j, m - i) == 1)
				k = k + 1;
		if(k==0)
			ini[m - i].entry = 1;
		max = 0;
		w = 0;
		if (!(ini[m - i].ct))
		{
			for (j = 0; j < ta.dim2(); j++)
				w = w + ta(m - i, j);
			w = w / 3;
		}
		else
			w = t;
		for (j = 0; j < G.dim2(); j++)
			if ((G(m - i, j) == 1) && (max < ini[j].pri))
				max = ini[j].pri;
		ini[m - i].pri = w + max;
	}
}

int find_biggest_pri(vector<task>&ini)
{
	unsigned int i;
	int max=0;
	for (i = 0; i < ini.size(); i++)
		if (ini[i].pri > ini[max].pri)
			max = i;
	return max;
}

// find the max in two numbers
int max2(int &m, int &n)
{
	if (m >= n)
		return m;
	else
		return n;
}

// if local schedule, return RTL
int d_rtl(task &vi, vector<task>&S, Matrix<int, 2>&G) 
{
	unsigned int i;
	unsigned int j;
	int max=0;
	if (S.size()!=0)
	{
		for (i = 0; i < G.dim2(); i++)
			if (G(i, vi.num - 1) == 1)
				for (j = 0; j < S.size(); j++)
					if ((S[j].num == i + 1)&&(max < max2(S[j].FTl, S[j].FTWR)))
						max = max2(S[j].FTl, S[j].FTWR);
	}
	return max;
}

// if cloud schedule, return RTWS
int d_rtws(task &vi, vector<task>&S, Matrix<int, 2>&G)
{
	unsigned int i;
	unsigned int j;
	int max=0;
	if (S.size()!=0)
	{
		for (i = 0; i < G.dim2(); i++)
			if (G(i, vi.num - 1) == 1)
				for (j = 0; j < S.size(); j++)
					if ((S[j].num == i + 1)&&(max < max2(S[j].FTl, S[j].FTWS)))
						max = max2(S[j].FTl, S[j].FTWS);
	}
	return max;
}

// if cloud schedule, return RTC
int d_rtc(task &vi, vector<task>&S, Matrix<int, 2>&G)
{
	unsigned int i;
	unsigned int j;
	int max=vi.FTWS;
	if (S.size()!=0)
	{
		for (i = 0; i < G.dim2(); i++)
			if (G(i, vi.num - 1) == 1)
				for (j = 0; j < S.size(); j++)
					if ((S[j].num == i + 1)&&(max < max2(vi.FTWS, S[j].FTC)))
						max = max2(vi.FTWS, S[j].FTC);
	}
	return max;
}

// if cloud schedule, return RTWR
int d_rtwr(task &vi)
{
	return vi.FTC;
}

// if local schedule, return the smallest finish time
int locals(task &vi, vector<task>&S, Matrix<int, 2>&G, Matrix<int, 2>&ta)
{
	vi.RTl = d_rtl(vi, S, G);
	unsigned int i;
	unsigned int j;
	int mint=INT_MAX;
	int ft;
	int max = 0; // find a local core's biggest finish time
	if (S.size()==0)
	{
		for (i = 0; i < ta.dim2(); i++)
		{
			ft = ta(vi.num - 1, i);
			if (mint > ft)
			{
				mint = ft;
				vi.chan = i;
			}
		}
		return mint;
	}
	for (i = 0; i < ta.dim2(); i++)
	{
		ft = vi.RTl + ta(vi.num - 1, i);
		max = 0;
		for (j = 0; j < S.size(); j++)
			if ((S[j].chan == i) && (max < S[j].FTl))
				max = S[j].FTl;
		if(max>vi.RTl)
			ft=max+ ta(vi.num - 1, i);
		if (mint > ft)
		{
			mint = ft;
			vi.chan = i;
		}
	}
	return mint;
}

// if cloud schedule, return the finish time
int clouds(task &vi, vector<task>&S, Matrix<int, 2>&G, int ts, int tc, int tr)
{
	vi.RTWS = d_rtws(vi, S, G);
	unsigned int i;
	int maxs = 0;
	int t;
	int maxc = 0;
	int maxr = 0;
	int ft;
	t = ts + tc + tr;
	if (S.size()==0)
	{
		vi.FTWS = ts;
		vi.RTC = ts;
		vi.FTC = ts + tc;
		vi.RTWR = ts+tc;
		return t;
	}
	for(i=0;i<S.size();i++)
		if (S[i].chan == 3)
			if (maxs < S[i].FTWS)
				maxs = S[i].FTWS;
	if (maxs > vi.RTWS)
		vi.FTWS = maxs + ts;
	else
		vi.FTWS = vi.RTWS + ts;
	vi.RTC = d_rtc(vi, S, G);
	for (i = 0; i < S.size(); i++)
		if (S[i].chan == 3)
			if (maxc < S[i].FTC)
				maxc = S[i].FTC;
	if (maxc > vi.RTC)
		vi.FTC = maxc + tc;
	else
		vi.FTC = vi.RTC + tc;
	vi.RTWR = d_rtwr(vi);
	for (i = 0; i < S.size(); i++)
		if (S[i].chan == 3)
			if (maxr < S[i].FTWR)
				maxr = S[i].FTWR;
	if (maxr > vi.RTWR)
		ft = maxr + tr;
	else
		ft = vi.RTWR + tr;
	return ft;
}

void initials(vector<task>&S, vector<task>&ini, Matrix<int, 2>&ta, Matrix<int, 2>&G, int ts, int tc, int tr)
{
	unsigned int i;
	int t;
	int maxp; // find the max priority in each iteration of ini
	int mint;  // find the minimum finish time of local
	int anot;  // perpare for another time (cloud)
	t = ts + tc + tr;
	for (i = 0; i < G.dim1(); i++)
	{
		maxp = find_biggest_pri(ini);
		if (!ini[maxp].ct)
		{
			mint = locals(ini[maxp], S, G, ta);
			anot = clouds(ini[maxp], S, G, ts, tc, tr);
			if (anot < mint)
			{
				ini[maxp].RTl = 0;
				ini[maxp].FTl = 0;
				ini[maxp].chan = 3;
				ini[maxp].FTWR = anot;
				ini[maxp].ST = anot - t;
			}
			else
			{
				ini[maxp].FTC = 0;
				ini[maxp].FTWS = 0;
				ini[maxp].RTWS = 0;
				ini[maxp].RTC = 0;
				ini[maxp].RTWR = 0;
				ini[maxp].FTWR = 0;
				ini[maxp].FTl = mint;
				ini[maxp].ST = mint - ta(ini[maxp].num - 1, ini[maxp].chan);
			}
		}
		else
		{
			ini[maxp].FTl = 0;
			ini[maxp].RTl = 0;
			ini[maxp].chan = 3;
			ini[maxp].FTWR= clouds(ini[maxp], S, G, ts, tc, tr);
			ini[maxp].ST = ini[maxp].FTWR - t;
		}
		S.push_back(ini[maxp]);
		ini.erase(ini.begin() + maxp);
	}
}

// return a task's finish time
int find_ft(task&vi)
{
	int max;
	max = max2(vi.FTl, vi.FTWR);
	return max;
 }

// print the sequence S
void prints(vector<task>&S)
{
	unsigned int i;
	int k,m;
	for (i = 0; i < S.size(); i++)
	{
		k = 1 + S[i].chan;
		m = find_ft(S[i]);
		cout << "Task" << S[i].num << ": ";
		switch (S[i].chan)
		{
		case 0:
			cout << "local core" << k << ", ";
			break;
		case 1:
			cout << "local core" << k << ", ";
			break;
		case 2:
			cout << "local core" << k << ", ";
			break;
		case 3:
			cout << "cloud" << ", ";
			break;
		default:
			break;
		}
		cout << "start time is: " << S[i].ST << ",finish time is: "<<m<<endl;
	}
}

// return the completion time of sequence S
int find_tcom(vector<task>&S)
{
	unsigned int i;
	int max=0;
	for (i = 0; i < S.size(); i++)
		if ((S[i].exit) && (max < find_ft(S[i])))
			max = find_ft(S[i]);
	return max;
}

// return the total energy of the sequence S
double find_en(vector<task>&S, int p1, int p2, int p3, double ps)
{
	unsigned int i;
	double ene=0;
	for (i = 0; i < S.size(); i++)
	{
		switch (S[i].chan)
		{
		case 0:
			ene = ene + p1 * (find_ft(S[i]) - S[i].ST);
			break;
		case 1:
			ene = ene + p2 * (find_ft(S[i]) - S[i].ST);
			break;
		case 2:
			ene = ene + p3 * (find_ft(S[i]) - S[i].ST);
			break;
		case 3:
			ene = ene + ps * (S[i].FTWS - S[i].ST);
			break;
		default:
			break;
		}
	}
	return ene;
}

//compute all the ready1 in a sequence
void get_ready1(vector<task>&S, Matrix<int, 2>&G)
{
	unsigned int i, j, k;
	int m;
	for (i = 0; i < S.size(); i++)
	{
		m = 0;
		for (j = 0; j < G.dim2(); j++)
			if (G(j, S[i].num-1) == 1)
				for (k = 0; k < S.size(); k++)
					if (S[k].num == j + 1)
						m = m + 1;
		S[i].ready1 = m;
	}
}

//compute all the ready2 in a sequence
void get_ready2(vector<task>&S)
{
	unsigned int i, j;
	int m;
	for (i = 0; i < S.size(); i++)
	{
		m = 0;
		for (j = 0; j < S.size(); j++)
			if ((S[i].chan == S[j].chan) && (S[j].ST < S[i].ST))
				m = m + 1;
		S[i].ready2 = m;
	}
}

// local schedule task vi whose local core is confirmed
int localse(task &vi, vector<task>&SN, Matrix<int, 2>&G, Matrix<int, 2>&ta)
{
	vi.RTl = d_rtl(vi, SN, G);
	unsigned int i;
	int ft;
	int max=0;
	if (SN.size()==0)
		ft = vi.RTl + ta(vi.num - 1, vi.chan);
	else
	{
		for (i = 0; i < SN.size(); i++)
			if ((SN[i].chan == vi.chan) && (max < SN[i].FTl))
				max = SN[i].FTl;
		if(max>vi.RTl)
			ft=max+ ta(vi.num - 1, vi.chan);
		else
			ft=vi.RTl+ ta(vi.num - 1, vi.chan);
	}
	return ft;
}

void kernel(vector<task>&S, vector<task>&SN, int ktar, task vtar, Matrix<int, 2>&G, Matrix<int, 2>&ta,int ts, int tc, int tr)
{
	unsigned int i;
	int m;
	int t;
	t = ts + tc + tr;
	vector<task>re;
	re = S;
	for (i = 0; i < re.size(); i++)
		if (vtar.num == re[i].num)
		{
			re[i].chan = ktar;
			if (ktar == 3)
			{
				re[i].FTl = 0;
				re[i].RTl = 0;
			}
		}
	while (re.size()!=0)
	{
		get_ready1(re, G);
		get_ready2(re);
		m = 0;
		while ((re[m].ready1!=0)&&(re[m].ready2 != 0))
			m = m + 1;
		if (re[m].chan == 3)
		{
			re[m].FTWR = clouds(re[m], SN, G, ts, tc, tr);
			re[m].ST = re[m].FTWR - t;
		}
		else
		{
			re[m].FTl = localse(re[m], SN, G, ta);
			re[m].ST = re[m].FTl - ta(re[m].num - 1, re[m].chan);
		}
		SN.push_back(re[m]);
		re.erase(re.begin() + m);
	}
}

void mcc(vector<task>&S, Matrix<int, 2>&G, Matrix<int, 2>&ta, int ts, int tc, int tr,int p1, int p2, int p3, double ps, int tmax)
{
	unsigned int i, j;
	int tcom;
	int tcom2;
	int a;
	double en;
	double en1;
	double en2;
	double ratio1=0;
	double ratio2;
	vector<task>SN;
	tcom = find_tcom(S);
	en = find_en(S, p1, p2, p3, ps);
	for (i = 0; i < S.size(); i++)
	{
		a = S[i].chan;
		if (S[i].chan != 3)
		{
			for (j = 0; j < 4; j++)
			{
				if (j != a)
				{
					SN.erase(SN.begin(), SN.end());
					en1 = find_en(S, p1, p2, p3, ps);
					kernel(S, SN, j, S[i], G, ta, ts, tc, tr);
					tcom2 = find_tcom(SN);
					en2 = find_en(SN, p1, p2, p3, ps);
					if ((en2 < en1) && (tcom >= tcom2))
						S = SN;
					else if ((en2 < en1) && (tcom2 <= tmax))
					{
						ratio2 = (en - en2) / (tcom2 - tcom);
						if (ratio2 > ratio1)
						{
							ratio1 = ratio2;
							S = SN;
						}
					}
				}
			}
		}
	}
}

void outerloop(vector<task>&S, Matrix<int, 2>&G, Matrix<int, 2>&ta, int ts, int tc, int tr, int p1, int p2, int p3, double ps, int tmax)
{
	double en;
	double en1=0;
	en = find_en(S, p1, p2, p3, ps);
	while (en1<en)
	{
		en= find_en(S, p1, p2, p3, ps);
		mcc(S, G, ta, ts, tc, tr, p1, p2, p3, ps, tmax);
		en1= find_en(S, p1, p2, p3, ps);
	}
}

int main()
{
	int N1 = 10;  // the number of tasks
	int K = 3;   // the number of local cores
	int ts1 = 3,tc1 = 1,tr1 = 1;
	unsigned int i;
	int t1;
	int tmax1 = 27;
	int p11 = 1;
	int p12 = 2;
	int p13 = 4;
        double ps1 = 0.5;
        double rt;
        clock_t start, end;
	t1 = ts1 + tc1 + tr1;
	Matrix<int, 2>G1(N1,N1);
	Matrix<int, 2>ta1(N1,K);
	vector<task>ini1(N1);
	vector<task>S1;
	G1(0, 1) = 1;
	G1(0, 2) = 1;
	G1(0, 3) = 1;
	G1(0, 4) = 1;
	G1(0, 5) = 1;
	G1(1, 7) = 1;
	G1(1, 8) = 1;
	G1(2, 6) = 1;
	G1(3, 7) = 1;
	G1(3, 8) = 1;
	G1(4, 8) = 1;
	G1(5, 7) = 1;
	G1(6, 9) = 1;
	G1(7, 9) = 1;
	G1(8, 9) = 1;
	ta1(0, 0) = 9;
	ta1(0, 1) = 7;
	ta1(0, 2) = 5;
	ta1(1, 0) = 8;
	ta1(1, 1) = 6;
	ta1(1, 2) = 5;
	ta1(2, 0) = 6;
	ta1(2, 1) = 5;
	ta1(2, 2) = 4;
	ta1(3, 0) = 7;
	ta1(3, 1) = 5;
	ta1(3, 2) = 3;
	ta1(4, 0) = 5;
	ta1(4, 1) = 4;
	ta1(4, 2) = 2;
	ta1(5, 0) = 7;
	ta1(5, 1) = 6;
	ta1(5, 2) = 4;
	ta1(6, 0) = 8;
	ta1(6, 1) = 5;
	ta1(6, 2) = 3;
	ta1(7, 0) = 6;
	ta1(7, 1) = 4;
	ta1(7, 2) = 2;
	ta1(8, 0) = 5;
	ta1(8, 1) = 3;
	ta1(8, 2) = 2;
	ta1(9, 0) = 7;
	ta1(9, 1) = 4;
	ta1(9, 2) = 2;
	primary(ini1, ta1, t1);
	prioritize(ini1, ta1, G1,t1);
	start = clock();
	initials(S1, ini1, ta1, G1, ts1, tc1, tr1);
	end = clock();
	cout << "Initial schedule:" << endl;
	prints(S1);
	rt = (double)(end - start) / (double)(CLOCKS_PER_SEC)*(double)(1000.000000);
	cout << " Now the total energy is: " << find_en(S1, p11, p12, p13, ps1)<<endl;
        cout << " Now the completion time is: " << find_tcom(S1) << endl;
	cout << "The running time of initial schedule of Graph 1 is "<<rt<<" ms"<< endl;
	start = clock();
	outerloop(S1, G1, ta1, ts1, tc1, tr1, p11, p12, p13, ps1, tmax1);
	end = clock();
	cout << "After Task Migration:" << endl;
	prints(S1);
	rt = (double)(end - start) / (double)(CLOCKS_PER_SEC)*(double)(1000.000000);
	cout << " Now the total energy is: " << find_en(S1, p11, p12, p13, ps1) << endl;
        cout << " Now the completion time is: " << find_tcom(S1) << endl;
	cout << "The running time of task migration of Graph 1 is "<<rt<<" ms"<< endl;
	cout << endl;
	int N2 = 8;  // the number of tasks
	int ts2 = 3, tc2 = 2, tr2 = 1;
	int t2;
	int tmax2=27 ;
	int p21=2;
	int p22=3;
	int p23=5;
	double ps2 = 1.5;
	t2 = ts2 + tc2 + tr2;
	Matrix<int, 2>G2(N2, N2);
	Matrix<int, 2>ta2(N2, K);
	vector<task>ini2(N2);
	vector<task>S2;
	G2(0, 2) = 1;
	G2(0, 3) = 1;
	G2(1, 3) = 1;
	G2(1, 4) = 1;
	G2(2, 5) = 1;
	G2(3, 5) = 1;
	G2(4, 5) = 1;
	G2(5, 6) = 1;
	G2(5, 7) = 1;
	ta2(0, 0) = 9;
	ta2(0, 1) = 8;
	ta2(0, 2) = 6;
	ta2(1, 0) = 8;
	ta2(1, 1) = 7;
	ta2(1, 2) = 6;
	ta2(2, 0) = 10;
	ta2(2, 1) = 9;
	ta2(2, 2) = 7;
	ta2(3, 0) = 6;
	ta2(3, 1) = 4;
	ta2(3, 2) = 2;
	ta2(4, 0) = 5;
	ta2(4, 1) = 3;
	ta2(4, 2) = 1;
	ta2(5, 0) = 7;
	ta2(5, 1) = 5;
	ta2(5, 2) = 4;
	ta2(6, 0) = 15;
	ta2(6, 1) = 11;
	ta2(6, 2) = 10;
	ta2(7, 0) = 6;
	ta2(7, 1) = 4;
	ta2(7, 2) = 1;
	primary(ini2, ta2, t2);
	prioritize(ini2, ta2, G2, t2);
	start = clock();
	initials(S2, ini2, ta2, G2, ts2, tc2, tr2);
	end = clock();
	cout << "Initial schedule:" << endl;
	prints(S2);
	rt = (double)(end - start) / (double)(CLOCKS_PER_SEC)*(double)(1000.000000);
	cout << " Now the total energy is: " << find_en(S2, p21, p22, p23, ps2) << endl;
	cout << " Now the completion time is: " << find_tcom(S2) << endl;
	cout << "The running time of initial schedule of Graph 2 is "<<rt<<" ms"<< endl;
	start = clock();
	outerloop(S2, G2, ta2, ts2, tc2, tr2, p21, p22, p23, ps2, tmax2);
	end = clock();
	cout << "After all of Task Migration:" << endl;
	prints(S2);
	rt = (double)(end - start) / (double)(CLOCKS_PER_SEC)*(double)(1000.000000);
	cout << " Now the total energy is: " << find_en(S2, p21, p22, p23, ps2) << endl;
	cout << " Now the completion time is: " << find_tcom(S2) << endl;
	cout << "The running time of initial schedule of Graph 2 is "<<rt<<" ms"<< endl;
	cout << endl;
}
