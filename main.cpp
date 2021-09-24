#include <bits/stdc++.h>

#include <chrono>
#include <iostream>
#include "HypervolumeIndicator.h"
using namespace std;
/*
 * r   ={0, 10000, -1}
 *     ={1000, 0, -1}
   a[2]={2,2,1};
   a[0]={1,2,2};
   a[1]={2,1,2};
 * */
double hv3Dv2(set<pair<double, pair<double, double> > > &points){
		if (points.empty()) return 0.0;
		// add the first point
		std::map<double, double> front2D;
		double inf = std::numeric_limits<double>::max();//inf can be more costly to work with!
		vector<double>artifitialx={-inf, 0, -inf};
		vector<double>artifitialy={0, -inf, -inf};
		front2D[artifitialx[0]] = artifitialx[1];
		front2D[artifitialy[0]] = artifitialy[1];
		double prev_x2 = points.begin()->first;
		double area = 0;//(refPoint[0] - points[0][0]) * (refPoint[1] - points[0][1]);
		double volume = 0.0;
		// process further points
		for( auto point:points){
			auto right= front2D.lower_bound(point.second.first);
			auto left=right;
			left--; //collineal points
			if (point.second.second >= left->second) continue;   // x is dominated note: strong dominance...
			// add chunk to volume
			volume += area * (point.first- prev_x2); //prev is the previous point that is not dominanted.
			while(right->second > point.second.second){
			       auto tmp=right;
                               right++;
			       area -= (right->first-tmp->first)*(left->second-tmp->second);
			       front2D.erase(tmp);
			}
			// add the new point
			area += (right->first- point.second.first) * (left->second - point.second.second);
			front2D[point.second.first] = point.second.second;
			prev_x2=point.first;
		}
		// add trailing chunk to volume
		volume += area * (-prev_x2);
		// return the result
		return volume;
}
double hv3D(vector<vector<double>> &points, vector<double> &refPoint){
		if (points.empty()) return 0.0;

		std::sort(points.begin(), points.end(),
					[] (vector<double> const& x, vector<double> const& y)
					{ return (x[2] < y[2]); }
				);
		// add the first point
		std::map<double, double> front2D;
		vector<double> const& x0 = points[0];
		front2D[x0[0]] = x0[1];
		double prev_x2 = x0[2];
		double area = (refPoint[0] - x0[0]) * (refPoint[1] - x0[1]);
		double volume = 0.0;

		// process further points
		for (size_t i=1; i< points.size(); i++)
		{
			vector<double> const& x = points[i];
			// check whether x is dominated and find "top" coordinate
			double t = refPoint[1];
			std::map<double, double>::iterator right = front2D.lower_bound(x[0]);
			std::map<double, double>::iterator left = right;
			if (right == front2D.end())
			{
				--left;
				t = left->second;
			}
			else
			{
				if (right->first == x[0])
				{
					t = left->second;
				}
				else if (left != front2D.begin())
				{
					--left;
					t = left->second;
				}
			}
			if (x[1] >= t) continue;   // x is dominated

			// add chunk to volume
			volume += area * (x[2] - prev_x2);

			// remove dominated points and corresponding areas
			while (right != front2D.end() && right->second >= x[1])
			{
				std::map<double, double>::iterator tmp = right;
				++right;
				const double r = (right == front2D.end()) ? refPoint[0] : right->first;
				area -= (r - tmp->first) * (t - tmp->second);
				front2D.erase(tmp);
			}

			// add the new point
			const double r = (right == front2D.end()) ? refPoint[0] : right->first;
			area += (r - x[0]) * (t - x[1]);
			front2D[x[0]] = x[1];

			// volume is processed up to here:
			prev_x2 = x[2];
		}

		// add trailing chunk to volume
		volume += area * (refPoint[2] - prev_x2);

		// return the result
		return volume;
}
vector<vector<double> > bnp3D(vector<vector<double> > inpoints, vector<double> &reference, int S){
   vector<vector<double> > outpoints;
   set< pair<double, pair<double, double > > >points;
   for(int t = 0; t < S; t++){
	pair<double, int> maxctr(0.0, -1);
	for(int i = 0; i < inpoints.size(); i++){
	   auto element=make_pair(inpoints[i][2]-reference[2], make_pair(inpoints[i][0]-reference[0], inpoints[i][1]-reference[1]));
	   auto pos=points.insert(element);
	   maxctr= max(maxctr, make_pair(hv3Dv2(points), i));
	   points.erase(pos.first);
        }
	   auto element=make_pair(inpoints[maxctr.second][2]-reference[2], make_pair(inpoints[maxctr.second][0]-reference[0], inpoints[maxctr.second][1]-reference[1]));
	   auto pos=points.insert(element);
       outpoints.push_back(inpoints[maxctr.second]); 
       iter_swap(inpoints.begin()+maxctr.second, inpoints.end()-1);
       inpoints.pop_back();
   }
   return outpoints;
}
vector<vector<double> > bnp2D(vector<vector<double> > inpoints, vector<double> &reference, int S){
   vector<vector<double> > outpoints;
   map<double, double> front2D;
   double inf = std::numeric_limits<double>::max();//inf can be more costly to work with!
   front2D[-inf]=0;
   front2D[0]=-inf;
   for(int t = 0; t < S; t++){
	pair<double, int> maxctr(-1.0, -1);
	for(int i = 0; i < inpoints.size(); i++){
	   auto point = inpoints[i];
	   for(int j = 0; j < point.size(); j++) point[j]-=reference[j];
	   auto right=front2D.lower_bound(point[0]);
	   auto left=right;
	   left--;
	   if(left->second < point[1])continue;
	   double area = (right->second- point[1]) * (left->first - point[0]);
	   cout << area<<endl;
	   maxctr= max(maxctr, make_pair(area, i));
          }
        front2D[inpoints[maxctr.second][0]-reference[0]] =inpoints[maxctr.second][1]-reference[1];
        outpoints.push_back(inpoints[maxctr.second]); 
        iter_swap(inpoints.begin()+maxctr.second, inpoints.end()-1);
        inpoints.pop_back();
   }
   return outpoints;
}
vector<vector<double> > lala(vector<vector<double> > &inpoints, vector<double> &reference, int S){
   vector<vector<double> > outpoints;
   set<pair<double, pair<double, double > > >points;
	for(int i = 0; i < inpoints.size(); i++){
	   auto element=make_pair(inpoints[i][2]-reference[2], make_pair(inpoints[i][0]-reference[0], inpoints[i][1]-reference[1]));
	   points.insert(element);
	}
	cout << hv3Dv2(points);
   return outpoints;
}
////////////////////////////////////////////////////////////7
vector<vector<double> >newpoints(int n, int m){
     srand(1);
     vector<vector<double> > points(n, vector<double>(m) ) ;
     for(int i = 0; i < n; i++){
   	  double x1= i/(double)n;
   	  double x2=rand()/(double)RAND_MAX;
	  if(m==2){
   	  points[i][0]=x1;
   	  points[i][1]=1.0-pow(x1, 0.2);
   	  }
	  else{
            points[i][0] = cos(M_PI*x1*0.5)*cos(M_PI*x2*0.5);
            points[i][1] = cos(M_PI*x1*0.5)*sin(M_PI*x2*0.5);
            points[i][2] = sin(M_PI*x1*0.5);
	  }
//     points[i][0] = x1*x2;
//     points[i][1] = x1*(1.0-x2);
//     points[i][2] = (1.0-x1);
     }
     return points;
}
ostream &operator<<(ostream &os, vector<vector<double> >&points){
   for(auto x:points){
     for(auto xi:x)
	os<< xi<<" ";
     os <<endl;
   }
   return os;
}
int main(){
using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
   auto a=newpoints(101, 2);
   vector<double> ref(2, 1.1);
//  cout << hv3D(a,ref);
  auto t1 = high_resolution_clock::now();
  for(int i = 0; i < 2500; i++){
   //auto b=bnp3D(a, ref, 100);
   auto b=bnp2D(a, ref, 100);
 //  cout << b<<endl;
  }
    auto t2 = high_resolution_clock::now();

    auto ms_int = duration_cast<milliseconds>(t2 - t1);
duration<double, std::milli> ms_double = t2 - t1;

  //  std::cout << ms_int.count() << "ms\n";
  //  std::cout << ms_double.count() << "ms";
//   cout <<b<<endl;
//   cout << hv3D(a, ref)<<endl;
//   cout << hv3Dv2(a, ref)<<endl;
}
