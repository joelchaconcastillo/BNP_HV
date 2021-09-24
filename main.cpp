#include <bits/stdc++.h>
#include "HypervolumeIndicator.h"
using namespace std;
/*
 * r   ={0, 10000, -1}
 *     ={1000, 0, -1}
   a[2]={2,2,1};
   a[0]={1,2,2};
   a[1]={2,1,2};
 * */
double hv3Dv2(vector<vector<double>> &points1, vector<double> &refPoint){
		if (points1.empty()) return 0.0;
//		std::sort(points.begin(), points.end(),
//					[] (vector<double> const& x, vector<double> const& y)
//					{ return (x[2] < y[2]); }
//				);
		map<double, pair<double, double > >points; //z,x,y,...
		for(auto i:points1)points[i[2]]=make_pair(i[0], i[1]);
		// add the first point
		std::map<double, double> front2D;
		vector<double>artifitialx={0, refPoint[1], -1};
		vector<double>artifitialy={refPoint[0], 0, -1};
		front2D[artifitialx[0]] = artifitialx[1];
		front2D[artifitialy[0]] = artifitialy[1];
//		front2D[points[0][0]]=points[0][1];
		//double prev_x2 = points[0][2];
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
		volume += area * (refPoint[2] - prev_x2);

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
void bnp3D(vector<vector<double> > &points, vector<double> &reference){

}
////////////////////////////////////////////////////////////7
int nobj=3, N;
vector<vector<double> >newpoints(int n, int m){
     srand(1);
     vector<vector<double> > points(n, vector<double>(m) ) ;
     for(int i = 0; i < n; i++){
   	  double x1= i/(double)n;
   	  double x2=rand()/(double)RAND_MAX;
   //	  points[i][0]=x1;
   //	  points[i][1]=1.0-pow(x1, 1.5);
//            points[i][0] = cos(M_PI*x1*0.5)*cos(M_PI*x2*0.5);
//            points[i][1] = cos(M_PI*x1*0.5)*sin(M_PI*x2*0.5);
//            points[i][2] = sin(M_PI*x1*0.5);
     points[i][0] = x1*x2;
     points[i][1] = x1*(1.0-x2);
     points[i][2] = (1.0-x1);
     }
     return points;
}
int main(){
   auto a=newpoints(10000, 3);
   vector<double> ref(3, 1.1);
 //  vector<vector<double> >a;
 //  a.push_back({1,0,0});
 //  a.push_back({0,10,0});
 //  a.push_back({0,0,10});
 //  a.push_back({10,0,1});
 //  a.push_back({0,0,10});
//   for(auto x:a){
//     for(auto xi:x)
//	cout << xi<<" ";
//     cout <<endl;
//   }
   cout << hv3D(a, ref)<<endl;
   cout << hv3Dv2(a, ref)<<endl;
  // cout << hv3Dv3(a, ref)<<endl;
}
