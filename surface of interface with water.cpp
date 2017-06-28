/*
created by Sedigheh Mahdavi for the paper
Mahdavi, Sedigheh, et al. "Computational structure analysis of biomacromolecule complexes by interface geometry." Computational biology and chemistry 47 (2013): 16-23.
Laboratory of Algorithms and Computational Geometry, Department of Mathematics and Computer Science,
Amirkabir University of Technology, Tehran, Iran
http://cg.aut.ac.ir/acg/ 
please contact the author through email:s_mahdavi@aut.ac.ir 
*/
///////*******************************************************************************************************************************
#include "stdafx.h"
#include <iostream>
#include <vector>
#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/range.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <fstream>
#include <CGAL/MP_Float.h>
#include <CGAL/long_double.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <cassert>
#include <CGAL/IO/Color.h>
#include <CGAL/PDB/Model.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Weighted_alpha_shape_euclidean_traits_3.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/double.h>
#include <CGAL/convex_hull_incremental_3.h>
#include <CGAL/Weighted_alpha_shape_euclidean_traits_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <list>
#include <algorithm>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/circulator.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gmpq.h>
#include <CGAL/enum.h>
#include <CGAL/intersections.h>
#include <boost/lexical_cast.hpp>
#include <cmath>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Two functors to compute the normals:  We assume the
// Simple_cartesian<double> Kernel here and use its global functions.
struct Facet_normal {
	template <class Facet>
	void operator()( Facet& f) {
		typename Facet::Halfedge_handle h = f.halfedge();
		typename Facet::Normal_3 normal = CGAL::cross_product(
			h->next()->vertex()->point() - h->vertex()->point(),
			h->next()->next()->vertex()->point()- h->next()->vertex()->point());
		f.normal() = normal / sqrt( normal * normal);
	}
};
struct Vertex_normal {
	template <class Vertex>
	void operator()( Vertex& v) {
		typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
		typedef typename Vertex::Halfedge_around_vertex_const_circulator Circ;
		Circ c = v.vertex_begin();
		Circ d = c;
		CGAL_For_all( c, d) {
			if ( ! c->is_border())
				normal = normal + c->facet()->normal();
		}
		v.normal() = normal / sqrt( normal * normal);
	}
};
// A redefined items class for the Polyhedron_3 with a refined vertex
// class that contains a member for the normal vector and a refined
// facet with a normal vector instead of the plane equation (this is
// an alternative solution instead of using Polyhedron_traits_with_normals_3).
template <class Refs, class T, class P, class Norm>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P> {
	Norm  norm;
public:
	My_vertex() {} // repeat mandatory constructors
	My_vertex( const P& pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt) {}
	typedef Norm Normal_3;
	Normal_3&       normal()       { return norm; }
	const Normal_3& normal() const { return norm; }
};
template <class Refs, class T, class Norm>
class My_facet : public CGAL::HalfedgeDS_face_base<Refs, T> {
	Norm  norm;
public:
	// no constructors to repeat, since only default constructor mandatory
	typedef Norm Normal_3;
	Normal_3&       normal()       { return norm; }
	const Normal_3& normal() const { return norm; }
};
struct My_items : public CGAL::Polyhedron_items_3 {
	template <class Refs, class Traits>
	struct Vertex_wrapper {
		typedef typename Traits::Point_3  Point;
		typedef typename Traits::Vector_3 Normal;
		typedef My_vertex<Refs, CGAL::Tag_true, Point, Normal> Vertex;
	};
	template <class Refs, class Traits>
	struct Face_wrapper {
		typedef typename Traits::Vector_3 Normal;
		typedef My_facet<Refs, CGAL::Tag_true, Normal> Face;
	};
};
using namespace std;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tie all types together and a small main function using it.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef  CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Weighted_alpha_shape_euclidean_traits_3<K> Gt;
typedef Gt::RT Weight;
typedef Gt::Bare_point Point;
typedef Gt::Weighted_point Weighted_point;
typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Rt;
typedef CGAL::Alpha_shape_3<Rt>        Alpha_shape_3;
typedef Rt::Vertex_handle Vertex_handle;
typedef Rt::Facet_iterator Facet_iterator2;
typedef Rt::Finite_facets_iterator Facet_iteratorrt;
typedef Rt::Finite_vertices_iterator  Finite_vertices_iterator;
typedef Rt::Finite_edges_iterator  Edge_iterator;
typedef Rt::Finite_cells_iterator  Finite_cells_iterator;
typedef Rt::Cell_circulator  	Cell_circulator; 
typedef Rt::All_cells_iterator   All_cells_iterator;
typedef Rt::Cell cell;
typedef Rt::Cell_handle     Cell_handle;
typedef K:: Point_3  	  Point_3 ;
typedef CGAL::Convex_hull_traits_3<K>             Traits1;
typedef K::Segment_3                              Segment_3;
typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;
typedef CGAL::Convex_hull_traits_3<K>   Convex_h_t;
typedef K::Triangle_3   	  Triangle_3; 
typedef  CGAL::Convex_hull_traits_3<K>::Plane_3      Plane_3  ; 
typedef Rt::Edge  	Edge; 
typedef CGAL::Convex_hull_d_traits_3<K>        Hull_traits_3;
typedef CGAL::Convex_hull_d< Hull_traits_3 >   Convex_hull_3;
typedef CGAL::Creator_uniform_3<double, Point_3>   Creator;
typedef CGAL::Polyhedron_3<K, My_items>           Polyhedron_3;
typedef Polyhedron_3::Vertex_iterator                    Vertex_iterator;
typedef Polyhedron_3::Facet::Plane_3  plan;
typedef Polyhedron_3::Facet     Facet ;
typedef Polyhedron_3::Halfedge  Halfedge;
typedef Polyhedron_3::Vertex  Vertex ;
typedef Polyhedron_3::Facet_iterator Facet_iterator1 ;
typedef Polyhedron_3::Facet_handle  Facet_handle;
typedef Polyhedron_3::Facet::Halfedge_around_facet_circulator CGAL_Halfedge_circulator;
//typedef for alphacomplex
typedef Alpha_shape_3::Cell_handle          Cell_handlealpha;
typedef Alpha_shape_3::Facet                Facetalpha;
typedef Alpha_shape_3::Edge                 Edgealpha;
typedef K::Vector_3                          Vector_3 ;
typedef Alpha_shape_3::Vertex_handle               Valpha;
//typedef for Delaunay Triangulation
typedef CGAL::Triangulation_3<K>      Triangulation;
typedef Triangulation::Edge  Edge1;
typedef Triangulation::Cell_handle    Cell_handleTr;
typedef Triangulation::Vertex_handle  Vertex_handleTr;
typedef Triangulation::Locate_type    Locate_typeTr;
typedef Triangulation::Point          PointTr;
typedef K::Tetrahedron_3  Tetrahedron_3;
typedef Triangulation::Finite_facets_iterator  Facet_iterator;
typedef Triangulation::Finite_cells_iterator  Finite_cells_iterator1;
typedef Rt::Facet_circulator   Facet_circulator;  
template <class Vector, class Range>
void insert(Vector &v, Range r) {
	v.insert(v.end(), r.begin(), r.end());
}
struct atoms
{
	char name;
	string id;
	double x,y,z,score;
	string res;
	char chain;
	double  r;
	string resno;
};
struct atomres
{
	Point_3 p;
	string res;
	string resno;
	char  chain;
};
struct residuenum
{
	string res1;
	string resno1;
	char  chain;
};
struct water
{
	double x,y,z;
	double r;
	char name;
	string resno;
	char chain;
};
struct interfaceatom 
{
	string res1;
	string resno1;
	char  chain;
	string nameatom;
	double x,y,z,score;
};
typedef vector<atoms> Vecatom;
typedef vector<Vecatom> Matatom;
typedef vector<cell> Veccell;
typedef vector<Veccell> Matcell;
typedef vector<residuenum>  resnum;
typedef Rt::Facet     facet ;
using namespace std;
/***************************************************************************main()*******************************************/
int main(int, char *[]) 
{
	char ch;
	std::vector<Weighted_point> P,Pia,Pib;
	std::vector<Cell_handle> cells;
	std::vector<Edge> iff,iff1,fedges ;//vector for edge interface after onionfilter
	int ci ;//count of infinite cell incident of edge
	int edgeinf=0;//Variable  for edge of convex hull by Regular_triangulation
	std::vector<K::Point_3> co1,co2;
	std::vector<K::Point_3> pointwater,iifchain;
	std::vector<K::Sphere_3> po1, po2;
	std::vector<K::Point_3> atom,atom1,atomi,atoma;
	K k;
	K::Sphere_3 sp;
	std::vector<K::Sphere_3> watersphere;
	std::vector<Weighted_point> waterweight,waterpoint;
	std::vector<Edge> edges,edgesf;//vairable is for Maintenance edge infinite iif
	ofstream pdb1 ("mol1.txt", ios::trunc);
	ofstream pdb2 ("mol2.txt", ios::trunc);
	ofstream pdb2c ("molc.txt", ios::trunc);
	ofstream pdbatom("interfaceatom.txt", ios::trunc);
	ofstream pdbatom1("interfaceatom1.txt", ios::trunc);
	ofstream pdbres("interfaceres.txt", ios::trunc);
	ofstream res("res.txt", ios::trunc);
	ofstream output("output.txt", ios::trunc); 
	ofstream numberres("numberres.xls", ios::trunc);
	std::ofstream pdbw ("molw.txt", std::ios::trunc);
	ofstream interfaceedge("interfaceedge.txt", ios::trunc); 
	ofstream allpoint("allpoint.txt", ios::trunc); 
	ofstream allsphere("allsphere.txt", ios::trunc); 
	vector<K::Point_3>  pointres;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//iterator for chains and get sphers of chains
	// double x,y,z;
	Matatom matrix;
	std::vector<atoms> protein;
	std::vector<char> chainv;
	std::vector<water> watermolecule;
	std::string a;
	std::vector<atomres> residue;
	std::cout<<"enter path file for example c:\\1si3.pdb "<<std::endl;
	std::cin>>a;
	std::cout<<std::endl;
	char * a1;
	int isize=a.length();
	a1=new char[isize+1];
	a.copy(a1,isize,0);
	a1[isize]='\0';
	std::ifstream pdb30(a1,std::ios::in);
	int cont_atom=0,cfile=0,u=0; 
	while (!pdb30.eof())
	{
		atoms  atomt;
		water  watert;
		cfile++;
		char line[256];
		pdb30.getline(line,256);
		string d=line;
		bool al=false,al1=false;
		string f="ATOM";
		string c,c1;
		int size=int(d.size());
		if (size>4)
		{
			for(int k=0;k<4;k++)
			{
				if (d[k]!=' ')
					c.insert(c.end(),d[k]);
			}
			if (f.compare(c)==0)
			{
				al=true;
				cont_atom++;
			}
		}
		string f1="HETATM";
		if (size>6)
		{
			for(int k=0;k<6;k++)
			{
				if (d[k]!=' ')
					c1.insert(c1.end(),d[k]);
			}
			if (f1.compare(c1)==0)
			{
				al1=true;
			}
		}
		//////////////////////
		/////print molecule
		int i1=0;
		int count_word1=0;
		bool b1=false;
		bool b2=false;
		while (i1<=(size-1))
		{
			if ((d[i1]==' ')&&(d[i1-1]!=' '))
			{
				string a;
				int m=i1-1;
				while (d[m]!=' ')
				{
					a.insert(a.end(),d[m]);
					if (m>0)
						m--;
					else 
						break;
				}//end while m;
				string c;
				for(int k=int(a.size())-1;k>=0;k--)
					c.insert(c.end(),a[k]);
				string f="COMPND";
				string f1="MOLECULE:";
				if (count_word1==0)
					if (f.compare(c)==0)
						b1=true;
				if (count_word1==2)
					if (f1.compare(c)==0)
						b2=true;
				if (b1&& b2)
				{
					//cout<<d<<endl;
					u=cfile;
					break;
				}
				count_word1++;
			}
			i1++;
		}

		/////////////////////
		if (al==true)
		{
			std::vector<string> words;
			int count_word=0;
			int i=0;
			while (i<=(size-1))
			{
				if ((d[i]==' ')&&(d[i-1]!=' '))
				{
					string a;
					int m=i-1;
					while (d[m]!=' ')
					{
						a.insert(a.end(),d[m]);
						if (m>0)
							m--;
						else 
							break;
					}//end while m;
					string c;
					for(int k=int(a.size())-1;k>=0;k--)
						c.insert(c.end(),a[k]);
					words.push_back(c);
					count_word++;
				}//end word
				i++;
			}//end while i line;
			///////////////condition for line 
			string t,t1,cy,cx;
			string ta,t1a;
			switch (count_word)
			{
			case 12:
				atomt.chain=(*(words.begin()+4))[0];
				atomt.x=boost::lexical_cast<double>(*(words.begin()+6));
				atomt.y=boost::lexical_cast<double>(*(words.begin()+7));
				atomt.z=boost::lexical_cast<double>(*(words.begin()+8));
				atomt.res=*(words.begin()+3);
				atomt.id=*(words.begin()+2);
				atomt.resno=*(words.begin()+5);
				atomt.name=(*(words.begin()+11))[0];
				break;
			case 11:
				cy=*(words.begin()+2);
				cx=*(words.begin()+4);
				if (int(cy.size())>3)
				{
					for(int k=0;k<3;k++)
						t.insert(t.end(),cy[k]);
					for(int k=3;k<=int(cy.size())-1;k++)
						t1.insert(t1.end(),cy[k]);
					atomt.chain=(*(words.begin()+3))[0];
					atomt.x=boost::lexical_cast<double>(*(words.begin()+5));
					atomt.y=boost::lexical_cast<double>(*(words.begin()+6));
					atomt.z=boost::lexical_cast<double>(*(words.begin()+7));
					atomt.res=t;
					atomt.resno=*(words.begin()+5);
					atomt.id=t1;
					atomt.name=(*(words.begin()+10))[0];
				}
				else 
					if (int(cx.size())>=4)
					{
						atomt.chain=(*(words.begin()+4))[0];
						atomt.x=boost::lexical_cast<double>(*(words.begin()+5));
						atomt.y=boost::lexical_cast<double>(*(words.begin()+6));
						atomt.z=boost::lexical_cast<double>(*(words.begin()+7));
						atomt.res=*(words.begin()+3);
						atomt.id=*(words.begin()+2);
						atomt.resno=*(words.begin()+5);
						atomt.name=(*(words.begin()+10))[0];
					}
					else
					{
						atomt.chain=(*(words.begin()+4))[0];
						atomt.x=boost::lexical_cast<double>(*(words.begin()+6));
						atomt.y=boost::lexical_cast<double>(*(words.begin()+7));
						atomt.z=boost::lexical_cast<double>(*(words.begin()+8));
						atomt.res=*(words.begin()+3);
						atomt.id=*(words.begin()+2);
						atomt.resno=*(words.begin()+5);
						atomt.name=(*(words.begin()+10))[0];
					}
					break;
			case 10:
				atomt.chain=(*(words.begin()+4))[0];
				atomt.x=boost::lexical_cast<double>(*(words.begin()+5));
				atomt.y=boost::lexical_cast<double>(*(words.begin()+6));
				atomt.z=boost::lexical_cast<double>(*(words.begin()+7));
				atomt.res=*(words.begin()+3);
				atomt.id=*(words.begin()+2);
				atomt.name=(*(words.begin()+9))[0];
				break;
			case 13:
				ta=*(words.begin()+2);
				t1a=*(words.begin()+3);
				int u=int(ta.size())-1;
				for(int k=0;k<u;k++)
					t1a.insert(t1a.end(),ta[k]);
				atomt.chain=(*(words.begin()+5))[0];
				atomt.x=boost::lexical_cast<double>(*(words.begin()+7));
				atomt.y=boost::lexical_cast<double>(*(words.begin()+8));
				atomt.z=boost::lexical_cast<double>(*(words.begin()+9));
				atomt.res=*(words.begin()+4);
				atomt.id=t1a;
				atomt.name=(*(words.begin()+12))[0];
				break;
			}
			//////////////////////////////////////
			////radius for atom/////////////////////////////////////////////////////////////////////////
			switch (atomt.name)
			{
			case 'C':
				atomt.r=2.89;
				break;
			case 'N':
				atomt.r=2.4025  ;
				break;
			case 'O':
				atomt.r= 2.3104;
				break;
			case 'S':
				atomt.r=3.4225;
				break;
			case 'P':
				atomt.r=3.24;
				break;
			}
			if (!(atomt.name=='H'))
				protein.push_back(atomt);
		}//end if al for atoms
		/////////////////////////////////////////////////////////////
		else
			if(al1==true) 
			{
				std::vector<string> words;
				int count_word=0;
				int i=0;
				while (i<=(size-1))
				{
					if ((d[i]==' ')&&(d[i-1]!=' '))
					{
						string a;
						int m=i-1;
						while (d[m]!=' ')
						{
							a.insert(a.end(),d[m]);
							if (m>0)
								m--;
							else 
								break;
						}//end while m;
						string c;
						for(int k=int(a.size())-1;k>=0;k--)
							c.insert(c.end(),a[k]);
						words.push_back(c);
						count_word++;
					}//end word
					i++;
				}//end while i line;
				///////////////condition for line 
				switch (count_word)
				{
				case 12:
					watert.chain=(*(words.begin()+4))[0];
					watert.x=boost::lexical_cast<double>(*(words.begin()+6));
					watert.y=boost::lexical_cast<double>(*(words.begin()+7));
					watert.z=boost::lexical_cast<double>(*(words.begin()+8));
					watert.name=(*(words.begin()+11))[0] ;
					watert.resno=(*(words.begin()+3));
					break;
				case 11:
					watert.chain=(*(words.begin()+3))[0];
					watert.x=boost::lexical_cast<double>(*(words.begin()+5));
					watert.y=boost::lexical_cast<double>(*(words.begin()+6));
					watert.z=boost::lexical_cast<double>(*(words.begin()+7));
					watert.name=(*(words.begin()+10))[0] ;
					watert.resno=(*(words.begin()+3));
					break;
				}// end of switch
				////radius for atom/////////////////////////////////////////////////////////////////////////
				watert.r=1.96;
				string b="HOH";
				if (watert.name=='O')
					watermolecule.push_back(watert);
			}
			//////////end if al1 for hetatm
	}//end while file;
	std::vector<water>::iterator iw=watermolecule.begin();
	while( iw!=watermolecule.end())
	{
		Point_3 t1= Point_3(iw->x,iw->y,iw->z);
		K::Sphere_3 t=K::Sphere_3(t1,iw->r);
		Weighted_point tw=Weighted_point(t.center(),t.squared_radius());
		waterweight.push_back(tw);
		pdbw<<t<<std::endl;
		iw++;
	}
	///////////////////////////////////////////////////////////////////
	std::cout<<"for RNA-protein and DNA-protein complexes ,please  enter number file 2 for chains in proteins: "<<std::endl;
	std::ifstream pdbw1("molw.txt");
	while (!pdbw1.eof())
	{
		K::Sphere_3  t;
		pdbw1 >> t;
		Weighted_point tw=Weighted_point(t.center(),t.squared_radius());
		watersphere.push_back(t);
		waterpoint.push_back(tw);
	}
	//////////////////////////////
	std::vector<atoms>::iterator ip=protein.begin();
	while (ip!=protein.end())
	{
		std::vector<char>::iterator ic=find(chainv.begin(),chainv.end(),ip->chain);
		if (ic==chainv.end())
			chainv.push_back(ip->chain);
		ip++;
	}
	cout<<cont_atom<<endl;
	cout<<chainv.size()<<endl;
	std::vector<char>::iterator ic=chainv.begin();
	while (ic!=chainv.end())
	{
		std::vector<atoms> at;
		std::vector<atoms>::iterator ip=protein.begin();
		while (ip!=protein.end())
		{
			if (ip->chain==*ic)
				at.push_back(*ip);
			ip++;
		}//end while
		matrix.push_back(at);
		ic++;
	}//end while chainv
	//write to mol1 and mol2
	std::vector<Vecatom>::iterator im=matrix.begin();
	while (im!=matrix.end())
	{
		std::vector<atoms> at1=*im;
		std::vector<atoms>::iterator ip=at1.begin();
		std::cout<<"please enter number file for this chain:     "<<ip->chain<<std::endl;
		std::cout<<at1.size()<<std::endl;
		output<<"please enter number file for this chain:     "<<ip->chain<<std::endl;
		output<<at1.size()<<std::endl;
		int number;
		int number1;
		std::cin>>number; 
		if (number==1)
		{
			while(ip!=at1.end())
			{
				Point_3 t1= Point_3(ip->x,ip->y,ip->z);
				K::Sphere_3 t=K::Sphere_3(t1,ip->r);
				Weighted_point tw=Weighted_point(t.center(),t.squared_radius());
				atomres rest;
				rest.chain=ip->chain;
				rest.res=ip->res;
				rest.resno=ip->resno;
				residue.push_back(rest);
				pdb1<<t<<std::endl;
				ip++;
			}
		}
		else if (number==2)
		{
			while(ip!=at1.end())
			{
				Point_3 t1= Point_3(ip->x,ip->y,ip->z);
				K::Sphere_3 t=K::Sphere_3(t1,ip->r);
				Weighted_point tw=Weighted_point(t.center(),t.squared_radius());
				atomres rest;
				rest.chain=ip->chain;
				rest.res=ip->res;
				rest.resno=ip->resno;
				residue.push_back(rest);
				pdb2<<t<<std::endl;
				ip++;
			}
		}
		im++;
	}
	//insert into iifchain
	std::ofstream pdbpvoid("pdbpvoid.txt", std::ios::trunc );
	//put atoms any molecoul in its vector   
	std::vector<K::Sphere_3> points1;
	std::ifstream pdb3("mol1.txt");
	K::Sphere_3   t;
	while (!pdb3.eof())
	{
		pdb3>> t;
		po1.push_back(t);
		co1.push_back(t.center());
		points1.push_back(t);
		allpoint<<t.center()<<std::endl;
		allsphere<<t<<std::endl;
		atomres r;
		pointres.push_back(t.center());

	}
	std::ifstream pdb4("mol2.txt");
	K::Sphere_3   t1;
	while (!pdb4.eof())
	{
		pdb4>> t1;
		po2.push_back(t1);
		co2.push_back(t1.center());
		points1.push_back(t1);
		allpoint<<t1.center()<<std::endl;
		allsphere<<t1<<std::endl;
		atomres r;
		pointres.push_back(t.center());

	}
	////////////////////////////////////
	vector<atomres>::iterator irrd=residue.begin();
	vector<K::Point_3>::iterator irrp=pointres.begin();
	while( irrd!=residue.end())
	{ 
		irrd->p=*irrp;
		irrp++;
		irrd++;
	}
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*******************************************************TR and find bicolor edge**********************************************************************/
	//transformation sphere points into Weighted points	  
	std::vector<K::Sphere_3>::iterator i=points1.begin();
	int numberofatom=0;
	while (i!=points1.end())
	{
		P.push_back(Weighted_point(i->center(),i->squared_radius()));
		i++;
		numberofatom++;
	}
	//////////// insert water on to P
	std::vector<Weighted_point>::iterator iw1=waterpoint.begin();
	while (iw1!=waterpoint.end())
	{
		P.push_back(*iw1);
		pointwater.push_back(iw1->point());
		K::Sphere_3 t=K::Sphere_3(iw1->point(),iw1->weight() );
		points1.push_back(t);
		//allpoint<<iw1->point()<<std::endl;
		allsphere<<t<<std::endl;
		iw1++;
	}

	//insert points into regular
	Rt T;
	T.insert (P.begin(), P.end());
	assert( T.is_valid() );
	assert( T.dimension() == 3 );
	Edge_iterator ei = T.finite_edges_begin();
	//travels of edge in rt
	std::ofstream pdb5("pdb1.txt", std::ios::trunc );//file output maintance facet in pw3d
	int kb=0; 
	for( ; ei != T.finite_edges_end(); ++ei) 
	{ 
		Segment_3 s;
		s=T.segment(*ei);
		//condition for interface edge
		if ((find(co1.begin(), co1.end(), s.source()) != co1.end()) && (find(co2.begin(), co2.end(), s.target()) != co2.end())||(find (co2.begin(), co2.end(), s.source()) != co2.end()) && (find (co1.begin(), co1.end(), s.target()) != co1.end()))
		{
			//obtain face dual of edge
			if (find(atomi.begin(),atomi.end(), s.source()) == atomi.end())

				atomi.push_back( s.source());

			if (find(atomi.begin(),atomi.end(), s.target()) == atomi.end())

				atomi.push_back( s.target());

			if (find(edgesf.begin(),edgesf.end(),*ei) == edgesf.end())
				edgesf.push_back(*ei);
		}
		////////inseert water interface atoms //////////////////////////////////////////////////////////////////////////////////
		if((find(pointwater.begin(), pointwater.end(),s.source())!= pointwater.end())&&(find(co1.begin(),co1.end(), s.target())!= co1.end()))
		{
			std::vector< Point_3>::iterator  i1=find(pointwater.begin(), pointwater.end(),s.source());
			bool t=false;
			Point_3 k;
			Edge_iterator ei1 = T.finite_edges_begin();
			for( ; ei1 != T.finite_edges_end(); ++ei1) 
			{ 
				Segment_3 s1;
				s1=T.segment(*ei1);
				if (s.source()==s1.source())
				{
					if (find(co2.begin(),co2.end(), s1.target())!= co2.end())
						t=true;
					k=s1.source();
				}
				if (s.source()==s1.target())
				{
					if (find(co2.begin(),co2.end(), s1.source())!= co2.end())
						t=true;
					k=s1.target();
				}
			}
			/////end for
			if (t) 
			{
				if (find(edgesf.begin(),edgesf.end(),*ei) == edgesf.end())
					edgesf.push_back(*ei);
				if (find(atomi.begin(),atomi.end(), s.source()) == atomi.end())
					atomi.push_back( s.source());
				if (find(atomi.begin(),atomi.end(), s.target()) == atomi.end())
					atomi.push_back( s.target());

				kb++;
			}
			//////////////////////////////end if 
		}///////////////////////////////////////////////////////////////////end if 
		if((find(pointwater.begin(), pointwater.end(),s.target())!= pointwater.end())&&(find(co1.begin(),co1.end(), s.source())!= co1.end()))
		{
			std::vector< Point_3>::iterator  i1=find(pointwater.begin(), pointwater.end(),s.target());
			bool t=false;
			Point_3 k;
			Edge_iterator ei1 = T.finite_edges_begin();
			for( ; ei1 != T.finite_edges_end(); ++ei1) 
			{ 
				Segment_3 s1;
				s1=T.segment(*ei1);
				if (s.target()==s1.source())
				{
					if (find(co2.begin(),co2.end(), s1.target())!= co2.end())
						t=true;
					k=s1.target();
				}

				if (s.target()==s1.target())
				{
					if (find(co2.begin(),co2.end(), s1.source())!= co2.end())
						t=true;
					k=s1.target();
				}
			}/////end for
			if (t) 
			{
				if (find(edgesf.begin(),edgesf.end(),*ei) == edgesf.end())
					edgesf.push_back(*ei);
				if (find(atomi.begin(),atomi.end(), s.source()) == atomi.end())
					atomi.push_back( s.source());
				if (find(atomi.begin(),atomi.end(), s.target()) == atomi.end())
					atomi.push_back( s.target());

				kb++;
			}
			//////////////////////////////end if 
		}///////////////////////////////////////////////////////////////////end if 
		////////insert water interface atoms //////////////////////////////////////////////////////////////////////////////////
		if((find(pointwater.begin(), pointwater.end(),s.source())!= pointwater.end())&&(find(co2.begin(),co2.end(), s.target())!= co2.end()))
		{
			std::vector< Point_3>::iterator  i1=find(pointwater.begin(), pointwater.end(),s.source());
			bool t=false;
			Point_3 k;
			Edge_iterator ei1 = T.finite_edges_begin();
			for( ; ei1 != T.finite_edges_end(); ++ei1) 
			{ 
				Segment_3 s1;
				s1=T.segment(*ei1);
				if (s.source()==s1.source())
				{
					if (find(co1.begin(),co1.end(), s1.target())!= co1.end())
						t=true;
					k=s1.source();
				}
				if (s.source()==s1.target())
				{
					if (find(co1.begin(),co1.end(), s1.source())!= co1.end())
						t=true;
					k=s1.target();
				}
			}
			/////end for
			if (t) 
			{
				if (find(edgesf.begin(),edgesf.end(),*ei) == edgesf.end())
					edgesf.push_back(*ei);
				if (find(atomi.begin(),atomi.end(), s.source()) == atomi.end())
					atomi.push_back( s.source());
				if (find(atomi.begin(),atomi.end(), s.target()) == atomi.end())
					atomi.push_back( s.target());

				kb++;
			}
			//////////////////////////////end if 
		}///////////////////////////////////////////////////////////////////end if 
		if((find(pointwater.begin(), pointwater.end(),s.target())!= pointwater.end())&&(find(co2.begin(),co2.end(), s.source())!= co2.end()))
		{
			std::vector< Point_3>::iterator  i1=find(pointwater.begin(), pointwater.end(),s.target());
			bool t=false;
			Point_3 k;
			Edge_iterator ei1 = T.finite_edges_begin();
			for( ; ei1 != T.finite_edges_end(); ++ei1) 
			{ 
				Segment_3 s1;
				s1=T.segment(*ei1);
				if (s.target()==s1.source())
				{
					if (find(co1.begin(),co1.end(), s1.target())!= co1.end())
						t=true;
					k=s1.target();
				}
				if (s.target()==s1.target())
				{
					if (find(co1.begin(),co1.end(), s1.source())!= co1.end())
						t=true;
					k=s1.target();
				}
			}
			/////end for
			if (t) 
			{
				if (find(edgesf.begin(),edgesf.end(),*ei) == edgesf.end())
					edgesf.push_back(*ei);
				if (find(atomi.begin(),atomi.end(), s.source()) == atomi.end())
					atomi.push_back( s.source());
				if (find(atomi.begin(),atomi.end(), s.target()) == atomi.end())
					atomi.push_back( s.target());

				kb++;
			}
			//////////////////////////////end if 
		}///////////////////////////////////////////////////////////////////end if*/ 
	}/////for 
		// atomi ,atoms of interface infinite and edgesf,edges of interface infinite
	////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////surface for molecule A
	std::vector<Point_3> atomia;
	std::vector<K::Sphere_3>::iterator isa1=po1.begin();
	while (isa1!=po1.end())
	{
		Pia.push_back(Weighted_point(isa1->center(),isa1->squared_radius()+2.8*sqrt(isa1->squared_radius())+1.96));
		isa1++;
	}
	//////////////////////////////////////////////////////////
	Alpha_shape_3  asa(Pia.begin(), Pia.end(),0, Alpha_shape_3::GENERAL);
	std::vector< Valpha>       verteiesa;
	asa.get_alpha_shape_vertices(std::back_inserter(verteiesa),
		Alpha_shape_3::SINGULAR);
	asa.get_alpha_shape_vertices(std::back_inserter(verteiesa),
		Alpha_shape_3::REGULAR);
	std::vector<K::Point_3> pointsura;
	std::vector< Valpha>::iterator  ipa4a=verteiesa.begin();
	int sa=0;  
	while (ipa4a!= verteiesa.end())
	{
		K::Point_3   k=(*ipa4a)->point();
		pointsura.push_back(k);
		ipa4a++;
		sa++;
	}
	////////////////////////////////////////surface for molecule B
	std::vector<Point_3> atomib;
	std::vector<K::Sphere_3>::iterator isb1=po2.begin();
	while (isb1!=po2.end())
	{
		Pib.push_back(Weighted_point(isb1->center(),isb1->squared_radius()+2.8*sqrt(isb1->squared_radius())+1.96));
		isb1++;
	}
	//////////////////////////////////////////////////////////
	Alpha_shape_3  asb(Pib.begin(), Pib.end(),0, Alpha_shape_3::GENERAL);
	std::vector< Valpha>       verteiesb;
	asb.get_alpha_shape_vertices(std::back_inserter(verteiesb),
		Alpha_shape_3::SINGULAR);
	asb.get_alpha_shape_vertices(std::back_inserter(verteiesb),
		Alpha_shape_3::REGULAR);
	std::vector<K::Point_3> pointsurb;
	std::vector< Valpha>::iterator  ipa4b=verteiesb.begin();
	int sb=0; 
	while (ipa4b!= verteiesb.end())
	{
		K::Point_3   k=(*ipa4b)->point();
		pointsurb.push_back(k);
		ipa4b++;
		sb++;
	}
	//////////////////////////////////  
	std::vector<K::Point_3>::iterator ik=atomi.begin();
	while (ik!=atomi.end())
	{

		if (find(pointwater.begin(),pointwater.end(),*ik)==pointwater.end())
		{
			if (((find(pointsura.begin(),pointsura.end(),*ik)==pointsura.end())&&(find(pointsurb.begin(),pointsurb.end(),*ik)==pointsurb.end())))
				atomi.erase(ik);
			else 
				ik++;
		}
		else 
			ik++;
	}
	////////////////////////////////////////////////////////////////////////////////////////
	std::vector<Edge>::iterator ealpha12=edgesf.begin();
	while ( ealpha12!=edgesf.end())
	{
		K::Segment_3 s1;
		s1=T.segment(*ealpha12);
		Point_3  p = s1.source();    
		Point_3  q =s1.target(); 
		if (((find(atomi.begin(),atomi.end(),p)==atomi.end())||(find(atomi.begin(),atomi.end(),q)==atomi.end())))
			edgesf.erase(ealpha12);
		else
			ealpha12++;
	}
	///////////////////////////////remove water moleculs none interface and other molecules
	std::vector<Edge>::iterator ealpha12w=edgesf.begin();
	std::vector<K::Point_3> pointwaterr;
	while (ealpha12w!=edgesf.end())
	{
		K::Segment_3 s1;
		s1=T.segment(*ealpha12w);
		Point_3  p = s1.source();    
		Point_3  q =s1.target();
		if(find(pointwater.begin(), pointwater.end(),p)!= pointwater.end()) 
		{
			if(find(co1.begin(), co1.end(),q)!= co1.end()) 
			{
				bool k=true;
				std::vector<Edge>::iterator ealpha12w1=edgesf.begin();
				while ( ealpha12w1!=edgesf.end())
				{
					K::Segment_3 s2;
					s2=T.segment(*ealpha12w1);
					Point_3  p = s2.source();    
					Point_3  q =s2.target();
					if (((p==s1.source())&&(find(co2.begin(), co2.end(),q)!= co2.end()))||((q==s1.source())&&(find(co2.begin(), co2.end(),p)!= co2.end())))
					{
						k=false;
						break;
					}
					ealpha12w1++;
				}//end while 
				if (k)
					if(find(pointwaterr.begin(), pointwaterr.end(),p)== pointwaterr.end()) 
						pointwaterr.push_back(p);
			}//end if 
			else
			{
				bool k=true;
				std::vector<Edge>::iterator ealpha12w1=edgesf.begin();
				while ( ealpha12w1!=edgesf.end())
				{
					K::Segment_3 s2;
					s2=T.segment(*ealpha12w1);
					Point_3  p = s2.source();    
					Point_3  q =s2.target();
					if (((p==s1.source())&&(find(co1.begin(), co1.end(),q)!= co1.end()))||((q==s1.source())&&(find(co1.begin(), co1.end(),p)!= co1.end())))
					{
						k=false;
						break;
					}
					ealpha12w1++;
				}//end while
				if (k)
					if(find(pointwaterr.begin(), pointwaterr.end(),p)== pointwaterr.end()) 
						pointwaterr.push_back(p);
			}
		}//endif 
		//////////////////////////////////
		if(find(pointwater.begin(), pointwater.end(),q)!= pointwater.end()) 
		{
			if(find(co1.begin(), co1.end(),p)!= co1.end()) 
			{
				bool k=true;
				std::vector<Edge>::iterator ealpha12w1=edgesf.begin();
				while ( ealpha12w1!=edgesf.end())
				{
					K::Segment_3 s2;
					s2=T.segment(*ealpha12w1);
					Point_3  p = s2.source();    
					Point_3  q =s2.target();
					if (((p==s1.target())&&(find(co2.begin(), co2.end(),q)!= co2.end()))||((q==s1.target())&&(find(co2.begin(), co2.end(),p)!= co2.end())))
					{
						k=false;
						break;
					}
					ealpha12w1++;
				}//end while 
				if (k)
					if(find(pointwaterr.begin(), pointwaterr.end(),q)== pointwaterr.end()) 
						pointwaterr.push_back(q);
			}//end if 
			else
			{
				bool k=true;
				std::vector<Edge>::iterator ealpha12w1=edgesf.begin();
				while ( ealpha12w1!=edgesf.end())
				{
					K::Segment_3 s2;
					s2=T.segment(*ealpha12w1);
					Point_3  p = s2.source();    
					Point_3  q =s2.target();
					if (((p==s1.target())&&(find(co1.begin(), co1.end(),q)!= co1.end()))||((q==s1.target())&&(find(co1.begin(), co1.end(),p)!= co1.end())))
					{
						k=false;
						break;
					}
					ealpha12w1++;
				}//end while
				if (k)
					if(find(pointwaterr.begin(), pointwaterr.end(),q)== pointwaterr.end()) 
						pointwaterr.push_back(q);
			}

		}//end if
		ealpha12w++;
	}
	/////////////////
	///////////////////////////////////////remove edges  /////////////////////////////
	std::vector<Edge>::iterator ealpha122=edgesf.begin();
	while (ealpha122!=edgesf.end())
	{
		K::Segment_3 s1;
		s1=T.segment(*ealpha122);
		Point_3  p = s1.source();    
		Point_3  q =s1.target();
		if ((find(pointwaterr.begin(), pointwaterr.end(),p)!= pointwaterr.end())||(find(pointwaterr.begin(), pointwaterr.end(),q)!= pointwaterr.end()))
			edgesf.erase(ealpha122);
		else
			ealpha122++;
	}
	///////////////////
	std::vector<K::Point_3> atomi1;
	std::vector<Edge>::iterator ealpha122a=edgesf.begin();
	while (ealpha122a!=edgesf.end())
	{
		K::Segment_3 s1;
		s1=T.segment(*ealpha122a);
		Point_3  p = s1.source();    
		Point_3  q =s1.target();
		if (find(atomi1.begin(), atomi1.end(),p)== atomi1.end())
			atomi1.push_back(p);
		if (find(atomi1.begin(), atomi1.end(),q)== atomi1.end())
			atomi1.push_back(q);
		ealpha122a++;
	}
	////////////
	///////////////////////////////
	std::vector<K::Point_3> pointwateri;
	/////////////////////////////////////////////////////////////////
	std::vector<Edge>::iterator ealpha1= edgesf.begin();
	int coutw=0,coutw1=0;
	while (ealpha1!= edgesf.end())
	{
		K::Segment_3 s1;
		s1=T.segment(*ealpha1);
		double h;
		double h1;
		double h2,r1,r2;
		Point_3  p = s1.source();   
		Point_3  q =s1.target(); 
		Point_3 p1=Point_3(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
		Point_3 q1=Point_3(CGAL::to_double(q.x()),CGAL::to_double(q.y()),CGAL::to_double(q.z()));
		h=sqrt(CGAL::to_double((p.x()-q.x())*(p.x()-q.x())+(p.y()-q.y())*(p.y()-q.y())+(p.z()-q.z())*(p.z()-q.z())));
		std::vector<K::Sphere_3>::iterator i1=points1.begin();
		h1=2.8;
		while (i1!=points1.end())
		{
			if (i1->center()==p) 
			{ 
				h1=sqrt(CGAL::to_double(i1->squared_radius()))+h1;
			}
			if (i1->center()==q)
			{
				h1=sqrt(CGAL::to_double(i1->squared_radius())) +h1;
			}
			i1++;
		}
		///////////condition for chains in interfce
		if ((h<=h1 ))
		{
			if ((find(atom1.begin(), atom1.end(), p) == atom1.end()))
			{
				atom1.push_back(p);
			}
			if ((find(atom1.begin(), atom1.end(), q) == atom1.end()))
			{
				atom1.push_back(q);
			}
			if (find(iff1.begin(),iff1.end(), *ealpha1) == iff1.end())
				iff1.push_back(*ealpha1);
			std::vector<Point_3>::iterator i3=find(pointwater.begin(), pointwater.end(),s1.target());
			if (i3!=pointwater.end())
			{
				if ((find(pointwateri.begin(), pointwateri.end(), s1.target()) == pointwateri.end()))
					pointwateri.push_back(s1.target());
			}
			std::vector<Point_3>::iterator i31=find(pointwater.begin(), pointwater.end(),s1.source());
			if (i31!=pointwater.end())
			{
				if ((find(pointwateri.begin(), pointwateri.end(), s1.source()) == pointwateri.end()))
					pointwateri.push_back(s1.source());
			}
		}
		ealpha1++;
	}
	
	///////////////////////////////remove water moleculs none interface and other molecules
	std::vector<Edge>::iterator ealpha12w1=iff1.begin();
	std::vector<K::Point_3> pointwaterr1;
	while (ealpha12w1!=iff1.end())
	{
		K::Segment_3 s1;
		s1=T.segment(*ealpha12w1);
		Point_3  p = s1.source();    
		Point_3  q =s1.target();
		if(find(pointwater.begin(), pointwater.end(),p)!= pointwater.end()) 
		{
			if(find(co1.begin(), co1.end(),q)!= co1.end()) 
			{
				bool k=true;
				std::vector<Edge>::iterator ealpha12w12=iff1.begin();
				while ( ealpha12w12!=iff1.end())
				{
					K::Segment_3 s2;
					s2=T.segment(*ealpha12w12);
					Point_3  p = s2.source();    
					Point_3  q =s2.target();
					if (((p==s1.source())&&(find(co2.begin(), co2.end(),q)!= co2.end()))||((q==s1.source())&&(find(co2.begin(), co2.end(),p)!= co2.end())))
					{
						k=false;
						break;
					}
					ealpha12w12++;
				}//end while 
				if (k)
					if(find(pointwaterr1.begin(), pointwaterr1.end(),p)== pointwaterr1.end()) 
						pointwaterr1.push_back(p);
			}//end if 
			else
			{
				bool k=true;
				std::vector<Edge>::iterator ealpha12w12=iff1.begin();
				while ( ealpha12w12!=iff1.end())
				{
					K::Segment_3 s2;
					s2=T.segment(*ealpha12w12);
					Point_3  p = s2.source();    
					Point_3  q =s2.target();
					if (((p==s1.source())&&(find(co1.begin(), co1.end(),q)!= co1.end()))||((q==s1.source())&&(find(co1.begin(), co1.end(),p)!= co1.end())))
					{
						k=false;
						break;
					}
					ealpha12w12++;
				}//end while
				if (k)
					if(find(pointwaterr1.begin(), pointwaterr1.end(),p)== pointwaterr1.end()) 
						pointwaterr1.push_back(p);
			}
		}//endif 
		//////////////////////////////////
		if(find(pointwater.begin(), pointwater.end(),q)!= pointwater.end()) 
		{
			if(find(co1.begin(), co1.end(),p)!= co1.end()) 
			{
				bool k=true;
				std::vector<Edge>::iterator ealpha12w12=iff1.begin();
				while ( ealpha12w12!=iff1.end())
				{
					K::Segment_3 s2;
					s2=T.segment(*ealpha12w12);
					Point_3  p = s2.source();    
					Point_3  q =s2.target();
					if (((p==s1.target())&&(find(co2.begin(), co2.end(),q)!= co2.end()))||((q==s1.target())&&(find(co2.begin(), co2.end(),p)!= co2.end())))
					{
						k=false;
						break;
					}
					ealpha12w12++;
				}//end while 
				if (k)
					if(find(pointwaterr1.begin(), pointwaterr1.end(),q)== pointwaterr1.end()) 
						pointwaterr1.push_back(q);
			}//end if 
			else
			{
				bool k=true;
				std::vector<Edge>::iterator ealpha12w12=iff1.begin();
				while ( ealpha12w12!=iff1.end())
				{
					K::Segment_3 s2;
					s2=T.segment(*ealpha12w12);
					Point_3  p = s2.source();    
					Point_3  q =s2.target();
					if (((p==s1.target())&&(find(co1.begin(), co1.end(),q)!= co1.end()))||((q==s1.target())&&(find(co1.begin(), co1.end(),p)!= co1.end())))
					{
						k=false;
						break;
					}
					ealpha12w12++;
				}//end while
				if (k)
					if(find(pointwaterr1.begin(), pointwaterr1.end(),q)== pointwaterr1.end()) 
						pointwaterr1.push_back(q);
			}

		}//end if
		ealpha12w1++;
	}
	/////////////////
	///////////////////////////////////////remove edges  /////////////////////////////
	std::vector<Edge>::iterator ealpha1223=iff1.begin();
	while (ealpha1223!=iff1.end())
	{
		K::Segment_3 s1;
		s1=T.segment(*ealpha1223);
		Point_3  p = s1.source();    
		Point_3  q =s1.target();
		if ((find(pointwaterr1.begin(), pointwaterr1.end(),p)!= pointwaterr1.end())||(find(pointwaterr1.begin(), pointwaterr1.end(),q)!= pointwaterr1.end()))
			iff1.erase(ealpha1223);
		else
			ealpha1223++;
	}
	///////////////////
	std::vector<K::Point_3> atomi12;
	std::vector<K::Point_3> wateri;
	std::vector<Edge>::iterator ealpha122a1=iff1.begin();
	while (ealpha122a1!=iff1.end())
	{
		K::Segment_3 s1;
		s1=T.segment(*ealpha122a1);
		interfaceedge<<s1<<std::endl;
		Point_3  p = s1.source();    
		Point_3  q =s1.target();
		if (find(atomi12.begin(), atomi12.end(),p)== atomi12.end())
			atomi12.push_back(p);
		if (find(atomi12.begin(), atomi12.end(),q)== atomi12.end())
			atomi12.push_back(q);
		if(find(pointwater.begin(), pointwater.end(),q)!= pointwater.end()) 
			if(find(wateri.begin(), wateri.end(),q)== wateri.end()) 
				wateri.push_back(q);
		if(find(pointwater.begin(), pointwater.end(),p)!= pointwater.end()) 
			if(find(wateri.begin(), wateri.end(),p)== wateri.end()) 
				wateri.push_back(p);
		ealpha122a1++;
	}
	////////////
	std::cout<<"number of interface atom  iif infinite  final "<<atomi12.size()<<std::endl;
	output<<"number of interface atom  iif infinite  final  "<<atomi12.size()<<std::endl;
	std::cout<<"number of interface facet iif  infinite final  "<<iff1.size()<<std::endl;
	output<<"number of interface facet iif  infinite final  "<<iff1.size()<<std::endl;
	std::cout<<"number of interface water  final  "<<wateri.size()<<std::endl;
	output<<"number of interface water  final  "<<wateri.size()<<std::endl;
	////////////////////////////////////////////////////////////////write interface atoms in a file
	ofstream pdbname("nameatomif.txt", ios::trunc);
	vector<Point_3>::iterator ir11= atomi12.begin();
	while (ir11!=atomi12.end())
	{
		vector<atoms>::iterator ip=protein.begin();
		while (ip!=protein.end())
		{
			atoms t=*ip;
			double x2,y2,z2;
			x2=t.x;
			y2=t.y;
			z2=t.z;
			K::Point_3 t1=*ir11;
			double x1,y1,z1;
			x1=CGAL::to_double(t1.x());
			y1=CGAL::to_double(t1.y());
			z1=CGAL::to_double(t1.z());
			if ((x1==x2)&&(y2==y1)&&(z1==z2))
				pdbname<<t.id<<"  "<<t.chain<<"  "<<t.res<<"  "<<t.resno<<endl;
			ip++;
		}
		ir11++;
	}
	///////////////////////residue in interface//////////////////////////////////////////////////////////////
	/////////////////////////////////////
	vector<residuenum>  resi;
	int gh=0;
	vector<Point_3>::iterator ir1= atomi12.begin();
	while (ir1!=atomi12.end())
	{
		vector<atomres>::iterator irr=residue.begin();
		while (irr!=residue.end())
		{
			Point_3 t1;
			string s;
			if (irr->p==*ir1)
			{
				residuenum f;
				f.res1=irr->res ;
				f.resno1= irr->resno   ;
				f.chain=irr->chain;
				bool n=false;
				vector<residuenum>::iterator ic=resi.begin();
				while (ic!=resi.end())
				{
					if ((f.res1==ic->res1)&&(f.resno1==ic->resno1)&&(f.chain==ic->chain))
					{
						n=true;
						break;
					}
					ic++;
				}
				if (!n)
				{
					resi.push_back(f);
				}
				break;
			}
			irr++;
		}
		ir1++;
	}
	cout<<"number of residues in interface :"<<resi.size()<<endl;
	output<<"number of residues in interface :"<<resi.size()<<endl;
	///////////////////////////////////////////////////////////////////////////////////////////////
	int num_Ala=0,num_Ile=0,num_Leu=0,num_Tyr=0,num_Trp=0,num_Thr=0,num_His=0,num_Ser=0,num_Gln=0,num_Asn=0,num_Lys=0,num_Glu=0,num_Asp=0,num_Arg=0,num_Val=0,num_Met=0,num_Gly=0,num_Phe=0,num_Pro=0,num_Cys=0,num_A=0,num_C=0,num_G=0,num_U=0,num_T=0,num_DA=0,num_DG=0,num_DT=0,num_DC=0;
	vector<residuenum>::iterator ic2=resi.begin();
	while (ic2!=resi.end())
	{
		string a=ic2->res1;	
		if (a=="AALA"||a=="BALA"||a=="ALA")
			num_Ala++;
		else if (a=="ILE"||a=="AILE"||a=="BILE")
			num_Ile++;
		else if (a=="LEU"||a=="ALEU"||a=="BLEU")
			num_Leu++;
		else if (a=="TYR"||a=="ATYR"||a=="BTYR")
			num_Tyr++;
		else if (a=="TRP"||a=="ATRP"||a=="BTRP")
			num_Trp++;
		else if (a=="THR"||a=="ATHR"||a=="BTHR")
			num_Thr++;
		else if (a=="HIS"||a=="AHIS"||a=="BHIS")
			num_His++;
		else if (a=="SER"||a=="ASER"||a=="BSER")
			num_Ser++;
		else if (a=="GLN"||a=="AGLN"||a=="BGLN")
			num_Gln++;
		else if (a=="ASN"||a=="AASN"||a=="BASN")
			num_Asn++;
		else if (a=="LYS"||a=="ALYS"||a=="BLYS")
			num_Lys++;
		else if (a=="GLU"||a=="AGLU"||a=="BGLU")
			num_Glu++;
		else if (a=="ASP"||a=="AASP"||a=="BASP")
			num_Asp++;
		else if (a=="ARG"||a=="AARG"||a=="BARG")
			num_Arg++;
		else if (a=="VAL"||a=="AVAL"||a=="BVAL")
			num_Val++;
		else if (a=="MET"||a=="AMET"||a=="BMET")
			num_Met++;
		else if (a=="GLY"||a=="AGLY"||a=="BGLY")
			num_Gly++;
		else if (a=="PHE"||a=="APHE"||a=="BPHE")
			num_Phe++;
		else if (a=="PRO"||a=="APRO"||a=="BPRO")
			num_Pro++;
		else if (a=="CYS"||a=="ACYS"||a=="BCYS")
			num_Cys++;
		else if (a=="A")
			num_A++;
		else if (a=="C")
			num_C++;
		else if (a=="G")
			num_G++;
		else if (a=="U")
			num_U++;
		else if (a=="DA")
			num_DA++;
		else if (a=="DG")
			num_DG++;
		else if (a=="DT")
			num_DT++;
		else if (a=="DC")
			num_DC++;
		else if (a=="T")
			num_T++;
		else 
			cout<<a<<endl;
		ic2++;
	}
	numberres<<"num_Ala= \t"<<num_Ala<<endl;
	numberres<<"num_Ile= \t"<<num_Ile<<endl;
	numberres<<"num_Leu= \t"<<num_Leu<<endl;
	numberres<<"num_Tyr= \t"<<num_Tyr<<endl;
	numberres<<"num_Trp= \t"<<num_Trp<<endl;
	numberres<<"num_Thr= \t"<<num_Thr<<endl;
	numberres<<"num_His= \t"<<num_His<<endl;
	numberres<<"num_Ser= \t"<<num_Ser<<endl;
	numberres<<"num_Gln= \t"<<num_Gln<<endl;
	numberres<<"num_Asn= \t"<<num_Asn<<endl;
	numberres<<"num_Lys= \t"<<num_Lys<<endl;
	numberres<<"num_Glu= \t"<<num_Glu<<endl;
	numberres<<"num_Asp= \t"<<num_Asp<<endl;
	numberres<<"num_Arg= \t"<<num_Arg<<endl;
	numberres<<"num_Val= \t"<<num_Val<<endl;
	numberres<<"num_Met= \t"<<num_Met<<endl;
	numberres<<"num_Gly= \t"<<num_Gly<<endl;
	numberres<<"num_Phe= \t"<<num_Phe<<endl;
	numberres<<"num_Pro= \t"<<num_Pro<<endl;
	numberres<<"num_Cys= \t"<<num_Cys<<endl;
	numberres<<"num_A= \t"<<num_A<<endl;
	numberres<<"num_C= \t"<<num_C<<endl;
	numberres<<"num_G= \t"<<num_G<<endl;
	numberres<<"num_U= \t"<<num_U<<endl;
	numberres<<"num_T= \t"<<num_T<<endl;
	numberres<<"num_DA= \t"<<num_DA<<endl;
	numberres<<"num_DG= \t"<<num_DG<<endl;
	numberres<<"num_DT= \t"<<num_DT<<endl;
	numberres<<"num_DC= \t"<<num_DC<<endl;
	///////////////////////////////////////////////////
	//////////////////////////////////////////volume//////////////////////////
	std::vector<K::Point_3>  atom11=atomi12;
	////////////////buird atom in interface or no ///////////////
	std::vector<K::Point_3>  atomb ;
	vector<K::Point_3>::iterator it2=atom11.begin();
	while (it2!=atom11.end())
	{
		bool t=true;
		Edge_iterator ei1 = T.finite_edges_begin();
		for( ; ei1 != T.finite_edges_end(); ++ei1) 
		{ 
			Segment_3 s;
			s=T.segment(*ei1);
			if (s.source()==*it2|| s.target()==*it2)
			{
				if ((find(atom11.begin(), atom11.end(),s.source()) == atom11.end()) ||(find(atom11.begin(), atom11.end(), s.target()) == atom11.end()))
				{
					t=false;
					break;
				}
			}
		}//end for
		if (t)
		{
			atomb.push_back(*it2);
		}
		//else
		it2++;
	}////end while it2
	//////////////////////////remove water in atomb/////////////////////////////////////////////
	vector<K::Point_3>::iterator iatomb1=atomb.begin();
	while (iatomb1!=atomb.end())
	{
		K::Point_3 t=*iatomb1;
		if (find(waterpoint.begin(),waterpoint.end(),t) != waterpoint.end())
			atomb.erase(iatomb1);
		else
			iatomb1++;
	}
	////////buired atom in only protein/////////////////////////////////////////////
	std::cout<<"if complex is protein-protein enter number 1 else enter number 2 "<<std::endl;
	int aa;
	std::cin>>aa;
	if (aa!=1)
	{
		vector<K::Point_3>::iterator iatomb1h=atomb.begin();
		while (iatomb1h!=atomb.end())
		{
			K::Point_3 t=*iatomb1h;
			if (find(co2.begin(),co2.end(),t) == co2.end())
				atomb.erase(iatomb1h);
			else
				iatomb1h++;
		}
	}
	///////////////////////////////////////////
	vector<K::Point_3>::iterator iatomb12=atomb.begin();
	while (iatomb12!=atomb.end())
	{
		K::Point_3 t=*iatomb12;
		pdbatom<<t<<endl;
		iatomb12++;
	}
	///////////////////////////////////////////////////
	//////////////////////////////////
	vector<Vertex_handle>    N;
	for(  Finite_vertices_iterator it = T.finite_vertices_begin();
		it != T.finite_vertices_end(); ++ it )
	{
		Point_3 d=it->point();
		if ((find(atomb.begin(), atomb.end(), d) != atomb.end()))
			N.push_back(it);
	} 
	///////////////////////////////////////////////////////////////////
	vector<Vertex_handle>::iterator it1=N.begin();
	double  v=0;
	while (it1!=N.end())
	{
		std::vector<K::Point_3>  r;
		std::vector < Cell_handle>   cell_cir;
		Finite_cells_iterator fc = T.finite_cells_begin();
		for( ; fc != T.finite_cells_end(); ++fc) 
		{
			Point_3  v11=fc->vertex(0)->point();
			Point_3 v2=fc->vertex(1)->point();
			Point_3 v3=fc->vertex(2)->point();
			Point_3 v4=fc->vertex(3)->point();
			if ((v11==(*it1)->point())||(v2==(*it1)->point())||(v3==(*it1)->point())||(v4==(*it1)->point()))
				cell_cir.push_back(fc);
		}
		std::vector< Cell_handle>::iterator ity=cell_cir.begin();
		while (ity!=cell_cir.end())
		{
			K::Point_3 t1;
			t1=T.dual(*ity);
			r.push_back(t1);   
			ity++;
		}
		///////////////create polyhedron
		//std::cout<<"size r  "<<r.size()<<std::endl;
		Polyhedron_3 Polyb;
		CGAL::convex_hull_3(r.begin(), r.end(),Polyb );
		std::for_each(  Polyb.facets_begin(),    Polyb.facets_end(),   Facet_normal());
		std::for_each(  Polyb.vertices_begin(),  Polyb.vertices_end(), Vertex_normal());
		std::vector<K::Point_3>  pointvertex;
		for ( Vertex_iterator icb = Polyb.vertices_begin(); icb!= Polyb.vertices_end(); ++icb)
			pointvertex.push_back(icb-> point());
		Triangulation Tb1;
		Tb1.insert (pointvertex.begin(), pointvertex.end());
		assert( Tb1.is_valid() );
		assert( Tb1.dimension() == 3 );
		Finite_cells_iterator1 fc1 = Tb1.finite_cells_begin();
		double v1=0;
		for( ; fc1 != Tb1.finite_cells_end(); ++fc1) 
		{
			//Tetrahedron_3 tet= Tb1.tetrahedron(*fc1);
			Point_3 v10=fc1->vertex(0)->point();
			Point_3 v2=fc1->vertex(1)->point();
			Point_3 v3=fc1->vertex(2)->point();
			Point_3 v4=fc1->vertex(3)->point();
			Tetrahedron_3 tet(v10,v2,v3,v4);
			v1=(CGAL::to_double(tet.volume()))+v1; 
		}
		v=v+v1;
		//////////////////////////////////////////
		it1++;
	}
	//////////////////////////////////////////////
	//////////////////////////////////////////////
	////////////////////calculate standar volume 
	double sv0=0;
	ifstream pdbif("interfaceatom.txt");
	vector<interfaceatom>  ifatom11;
	while (!pdbif.eof())
	{
		double x1,y1,z1;
		char line[400];
		pdbif.getline(line,400);
		string d=line;
		int size=int(d.size());
		vector<string> words;
		int count_word=0;
		int i=1;
		while (i<=(size))
		{
			if (d[1]==' ')
				break;
			if (((d[i]==' ')&&(d[i-1]!=' '))||i==size)
			{
				string a;
				int m;
				m=i-1;
				while (d[m]!=' ')
				{
					a.insert(a.end(),d[m]);
					if (m>0)
						m--;
					else 
						break;
				}//end while m;
				string c;
				for(int k=int(a.size())-1;k>=0;k--)
					c.insert(c.end(),a[k]); 
				words.push_back(c);
				count_word++;
			}//end word
			i++;
		}//end while i line;
		if (!words.empty())
		{
			x1=boost::lexical_cast<double>(*(words.begin()));
			y1=boost::lexical_cast<double>(*(words.begin()+1));
			z1=boost::lexical_cast<double>(*(words.begin()+2));
		}
		////////////compare by all of atoms in complex
		vector<atoms>::iterator ip=protein.begin();
		while (ip!=protein.end())
		{
			atoms t=*ip;
			double x2,y2,z2;
			x2=t.x;
			y2=t.y;
			z2=t.z;
			if ((x1==x2)&&(y2==y1)&&(z1==z2))
			{
				interfaceatom  it;
				it.res1=t.res;
				it.resno1=t.resno;
				it.chain=t.chain;
				it.nameatom=t.id;
				it.x=t.x;
				it.y=t.y;
				it.z=t.z;
				bool y=true;
				vector<interfaceatom >::iterator is=ifatom11.begin();
				while(is!=ifatom11.end())
				{
					interfaceatom t=*is;
					if((t.x==it.x)&&(t.y==it.y)&&(t.z==it.z))
						y=false;
					is++;
				}
				if (y)
					ifatom11.push_back(it);
			}
			ip++;
		}
	}///end while file
	//////////////////////////////////
	std::cout<<ifatom11.size()<<std::endl;
	output<<ifatom11.size()<<std::endl;
	vector<interfaceatom>::iterator   iatom1=ifatom11.begin();
	while (iatom1!=ifatom11.end())
	{
		interfaceatom  t=*iatom1;
		double v0=0;
		if (t.res1=="GLY")
		{
			if  (t.nameatom=="N")
				v0=14.5;
			else if (t.nameatom=="CA")
				v0=23.3;
			else if (t.nameatom=="C")
				v0=9.5;
			else if (t.nameatom=="O")  
				v0=16.2;
		}//end GLY
		else if (t.res1=="ALA")
		{
			if  (t.nameatom=="N")
				v0=13.8;
			else if (t.nameatom=="CA")
				v0=14;
			else if (t.nameatom=="C")
				v0=8.8;
			else if (t.nameatom=="O")  
				v0=16;
			else if (t.nameatom=="CB")  
				v0=36.9;
		}//end ALA
		else if (t.res1=="VAL")
		{
			if  (t.nameatom=="N")
				v0=13.6;
			else if (t.nameatom=="CA")
				v0=13.1;
			else if (t.nameatom=="C")
				v0=8.5;
			else if (t.nameatom=="O")  
				v0=16;
			else if (t.nameatom=="CB")  
				v0=14.5;
			else if (t.nameatom=="CG1")  
				v0=36.1;
			else if (t.nameatom=="CG2")  
				v0=36;
		}//end vAL
		else if (t.res1=="LEU")
		{
			if  (t.nameatom=="N")
				v0=13.7;
			else if (t.nameatom=="CA")
				v0=13.1;
			else if (t.nameatom=="C")
				v0=8.7;
			else if (t.nameatom=="O")  
				v0=16;
			else if (t.nameatom=="CB")  
				v0=22.8;
			else if (t.nameatom=="CG")  
				v0=14.6;
			else if (t.nameatom=="CD1")  
				v0=37.4;
			else if (t.nameatom=="CD2")  
				v0=36.8;
		}//end leu
		else if (t.res1=="ILE")
		{
			if  (t.nameatom=="N")
				v0=13.6;
			else if (t.nameatom=="CA")
				v0=12.9;
			else if (t.nameatom=="C")
				v0=8.5;
			else if (t.nameatom=="O")  
				v0=16.1;
			else if (t.nameatom=="CB")  
				v0=14.2;
			else if (t.nameatom=="CG1")  
				v0=24.1;
			else if (t.nameatom=="CG2")  
				v0=35.7;
			else if (t.nameatom=="CD1")  
				v0=38.5;
		}//end ile
		else if (t.res1=="MET")
		{
			if  (t.nameatom=="N")
				v0=13.5;
			else if (t.nameatom=="CA")
				v0=13.2;
			else if (t.nameatom=="C")
				v0=8.8;
			else if (t.nameatom=="O")  
				v0=15.9;
			else if (t.nameatom=="CB")  
				v0=23.5;
			else if (t.nameatom=="CG")  
				v0=24.0;
			else if (t.nameatom=="SD")  
				v0=29.6;
			else if (t.nameatom=="CE")  
				v0=36.8;
		}//end MET
		else if (t.res1=="PRO")
		{
			if  (t.nameatom=="N")
				v0=8.5;
			else if (t.nameatom=="CA")
				v0=13.9;
			else if (t.nameatom=="C")
				v0=8.7;
			else if (t.nameatom=="O")  
				v0=16.3;
			else if (t.nameatom=="CB")  
				v0=25.3;
			else if (t.nameatom=="CG")  
				v0=25.8;
			else if (t.nameatom=="CD")  
				v0=23.6;
		}//end pro
		else if (t.res1=="HIS")
		{
			if  (t.nameatom=="N")
				v0=13.6;
			else if (t.nameatom=="CA")
				v0=13.3;
			else if (t.nameatom=="C")
				v0=8.8;
			else if (t.nameatom=="O")  
				v0=16.2;
			else if (t.nameatom=="CB")  
				v0=23.2;
			else if (t.nameatom=="CG")  
				v0=9.9;
			else if (t.nameatom=="ND1")  
				v0=15.1;
			else if (t.nameatom=="CD2")  
				v0=21.3;
			else if (t.nameatom=="CE1")  
				v0=20.8;
			else if (t.nameatom=="NE2")  
				v0=14.9;
		}//end his
		else if (t.res1=="PHE")
		{
			if  (t.nameatom=="N")
				v0=13.6;
			else if (t.nameatom=="CA")
				v0=13.4;
			else if (t.nameatom=="C")
				v0=8.6;
			else if (t.nameatom=="O")  
				v0=16.0;
			else if (t.nameatom=="CB")  
				v0=23.6;
			else if (t.nameatom=="CG")  
				v0=9.7;
			else if (t.nameatom=="CD1")  
				v0=20.3;
			else if (t.nameatom=="CD2")  
				v0=21.0;
			else if (t.nameatom=="CE1")  
				v0=21.4;
			else if (t.nameatom=="CE2")  
				v0=21.5;
			else if (t.nameatom=="CZ")  
				v0=21.6;
		}//end phe
		else if (t.res1=="TYR")
		{
			if  (t.nameatom=="N")
				v0=13.5;
			else if (t.nameatom=="CA")
				v0=13.3;
			else if (t.nameatom=="C")
				v0=8.8;
			else if (t.nameatom=="O")  
				v0=16.0;
			else if (t.nameatom=="CB")  
				v0=23.4;
			else if (t.nameatom=="CG")  
				v0=9.6;
			else if (t.nameatom=="CD1")  
				v0=20.0;
			else if (t.nameatom=="CD2")  
				v0=20.8;
			else if (t.nameatom=="CE1")  
				v0=20.4;
			else if (t.nameatom=="CE2")  
				v0=20.3;
			else if (t.nameatom=="CZ")  
				v0=9.8;
			else if (t.nameatom=="OH")  
				v0=18.8;
		}//end TYR
		else if (t.res1=="TRP")
		{
			if  (t.nameatom=="N")
				v0=13.7;
			else if (t.nameatom=="CA")
				v0=13.2;
			else if (t.nameatom=="C")
				v0=8.7;
			else if (t.nameatom=="O")  
				v0=16.0;
			else if (t.nameatom=="CB")  
				v0=24.0;
			else if (t.nameatom=="CG")  
				v0=9.8;
			else if (t.nameatom=="CD1")  
				v0=20.5;
			else if (t.nameatom=="NE1")  
				v0=16.7;
			else if (t.nameatom=="CD2")  
				v0=10.1;
			else if (t.nameatom=="CE2")  
				v0=9.7;
			else if (t.nameatom=="CE3")  
				v0=20.3;
			else if (t.nameatom=="CZ3")  
				v0=21.4;
			else if (t.nameatom=="CZ2")  
				v0=20.7;
			else if (t.nameatom=="CH2")  
				v0=21.1;
		}//end TRp
		else if (t.res1=="CYH")
		{
			if  (t.nameatom=="N")
				v0=14;
			else if (t.nameatom=="CA")
				v0=13.4;
			else if (t.nameatom=="C")
				v0=8.8;
			else if (t.nameatom=="O")  
				v0=16.3;
			else if (t.nameatom=="CB")  
				v0=23.1;
			else if (t.nameatom=="SG")  
				v0=35.3;
		}//end cyh
		////////////////
		else if (t.res1=="CYS")
		{
			if  (t.nameatom=="N")
				v0=13.8;
			else if (t.nameatom=="CA")
				v0=13.2;
			else if (t.nameatom=="C")
				v0=8.7;
			else if (t.nameatom=="O")  
				v0=16.3;
			else if (t.nameatom=="CB")  
				v0=23.4;
			else if (t.nameatom=="SG")  
				v0=27.7;
		}//end cys
		else if (t.res1=="SER")
		{
			if  (t.nameatom=="N")
				v0=13.9;
			else if (t.nameatom=="CA")
				v0=13.4;
			else if (t.nameatom=="C")
				v0=8.8;
			else if (t.nameatom=="O")  
				v0=16.1;
			else if (t.nameatom=="CB")  
				v0=23.7;
			else if (t.nameatom=="OG")  
				v0=18.9;
		}//end ser
		else if (t.res1=="THR")
		{
			if  (t.nameatom=="N")
				v0=13.6;
			else if (t.nameatom=="CA")
				v0=13.2;
			else if (t.nameatom=="C")
				v0=8.7;
			else if (t.nameatom=="O")  
				v0=16.0;
			else if (t.nameatom=="CB")  
				v0=14.8;
			else if (t.nameatom=="OG1")  
				v0=17.9;
			else if (t.nameatom=="CG2")  
				v0=35.8;
		}//end thr
		else if (t.res1=="ASN")
		{
			if  (t.nameatom=="N")
				v0=13.9;
			else if (t.nameatom=="CA")
				v0=13.1;
			else if (t.nameatom=="C")
				v0=8.9;
			else if (t.nameatom=="O")  
				v0=16.3;
			else if (t.nameatom=="CB")  
				v0=23.4;
			else if (t.nameatom=="CG")  
				v0=9.5;
			else if (t.nameatom=="OD1")  
				v0=16.8;
			else if (t.nameatom=="ND2")  
				v0=23.3;
		}//end asn
		else if (t.res1=="GLN")
		{
			if  (t.nameatom=="N")
				v0=13.7;
			else if (t.nameatom=="CA")
				v0=13.3;
			else if (t.nameatom=="C")
				v0=8.7;
			else if (t.nameatom=="O")  
				v0=15.9;
			else if (t.nameatom=="CB")  
				v0=23.2;
			else if (t.nameatom=="CG")  
				v0=23.6;
			else if (t.nameatom=="CD")  
				v0=9.7;
			else if (t.nameatom=="OE1")  
				v0=17.1;
			else if (t.nameatom=="NE2")  
				v0=23.6;
		}//end gln
		else if (t.res1=="ASP")
		{
			if  (t.nameatom=="N")
				v0=13.9;
			else if (t.nameatom=="CA")
				v0=13.2;
			else if (t.nameatom=="C")
				v0=8.7;
			else if (t.nameatom=="O")  
				v0=16.3;
			else if (t.nameatom=="CB")  
				v0=23.0;
			else if (t.nameatom=="CG")  
				v0=9.2;
			else if (t.nameatom=="OD1")  
				v0=15.3;
			else if (t.nameatom=="OD2")  
				v0=15.9;
		}//end asp
		else if (t.res1=="GLU")
		{
			if  (t.nameatom=="N")
				v0=13.7;
			else if (t.nameatom=="CA")
				v0=13.4;
			else if (t.nameatom=="C")
				v0=8.6;
			else if (t.nameatom=="O")  
				v0=15.7;
			else if (t.nameatom=="CB")  
				v0=23.3;
			else if (t.nameatom=="CG")  
				v0=23.0;
			else if (t.nameatom=="CD")  
				v0=9.4;
			else if (t.nameatom=="OE1")  
				v0=15.4;
			else if (t.nameatom=="OE2")  
				v0=15.8;
		}//end glu
		else if (t.res1=="LYS")
		{
			if  (t.nameatom=="N")
				v0=13.3;
			else if (t.nameatom=="CA")
				v0=13.4;
			else if (t.nameatom=="C")
				v0=8.7;
			else if (t.nameatom=="O")  
				v0=16.2;
			else if (t.nameatom=="CB")  
				v0=23.0;
			else if (t.nameatom=="CG")  
				v0=23.7;
			else if (t.nameatom=="CD")  
				v0=24.1;
			else if (t.nameatom=="CE")  
				v0=23.3;
			else if (t.nameatom=="NZ")  
				v0=21.7;
		}//end lys
		else if (t.res1=="ARG")
		{
			if  (t.nameatom=="N")
				v0=13.6;
			else if (t.nameatom=="CA")
				v0=13.4;
			else if (t.nameatom=="C")
				v0=8.8;
			else if (t.nameatom=="O")  
				v0=16.0;
			else if (t.nameatom=="CB")  
				v0=22.9;
			else if (t.nameatom=="CD")  
				v0=23.3;
			else if (t.nameatom=="CG")  
				v0=23.2;
			else if (t.nameatom=="NE")  
				v0=15.2;
			else if (t.nameatom=="CZ")  
				v0=9.6;
			else if (t.nameatom=="NH1")  
				v0=22.3;
			else if (t.nameatom=="NH2")  
				v0=23.2;
		}//end arg
		sv0=sv0+v0;
		iatom1++;
	}
	double indexv= v/sv0 ;
	std::cout<<"index v==   "<<indexv<<std::endl;
	output<<"index v==   "<<indexv<<std::endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////
	std::cin>>ch ;
	return 0;
}

















