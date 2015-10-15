/**
            Implementation of Birch Algorithm specific for an indexing module of speech feature vectors.
			Code written by:  Murali Raghu Babu B 
			Computer Science and Engineering Department
			Indian Institute of Technology Guwahati
			India
			Year: 2014
			Refer to this paper for the Birch Algorithm http://www.cs.sfu.ca/CourseCentral/459/han/papers/zhang96.pdf
			This code is of the Birch Algorithm along with the merging requirement mentioned to avoid bad clustering in the presence of skewed data.
			For any queries feel free to get back to me at muraliraghubabu1994@gmail.com			

**/

//including the necessary header files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <limits>
#include <stdlib.h>
#include <iomanip>

using namespace std;

	class cfpoint;						//pre-declaration of the classes to avoid errors and warnings
	class cfentry;
	class cfnode;
	class cftree;

	int dim=0;							//global variable - dimesnion of the feature vectors
	int b=0;							//no of entries in the non leaf nodes after the birch tree
	int l=0;							//no of entries in the leaf nodes
	double threshold = 0;				//threshold which the leaf nodes should not violate
	cfnode* leafliststart=NULL;			//pointer to the starting leaf node
	cfentry* rootsplitentry=NULL;		//the entry which 
	int totalentries=0;					//total no of entries in the tree
	int totalnodes=0;					//total no of nodes in the tree
	int totalpoint=0;					//total no of data points


    double dtw(cfentry* e1, cfentry* e2);			//calculates the dtw distance between the two feature vectors

    double euclidean(double *sum1, double* sum2)	//calculates the euclidean distance between the two vectors
    {
    int i=0; double distance=0;
    while(i<dim)
    {
        distance+=(sum1[i]-sum2[i])*(sum1[i]-sum2[i]);
        i++;
    }
    return sqrt(distance);				//returns the distance between the two vectors
	}

	class cfpoint									//class for each point
	{
	private:
		int id;					//the id no of the feature vector
		double *sumx;			//the actual data point
		string filen;			//the filename associated with the feature vector
	public:
		friend class cfentry;			//friend classes to utilize the functions here
		friend class cfnode;
		friend class cftree;
		cfpoint(int i, string s,double *sum)		//constructor to initialize the point
		{
			id=i;
			filen=s;
			sumx=new double[dim];
			i=0;
			while(i<dim)
			{
				sumx[i]=sum[i];
				i++;
			}
		}
		cfpoint(double *sum)						//constructor to initialize the point - constructor overloading
		{
			int i=0;
			sumx=new double[dim]; id=0; filen="";
			while(i<dim)
			{
				sumx[i]=sum[i];
				i++;
			}
		}
		cfpoint()									//basic constructor
		{
			id=0;
			sumx=new double[dim];
			filen="";
		}
		double* getxsum()							//returns the feature vector
		{
			return sumx;
		}
		double getsumx(int i)						//returns the value of the feature vector at that dimension
		{
			return sumx[i];
		}
		string getfilen()							//returns the filename of that vector
		{
			return filen;
		}
		int getid()									//returns the id of this feature vector
		{
			return id;
		}
		void setsumx(int i, double a)				//changing the value of that feature vector at that dimension
		{
			sumx[i]=a;
		}
	};

	class cfentry										//class for the entry in each node
	{
	private:
		cfnode* child;								//pointer to the child of the node it is in
		cfnode* nodeparent;							//pointer to the parent of the node it is in
		int *id;									//pointer to the id vector
		long int n, maxval;							//no of points in this entry
		bool inleaf;								//whether this entry is in a leafnode or not
		cfpoint **comparepoints;					//pointer to the set of cfpoints that represents this entry - used for distance calculation
		int ncompare;								//no of the compare points
		int maxcmp;									//max compare points

	public:
		friend class cfpoint;											//friend classes to utilize the functions in it
		friend class cfnode;
		friend class cftree;
		friend void searchtimitquery(cfentry* query, cfnode* root);		//friend function to utilize the variables
		cfentry()														//constructor to initilize it
		{
			child=NULL;													//initializing all the parameters in the entry.
			nodeparent=NULL;
			id = new int[1];
			n=0;
			maxval=1;
			inleaf=true;
			comparepoints = new cfpoint*[1];
			ncompare=0;
			maxcmp=1;
		}
		cfnode* getchild()												//get the child node
		{
			return child;
		}
		cfnode* getnodeparent()											//get the parent node
		{
			return nodeparent;
		}
		cfpoint* getcomparepoints(int i)								// get the respective compare point
		{
			return comparepoints[i];
		}
		int getncompare()												// get the ncompare value
		{
			return ncompare;
		}
		int getid(int i)												//get the id of the ith point
		{
			return id[i];
		}
		int getn()														//get the no of points in this entry			
		{
			return n;
		}
		bool getinleaf()												//get to know if in leaf or not
		{
			return inleaf;
		}
		void addthispoint(cfpoint *p, int q)			//adding the cfpoint in that entry.
		{
			if(p)										// if it is not null
			{
				comparepoints[ncompare]=p;
				ncompare++;
				int i=0;
				if(ncompare==maxcmp)
				{
					cfpoint **newp = new cfpoint* [2*maxcmp];
					maxcmp = 2*maxcmp;

					i=0;
					while(i<ncompare)
					{
						newp[i]=comparepoints[i];
						i++;
					}
					free(comparepoints);
					comparepoints=newp;
				}

				if(n==0)
				{
					id[n]=q;
					n++;
					if(n==maxval)
					{
						int *newid= new int[2*maxval];
						maxval=2*maxval;
						i=0;
						while(i<n)
						{
							newid[i]=id[i];
							i++;
						}
						free(id);
						id=newid;
					}
				}

			}
		 }
		bool checkthreshold(cfentry *e)					//checks the threshold radius to see if the cluster is not violating the threshold constraint
		{

			if(getupdatedradius(e)>threshold)
				{   return false;   }
		return true;
		}
		double getupdatedradius(cfentry *e)					//get the distance between the current entry and the entry e
		{
			 if(ncompare!=0)
			{
					double radius=0;
					radius=dtw(this,e);
				return radius;
			}
			else
			{
				return 0.0;
			}
		}
		void modifycomparepoints(cfentry *e1, cfentry *e2)		//modify the compare points by combining these two entries
		{
			int a=e2->getncompare();
			int b=e1->getncompare();

			if(a!=0&&b!=0) {
			double matrix[a][b]; int i=0, j=0; double dmin=0; double d1, d2, d3;
			double row[a][b];

			while(j<b)
			{
				i=0;
				while(i<a)
				{
					d1=std::numeric_limits<double>::max();
					if(i>=0&&j-1>=0)
					{
						d1=matrix[i][j-1];
					}

					d2=std::numeric_limits<double>::max();
					if(i-1>=0&&j-1>=0)
					{
						d2=matrix[i-1][j-1];
					}

					d3=std::numeric_limits<double>::max();
					if(i-2>=0&&j-1>=0)
					{
						d3=matrix[i-2][j-1];
					}


					if(d1<=d2&&d1<=d3)
					{
						dmin=d1;
						if(d1!=std::numeric_limits<double>::max())
						{
						row[i][j]=i;
						//column[i][j-1]=j-1;
						}
						else
						{
							dmin=0;
							row[i][j]=-1;
						   // column[i][j]=-1;
						}
					}
					else if(d2<=d1&&d2<=d3)
					{
						dmin=d2;
						if(d2!=std::numeric_limits<double>::max())
						{
						row[i][j]=i-1;
					   // column[i][j]=j-1;
						}
						else
						{
							dmin=0;
							row[i][j]=-1;
						   // column[i][j]=-1;
						}
					}
					else
					{
						dmin=d3;
						if(d3!=std::numeric_limits<double>::max())
						{
						row[i][j]=i-2;
						//column[i][j]=j-1;
						}
						else
						{
							dmin=0;
							row[i][j]=-1;
							//column[i][j]=-1;
						}
					}


					matrix[i][j]=euclidean(e1->getcomparepoints(j)->getxsum(),e2->getcomparepoints(i)->getxsum())+ dmin;
					i++;
				}
				j++;
			}


			i=a-1; j=b-1; int k=0;
			while(j>=0)
			{
				k=0;
				while(k<dim)
				{
					e1->comparepoints[j]->setsumx(k,(e1->getcomparepoints(j)->getsumx(k) + e2->getcomparepoints(i)->getsumx(k))/2);
					k++;
				}

				if(row[i][j]!=-1) i=row[i][j];
				j--;
			}
			}

		}
		bool addentry(cfentry *e)								//adding entry in this entry
		{
			if(checkthreshold(e))
			{
				int i=0,j=0;

				if(ncompare==0)
				{
					i=0;
					while(i<e->ncompare)
					{
						comparepoints[ncompare]=e->comparepoints[i];
						ncompare++;
						if(ncompare==maxcmp)
						{
							cfpoint **newcompare;
							maxcmp=2*maxcmp;
							newcompare=new cfpoint *[maxcmp];
							j=0;
							while(j<ncompare)
							{
								newcompare[j]=comparepoints[j];
								j++;
							}
							free(comparepoints);
							comparepoints=newcompare;
						}
						i++;
					}
				}
				else
				{
					modifycomparepoints(this, e);
				}
				i=0;
				while(i<e->n)
				{
					id[n]=e->id[i];
					n++;
					if(n==maxval)
					{
						int *newid;
						maxval=2*maxval;
						newid = new int[maxval] ;
						j=0;
						while(j<n)
						{
							newid[j]=id[j]; j++;
						}
						free(id);
						id=newid;
					}
					i++;
				}



			return true;
			}
			return false;
		}
		void updateentry(cfentry *e);					//defined below
		void updateentrysame(cfentry *e);
		void updateentryfull();
	};

    double dtw(cfentry* e1, cfentry* e2)				//returns the dynamic time warping distance between these two entries.
    {
        int a=e2->getncompare();
        int b=e1->getncompare();


        double matrix[a][b]; int i=0, j=0; double dmin=0; double d1, d2, d3;
        double row[a][b];


        while(j<b)
        {
            i=0;
            while(i<a)
            {
                d1=std::numeric_limits<double>::max();
                if(i>=0&&j-1>=0)
                {
                    d1=matrix[i][j-1];
                }

                d2=std::numeric_limits<double>::max();
                if(i-1>=0&&j-1>=0)
                {
                    d2=matrix[i-1][j-1];
                }

                d3=std::numeric_limits<double>::max();
                if(i-2>=0&&j-1>=0)
                {
                    d3=matrix[i-2][j-1];
                }


                if(d1<=d2&&d1<=d3)
                {
                    dmin=d1;
                    if(d1!=std::numeric_limits<double>::max())
                    {
                    row[i][j]=i;
                    //column[i][j-1]=j-1;
                    }
                    else
                    {
                        dmin=0;
                        row[i][j]=-1;
                       // column[i][j]=-1;
                    }
                }
                else if(d2<=d1&&d2<=d3)
                {
                    dmin=d2;
                    if(d2!=std::numeric_limits<double>::max())
                    {
                    row[i][j]=i-1;
                   // column[i][j]=j-1;
                    }
                    else
                    {
                        dmin=0;
                        row[i][j]=-1;
                       // column[i][j]=-1;
                    }
                }
                else
                {
                    dmin=d3;
                    if(d3!=std::numeric_limits<double>::max())
                    {
                    row[i][j]=i-2;
                    //column[i][j]=j-1;
                    }
                    else
                    {
                        dmin=0;
                        row[i][j]=-1;
                        //column[i][j]=-1;
                    }
                }


                matrix[i][j]=euclidean(e1->getcomparepoints(j)->getxsum(),e2->getcomparepoints(i)->getxsum())+ dmin;
                i++;
            }
            j++;
        }
        //row[a-1][b-1]=a-1;
        return matrix[a-1][b-1];

    }


	class cfentrypair								//class for an entry pair basically useful for comparing for the merging refinement
	{

	private:
		cfentry *e1;
		int index1;
		int index2;
		cfentry *e2;
	public:
		friend class cfpoint;
		friend class cfentry;
		friend class cfnode;
		friend class cftree;
		cfentrypair()
		{
			e1=NULL;
			e2=NULL;
		}
		cfentrypair(cfentry *e1, cfentry *e2, int i1,int i2)
		{
			this->e1=e1;
			this->e2=e2;
			index1=i1;
			index2=i2;
		}
		bool equals(cfentrypair *e)
		{
			if((this->e1==e->e1&&this->e2==e->e2&&index1==e->index1&&index2==e->index2))
				return true;
			if((this->e1==e->e2&&this->e2==e->e1&&index1==e->index2&&index2==e->index1))
				return true;
		return false;
		}
	};


	class cfnode
	{

	private:
		cfentry **entries;
		bool isleaf;
		int countentries;
		int maxentries;
		cfentry* parent;
		cfnode* previous;
		cfnode* next;

	public:
		friend class cfpoint;
		friend class cfentry;
		friend class cftree;
		friend cfentry* getclosestentry(cfentry* p, cfnode* x);
		friend void searchtimitquery(cfentry* query, cfnode* root);
		cfnode()
		{
			entries = new cfentry *[b];
			isleaf=true;
			countentries=0;
			maxentries=b;
			next=NULL;
			parent=NULL;
			previous=NULL;
		}
		cfnode* getnext()
		{
			return next;
		}
		cfnode* getprevious()
		{
			return previous;
		}
		cfentry* getparent()
		{
			return parent;
		}
		bool getisleaf()
		{
			return isleaf;
		}
		cfentry* getentries(int i)
		{
			return entries[i];
		}
		int getcountentries()
		{
			return countentries;
		}
		int getmaxentries()
		{
			return maxentries;
		}
		cfentry* getentry(int i)
		{
			return entries[i];
		}
		bool addentry(cfentry *e)
		{
			if(countentries!=0)
			{
				cfentry *closestentry;
				closestentry=this->getclosestentry(e);

				if(closestentry->addentry(e))
					{
						cfentry *x=parent;
						while(x!=NULL)
						{
							x->updateentry(e);
							x=x->nodeparent->parent;

						}
					return true;
					}
				else
				{

					if(countentries<maxentries)
					{
						entries[countentries]=e;
						e->nodeparent=this;
						e->inleaf=this->isleaf;
						countentries++;
						cfentry *x=parent;
						while(x!=NULL)
						{
							x->updateentry(e);
							x=x->nodeparent->parent;
						}
					return true;
					}
					else
					{
						if(parent)
						{
							if(split(e))
							{   parent=NULL;
								previous=NULL;
								next=NULL;
								/// here something free(this) like that;;;;; ///////////////////
							return true;
							}

							else {
								return false;
								}
						}
						else
						{
							rootsplitentry=e;
							return false;
						}
					}
				}
			}
			else
			{
				entries[countentries]=e;
				e->nodeparent=this;
				e->inleaf=this->isleaf;
				countentries++;
				cfentry *x=parent;
				while(x!=NULL)
					{
						x->updateentry(e);
						x=x->nodeparent->parent;
					}
			return true;
			}

		return false;
		}
		cfentry* getclosestentry(cfentry *e)
		{
				int i=0, t=0; double dist=0;

				dist=dtw(e,entries[i]);

				i++;
				while(i<countentries)
				{
					if(dist>dtw(e,entries[i]))
					{   dist=dtw(e,entries[i]);
						t=i;
					}
					i++;
				}
				if(countentries==0) return NULL;

		return entries[t];
		}
		cfentrypair* getfarthestentries()
		{
			int i=0,j=1; double maxdist=0;
			maxdist=dtw(entries[i],entries[j]);
			int e1=0,e2=1;
			while(i<countentries)
			{
				j=i+1;
				while(j<countentries)
				{
					if(maxdist<dtw(entries[i],entries[j]))
					{
						maxdist=dtw(entries[i],entries[j]);
						e1=i;
						e2=j;
					}
					j++;
				}
				i++;
			}
			cfentrypair* farthestpair=new cfentrypair(entries[e1],entries[e2],e1,e2);

		return farthestpair;
		}
		bool split(cfentry *entry)
		{
			cfentrypair* farthestpair= getfarthestentries();

			cfentry * newe1=new cfentry;
			cfnode* node1=new cfnode;
			cfentry* newe2=new cfentry;
			cfnode* node2=new cfnode;
			newe1->child=node1;
			node1->parent=newe1;
			newe2->child=node2;
			node2->parent=newe2;
			node1->maxentries=maxentries;
			node2->maxentries=maxentries;
			node1->entries = new cfentry *[maxentries];
			node2->entries = new cfentry *[maxentries];
			node1->isleaf=isleaf;
			node2->isleaf=isleaf;
			newe1->inleaf=false;
			newe2->inleaf=false;


			if(isleaf)
			{
				node1->previous=previous;
				if(previous)
					previous->next=node1;
				else
					leafliststart=node1;

				node1->next=node2;
				node2->next=next;
				node2->previous=node1;
				if(next)
					next->previous=node2;

			}


			node1->countentries=0; node2->countentries=0;

			redistributeentries(farthestpair,newe1,newe2, node1, node2,entry);

			newe1->child=node1;
			node1->parent=newe1;
			newe2->child=node2;
			node2->parent=newe2;

			if(parent)
			{

				cfnode* earlynode=parent->nodeparent;
				cfentry *removeentry=parent;

				int hd=replaceentryearlynode(earlynode,removeentry,newe1);
				if(addentryearlynode(earlynode,newe2,newe1,hd))
					{ /// freenode=this; /// free node this look here....////////////
						return true;
					}
				else
					{   return false; }

			}
			else
			{
					rootsplitentry=entry;
					return false;
			}
		}
		void redistributeentries(cfentrypair* farthestpair,cfentry *newe1, cfentry* newe2,cfnode* node1, cfnode* node2, cfentry* entry)
		{
			int i=node1->countentries; int j=node2->countentries;
			node1->entries[i]=farthestpair->e1;
			node2->entries[j]=farthestpair->e2;
			node1->entries[i]->nodeparent=node1;
			node2->entries[j]->nodeparent=node2;
			node1->entries[i]->inleaf=node1->isleaf;
			node2->entries[j]->inleaf=node2->isleaf;

			newe1->child=node1;
			node1->parent=newe1;
			newe2->child=node2;
			node2->parent=newe2;

			int index1=farthestpair->index1;
			int index2=farthestpair->index2;
			i++; node1->countentries=i;
			j++; node2->countentries=j;

			newe1->updateentrysame(farthestpair->e1);

			newe2->updateentrysame(farthestpair->e2);

			int k=0;

	/// ///////////////////// here l cannot be equal to 1 if you want it to be 1 change the code here..!!!

			while(k<countentries)
			{
				if(k!=index1&&k!=index2)
				{
					if(dtw(newe1,entries[k])<=dtw(newe2,entries[k]))
					   {

							if(i<node1->maxentries)
						  {
							node1->entries[i]=entries[k];
							entries[k]->inleaf=node1->isleaf;
							node1->entries[i]->nodeparent=node1;
							i++;
							node1->countentries=i;
							newe1->updateentry(entries[k]);
							//new cfpoint**
						  }
						  else
						  {
							  node2->entries[j]=entries[k];
							entries[k]->inleaf=node2->isleaf;
							node2->entries[j]->nodeparent=node2;
							j++;
							node2->countentries=j;
							newe2->updateentry(entries[k]);
						  }
					   }
					   else
					   {
							if(j<node2->maxentries)
							{

							node2->entries[j]=entries[k];
							entries[k]->inleaf=node2->isleaf;
							node2->entries[j]->nodeparent=node2;
							j++;
							node2->countentries=j;
							newe2->updateentry(entries[k]);
						   }
						   else
						   {
							   node1->entries[i]=entries[k];
							entries[k]->inleaf=node1->isleaf;
							node1->entries[i]->nodeparent=node1;
							i++;
							node1->countentries=i;
							newe1->updateentry(entries[k]);
						   }
					   }
				}

				k++;
			}

			if(dtw(newe1,entry)<=dtw(newe2,entry))
					   {

						if(i<node1->maxentries)
						  {
							node1->entries[i]=entry;
							entry->inleaf=node1->isleaf;
							node1->entries[i]->nodeparent=node1;
							i++;
							node1->countentries=i;
							newe1->updateentry(entry);
						  }
						  else
						  {
								node2->entries[j]=entry;
							entry->inleaf=node2->isleaf;
							node2->entries[j]->nodeparent=node2;
							j++;
							node2->countentries=j;
							newe2->updateentry(entry);
						  }


					   }
					   else
					   {
							if(j<node2->maxentries)
							{
								node2->entries[j]=entry;
							entry->inleaf=node2->isleaf;
							node2->entries[j]->nodeparent=node2;
							j++;
							node2->countentries=j;
							newe2->updateentry(entry);
						   }
						   else
						   {
							   node1->entries[i]=entry;
							entry->inleaf=node1->isleaf;
							node1->entries[i]->nodeparent=node1;
							i++;
							node1->countentries=i;
							newe1->updateentry(entry);
						   }

					   }

			node1->countentries=i;
			node2->countentries=j;


			i=0;
			while(i<node1->countentries)
			{
				node1->entries[i]->nodeparent=node1;
				node1->entries[i]->inleaf=node1->isleaf;
				i++;
			}
			j=0;
			while(j<node2->countentries)
			{
				node2->entries[j]->nodeparent=node2;
				node2->entries[j]->inleaf=node2->isleaf;
				j++;
			}

			newe1->child=node1;
			node1->parent=newe1;
			newe2->child=node2;
			node2->parent=newe2;

		}
		int replaceentryearlynode(cfnode* earlynode,cfentry *removeentry, cfentry *newentry)
		{
			int i=0; int hd=0;
			while(i<earlynode->countentries)
			{
				if(earlynode->entries[i]==removeentry)
					{
						earlynode->entries[i]=newentry;
						earlynode->entries[i]->nodeparent=earlynode;
						earlynode->entries[i]->inleaf=earlynode->isleaf;
						hd=i;
						removeentry->child=NULL; removeentry->nodeparent=NULL;
						//free(removeentry->child);

						cfentry *x=earlynode->parent;
						while(x!=NULL)
						{
							x->updateentryfull();
							x=x->nodeparent->parent;
						}
						break;
					}
				i++;
			}
		return hd;
		}
		bool addentryearlynode(cfnode* earlynode,cfentry* newe2,cfentry* newe1, int hd)
		{

			if(earlynode->countentries==earlynode->maxentries)
			{
				if(earlynode->split(newe2))
					return true;
				else
					return false;
			}
			else
			{
					int s=earlynode->countentries;
					earlynode->entries[s]=newe2;
					newe2->nodeparent=earlynode;
					newe2->inleaf=earlynode->isleaf;
					earlynode->countentries++;

					cfentry *x=earlynode->parent;
					 while(x!=NULL)
					 {
						 x->updateentryfull();
						 x=x->nodeparent->parent;
					 }

			return true;
			}
		}

	};

    void cfentry::updateentry(cfentry* e)					//just update this entries with the new entry into a new entry
    {
        if(!(inleaf))
            {
                modifycomparepoints(this, e);
            }
            else
            {
                cout<<" CHECK ONCE AGAIN HERE "<<endl;
            }
    }

    void cfentry::updateentrysame(cfentry* e)				//updates the entry with the given entry where points are already present
    {
        int ui=0; int j=0;
        ncompare=0;
        free(comparepoints);
        maxcmp=e->ncompare+1;
        comparepoints = new cfpoint*[maxcmp];
        while(ui<e->getncompare())
        {
            comparepoints[ncompare]=new cfpoint(e->getcomparepoints(ui)->getxsum());

            ncompare++;

            if(ncompare==maxcmp)
                    {
                        cfpoint **newcompare;
                        maxcmp=2*maxcmp;
                        newcompare=new cfpoint *[maxcmp];
                        j=0;
                        while(j<ncompare)
                        {
                            newcompare[j]=comparepoints[j];
                            j++;
                        }
                        free(comparepoints);
                        comparepoints=newcompare;
                    }

            ui++;
        }
    }

    void cfentry::updateentryfull()							//updating the entries of both the node and its children
    {
        int i=0;
        updateentrysame(child->entries[i]);
        i++;
        while(i<child->countentries)
        {
            updateentry(child->entries[i]);
            i++;
        }
    }


	class cftree											//class for the entire clustering feature tree
	{
	private:
		cfnode *root;
	public:
		friend class cfpoint;
		friend class cfentry;
		friend class cfnode;
		cftree()							//constructor
		{
			root=NULL;
		}
		cfnode* getroot()					//get the root node
		{
			return root;
		}
		bool addnewpoint(cfentry *p)		//add this new entry to the node
		{
		   if(root==NULL)
		   {
			   root=new cfnode;
			   root->countentries=0;
			   root->maxentries=l;
			   root->isleaf=true;
			   if(l!=b) root->entries= new cfentry* [l];
			   root->next=NULL; root->previous=NULL; root->parent=NULL;

			   if(root->addentry(p))
			   {
					leafliststart=root;
					checkmemorysize();

					return true;
			   }
		   }
		   else if(root->isleaf)
			{

			   if(root->addentry(p))
				{
					checkmemorysize();
					return true;
				}
			   else
			   {
					if(rootsplitentry!=NULL)
					{
						rootsplit();
						checkmemorysize();
						return true;
					}
			   }
			}
			else
			{
				cfnode *x;
				x=root; cfentry* y;

				while(!(x->isleaf))
				{
					y=getclosestentry(p,x);
					x=y->child;
				}

			   if(x->addentry(p))
			   {
				   checkmemorysize();
				   return true;
			   }
			   else
			   {
				   if(rootsplitentry!=NULL)
					{
						rootsplit();
						checkmemorysize();

					return true;
					}
			   }

			}
		return false;
		}
		cfentry* getclosestentry(cfentry* p, cfnode* x)	//get the closest entry in that node to the entry
		{
			 if(x->getcountentries()!=0)
				{
					int i=0, t=0; double dist=0;
					dist=dtw(p,x->entries[i]);
					i++;
					while(i<x->countentries)
					{
						if(dist>dtw(p,x->entries[i]))
						{   dist=dtw(p,x->entries[i]);
							t=i;
						}
						i++;
					}
				return x->entries[t];
				}
			 else
				return NULL;
		}
		void rootsplit()					//split the root node in the birch tree
		{
			cfentrypair* farthestpair= root->getfarthestentries();
			cfentry * newe1=new cfentry;
			cfnode* node1=new cfnode;
			cfentry* newe2=new cfentry;
			cfnode* node2=new cfnode;

			newe1->child=node1;
			node1->parent=newe1;
			newe2->child=node2;
			node2->parent=newe2;

			node1->maxentries=root->maxentries;
			node2->maxentries=root->maxentries;
			node1->entries = new cfentry *[root->maxentries];
			node2->entries = new cfentry *[root->maxentries];

			node1->isleaf=root->isleaf;
			node2->isleaf=root->isleaf;
			newe1->inleaf=false;
			newe2->inleaf=false;
			int i=0;

			if(root->isleaf)
			{
				node1->previous=root->previous;
				if(root->previous)
					root->previous->next=node1;
				else
					leafliststart=node1;
				node1->next=node2;
				node2->next=root->next;
				node2->previous=node1;
				if(root->next)
					root->next->previous=node2;
				root->next=NULL;
				root->previous=NULL;
			}

				root->redistributeentries(farthestpair,newe1,newe2,node1,node2,rootsplitentry);

				rootsplitentry=NULL;
				cfnode* newroot= new cfnode;
				newroot->isleaf=false;

				i=0;

				newroot->entries[i]=newe1;
				newroot->entries[i]->nodeparent=newroot;

				i++;
				newroot->entries[i]=newe2;
				newroot->entries[i]->nodeparent=newroot;

				i++;
				newroot->countentries=i;

				newe1->child=node1;
				node1->parent=newe1;
				newe2->child=node2;
				node2->parent=newe2;
				i=0;
				while(i<newroot->countentries)
				{
					newroot->entries[i]->nodeparent=newroot;
					newroot->entries[i]->inleaf=false;
					i++;
				}
				root->parent=NULL; root->previous=NULL; root->next=NULL;
				root=newroot;

		}
		void checkmemorysize()				//checks memory and when to split
		{
		   cfnode* k; int n=0;
			/**k=leafliststart; int co=0;

			while(k!=NULL)
			{
				n++; co=co+k->countentries; k=k->next;
			}

			if(n>100)
			{
				cout<<" REBUILDING, N: "<<n<<" CO: "<<co<<endl;
				updatethreshold();
				root=rebuildtree(co);
				//cout<<" Re-checking "<<endl;
				//checkmemorysize();
			}
			**/
			k=leafliststart; int co=0;

			while(k!=NULL)
			{
				n++; co=co+k->countentries; k=k->next;
			}
			 k=root;

			totalnodes=0;
			totalentries=0;
			totalpoint=0;
			totalsize(k);

			 n=0;
			k=leafliststart;
			while(k!=NULL)
			{
				 k=k->next;
				n++;
			}

			long memory= n*(4*l + 32) + (totalnodes-n)*(4*b + 32) + totalentries*(8*dim*2 + 36+12) + totalpoint*8;
			if(memory>(1.5*1024*1024))
			{
			   // cout<<" MEMORY : "<<memory<<" GIVEN : 1,572,864"<<" total nodes, leaf, entries : "<<totalnodes<<"  "<<n<<"  "<<totalentries<<endl;
				updatethreshold();
				//cout<<"REBUILDING TREE "<<endl;
				cout<<"       "<<n<<"  "<<co<<endl;
				root=rebuildtree(co);
			}

		}
		void totalfree(cfnode* l)			//free an entire node
		{
			int i=0;
			cfnode *n; cfnode *x;
			n=leafliststart;
			while(n!=NULL)
			{
				x=n->next;
				n->parent->child=NULL;
				free(n);
				n=x;
			}

			if(l)
			{
				if(!(l->isleaf))
				{
					i=0;
					while(i<l->countentries)
					{
						if(l->entries[i]->child)
							{
									totalfree(l->entries[i]->child);
							}
						else ;
						free(l->entries[i]);

						i++;
					}

					free(l->entries);
					free(l);
				}
				else
				{
					l->parent->child=NULL;
					free(l);
				}
			}

		}
		void totalsize(cfnode* l)			//calculate the size of the node
		{
			int i=0;
			totalnodes=totalnodes+1;
			totalentries=totalentries+l->countentries;
			while(i<l->countentries)
			{

				if(l->entries[i]->child)
					totalsize(l->entries[i]->child);
				else
					totalpoint=totalpoint + l->entries[i]->getn();
				i++;
			}
		}
		void updatethreshold()				//update the threshold based on the heuristic given in the paper
		{

			double avgdist=0, newthreshold=0; int n=0;
			cfnode* l; int i=0;
			l=leafliststart;
			while(l!=NULL)
			{
				i=0;
				while(i<l->countentries)
				{
					avgdist+=getclosestdistance(l,i);
					if(getclosestdistance(l,i)!=0) n++;
					i++;
				}
				l=l->next;
			}
			l=leafliststart; i=0;
			while(l!=NULL)
			{
				l=l->next; i++;
			}

			if(n>0)
				newthreshold = avgdist/n;

			if(newthreshold <= threshold)
			   { if(n!=0) newthreshold=threshold+ threshold/n ;
				  else newthreshold=2*threshold;
			   cout<<" threshold less. so updated other way "<<endl; }

			threshold=newthreshold;

			cout<<" threshold : "<<setprecision(30)<<threshold<<endl;

		}
		double getclosestdistance(cfnode*l, int i)		//get the distance between closest entries in the node for threshold update
		{
			int j=0; double mindist=0;
			mindist=dtw(l->entries[i],l->entries[j]);
			j++;
			if(l->countentries>1)
			{
				if(mindist==0)
					{  mindist=dtw(l->entries[i],l->entries[j]);
					j++;
					}
				while(j<l->countentries)
				{   if(j!=i)
					{
						if(mindist>dtw(l->entries[i],l->entries[j]))
								mindist=dtw(l->entries[i],l->entries[j]);
					}
					j++;
				}

			return mindist;
			}
			else
				return 0;
		}
		cfnode* rebuildtree(int co)					//rebuild the tree after updating the threshold
		{
			cftree* newt= new cftree;
			cfnode* l=leafliststart; int i=0; int h=0;

			cfentry **entries; entries=new cfentry *[co];
			while(l!=NULL)
			{
				i=0;
				while(i<l->countentries)
				{
					entries[h]=l->getentries(i);
					h++;
					i++;
				}
				l=l->next;
			}

			totalfree(root);
				i=0;
				while(i<co)
				{
					h=0;

					cfentry* e;
					e=entries[i];
					newt->addnewpoint(e);

					i++;
				}
				if(i==co) cout<<" FINAL UPBUILD DONE...!!!   "<<i<<endl;

			return newt->getroot();
		}

	};

    ofstream foutput;										//outputs the necessary files for that query
    void getfiles(cfentry *y)
    {
        int i=0;
        while(i<y->getn())
        {
            foutput<<"    "<<y->getid(i)<<endl;
            i++;
        }
    }
    cfentry* getclosestentry(cfentry* p, cfnode* x)			//get the closest entry for that entry in that node
    {
         if(x->getcountentries()!=0)
            {
                int i=0, t=0; double dist=0;
                dist=dtw(p,x->entries[i]);
                i++;
                while(i<x->countentries)
                {
                    if(dist>dtw(p,x->entries[i]))
                    {   dist=dtw(p,x->entries[i]);
                        t=i;
                    }
                    i++;
                }
            return x->entries[t];
            }
         else
            return NULL;
    }
    void searchtimitquery(cfentry* query, cfnode* root)		//searches the tree for a given query and outputs the necessary file names
    {
            cfnode *x;
            x=root; cfentry* y;

            while(!(x->isleaf))
            {
                y=getclosestentry(query,x);
                x=y->child;
            }
            y=getclosestentry(query,x);
            foutput<<" "<<y->getn()<<endl;
            getfiles(y);
    }


int main()				//main function builds the indexing module and outputs for their queries
{


    int n; ifstream fin;
    cout<<" ===================================================================="<<endl;
    cout<<"                                                                     "<<endl;
    cout<<"                        IIIT HYDERABAD                               "<<endl;
    cout<<"              Language Technologies Research Center                  "<<endl;
    cout<<"                                                                     "<<endl;
    cout<<"                  IMPLEMENTING BIRCH ALGORITHM                "<<endl;
    cout<<"                                                                     "<<endl;
    cout<<" ===================================================================="<<endl;
    cout<<" Taking inputs here..."<<endl;
    cout<<" enter dimension value: ";
    cin>>dim;
    cout<<" enter the value of b (no . of entries in non-leaf node) : ";
    cin>>b;
    cout<<" enter the value of l (no . of entries in leaf node)  : ";
    cin>>l;
    cout<<" enter threshold value for violation (typically, initially zero): ";
    cin>>threshold;
    cout<<" ===================================================================="<<endl;
    cout<<" TAKING INPUT FROM FILE "<<endl;
    cout<<" enter yes or no :  "<<endl; string c; 

    fin.open("train.txt"); double *p; p=new double[dim]; string phone;
    int i=0,k=0; long int q=0; string s; long int jc=0;

    cftree *tree1= new cftree;
    while(!fin.eof())
    {

        fin>>n; q++; jc++;
        cfentry* entry = new cfentry[1];
        i=0;

        while(i<n)
        {
            k=0;
            while(k<dim)
            {
                fin>>p[k]; k++;
            }
            fin>>s; jc++;			//s contains the filename
									//p contains the feature vector
            cfpoint* point =new cfpoint(q,s,p);		//creating new cfpoint
            entry->addthispoint(point,q);

            i++;
        }

        if(tree1->addnewpoint(entry)) cout<<"  SUCCESS "<<q<<"  line:-  "<<jc<<"  phone :- "<<s<<endl;

        cfnode* k=leafliststart; int co=0;
        n=0;
        while(k!=NULL)
        {
            n++; co=co+k->getcountentries(); k=k->getnext();
        }

        cout<<" NODES : "<<n<<" TOTAL ENTRIES : "<<co<<endl;
    }

    cout<<" Finished reading input till here "<<endl;
    cout<<" The Clustering Feature Tree is ready...!! "<<endl;
    cout<<" shall we continue... ? ";
    char cc; cin>>cc;

    fin.close();						//tests the output queries from the files
    fin.open("test.txt");
    foutput.open("output.txt");
    q=0; jc=0;

    while(!fin.eof())
    {
        fin>>n; q++; jc++; 
        cfentry *entry = new cfentry;
        i=0;
        while(i<n)
        {
            k=0;
            while(k<dim)
            {
                fin>>p[k]; k++;
            }
            fin>>s; jc++;
            cfpoint* point = new cfpoint(q,s,p);
            entry->addthispoint(point,q);
            i++;
        }
        foutput<<q<<endl;
        searchtimitquery(entry,tree1->getroot());
        cout<<"  SUCCESS "<<q<<"  line:-  "<<jc<<endl;
    }

    fin.close();
    foutput.close();

    return 0;
}



/**

The code has been written as a part of my internship by Murali Raghu Babu B, IIT Guwahati
Contact me at muraliraghubabu1994@gmail.com or +91 8473994919



**/




